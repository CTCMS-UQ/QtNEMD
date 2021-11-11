from PyQt5 import QtWidgets, QtCore
import pyqtgraph as pg
import sys 
import os
import InputManager
import LammpsDriver

sys.path.append("GUI-resources")
from ui_mainwindow import Ui_MainWindow

class MainWindow(QtWidgets.QMainWindow):
    # Custom signal to communicate with real-time plot widgets
    # This needs to go here and not in the constructor. See: 
    #https://stackoverflow.com/questions/2970312/pyqt4-qtcore-pyqtsignal-object-has-no-attribute-connect
    # for an explanation
    timestep_update = QtCore.pyqtSignal(int)

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        # Set the OMP_NUM_THREADS environment variable so LAMMPS doesn't complain
        os.environ['OMP_NUM_THREADS'] = str(1)
        self.md = LammpsDriver.MDInterface()
        # Timestep
        self.tau = 0

        # Keep track of the open plotting windows
        self.window_list = []

        # Register the signal handlers
        self.register_handlers()
        
        # Start with the pause button disabled until we start the simulation
        self.ui.pause_button.setEnabled(False)

        # Initialise the simulation to its default values
        self.params = InputManager.InputManager()
        # Print the input values to the QTextBrowser widget
        self.param_str = self.params.format_params()
        self.ui.input_textbrowser.setPlainText(self.param_str)
        self.initialise_simulation()

    def register_handlers(self):
        # Set up all signal handlers not already setup in Qt Designer
        self.ui.start_button.clicked.connect(self.start_sim)
        self.ui.pause_button.clicked.connect(self.pause_sim)
        self.ui.restart_button.clicked.connect(self.restart_sim)
       
        # Add signals so spin-boxes change the underlying simulation parameters
        self.ui.lj_eps_spinbox.editingFinished.connect(self.update_parameters)
        self.ui.lj_sigma_spinbox.editingFinished.connect(self.update_parameters)
        self.ui.field_spinbox.editingFinished.connect(self.update_parameters)
        #self.ui.e0_spinbox.editingFinished.connect(self.update_parameters)
        self.ui.temp_spinbox.editingFinished.connect(self.update_parameters)
        self.ui.density_spinbox.editingFinished.connect(self.update_parameters)

        # Checkbox to toggle NEMD field
        self.ui.nemd_checkbox.stateChanged.connect(self.toggle_ne_field)

        # Now initialise (but don't start) a timer to update the plot
        self.sim_timer = QtCore.QTimer()
        self.sim_timer.setInterval(8)

        # Now do a longer timer for updating the UI elements (but not the plot)
        self.gui_timer = QtCore.QTimer()
        self.gui_timer.setInterval(60)
        # Connect the timer's "timeout" (finished) event to our update function
        self.sim_timer.timeout.connect(self.update_plot_data)
        self.gui_timer.timeout.connect(self.update_GUI_elements)

        # Now connect the file editing dialog to our open and close menu buttons
        self.ui.open_input_file.triggered.connect(self.open_input_file)
        self.ui.save_input_file.triggered.connect(self.save_to_file)
        self.ui.plot_action.triggered.connect(self.open_new_plot)
        self.ui.actionQuit.triggered.connect(self.clean_exit)
        finish = QtWidgets.QAction("Quit", self)
        finish.triggered.connect(self.closeEvent)

    def initialise_simulation(self):

        self.md.setup()
        # Only plot the fluid particles for now
        coords = self.lmp.numpy.extract_atom("x")
        x = coords[:,0]
        y = coords[:,1]

        # Fix the X and Y ranges so they don't constantly shift throughout the simulation
        self.ui.plot_window.clear()
        self.ui.plot_window.setXRange(0, self.lmp.get_thermo("lx"))
        self.ui.plot_window.setYRange(0, self.lmp.get_thermo("ly"))
        self.ui.plot_window.setBackground('w')
        self.pos_data = self.ui.plot_window.plot(x, y, pen=None, symbol = 'o')

        # Now do the g(2) radial-distribution function
        g2_compute = self.lmp.numpy.extract_compute("g2", lammps.LMP_STYLE_GLOBAL, lammps.LMP_TYPE_ARRAY)
        r = g2_compute[:,0]
        g2 = g2_compute[:,1]
        self.ui.g2_window.setXRange(0, max(r)+1)
        self.ui.g2_window.setYRange(0, max(g2)+1)
        self.ui.g2_window.setBackground('w')
        self.g2_data = self.ui.g2_window.plot(r, g2, color='k')
        

        # Initialise the N, V and T labels
        self.ui.npart_label.setText(f"N particles = {npart}")
        area = self.lmp.get_thermo("lx")*self.lmp.get_thermo("ly")
        self.ui.volume_label.setText(f"Area = {area:.2f}")
        self.ui.temperature_label.setText(f"Temperature = {self.params.temp}")

        # Finally, set the simulation controls to the correct value
        self.ui.lj_eps_spinbox.setValue(self.params.eps)
        self.ui.lj_sigma_spinbox.setValue(self.params.sigma)
        self.ui.field_spinbox.setValue(self.params.flowrate)
        self.ui.temp_spinbox.setValue(self.params.temp)
        # Need to get LAMMPS to compute the kinetic energy
        #self.ui.e0_spinbox.setValue(TTCF.inener.e0)
        self.ui.density_spinbox.setValue(self.params.reduced_density)

    def update_parameters(self):
        # Get the widget which sent this signal, as well as its new value
        sender = self.sender()
        value = sender.value()

        # Now change the appropriate simulation parameter
        if sender == self.ui.lj_eps_spinbox:
            self.params.eps = value
            self.params.update_parameters(self.lmp)

        if sender == self.ui.lj_sigma_spinbox:
            self.params.sigma = value
            self.params.update_parameters(self.lmp)

        elif sender == self.ui.field_spinbox:
            self.params.flowrate = value
            self.params.update_parameters(self.lmp)

        # These spinboxes control initial parameters, and require the simulation to be restarted after
        # changing
        elif sender == self.ui.temp_spinbox:
            self.params.temp = value
            self.params.reset_and_update(self.lmp)

        elif sender == self.ui.density_spinbox:
            self.params.reduced_density = value
            self.params.reset_and_update(self.lmp)

        else:
            print("Unknown sender")
            pass
        ## Finally, update the input values in the QTextBrowser widget
        self.param_str = self.params.format_params()
        self.ui.input_textbrowser.setPlainText(self.param_str)
        self.ui.input_textbrowser.repaint()

    def toggle_ne_field(self, state):
        # Toggles the nonequilibrium field (on or off) based on the status of nemd_checkbox
        if state == QtCore.Qt.Checked:
            self.params.toggle_nemd(self.lmp)
        else:
            self.params.toggle_nemd(self.lmp)

        ## Finally, update the input values in the QTextBrowser widget
        self.param_str = self.params.format_params()
        self.ui.input_textbrowser.setPlainText(self.param_str)
        self.ui.input_textbrowser.repaint()


    ################################# Plotting routines ################################
    def update_plot_data(self):
        # First, run an MD timestep
        #self.lmp.command(f"run 1 pre no post no")
        self.lmp.command(f"run 1")

        # Get the current timestep
        self.tau = self.lmp.extract_global("ntimestep")

        # Only plot the fluid particles for now
        x = coords[:,0]
        y = coords[:,1]
        self.pos_data.setData(x, y)  # Update the data.

        # Send a signal that we've moved forward a timestep. This is currently useless, but will get
        # used to synchronise other real-time plots
        self.timestep_update.emit(self.tau)

        # Now do the g(2) radial-distribution function
        g2_compute = self.lmp.numpy.extract_compute("g2", lammps.LMP_STYLE_GLOBAL, lammps.LMP_TYPE_ARRAY)
        r = g2_compute[:,0]
        g2 = g2_compute[:,1]
        self.ui.g2_window.setXRange(0, max(r)+1)
        self.ui.g2_window.setYRange(0, max(g2)+1)
        self.g2_data.setData(r, g2)

    def update_GUI_elements(self):
        # Update the temperature and volume labels
        temp = self.lmp.get_thermo("temp")
        vol = self.lmp.get_thermo("vol")
        self.ui.temperature_label.setText(f"Temperature = {temp:.2f}")
        self.ui.volume_label.setText(f"Volume = {vol:.2f}")

    def start_sim(self):
        """ Start the simulation.

            The timer has already been initialised and linked to the update function, so we only need to
            start the timer here."""
        self.sim_timer.start()
        self.gui_timer.start()

        self.ui.start_button.setEnabled(False)
        self.ui.pause_button.setEnabled(True)

        # Also want to disable the Npart spinbox, since it makes no sense to change the particle number
        # while the simulation is running
        self.ui.temp_spinbox.setEnabled(False)
        #self.ui.e0_spinbox.setEnabled(False)
        #self.ui.field_spinbox.setEnabled(False)
        self.ui.density_spinbox.setEnabled(False)
        self.ui.lj_eps_spinbox.setEnabled(False)
        self.ui.lj_sigma_spinbox.setEnabled(False)
        
        # Finally, run an MD timestep
        self.lmp.command(f"run 1 ")
        
    def pause_sim(self):
        """ Pause the simulation.
            
            This is simplest to achieve by simply stopping the timer temporarily, so the plot stops
            updating. It will start back up again when the timer is restarted."""
        if self.sim_timer.isActive():
            self.sim_timer.stop()
            self.gui_timer.stop()

        self.ui.start_button.setEnabled(True)
        self.ui.pause_button.setEnabled(False)

    def restart_sim(self, sender = None):
        """ Restart the simulation by stopping the timer and reinitialising parameters."""
        if self.sim_timer.isActive():
            self.sim_timer.stop()
            self.gui_timer.stop()
        self.tau = 0

        self.params.reset_and_update(self.lmp)

        # Calculate thermodynamic variables
        self.lmp.command("run 0")

        # Only plot the fluid particles for now
        coords = self.lmp.numpy.extract_atom("x")
        x = coords[:,0]
        y = coords[:,1]

        # Fix the X and Y ranges so they don't constantly shift throughout the simulation
        self.ui.plot_window.clear()
        self.ui.plot_window.setXRange(0, self.lmp.get_thermo("lx"))
        self.ui.plot_window.setYRange(0, self.lmp.get_thermo("ly"))

        self.ui.plot_window.setBackground('w')
        self.pos_data = self.ui.plot_window.plot(x, y, pen=None, symbol = 'o')

        # Now do the g(2) radial-distribution function
        g2_compute = self.lmp.numpy.extract_compute("g2", lammps.LMP_STYLE_GLOBAL, lammps.LMP_TYPE_ARRAY)
        r = g2_compute[:,0]
        g2 = g2_compute[:,1]
        self.ui.g2_window.clear()
        self.g2_data = self.ui.g2_window.plot(r, g2)

        # Re-enable buttons which can't be changed while the simulation is running
        self.ui.temp_spinbox.setEnabled(True)
        #self.ui.e0_spinbox.setEnabled(True)
        #self.ui.field_spinbox.setEnabled(True)
        self.ui.density_spinbox.setEnabled(True)
        self.ui.start_button.setEnabled(True)
        self.ui.lj_eps_spinbox.setEnabled(True)
        self.ui.lj_sigma_spinbox.setEnabled(True)

        self.initialise_simulation()

    ######################### I/O Control routines ####################################
    def open_input_file(self):
        pass
        #self.sim_timer.stop()
        #input_file, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Open File", "", "Input Files (*.in)")
        #if input_file:
        #    self.params.read_from_file(input_file)
        #    # Update the input values in the QTextBrowser widget
        #    self.param_str = self.params.format_params()
        #    self.ui.input_textbrowser.setPlainText(self.param_str)
        #    self.ui.input_textbrowser.repaint()
        #    self.restart_sim()
        #    self.initialise_simulation()
            
    def save_to_file(self):
        self.sim_timer.stop()
        self.gui_timer.stop()
        output_file, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save File", "", "Input Files (*.in)")
        if output_file:
            out_string = self.params.format_params()
            with open(output_file, 'w') as ofp:
                ofp.write(out_string)

    def open_new_plot(self):
        pass
        #new_plot = PlotManager.PlotManager()
        #new_plot.show()
        #self.window_list.append(new_plot)

        ## Now connect signals so the main window can communicate with the floating plot widget
        #self.timestep_update.connect(new_plot.update)

    ############################# Clean exit ################################
    def clean_exit(self):
        # Close all open windows when the main window is closed
        for window in self.window_list:
            window.close()
        
        self.close()

    def closeEvent(self, event):
        # Close all open windows when the main window is closed
        for window in self.window_list:
            window.close()
        
        event.accept()
