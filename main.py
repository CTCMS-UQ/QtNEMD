#!/usr/bin/python3
# Import the MD Fortran routines first
import TTCF
import InputManager
import PlotManager
import sys 

from PyQt5 import QtWidgets, QtCore, QtGui
import pyqtgraph as pg

from threading import Thread

sys.path.append("GUI-resources")
from ui_mainwindow import Ui_MainWindow

VERSION_NUMBER=0.01

def thread_func(nsteps, iflag):
    TTCF.md(nsteps, iflag)

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
        self.ui.lj_spinbox.editingFinished.connect(self.update_parameters)
        self.ui.field_spinbox.editingFinished.connect(self.update_parameters)
        self.ui.e0_spinbox.editingFinished.connect(self.update_parameters)
        self.ui.tr_spinbox.editingFinished.connect(self.update_parameters)
        self.ui.density_spinbox.editingFinished.connect(self.update_parameters)

        # Checkbox to toggle NEMD field
        self.ui.nemd_checkbox.stateChanged.connect(self.toggle_ne_field)

        # Now initialise (but don't start) a timer to update the plot
        self.sim_timer = QtCore.QTimer()
        self.sim_timer.setInterval(8)
        # Connect the timer's "timeout" (finished) event to our update function
        self.sim_timer.timeout.connect(self.update_plot_data)

        # Now connect the file editing dialog to our open and close menu buttons
        self.ui.open_input_file.triggered.connect(self.open_input_file)
        self.ui.save_input_file.triggered.connect(self.save_to_file)
        self.ui.plot_action.triggered.connect(self.open_new_plot)
        self.ui.actionQuit.triggered.connect(self.clean_exit)
        finish = QtWidgets.QAction("Quit", self)
        finish.triggered.connect(self.closeEvent)

    def initialise_simulation(self):
        TTCF.setup() 

        # Only plot the fluid particles for now
        self.x = TTCF.coord.x[:TTCF.nopart.npart]
        self.y = TTCF.coord.y[:TTCF.nopart.npart]
        
        # Spin up a thread which will do MD steps
        self.md_thread = Thread(target = thread_func, args=(1,0))

        # Fix the X and Y ranges so they don't constantly shift throughout the simulation
        self.ui.plot_window.clear()
        self.ui.plot_window.setXRange(0, TTCF.parm.cubex)
        self.ui.plot_window.setYRange(0, TTCF.parm.cubey)

        self.ui.plot_window.setBackground('w')
        self.data =  self.ui.plot_window.plot(self.x, self.y, pen=None, symbol = 'o')

        # Initialise the N, V and T labels
        self.ui.npart_label.setText(f"N particles = {TTCF.nopart.npart}")
        volume = TTCF.nopart.npart / TTCF.parm.drf
        self.ui.volume_label.setText(f"Volume = {volume:.2f}")
        self.ui.temperature_label.setText(f"Temperature = {TTCF.averg.temp:.2f}")

        # Finally, set the simulation controls to the correct value
        self.ui.lj_spinbox.setValue(TTCF.parm.kf)
        self.ui.field_spinbox.setValue(TTCF.parm.fe0)
        self.ui.tr_spinbox.setValue(TTCF.inener.tr)
        self.ui.e0_spinbox.setValue(TTCF.inener.e0)
        self.ui.density_spinbox.setValue(TTCF.parm.drf)

    def update_parameters(self):
        # Get the widget which sent this signal, as well as its new value
        sender = self.sender()
        value = sender.value()

        # Now change the appropriate simulation parameter
        if sender == self.ui.lj_spinbox:
            TTCF.parm.kf = value
            self.params.kf = value

        elif sender == self.ui.field_spinbox:
            TTCF.parm.fe0 = value
            self.params.fe0 = value

        # These spinboxes control initial parameters, and require the simulation to be restarted after
        # changing
        elif sender == self.ui.tr_spinbox:
            TTCF.inener.tr = value
            self.params.tr = value
            self.restart_sim(sender)

        elif sender == self.ui.e0_spinbox:
            TTCF.inener.e0 = value
            self.params.e0 = value
            self.restart_sim(sender)

        elif sender == self.ui.density_spinbox:
            TTCF.parm.drf = value
            self.params.drf = value
            self.restart_sim(sender)

        else:
            print("Unknown sender")
            pass
        # Finally, update the input values in the QTextBrowser widget
        self.param_str = self.params.format_params()
        self.ui.input_textbrowser.setPlainText(self.param_str)
        self.ui.input_textbrowser.repaint()

    def toggle_ne_field(self, state):
        # Toggles the nonequilibrium field (on or off) based on the status of nemd_checkbox
        if state == QtCore.Qt.Checked:
            self.params.do_nemd = True
        else:
            self.params.do_nemd = False

    ################################# Plotting routines ################################
    def update_plot_data(self):
        # Get the number of steps to run our MD for and run
        plot_step = TTCF.iparm.nplot

        self.tau += plot_step

        # Set up field strength based on whether or not we're doing NEMD
        if self.params.do_nemd:
            TTCF.parm.field = self.params.fe0
            iflag = 1
        else:
            TTCF.parm.field = 0.0
            iflag = 0
            
        # join() the MD thread since we can't update until it's finished the MD step
        self.md_thread.join()

        # Only plot the fluid particles for now
        self.x = TTCF.coord.x[:TTCF.nopart.npart]
        self.y = TTCF.coord.y[:TTCF.nopart.npart]

        self.timestep_update.emit(self.tau)
        # Restart the MD thread in the background while we update the plot
        self.md_thread = Thread(target = thread_func, args=(1,0))
        self.md_thread.start()
                                                                               
        self.data.setData(self.x, self.y)  # Update the data.

        # Update the temperature label
        self.ui.temperature_label.setText(f"Temperature = {TTCF.averg.temp:.2f}")
        volume = TTCF.nopart.npart / TTCF.parm.drf
        self.ui.volume_label.setText(f"Volume = {volume:.2f}")

        # Send the signal to update any open real-time plot windows

    def start_sim(self):
        """ Start the simulation.

            The timer has already been initialised and linked to the update function, so we only need to
            start the timer here."""
        self.sim_timer.start()

        self.ui.start_button.setEnabled(False)
        self.ui.pause_button.setEnabled(True)

        # Also want to disable the Npart spinbox, since it makes no sense to change the particle number
        # while the simulation is running
        self.ui.tr_spinbox.setEnabled(False)
        self.ui.e0_spinbox.setEnabled(False)
        
        # Finally, start the background computational (MD) thread. Check if it's running first so we
        # don't get conflicts
        try:
            self.md_thread.start()
        except RuntimeError:
            self.md_thread = Thread(target = thread_func, args=(1,0))

    def pause_sim(self):
        """ Pause the simulation.
            
            This is simplest to achieve by simply stopping the timer temporarily, so the plot stops
            updating. It will start back up again when the timer is restarted."""
        if self.sim_timer.isActive():
            self.sim_timer.stop()

        self.ui.start_button.setEnabled(True)
        self.ui.pause_button.setEnabled(False)

        # Join the worker thread so it isn't left hanging
        if(self.md_thread.is_alive()):
            self.md_thread.join()
            self.md_thread = Thread(target = thread_func, args=(1,0))

    def restart_sim(self, sender = None):
        """ Restart the simulation by stopping the timer and reinitialising parameters."""
        if self.sim_timer.isActive():
            self.sim_timer.stop()
        self.tau = 0

        # Join the worker thread so it isn't left hanging
        if(self.md_thread.is_alive()):
            self.md_thread.join()
            self.md_thread = Thread(target = thread_func, args=(1,0))

        TTCF.setup() 

        # Only plot the fluid particles for now
        self.x = TTCF.coord.x[:TTCF.nopart.npart]
        self.y = TTCF.coord.y[:TTCF.nopart.npart]
        
        # Fix the X and Y ranges so they don't constantly shift throughout the simulation
        self.ui.plot_window.setXRange(0, TTCF.parm.cubex)
        self.ui.plot_window.setYRange(0, TTCF.parm.cubey)

        self.ui.plot_window.setBackground('w')
        self.data.setData(self.x, self.y, pen=None, symbol = 'o')

        # Re-enable buttons which can't be changed while the simulation is running
        self.ui.tr_spinbox.setEnabled(True)
        self.ui.e0_spinbox.setEnabled(True)
        self.ui.start_button.setEnabled(True)
        self.initialise_simulation()

    ######################### I/O Control routines ####################################
    def open_input_file(self):
        self.sim_timer.stop()
        input_file, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Open File", "", "Input Files (*.in)")
        if input_file:
            self.params.read_from_file(input_file)
            # Update the input values in the QTextBrowser widget
            self.param_str = self.params.format_params()
            self.ui.input_textbrowser.setPlainText(self.param_str)
            self.ui.input_textbrowser.repaint()
            self.restart_sim()
            self.initialise_simulation()
            
    def save_to_file(self):
        self.sim_timer.stop()
        output_file, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save File", "", "Input Files (*.in)")
        if output_file:
            out_string = self.params.format_params()
            with open(output_file, 'w') as ofp:
                ofp.write(out_string)

    def open_new_plot(self):
        new_plot = PlotManager.PlotManager()
        new_plot.show()
        self.window_list.append(new_plot)

        # Now connect signals so the main window can communicate with the floating plot widget
        self.timestep_update.connect(new_plot.update)

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

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)

    window = MainWindow()
    window.show()

    sys.exit(app.exec_())
