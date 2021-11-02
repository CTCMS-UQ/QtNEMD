#!/usr/bin/python
# Widget and logic to control real-time plotting
import sys
import TTCF
import numpy as np
from PyQt5 import QtWidgets, QtCore, QtGui
import pyqtgraph as pg

sys.path.append("GUI-resources")
from ui_floating_plot import Ui_PlotWidgetWindow

class PlotManager(QtWidgets.QMainWindow):
    
    def __init__(self, *args, **kwargs):
        super(PlotManager, self).__init__(*args, **kwargs)
        self.ui = Ui_PlotWidgetWindow()
        self.ui.setupUi(self)

        self.setup_plot()
        
        # Handle shut-downs properly
        self.ui.actionQuit.triggered.connect(self.close)

    def setup_plot(self):
        # Setup the basic plotting details for tracking TTCF
        self.x = []
        self.y = []
        self.z = []
        self.t = []

        # Array to hold the data to plot. Could be multidimensional depending on how many components (x,
        # y, z) the user requests
        self.to_plot = []
        
        self.ui.plot_window.setBackground('w')
 
        # Enable automatic scaling of the x-axis (time) and add a legend
        vbox = self.ui.plot_window.getViewBox()
        vbox.enableAutoRange(axis=vbox.XAxis)
        vbox.enableAutoRange(axis=vbox.YAxis)
        self.ui.plot_window.addLegend()

        # Now register handlers to switch between quantities
        self.ui.msd_button.clicked.connect(self.switch_to_msd)
        self.ui.ttcf_button.clicked.connect(self.switch_to_ttcf)

        # Do the initial setup
        if self.ui.msd_button.isChecked():
            self.plotter = MSDPlotter(parent=self)

        elif self.ui.ttcf_button.isChecked():
            self.switch_to_tcf()
        
    def switch_to_msd(self):
        self.plotter = MSDPlotter()

    def switch_to_ttcf(self):
        self.plotter = TCFPlotter()

    ############################### Custom slot for inter-widget communication ######################
    @QtCore.pyqtSlot(int)
    def update(self, tau):
        """ Updates the plot when the main widget sends the 'timestep_update(int)' signal."""

        self.ui.plot_window.clear()
        
        # Decide which quantities and coordinates to plot based on which checkboxes are ticked
        self.plotter.update(tau)
        
        if self.ui.x_checkbox.isChecked():
            self.ui.plot_window.plot(self.plotter.t, self.plotter.x, name='x')

        if self.ui.y_checkbox.isChecked():
            pen = pg.mkPen(color=(0, 0, 0))
            self.ui.plot_window.plot(self.plotter.t, self.plotter.y, name='y', pen=pen)

        if self.ui.z_checkbox.isChecked():
            pen = pg.mkPen(color=(0, 0, 0), style=QtCore.Qt.DashLine)
            self.ui.plot_window.plot(self.plotter.t, self.plotter.z, name='z', pen=pen)

class GenericPlotter:
    """ Helper class to handle plotting of various quantities."""
    def __init__(self, parent):
        self.parent = parent

        # Get the initial x, y, z and t arrays to empty lists. They'll be updated once we start
        # plotting
        self.x = []
        self.y = []
        self.z = []
        self.t = []

    def update(self):
        pass

class MSDPlotter(GenericPlotter):
    """ Helper class to plot mean-squared displacement."""
    
    def update(self, tau):
        # First, update the timestep
        self.t.append(tau*TTCF.simul.delta)

        # Now get the latest values of the x,y,z components from the TTCF module
        # Need to update all three arrays so they maintain the same shape. Otherwise, we can't switch 
        # between arrays once the plot has started
        self.x.append(float(TTCF.count.msdx))
        self.y.append(float(TTCF.count.msdy))
        self.z.append(float(TTCF.count.msdz)) 

class TCFPlotter:
    """ Helper class to handle plotting mean-squared displacement."""
    def __init__(self, parent):
        self.parent = parent

        # Get the initial x, y, z and t arrays to empty lists. They'll be updated once we start
        # plotting
        self.x = []
        self.y = []
        self.z = []
        self.t = []

    def update(self, tau):
        # First, update the timestep
        #self.t.append(tau*TTCF.simul.delta)

        # Now get the latest values of the x,y,z components from the TTCF module
        # Need to update all three arrays so they maintain the same shape. Otherwise, we can't switch 
        # between arrays once the plot has started
        #self.x.append(float(TTCF.count.tcfx))
        #self.y.append(float(TTCF.count.tcfy))
        #self.z.append(float(TTCF.count.tcfz)) 
        pass
