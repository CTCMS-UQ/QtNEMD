# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'GUI-resources/mainwindow.ui'
#
# Created by: PyQt5 UI code generator 5.15.6
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1454, 727)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(MainWindow.sizePolicy().hasHeightForWidth())
        MainWindow.setSizePolicy(sizePolicy)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout(self.centralwidget)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout()
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.plot_window = PlotWidget(self.centralwidget)
        self.plot_window.setObjectName("plot_window")
        self.verticalLayout_2.addWidget(self.plot_window)
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.temperature_label = QtWidgets.QLabel(self.centralwidget)
        self.temperature_label.setObjectName("temperature_label")
        self.horizontalLayout_3.addWidget(self.temperature_label)
        self.volume_label = QtWidgets.QLabel(self.centralwidget)
        self.volume_label.setObjectName("volume_label")
        self.horizontalLayout_3.addWidget(self.volume_label)
        self.npart_label = QtWidgets.QLabel(self.centralwidget)
        self.npart_label.setObjectName("npart_label")
        self.horizontalLayout_3.addWidget(self.npart_label)
        self.verticalLayout_2.addLayout(self.horizontalLayout_3)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.start_button = QtWidgets.QPushButton(self.centralwidget)
        self.start_button.setObjectName("start_button")
        self.horizontalLayout.addWidget(self.start_button)
        self.pause_button = QtWidgets.QPushButton(self.centralwidget)
        self.pause_button.setObjectName("pause_button")
        self.horizontalLayout.addWidget(self.pause_button)
        self.restart_button = QtWidgets.QPushButton(self.centralwidget)
        self.restart_button.setObjectName("restart_button")
        self.horizontalLayout.addWidget(self.restart_button)
        self.verticalLayout_2.addLayout(self.horizontalLayout)
        self.horizontalLayout_2.addLayout(self.verticalLayout_2)
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalLayout_6 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_6.setObjectName("horizontalLayout_6")
        self.nemd_checkbox = QtWidgets.QCheckBox(self.centralwidget)
        self.nemd_checkbox.setObjectName("nemd_checkbox")
        self.horizontalLayout_6.addWidget(self.nemd_checkbox)
        self.twod_button = QtWidgets.QRadioButton(self.centralwidget)
        self.twod_button.setChecked(True)
        self.twod_button.setObjectName("twod_button")
        self.dimensionButtonGroup = QtWidgets.QButtonGroup(MainWindow)
        self.dimensionButtonGroup.setObjectName("dimensionButtonGroup")
        self.dimensionButtonGroup.addButton(self.twod_button)
        self.horizontalLayout_6.addWidget(self.twod_button)
        self.threed_button = QtWidgets.QRadioButton(self.centralwidget)
        self.threed_button.setObjectName("threed_button")
        self.dimensionButtonGroup.addButton(self.threed_button)
        self.horizontalLayout_6.addWidget(self.threed_button)
        self.verticalLayout.addLayout(self.horizontalLayout_6)
        self.horizontalLayout_4 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_4.setObjectName("horizontalLayout_4")
        self.formLayout = QtWidgets.QFormLayout()
        self.formLayout.setObjectName("formLayout")
        self.field_spinbox = QtWidgets.QDoubleSpinBox(self.centralwidget)
        self.field_spinbox.setMinimumSize(QtCore.QSize(129, 0))
        self.field_spinbox.setStatusTip("")
        self.field_spinbox.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.field_spinbox.setDecimals(2)
        self.field_spinbox.setMaximum(10.0)
        self.field_spinbox.setSingleStep(0.01)
        self.field_spinbox.setObjectName("field_spinbox")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.field_spinbox)
        self.field_label = QtWidgets.QLabel(self.centralwidget)
        self.field_label.setObjectName("field_label")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.field_label)
        self.lj_eps_spinbox = QtWidgets.QDoubleSpinBox(self.centralwidget)
        self.lj_eps_spinbox.setMinimumSize(QtCore.QSize(129, 0))
        self.lj_eps_spinbox.setStatusTip("")
        self.lj_eps_spinbox.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.lj_eps_spinbox.setDecimals(2)
        self.lj_eps_spinbox.setSingleStep(0.01)
        self.lj_eps_spinbox.setObjectName("lj_eps_spinbox")
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.lj_eps_spinbox)
        self.lj_eps_label = QtWidgets.QLabel(self.centralwidget)
        self.lj_eps_label.setObjectName("lj_eps_label")
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.lj_eps_label)
        self.horizontalLayout_4.addLayout(self.formLayout)
        self.formLayout_2 = QtWidgets.QFormLayout()
        self.formLayout_2.setObjectName("formLayout_2")
        self.density_spinbox = QtWidgets.QDoubleSpinBox(self.centralwidget)
        self.density_spinbox.setMinimumSize(QtCore.QSize(129, 0))
        self.density_spinbox.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.density_spinbox.setMaximum(10.0)
        self.density_spinbox.setSingleStep(0.01)
        self.density_spinbox.setObjectName("density_spinbox")
        self.formLayout_2.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.density_spinbox)
        self.density_label = QtWidgets.QLabel(self.centralwidget)
        self.density_label.setObjectName("density_label")
        self.formLayout_2.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.density_label)
        self.temp_spinbox = QtWidgets.QDoubleSpinBox(self.centralwidget)
        self.temp_spinbox.setMinimumSize(QtCore.QSize(129, 0))
        self.temp_spinbox.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.temp_spinbox.setDecimals(2)
        self.temp_spinbox.setMinimum(0.01)
        self.temp_spinbox.setMaximum(10.0)
        self.temp_spinbox.setSingleStep(0.01)
        self.temp_spinbox.setObjectName("temp_spinbox")
        self.formLayout_2.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.temp_spinbox)
        self.tr_label = QtWidgets.QLabel(self.centralwidget)
        self.tr_label.setObjectName("tr_label")
        self.formLayout_2.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.tr_label)
        self.horizontalLayout_4.addLayout(self.formLayout_2)
        self.verticalLayout.addLayout(self.horizontalLayout_4)
        self.label = QtWidgets.QLabel(self.centralwidget)
        self.label.setAlignment(QtCore.Qt.AlignCenter)
        self.label.setObjectName("label")
        self.verticalLayout.addWidget(self.label)
        self.g2_window = PlotWidget(self.centralwidget)
        self.g2_window.setObjectName("g2_window")
        self.verticalLayout.addWidget(self.g2_window)
        self.horizontalLayout_2.addLayout(self.verticalLayout)
        self.verticalLayout_3 = QtWidgets.QVBoxLayout()
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.textbrowser_label = QtWidgets.QLabel(self.centralwidget)
        self.textbrowser_label.setAlignment(QtCore.Qt.AlignCenter)
        self.textbrowser_label.setObjectName("textbrowser_label")
        self.verticalLayout_3.addWidget(self.textbrowser_label)
        self.input_textbrowser = QtWidgets.QTextBrowser(self.centralwidget)
        self.input_textbrowser.setMinimumSize(QtCore.QSize(250, 271))
        self.input_textbrowser.setLineWrapMode(QtWidgets.QTextEdit.NoWrap)
        self.input_textbrowser.setObjectName("input_textbrowser")
        self.verticalLayout_3.addWidget(self.input_textbrowser)
        self.horizontalLayout_2.addLayout(self.verticalLayout_3)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1454, 29))
        self.menubar.setObjectName("menubar")
        self.menuView = QtWidgets.QMenu(self.menubar)
        self.menuView.setObjectName("menuView")
        self.menuFile = QtWidgets.QMenu(self.menubar)
        self.menuFile.setObjectName("menuFile")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.plot_action = QtWidgets.QAction(MainWindow)
        self.plot_action.setObjectName("plot_action")
        self.open_input_file = QtWidgets.QAction(MainWindow)
        self.open_input_file.setObjectName("open_input_file")
        self.save_input_file = QtWidgets.QAction(MainWindow)
        self.save_input_file.setObjectName("save_input_file")
        self.actionQuit = QtWidgets.QAction(MainWindow)
        self.actionQuit.setObjectName("actionQuit")
        self.actionAbout = QtWidgets.QAction(MainWindow)
        self.actionAbout.setObjectName("actionAbout")
        self.actionStructure_plots = QtWidgets.QAction(MainWindow)
        self.actionStructure_plots.setObjectName("actionStructure_plots")
        self.menuView.addAction(self.plot_action)
        self.menuView.addAction(self.actionStructure_plots)
        self.menuFile.addAction(self.open_input_file)
        self.menuFile.addAction(self.save_input_file)
        self.menuFile.addAction(self.actionQuit)
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuView.menuAction())

        self.retranslateUi(MainWindow)
        self.actionQuit.triggered.connect(MainWindow.close) # type: ignore
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "QtNEMD"))
        self.temperature_label.setText(_translate("MainWindow", "Temperature = "))
        self.volume_label.setText(_translate("MainWindow", "Volume = "))
        self.npart_label.setText(_translate("MainWindow", "N particles = "))
        self.start_button.setText(_translate("MainWindow", "Start"))
        self.pause_button.setText(_translate("MainWindow", "Pause"))
        self.restart_button.setText(_translate("MainWindow", "Reset"))
        self.nemd_checkbox.setText(_translate("MainWindow", "Enable NEMD"))
        self.twod_button.setText(_translate("MainWindow", "2D"))
        self.threed_button.setText(_translate("MainWindow", "3D (Colour = depth)"))
        self.field_spinbox.setToolTip(_translate("MainWindow", "Field strength"))
        self.field_spinbox.setWhatsThis(_translate("MainWindow", "Spin box input for field strength"))
        self.field_spinbox.setAccessibleName(_translate("MainWindow", "Field strength"))
        self.field_spinbox.setAccessibleDescription(_translate("MainWindow", "Spin box input for field strength"))
        self.field_label.setText(_translate("MainWindow", "Field strength"))
        self.lj_eps_spinbox.setToolTip(_translate("MainWindow", "LJ strength"))
        self.lj_eps_spinbox.setWhatsThis(_translate("MainWindow", "Spin box input for LJ repulsion strength"))
        self.lj_eps_spinbox.setAccessibleName(_translate("MainWindow", "LJ strength"))
        self.lj_eps_spinbox.setAccessibleDescription(_translate("MainWindow", "Spin box input for LJ repulsion strength"))
        self.lj_eps_label.setText(_translate("MainWindow", "<html><head/><body><p>WCA ε</p></body></html>"))
        self.density_label.setText(_translate("MainWindow", "Fluid density"))
        self.temp_spinbox.setStatusTip(_translate("MainWindow", "Set the temperature"))
        self.temp_spinbox.setWhatsThis(_translate("MainWindow", "Spinbox to set the temperature of the simulation"))
        self.temp_spinbox.setAccessibleName(_translate("MainWindow", "Temperature control"))
        self.temp_spinbox.setAccessibleDescription(_translate("MainWindow", "Spinbox to set the temperature of the simulation (can only be used while the smulation is stopped)"))
        self.tr_label.setText(_translate("MainWindow", "Temperature"))
        self.label.setText(_translate("MainWindow", "<html><head/><body><p>g<span style=\" vertical-align:super;\">(2)</span>(r) structure factor</p></body></html>"))
        self.textbrowser_label.setText(_translate("MainWindow", "LAMMPS Input file (read only)"))
        self.input_textbrowser.setToolTip(_translate("MainWindow", "Input parameters for simulation"))
        self.input_textbrowser.setWhatsThis(_translate("MainWindow", "Input parameters for simulation in input file format"))
        self.menuView.setTitle(_translate("MainWindow", "View"))
        self.menuFile.setTitle(_translate("MainWindow", "File"))
        self.plot_action.setText(_translate("MainWindow", "New real-time plot"))
        self.plot_action.setToolTip(_translate("MainWindow", "Plotting dialog"))
        self.plot_action.setWhatsThis(_translate("MainWindow", "Select quantities to plot in real time"))
        self.plot_action.setShortcut(_translate("MainWindow", "Ctrl+P"))
        self.open_input_file.setText(_translate("MainWindow", "Read parameters from file"))
        self.open_input_file.setShortcut(_translate("MainWindow", "Ctrl+O"))
        self.save_input_file.setText(_translate("MainWindow", "Save parameters to file"))
        self.save_input_file.setShortcut(_translate("MainWindow", "Ctrl+S"))
        self.actionQuit.setText(_translate("MainWindow", "Quit"))
        self.actionQuit.setShortcut(_translate("MainWindow", "Ctrl+Q"))
        self.actionAbout.setText(_translate("MainWindow", "About"))
        self.actionStructure_plots.setText(_translate("MainWindow", "Structure plots"))
from pyqtgraph import PlotWidget
