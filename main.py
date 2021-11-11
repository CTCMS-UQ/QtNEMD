#!/usr/bin/python3
from PyQt5 import QtWidgets
import sys 
import MainWindow


VERSION_NUMBER=0.02

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)

    window = MainWindow.MainWindow()
    window.show()

    sys.exit(app.exec_())
