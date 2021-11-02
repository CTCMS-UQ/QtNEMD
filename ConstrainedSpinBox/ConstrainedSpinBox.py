#!/usr/bin/python
import sys
import os

from PyQt5.QtWidgets import (QApplication, QWidget, QSpinBox, QLabel)
from PyQt5.QtCore import Qt

class ConstrainedSpinBox(QSpinBox):
    """ Custom subclass of QSpinBox which is constrained to only take values from a finite list.
    
        The list of allowed values can be non-contiguous, but is currently hard-coded.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Keep track of the index we're currently at in the list of allowed values. This will be
        # incremented or decremented by stepBy
        self.index = 0

        # Make a list of all allowed values for NPart = 4*N^3
        self.allowed_values = [(4*N**3) for N in range(2,10)]
        self.allowed_values.sort()
        self.setRange(self.allowed_values[0], self.allowed_values[-1])
    
    def stepBy(self, steps):
        """ Overrides the base class's method to step through the list of allowed values."""

        # Clamp the step size to stay within the bounds of the list of allowed values
        raw_index = self.index + steps
        clamped_index = max(0, min(raw_index, len(self.allowed_values)-1))
        self.index = clamped_index
        self.setValue(self.allowed_values[self.index])

    def setValue(self, val):
        """ Overrides the base class's method to make sure that the index keeps pace with new values."""
        super().setValue(val)
        
        # Find the index of the new value
        try:
            self.index = self.allowed_values.index(val)
        except ValueError:
            print("Error: invalid value passed to ConstrainedSpinbox")
            #raise ValueError
