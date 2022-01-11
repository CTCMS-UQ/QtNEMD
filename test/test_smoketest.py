#!/usr/bin/python3
import unittest
import sys
sys.path.append("../")
import FortranDriver
import numpy as np

class DriverSmokeTest(unittest.TestCase):
    """ Set of "smoke tests" to check the driver for the Fortran backend runs at all.
    
    These tests are somewhat unsophisticated and are not strictly unit tests, 
    but provide value because if they fail then there's no point in continuing
    to test."""
    def setUp(self):
        """ Initialise the Fortran backend."""
        self.md = FortranDriver.MDInterface()
        self.md.setup()
    
    def test_md_short(self):
        """ Test a single step of the MD routine. """
        x_start, y_start, z_start = self.md.x.copy(), self.md.y.copy(), self.md.z.copy()
        self.md.run(1)

        # Particle position should be different to when we start
        x_end, y_end, z_end = self.md.x.copy(), self.md.y.copy(), self.md.z.copy()
        self.assertFalse(np.all(x_start == x_end))
        self.assertFalse(np.all(y_start == y_end))
        self.assertFalse(np.all(z_start == z_end))

        # Particles should not be outside the box bounds (0 < x < md.box_bounds[0])
        # (x,y,z) >= 0
        self.assertTrue(np.all(self.md.x >= 0))
        self.assertTrue(np.all(self.md.y >= 0))
        self.assertTrue(np.all(self.md.z >= 0))
        # (x,y,z) < (xmax,ymax,zmax)
        self.assertTrue(np.all(self.md.x < self.md.box_bounds[0]))
        self.assertTrue(np.all(self.md.y < self.md.box_bounds[0]))
        self.assertTrue(np.all(self.md.z < self.md.box_bounds[0]))
        # Basic thermodynamic properties
        self.assertGreater(self.md.temp, 0.0)

    def test_md_long(self):
        """ Test many steps of the MD routine"""
        x_start, y_start, z_start = self.md.x.copy(), self.md.y.copy(), self.md.z.copy()
        self.md.run(1000)
        x_end, y_end, z_end = self.md.x.copy(), self.md.y.copy(), self.md.z.copy()
        self.assertFalse(np.all(x_start == x_end))
        self.assertFalse(np.all(y_start == y_end))
        self.assertFalse(np.all(z_start == z_end))

        # Particles should not be outside the box bounds (0 < x < md.box_bounds[0])
        # (x,y,z) >= 0
        self.assertTrue(np.all(self.md.x >= 0))
        self.assertTrue(np.all(self.md.y >= 0))
        self.assertTrue(np.all(self.md.z >= 0))
        # (x,y,z) < (xmax,ymax,zmax)
        self.assertTrue(np.all(self.md.x < self.md.box_bounds[0]))
        self.assertTrue(np.all(self.md.y < self.md.box_bounds[0]))
        self.assertTrue(np.all(self.md.z < self.md.box_bounds[0]))

        # Basic thermodynamic properties
        self.assertGreater(self.md.temp, 0.0)

if __name__ == '__main__':
    unittest.main()
