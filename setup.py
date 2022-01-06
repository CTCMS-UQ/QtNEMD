import sys
import subprocess
from setuptools import setup
from setuptools import find_packages
from setuptools.command.build_ext import build_ext

class Build(build_ext):
     """Build the Fortran backend. """
     def run(self):
         build_command = ["make"]
         if subprocess.call(build_command) != 0:
             sys.exit(-1)
         build_ext.run(self)

setup(
    name='QtNEMD', #project name
    version='0.0.4-alpha',
    description='Graphical frontend for non-equilibrium molecular-dynamics',
    #url
    author='Emily Kahl',
    author_email='e.kahl@uq.edu.au', 
    license='GPL-3.0',
    packages=find_packages(),
    install_requires=['numpy','PyQt5','pyqtgraph'],
    cmdclass={
         'build_ext': Build,
     }
)
