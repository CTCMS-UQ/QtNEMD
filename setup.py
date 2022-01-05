from setuptools import setup
from setuptools import find_packages

setup(
    name='QtNEMD', #project name
    version='0.0.4-alpha',
    description='Graphical frontend for non-equilibrium molecular-dynamics',
    #url
    author='Emily Kahl',
    author_email='e.kahl@uq.edu.au', 
    license='GPL-3.0',
    packages=find_packages(),
    install_requires=['numpy','PyQt5','pyqtgraph']
)
