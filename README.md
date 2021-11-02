# QtNEMD - A friendly python front-end for non-equilibrium molecular dynamics
QtNEMD is a python application which provides a 
friendly graphical front-end for molecular dynamics calculations. It uses 
[numpy's f2py utility](https://numpy.org/doc/stable/f2py/) to wrap the
Fortran 77 codes of the Bernhardt group and provides real-time visualisation of
particle positions in two or three dimensions.

This project is currently a work in progress: a roadmap of future features can
be found in `Doc/roadmap.md`.

## Contents
This repository contains the following files:
- `main.py`: contains the functions to set-up, run and visualise the MD
  simulations,
- `Walls_SSGK.f`, `SSGK.inc`: Contains the back-end numerical MD routines.
- `templace_SSGK.in`: Sample input file,
- `driver.f`: set of helper functions for interfacing with the MD routines
  (setting up variables, common-data, etc).

## Compiling
A makefile has been provided for convenience and to document the necessary f2py
commands to compile the python module. Simply do `make` in the root directory
of this repository to compile the `TTCF` module.

First, the makefile generates a `signature file` (`TTCF.pyf`) telling `f2py` 
what subroutines and common variables are contained in the Fortran back-end 
files. Next, it uses the signature file plus the `gfortran` compiler to
generate a python module containing the back-end numerical routines. Python
modules compiled for one system may not work on others (e.g. with a different
operating system), so it's a good idea to do a fresh compile on every machine.

The UI elements and layout are defined in the `.ui` files in `GUI-resources`, 
which must be converted into Qt python modules before they can be used. To do
this, either run `make all` (which is the default target) or `make gui` to only
build the GUI files.

## Usage
To run the main program:
`./main.py`

## Using the MD routines in custom scripts
The Fortran variables and routines are compiled into a single python module
called `TTCF`. All that needs to be done to incorporate them is to add `import
TTCF` to python scripts (or do it in the interpreter).

The `TTCF` module contains three `helper` functions for setting up and running
calculations:
- `main(input_file)`: sets up and runs the MD calculation described in the
  input file,
- `setup(input_file)`: sets up shared variables and initialises particle
  positions but does not run and MD calculations,
- `teardown()`: post-processes simulation trajectories, writes output to files
  and ensures all temporary files are closed.
