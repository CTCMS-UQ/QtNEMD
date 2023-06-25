# QtNEMD - A friendly python front-end for non-equilibrium molecular dynamics
QtNEMD is a python application which provides a friendly graphical front-end in Qt for molecular 
dynamics calculations using software for developed by the Bernhardt group at the University of
Queensland.

This project is currently a work in progress, bug reports and feature requests welcome.

# Installation
## Prerequisites
You'll need the following software packages installed on your computer:
  1) Python 3
  2) The `pip` python package manager (usually bundled with your python distribution)
  3) Python virtual environment manager `venv`
  4) `make` (the software has been tested on Gnu Make, but other dialects may work)
  5) The `g++` compiler
  6) CMake

## Automatic installation
There is an `install.sh` script in the top-level directory which *should* automatically install all
dependencies (provided you've already installed the pre-requisites). This script as been tested in
Linux (Fedora/CentOS and Ubuntu) and Windows Subsystem for Linux (running Ubuntu). It
does not currently work on MacOS on Apple Silicon (e.g. the M1), as PyQT5 has not (yet) been ported
to this platform. To run it, do:

```
bash ./install.sh
```

This will download all python dependencies and will download and build LAMMPS. LAMMPS is a very
large codebase, so this step may take some time on slow internet connections.

## Manual installation
This software depends on the following python modules:
  1) `numpy`
  2) `PyQt5` graphical user interface library
  3) `pyqtgraph`

These packages can be automatically installed by `pip` via the following command:

```
pip install -r requirements.txt
```

You'll also need to install the [LAMMPS molecular dynamics
software](https://github.com/lammps/lammps) and compile it with [its Python 
interface](https://docs.lammps.org/Python_head.html). Follow the instructions in the [LAMMPS
documentation](https://docs.lammps.org/Build.html)

# Running QtNEMD
Run the code by doing:

```
chmod +x main.py
./main.py
```

Again, this code has been tested on Fedora and Ubuntu Linux and Windows Subsystem for Linux. It
does not currently work on MacOS on Apple Silicon (e.g. the M1).
