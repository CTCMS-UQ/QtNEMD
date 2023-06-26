#!/bin/bash

# Install python dependencies
python3 -m pip install --upgrade pip
pip3 install -r requirements.txt

# Now download and build LAMMPS from the CTCMS branch
git clone https://github.com/CTCMS-UQ/lammps.git
cd lammps
git checkout sprint-develop

mkdir -p build
# Basic CPU
cmake -B build ./cmake -C ./cmake/presets/basic.cmake -C ./cmake/presets/nolib.cmake \
-D CMAKE_CXX_COMPILER=g++ \
-D BUILD_SHARED_LIBS=yes \
-D PKG_MOLECULE=yes \
-D PKG_MOL-SLLOD=yes \
-D BUILD_MPI=no -D CMAKE_BUILD_TYPE=RelWithDebInfo

cd build
make -j $(nproc)
make install-python
