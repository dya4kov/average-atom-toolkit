Simulations of materials in extremes require high-quality wide-range equations of state. Average Atom Tools are designed to provide electron properties at high temperatures: energy spectrum, electron density, self-consistent potential, and equation of state. The approach is based on the Finite Temperature Thomas-Fermi atom model with its extensions and corrections.

## Setup instructions

#### Download project and dependencies

1. git clone this project.
2. cd to project directory
3. mkdir **external**
4. cd **external**
5. git clone https://github.com/pybind/pybind11.git

#### Building with CMake:

1. mkdir /path/where/to/build
2. cd /path/where/to/build
3. cmake -DCMAKE_INSTALL_PREFIX=/path/where/to/install /path/to/project
4. make install

#### Don't forget to

set LD_LIBRARY_PATH=/path/where/to/install/lib, PYTHONPATH=/path/where/to/install/python