         _______      _______      _______     _________
        (  ____ )    (  ___  )    (  ____ \    \__   __/
        | (    )|    | (   ) |    | (    \/       ) (
        | (____)|    | (___) |    | (__           | |
        |  _____)    |  ___  |    |  __)          | |
        | (          | (   ) |    | (             | |
        | )          | )   ( |    | )          ___) (___
        |/           |/     \|    |/           \_______/
        Projected    Average      Force        Integrator

## [Main Page](README.md)

## [Getting Started](TUTORIAL.md)

# Installation

## Install `cmake`
- On a cluster, try `module load cmake`
- On Linux, try e.g. `[apt/yum] install cmake`
- On OSX, try e.g.  `[conda/brew] install cmake`
- Alternatively [download](https://cmake.org/download/) and install `cmake` manually

## Set local install environment
- This is where `LAMMPS` will be installed, and where `PAFI` will look
```bash
export PREFIX=${HOME}/.local # example value
```

## Build `LAMMPS`

*For LAMMPS version older than 28 July 2021, or to statically link with traditional make, please follow [these instructions](STATIC_MAKE.md)*

1. [Download](https://lammps.sandia.gov/download.html) a tarball or clone into `LAMMPS` source code
```bash
git clone https://github.com/lammps/lammps.git
```

2. Go to root of `LAMMPS` distribution
```bash
# go to root of distribution
cd /path/to/lammps
```

2. Create a file `my_options.cmake`:
```cmake

# set installation location
set(CMAKE_INSTALL_PREFIX "$ENV{PREFIX}")

# enforce c++11 standards
set(CCFLAGS -g -O3 -std=c++11)

# compile a binary and shared library
set(BUILD_SHARED_LIBS ON CACHE BOOL "" FORCE)

# allow error messages (very useful)
set(LAMMPS_EXCEPTIONS ON CACHE BOOL "" FORCE)

# minimal packages to run example (MANYBODY, EXTRA-FIX) and generate new pathway (REPLICA for "fix neb")
set(ALL_PACKAGES MANYBODY EXTRA-FIX REPLICA)

foreach(PKG ${ALL_PACKAGES})
  set(PKG_${PKG} ON CACHE BOOL "" FORCE)
endforeach()
```

Build `LAMMPS`:
```bash
# create build folder
mkdir build
cd build

# configure LAMMPS compilation
cmake -C ../my_options.cmake ../cmake

# compile LAMMPS
cmake --build .

# install LAMMPS into $PREFIX
cmake --install .
```

## Compile `PAFI`

1. Specify compiler in CMakeLists.txt:
```cmake
   set(CMAKE_CXX_COMPILER path/to/mpic++)
```

2. Make pafi build folder, run `cmake`, ensuring `PREFIX` is in your environment
```bash
   export PREFIX=${HOME}/.local # if in different shell to LAMMPS compilation
   mkdir build
   cd build
   cmake ..
   make # or try make -j4 for parallel make using 4 cores
```
