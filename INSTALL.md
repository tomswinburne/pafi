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
- Either [download](https://cmake.org/download/) or try e.g.
```bash
module load cmake # linux cluster
conda install cmake # unix
apt install cmake # ubuntu
brew install cmake # osx
... # other options available!

```
## Set local install location
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

2. Configure `LAMMPS` with at least `EXTRA-FIX` package. Example options in `lammps_options.cmake` or see `path/to/lammps/cmake/presets`

3. Build `LAMMPS`:
```bash
# create build folder
mkdir build
cd build

# configure LAMMPS compilation
cmake -C /path/to/lammps_options.cmake ../cmake

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
