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
*For LAMMPS version older than 28 July 2021, or to statically link with traditional make, please follow [these instructions](STATIC_MAKE.md)*

- Install `cmake` by [download](https://cmake.org/download/) or try e.g.
`[conda/apt/brew] install cmake` or `module load cmake`

- Set local install location
```bash
export PREFIX=${HOME}/.local # example value
```
## Build LAMMPS
- [Download](https://lammps.sandia.gov/download.html) or `git clone https://github.com/lammps/lammps.git`

- Configure `LAMMPS` with `EXTRA-FIX` package. See `lammps_options.cmake` in `PAFI` or `path/to/lammps/cmake/presets`, then build:
```bash
cd /path/to/lammps
mkdir build
cd build
cmake -C /path/to/pafi/lammps_options.cmake ../cmake
cmake --build .
cmake --install .
```

## Build PAFI

- Specify compiler in `CMakeLists.txt`:
```cmake
   set(CMAKE_CXX_COMPILER path/to/mpic++)
```

- Make pafi build folder, run `cmake`, ensuring `PREFIX` is in your environment
```bash
cd path/to/pafi
export PREFIX=${HOME}/.local # if in different shell to LAMMPS compilation
mkdir build
cd build
cmake ..
make # or try make -j4 for parallel make using 4 cores
```
