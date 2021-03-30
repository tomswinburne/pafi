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

## Compile `LAMMPS` with `USER-MISC` package
1. `PAFI` is now integrated into `LAMMPS` as part of the `USER-MISC` package.
You can [download](https://lammps.sandia.gov/download.html) a tarball from the `LAMMPS`
website or clone the public repository with
```bash
git clone https://github.com/lammps/lammps.git
```

2. Install `USER-MISC` and any packages you desire (e.g. replica for `NEB`)
```bash
cd /path/to/lammps/src
make yes-user-misc
make yes-replica # for NEB calculation
make yes-package # (e.g. manybody for EAM potentials etc)
```

3. In the appropriate Makefile add `-std=c++11` to `CCFLAGS` and `LINKFLAGS` and
add `-DLAMMPS_EXCEPTIONS` to `LMP_INC` to allow `PAFI` to read `LAMMPS` error messages.
This is very useful when running your own simulations. For `src/MAKE/Makefile.mpi` this reads
 ```make
CCFLAGS =	-g -O3 -std=c++11
LINKFLAGS =	-g -O3 -std=c++11
LMP_INC =	-DLAMMPS_GZIP -DLAMMPS_MEMALIGN=64  -DLAMMPS_EXCEPTIONS
```

4. Compile static library and binary Consult [LAMMPS documentation](http://lammps.sandia.gov/doc/Section_start.html) for details
```bash
   make mpi mode=lib # liblammps_mpi.a library for pafi
   make mpi # lmp_mpi binary for running initial NEB calculation if desired
```

4. Copy library to your local lib/ and headers to local include/, at e.g. ${HOME}/.local
```bash
  export PREFIX=${HOME}/.local # example value
  cp liblammps_mpi.a ${PREFIX}/lib/liblammps_mpi.a
  mkdir ${PREFIX}/include/lammps
  cp *.h ${PREFIX}/include/lammps/
```


## Compile `PAFI`
0. `PAFI` requires `cmake` to compile:
- On a cluster, try `module load cmake`
- On Linux, try `[apt/yum] install cmake`
- Alternatively [download](https://cmake.org/download/) and install `cmake` manually

*Technical point: `LAMMPS` can also be built with `cmake` . However, this is causes
[complications](https://lammps.sandia.gov/doc/Build_link.html) with static linking.*

1. Specify compiler in CMakeLists.txt:
```make
   set(CMAKE_CXX_COMPILER path/to/mpic++)
```

2. Make pafi build folder, run cmake and make
```bash
   export PREFIX=${HOME}/.local # if in different shell to LAMMPS compilation
   mkdir build
   cd build
   cmake ..
   make # or try make -j4 for parallel make using 4 cores
```
