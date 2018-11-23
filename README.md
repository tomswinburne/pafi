         _______      _______      _______     _________
        (  ____ )    (  ___  )    (  ____ \    \__   __/
        | (    )|    | (   ) |    | (    \/       ) (
        | (____)|    | (___) |    | (__           | |
        |  _____)    |  ___  |    |  __)          | |
        | (          | (   ) |    | (             | |
        | )          | )   ( |    | )          ___) (___
        |/           |/     \|    |/           \_______/
        Projected    Average      Force        Integrator


v0.4 :copyright: TD Swinburne and M-C Marinica 2018 MIT License



Beta version of code used in [this paper](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.120.135503)
> *Unsupervised Calculation of Free Energy Barriers in Large Crystalline Systems*   
> T.D. Swinburne and M.-C. Marinica, Physical Review Letters 120 (13), 135503, 2018



## Patching and compiling LAMMPS

1. [Download LAMMPS](http://lammps.sandia.gov/download.html) and read the [installation instructions](http://lammps.sandia.gov/doc/Section_start.html)

2. Before installing any packages, in the LAMMPS root directory (parent directory of src folder) run the pafi patch, adding the USER-PAFI package and support for PAFI atom styles
```
cd /path/to/lammps_repo/
patch -p0 < /path/to/user-pafi_lammps_Aug18.patch # modifies /path/to/lammps_repo/src folder
```

3. Install USER-PAFI and any packages you desire (replica for NEB)
```
cd src
make yes-user-pafi
make yes-replica # for NEB calculation
make yes-package_name
```

3. Compile as a static library (and optionally binary initial NEB calculation)
```
   make mpi mode=lib # liblammps_mpi.a library for pafi
   make mpi # lmp_mpi binary for running initial NEB calculation if desired
```
4. Copy library to your local lib/ and headers to local include/, at e.g. ${HOME}/.local
```
   cp liblammps_mpi.a ${HOME}/.local/lib
   mkdir ${HOME}/.local/include/lammps
   cp *.h ${HOME}/.local/include/lammps/
```

## Compiling Boost C++ libraries
1. Download tarball from https://www.boost.org/users/download/
2. Bootstrap with the local install prefix
```
   ./bootstrap.sh --prefix=${HOME}/.local
```
3. Build with gcc and c++11
```
   ./b2 toolset=gcc cxxflags=-std=c++11
```
4. Install
```
   ./b2 install
```

## Compiling PAFI with cmake (recommended)
1. Download and install cmake from https://cmake.org/download/
2. Specify environment variables in CMakeLists.txt:
```
   set(CMAKE_INCLUDE_PATH ${HOME}/.local/include)
   set(CMAKE_LIBRARY_PATH ${HOME}/.local/lib)
   set(CMAKE_CXX_COMPILER ${HOME}/.local/bin/mpicxx)
   set(CMAKE_C_COMPILER ${HOME}/.local/bin/mpicc)
```
3. Make pafi build folder and run cmake
```
   mkdir build
   cd build
   cmake ..
```
## Compiling PAFI with make (doesn't run configure)
1. Specify environment variables with PREFIX
```
   PREFIX=${HOME}/.local/
```
2. Make pafi build folder and run make
```
   mkdir build
   cd build
   make -f../makefile
```

## Calculation of free energy barrier between states

1. First set up a LAMMPS neb calculation as described [here](http://lammps.sandia.gov/doc/neb.html)

2. In the LAMMPS script, use the "write_data" command to dump all the NEB knots, i.e.
```
variable u uloop N_NEB_IMAGES

neb etol ftol N1 N2 Nevery file-style arg keyword

write_data neb_knot_file.$u
```
3. Move the knots to folder and include and configure the configuratio xml file, as shown in the examples

4. Move run.sh to run directory and change path for pafi binary soft link

5. Run PAFI as e.g.
```
mpirun -np NPROCS ./pafi
```
## Output

TBA
