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

Please cite the above when publishing results using PAFI

## Patch and compile LAMMPS

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
## Compile Eigen C++ library
0. Eigen is only used in Boundary.cpp to invert one 3x3 matrix for periodic boundary conditions with triclinic matricies. Straightforward to replace with your own code! (Trivial when cell is orthorhombic)
1. Download and extract tarball from http://eigen.tuxfamily.org/index.php
2. In root of distribution, make build directory and compile
```
mkdir build
cd build
FC="nofortran" cmake -DCMAKE_INSTALL_PREFIX=${HOME}/.local ../
```
3. Install
```
make install
```
## Compile Boost C++ library
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

## Compile PAFI
1. Download and install cmake from https://cmake.org/download/
2. Specify environment variables in CMakeLists.txt:
```
   set(HOME your/home/path)
   set(CMAKE_INCLUDE_PATH ${HOME}/.local/include)
   set(CMAKE_LIBRARY_PATH ${HOME}/.local/lib)
   set(CMAKE_CXX_COMPILER path/to/mpic++)
```
3. Make pafi build folder, run cmake and make
```
   mkdir build
   cd build
   cmake ..
   make
```

## Calculation of free energy barrier between states

0. Tarball in example folder has premade NEB calculation (SIA in EAM-Fe) for testing

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

1. PAFI will try to create a folder with a name given by the DumpFolder parameter in config.xml. If it already exists, PAFI will create the first available directory named DumpFolder_i, where i is an integer less than 20.

2. In DumpFolder PAFI will create subfolders for each temperature in the run. In each subfolder there will be a file `dev_r_T.dat` with the ensemble average and variance pathway deviation from each hyperplane and a file `free_energy_profile_T` that has the integrated FEP.

## TODO
1. Restart files from pathway deviations
2. Smoothed spline interpolation for more general reference pathways
3. tbc....
