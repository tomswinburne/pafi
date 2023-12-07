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

To run `PAFI`, you need the `LAMMPS` library. `PAFI` can be built using either the static `*.a` or shared `*.so` libraries of `LAMMPS`. Therefore, you should first build `LAMMPS` and its libraries, and then compile `PAFI`. Here, we detail only the procedure based on `cmake`. 

*For LAMMPS version older than 28 July 2021, or to link with traditional make, please follow [these instructions](STATIC_MAKE.md)*

For `PAFI` installation using `cmake`  please follow the following steps. 


**Note:**  Be sure that `cmake` is intalled on your system. Install `cmake` by [download](https://cmake.org/download/) or try e.g.
`[conda/apt/brew] install cmake` or `module load cmake`. Attention:  the version should be newer than `3.20.0`.

1. **Set the local installation directory** of `LAMMPS` using the `$PREFIX` environment variable. This floder will also be where PAFI accesses the necessary library includes. 
   In bash a variable can be easily set in the shell:

   ```bash
   export PREFIX=${HOME}/.local # example value
   ```
2. **Build LAMMPS**

  * 2.1 Download [here](https://lammps.sandia.gov/download.html) or clone by `git clone https://github.com/lammps/lammps.git`

  * 2.2 Configure `LAMMPS` with the `EXTRA-FIX` package and any others you like for the desired force field. We usually use `lammps_options.cmake` in `PAFI` to set up the 
    necessary packages.  Otherwise, for the same result, you can use the `-D` option in `cmake` as shown in the example below.  

    ```bash
    cd /path/to/lammps
    mkdir build
    cd build
    cmake  -DLAMMPS_EXCEPTIONS=on -DLAMMPS_LIB_MPI=on -DLAMMPS_MEMALIGN=64   \ 
           -D PKG_REPLICA=ON -D PKG_MANYBODY=0N -D PKG_EXTRA-FIX=ON  \
           -C /path/to/pafi/lammps_options.cmake -C ../cmake/preset/my.cmake   ../cmake
    ```
    In the above example, the packages `EXTRA-FIX`, `MANYBODY`, or `REPLICA` are set to `ON` directly by CMake. You can add any package you want in this way. 
    If you inspect the contents of `/path/to/pafi/lammps_options.cmake`, you will find other methods to set desired packages to `ON`.

  * 2.3 In the above example, the extra files include:

    *   1.3a The file `/path/to/pafi/lammps_options.cmake` can be used to specify additional LAMMPS compilation options. For instance, the `LAMMPS_EXCEPTIONS` option can be set in this file, as well as `C++` flags, etc.
    *   1.3b The file `../cmake/presets/my.cmake` can be used for custom compiler choices. In `../cmake/presets`, there are many example files provided by LAMMPS developers. As an example, here is my choice to use GNU `g++` and `gcc` for `C++` and `C` compilers, and `ifort` for the Fortran compiler.

    ```bash
    # preset that will explicitly request gcc/g++ compilers with support for MPI and OpenMP
    set(CMAKE_CXX_COMPILER "g++" CACHE STRING "" FORCE)
    set(CMAKE_C_COMPILER "gcc" CACHE STRING "" FORCE)
    set(CMAKE_Fortran_COMPILER "ifort" CACHE STRING "" FORCE)
    set(CMAKE_Fortran_FLAGS_DEBUG "-Wall -Wextra -g" CACHE STRING "" FORCE)
    set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-Wall -Wextra -g -O2 -DNDEBUG" CACHE STRING "" FORCE)
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "" FORCE)
    set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Og -g" CACHE STRING "" FORCE)
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -O2 -DNDEBUG" CACHE STRING "" FORCE)
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "" FORCE)
    set(MPI_CXX "g++" CACHE STRING "" FORCE)
    set(MPI_CXX_COMPILER "mpicxx" CACHE STRING "" FORCE)
    set(MPI_C "gcc" CACHE STRING "" FORCE)
    set(MPI_C_COMPILER "mpicc" CACHE STRING "" FORCE)
    set(CMAKE_C_FLAGS_DEBUG "-Wall -Og -g" CACHE STRING "" FORCE)
    set(CMAKE_C_FLAGS_RELWITHDEBINFO "-g -O2 -DNDEBUG" CACHE STRING "" FORCE)
    set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "" FORCE)
    set(MPI_Fortran "ifort" CACHE STRING "" FORCE)
    set(MPI_Fortran_COMPILER "mpif90" CACHE STRING "" FORCE)
    unset(HAVE_OMP_H_INCLUDE CACHE)
    ```
      
  * 2.4 Build `LAMMPS` on 5 processors (5 is an example :) ), create the shared `*lammps*.so` library, and install it. : 
    ```bash
    cmake --build . -j5
    cmake --install .
    ```
    If you have more than 5 processors, you can specify more by using `-j no_of_procs`, which will make the compilation faster.
    At the end the libraries used then by `PAFI` will be installed in `${PREFIX}/lib` and the includes in `${PREFIX}/include`.
    Now the directory `${PREFIX}/lib` sould contains `*.so` `LAMMPS` library. Here is mine: 
    ```bash
    [pafi@irene195 lammps.git]$ ll ${PREFIX}/lib/
    -rw-r----- 1 pafi users       994 Dec  7 10:12 liblammps.pc
    -rwxr-x--- 1 pafi users  88719736 Dec  7 10:12 liblammps.so
    -rwxr-x--- 1 pafi users  88719736 Dec  7 10:12 liblammps.so.0
    ```
    

3. **Build PAFI**
   
   `PAFI` compilation uses its own `cmake`. For example you can specify the compiler in `path/to/pafi/CMakeLists.txt`:
    ```cmake
       set(CMAKE_CXX_COMPILER path/to/mpic++)
    ```

   Make pafi build folder, run `cmake`, ensuring `PREFIX` is in your environment
   ```bash
   cd path/to/pafi
   mkdir build
   cd build
   cmake ..
   make -j5 
   ```
4.  **Test your instalation**
   This is simple. Go into `/path/to/pafi/example/pathway_test/` and then `./run.sh > out`. If the output of the 
   `tail -10 out` is something like that:
    ```bash
      0.8        -3.77761e-05        -7.48921e-05            0.855823                   0         1.02009e-15             1.00008         1.92701e-05                   1
      # 1/1 2/2 8
      0.9        -7.91794e-05        -0.000158042            0.577326                   0         1.17199e-15             1.00035         1.99967e-05                   1
      # 1/1 2/2 8
        1        -1.74119e-10         -2.7142e-10         1.26246e-09                   0         6.10929e-17                   1         1.72426e-13                   1
    
     T=0K run complete, integrating FEP....
    
     Integration complete
    ``` 
   You have done it! Now `PAFI` is ready for many HPC projects !!!! 