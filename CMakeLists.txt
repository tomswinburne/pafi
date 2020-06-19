
cmake_minimum_required(VERSION 3.5)
set(PREFIX /home/tomswinburne/.local) # example value
set(CMAKE_INCLUDE_PATH ${PREFIX}/include) # location of copied lammps header directory (lammps/*.h)
set(CMAKE_LIBRARY_PATH ${PREFIX}/lib) # location of liblammps_mpi.a
set(CMAKE_CXX_COMPILER /usr/bin/mpicxx) # c++11 compatible compiler
set(CMAKE_CXX_COMPILER /usr/bin/mpicc) # c++11 compatible compiler

project (pafi)

include_directories(core)
include_directories(rapidxml)

include_directories(${CMAKE_INCLUDE_PATH})
link_directories(${CMAKE_LIBRARY_PATH})


set(USE_MILADY 0)



if(USE_MILADY)
  set(MILADY_PATH ${PREFIX}/../Code/milady_lammps/lib/milady)
  include_directories(${MILADY_PATH})
  link_directories(${MILADY_PATH})
endif(USE_MILADY)

# NEWFIX : only enable if compiled with gatherscatter LAMMPS branch
# VERBOSE : debugging output
add_definitions(-std=c++11 -O3 -DLAMMPS_EXCEPTIONS -DNEWFIX) #-DVERBOSE)

file(GLOB src_files "core/*")

add_executable(pafi ${src_files} pafi/pafi.cpp pafi/pafi.hpp)
add_executable(pafi-lammps-path ${src_files} pafi/pafi-lammps-path.cpp pafi/pafi.hpp)


if(USE_MILADY)
  target_link_libraries(pafi lammps_mpi gfortran lapack milady)
else()
  target_link_libraries(pafi lammps_mpi)# gfortran lapack milady)
endif(USE_MILADY)

target_link_libraries(pafi-lammps-path lammps_mpi)
