
cmake_minimum_required(VERSION 3.5)
set(DEBUG 0)

set(PREFIX "$ENV{PREFIX}") # example value
set(CMAKE_CXX_COMPILER /usr/local/bin/mpicxx) # example value
set(CMAKE_C_COMPILER /usr/local/bin/mpicc) # example value
set(CMAKE_INCLUDE_PATH ${PREFIX}/include) # location of copied lammps header directory (lammps/*.h)
set(CMAKE_LIBRARY_PATH ${PREFIX}/lib) # location of liblammps_mpi.a

project (pafi)

include_directories(core)
include_directories(rapidxml)

include_directories(${CMAKE_INCLUDE_PATH})
link_directories(${CMAKE_LIBRARY_PATH})

# include VERBOSE for debugging output
add_definitions(-std=c++11 -O3 -DLAMMPS_EXCEPTIONS -DLAMMPS_LIB_MPI)# -DVERBOSE)
if(DEBUG)
	add_definitions(-g -DVERBOSE)
endif(DEBUG)
file(GLOB src_files "core/*")

file(GLOB ti_files "pafi-ti/*")

add_executable(pafi ${src_files} pafi/pafi.cpp pafi/pafi.hpp)
add_executable(pafi-ti ${src_files} pafi/pafi.cpp ${ti_files})

add_executable(pafi-lammps-preparation ${src_files} pafi/pafi-lammps-preparation.cpp pafi/pafi.hpp)
add_executable(pafi-path-test ${src_files} pafi/pafi-path-test.cpp pafi/pafi.hpp)


target_link_libraries(pafi lammps_mpi)
target_link_libraries(pafi-ti lammps_mpi)
target_link_libraries(pafi-lammps-preparation lammps_mpi)
target_link_libraries(pafi-path-test lammps_mpi)
