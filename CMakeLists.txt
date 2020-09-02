
cmake_minimum_required(VERSION 3.5)
set(PREFIX $ENV{HOME}/.local) # example value - this should match value in lammps.cmake
set(CMAKE_INCLUDE_PATH ${PREFIX}/include) # location of copied lammps header directory (lammps/*.h)
set(CMAKE_LIBRARY_PATH ${PREFIX}/lib) # location of liblammps_mpi.a
set(CMAKE_CXX_COMPILER /usr/bin/env mpic++) # c++11 compatible compiler

project (pafi)

include_directories(core)
include_directories(rapidxml)

include_directories(${CMAKE_INCLUDE_PATH})
link_directories(${CMAKE_LIBRARY_PATH})


# VERBOSE : debugging output
add_definitions(-std=c++11 -O3 -DVERBOSE -DLAMMPS_EXCEPTIONS)
file(GLOB src_files "core/*")

add_executable(pafi ${src_files} pafi/pafi.cpp pafi/pafi.hpp)
add_executable(pafi-lammps-path ${src_files} pafi/pafi-lammps-path.cpp pafi/pafi.hpp)


target_link_libraries(pafi lammps)

target_link_libraries(pafi-lammps-path lammps)
