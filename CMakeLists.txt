set(CMAKE_INCLUDE_PATH /home/tomswinburne/.local/include) # location of copied lammps header directory (lammps/*.h)
set(CMAKE_LIBRARY_PATH /home/tomswinburne/.local/lib) # location of liblammps_mpi.a
set(CMAKE_CXX_COMPILER /usr/bin/mpicxx) # c++11 compatible compiler

cmake_minimum_required(VERSION 2.8.9)
project (pafi)

include_directories(core)
include_directories(rapidxml)

include_directories(${CMAKE_INCLUDE_PATH})
link_directories(${CMAKE_LIBRARY_PATH})

add_definitions(-std=c++11 -O3 -g)

file(GLOB src_files "core/*")

add_executable(pafi ${src_files} pafi/pafi.cpp pafi/pafi.hpp)

target_link_libraries(pafi lammps_mpi)
