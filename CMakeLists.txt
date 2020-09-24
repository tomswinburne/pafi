
cmake_minimum_required(VERSION 3.5)
set(PREFIX "$ENV{PREFIX}") # example value
set(CMAKE_CXX_COMPILER /usr/bin/env mpic++) # example value
set(CMAKE_INCLUDE_PATH ${PREFIX}/include) # location of copied lammps header directory (lammps/*.h)
set(CMAKE_LIBRARY_PATH ${PREFIX}/lib) # location of liblammps_mpi.a

project (pafi)

include_directories(core)
include_directories(rapidxml)

include_directories(${CMAKE_INCLUDE_PATH})
link_directories(${CMAKE_LIBRARY_PATH})

add_definitions(-std=c++11 -O3 -DLAMMPS_EXCEPTIONS) # -DVERBOSE # include for debugging output

file(GLOB src_files "core/*")

add_executable(pafi ${src_files} pafi/pafi.cpp pafi/pafi.hpp)
add_executable(pafi-lammps-path ${src_files} pafi/pafi-lammps-path.cpp pafi/pafi.hpp)


target_link_libraries(pafi lammps)

target_link_libraries(pafi-lammps-path lammps)
