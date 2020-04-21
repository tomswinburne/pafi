set(PREFIX /home/tomswinburne/.local) # example value
set(CMAKE_INCLUDE_PATH ${PREFIX}/newinclude) # location of copied lammps header directory (lammps/*.h)
set(CMAKE_LIBRARY_PATH ${PREFIX}/lib) # location of liblammps_mpi.a
set(CMAKE_CXX_COMPILER /usr/bin/mpicxx) # c++11 compatible compiler

cmake_minimum_required(VERSION 2.8.9)
project (pafi)



include_directories(core)
include_directories(rapidxml)

include_directories(${CMAKE_INCLUDE_PATH})
link_directories(${CMAKE_LIBRARY_PATH})

add_definitions(-std=c++11 -O3 -DLAMMPS_EXCEPTIONS -DNEWFIX)#-DVERBOSE)

file(GLOB src_files "core/*")

add_executable(pafi ${src_files} pafi/pafi.cpp pafi/pafi.hpp)
add_executable(pafi-lammps-path ${src_files} pafi/pafi-lammps-path.cpp pafi/pafi.hpp)

target_link_libraries(pafi lammps_mpigs)
target_link_libraries(pafi-lammps-path lammps_mpigs)
