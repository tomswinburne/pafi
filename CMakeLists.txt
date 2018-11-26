set(CMAKE_INCLUDE_PATH /home/tomswinburne/.local/include)
set(CMAKE_LIBRARY_PATH /home/tomswinburne/.local/lib)
set(CMAKE_CXX_COMPILER /home/tomswinburne/miniconda2/bin/mpicxx)
set(CMAKE_C_COMPILER /home/tomswinburne/miniconda2/bin/mpicc)

cmake_minimum_required(VERSION 2.8.9)
project (pafi)

include_directories(core)
include_directories(${CMAKE_INCLUDE_PATH})
link_directories(${CMAKE_LIBRARY_PATH})

add_definitions(-std=c++11 -O3 -g -DBOOST_ALL_DYN_LINK)

file(GLOB src_files "core/*")

set(pafi_src_files ${src_files} pafi/pafi.cpp pafi/pafi.hpp)
set(prepare_src_files ${src_files} prepare/prepare.cpp prepare/prepare.hpp)
add_executable(pafi ${pafi_src_files})
target_link_libraries ( pafi pthread lammps_mpi boost_random boost_system boost_filesystem)

# For later implemenation-
# boost_timer boost_serialization
# boost_wserialization boost_log boost_log_setup
