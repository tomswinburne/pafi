
cmake_minimum_required(VERSION 3.5)
set(DEBUG 0)

if(DEBUG)
	message("DEBUG=1 MAKING VERBOSE BINARIES!")
	add_definitions(-g -DVERBOSE)
endif(DEBUG)


set(PREFIX "$ENV{PREFIX}") # example value
set(CMAKE_CXX_COMPILER mpicxx) # example value
set(CMAKE_C_COMPILER mpicc) # example value
set(CMAKE_INCLUDE_PATH ${PREFIX}/include) # location of copied lammps header directory (lammps/*.h)
set(CMAKE_LIBRARY_PATH ${PREFIX}/lib) # location of liblammps_mpi.a

# Just to find lammps for now
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules)



project (pafi)

include_directories(core)
include_directories(rapidxml)

include_directories(${CMAKE_INCLUDE_PATH})
link_directories(${CMAKE_LIBRARY_PATH})

find_package(LAMMPS REQUIRED)
include_directories(${LAMMPS_INCLUDE_DIRS})

# include VERBOSE for debugging output
add_definitions(-std=c++11 -O3 -DLAMMPS_EXCEPTIONS -DLAMMPS_LIB_MPI)# -DVERBOSE)
file(GLOB src_files "core/*")

if(DEBUG)
	add_executable(pafi-verbose ${src_files} pafi/pafi.cpp pafi/pafi.hpp)
	add_executable(pafi-lammps-preparation-verbose ${src_files} utils/pafi-lammps-preparation.cpp utils/pafi-lammps-preparation.hpp)
	add_executable(pafi-path-test-verbose ${src_files} utils/pafi-path-test.cpp utils/pafi-path-test.hpp)
	target_link_libraries(pafi-verbose ${LAMMPS_LIBRARIES})
	target_link_libraries(pafi-lammps-preparation-verbose ${LAMMPS_LIBRARIES})
	target_link_libraries(pafi-path-test-verbose ${LAMMPS_LIBRARIES})
else(DEBUG)
	add_executable(pafi ${src_files} pafi/pafi.cpp pafi/pafi.hpp)
	add_executable(pafi-lammps-preparation ${src_files} utils/pafi-lammps-preparation.cpp utils/pafi-lammps-preparation.hpp)
	add_executable(pafi-path-test ${src_files} utils/pafi-path-test.cpp utils/pafi-path-test.hpp)
	target_link_libraries(pafi ${LAMMPS_LIBRARIES})
	target_link_libraries(pafi-lammps-preparation ${LAMMPS_LIBRARIES})
	target_link_libraries(pafi-path-test ${LAMMPS_LIBRARIES})
endif(DEBUG)
