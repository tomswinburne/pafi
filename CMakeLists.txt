
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
set(LAMMPS_INCLUDE_DIR ${PREFIX}/include) # location of copied lammps header directory (lammps/*.h)
set(CMAKE_LIBRARY_PATH ${PREFIX}/lib) # location of liblammps_mpi.a
set(LAMMPS_LIBRARY_DIR ${PREFIX}/lib) # location of liblammps_mpi.a


# Just to find lammps for now
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules)



project (pafi)

include_directories(core)
include_directories(rapidxml)

include_directories(${CMAKE_INCLUDE_PATH})
link_directories(${CMAKE_LIBRARY_PATH})

#find_package(LAMMPS REQUIRED)
#include_directories(${LAMMPS_INCLUDE_DIRS})

set(LAMMPS_LIBRARIES "${LAMMPS_LIBRARY_DIR}/liblammps.so")
# Create an imported target (if necessary)
add_library(LAMMPS::Library SHARED IMPORTED)
set_target_properties(LAMMPS::Library PROPERTIES
  IMPORTED_LOCATION "${LAMMPS_LIBRARIES}"
  INTERFACE_INCLUDE_DIRECTORIES "${LAMMPS_INCLUDE_DIR}"
)

#ww # Specify the static library file
#ww set(LAMMPS_LIBRARIES "${LAMMPS_LIBRARY_DIR}/liblammps.a")
#ww 
#ww # Create an imported target
#ww add_library(LAMMPS::Library STATIC IMPORTED)
#ww set_target_properties(LAMMPS::Library PROPERTIES
#ww   IMPORTED_LOCATION "${LAMMPS_LIBRARIES}"
#ww   INTERFACE_INCLUDE_DIRECTORIES "${LAMMPS_INCLUDE_DIR}"
#ww )

#ww # Specify the static library file
#ww set(LAMMPS_LIBRARIES "${LAMMPS_LIBRARY_DIR}/liblammps.a")
#ww 
#ww # Create an imported target
#ww add_library(LAMMPS::Library STATIC IMPORTED)
#ww set_target_properties(LAMMPS::Library PROPERTIES
#ww   IMPORTED_LOCATION "${LAMMPS_LIBRARIES}"
#ww   INTERFACE_INCLUDE_DIRECTORIES "${LAMMPS_INCLUDE_DIR}"
#ww )

message( STATUS "*** !!!!!! PREFIX            : " ${PREFIX} )
message( STATUS "*** !!!!!! CMAKE_INCLUDE_PATH: " ${CMAKE_INCLUDE_PATH})
message( STATUS "*** !!!!!! LAMMPS_INCLUDE_DIR: " ${LAMMPS_INCLUDE_DIR})
message( STATUS "*** !!!!!! CMAKE_LIBRARY_PATH: " ${CMAKE_LIBRARY_PATH})
message( STATUS "*** !!!!!! LAMMPS_LIBRARY_DIR: " ${LAMMPS_LIBRARY_DIR})
message( STATUS "*** !!!!!! LAMMPS_LIBRARIES: " ${LAMMPS_LIBRARIES})


# include VERBOSE for debugging output
add_definitions(-std=c++11 -O3 -DLAMMPS_EXCEPTIONS -DLAMMPS_LIB_MPI) -DLAMMPS_MEMALIGN=64# -DVERBOSE)
#dd_definitions(-std=c++11 -O3 -DLAMMPS_EXCEPTIONS -DLAMMPS_LIB_MPI -DLAMMPS_MEMALIGN=64 -DMPICH_SKIP_MPICXX -DOMPI_SKIP_MPICXX=1 )# -DVERBOSE)

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
	#intel target_link_libraries(pafi ${LAMMPS_LIBRARIES}  imf c iomp5 pthread m dl)
	target_link_libraries(pafi-lammps-preparation ${LAMMPS_LIBRARIES})
	target_link_libraries(pafi-path-test ${LAMMPS_LIBRARIES})
endif(DEBUG)
