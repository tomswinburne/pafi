# preset that turns on just a few, frequently used packages
# Example LAMMPS configuration using cmake build

# PREFIX should match the value in PAFI's CMakeLists.txt. (Could be set by environment variable, i.e. $ENV{PREFIX})
set(PREFIX "$ENV{HOME}/.local")

# Desired packages - here the minimial to run NEB and PAFI with standard EAM potentials
set(ALL_PACKAGES MANYBODY USER-MISC REPLICA)


# Following should be left alone unless you know what you're doing....
set(LAMMPS_EXCEPTIONS ON CACHE BOOL "" FORCE)
set(BUILD_LIB ON CACHE BOOL "" FORCE)
set(BUILD_SHARED_LIBS ON CACHE BOOL "" FORCE)
set(CXX_LINK_FLAGS "-std=c++11" ON CACHE STRING "" FORCE)
set(CMAKE_INSTALL_PREFIX ${PREFIX} CACHE PATH "Custom Install Path" FORCE )

set(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT OFF CACHE BOOL "" FORCE)

foreach(PKG ${ALL_PACKAGES})
  set(PKG_${PKG} ON CACHE BOOL "" FORCE)
endforeach()
