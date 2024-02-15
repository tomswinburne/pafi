# FOR LAMMPS CMAKE BUILD

# set installation location
set(CMAKE_INSTALL_PREFIX "$ENV{PREFIX}" CACHE BOOL "" FORCE)
set(CMAKE_BUILD_PARALLEL_LEVEL 4) 
#set(Python_FIND_FRAMEWORK LAST)

#set(PYTHON_EXECUTABLE /path/to/python CACHE BOOL "" FORCE) 

# enforce c++11 standards
set(CCFLAGS -g -O3 -std=c++11)

# compile a binary and shared library
set(BUILD_SHARED_LIBS ON CACHE BOOL "" FORCE)

# allow error messages (very useful)
set(LAMMPS_EXCEPTIONS ON CACHE BOOL "" FORCE)

# minimal packages to run example (MANYBODY, EXTRA-FIX) and generate new pathway (REPLICA for "fix neb")
set(ALL_PACKAGES MANYBODY EXTRA-FIX REPLICA ML-SNAP MOLECULE PHONON)

foreach(PKG ${ALL_PACKAGES})
  set(PKG_${PKG} ON CACHE BOOL "" FORCE)
endforeach()
