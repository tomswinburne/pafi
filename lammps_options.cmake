# FOR LAMMPS CMAKE BUILD

# set installation location
set(CMAKE_INSTALL_PREFIX "$ENV{PREFIX}" CACHE BOOL "" FORCE)

# enforce c++11 standards
#set(CCFLAGS -g -O3 -std=c++11)
set(CMAKE_CXX_STANDARD 11) # or a newer standard like 14, 17, 20
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CCFLAGS -std=c++11 -O3  -DLAMMPS_EXCEPTIONS -DLAMMPS_LIB_MPI -DLAMMPS_MEMALIGN=64 )

# compile a binary and shared library
set(BUILD_SHARED_LIBS ON CACHE BOOL "" FORCE)

# allow error messages (very useful)
set(LAMMPS_EXCEPTIONS ON CACHE BOOL "" FORCE)

# minimal packages to run example (MANYBODY, EXTRA-FIX) and generate new pathway (REPLICA for "fix neb")
set(ALL_PACKAGES MANYBODY EXTRA-FIX REPLICA)

foreach(PKG ${ALL_PACKAGES})
  set(PKG_${PKG} ON CACHE BOOL "" FORCE)
endforeach()
