# - Find lammps
# Find the native LAMMPS headers and libraries.
#
#  LAMMPS_INCLUDE_DIRS - where to find library.h, etc.
#  LAMMPS_LIBRARIES    - List of libraries when using lammps.
#  LAMMPS_FOUND        - True if lammps found.
#

find_path(LAMMPS_INCLUDE_DIR lammps/library.h)

find_library(LAMMPS_LIBRARY NAMES lammps_mpi lammps)

set(LAMMPS_LIBRARIES ${LAMMPS_LIBRARY} )
set(LAMMPS_INCLUDE_DIRS ${LAMMPS_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LAMMPS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(LAMMPS DEFAULT_MSG LAMMPS_LIBRARY LAMMPS_INCLUDE_DIR )

mark_as_advanced(LAMMPS_INCLUDE_DIR LAMMPS_LIBRARY )
