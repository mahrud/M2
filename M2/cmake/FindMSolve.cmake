# Try to find the msolve libraries
# See https://msolve.lip6.fr/
#
# This file sets up msolve for CMake. Once done this will define
#  MSOLVE_FOUND		- system has the MSOLVE library with correct version
#  MSOLVE_INCLUDE_DIR	- the MSOLVE include directory
#  MSOLVE_LIBRARIES	- the MSOLVE library
#
# Copyright (c) 2021, Mahrud Sayrafi, <mahrud@umn.edu>
#
# Redistribution and use is allowed according to the terms of the BSD license.

if(NOT MSOLVE_FOUND)
  set(MSOLVE_INCLUDE_DIR NOTFOUND)
  set(MSOLVE_LIBRARIES NOTFOUND)

  # search first if an MSolveConfig.cmake is available in the system,
  # if successful this would set MSOLVE_INCLUDE_DIR and the rest of
  # the script will work as usual
  find_package(MSolve ${MSolve_FIND_VERSION} NO_MODULE QUIET)

  if(NOT MSOLVE_INCLUDE_DIR)
    find_path(MSOLVE_INCLUDE_DIR NAMES msolve/msolve.h neogb/neogb.h
      HINTS ENV MSOLVEDIR
      PATHS ${INCLUDE_INSTALL_DIR} ${CMAKE_INSTALL_PREFIX}/include
      PATH_SUFFIXES msolve
      )
  endif()

  # msolve installs two libraries: libmsolve holds the solver and the rational
  # (multi-modular) Groebner routines, libneogb holds the F4 engine used over
  # finite fields. Both are needed, libmsolve first since it depends on libneogb.
  if(NOT MSOLVE_LIBRARY)
    find_library(MSOLVE_LIBRARY NAMES msolve
      HINTS ENV MSOLVEDIR
      PATHS ${LIB_INSTALL_DIR} ${CMAKE_INSTALL_PREFIX}/lib
      )
  endif()

  if(NOT NEOGB_LIBRARY)
    find_library(NEOGB_LIBRARY NAMES neogb
      HINTS ENV MSOLVEDIR
      PATHS ${LIB_INSTALL_DIR} ${CMAKE_INSTALL_PREFIX}/lib
      )
  endif()

  if(MSOLVE_LIBRARY AND NEOGB_LIBRARY)
    set(MSOLVE_LIBRARIES ${MSOLVE_LIBRARY} ${NEOGB_LIBRARY})
  endif()

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(MSolve DEFAULT_MSG
    MSOLVE_INCLUDE_DIR MSOLVE_LIBRARY NEOGB_LIBRARY)

  mark_as_advanced(MSOLVE_INCLUDE_DIR MSOLVE_LIBRARY NEOGB_LIBRARY MSOLVE_LIBRARIES)
endif()
