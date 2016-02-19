#.rst:
# FindLibIomp
# --------
#
# Find libiomp Intel OpenMP runtime library
#
# The module defines the following variables:
#
# ::
#
#   LIBIOMP_FOUND - true if system has libiomp
#   LIBIOMP_INCLUDE_DIR - libiomp include directory
#   LIBIOMP_LIBRARIES - libraries needed to use libiomp
#   LIBIOMP_DEFINITIONS - compiler switches required for using libiomp
#   LIBIOMP_VERSION_STRING - version of libiomp found
#
set(CMAKE_REQUIRED_QUIET ${LIBIOMP_FIND_QUIETLY})

if(NOT CMAKE_REQUIRED_QUIET)
  message(STATUS "Looking for Intel OpenMP runtime library (libiomp)")
endif()

find_path(LIBIOMP_INCLUDE_DIR NAMES omp.h PATH_SUFFIXES include libiomp include/libiomp)
find_library(LIBIOMP_LIBRARIES NAMES iomp5 libiomp5 PATH_SUFFIXES lib lib64)
get_filename_component(LIBIOMP_LIBRARIES_DIR ${LIBIOMP_LIBRARIES} DIRECTORY)

if(LIBIOMP_INCLUDE_DIR AND EXISTS "${LIBIOMP_INCLUDE_DIR}/omp.h")
  # The matching code within this block works, but could probably be simplified
  file(STRINGS "${LIBIOMP_INCLUDE_DIR}/omp.h" kmp_major_version_str
    REGEX "#[\t ]+define[\t ]+KMP_VERSION_MAJOR+[ \t]+[0-9]+")
  file(STRINGS "${LIBIOMP_INCLUDE_DIR}/omp.h" kmp_minor_version_str
    REGEX "#[\t ]+define[\t ]+KMP_VERSION_MINOR+[ \t]+[0-9]+")
    file(STRINGS "${LIBIOMP_INCLUDE_DIR}/omp.h" kmp_build_version_str
    REGEX "#[\t ]+define[\t ]+KMP_VERSION_BUILD+[ \t]+[0-9]+")
  
  string(REGEX REPLACE "#[\t ]+define[\t ]+KMP_VERSION_MAJOR[\t ]+([0-9]+)" "\\1" KMP_VERSION_MAJOR "${kmp_major_version_str}")
  string(REGEX REPLACE "#[\t ]+define[\t ]+KMP_VERSION_MINOR[\t ]+([0-9]+)" "\\1" KMP_VERSION_MINOR "${kmp_minor_version_str}")
  string(REGEX REPLACE "#[\t ]+define[\t ]+KMP_VERSION_BUILD[\t ]+([0-9]+)" "\\1" KMP_VERSION_BUILD "${kmp_build_version_str}")

  set(LIBIOMP_VERSION_STRING "${KMP_VERSION_MAJOR}.${KMP_VERSION_MINOR}.${KMP_VERSION_BUILD}")
  unset(kmp_version_str)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LibIomp
  REQUIRED_VARS LIBIOMP_INCLUDE_DIR LIBIOMP_LIBRARIES LIBIOMP_LIBRARIES_DIR
  VERSION_VAR LIBIOMP_VERSION_STRING)

if(LIBIOMP_FOUND)
  message(STATUS "Looking for Intel OpenMP runtime library (libiomp) -- found")
  message(STATUS "Intel OpenMP runtime library (libiomp) version: ${LIBIOMP_VERSION_STRING}")
else()
  message(STATUS "Looking for Intel OpenMP runtime library (libiomp) -- not found")
endif()
  
mark_as_advanced(LIBIOMP_INCLUDE_DIR LIBIOMP_LIBRARIES)
# FindLibIomp.cmake ends here
