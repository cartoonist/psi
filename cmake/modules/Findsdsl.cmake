# Find 'sdsl-lite' library.
#
# This set the following variables:
#   - sdsl_FOUND
#   - sdsl_VERSION
#   - sdsl_INCLUDE_DIRS
#   - sdsl_LIBRARIES
#
# and the following imported targets:
#   - sdsl::sdsl

if(sdsl_INCLUDE_DIRS)
  set(sdsl_FIND_QUIETLY TRUE)
else()
  # Try pkg-config, first.
  find_package(PkgConfig QUIET)
  pkg_check_modules(sdsl QUIET sdsl-lite>=2.1.0)
  # If sdsl_INCLUDE_DIRS is not set, this searches for the header/library file.
  find_path(sdsl_INCLUDE_DIRS sdsl/config.hpp)
  find_library(sdsl_LIBRARIES sdsl)
endif(sdsl_INCLUDE_DIRS)

## handle the QUIETLY and REQUIRED arguments and set sdsl_FOUND to TRUE if
## all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(sdsl DEFAULT_MSG sdsl_INCLUDE_DIRS sdsl_LIBRARIES)

mark_as_advanced(sdsl_FOUND sdsl_VERSION sdsl_INCLUDE_DIRS sdsl_LIBRARIES)

# Define `sdsl::sdsl` imported target
if(sdsl_FOUND AND NOT TARGET sdsl::sdsl)
  add_library(sdsl::sdsl INTERFACE IMPORTED)
  set_target_properties(sdsl::sdsl PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${sdsl_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES "${sdsl_LIBRARIES}")
endif()
