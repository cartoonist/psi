cmake_minimum_required(VERSION 3.10)
project(psi VERSION 0.0.2 LANGUAGES CXX)

# Options
option(BUILD_PSIKT "Build PSI seeding/indexing command-line tool (`psikt`)" ON)
option(BUILD_AUX_TOOLS "Build auxiliary tools" OFF)
option(BUILD_TESTING "Build test programs" OFF)
option(BUILD_DOC "Build documentation" OFF)

# Adding CMake module path.
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)
# Finding dependencies.
find_package(ZLIB REQUIRED)   # required by SeqAn
find_package(BZip2 REQUIRED)  # required by SeqAn
find_package(OpenMP)          # required by SeqAn
find_package(SeqAn QUIET)
find_package(sdsl REQUIRED)
find_package(gum QUIET)
# Include external dependencies
include(${PROJECT_SOURCE_DIR}/ext/external.cmake)

# Force CMake to use `-I` instead of `-isystem`
set(CMAKE_NO_SYSTEM_FROM_IMPORTED TRUE)

# Load `BuildDoxygen` module
if(BUILD_DOC)
  include(BuildDoxygen)
endif(BUILD_DOC)
# Load `CheckGitRevision` module
include(CheckGitRevision)

# Defining the configure header file
configure_file(include/psi/config.hpp.in ${PROJECT_SOURCE_DIR}/include/psi/config.hpp @ONLY)
# Setting header files
file(GLOB HEADER_FILES "${PROJECT_SOURCE_DIR}/include/psi/*.hpp")
# Creating an INTERFACE library
add_library(psi INTERFACE)
# Defining target sources
target_sources(psi INTERFACE "$<BUILD_INTERFACE:${HEADER_FILES}>")
# Defining include directories
target_include_directories(psi
  INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  SYSTEM INTERFACE $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
  INTERFACE gum::gum
  INTERFACE sdsl::sdsl
  INTERFACE SeqAn::SeqAn)
# Defining link libraries
target_link_libraries(psi
  INTERFACE gum::gum
  INTERFACE sdsl::sdsl
  INTERFACE SeqAn::SeqAn
  INTERFACE "atomic")
# Define `PSI_DEBUG`
if(CMAKE_BUILD_TYPE STREQUAL "DEBUG")
  target_compile_definitions(psi INTERFACE PSI_DEBUG)
endif()
# Use C++17
target_compile_features(psi INTERFACE cxx_std_17)
# Defining `psi::psi` alias
add_library(psi::psi ALIAS psi)
# Install targets
install(DIRECTORY include/psi DESTINATION include PATTERN "*.hpp.in" EXCLUDE)

# Adding `psikt` submodule
if(BUILD_PSIKT)
  add_subdirectory(src)
endif(BUILD_PSIKT)

# Adding auxiliary tools submodule
if(BUILD_AUX_TOOLS)
  add_subdirectory(tools)
endif(BUILD_AUX_TOOLS)

# Adding test submodule
if(BUILD_TESTING)
  add_subdirectory(test)
endif(BUILD_TESTING)

# Registering tests.
include(CTest)
add_test(NAME TestUtils COMMAND Tests "[utils]")
add_test(NAME TestVarGraph COMMAND Tests "[graph]")
add_test(NAME TestStat COMMAND Tests "[stat]")
add_test(NAME TestSequence COMMAND Tests "[sequence]")
add_test(NAME TestTraverser COMMAND Tests "[traverser]")
add_test(NAME TestIndexIter COMMAND Tests "[index][iterator]")
add_test(NAME TestFMIndex COMMAND Tests "[fmindex]")
add_test(NAME TestPath COMMAND Tests "[path]")
add_test(NAME TestPathSet COMMAND Tests "[pathset]")
add_test(NAME TestPathIndex COMMAND Tests "[pathindex]")
add_test(NAME TestMapper COMMAND Tests "[mapper]")
