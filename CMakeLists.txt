cmake_minimum_required(VERSION 3.10)
project(psi VERSION 0.9.1 LANGUAGES CXX)

# Options
option(BUILD_PSIKT "Build PSI seeding/indexing command-line tool (`psikt`)" ON)
option(BUILD_PSI_AUX_TOOLS "Build PSI auxiliary tools" OFF)
option(BUILD_TESTING "Build test programs" OFF)
option(BUILD_DOC "Build documentation" OFF)

# Include external modules
include(GNUInstallDirs)
# Adding CMake module path.
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)
# Finding dependencies.
find_package(ZLIB REQUIRED)   # required by SeqAn
find_package(BZip2 REQUIRED)  # required by SeqAn
find_package(OpenMP)          # required by SeqAn
find_package(SeqAn QUIET)
find_package(gum QUIET)
find_package(kseq++ QUIET)
find_package(pairg QUIET)
# Include external dependencies
include(${PROJECT_SOURCE_DIR}/ext/external.cmake)

# Load `BuildDoxygen` module
if(BUILD_DOC)
  include(BuildDoxygen)
endif(BUILD_DOC)

# Load `CheckGitRevision` module
include(CheckGitRevision)

# Generating the configure header file
if(GIT_REVISION)
  set(PSI_GIT_REVISION ${GIT_REVISION})
endif(GIT_REVISION)
if(GIT_COMMIT_DATE)
  set(PSI_GIT_COMMIT_DATE ${GIT_COMMIT_DATE})
endif(GIT_COMMIT_DATE)
configure_file(include/psi/config.hpp.in ${PROJECT_SOURCE_DIR}/include/psi/config.hpp @ONLY)
# Setting header files
file(GLOB HEADER_FILES "${PROJECT_SOURCE_DIR}/include/psi/*.hpp")
# Creating an INTERFACE library
add_library(psi INTERFACE)
add_library(psi_novg INTERFACE)
# Defining target sources
target_sources(psi INTERFACE "$<BUILD_INTERFACE:${HEADER_FILES}>")
target_sources(psi_novg INTERFACE "$<BUILD_INTERFACE:${HEADER_FILES}>")
# Defining include directories
target_include_directories(psi
  INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>;$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>
  INTERFACE gum::gum
  INTERFACE kseq++::kseq++
  INTERFACE SeqAn::SeqAn)
target_include_directories(psi_novg
  INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  SYSTEM INTERFACE $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
  INTERFACE gum::gum_novg
  INTERFACE kseq++::kseq++
  INTERFACE SeqAn::SeqAn)
# Defining link libraries
target_link_libraries(psi
  INTERFACE gum::gum
  INTERFACE kseq++::kseq++
  INTERFACE SeqAn::SeqAn
  INTERFACE pairg::libpairg
  INTERFACE "atomic")
target_link_libraries(psi_novg
  INTERFACE gum::gum_novg
  INTERFACE kseq++::kseq++
  INTERFACE SeqAn::SeqAn
  INTERFACE pairg::libpairg
  INTERFACE "atomic")
# Define `PSI_DEBUG`
if(CMAKE_BUILD_TYPE STREQUAL "DEBUG")
  target_compile_definitions(psi INTERFACE PSI_DEBUG)
  target_compile_definitions(psi_novg INTERFACE PSI_DEBUG)
endif()
# Define `PSI_STATS`
if(PSI_STATS)
  target_compile_definitions(psi INTERFACE PSI_STATS)
  target_compile_definitions(psi_novg INTERFACE PSI_STATS)
endif(PSI_STATS)
# Use C++17
target_compile_features(psi INTERFACE cxx_std_17)
target_compile_features(psi_novg INTERFACE cxx_std_17)
# Defining `psi::psi` alias
add_library(psi::psi ALIAS psi)
add_library(psi::psi_novg ALIAS psi_novg)
# Generating pkgconfig file
configure_file(psi.pc.in ${PROJECT_SOURCE_DIR}/psi.pc @ONLY)
# Install targets
install(FILES ${HEADER_FILES} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/psi)
install(FILES psi.pc DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

# Adding `psikt` submodule
if(BUILD_PSIKT)
  add_subdirectory(src)
endif(BUILD_PSIKT)

# Adding auxiliary tools submodule
if(BUILD_PSI_AUX_TOOLS)
  add_subdirectory(tools)
endif(BUILD_PSI_AUX_TOOLS)

# Adding test submodule
if(BUILD_TESTING)
  add_subdirectory(test)
endif(BUILD_TESTING)

# Registering tests.
include(CTest)
add_test(NAME TestUtils COMMAND psi-tests "[utils]")
add_test(NAME TestGraph COMMAND psi-tests "[graph]")
add_test(NAME TestCRSMatrix COMMAND psi-tests "[crsmatrix]")
add_test(NAME TestGraphIter COMMAND psi-tests "[graph][iterator]")
add_test(NAME TestStats COMMAND psi-tests "[stats]")
add_test(NAME TestSequence COMMAND psi-tests "[sequence]")
add_test(NAME TestTraverser COMMAND psi-tests "[traverser]")
add_test(NAME TestIndexIter COMMAND psi-tests "[index][iterator]")
add_test(NAME TestFMIndex COMMAND psi-tests "[fmindex]")
add_test(NAME TestPath COMMAND psi-tests "[path]")
add_test(NAME TestPathSet COMMAND psi-tests "[pathset]")
add_test(NAME TestPathIndex COMMAND psi-tests "[pathindex]")
add_test(NAME TestSeedFinder COMMAND psi-tests "[seedfinder]")
