cmake_minimum_required(VERSION 3.10)
project(psi VERSION 0.9.1 LANGUAGES CXX)

# Options
option(BUILD_PSIKT "Build PSI seeding/indexing command-line tool (`psikt`)" ON)
option(BUILD_PSI_AUX_TOOLS "Build PSI auxiliary tools" OFF)
option(BUILD_TESTING "Build test programs" OFF)
option(BUILD_DOC "Build documentation" OFF)
option(USE_BUNDLED_SEQAN "Use bundled SeqAn library" OFF)
option(USE_BUNDLED_GUM "Use bundled gum library" OFF)
option(USE_BUNDLED_KSEQPP "Use bundled kseq++ library" OFF)
option(USE_BUNDLED_PAIRG "Use bundled PairG library" OFF)
option(USE_BUNDLED_ALL "Use all bundled dependencies" OFF)

if(USE_BUNDLED_ALL)
  set(USE_BUNDLED_SEQAN on)
  set(USE_BUNDLED_GUM on)
  set(USE_BUNDLED_KSEQPP on)
  set(USE_BUNDLED_PAIRG on)
endif(USE_BUNDLED_ALL)

# Include external modules
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
# Adding CMake module path.
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)
# Finding dependencies.
find_package(ZLIB REQUIRED)   # required by SeqAn
find_package(BZip2 REQUIRED)  # required by SeqAn
find_package(OpenMP)          # required by SeqAn
# Bundled dependencies
find_package(SeqAn QUIET)
find_package(gum 0.1.1 QUIET)
find_package(kseq++ 0.2.3 QUIET)
find_package(pairg 0.2.0 QUIET)
# Handle bundled dependencies
include(${PROJECT_SOURCE_DIR}/ext/external.cmake)

# Load `BuildDoxygen` module
if(BUILD_DOC)
  include(BuildDoxygen)
endif(BUILD_DOC)

# Load `CheckGitRevision` module
include(CheckGitRevision)

# Creating the INTERFACE library
add_library(psi INTERFACE)
# Defining include directories
target_include_directories(psi
  INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>;$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>
  INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/gum::gum>;$<INSTALL_INTERFACE:gum::gum>
  INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/kseq++::kseq++>;$<INSTALL_INTERFACE:kseq++::kseq++>
  INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/SeqAn::SeqAn>;$<INSTALL_INTERFACE:SeqAn::SeqAn>
  INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/pairg::libpairg>;$<INSTALL_INTERFACE:pairg::libpairg>)
# Defining link libraries
target_link_libraries(psi
  INTERFACE $<BUILD_INTERFACE:gum::gum>;$<INSTALL_INTERFACE:gum::gum>
  INTERFACE $<BUILD_INTERFACE:kseq++::kseq++>;$<INSTALL_INTERFACE:kseq++::kseq++>
  INTERFACE $<BUILD_INTERFACE:SeqAn::SeqAn>;$<INSTALL_INTERFACE:SeqAn::SeqAn>
  INTERFACE $<BUILD_INTERFACE:pairg::libpairg>;$<INSTALL_INTERFACE:pairg::libpairg>
  INTERFACE "atomic")
# Define `PSI_DEBUG`
if(CMAKE_BUILD_TYPE STREQUAL "DEBUG")
  target_compile_definitions(psi INTERFACE PSI_DEBUG)
endif()
# Define `PSI_STATS`
if(PSI_STATS)
  target_compile_definitions(psi INTERFACE PSI_STATS)
endif(PSI_STATS)
# Use C++17
target_compile_features(psi INTERFACE cxx_std_17)
# Generating the configure header file
if(GIT_REVISION)
  set(PSI_GIT_REVISION ${GIT_REVISION})
endif(GIT_REVISION)
if(GIT_COMMIT_DATE)
  set(PSI_GIT_COMMIT_DATE ${GIT_COMMIT_DATE})
endif(GIT_COMMIT_DATE)
configure_file(include/psi/config.hpp.in ${PROJECT_SOURCE_DIR}/include/psi/config.hpp @ONLY)
# Generating pkgconfig file
configure_file(psi.pc.in ${PROJECT_SOURCE_DIR}/psi.pc @ONLY)
# Setting header files
file(GLOB HEADER_FILES RELATIVE "${PROJECT_SOURCE_DIR}/include" "${PROJECT_SOURCE_DIR}/include/psi/*.hpp")
list(TRANSFORM HEADER_FILES PREPEND "${PROJECT_SOURCE_DIR}/include/" OUTPUT_VARIABLE BUILD_HEADER_FILES)
list(TRANSFORM HEADER_FILES PREPEND "${CMAKE_INSTALL_FULL_INCLUDEDIR}/" OUTPUT_VARIABLE INSTALL_HEADER_FILES)
# Defining target sources
target_sources(psi
  INTERFACE "$<BUILD_INTERFACE:${BUILD_HEADER_FILES}>;$<INSTALL_INTERFACE:${INSTALL_HEADER_FILES}>")
# Defining `psi::psi` alias
add_library(psi::psi ALIAS psi)
# Install targets
install(FILES ${BUILD_HEADER_FILES} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/psi)
install(FILES psi.pc DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
# Exporting targets
install(TARGETS psi EXPORT psi-targets)
install(EXPORT psi-targets NAMESPACE psi:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/psi)
install(FILES ${PROJECT_SOURCE_DIR}/cmake/modules/SeqAnTarget.cmake DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/psi)
# Creating the package configuration file
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/psi-config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/psi-config.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/psi)
# Generating the version file for the configuration file
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/psi-config-version.cmake"
  VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}"
  COMPATIBILITY AnyNewerVersion)
# Install generated configuration files
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/psi-config.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/psi-config-version.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/psi)

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
