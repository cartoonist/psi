# When SeqAn is not found
if(NOT SeqAn_FOUND)
  if(NOT USE_BUNDLED_SEQAN)
    message(FATAL_ERROR "SeqAn library not found. "
      "Pass in `-DUSE_BUNDLED_SEQAN=on` when running cmake to use the bundled version. "
      "It will be installed alongside the library.")
  endif()
  message(STATUS "Using bundled SeqAn library")
  set(SEQAN_SOURCE_DIR ${PROJECT_SOURCE_DIR}/ext/seqan)
  set(SEQAN_BUILD_SYSTEM "SEQAN_RELEASE_LIBRARY" CACHE STRING "Only install SeqAn library")
  set(SEQAN_NO_DOX ON CACHE BOOL "No dox for SeqAn")
  execute_process(
    COMMAND git submodule update --init --recursive -- ${SEQAN_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  set(CMAKE_INCLUDE_PATH ${SEQAN_SOURCE_DIR}/include ${CMAKE_INCLUDE_PATH})
  set(CMAKE_PREFIX_PATH ${SEQAN_SOURCE_DIR}/util/cmake ${CMAKE_PREFIX_PATH})
  set(CMAKE_MODULE_PATH ${SEQAN_SOURCE_DIR}/util/cmake ${CMAKE_MODULE_PATH})
  find_package(SeqAn REQUIRED CONFIG)
endif()

include(SeqAnTarget)

# When `kseq++` is not found
if(NOT TARGET kseq++::kseq++)
  if(NOT USE_BUNDLED_KSEQPP)
    message(FATAL_ERROR "kseq++ library not found. "
      "Pass in `-DUSE_BUNDLED_KSEQPP=on` when running cmake to use the bundled version. "
      "It will be installed alongside the library.")
  endif()
  message(STATUS "Using bundled kseq++ library")
  set(KSEQPP_SOURCE_DIR ${PROJECT_SOURCE_DIR}/ext/kseq++)
  execute_process(
    COMMAND git submodule update --init --recursive -- ${KSEQPP_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  set(BUILD_TESTING_SAVED "${BUILD_TESTING}")
  set(BUILD_TESTING OFF)
  add_subdirectory(${KSEQPP_SOURCE_DIR})
  set(BUILD_TESTING "${BUILD_TESTING_SAVED}")
endif()

# When `diverg` is not found
# NOTE: PSI uses gum library provided transitively by DiVerG.
if(NOT TARGET diverg::diverg)
  if(NOT USE_BUNDLED_DIVERG)
    message(FATAL_ERROR "DiVerG library not found. "
      "Pass in `-DUSE_BUNDLED_DIVERG=on` when running cmake to use the bundled version. "
      "It will be installed alongside the library.")
  endif()
  message(STATUS "Using bundled DiVerG library")
  set(DIVERG_SOURCE_DIR ${PROJECT_SOURCE_DIR}/ext/diverg)
  execute_process(
    COMMAND git submodule update --init --recursive -- ${DIVERG_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  set(BUILD_TESTING_SAVED "${BUILD_TESTING}")
  set(BUILD_TESTING OFF)
  # Always build DiVerG with ETI: the header-only variant might fail under nvcc.
  set(DIVERG_ETI ON)
  # PSI does not use DiVerG's KokkosKernels-backed utilities, so keep them (and the
  # KokkosKernels dependency) out of the build.
  set(DIVERG_ENABLE_UTILS OFF)
  # Request the gum features PSI needs so DiVerG builds its bundled gum with them.
  set(GUM_WITH_VG ON)
  set(GUM_WITH_VGIO OFF)
  set(GUM_WITH_HG ON)
  set(GUM_WITH_BDSG OFF)
  # `DIVERG_ENABLE_OPENMP`/`DIVERG_ENABLE_CUDA`/`DIVERG_STATS` are forwarded via the
  # PSI cache variables of the same name. DiVerG handles its Kokkos and gum deps.
  add_subdirectory(${DIVERG_SOURCE_DIR})
  set(BUILD_TESTING "${BUILD_TESTING_SAVED}")
endif()
