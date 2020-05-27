# When SeqAn is not found
if(NOT SeqAn_FOUND)
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

# If found, but `SeqAn::SeqAn` target is not defined
if(NOT TARGET SeqAn::SeqAn)
  add_library(SeqAn::SeqAn INTERFACE IMPORTED)
  #string(REGEX REPLACE "-D" "" SEQAN_DEFINITIONS_NOD "${SEQAN_DEFINITIONS}")
  set_target_properties(SeqAn::SeqAn PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${SEQAN_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES "${SEQAN_LIBRARIES}")
    #INTERFACE_COMPILE_DEFINITIONS "${SEQAN_DEFINITIONS_NOD}")
endif()

# When `gum` is not found
if(NOT TARGET gum::gum)
  set(GUM_SOURCE_DIR ${PROJECT_SOURCE_DIR}/ext/gum)
  execute_process(
    COMMAND git submodule update --init --recursive -- ${GUM_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  add_subdirectory(${GUM_SOURCE_DIR})
endif()

# When `kseq++` is not found
if(NOT TARGET kseq++::kseq++)
  set(KSEQPP_SOURCE_DIR ${PROJECT_SOURCE_DIR}/ext/kseq++)
  execute_process(
    COMMAND git submodule update --init --recursive -- ${KSEQPP_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  add_subdirectory(${KSEQPP_SOURCE_DIR})
endif()
