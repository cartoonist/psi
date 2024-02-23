# If found, but `SeqAn::SeqAn` target is not defined
if(NOT TARGET SeqAn::SeqAn)
  add_library(SeqAn::SeqAn INTERFACE IMPORTED)
  #string(REGEX REPLACE "-D" "" SEQAN_DEFINITIONS_NOD "${SEQAN_DEFINITIONS}")
  set_target_properties(SeqAn::SeqAn PROPERTIES
    #INTERFACE_COMPILE_DEFINITIONS "${SEQAN_DEFINITIONS_NOD}",
    INTERFACE_INCLUDE_DIRECTORIES "${SEQAN_INCLUDE_DIRS}")
  # NOTE: In order to handle debug/optimized keywords, `${SEQAN_LIBRARIES}`
  # should be passed to `target_link_libraries` (with no quotation).
  target_link_libraries(SeqAn::SeqAn INTERFACE ${SEQAN_LIBRARIES})

  message(STATUS "Imported SeqAn target (version ${SEQAN_VERSION_STRING})")
endif()
