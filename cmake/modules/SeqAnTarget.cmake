# If found, but `SeqAn::SeqAn` target is not defined
if(NOT TARGET SeqAn::SeqAn)
  add_library(SeqAn::SeqAn INTERFACE IMPORTED)
  #string(REGEX REPLACE "-D" "" SEQAN_DEFINITIONS_NOD "${SEQAN_DEFINITIONS}")
  set_target_properties(SeqAn::SeqAn PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${SEQAN_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES "${SEQAN_LIBRARIES}")
  #INTERFACE_COMPILE_DEFINITIONS "${SEQAN_DEFINITIONS_NOD}")

  message(STATUS "Imported SeqAn target (version ${SEQAN_VERSION_STRING})")
endif()
