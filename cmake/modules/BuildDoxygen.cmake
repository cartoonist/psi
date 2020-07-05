# Source: https://vicrucann.github.io/tutorials/quick-cmake-doxygen

# Check if Doxygen is installed
find_package(Doxygen)

if(Doxygen_FOUND)
  set(DOXYGEN_FILE ${PROJECT_SOURCE_DIR}/Doxyfile)
  add_custom_target(doc
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_FILE}
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMENT "Generating API documentation with Doxygen"
    VERBATIM)
else(Doxygen_FOUND)
  message("Doxygen need to be installed to generate the doxygen documentation")
endif(Doxygen_FOUND)
