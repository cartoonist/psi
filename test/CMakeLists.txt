# Enabling `assert` by setting the build type to 'Debug' for test module
set(CMAKE_BUILD_TYPE "Debug")

# Prepare "Catch" library for other executables
if(NOT TARGET Catch2::Catch2)
  set(CATCH_SOURCE_DIR ext/catch2)
  execute_process(COMMAND git submodule update --init -- ${CATCH_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  add_subdirectory(${CATCH_SOURCE_DIR})
endif()

# Defining the configure header file
configure_file(src/test_config.hpp.in ${PROJECT_BINARY_DIR}/test/src/test_config.hpp @ONLY)

# Defining target psi-tests
file(GLOB TEST_SOURCES "src/*.cpp")
add_executable(psi-tests ${TEST_SOURCES})
target_compile_options(psi-tests PRIVATE -g -Wall -Wpedantic -Werror)
target_include_directories(psi-tests
  PRIVATE ${PROJECT_BINARY_DIR}/test/src
  PRIVATE psi::psi)
target_link_libraries(psi-tests
  PRIVATE Catch2::Catch2
  PRIVATE psi::psi)
