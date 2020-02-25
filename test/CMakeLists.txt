# Enabling `assert` by setting the build type to 'Debug' for test module
set(CMAKE_BUILD_TYPE "Debug")

# Prepare "Catch" library for other executables
set(CATCH_SOURCE_DIR ext/catch2)
execute_process(COMMAND git submodule update --init -- ${CATCH_SOURCE_DIR}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_subdirectory(${CATCH_SOURCE_DIR})

# Defining the configure header file
configure_file(src/test_config.hpp.in ${PROJECT_BINARY_DIR}/test/src/test_config.hpp @ONLY)

# Defining target Tests
file(GLOB TEST_SOURCES "src/*.cpp")
add_executable(Tests ${TEST_SOURCES})
target_compile_options(Tests PRIVATE -g -Wall -Wpedantic -Werror)
target_include_directories(Tests
  PRIVATE ${PROJECT_BINARY_DIR}/test/src
  PRIVATE psi::psi
  PRIVATE VGio::VGio)
target_link_libraries(Tests
  PRIVATE psi::psi
  PRIVATE VGio::VGio)
