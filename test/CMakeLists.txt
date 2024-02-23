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

# Default compile option
set(PSI_TESTS_DEFAULT_CXXOPS -g -Wall -Wpedantic -Werror)
# Default link option
set(PSI_TESTS_DEFAULT_LDOPS -Wl,-ld_classic)

if (EASY_PSI_EASY)
  list(REMOVE_ITEM PSI_TESTS_DEFAULT_CXXOPS -Werror)
endif(EASY_PSI_EASY)

# Defining target psi-tests
file(GLOB TEST_SOURCES "src/*.cpp")
add_executable(psi-tests ${TEST_SOURCES})
target_compile_options(psi-tests PRIVATE ${PSI_TESTS_DEFAULT_CXXOPS})
target_link_options(psi-tests PRIVATE ${PSI_TESTS_DEFAULT_LDOPS})
target_include_directories(psi-tests
  PRIVATE ${PROJECT_BINARY_DIR}/test/src
  PRIVATE psi::psi
  PRIVATE psi::minivgio)
target_link_libraries(psi-tests
  PRIVATE Catch2::Catch2
  PRIVATE psi::psi
  PRIVATE psi::minivgio)
