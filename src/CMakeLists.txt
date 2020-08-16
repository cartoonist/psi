# Setting the build type if it is not already set
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

# Prepare "spdlog" dependency
if(NOT TARGET spdlog::spdlog)
  set(SPDLOG_SOURCE_DIR ext/spdlog)
  execute_process(COMMAND git submodule update --init -- ${SPDLOG_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  add_subdirectory(${SPDLOG_SOURCE_DIR} EXCLUDE_FROM_ALL)
endif()

# Defining target 'psikt': PSI indexing/seeding tool
set(SOURCES "psikt.cpp")
add_executable(psikt ${SOURCES})
target_compile_options(psikt PRIVATE -g -Wall -Wpedantic -Werror)
target_include_directories(psikt
  PRIVATE ${PROJECT_BINARY_DIR}/src
  PRIVATE psi::psi
  PRIVATE spdlog::spdlog_header_only)
target_link_libraries(psikt
  PRIVATE psi::psi
  PRIVATE spdlog::spdlog_header_only)
# Install targets
install(TARGETS psikt DESTINATION bin)
