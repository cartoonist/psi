# Setting the build type if it is not already set
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

# Prepare "spdlog" dependency
if(NOT TARGET spdlog::spdlog)
  set(SPDLOG_SOURCE_DIR ${PROJECT_SOURCE_DIR}/src/ext/spdlog)
  execute_process(COMMAND git submodule update --init -- ${SPDLOG_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  add_subdirectory(${SPDLOG_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/ext/spdlog EXCLUDE_FROM_ALL)
endif()

# Prepare "cxxopts" dependency
if(NOT TARGET cxxopts::cxxopts)
  set(CXXOPTS_SOURCE_DIR ext/cxxopts)
  execute_process(COMMAND git submodule update --init -- ${CXXOPTS_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  add_subdirectory(${CXXOPTS_SOURCE_DIR} EXCLUDE_FROM_ALL)
endif()

# Defining target 'dindexctl': PSI distance index controller tool
set(DINDEXCTL_SOURCES "src/dindexctl.cpp")
add_executable(dindexctl ${DINDEXCTL_SOURCES})
target_compile_options(dindexctl PRIVATE -g -Wall -Wpedantic -Werror)
target_include_directories(dindexctl
  PRIVATE psi::psi
  PRIVATE cxxopts::cxxopts
  PRIVATE spdlog::spdlog_header_only)
target_link_libraries(dindexctl
  PRIVATE psi::psi
  PRIVATE cxxopts::cxxopts
  PRIVATE spdlog::spdlog_header_only)
# Install targets
install(TARGETS dindexctl DESTINATION bin)
