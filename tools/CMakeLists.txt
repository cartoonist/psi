# Include external modules
include(GNUInstallDirs)

if(CMAKE_VERSION VERSION_LESS "3.18.0")
  set(PSI_LD_HAS_LD_CLASSIC FALSE)
else()
  include(CheckLinkerFlag)
  check_linker_flag(CXX -Wl,-ld_classic PSI_LD_HAS_LD_CLASSIC)
endif()

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

# Default compile option
set(PSI_TOOLS_DEFAULT_CXXOPS -g -Wall -Wpedantic -Werror)
# Default link option
set(PSI_TOOLS_DEFAULT_LDOPS)
if (PSI_LD_HAS_LD_CLASSIC)
  list(APPEND PSI_TOOLS_DEFAULT_LDOPS -Wl,-ld_classic)
endif (PSI_LD_HAS_LD_CLASSIC)

if (EASY_PSI_EASY)
  list(REMOVE_ITEM PSI_TOOLS_DEFAULT_CXXOPS -Werror)
endif(EASY_PSI_EASY)

# Defining target 'ggsim': PSI read simulation tool
set(GGSIM_SOURCES "src/ggsim.cpp")
add_executable(ggsim ${GGSIM_SOURCES})
target_compile_options(ggsim PRIVATE ${PSI_TOOLS_DEFAULT_CXXOPS})
target_link_options(ggsim PRIVATE ${PSI_TOOLS_DEFAULT_LDOPS})
target_include_directories(ggsim
  PRIVATE psi::psi
  PRIVATE cxxopts::cxxopts
  PRIVATE spdlog::spdlog_header_only
  PRIVATE psi::minivgio)
target_link_libraries(ggsim
  PRIVATE psi::psi
  PRIVATE cxxopts::cxxopts
  PRIVATE spdlog::spdlog_header_only
  PRIVATE psi::minivgio)
# Install targets
install(TARGETS ggsim DESTINATION ${CMAKE_INSTALL_BINDIR})

# Defining target 'dindexctl': PSI distance index controller tool
set(DINDEXCTL_SOURCES "src/dindexctl.cpp")
add_executable(dindexctl ${DINDEXCTL_SOURCES})
target_compile_options(dindexctl PRIVATE ${PSI_TOOLS_DEFAULT_CXXOPS})
target_link_options(dindexctl PRIVATE ${PSI_TOOLS_DEFAULT_LDOPS})
target_include_directories(dindexctl
  PRIVATE psi::psi
  PRIVATE cxxopts::cxxopts
  PRIVATE spdlog::spdlog_header_only
  PRIVATE psi::minivgio)
target_link_libraries(dindexctl
  PRIVATE psi::psi
  PRIVATE cxxopts::cxxopts
  PRIVATE spdlog::spdlog_header_only
  PRIVATE psi::minivgio)
# Install targets
install(TARGETS dindexctl DESTINATION ${CMAKE_INSTALL_BINDIR})

# Defining target 'pindexctl': PSI path index controller tool
set(PINDEXCTL_SOURCES "src/pindexctl.cpp")
add_executable(pindexctl ${PINDEXCTL_SOURCES})
target_compile_options(pindexctl PRIVATE ${PSI_TOOLS_DEFAULT_CXXOPS})
target_link_options(pindexctl PRIVATE ${PSI_TOOLS_DEFAULT_LDOPS})
target_include_directories(pindexctl
  PRIVATE psi::psi
  PRIVATE cxxopts::cxxopts
  PRIVATE spdlog::spdlog_header_only
  PRIVATE psi::minivgio)
target_link_libraries(pindexctl
  PRIVATE psi::psi
  PRIVATE cxxopts::cxxopts
  PRIVATE spdlog::spdlog_header_only
  PRIVATE psi::minivgio)
# Install targets
install(TARGETS pindexctl DESTINATION ${CMAKE_INSTALL_BINDIR})

# Defining target 'loci_stats': PSI starting loci statistics reporting tool
set(LOCI_STATS_SOURCES "src/loci_stats.cpp")
add_executable(loci_stats ${LOCI_STATS_SOURCES})
target_compile_options(loci_stats PRIVATE ${PSI_TOOLS_DEFAULT_CXXOPS})
target_link_options(loci_stats PRIVATE ${PSI_TOOLS_DEFAULT_LDOPS})
target_include_directories(loci_stats
  PRIVATE psi::psi
  PRIVATE cxxopts::cxxopts
  PRIVATE spdlog::spdlog_header_only
  PRIVATE psi::minivgio)
target_link_libraries(loci_stats
  PRIVATE psi::psi
  PRIVATE cxxopts::cxxopts
  PRIVATE spdlog::spdlog_header_only
  PRIVATE psi::minivgio)
# Install targets
install(TARGETS loci_stats DESTINATION ${CMAKE_INSTALL_BINDIR})

# Defining target 'sloci': Starting loci inspection tool
set(SLOCI_SOURCES "src/sloci.cpp")
add_executable(sloci ${SLOCI_SOURCES})
target_compile_options(sloci PRIVATE ${PSI_TOOLS_DEFAULT_CXXOPS})
target_link_options(sloci PRIVATE ${PSI_TOOLS_DEFAULT_LDOPS})
target_include_directories(sloci
  PRIVATE psi::psi
  PRIVATE cxxopts::cxxopts
  PRIVATE spdlog::spdlog_header_only
  PRIVATE psi::minivgio)
target_link_libraries(sloci
  PRIVATE psi::psi
  PRIVATE cxxopts::cxxopts
  PRIVATE spdlog::spdlog_header_only
  PRIVATE psi::minivgio)
# Install targets
install(TARGETS sloci DESTINATION ${CMAKE_INSTALL_BINDIR})

# Defining target 'alice': PSI ALignment InspeCtor and analysEr
set(ALICE_SOURCES "src/alice.cpp")
add_executable(alice ${ALICE_SOURCES})
target_compile_options(alice PRIVATE ${PSI_TOOLS_DEFAULT_CXXOPS})
target_link_options(alice PRIVATE ${PSI_TOOLS_DEFAULT_LDOPS})
target_include_directories(alice
  PRIVATE psi::psi
  PRIVATE cxxopts::cxxopts
  PRIVATE spdlog::spdlog_header_only
  PRIVATE psi::minivgio)
target_link_libraries(alice
  PRIVATE psi::psi
  PRIVATE cxxopts::cxxopts
  PRIVATE spdlog::spdlog_header_only
  PRIVATE psi::minivgio)
# Install targets
install(TARGETS alice DESTINATION ${CMAKE_INSTALL_BINDIR})
