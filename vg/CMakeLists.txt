# Setting the build type if it is not already set
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

if(TARGET psi::minivgio)
  message(DEBUG "Skipping listfile as `psi::minivgio` is already defined")
  return()
endif()

# Generate vg protobuf source and header
protobuf_generate_cpp(VGPROTO_SRCS VGPROTO_HDRS vg.proto)

add_library(minivgio STATIC ${VGPROTO_SRCS} ${VGPROTO_HDRS})
target_compile_options(minivgio PRIVATE -g -Wall -Wpedantic -Werror)
target_include_directories(minivgio
  PUBLIC ${PROJECT_SOURCE_DIR}  # vg/stream.hpp
  PUBLIC ${PROJECT_BINARY_DIR}  # vg/vg.pb.*
  PUBLIC protobuf::libprotobuf)
target_link_libraries(minivgio
  PUBLIC protobuf::libprotobuf)

add_library(psi::minivgio ALIAS minivgio)
