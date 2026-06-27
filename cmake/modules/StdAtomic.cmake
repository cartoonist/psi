# Check if std::atomic needs -latomic
# Source: https://github.com/ccache/ccache/blob/master/cmake/StdAtomic.cmake

include(CheckCXXSourceCompiles)

set(
  check_std_atomic_source_code
  [=[
    #include <atomic>
    int main()
    {
      std::atomic<long long> x;
      ++x;
      (void)x.load();
      return 0;
    }
  ]=])

check_cxx_source_compiles("${check_std_atomic_source_code}" std_atomic_without_libatomic)

if(NOT std_atomic_without_libatomic)
  set(CMAKE_REQUIRED_LIBRARIES atomic)
  check_cxx_source_compiles("${check_std_atomic_source_code}" std_atomic_with_libatomic)
  set(CMAKE_REQUIRED_LIBRARIES)
endif()

function(target_link_atomic A_TARGET SCOPE)
  if(NOT std_atomic_without_libatomic)
    if(NOT std_atomic_with_libatomic)
      message(FATAL_ERROR "Toolchain doesn't support std::atomic with nor without -latomic")
    else()
      message(VERBOSE "Toolchain needs -latomic")
      target_link_libraries(${A_TARGET} ${SCOPE} atomic)
    endif()
  else()
    message(VERBOSE "Toolchain does not need -latomic")
  endif()
endfunction()
