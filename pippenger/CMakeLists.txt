cmake_minimum_required(VERSION 3.20)
project(multi_instance_pippenger VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(FetchContent)
if(MSVC AND CMAKE_SIZEOF_VOID_P EQUAL 8)
  get_filename_component(_msvc_bin "${CMAKE_CXX_COMPILER}" DIRECTORY)
  set(ML64 "${_msvc_bin}/ml64.exe" CACHE FILEPATH "MSVC x64 assembler" FORCE)
endif()
set(MCL_BUILD_SAMPLE OFF CACHE BOOL "" FORCE)
set(MCL_BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(MCL_TEST_WITH_GMP OFF CACHE BOOL "" FORCE)
set(MCL_USE_GMP OFF CACHE BOOL "" FORCE)
if(MSVC)
  set(MCL_MSVC_RUNTIME_DLL ON CACHE BOOL "" FORCE)
endif()

FetchContent_Declare(mcl
  URL https://github.com/herumi/mcl/archive/refs/tags/v3.00.zip
  URL_HASH SHA256=8298dcc4fe8619bd6b19597b211300e4503a14838d1ae3906c44a25242280ccb
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
FetchContent_MakeAvailable(mcl)

if(WIN32 AND TARGET mcl)
  set_target_properties(mcl PROPERTIES OUTPUT_NAME mcl_shared)
endif()

add_library(pippenger
  src/initialization.cpp
  src/pippenger.cpp
  src/pippenger_core.cpp
  src/pippenger_validation.cpp)
target_include_directories(pippenger PUBLIC include PRIVATE src)
target_link_libraries(pippenger PUBLIC mcl::mcl_st)
target_compile_features(pippenger PUBLIC cxx_std_20)

add_executable(run_pippenger bench/run_pippenger.cpp)
target_link_libraries(run_pippenger PRIVATE pippenger)

add_library(pinkas
  src/pinkas_validation.cpp
  src/protocol.cpp
  src/protocol_math.cpp
  src/prover.cpp
  src/serialization.cpp
  src/transcript.cpp
  src/verifier.cpp)
target_include_directories(pinkas PUBLIC include PRIVATE src)
target_link_libraries(pinkas PUBLIC pippenger)
target_compile_features(pinkas PUBLIC cxx_std_20)

add_executable(run_pinkas bench/run_pinkas.cpp)
target_link_libraries(run_pinkas PRIVATE pinkas)

include(CTest)
if(BUILD_TESTING)
  add_executable(pippenger_tests tests/pippenger_tests.cpp)
  target_include_directories(pippenger_tests PRIVATE src)
  target_link_libraries(pippenger_tests PRIVATE pinkas)
endif()

function(enable_strict_warnings target)
  if(MSVC)
    target_compile_options(${target} PRIVATE /W4 /permissive-)
  else()
    target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
  endif()
endfunction()

enable_strict_warnings(pippenger)
enable_strict_warnings(run_pippenger)
enable_strict_warnings(pinkas)
enable_strict_warnings(run_pinkas)
if(BUILD_TESTING)
  enable_strict_warnings(pippenger_tests)
endif()

if(BUILD_TESTING)
  add_test(NAME pippenger_tests COMMAND pippenger_tests)
endif()
