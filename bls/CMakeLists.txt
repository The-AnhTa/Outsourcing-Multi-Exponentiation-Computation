cmake_minimum_required(VERSION 3.20)
project(bls_agg_bf VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(BLS_BUILD_BENCHMARKS "Build the BLS aggregate benchmark" ON)

function(bls_enable_strict_warnings target)
  if(MSVC)
    target_compile_options(${target} PRIVATE /W4 /permissive-)
  else()
    target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
  endif()
endfunction()

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

add_library(bls_agg_bf
  src/internal/crypto.cpp
  src/internal/validation.cpp
  src/setup.cpp
  src/transcript.cpp
  src/protocol.cpp
  src/serialization.cpp)
target_include_directories(bls_agg_bf PUBLIC include PRIVATE src)
target_link_libraries(bls_agg_bf PUBLIC mcl::mcl_st)
target_compile_features(bls_agg_bf PUBLIC cxx_std_20)
bls_enable_strict_warnings(bls_agg_bf)

if(BLS_BUILD_BENCHMARKS)
  add_executable(run_bls bench/run_bls.cpp)
  target_link_libraries(run_bls PRIVATE bls_agg_bf)
  bls_enable_strict_warnings(run_bls)
endif()

include(CTest)
if(BUILD_TESTING)
  add_executable(bls_tests tests/bls_tests.cpp)
  target_link_libraries(bls_tests PRIVATE bls_agg_bf)
  bls_enable_strict_warnings(bls_tests)
  add_test(NAME bls_tests COMMAND bls_tests)
endif()
