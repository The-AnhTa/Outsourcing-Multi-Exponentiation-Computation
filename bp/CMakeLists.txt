cmake_minimum_required(VERSION 3.20)
project(bp_ipa VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(BP_BUILD_BENCHMARKS "Build BP protocol benchmark executables" ON)

function(bp_enable_strict_warnings target)
  if(MSVC)
    target_compile_options(${target} PRIVATE /W4 /permissive-)
  else()
    target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
  endif()
endfunction()

include(FetchContent)
if(MSVC AND CMAKE_SIZEOF_VOID_P EQUAL 8)
  get_filename_component(_bp_msvc_bin "${CMAKE_CXX_COMPILER}" DIRECTORY)
  set(ML64 "${_bp_msvc_bin}/ml64.exe" CACHE FILEPATH "MSVC x64 assembler" FORCE)
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

add_library(bp_ipa
  src/internal/protocol_utils.cpp
  src/internal/protocol_validation.cpp
  src/internal/bp_transcript.cpp
  src/bp.cpp
  src/hp_transcript_internal.cpp
  src/hp_vme_internal.cpp
  src/helper_prover.cpp
  src/helper_verifier.cpp)
target_include_directories(bp_ipa PUBLIC include PRIVATE src)
target_link_libraries(bp_ipa PUBLIC mcl::mcl_st)
target_compile_features(bp_ipa PUBLIC cxx_std_20)
bp_enable_strict_warnings(bp_ipa)

if(BP_BUILD_BENCHMARKS)
  add_executable(run_bp bench/run_bp.cpp)
  target_link_libraries(run_bp PRIVATE bp_ipa)
  bp_enable_strict_warnings(run_bp)

  add_executable(run_helper_prover bench/run_helper_prover.cpp)
  target_link_libraries(run_helper_prover PRIVATE bp_ipa)
  bp_enable_strict_warnings(run_helper_prover)

  add_executable(run_helper_verifier bench/run_helper_verifier.cpp)
  target_link_libraries(run_helper_verifier PRIVATE bp_ipa)
  bp_enable_strict_warnings(run_helper_verifier)
endif()

include(CTest)
if(BUILD_TESTING)
  add_executable(bp_tests tests/bp_tests.cpp)
  target_include_directories(bp_tests PRIVATE src)
  target_link_libraries(bp_tests PRIVATE bp_ipa)
  bp_enable_strict_warnings(bp_tests)
  add_test(NAME bp_tests COMMAND bp_tests)
endif()
