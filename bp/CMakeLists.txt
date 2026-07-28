cmake_minimum_required(VERSION 3.20)
project(bp_ipa LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

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
  src/bp.cpp
  src/hp_transcript_internal.cpp
  src/hp_vme_internal.cpp
  src/helper_prover.cpp
  src/helper_verifier.cpp)
target_include_directories(bp_ipa PUBLIC include)
target_link_libraries(bp_ipa PUBLIC mcl::mcl_st)
target_compile_features(bp_ipa PUBLIC cxx_std_20)
if(MSVC)
  target_compile_options(bp_ipa PRIVATE /W4 /permissive-)
else()
  target_compile_options(bp_ipa PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_executable(run_bp bench/run_bp.cpp)
target_link_libraries(run_bp PRIVATE bp_ipa)
add_executable(run_helper_prover bench/run_helper_prover.cpp)
target_link_libraries(run_helper_prover PRIVATE bp_ipa)
add_executable(run_helper_verifier bench/run_helper_verifier.cpp)
target_link_libraries(run_helper_verifier PRIVATE bp_ipa)
