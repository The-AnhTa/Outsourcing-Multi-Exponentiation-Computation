cmake_minimum_required(VERSION 3.20)
project(vpip_bf LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(FetchContent)
if(MSVC AND CMAKE_SIZEOF_VOID_P EQUAL 8)
  get_filename_component(_vme_msvc_bin "${CMAKE_CXX_COMPILER}" DIRECTORY)
  set(ML64 "${_vme_msvc_bin}/ml64.exe" CACHE FILEPATH "MSVC x64 assembler" FORCE)
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

add_library(vpip_bf
  src/group_utils.cpp src/transcript.cpp src/setup.cpp src/phase1.cpp src/phase2.cpp
  src/proof.cpp src/protocol.cpp src/verify_reference.cpp)
target_sources(vpip_bf PRIVATE src/symbolic_gt.cpp src/batch_inversion.cpp
  src/gt_multiexp.cpp src/verify_deferred.cpp src/serialization.cpp)
target_include_directories(vpip_bf PUBLIC include)
target_link_libraries(vpip_bf PUBLIC mcl::mcl_st)
target_compile_features(vpip_bf PUBLIC cxx_std_20)
if(MSVC)
  target_compile_options(vpip_bf PRIVATE /W4 /permissive-)
else()
  target_compile_options(vpip_bf PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_executable(run_vpip bench/run_vpip.cpp)
target_link_libraries(run_vpip PRIVATE vpip_bf)
