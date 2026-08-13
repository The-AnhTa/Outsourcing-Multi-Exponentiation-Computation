cmake_minimum_required(VERSION 3.20)
project(rexp_bn254 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT TARGET mcl::mcl_st)
  include(FetchContent)
  set(MCL_BUILD_SAMPLE OFF CACHE BOOL "" FORCE)
  set(MCL_BUILD_TESTING OFF CACHE BOOL "" FORCE)
  set(MCL_TEST_WITH_GMP OFF CACHE BOOL "" FORCE)
  set(MCL_USE_GMP OFF CACHE BOOL "" FORCE)
  if(MSVC)
    set(MCL_MSVC_RUNTIME_DLL ON CACHE BOOL "" FORCE)
  endif()
  FetchContent_Declare(mcl URL https://github.com/herumi/mcl/archive/refs/tags/v3.00.zip
    URL_HASH SHA256=8298dcc4fe8619bd6b19597b211300e4503a14838d1ae3906c44a25242280ccb
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
  FetchContent_MakeAvailable(mcl)
endif()
if(WIN32 AND TARGET mcl)
  set_target_properties(mcl PROPERTIES OUTPUT_NAME mcl_shared)
endif()
add_library(rexp
  src/internal/crypto.cpp
  src/internal/dory_transcript.cpp
  src/internal/rexp_transcript.cpp
  src/dory_setup.cpp
  src/dory_prove.cpp
  src/dory_serialization.cpp
  src/dory.cpp
  src/rexp_setup.cpp
  src/rexp_prove.cpp
  src/rexp_serialization.cpp
  src/rexp_validation.cpp
  src/rexp.cpp)
target_include_directories(rexp PUBLIC include)
target_link_libraries(rexp PUBLIC mcl::mcl_st)
if(MSVC)
  target_compile_options(rexp PRIVATE /W4 /permissive- /EHsc)
else()
  target_compile_options(rexp PRIVATE -Wall -Wextra -Wpedantic)
endif()
add_executable(run_rexp bench/run_rexp.cpp)
target_link_libraries(run_rexp PRIVATE rexp)
if(MSVC)
  target_compile_options(run_rexp PRIVATE /W4 /permissive- /EHsc)
else()
  target_compile_options(run_rexp PRIVATE -Wall -Wextra -Wpedantic)
endif()

include(CTest)
if(BUILD_TESTING)
  add_library(rexp_header_checks OBJECT
    tests/headers/metrics.cpp
    tests/headers/proof.cpp
    tests/headers/prove.cpp
    tests/headers/rexp.cpp
    tests/headers/serialization.cpp
    tests/headers/setup.cpp
    tests/headers/verify.cpp)
  target_include_directories(rexp_header_checks PRIVATE include)
  target_link_libraries(rexp_header_checks PRIVATE mcl::mcl_st)
  add_executable(rexp_tests tests/rexp_tests.cpp)
  target_link_libraries(rexp_tests PRIVATE rexp)
  if(MSVC)
    target_compile_options(rexp_header_checks PRIVATE /W4 /permissive- /EHsc)
    target_compile_options(rexp_tests PRIVATE /W4 /permissive- /EHsc)
  else()
    target_compile_options(rexp_header_checks PRIVATE -Wall -Wextra -Wpedantic)
    target_compile_options(rexp_tests PRIVATE -Wall -Wextra -Wpedantic)
  endif()
  add_test(NAME rexp_tests COMMAND rexp_tests)
endif()
