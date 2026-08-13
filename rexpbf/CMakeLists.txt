cmake_minimum_required(VERSION 3.20)
project(rexpbf_bn254 VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(FetchContent)
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

add_library(rexpbf
  src/internal/crypto.cpp
  src/internal/protocol.cpp
  src/serialization.cpp
  src/pairing.cpp
  src/setup.cpp
  src/transcript.cpp
  src/proof.cpp
  src/prove.cpp
  src/gt_multiexp.cpp
  src/validation.cpp
  src/verify_reference.cpp
  src/internal/symbolic.cpp
  src/verify.cpp)
target_include_directories(rexpbf PUBLIC include)
target_compile_features(rexpbf PUBLIC cxx_std_20)
target_link_libraries(rexpbf PUBLIC mcl::mcl_st)
if(MSVC)
  target_compile_options(rexpbf PRIVATE /W4 /permissive- /EHsc)
else()
  target_compile_options(rexpbf PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_executable(run_rexpbf bench/run_rexpbf.cpp)
target_link_libraries(run_rexpbf PRIVATE rexpbf)
if(MSVC)
  target_compile_options(run_rexpbf PRIVATE /W4 /permissive- /EHsc)
endif()

include(CTest)
if(BUILD_TESTING)
  add_executable(rexpbf_tests tests/rexpbf_tests.cpp)
  target_link_libraries(rexpbf_tests PRIVATE rexpbf)
  if(MSVC)
    target_compile_options(rexpbf_tests PRIVATE /W4 /permissive- /EHsc)
  else()
    target_compile_options(rexpbf_tests PRIVATE -Wall -Wextra -Wpedantic)
  endif()
  add_test(NAME rexpbf_tests COMMAND rexpbf_tests)
endif()
