cmake_minimum_required(VERSION 3.20)
project(multi_instance_pippenger LANGUAGES CXX)

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

add_library(pippenger src/pippenger.cpp)
target_include_directories(pippenger PUBLIC include)
target_link_libraries(pippenger PUBLIC mcl::mcl_st)

add_executable(run_pippenger run_pippenger.cpp)
target_link_libraries(run_pippenger PRIVATE pippenger)

add_library(pinkas src/pinkas.cpp)
target_include_directories(pinkas PUBLIC include)
target_link_libraries(pinkas PUBLIC pippenger)

add_executable(run_pinkas run_pinkas.cpp)
target_link_libraries(run_pinkas PRIVATE pinkas)

if(MSVC)
  target_compile_options(pippenger PRIVATE /W4 /permissive-)
  target_compile_options(run_pippenger PRIVATE /W4 /permissive-)
  target_compile_options(pinkas PRIVATE /W4 /permissive-)
  target_compile_options(run_pinkas PRIVATE /W4 /permissive-)
else()
  target_compile_options(pippenger PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(run_pippenger PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(pinkas PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(run_pinkas PRIVATE -Wall -Wextra -Wpedantic)
endif()
