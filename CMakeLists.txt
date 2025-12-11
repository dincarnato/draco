cmake_minimum_required(VERSION 3.8)
project(draco VERSION 1.3.0)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")
subdirs(src test)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "10")
    message(FATAL_ERROR "Insufficient gcc version. Minimal supported version for GCC is 10.")
  endif()
endif()

if(POLICY CMP0167)
  cmake_policy(SET CMP0167 OLD)
endif()

include(FetchContent)

find_package(Boost 1.66 REQUIRED)
find_package(Armadillo REQUIRED)
find_package(OpenMP REQUIRED)

option(ARMA_NO_WRAPPER "Do not use Armadillo wrappers. Increases speed." OFF)
option(LINK_TIME_OPTIMIZATIONS "Optimize at link-time. Only set for Release building." OFF)
option(NATIVE_BUILD "Build for native arch. This will probabily let everything crash if run on other computers" OFF)
option(USE_ARMA_OPENMP "Uses Armadillo OpenMP support" OFF)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  option(USE_LIBCXX "Use libc++ instead of libstdc++" OFF)
  if (USE_LIBCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
    link_libraries(c++abi)
  else()
    link_libraries(stdc++fs)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  link_libraries(stdc++fs)
endif()

find_package(PkgConfig REQUIRED)

FIND_PACKAGE(dlib)
if(NOT dlib_FOUND)
    message(STATUS "dlib not found using CMake packages. Trying pkgconfig.")
    pkg_check_modules(dlib REQUIRED dlib-1)
    set(dlib_USES_PKGCONFIG ON)
endif()

pkg_check_modules(TBB tbb REQUIRED)

FetchContent_Declare(
  fmt
  GIT_REPOSITORY https://github.com/fmtlib/fmt
  GIT_TAG        12.0.0)

message(STATUS "Fetching fmt library")
FetchContent_MakeAvailable(fmt)

if(ARMA_NO_WRAPPER)
    add_definitions(-DARMA_DONT_USE_WRAPPER)
    find_library(OPENBLAS_LIB openblas)
    if(NOT OPENBLAS_LIB)
        message(SEND_ERROR "Armadillo wrappers not set OpenBlas is not found. Please provide the OpenBlas library or unset ARMA_NO_WRAPPER")
    else()
        message(STATUS "Openblas library found: ${OPENBLAS_LIB}")
        set(ARMADILLO_LIBRARIES "-lopenblas")
    endif()
endif()

if(NOT USE_ARMA_OPENMP)
    add_definitions(-DARMA_DONT_USE_OPENMP)
endif()

if(LINK_TIME_OPTIMIZATIONS)
    set_property(GLOBAL PROPERTY INTERPROCEDURAL_OPTIMIZATION True)

    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -flto")
    endif()
endif()

if(NATIVE_BUILD)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

enable_testing()

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED 23)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-missing-braces ${OpenMP_CXX_FLAGS}")
