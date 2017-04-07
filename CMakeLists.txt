cmake_minimum_required(VERSION 2.6)
cmake_policy(SET CMP0012 NEW)

project (surf)

if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "Release")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/)
include(FindOptionalPackage)

find_optional_package(OpenMP)

if (NOT DISABLE_OPENMP AND OPENMP_FOUND)
    add_definitions(-DOPENMP)
endif()

find_optional_package(GSL)

if (NOT DISABLE_GSL AND GSL_FOUND)
    add_definitions(-DHAVE_GSL)
endif()

find_package(PkgConfig REQUIRED)

option(ENABLE_NLOPT "Use NLOPT optimization library" OFF)

if(ENABLE_NLOPT)
    pkg_search_module(NLOPT nlopt REQUIRED)
    add_definitions(-DHAVE_NLOPT)
endif()

option(ENABLE_EINSPLINE "Use EINSPLINE interpolation library" OFF)

if(ENABLE_EINSPLINE)
    pkg_search_module(EINSPLINE einspline REQUIRED)
    add_definitions(-DHAVE_EINSPLINE)
endif()

add_subdirectory(src)
