cmake_minimum_required(VERSION 2.6)
cmake_policy(SET CMP0012 NEW)

project (surf)

if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "Release")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/)
include(FindOptionalPackage)
include(FindOptionalPkgPackage)

find_optional_package(OpenMP)

if (NOT DISABLE_OPENMP AND OPENMP_FOUND)
    add_definitions(-DOPENMP)
endif()

find_optional_package(GSL)

if (NOT DISABLE_GSL AND GSL_FOUND)
    add_definitions(-DHAVE_GSL)
endif()

find_package(PkgConfig REQUIRED)

find_optional_pkg_package(NLOPT)

if (NOT DISABLE_NLOPT AND NLOPT_FOUND)
    add_definitions(-DHAVE_NLOPT)
endif()

find_optional_pkg_package(EINSPLINE)

if (NOT DISABLE_EINSPLINE AND EINSPLINE_FOUND)
    add_definitions(-DHAVE_EINSPLINE)
endif()

add_subdirectory(src)
