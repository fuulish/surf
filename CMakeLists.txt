cmake_minimum_required(VERSION 2.6)
cmake_policy(SET CMP0012 NEW)

project (surf)

option(ENABLE_OPENMP "Use OpenMP parallelization" ON)
find_package(OpenMP)
if (OPENMP_FOUND AND ENABLE_OPENMP)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    add_definitions(-DOPENMP)
endif()

option(ENABLE_GSL "Use GNU Scientific Library" ON)
find_package (GSL)
if (GSL_FOUND  AND ENABLE_GSL)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${GSL_C_FLAGS}")
    add_definitions(-DHAVE_GSL)
endif()

find_package(PkgConfig REQUIRED)

option(ENABLE_NLOPT "Use NLOPT library" OFF)
pkg_search_module(NLOPT nlopt)
if (NLOPT_FOUND AND ENABLE_NLOPT)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${NLOPT_C_FLAGS}")
    add_definitions(-DHAVE_NLOPT)
endif()

option(ENABLE_EINSPLINE "Use EINSPLINE library" OFF)
pkg_search_module(EINSPLINE einspline)
if (EINSPLINE_FOUND AND ENABLE_EINSPLINE)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${EINSPLINE_C_FLAGS}")
    add_definitions(-DHAVE_EINSPLINE)
endif()

set(CMAKE_C_FLAGS "-O2")

add_subdirectory(src)
