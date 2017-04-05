cmake_minimum_required(VERSION 2.6)

project (surf)

option(ENABLE_OPENMP "Use OpenMP parallelization" ON)
find_package(OpenMP)
if (OPENMP_FOUND AND ENABLE_OPENMP)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    add_definitions(-DOPENMP)
endif()

option(ENABLE_OPTSURF "Use optimized surface calculation" OFF)
if (ENABLE_OPENMP)
    add_definitions(-DOPTSURF)
endif()

add_subdirectory(src)
