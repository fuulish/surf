cmake_minimum_required(VERSION 2.6)

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
    add_definitions(-DHAVE_GSL)
endif()

add_subdirectory(src)
