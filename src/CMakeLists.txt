add_subdirectory(xdrfile)
add_definitions(-DHAVE_XDRFILE)

include_directories(${CMAKE_SOURCE_DIR}/include)

file(GLOB sources *.c)

add_executable(surf ${sources})

target_link_libraries(surf xdrfile m)

if (OPENMP_FOUND AND NOT DISABLE_OPENMP)
    target_link_libraries(surf ${OpenMP_LIBRARIES})
    target_include_directories(surf PUBLIC ${OpenMP_INCLUDE_DIRS})
endif()

if (GSL_FOUND AND NOT DISABLE_GSL)
    target_link_libraries(surf ${GSL_LIBRARIES})
    target_include_directories(surf PUBLIC ${GSL_INCLUDE_DIRS})
endif()

if (NLOPT_FOUND AND NOT DISABLE_NLOPT)
    target_link_libraries(surf ${NLOPT_LIBRARIES})
    target_include_directories(surf PUBLIC ${NLOPT_INCLUDE_DIRS})
endif()

if (EINSPLINE_FOUND AND NOT DISABLE_EINSPLINE)
    target_link_libraries(surf ${EINSPLINE_LIBRARIES})
    target_include_directories(surf PUBLIC ${EINSPLINE_INCLUDE_DIRS})
endif()

install(TARGETS surf DESTINATION bin)
