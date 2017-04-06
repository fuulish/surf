add_subdirectory(xdrfile)
add_definitions(-DHAVE_XDRFILE)

include_directories(${CMAKE_SOURCE_DIR}/include)

file(GLOB sources *.c)

add_executable(surf ${sources})

target_link_libraries(surf xdrfile m)

if (NLOPT_FOUND)
    target_link_libraries(surf ${NLOPT_LIBRARIES})
    target_include_directories(surf PUBLIC ${NLOPT_INCLUDE_DIRS})
    target_compile_options(surf PUBLIC ${NLOPT_CFLAGS_OTHER})
endif()

if (EINSPLINE_FOUND)
    target_link_libraries(surf ${EINSPLINE_LIBRARIES})
    target_include_directories(surf PUBLIC ${EINSPLINE_INCLUDE_DIRS})
    target_compile_options(surf PUBLIC ${EINSPLINE_CFLAGS_OTHER})
endif()

install(TARGETS surf DESTINATION bin)
