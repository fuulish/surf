add_subdirectory(xdrfile)
add_definitions(-DHAVE_XDRFILE)

include_directories(${CMAKE_SOURCE_DIR}/include)

file(GLOB sources *.c)

add_executable(surf ${sources})

target_link_libraries(surf xdrfile m)

install(TARGETS surf DESTINATION bin)
