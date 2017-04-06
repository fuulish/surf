add_subdirectory(xdrfile)

file(GLOB sources *.c)

add_executable(surf ${sources})

target_link_libraries(surf xdrfile m)

install(TARGETS surf DESTINATION bin)
