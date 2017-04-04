file(GLOB sources *.c)

add_executable(surf ${sources})

target_link_libraries(surf m)
