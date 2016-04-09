cmake_minimum_required(VERSION 2.8)

enable_testing()

add_subdirectory(src)
add_subdirectory(test)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "/usr/share/cmake-${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}/Modules/")

message(STATUS "CMAKE_MODULE_PATH: " ${CMAKE_MODULE_PATH})

