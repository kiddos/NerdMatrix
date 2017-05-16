cmake_minimum_required(VERSION 2.8)

cmake_policy(SET CMP0048 NEW)
project(NerdMatrix VERSION "0.1.0" LANGUAGES CXX)

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

enable_testing()
add_subdirectory(src)
add_subdirectory(test)
