cmake_minimum_required(VERSION 2.8)

cmake_policy(SET CMP0048 NEW)
project(NerdMatrix VERSION "0.1.0" LANGUAGES CXX)

enable_testing()
add_subdirectory(src)
add_subdirectory(test)
