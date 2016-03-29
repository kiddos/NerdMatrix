set(CMAKE_CXX_FLAGS "-g -Wall")

add_executable(test1 mattest1.cpp common.cpp)
add_executable(test2 mattest2.cpp common.cpp)
add_executable(test3 mattestnn.cpp)
target_link_libraries(test2 armadillo)

add_test(NAME constructortest COMMAND ./test1)
add_test(NAME optest COMMAND ./test2)
