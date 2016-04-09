set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -O3")

include (CheckIncludeFiles)

check_include_files(omp.h HAVE_OPENMP)
if (HAVE_OPENMP)
  message(STATUS "using OpenMP")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -DPARALLEL")
endif (HAVE_OPENMP)

add_executable(test1 mattest1.cpp common.cpp)
add_executable(test2 mattest2.cpp common.cpp)
add_executable(test3 mattest3.cpp)
add_executable(test4 mattest4.cpp)
add_executable(testnn mattestnn.cpp)

target_link_libraries(test2 armadillo)

add_test(NAME constructortest COMMAND ./test1)
add_test(NAME optest COMMAND ./test2)
add_test(NAME rowoptest COMMAND ./test3)
add_test(NAME optest2 COMMAND ./test4)

