set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -std=c++11 -Wall -O2")

include (CheckIncludeFiles)

check_include_files(omp.h HAVE_OPENMP)
if (HAVE_OPENMP)
  message(STATUS "using OpenMP")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -DPARALLEL")
endif (HAVE_OPENMP)

add_executable(test1 mattest1.cpp)
add_executable(test2 mattest2.cpp)
add_executable(test3 mattest3.cpp)
add_executable(test4 mattest4.cpp)
add_executable(test5 mattest5.cpp)
add_executable(test6 mattest6.cpp)
add_executable(testnn mattestnn.cpp)
add_executable(testnn2 mattestnn2.cpp)

target_link_libraries(test2 armadillo)
target_link_libraries(testnn2 opencv_core opencv_highgui opencv_imgproc)

add_test(NAME constructortest COMMAND ./test1)
add_test(NAME optest COMMAND ./test2)
add_test(NAME rowoptest COMMAND ./test3)
add_test(NAME optest2 COMMAND ./test4)
add_test(NAME optest3 COMMAND ./test5)
add_test(NAME statictest COMMAND ./test6)

