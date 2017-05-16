set(CMAKE_CXX_STANDARD 11)
set(BUILD_SHARED_LIBS ON)

find_library(OPENBLAS_LIBRARIES openblas)
if (OPENBLAS_LIBRARIES)
  set(CMAKE_INCLUDE_CURRENT_DIR ON)

  file(GLOB_RECURSE NerdMatrixSource NerdMatrix/*.cc NerdMatrix/*.h)
  add_library(NerdMatrix ${NerdMatrixSource})
  message(STATUS "Linking ${BLAS_LIBRARIES}")
  target_link_libraries(NerdMatrix ${OPENBLAS_LIBRARIES})
endif ()
