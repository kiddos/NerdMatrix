#ifndef MAT_OP_H
#define MAT_OP_H

#include "mat.h"

template<class T> core::mat<T> operator+ (T v, core::mat<T> m) {
  core::mat<T> newmat(m.nrows, m.ncols);
  for (uint32_t i = 0 ; i < m.nrows * m.ncols ; ++i) {
    newmat.data[i] = v + m.data[i];
  }
  return newmat;
}

template<class T> core::mat<T> operator- (T v, core::mat<T> m) {
  core::mat<T> newmat(m.nrows, m.ncols);
  for (uint32_t i = 0 ; i < m.nrows * m.ncols ; ++i) {
    newmat.data[i] = v - m.data[i];
  }
  return newmat;
}

template<class T> core::mat<T> operator* (T v, core::mat<T> m) {
  core::mat<T> newmat(m.nrows, m.ncols);
  for (uint32_t i = 0 ; i < m.nrows * m.ncols ; ++i) {
    newmat.data[i] = v * m.data[i];
  }
  return newmat;
}

template<class T> core::mat<T> operator/ (T v, core::mat<T> m) {
  core::mat<T> newmat(m.nrows, m.ncols);
  for (uint32_t i = 0 ; i < m.nrows * m.ncols ; ++i) {
    newmat.data[i] = v / m.data[i];
  }
  return newmat;
}

#endif /* end of include guard: MAT_OP_H */
