#ifndef MAT_IO_H
#define MAT_IO_H

#include <iostream>
#include <fstream>
#include "mat.h"

template <class T> std::ostream &operator<< (std::ostream &out, core::mat<T> m) {
  for (uint32_t i = 0 ; i < m.nrows ; ++i) {
    for (uint32_t j = 0 ; j < m.ncols ; ++j) {
      out << m.data[i * m.ncols + j] << " ";
    }
    out << std::endl;
  }
  return out;
}

#endif /* end of include guard: MAT_IO_H */
