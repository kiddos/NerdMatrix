#ifndef MAT_STATIC_H
#define MAT_STATIC_H

#include <stdlib.h>
#include <math.h>
#include "mat.h"

namespace core {
  mat<double> randn(const uint32_t nrows, const uint32_t ncols) {
    mat<double> ran(nrows, ncols);
    double sum = 0, dev = 0;
    for (uint32_t i = 0 ; i < nrows * ncols ; ++i) {
      ran.data[i] = rand() % 10000 / 10000.0;
      sum += ran.data[i];
    }
    for (uint32_t i = 0 ; i < nrows * ncols ; ++i) {
      dev += pow(ran.data[i] - sum, 2);
    }
    dev = sqrt(dev / (nrows * ncols - 1));
    for (uint32_t i = 0 ; i < nrows * ncols ; ++i) {
      ran.data[i] /= dev;
    }
    return ran;
  }

  mat<int> randu(const uint32_t nrows, const uint32_t ncols, int range) {
    mat<int> ran(nrows, ncols);
    for (uint32_t i = 0 ; i < nrows * ncols ; ++i) {
      ran.data[i] = rand() % range;
    }
    return ran;
  }

  mat<double> randf(const uint32_t nrows, const uint32_t ncols) {
    mat<double> ran(nrows, ncols);
    for (uint32_t i = 0 ; i < nrows * ncols ; ++i) {
      ran.data[i] = rand() % 10000 / 10000.0;
    }
    return ran;
  }
}

#endif /* end of include guard: MAT_STATIC_H */
