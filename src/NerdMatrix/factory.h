#ifndef FACTORY_H
#define FACTORY_H

#include <chrono>
#include <random>
#include "NerdMatrix/except/nerdmatrix_exception.h"
#include "NerdMatrix/matrix.h"

namespace nerd {
namespace matrix {
namespace factory {

template <typename T>
NerdMatrix<T> Zeros(int rows, int cols) {
  return NerdMatrix<T>(rows, cols);
}

template <typename T>
NerdMatrix<T> Ones(int rows, int cols) {
  NerdMatrix<T> m(rows, cols);
  int n = rows * cols;
  for (int i = 0; i < n; ++i) {
    m[i] = static_cast<T>(1);
  }
  return m;
}

template <typename T, typename U>
NerdMatrix<T> Const(int rows, int cols, U val) {
  NerdMatrix<T> m(rows, cols);
  int n = rows * cols;
  for (int i = 0; i < n; ++i) {
    m[i] = static_cast<T>(val);
  }
  return m;
}

template <typename T>
NerdMatrix<T> Identity(int n) {
  NerdMatrix<T> m(n, n);
  for (int i = 0; i < n; ++i) {
    m[i * n + i] = static_cast<T>(1);
  }
  return m;
}

template <typename T>
NerdMatrix<T> Randn(int rows, int cols, double mean, double variance) {
  std::mt19937 gen(std::chrono::system_clock::now().time_since_epoch().count());
  std::normal_distribution<double> dist(mean, variance);

  NerdMatrix<T> m(rows, cols);
  int n = rows * cols;
  for (int i = 0; i < n; ++i) {
    m[i] = static_cast<T>(dist(gen));
  }
  return m;
}

template <typename T>
NerdMatrix<T> Randn(int rows, int cols, double mean, double variance,
                    int seeds) {
  std::mt19937 gen(seeds);
  std::normal_distribution<double> dist(mean, variance);

  NerdMatrix<T> m(rows, cols);
  int n = rows * cols;
  for (int i = 0; i < n; ++i) {
    m[i] = static_cast<T>(dist(gen));
  }
  return m;
}

template <typename T>
NerdMatrix<T> Randu(int rows, int cols, double minval, double maxval) {
  std::mt19937 gen(std::chrono::system_clock::now().time_since_epoch().count());
  std::uniform_real_distribution<double> dist(minval, maxval);

  NerdMatrix<T> m(rows, cols);
  int n = rows * cols;
  for (int i = 0; i < n; ++i) {
    m[i] = static_cast<T>(dist(gen));
  }
  return m;
}

template <typename T>
NerdMatrix<T> Randu(int rows, int cols, double minval, double maxval,
                    int seeds) {
  std::mt19937 gen(seeds);
  std::uniform_real_distribution<double> dist(minval, maxval);

  NerdMatrix<T> m(rows, cols);
  int n = rows * cols;
  for (int i = 0; i < n; ++i) {
    m[i] = static_cast<T>(dist(gen));
  }
  return m;
}

} /* end of factory namespace */
} /* end of matrix namespace */
} /* end of nerd namespace */

#endif /* end of include guard: FACTORY_H */
