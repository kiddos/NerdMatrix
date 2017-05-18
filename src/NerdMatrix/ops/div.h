#ifndef DIV_H
#define DIV_H

#include "NerdMatrix/except/nerdmatrix_illegal_ops.h"
#include "NerdMatrix/matrix.h"
#include "NerdMatrix/wrapper/blas_wrapper.h"

namespace nerd {
namespace matrix {
namespace ops {

template <typename T>
NerdMatrix<T> operator/(const NerdMatrix<T>& m1, const NerdMatrix<T>& m2) {
  if (m1.rows() != m2.rows() || m1.cols() != m2.cols()) {
    throw except::IllegalOpsException(
        "Matrix Element-wise Division: Input Size does not match.");
  } else {
    NerdMatrix<T> result(m1.rows(), m1.cols());
    int n = m1.rows() * m1.cols();
    for (int i = 0; i < n; ++i) {
      result[i] = m1.data()[i] / m2.data()[i];
    }
    return result;
  }
}

template <typename T, typename U>
NerdMatrix<T> operator/(const NerdMatrix<T>& m, U val) {
  NerdMatrix<T> result(m.rows(), m.cols());
  int n = m.rows() * m.cols();
  T v = static_cast<T>(val);
  for (int i = 0; i < n; ++i) {
    result[i] = m.data()[i] / v;
  }
  return result;
}

template <typename T, typename U>
void Div(const NerdMatrix<T>& m, U u, NerdMatrix<T>& result) {
  result = m / u;
}

} /* end of ops namespace */
} /* end of matrix namespace */
} /* end of nerd namespace */

#endif /* end of include guard: DIV_H */
