#ifndef SUB_H
#define SUB_H

#include "NerdMatrix/except/nerdmatrix_illegal_ops.h"
#include "NerdMatrix/matrix.h"
#include "NerdMatrix/wrapper/blas_wrapper.h"

namespace nerd {
namespace matrix {
namespace ops {

template <typename T>
void Sub(const NerdMatrix<T>& m1, const NerdMatrix<T>& m2,
         NerdMatrix<T>& result) {
  if (m1.rows() != m2.rows() || m1.cols() != m2.cols()) {
    throw except::IllegalOpsException(
        "Matrix Subtract: Input Size does not match.");
  } else {
    result = NerdMatrix<float>(m1.rows(), m1.cols());
    int n = m1.rows() * m1.cols();
    for (int i = 0; i < n; ++i) {
      result[i] = m1.data()[i] - m2.data()[i];
    }
  }
}

template <typename T>
NerdMatrix<T> operator-(const NerdMatrix<T>& m1, const NerdMatrix<T>& m2) {
  if (m1.rows() != m2.rows() || m1.cols() != m2.cols()) {
    throw except::IllegalOpsException(
        "Matrix Subtract: Input Size does not match.");
  } else {
    NerdMatrix<T> result(m1.rows(), m1.cols());
    int n = m1.rows() * m1.cols();
    for (int i = 0; i < n; ++i) {
      result[i] = m1.data()[i] - m2.data()[i];
    }
    return result;
  }
}

} /* end of ops namespace */
} /* end of matrix namespace */
} /* end of nerd namespace */

#endif /* end of include guard: SUB_H */
