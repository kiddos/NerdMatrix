#ifndef ADD_H
#define ADD_H

#include "NerdMatrix/matrix.h"
#include "NerdMatrix/wrapper/blas_wrapper.h"
#include "NerdMatrix/except/nerdmatrix_illegal_ops.h"

namespace nerd {
namespace matrix {
namespace ops {

template <typename T>
void Add(const NerdMatrix<T>& m1, const NerdMatrix<T>& m2,
         NerdMatrix<T>& result) {
  if (m1.rows() != m2.rows() || m1.cols() != m2.cols()) {
    throw except::IllegalOpsException("Matrix Add: Input Size does not match.");
  } else {
    result = NerdMatrix<float>(m1.rows(), m1.cols());
    int n = m1.rows() * m1.cols();
    T alpha = static_cast<T>(1);
    wrapper::Axpy(n, alpha, m1.data(), &result[0]);
    wrapper::Axpy(n, alpha, m2.data(), &result[0]);
  }
}

template <typename T>
NerdMatrix<T> operator+(const NerdMatrix<T>& m1, const NerdMatrix<T> m2) {
  if (m1.rows() != m2.rows() || m1.cols() != m2.cols()) {
    throw except::IllegalOpsException("Matrix Add: Input Size does not match.");
  } else {
    NerdMatrix<T> result(m1.rows(), m1.cols());
    int n = m1.rows() * m1.cols();
    T alpha = static_cast<T>(1);
    wrapper::Axpy(n, alpha, m1.data(), &result[0]);
    wrapper::Axpy(n, alpha, m2.data(), &result[0]);
    return result;
  }
}

} /* end of ops namespace */
} /* end of matrix namespace */
} /* end of nerd namespace */

#endif /* end of include guard: ADD_H */
