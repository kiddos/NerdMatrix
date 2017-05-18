#ifndef MATMUL_H
#define MATMUL_H

#include "NerdMatrix/except/nerdmatrix_illegal_ops.h"
#include "NerdMatrix/matrix.h"
#include "NerdMatrix/wrapper/blas_wrapper.h"

namespace nerd {
namespace matrix {
namespace ops {

template <typename T>
void MatMul(const NerdMatrix<T>& a, const NerdMatrix<T>& b, NerdMatrix<T>& c) {
  if (a.cols() != b.rows()) {
    throw except::IllegalOpsException(
        "Matrix Multiplication: Input Size does not match.");
  } else {
    int m = a.rows();
    int k = a.cols();
    int n = b.cols();
    c = NerdMatrix<T>(m, n);

    T alpha = static_cast<T>(1);
    T beta = static_cast<T>(0);
    wrapper::Gemm(m, n, k, alpha, a.data(), b.data(), beta, &c[0]);
  }
}

template <typename T>
NerdMatrix<T> operator%(const NerdMatrix<T>& a, const NerdMatrix<T>& b) {
  NerdMatrix<T> c;
  MatMul(a, b, c);
  return c;
}

} /* end of ops namespace */
} /* end of matrix namespace */
} /* end of nerd namespace */

#endif /* end of include guard: MATMUL_H */
