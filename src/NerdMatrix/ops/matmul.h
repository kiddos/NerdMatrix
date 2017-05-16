#ifndef MATMUL_H
#define MATMUL_H

#include "NerdMatrix/except/nerdmatrix_illegal_ops.h"
#include "NerdMatrix/matrix.h"
#include "NerdMatrix/wrapper/blas_wrapper.h"

namespace nerd {
namespace matrix {
namespace ops {

template <typename T>
void MatMul(const NerdMatrix<T>& m1, const NerdMatrix<T>& m2,
            NerdMatrix<T>& result) {
  if (m1.cols() != m2.rows()) {
    throw except::IllegalOpsException(
        "Matrix Multiplication: Input Size does not match.");
  } else {
  }
}

} /* end of ops namespace */
} /* end of matrix namespace */
} /* end of nerd namespace */

#endif /* end of include guard: MATMUL_H */
