#ifndef CONCAT_H
#define CONCAT_H

#include "NerdMatrix/except/nerdmatrix_illegal_ops.h"
#include "NerdMatrix/matrix.h"
#include "NerdMatrix/wrapper/blas_wrapper.h"

namespace nerd {
namespace matrix {
namespace ops {

template <typename T>
void Concat(const NerdMatrix<T>& a, const NerdMatrix<T>& b, NerdMatrix<T>& c,
            int axis) {
  if (axis == 0) {
    if (a.cols() != b.cols()) {
      throw except::IllegalOpsException("Matrix Concat: Invalid size.");
    }

    int colsize = a.cols();
    int rowsize = a.rows() + b.rows();
    c = NerdMatrix<T>(rowsize, colsize);
    wrapper::Copy(a.rows() * a.cols(), a.data(), &c[0]);
    wrapper::Copy(b.rows() * b.cols(), b.data(), &c[a.rows() * a.cols()]);
  } else if (axis == 1) {
    if (a.rows() != b.rows()) {
      throw except::IllegalOpsException("Matrix Concat: Invalid size.");
    }

    int colsize = a.cols() + b.cols();
    int rowsize = a.rows();
    c = NerdMatrix<T>(rowsize, colsize);
    for (int i = 0; i < rowsize; ++i) {
      for (int j = 0; j < a.cols(); ++j) {
        c[i * colsize + j] = a.data()[i * a.cols() + j];
      }
      int offset = a.cols();
      for (int j = 0; j < b.cols(); ++j) {
        c[i * colsize + j + offset] = b.data()[i * b.cols() + j];
      }
    }
  } else {
    throw except::IllegalOpsException("Matrix Concat: Invalid axis.");
  }
}

} /* end of ops namespace */
} /* end of matrix namespace */
} /* end of nerd namespace */

#endif /* end of include guard: CONCAT_H */
