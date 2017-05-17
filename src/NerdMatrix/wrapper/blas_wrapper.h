#ifndef BLAS_WRAPPER_H
#define BLAS_WRAPPER_H

#include <complex>

namespace nerd {
namespace matrix {
namespace wrapper {

template <typename T>
void Copy(int n, const T* const src, T* dst);

template <typename T>
void Axpy(int n, T alpha, const T* const src, T* dst);

template <typename T>
void Gemm(int m, int n, int k, T alpha, const T* const a, const T* const b,
          T beta, T* dst);

} /* end of wrapper namespace */
} /* end of matrix namespace */
} /* end of nerd namespace */

#endif /* end of include guard: BLAS_WRAPPER_H */
