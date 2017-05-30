#ifndef BLAS_WRAPPER_H
#define BLAS_WRAPPER_H

namespace nerd {
namespace matrix {
namespace wrapper {

template <typename T>
void Copy(int n, const T* const src, T* dst);

template <typename T>
void Axpy(int n, T alpha, const T* const src, T* dst);

template <typename T>
void Swap(int n, T* x, T* y);

template <typename T>
void Scale(int n, T alpha, T* x);

template <typename T>
T Dot(int n, const T* x, const T* y);

template <typename T>
void Gemm(int m, int n, int k, T alpha, const T* const a, const T* const b,
          T beta, T* dst);

} /* end of wrapper namespace */
} /* end of matrix namespace */
} /* end of nerd namespace */

#endif /* end of include guard: BLAS_WRAPPER_H */
