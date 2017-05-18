#include "NerdMatrix/wrapper/blas_wrapper.h"
#include <complex>
#include <cblas.h>

namespace nerd {
namespace matrix {
namespace wrapper {

template <>
void Copy<float>(int n, const float* const src, float* dst) {
  cblas_scopy(n, src, 1, dst, 1);
}

template <>
void Copy<double>(int n, const double* const src, double* dst) {
  cblas_dcopy(n, src, 1, dst, 1);
}

template <>
void Copy<std::complex<float>>(int n, const std::complex<float>* const src,
                               std::complex<float>* dst) {
  cblas_ccopy(n, src, 1, dst, 1);
}

template <>
void Copy<std::complex<double>>(int n, const std::complex<double>* const src,
                                std::complex<double>* dst) {
  cblas_zcopy(n, src, 1, dst, 1);
}

template <>
void Axpy(int n, float alpha, const float* const src, float* dst) {
  cblas_saxpy(n, alpha, src, 1, dst, 1);
}

template <>
void Axpy(int n, double alpha, const double* const src, double* dst) {
  cblas_daxpy(n, alpha, src, 1, dst, 1);
}

template <>
void Axpy(int n, std::complex<float> alpha,
          const std::complex<float>* const src, std::complex<float>* dst) {
  cblas_caxpy(n, &alpha, src, 1, dst, 1);
}

template <>
void Axpy(int n, std::complex<double> alpha,
          const std::complex<double>* const src, std::complex<double>* dst) {
  cblas_zaxpy(n, &alpha, src, 1, dst, 1);
}

template <>
void Gemm(int m, int n, int k, float alpha, const float* const a,
          const float* const b, float beta, float* c) {
  cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, a, k,
              b, n, beta, c, n);
}

template <>
void Gemm(int m, int n, int k, double alpha, const double* const a,
          const double* const b, double beta, double* c) {
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, a, k,
              b, n, beta, c, n);
}

template <>
void Gemm(int m, int n, int k, std::complex<float> alpha,
          const std::complex<float>* const a,
          const std::complex<float>* const b, std::complex<float> beta,
          std::complex<float>* c) {
  cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, &alpha, a, k,
              b, n, &beta, c, n);
}

template <>
void Gemm(int m, int n, int k, std::complex<double> alpha,
          const std::complex<double>* const a,
          const std::complex<double>* const b, std::complex<double> beta,
          std::complex<double>* c) {
  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, &alpha, a, k,
              b, n, &beta, c, n);
}

} /* end of wrapper namespace */
} /* end of matrix namespace */
} /* end of nerd namespace */
