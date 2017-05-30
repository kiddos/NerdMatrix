#include "NerdMatrix/wrapper/blas_wrapper.h"
#include <complex>
#include <cblas.h>

namespace nerd {
namespace matrix {
namespace wrapper {

typedef std::complex<float> complexf;
typedef std::complex<double> complexd;

template <>
void Copy(int n, const float* const src, float* dst) {
  cblas_scopy(n, src, 1, dst, 1);
}

template <>
void Copy(int n, const double* const src, double* dst) {
  cblas_dcopy(n, src, 1, dst, 1);
}

template <>
void Copy(int n, const complexf* const src,
                               complexf* dst) {
  cblas_ccopy(n, src, 1, dst, 1);
}

template <>
void Copy(int n, const std::complex<double>* const src,
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
void Axpy(int n, complexf alpha,
          const complexf* const src, complexf* dst) {
  cblas_caxpy(n, &alpha, src, 1, dst, 1);
}

template <>
void Axpy(int n, std::complex<double> alpha,
          const std::complex<double>* const src, std::complex<double>* dst) {
  cblas_zaxpy(n, &alpha, src, 1, dst, 1);
}

template <>
void Swap(int n, float* x, float* y) {
  cblas_sswap(n, x, 1, y, 1);
}

template <>
void Swap(int n, double* x, double* y) {
  cblas_dswap(n, x, 1, y, 1);
}

template <>
void Swap(int n, complexf* x, complexf* y) {
  cblas_cswap(n, x, 1, y, 1);
}

template <>
void Swap(int n, complexd* x, complexd* y) {
  cblas_zswap(n, x, 1, y, 1);
}

template <>
void Scale(int n, float alpha, float* x) {
  cblas_sscal(n, alpha, x, 1);
}

template <>
void Scale(int n, double alpha, double* x) {
  cblas_dscal(n, alpha, x, 1);
}

template <>
void Scale(int n, complexf alpha, complexf* x) {
  cblas_cscal(n, &alpha, x, 1);
}

template <>
void Scale(int n, complexd alpha, complexd* x) {
  cblas_zscal(n, &alpha, x, 1);
}

template <>
float Dot(int n, const float* x, const float* y) {
  return cblas_sdot(n, x, 1, y, 1);
}

template <>
double Dot(int n, const double* x, const double* y) {
  return cblas_ddot(n, x, 1, y, 1);
}

template <>
complexf Dot(int n, const complexf* x, const complexf* y) {
  complexf c;
  cblas_cdotc_sub(n, x, 1, y, 1, &c);
  return c;
}

template <>
complexd Dot(int n, const complexd* x, const complexd* y) {
  complexd c;
  cblas_zdotc_sub(n, x, 1, y, 1, &c);
  return c;
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
void Gemm(int m, int n, int k, complexf alpha,
          const complexf* const a,
          const complexf* const b, complexf beta,
          complexf* c) {
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
