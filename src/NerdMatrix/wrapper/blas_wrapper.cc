#include "NerdMatrix/wrapper/blas_wrapper.h"
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
  cblas_zcopy(n, src, 1, dst, 1);
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
  cblas_zaxpy(n, &alpha, src, 1, dst, 1);
}

template <>
void Axpy(int n, std::complex<double> alpha,
          const std::complex<double>* const src, std::complex<double>* dst) {
  cblas_zaxpy(n, &alpha, src, 1, dst, 1);
}

} /* end of wrapper namespace */
} /* end of matrix namespace */
} /* end of nerd namespace */
