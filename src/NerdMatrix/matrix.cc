#include "NerdMatrix/matrix.h"
#include "NerdMatrix/except/nerdmatrix_exception.h"
#include "NerdMatrix/wrapper/blas_wrapper.h"

namespace nerd {
namespace matrix {

template <typename T>
NerdMatrix<T>::NerdMatrix() : rows_(0), cols_(0), data_(nullptr) {}

template <typename T>
NerdMatrix<T>::NerdMatrix(int rows, int cols)
    : rows_(rows), cols_(cols), data_(nullptr) {
  if (rows_ > 0 && cols_ > 0) {
    data_ = new T[rows * cols]{0};
    if (data_ == nullptr) {
      throw except::NerdMatrixException("System runs out of memory");
    }
  } else {
    throw except::NerdMatrixException(
        "Constructor With Invalid rows or column size");
  }
}

template <typename T>
NerdMatrix<T>::NerdMatrix(const NerdMatrix& m)
    : rows_(m.rows_), cols_(m.cols_), data_(nullptr) {
  if (rows_ > 0 && cols_ > 0) {
    data_ = new T[rows_ * cols_];
    wrapper::Copy(rows_ * cols_, m.data_, data_);
  }
}

template <typename T>
NerdMatrix<T>::NerdMatrix(NerdMatrix&& m)
    : rows_(m.rows_), cols_(m.cols_), data_(m.data_) {
  m.data_ = nullptr;
}

template <typename T>
NerdMatrix<T>& NerdMatrix<T>::operator=(const NerdMatrix& m) {
  int n = rows_ * cols_;
  rows_ = m.rows_;
  cols_ = m.cols_;
  if (n < m.rows_ * m.cols_) {
    delete[] data_;
    data_ = new T[rows_ * cols_];
    wrapper::Copy(rows_ * cols_, m.data_, data_);
  } else {
    wrapper::Copy(rows_ * cols_, m.data_, data_);
  }
  return *this;
}

template <typename T>
NerdMatrix<T>& NerdMatrix<T>::operator=(NerdMatrix&& m) {
  rows_ = m.rows_;
  cols_ = m.cols_;
  if (data_) delete[] data_;
  data_ = m.data_;
  m.data_ = nullptr;
  return *this;
}

template <typename T>
NerdMatrix<T>::~NerdMatrix() {
  if (data_) delete[] data_;
}

template class NerdMatrix<float>;
template class NerdMatrix<double>;
template class NerdMatrix<std::complex<float>>;
template class NerdMatrix<std::complex<double>>;

} /* end of matrix namespace */
} /* end of nerd namespace */
