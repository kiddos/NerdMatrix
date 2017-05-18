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
    data_ = new T[rows_ * cols_]{0};
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
    data_ = new T[rows_ * cols_]{0};
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

template <typename T>
NerdMatrix<T> NerdMatrix<T>::GetRow(int i) {
  if (i < 0 || i >= rows_) {
    throw except::NerdMatrixException("NerdMatrix GetRows: Invalid index.");
  }

  NerdMatrix<T> row(1, cols_);
  // wrapper::Copy(cols_, &data_[i * cols_], row.data_);
  for (int j = 0; j < cols_; ++j) {
    row.data_[j] = data_[i * cols_ + j];
  }
  return row;
}

template <typename T>
NerdMatrix<T> NerdMatrix<T>::GetCol(int j) {
  if (j < 0 || j >= cols_) {
    throw except::NerdMatrixException("NerdMatrix GetRows: Invalid index.");
  }

  NerdMatrix<T> col(cols_, 1);
  for (int i = 0; i < rows_; ++i) {
    col.data_[i] = data_[i * cols_ + j];
  }
  return col;
}

template <typename T>
NerdMatrix<T> NerdMatrix<T>::GetRows(int i1, int i2) {
  if (i1 < 0 || i2 > rows_) {
    throw except::NerdMatrixException(
        "NerdMatrix GetRows: Invalid starting index.");
  }
  if (i2 <= i1) {
    throw except::NerdMatrixException("NerdMatrix GetRows: Invalid size.");
  }

  int nrows = i2 - i1;
  NerdMatrix<T> rows(nrows, cols_);
  for (int i = 0; i < nrows; ++i) {
    for (int j = 0; j < cols_; ++j) {
      rows.data_[i * cols_ + j] = data_[(i1 + i) * cols_ + j];
    }
  }
  return rows;
}

template <typename T>
NerdMatrix<T> NerdMatrix<T>::GetCols(int j1, int j2) {
  if (j1 < 0 || j2 > cols_) {
    throw except::NerdMatrixException(
        "NerdMatrix GetRows: Invalid starting index.");
  }
  if (j2 <= j1) {
    throw except::NerdMatrixException("NerdMatrix GetRows: Invalid size.");
  }

  int ncols = j2 - j1;
  NerdMatrix<T> cols(rows_, ncols);
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < ncols; ++j) {
      cols.data_[i * ncols + j] = data_[i * cols_ + j + j1];
    }
  }
  return cols;
}

template <typename T>
NerdMatrix<T> NerdMatrix<T>::SubMatrix(int i1, int i2, int j1, int j2) {
  if (i1 < 0 || i2 > rows_) {
    throw except::NerdMatrixException(
        "NerdMatrix GetRows: Invalid starting index.");
  }
  if (i2 <= i1) {
    throw except::NerdMatrixException("NerdMatrix GetRows: Invalid size.");
  }

  if (j1 < 0 || j2 > cols_) {
    throw except::NerdMatrixException(
        "NerdMatrix GetRows: Invalid starting index.");
  }
  if (j2 <= j1) {
    throw except::NerdMatrixException("NerdMatrix GetRows: Invalid size.");
  }

  int nrows = i2 - i1;
  int ncols = j2 - j1;
  NerdMatrix<T> sub(nrows, ncols);
  for (int i = 0; i < nrows; ++i) {
    for (int j = 0; j < ncols; ++j) {
      sub[i * ncols + j] = data_[(i + i1) * cols_ + j + j1];
    }
  }
  return sub;
}

template class NerdMatrix<float>;
template class NerdMatrix<double>;
template class NerdMatrix<std::complex<float>>;
template class NerdMatrix<std::complex<double>>;

} /* end of matrix namespace */
} /* end of nerd namespace */
