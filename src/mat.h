#ifndef MAT_H
#define MAT_H

#ifdef PARALLEL
#include <omp.h>
#endif

typedef unsigned int uint32_t;

namespace core {

#ifndef nullptr
#define nullptr ((T *) 0)
#endif

template <class T> class mat {
 public:
  mat();
  mat(const mat &m);
  mat(uint32_t nrows, uint32_t ncols);
  mat(uint32_t nrows, uint32_t ncols, T *data);
  mat(uint32_t nrows, uint32_t ncols, T **data);
  ~mat();
  void operator=(T d);
  mat& operator=(const mat &m);

  // fundamental operation
  mat add(const mat &m) const;
  mat subtract(const mat &m) const;
  mat multiply(const mat &m) const;
  T sum() const;
  mat sum(uint32_t dim) const;
  // row operation
  void addrow(uint32_t i, T v);
  void subtractrow(uint32_t i, T v);
  void multiplyrow(uint32_t i, T v);
  void dividerow(uint32_t i, T v);
  void addcol(uint32_t i, T v);
  void subtractcol(uint32_t i, T v);
  void multiplycol(uint32_t i, T v);
  void dividecol(uint32_t i, T v);
  mat appendrow(T v) const;
  mat insertrow(T v) const;
  mat appendcol(T v) const;
  mat insertcol(T v) const;
  // math operation
  mat f(T (*g)(T)) const;
  void func(T (*g)(T));
  // transpose and inverse
  mat t() const; // transpose
  mat inv() const; // inverse
  void transpose();
  void inverse();

  // operator
  mat operator+(const mat &m) const;
  mat operator-(const mat &m) const;
  mat operator*(const mat &m) const;
  mat operator%(const mat &m) const;

  mat operator-() const;
  mat operator+(T d) const;
  mat operator-(T d) const;
  mat operator*(T d) const;
  mat operator/(T d) const;

  T& operator()(uint32_t i, uint32_t j);
  T& operator[](uint32_t i);

  // static functions
  static mat zeros(uint32_t nrows, uint32_t ncols);
  static mat ones(uint32_t nrows, uint32_t ncols);
  static mat eyes(uint32_t dim);

  uint32_t nrows, ncols;
  T *data;
};

// implementation
template <class T> mat<T>::mat() {
  nrows = ncols = 0;
  data = nullptr;
}
template <class T> mat<T>::mat(const mat<T> &m) {
  nrows = m.nrows;
  ncols = m.ncols;
  data = new T[nrows * ncols];

  for (uint32_t i = 0 ; i < nrows * ncols ; i ++) {
    data[i] = m.data[i];
  }
}
template <class T> mat<T>::mat(uint32_t nrows, uint32_t ncols) {
  this->nrows = nrows;
  this->ncols = ncols;
  this->data = new T[nrows * ncols];

  for (uint32_t i = 0 ; i < nrows * ncols ; ++i) {
    data[i] = 0;
  }
}
template <class T> mat<T>::mat(uint32_t nrows, uint32_t ncols, T *data) {
  this->nrows = nrows;
  this->ncols = ncols;
  this->data = new T[nrows * ncols];

  for (uint32_t i = 0 ; i < nrows * ncols; i ++) {
    this->data[i] = data[i];
  }
}
template <class T> mat<T>::mat(uint32_t nrows, uint32_t ncols, T **data) {
  this->nrows = nrows;
  this->ncols = ncols;
  this->data = new T[nrows * ncols];

  for (uint32_t i = 0 ; i < nrows ; i ++) {
    for (uint32_t j = 0 ; j < ncols ; j ++) {
      this->data[i * ncols + j] = data[i][j];
    }
  }
}
template <class T> mat<T>::~mat() {
  if (data != nullptr)
    delete[] data;
}
// fundamental operation
template <class T> mat<T> mat<T>::add(const mat<T> &m) const {
  if (ncols != m.ncols || nrows != m.nrows) {
    return mat();
  } else {
    mat newmat(nrows, ncols);

    for (uint32_t i = 0 ; i < nrows * ncols ; ++i) {
      newmat.data[i] = data[i] + m.data[i];
    }
    return newmat;
  }
}
template <class T> mat<T> mat<T>::subtract(const mat<T> &m) const {
  if (this->ncols != m.ncols || this->nrows != m.nrows) {
    return mat();
  } else {
    mat newmat(nrows, ncols);

    for (uint32_t i = 0 ; i < nrows * ncols ; ++i) {
      newmat.data[i] = data[i] - m.data[i];
    }
    return newmat;
  }
}
template <class T> mat<T> mat<T>::multiply(const mat<T> &m) const {
  if (this->ncols != m.nrows) {
    return mat();
  } else {
    mat newmat(nrows, m.ncols);

    for (uint32_t i = 0 ; i < nrows * m.ncols ; ++i) {
      newmat.data[i] = 0;
      for (uint32_t k = 0 ; k < ncols ; ++k) {
        newmat.data[i] += data[(i / m.ncols) * ncols + k] *
            m.data[k * m.ncols + (i % m.ncols)];
      }
    }
    return newmat;
  }
}
template <class T> T mat<T>::sum() const {
  mat sum = mat(1, 1);

  for (uint32_t i = 0 ; i < this->nrows ; ++i) {
    for (uint32_t j = 0 ; j < this->ncols ; ++j) {
      sum.data[0] += data[i * this->ncols + j];
    }
  }
  return sum.data[0];
}
template <class T> mat<T> mat<T>::sum(uint32_t dim) const {
  mat sum;
  switch (dim) {
    case 0:
      sum = mat(1, 1);

      for (uint32_t i = 0 ; i < this->nrows ; ++i) {
        for (uint32_t j = 0 ; j < this->ncols ; ++j) {
          sum.data[0] += data[i * this->ncols + j];
        }
      }
      return sum;
      break;
    case 1:
      sum = mat(1, this->ncols);
      for (uint32_t j = 0 ; j < this->ncols ; ++j) {
        for (uint32_t i = 0 ; i < this->nrows ; ++i) {
          sum.data[j] += data[i * this->ncols + j];
        }
      }
      return sum;
      break;
    case 2:
      sum = mat(this->nrows, 1);
      for (uint32_t i = 0 ; i < this->nrows ; ++i) {
        for (uint32_t j = 0 ; j < this->ncols ; ++j) {
          sum.data[i] += data[i * this->ncols + j];
        }
      }
      return sum;
      break;
  }
  return sum;
}
// row operation on self
template <class T> void mat<T>::addrow(uint32_t i, T v) {
  const uint32_t rowindex = i * ncols;
  for (uint32_t k = 0 ; k < ncols ; ++k) {
    data[rowindex + k] += v;
  }
}
template <class T> void mat<T>::subtractrow(uint32_t i, T v) {
  const uint32_t rowindex = i * ncols;
  for (uint32_t k = 0 ; k < ncols ; ++k) {
    data[rowindex + k] -= v;
  }
}
template <class T> void mat<T>::multiplyrow(uint32_t i, T v) {
  const uint32_t rowindex = i * ncols;
  for (uint32_t k = 0 ; k < ncols ; ++k) {
    data[rowindex + k] *= v;
  }
}
template <class T> void mat<T>::dividerow(uint32_t i, T v) {
  const uint32_t rowindex = i * ncols;
  for (uint32_t k = 0 ; k < ncols ; ++k) {
    data[rowindex + k] /= v;
  }
}
template <class T> void mat<T>::addcol(uint32_t i, T v) {
  for (uint32_t k = 0 ; k < nrows ; ++k) {
    data[k * ncols + i] += v;
  }
}
template <class T> void mat<T>::subtractcol(uint32_t i, T v) {
  for (uint32_t k = 0 ; k < nrows ; ++k) {
    data[k * ncols + i] -= v;
  }
}
template <class T> void mat<T>::multiplycol(uint32_t i, T v) {
  for (uint32_t k = 0 ; k < nrows ; ++k) {
    data[k * ncols + i] *= v;
  }
}
template <class T> void mat<T>::dividecol(uint32_t i, T v) {
  for (uint32_t k = 0 ; k < nrows ; ++k) {
    data[k * ncols + i] /= v;
  }
}
// row operation
template <class T> mat<T> mat<T>::insertrow(T v) const {
  mat<T> newmat(this->nrows + 1, this->ncols);
  for (uint32_t i = 0 ; i < newmat.nrows ; ++i) {
    for (uint32_t j = 0 ; j < newmat.ncols ; ++j) {
      if (i == 0) {
        newmat.data[i * newmat.ncols + j] = v;
      } else {
        newmat.data[i * newmat.ncols + j] = this->data[(i-1) * this->ncols + j];
      }
    }
  }
  return newmat;
}
template <class T> mat<T> mat<T>::appendrow(T v) const {
  mat<T> newmat(this->nrows + 1, this->ncols);
  for (uint32_t i = 0 ; i < newmat.nrows ; ++i) {
    for (uint32_t j = 0 ; j < newmat.ncols ; ++j) {
      if (i == newmat.nrows - 1) {
        newmat.data[i * newmat.ncols + j] = v;
      } else {
        newmat.data[i * newmat.ncols + j] = this->data[(i-1) * this->ncols + j];
      }
    }
  }
  return newmat;
}
template <class T> mat<T> mat<T>::insertcol(T v) const {
  mat<T> newmat(this->nrows, this->ncols + 1);
  for (uint32_t i = 0 ; i < newmat.nrows ; ++ i) {
    newmat.data[i * newmat.ncols] = v;
    for (uint32_t j = 1 ; j < newmat.ncols ; ++ j) {
      newmat.data[i * newmat.ncols + j] = this->data[i * this->ncols + (j-1)];
    }
  }
  return newmat;
}
template <class T> mat<T> mat<T>::appendcol(T v) const {
  mat<T> newmat(this->nrows, this->ncols + 1);
  for (uint32_t i = 0 ; i < newmat.nrows ; ++ i) {
    for (uint32_t j = 0 ; j < newmat.ncols - 1 ; ++ j) {
      newmat.data[i * newmat.ncols + j] = this->data[i * this->ncols + (j-1)];
    }
    newmat.data[i * newmat.ncols + newmat.ncols - 1] = v;
  }
  return newmat;
}
// math operation
template <class T> mat<T> mat<T>::f(T (*g)(T)) const {
  mat<T> m(nrows, ncols);
  for (uint32_t i = 0 ; i < nrows ; ++i) {
    for (uint32_t j = 0 ; j < ncols ; ++j) {
      m.data[i * ncols + j] = g(data[i * ncols + j]);
    }
  }
  return m;
}
template <class T> void mat<T>::func(T (*g)(T)) {
  for (uint32_t i = 0 ; i < nrows ; ++i) {
    for (uint32_t j = 0 ; j < ncols ; ++j) {
      data[i * ncols + j] = g(data[i * ncols + j]);
    }
  }
}
// transpose and inverse
template <class T> mat<T> mat<T>::t() const {
  mat m(this->ncols, this->nrows);
  for (uint32_t i = 0 ; i < m.nrows ; ++i) {
    for (uint32_t j = 0 ; j < m.ncols ; ++j) {
      m.data[i * m.ncols + j] = this->data[j * this->ncols + i];
    }
  }
  return m;
}
template <class T> void mat<T>::transpose() {
  T *newdata = new T[nrows * ncols];
  for (uint32_t i = 0 ; i < nrows ; ++i) {
    for (uint32_t j = 0 ; j < ncols ; ++j) {
      newdata[j * nrows + i] = data[i * ncols + j];
    }
  }
  if (data) delete[] data;
  data = newdata;
}
// TODO need to fix inverse function rounding trailing
template <class T> mat<T> mat<T>::inv() const {
  if (nrows != ncols) return mat();

  mat m = eyes(nrows);
  mat temp(nrows, ncols, data);

  for (uint32_t k = 0 ; k < ncols ; ++k) {
    uint32_t nonzerorow = k;
    for (uint32_t i = k ; i < nrows ; ++i) {
      if (temp.data[i * ncols + k] != 0) {
        nonzerorow = i;
        m.dividerow(nonzerorow, temp.data[nonzerorow * ncols + k]);
        temp.dividerow(nonzerorow, temp.data[nonzerorow * ncols + k]);
        break;
      }
    }

    if (nonzerorow != k) {
      for (uint32_t j = 0 ; j < ncols ; ++j) {
        T t = temp.data[k * ncols + j];
        temp.data[k * ncols + j] = temp.data[nonzerorow * ncols + j];
        temp.data[nonzerorow * ncols + j] = t;

        t = m.data[k * ncols + j];
        m.data[k * ncols + j] = m.data[nonzerorow * ncols + j];
        m.data[nonzerorow * ncols + j] = t;
      }
    }

    for (uint32_t i = 0 ; i < nrows ; ++i) {
      if (i != k) {
        for (uint32_t j = 0 ; j < ncols ; ++j) {
          m.data[i * ncols + j] -= temp.data[i * ncols + k] * m.data[k * ncols + j];
          if (j != k) {
            temp.data[i * ncols + j] -= temp.data[i * ncols + k] *
                    temp.data[k * ncols + j];
          }
        }
        temp.data[i * ncols + k] = 0;
      }
    }
  }

  return m;
}
template <class T> void mat<T>::inverse() {
  if (nrows != ncols) return;

  // TODO check inversable
  mat m = eyes(nrows);

  for (uint32_t k = 0 ; k < ncols ; ++k) {
    uint32_t nonzerorow = k;
    for (uint32_t i = k ; i < nrows ; ++i) {
      if (data[i * ncols + k] != 0) {
        nonzerorow = i;
        m.dividerow(nonzerorow, data[nonzerorow * ncols + k]);
        dividerow(nonzerorow, data[nonzerorow * ncols + k]);
        break;
      }
    }

    if (nonzerorow != k) {
      for (uint32_t j = 0 ; j < ncols ; ++j) {
        T temp = data[k * ncols + j];
        data[k * ncols + j] = data[nonzerorow * ncols + j];
        data[nonzerorow * ncols + j] = temp;

        temp = m.data[k * ncols + j];
        m.data[k * ncols + j] = m.data[nonzerorow * ncols + j];
        m.data[nonzerorow * ncols + j] = temp;
      }
    }

    for (uint32_t i = 0 ; i < nrows ; ++i) {
      if (i != k) {
        for (uint32_t j = 0 ; j < ncols ; ++j) {
          m.data[i * ncols + j] -= data[i * ncols + k] * m.data[k * ncols + j];
          if (j != k) {
            data[i * ncols + j] -= data[i * ncols + k] * data[k * ncols + j];
          }
        }
        data[i * ncols + k] = 0;
      }
    }
  }

  // copy back
  for (uint32_t i = 0 ; i < nrows; ++ i) {
    for (uint32_t j = 0 ; j < ncols ; ++j) {
      data[i * ncols + j] = m.data[i * ncols + j];
    }
  }
}
// operators
template <class T> void mat<T>::operator=(T d) {
  // delete old data
  if (data) delete[] data;

  nrows = 1;
  ncols = 1;
  data = new T[1];
  data[0] = d;
}
template <class T> mat<T>& mat<T>::operator=(const mat<T> &m) {
  // delete old data
  if (data) delete[] data;

  nrows = m.nrows;
  ncols = m.ncols;
  data = new T[m.nrows * m.ncols];
  for (uint32_t i = 0 ; i < nrows * ncols ; ++i) {
    data[i] = m.data[i];
  }
  return *this;
}
template <class T> mat<T> mat<T>::operator+(const mat<T> &m) const {
  return this->add(m);
}
template <class T> mat<T> mat<T>::operator-(const mat<T> &m) const {
  return this->subtract(m);
}
template <class T> mat<T> mat<T>::operator*(const mat<T> &m) const {
  return this->multiply(m);
}
template <class T> mat<T> mat<T>::operator%(const mat<T> &m) const {
  if (this->ncols != m.ncols || this->nrows != m.nrows) {
    return mat();
  } else {
    mat newMat(this->nrows, this->ncols);

    for (uint32_t i = 0 ; i < nrows ; ++i) {
      for (uint32_t j = 0 ; j < ncols ; ++j) {
        const uint32_t index = i * ncols + j;
        newMat.data[index] = this->data[index] * m.data[index];
      }
    }
    return newMat;
  }
}
template <class T> mat<T> mat<T>::operator-() const {
  mat m(this->nrows, this->ncols);
  for (uint32_t i = 0 ; i < m.nrows ; ++ i) {
    for (uint32_t j = 0 ; j < m.ncols ; ++ j) {
      m.data[i * m.ncols + j] = - data[i * m.ncols + j];
    }
  }
  return m;
}
template <class T> mat<T> mat<T>::operator+(T d) const {
  mat m(nrows, ncols);
  for (uint32_t i = 0 ; i < m.nrows ; ++ i) {
    for (uint32_t j = 0 ; j < m.ncols ; ++ j) {
      m.data[i * m.ncols + j] = this->data[i * this->ncols + j] + d;
    }
  }
  return m;
}
template <class T> mat<T> mat<T>::operator-(T d) const {
  mat m(this->nrows, this->ncols);
  for (uint32_t i = 0 ; i < m.nrows ; ++ i) {
    for (uint32_t j = 0 ; j < m.ncols ; ++ j) {
      m.data[i * m.ncols + j] = this->data[i * this->ncols + j] - d;
    }
  }
  return m;
}
template <class T> mat<T> mat<T>::operator*(T d) const {
  mat m(this->nrows, this->ncols);
  for (uint32_t i = 0 ; i < m.nrows ; ++ i) {
    for (uint32_t j = 0 ; j < m.ncols ; ++ j) {
      m.data[i * m.ncols + j] = this->data[i * this->ncols + j] * d;
    }
  }
  return m;
}
template <class T> mat<T> mat<T>::operator/(T d) const {
  mat m(nrows, ncols);
  for (uint32_t i = 0 ; i < m.nrows ; ++ i) {
    for (uint32_t j = 0 ; j < m.ncols ; ++ j) {
      m.data[i * m.ncols + j] = this->data[i * this->ncols + j] / d;
    }
  }
  return m;
}
template <class T> T& mat<T>::operator()(uint32_t i, uint32_t j) {
  return data[i * ncols + j];
}
template <class T> T& mat<T>::operator[](uint32_t i) {
  return data[i];
}
// static functions
template <class T> mat<T> mat<T>::zeros(uint32_t nrows, uint32_t ncols) {
  return mat(nrows, ncols);
}
template <class T> mat<T> mat<T>::ones(uint32_t nrows, uint32_t ncols) {
  mat m(nrows, ncols);
  for (uint32_t i = 0 ; i < nrows ; ++ i) {
    for (uint32_t j = 0 ; j < ncols ; ++ j) {
      m.data[i * ncols + j] = 1;
    }
  }
  return m;
}
template <class T> mat<T> mat<T>::eyes(uint32_t dim) {
  mat m(dim, dim);
  for (uint32_t i = 0 ; i < dim ; i ++) {
    m.data[i * m.ncols + i] = 1;
  }
  return m;
}

typedef mat<double> Mat;
#undef nullptr
} // end of core namesapce

#endif /* end of include guard: MAT_H */
