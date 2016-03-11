#ifndef MAT_H
#define MAT_H

typedef unsigned int uint32_t;

template <class T>
class mat {
public:
  mat();
  mat(const mat &m);
  mat(uint32_t nrows, uint32_t ncols);
  mat(uint32_t nrows, uint32_t ncols, T *data);
  mat(uint32_t nrows, uint32_t ncols, T **data);
  virtual ~mat();
  mat add(mat m);
  mat subtract(mat m);
  mat multiply(mat m);
  mat sum(uint32_t dim);

  virtual void operator=(T d);
  virtual void operator=(const mat &m);
  mat operator+(const mat &m) const;
  mat operator-(const mat &m) const;
  mat operator*(const mat &m) const;

  mat operator+(T d) const;
  mat operator-(T d) const;
  mat operator*(T d) const;
  mat operator/(T d) const;

  T operator()(uint32_t i, uint32_t j) const;

  // static functions
  static mat zero(uint32_t nrows, uint32_t ncols);
  static mat ones(uint32_t nrows, uint32_t ncols);
  static mat eye(uint32_t dim);

  // operation
  void addrow(const mat &m);
  void addcol(const mat &m);
  void transpose();
  mat t();

  uint32_t nrows, ncols;
  T *data;
};

#define nullptr ((T *) 0)

template <class T> mat<T>::mat() {
  nrows = ncols = 0;
  data = nullptr;
}

template <class T> mat<T>::mat(const mat<T> &m) {
  nrows = m.nrows;
  ncols = m.ncols;
  data = new T[nrows * ncols];

  for (uint32_t i = 0 ; i < nrows ; i ++) {
    for (uint32_t j = 0 ; j < ncols ; j ++) {
      data[i * ncols + j] = m.data[i * ncols + j];
    }
  }
}

template <class T> mat<T>::mat(uint32_t nrows, uint32_t ncols) {
  this->nrows = nrows;
  this->ncols = ncols;
  this->data = new T[nrows * ncols];

  for (uint32_t i = 0 ; i < nrows ; i ++) {
    for (uint32_t j = 0 ; j < ncols ; j ++) {
      data[i * ncols + j] = 0;
    }
  }
}

template <class T> mat<T>::mat(uint32_t nrows, uint32_t ncols, T *data) {
  this->nrows = nrows;
  this->ncols = ncols;
  this->data = new T[nrows * ncols];

  for (uint32_t i = 0 ; i < nrows ; i ++) {
    for (uint32_t j = 0 ; j < ncols ; j ++) {
      this->data[i * ncols + j] = data[i * ncols + j];
    }
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
    delete data;
}

template <class T> mat<T> mat<T>::add(mat<T> m) {
  if (ncols != m.ncols || nrows != m.nrows) {
    return mat();
  } else {
    mat newMat(this->nrows, this->ncols);
    const uint32_t ncols = this->ncols, nrows = this->nrows;
    for (uint32_t i = 0 ; i < nrows ; i ++) {
      for (uint32_t j = 0 ; j < ncols ; j ++) {
        const uint32_t index = i * ncols + j;
        newMat.data[index] = this->data[index] + m.data[index];
      }
    }
    return newMat;
  }
}

template <class T> mat<T> mat<T>::subtract(mat<T> m) {
  if (this->ncols != m.ncols || this->nrows != m.nrows) {
    return mat();
  } else {
    mat newMat(this->nrows, this->ncols);
    const uint32_t ncols = this->ncols, nrows = this->nrows;
    for (uint32_t i = 0 ; i < nrows ; i ++) {
      for (uint32_t j = 0 ; j < ncols ; j ++) {
        const uint32_t index = i * ncols + j;
        newMat.data[index] = this->data[index] - m.data[index];
      }
    }
    return newMat;
  }
}

template <class T> mat<T> mat<T>::multiply(mat<T> m) {
  if (this->ncols != m.nrows) {
    return mat();
  } else {
    mat newMat(this->nrows, m.ncols);
    for (uint32_t i = 0 ; i < this->nrows ; i ++) {
      for (uint32_t j = 0 ; j < m.ncols ; j ++) {
        const uint32_t index = i * m.ncols + j;
        newMat.data[index] = 0;
        for (uint32_t k = 0 ; k < this->ncols ; k ++) {
          newMat.data[index] += this->data[i * this->ncols + k] *
              m.data[k * m.ncols + j];
        }
      }
    }
    return newMat;
  }
}

template <class T> mat<T> mat<T>::sum(uint32_t dim) {
  mat sum;
  switch (dim) {
    sum = mat(1, 1);
    case 0:
      for (uint32_t i = 0 ; i < this->nrows ; i ++) {
        for (uint32_t j = 0 ; j < this->ncols ; j ++) {
          sum.data[0] += data[i * this->ncols + j];
        }
      }
      return sum;
      break;
    case 1:
      sum = mat(1, this->ncols);
      for (uint32_t j = 0 ; j < this->ncols ; j ++) {
        for (uint32_t i = 0 ; i < this->nrows ; i ++) {
          sum.data[j] += data[i * this->ncols + j];
        }
      }
      return sum;
      break;
    case 2:
      sum = mat(this->nrows, 1);
      for (uint32_t i = 0 ; i < this->nrows ; i ++) {
        for (uint32_t j = 0 ; j < this->ncols ; j ++) {
          sum.data[i] += data[i * this->ncols + j];
        }
      }
      return sum;
      break;
  }
  return sum;
}

template <class T> void mat<T>::operator=(T d) {
  nrows = 1;
  ncols = 1;
  delete data;
  data = new T[1];
  data[0] = d;
}

template <class T> void mat<T>::operator=(const mat<T> &m) {
  nrows = m.nrows;
  ncols = m.ncols;
  delete data;
  data = new T[m.nrows * m.ncols];

  for (uint32_t i = 0 ; i < nrows ; i ++) {
    for (uint32_t j = 0 ; j < ncols ; j ++) {
      const uint32_t index = i * ncols + j;
      data[index] = m.data[index];
    }
  }
}

template <class T> mat<T> mat<T>::operator+(const mat<T> &m) const {
  return add(m);
}

template <class T> mat<T> mat<T>::operator-(const mat<T> &m) const {
  return subtract(m);
}

template <class T> mat<T> mat<T>::operator*(const mat<T> &m) const {
  return multiply(m);
}

template <class T> T mat<T>::operator()(uint32_t i, uint32_t j) const {
  return data[i * ncols + j];
}

template <class T> mat<T> mat<T>::zero(uint32_t nrows, uint32_t ncols) {
  return mat(nrows, ncols);
}

template <class T> mat<T> mat<T>::eye(uint32_t dim) {
  mat m(dim, dim);
  for (uint32_t i = 0 ; i < dim ; i ++) {
    m.data[i * m.ncols + i] = 1;
  }
  return m;
}

#undef nullptr

#endif /* end of include guard: MAT_H */
