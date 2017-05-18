#ifndef MATRIX_H
#define MATRIX_H

#include <complex>

namespace nerd {
namespace matrix {

template <typename T>
class NerdMatrix {
 public:
  NerdMatrix();
  NerdMatrix(int rows, int cols);
  NerdMatrix(const NerdMatrix& m);
  NerdMatrix(NerdMatrix&& m);
  NerdMatrix& operator=(const NerdMatrix& m);
  NerdMatrix& operator=(NerdMatrix&& m);
  ~NerdMatrix();

  int rows() const { return rows_; }
  int cols() const { return cols_; }
  T* data() const { return data_; }
  T& operator[](int index) { return data_[index]; }

  NerdMatrix<T> GetRow(int i);
  NerdMatrix<T> GetCol(int j);
  NerdMatrix<T> GetRows(int i1, int i2);
  NerdMatrix<T> GetCols(int j1, int j2);
  NerdMatrix<T> SubMatrix(int i1, int i2, int j1, int j2);

 private:
  int rows_, cols_;
  T* data_;
};

typedef NerdMatrix<float> NerdMatrixf;
typedef NerdMatrix<double> NerdMatrixd;
typedef NerdMatrix<std::complex<float>> NerdMatrixc;
typedef NerdMatrix<std::complex<double>> NerdMatrixz;

} /* end of matrix namespace */
} /* end of nerd namespace */

#endif /* end of include guard: MATRIX_H */
