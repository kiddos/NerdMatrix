#ifndef MATRIX_H
#define MATRIX_H

namespace nerd {
namespace matrix {

template <typename T>
class NerdMatrix {
 public:
  NerdMatrix();
  NerdMatrix(int rows, int cols);
  NerdMatrix(const NerdMatrix& NerdMatrix);
  NerdMatrix(NerdMatrix&& NerdMatrix);
  NerdMatrix& operator=(const NerdMatrix& NerdMatrix);
  NerdMatrix& operator=(NerdMatrix&& NerdMatrix);
  ~NerdMatrix();

  int rows() const { return rows_; }
  int cols() const { return cols_; }
  T* data() const { return data_; }
  T& operator[](int index) { return data_[index]; }

 private:
  int rows_, cols_;
  T* data_;
};

} /* end of matrix namespace */
} /* end of nerd namespace */

#endif /* end of include guard: MATRIX_H */