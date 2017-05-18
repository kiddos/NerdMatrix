#ifndef IO_H
#define IO_H

#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include "NerdMatrix/matrix.h"

namespace nerd {
namespace matrix {
namespace ops {

template <typename T>
std::ostream& operator<<(std::ostream& os, const NerdMatrix<T>& m) {
  for (int i = 0; i < m.rows(); ++i) {
    for (int j = 0; j < m.cols(); ++j) {
      os << std::setw(10) << m[i * m.cols() + j];
    }
  }
  return os;
}

template <typename T>
class MatrixFile {
 public:
  MatrixFile(const std::string& filepath);

  void Save(const NerdMatrix<T>& m);
  std::vector<NerdMatrix<T>> Load();

 private:
  std::string filepath_;
};

} /* end of ops namespace */
} /* end of matrix namespace */
} /* end of nerd namespace */

#endif /* end of include guard: IO_H */
