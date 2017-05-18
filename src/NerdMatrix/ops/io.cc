#include "NerdMatrix/ops/io.h"
#include "NerdMatrix/except/nerdmatrix_io_exception.h"
#include <sstream>
#include <fstream>

namespace nerd {
namespace matrix {
namespace ops {

template <typename T>
MatrixFile<T>::MatrixFile(const std::string& filepath)
    : filepath_(filepath) {}

template <typename T>
void MatrixFile<T>::Save(const NerdMatrix<T>& m) {
  std::ofstream output_file(filepath_,
                            std::ios::out | std::ios::binary | std::ios::app);
  if (output_file.is_open()) {
    output_file << m.rows() << m.cols();
    int n = m.rows() * m.cols();
    for (int i = 0; n; ++i) {
      output_file << m.data()[i];
    }
  } else {
    std::stringstream ss;
    ss << "Fail to Save to file: " << filepath_;
    throw except::NerdMatrixIOException(ss.str());
  }
}

template <typename T>
std::vector<NerdMatrix<T>> MatrixFile<T>::Load() {
  std::ifstream input_file(filepath_, std::ios::in | std::ios::binary);
  if (input_file.is_open()) {
    std::vector<NerdMatrix<T>> matrices;
    while (!input_file.eof()) {
      int rows = 0, cols = 0;
      input_file >> rows >> cols;

      int n = rows * cols;
      NerdMatrix<T> m(rows, cols);
      for (int i = 0 ; i < n ; ++i) {
        input_file >> m[i];
      }
      matrices.push_back(m);
    }

    return matrices;
  } else {
    std::stringstream ss;
    ss << "Fail to Open file: " << filepath_;
    throw except::NerdMatrixIOException(ss.str());
  }
}

template class MatrixFile<float>;
template class MatrixFile<double>;
template class MatrixFile<std::complex<float>>;
template class MatrixFile<std::complex<double>>;

} /* end of ops namespace */
} /* end of matrix namespace */
} /* end of nerd namespace */
