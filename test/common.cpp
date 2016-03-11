#include "common.h"

using std::cout;
using std::endl;

void print_mat(mat<int> m) {
  for (uint32_t i = 0 ; i < m.nrows ; i ++) {
    for (uint32_t j = 0 ; j < m.ncols ; j ++) {
      cout << m.data[i * m.ncols + j] << " ";
    }
    cout << endl;
  }
}

void print_mat(mat<double> m) {
  for (uint32_t i = 0 ; i < m.nrows ; i ++) {
    for (uint32_t j = 0 ; j < m.ncols ; j ++) {
      cout << m.data[i * m.ncols + j] << " ";
    }
    cout << endl;
  }
}

void print_mat(mat<float> m) {
  for (uint32_t i = 0 ; i < m.nrows ; i ++) {
    for (uint32_t j = 0 ; j < m.ncols ; j ++) {
      cout << m.data[i * m.ncols + j] << " ";
    }
    cout << endl;
  }
}

void print_mat(mat<uint64_t> m) {
  for (uint32_t i = 0 ; i < m.nrows ; i ++) {
    for (uint32_t j = 0 ; j < m.ncols ; j ++) {
      cout << m.data[i * m.ncols + j] << " ";
    }
    cout << endl;
  }
}

