#include <gtest/gtest.h>
#include "NerdMatrix/except/nerdmatrix_illegal_ops.h"
#include "NerdMatrix/matrix.h"
#include "NerdMatrix/ops/matmul.h"

using std::complex;
using nerd::matrix::NerdMatrix;

template <typename T>
class TestNerdMatrixMatMul : public ::testing::Test {
 public:
  int size_;

  void TestMatMul() {
    int m = size_;
    int n = size_;
    int k = size_;
    NerdMatrix<T> a(m, k);
    NerdMatrix<T> b(k, n);
    NerdMatrix<T> c(m, n);
    nerd::matrix::ops::MatMul(a, b, c);

    for (int i = 0; i < m; ++i) {
      for (int j = 0; j < n; ++j) {
        T value = T{0};
        for (int l = 0; l < k; ++l) {
          value += a[i * k + l] * b[l * n + j];
        }
        EXPECT_LT(std::abs(c[i * n + j] - value), 1e-4);
      }
    }

    b = NerdMatrix<T>(k + 1, n);
    EXPECT_THROW(nerd::matrix::ops::MatMul(a, b, c),
                 nerd::except::IllegalOpsException);
  }

  void TestMatMulOperator() {
    using namespace nerd::matrix::ops;

    int m = size_;
    int n = size_;
    int k = size_;
    NerdMatrix<T> a(m, k);
    NerdMatrix<T> b(k, n);
    NerdMatrix<T> c(m, n);

    c = a % b;
    for (int i = 0; i < m; ++i) {
      for (int j = 0; j < n; ++j) {
        T value = T{0};
        for (int l = 0; l < k; ++l) {
          value += a[i * k + l] * b[l * n + j];
        }
        EXPECT_LT(std::abs(c[i * n + j] - value), 1e-4);
      }
    }

    b = NerdMatrix<T>(k + 1, n);
    EXPECT_THROW(c = a % b, nerd::except::IllegalOpsException);
  }
};

typedef ::testing::Types<float, double, complex<float>, complex<double>>
    TestTypes;
TYPED_TEST_CASE(TestNerdMatrixMatMul, TestTypes);

TYPED_TEST(TestNerdMatrixMatMul, MatMul) {
  this->size_ = 256;
  this->TestMatMul();
}

TYPED_TEST(TestNerdMatrixMatMul, MatMulOperator) {
  this->size_ = 256;
  this->TestMatMulOperator();
}
