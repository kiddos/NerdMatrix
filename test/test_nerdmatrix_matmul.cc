#include <gtest/gtest.h>
#include <complex>
#include "NerdMatrix/except/nerdmatrix_illegal_ops.h"
#include "NerdMatrix/matrix.h"
#include "NerdMatrix/ops/matmul.h"

using std::complex;
using nerd::matrix::NerdMatrix;

template <typename T>
class TestNerdMatrixMatMul : public ::testing::Test {
 public:
  void TestMatMul() {
    int m = 1024;
    int n = 1024;
    int k = 1024;
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
};

typedef ::testing::Types<float, double, complex<float>, complex<double>>
    TestTypes;
TYPED_TEST_CASE(TestNerdMatrixMatMul, TestTypes);

TYPED_TEST(TestNerdMatrixMatMul, MatMul) { this->TestMatMul(); }
