#include <gtest/gtest.h>
#include <chrono>
#include <random>
#include "NerdMatrix/except/nerdmatrix_exception.h"
#include "NerdMatrix/matrix.h"
#include "NerdMatrix/ops/add.h"
#include "NerdMatrix/ops/sub.h"
#include "NerdMatrix/ops/mul.h"
#include "NerdMatrix/ops/div.h"

using std::complex;
using nerd::matrix::NerdMatrix;
using namespace nerd::matrix::ops;

template <typename T>
void RandomInit(NerdMatrix<T>& m) {
  std::mt19937 gen(std::chrono::system_clock::now().time_since_epoch().count());
  std::normal_distribution<double> dist(0, 1.0f);
  for (int i = 0; i < m.rows() * m.cols(); ++i) {
    m[i] = static_cast<T>(dist(gen));
  }
}

template <typename T>
class TestNerdMatrixElementwiseOps : public ::testing::Test {
 public:
  void TestAdd() {
    int size = 1024;
    NerdMatrix<T> a(size, size);
    NerdMatrix<T> b(size, size);
    RandomInit(a);
    RandomInit(b);

    NerdMatrix<T> c = a + b;

    for (int i = 0 ; i < size ; ++i) {
      EXPECT_EQ(c[i], a[i] + b[i]);
    }

    NerdMatrix<T> d(size, size -1);
    EXPECT_THROW(c = c + d, nerd::except::IllegalOpsException);
  }

  void TestSub() {
    int size = 1024;
    NerdMatrix<T> a(size, size);
    NerdMatrix<T> b(size, size);
    RandomInit(a);
    RandomInit(b);

    NerdMatrix<T> c = a - b;

    for (int i = 0 ; i < size ; ++i) {
      EXPECT_EQ(c[i], a[i] - b[i]);
    }

    NerdMatrix<T> d(size, size -1);
    EXPECT_THROW(c = c - d, nerd::except::IllegalOpsException);
  }

  void TestMul() {
    int size = 1024;
    NerdMatrix<T> a(size, size);
    NerdMatrix<T> b(size, size);
    RandomInit(a);
    RandomInit(b);

    NerdMatrix<T> c = a * b;

    for (int i = 0 ; i < size ; ++i) {
      EXPECT_EQ(c[i], a[i] * b[i]);
    }

    NerdMatrix<T> d(size, size -1);
    EXPECT_THROW(c = c * d, nerd::except::IllegalOpsException);
  }

  void TestDiv() {
    int size = 1024;
    NerdMatrix<T> a(size, size);
    NerdMatrix<T> b(size, size);
    RandomInit(a);
    RandomInit(b);

    NerdMatrix<T> c = a / b;

    for (int i = 0 ; i < size ; ++i) {
      EXPECT_EQ(c[i], a[i] / b[i]);
    }

    NerdMatrix<T> d(size, size -1);
    EXPECT_THROW(c = c / d, nerd::except::IllegalOpsException);
  }
};

typedef ::testing::Types<float, double, complex<float>, complex<double>> TestTypes;
TYPED_TEST_CASE(TestNerdMatrixElementwiseOps, TestTypes);

TYPED_TEST(TestNerdMatrixElementwiseOps, Add) {
  this->TestAdd();
}

TYPED_TEST(TestNerdMatrixElementwiseOps, Sub) {
  this->TestSub();
}

TYPED_TEST(TestNerdMatrixElementwiseOps, Mul) {
  this->TestMul();
}

TYPED_TEST(TestNerdMatrixElementwiseOps, Div) {
  this->TestDiv();
}
