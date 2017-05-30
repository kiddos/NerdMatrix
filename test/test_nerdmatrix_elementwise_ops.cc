#include <gtest/gtest.h>
#include <chrono>
#include <random>
#include "NerdMatrix/except/nerdmatrix_exception.h"
#include "NerdMatrix/factory.h"
#include "NerdMatrix/matrix.h"
#include "NerdMatrix/ops/add.h"
#include "NerdMatrix/ops/div.h"
#include "NerdMatrix/ops/mul.h"
#include "NerdMatrix/ops/sub.h"

using std::complex;
using nerd::matrix::NerdMatrix;
using namespace nerd::matrix::factory;
using namespace nerd::matrix::ops;

template <typename T>
class TestNerdMatrixElementwiseOps : public ::testing::Test {
 public:
  int size_;

  void TestAdd() {
    NerdMatrix<T> a = Randn<T>(size_, size_, 0, 1);
    NerdMatrix<T> b = Randn<T>(size_, size_, 0, 1);
    NerdMatrix<T> c = a + b;
    for (int i = 0; i < size_; ++i) {
      EXPECT_EQ(c[i], a[i] + b[i]);
    }

    // function call test
    Add(a, b, c);
    for (int i = 0; i < size_; ++i) {
      EXPECT_EQ(c[i], a[i] + b[i]);
    }

    // single value
    c = a + 1;
    for (int i = 0; i < size_; ++i) {
      EXPECT_EQ(c[i], a[i] + static_cast<T>(1));
    }

    // add equal
    NerdMatrix<T> temp = c;
    c += a;
    for (int i = 0; i < size_; ++i) {
      EXPECT_EQ(c[i], temp[i] + a[i]);
    }

    NerdMatrix<T> d(size_, size_ - 1);
    EXPECT_THROW(c = c + d, nerd::except::IllegalOpsException);
  }

  void TestSub() {
    NerdMatrix<T> a = Randn<T>(size_, size_, 0, 1);
    NerdMatrix<T> b = Randn<T>(size_, size_, 0, 1);
    NerdMatrix<T> c = a - b;

    // function call test
    Sub(a, b, c);
    for (int i = 0; i < size_; ++i) {
      EXPECT_EQ(c[i], a[i] - b[i]);
    }

    // single value
    c = a - 1;
    for (int i = 0; i < size_; ++i) {
      EXPECT_EQ(c[i], a[i] - static_cast<T>(1));
    }

    // sub equal
    NerdMatrix<T> temp = c;
    c -= a;
    for (int i = 0; i < size_; ++i) {
      EXPECT_EQ(c[i], temp[i] - a[i]);
    }

    NerdMatrix<T> d(size_, size_ - 1);
    EXPECT_THROW(c = c - d, nerd::except::IllegalOpsException);
  }

  void TestMul() {
    NerdMatrix<T> a = Randn<T>(size_, size_, 0, 1);
    NerdMatrix<T> b = Randn<T>(size_, size_, 0, 1);
    NerdMatrix<T> c = a * b;

    // function call test
    Mul(a, b, c);
    for (int i = 0; i < size_; ++i) {
      EXPECT_EQ(c[i], a[i] * b[i]);
    }

    // single value
    c = a * 10;
    for (int i = 0; i < size_; ++i) {
      EXPECT_EQ(c[i], a[i] * static_cast<T>(10));
    }

    // mul equal
    NerdMatrix<T> temp = c;
    c *= a;
    for (int i = 0; i < size_; ++i) {
      EXPECT_EQ(c[i], temp[i] * a[i]);
    }

    NerdMatrix<T> d(size_, size_ - 1);
    EXPECT_THROW(c = c * d, nerd::except::IllegalOpsException);
  }

  void TestDiv() {
    NerdMatrix<T> a = Randn<T>(size_, size_, 0, 1);
    NerdMatrix<T> b = Randn<T>(size_, size_, 0, 1);
    NerdMatrix<T> c = a / b;

    // function call test
    Div(a, b, c);
    for (int i = 0; i < size_; ++i) {
      EXPECT_EQ(c[i], a[i] / b[i]);
    }

    // single value
    c = a / 10;
    for (int i = 0; i < size_; ++i) {
      EXPECT_EQ(c[i], a[i] / static_cast<T>(10));
    }

    // div equal
    NerdMatrix<T> temp = c;
    c /= a;
    for (int i = 0; i < size_; ++i) {
      EXPECT_EQ(c[i], temp[i] / a[i]);
    }

    NerdMatrix<T> d(size_, size_ - 1);
    EXPECT_THROW(c = c / d, nerd::except::IllegalOpsException);
  }
};

typedef ::testing::Types<float, double, complex<float>, complex<double>>
    TestTypes;
TYPED_TEST_CASE(TestNerdMatrixElementwiseOps, TestTypes);

TYPED_TEST(TestNerdMatrixElementwiseOps, Add) {
  this->size_ = 1024;
  this->TestAdd();
}

TYPED_TEST(TestNerdMatrixElementwiseOps, Sub) {
  this->size_ = 1024;
  this->TestSub();
}

TYPED_TEST(TestNerdMatrixElementwiseOps, Mul) {
  this->size_ = 1024;
  this->TestMul();
}

TYPED_TEST(TestNerdMatrixElementwiseOps, Div) {
  this->size_ = 1024;
  this->TestDiv();
}
