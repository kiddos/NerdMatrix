#include <gtest/gtest.h>
#include "NerdMatrix/except/nerdmatrix_exception.h"
#include "NerdMatrix/factory.h"
#include "NerdMatrix/matrix.h"

using std::complex;
using nerd::matrix::NerdMatrix;

template <typename T>
class TestNerdMatrixFactory : public ::testing::Test {
 public:
  int size_;

  void TestZeros() {
    NerdMatrix<T> m = nerd::matrix::factory::Zeros<T>(size_, size_);
    int n = m.rows() * m.cols();
    for (int i = 0; i < n; ++i) {
      EXPECT_EQ(m[i], static_cast<T>(0));
    }

    EXPECT_THROW(nerd::matrix::factory::Zeros<T>(0, size_),
                 nerd::except::NerdMatrixException);
    EXPECT_THROW(nerd::matrix::factory::Zeros<T>(size_, 0),
                 nerd::except::NerdMatrixException);
    EXPECT_THROW(nerd::matrix::factory::Zeros<T>(-1, size_),
                 nerd::except::NerdMatrixException);
    EXPECT_THROW(nerd::matrix::factory::Zeros<T>(size_, -1),
                 nerd::except::NerdMatrixException);
  }

  void TestOnes() {
    NerdMatrix<T> m = nerd::matrix::factory::Ones<T>(size_, size_);
    int n = m.rows() * m.cols();
    for (int i = 0; i < n; ++i) {
      EXPECT_EQ(m[i], static_cast<T>(1));
    }

    EXPECT_THROW(nerd::matrix::factory::Ones<T>(0, size_),
                 nerd::except::NerdMatrixException);
    EXPECT_THROW(nerd::matrix::factory::Ones<T>(size_, 0),
                 nerd::except::NerdMatrixException);
    EXPECT_THROW(nerd::matrix::factory::Ones<T>(-1, size_),
                 nerd::except::NerdMatrixException);
    EXPECT_THROW(nerd::matrix::factory::Ones<T>(size_, -1),
                 nerd::except::NerdMatrixException);
  }

  void TestConst() {
    std::mt19937 gen(
        std::chrono::system_clock::now().time_since_epoch().count());
    std::normal_distribution<float> dist(0.0f, 1.0f);

    T val = static_cast<T>(dist(gen));
    NerdMatrix<T> m = nerd::matrix::factory::Const<T>(size_, size_, val);
    int n = m.rows() * m.cols();
    for (int i = 0; i < n; ++i) {
      EXPECT_EQ(m[i], val);
    }

    EXPECT_THROW(nerd::matrix::factory::Const<T>(0, size_, val),
                 nerd::except::NerdMatrixException);
    EXPECT_THROW(nerd::matrix::factory::Const<T>(size_, 0, val),
                 nerd::except::NerdMatrixException);
    EXPECT_THROW(nerd::matrix::factory::Const<T>(-1, size_, val),
                 nerd::except::NerdMatrixException);
    EXPECT_THROW(nerd::matrix::factory::Const<T>(size_, -1, val),
                 nerd::except::NerdMatrixException);
  }

  void TestIdentity() {
    NerdMatrix<T> m = nerd::matrix::factory::Identity<T>(size_);
    EXPECT_EQ(m.rows(), m.cols());

    for (int i = 0; i < m.rows(); ++i) {
      for (int j = 0; j < m.cols(); ++j) {
        if (i == j) {
          EXPECT_EQ(m[i * m.cols() + j], static_cast<T>(1));
        } else {
          EXPECT_EQ(m[i * m.cols() + j], static_cast<T>(0));
        }
      }
    }

    EXPECT_THROW(nerd::matrix::factory::Identity<T>(0),
                 nerd::except::NerdMatrixException);
    EXPECT_THROW(nerd::matrix::factory::Identity<T>(-1),
                 nerd::except::NerdMatrixException);
  }
};

typedef ::testing::Types<float, double, complex<float>, complex<double>>
    TestTypes;
TYPED_TEST_CASE(TestNerdMatrixFactory, TestTypes);

TYPED_TEST(TestNerdMatrixFactory, Zeros) {
  this->size_ = 4096;
  this->TestZeros();
}

TYPED_TEST(TestNerdMatrixFactory, Ones) {
  this->size_ = 4096;
  this->TestOnes();
}

TYPED_TEST(TestNerdMatrixFactory, Const) {
  this->size_ = 4096;
  this->TestConst();
}

TYPED_TEST(TestNerdMatrixFactory, Identity) {
  this->size_ = 4096;
  this->TestIdentity();
}
