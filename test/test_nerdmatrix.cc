#include <gtest/gtest.h>
#include <chrono>
#include <random>
#include "NerdMatrix/except/nerdmatrix_exception.h"
#include "NerdMatrix/matrix.h"

using std::complex;
using nerd::matrix::NerdMatrix;

template <typename T>
void RandomInitialize(NerdMatrix<T>& m) {
  std::mt19937 gen(std::chrono::system_clock::now().time_since_epoch().count());
  std::normal_distribution<float> dist(0, 1.0f);
  int n = m.rows() * m.cols();
  for (int i = 0; i < n; ++i) {
    m[i] = static_cast<T>(dist(gen));
  }
}

template <typename T>
class TestNerdMatrix : public ::testing::Test {
 public:
  int size_;

  void TestDefaultConstructor() {
    NerdMatrix<T> m;

    EXPECT_EQ(m.rows(), 0);
    EXPECT_EQ(m.cols(), 0);
    EXPECT_EQ(m.data(), nullptr);
  }

  void TestConstructor() {
    NerdMatrix<T> small(1, 1);
    EXPECT_EQ(small.rows(), 1);
    EXPECT_EQ(small.cols(), 1);
    EXPECT_NE(small.data(), nullptr);

    NerdMatrix<T> medium(size_, size_);
    EXPECT_EQ(medium.rows(), size_);
    EXPECT_EQ(medium.cols(), size_);
    EXPECT_NE(medium.data(), nullptr);

    for (int i = 0; i < size_ * size_; ++i) {
      EXPECT_EQ(medium[i], T{0});
    }

    // exception cases
    EXPECT_THROW(NerdMatrix<T> empty1(0, 1), nerd::except::NerdMatrixException);
    EXPECT_THROW(NerdMatrix<T> empty2(1, 0), nerd::except::NerdMatrixException);
    EXPECT_THROW(NerdMatrix<T> empty3(0, 0), nerd::except::NerdMatrixException);
  }

  void TestCopyConstructor() {
    NerdMatrix<T> m1(size_, size_);
    RandomInitialize(m1);

    NerdMatrix<T> m2 = m1;

    EXPECT_EQ(m1.rows(), m2.rows());
    EXPECT_EQ(m1.cols(), m2.cols());
    EXPECT_NE(m1.data(), m2.data());

    // check if data is copy
    for (int i = 0; i < size_ * size_; ++i) {
      EXPECT_EQ(m1[i], m2[i]);
    }
  }

  void TestMoveConstructor() {
    NerdMatrix<T> m(std::move(NerdMatrix<T>(size_, size_)));

    EXPECT_EQ(m.rows(), size_);
    EXPECT_EQ(m.cols(), size_);
    EXPECT_NE(m.data(), nullptr);

    EXPECT_THROW(NerdMatrix<T> empty(std::move(NerdMatrix<T>(0, 0))),
                 nerd::except::NerdMatrixException);
  }

  void TestCopyOperator() {
    NerdMatrix<T> m1(size_, size_);
    RandomInitialize(m1);

    NerdMatrix<T> m2(size_, size_);
    m2 = m1;
    EXPECT_EQ(m1.rows(), m2.rows());
    EXPECT_EQ(m1.cols(), m2.cols());
    EXPECT_NE(m1.data(), m2.data());

    // check if data is copy
    for (int i = 0; i < size_ * size_; ++i) {
      EXPECT_EQ(m1[i], m2[i]);
    }
  }

  void TestMoveOperator() {
    NerdMatrix<T> m(size_, size_);
    RandomInitialize(m);

    m = std::move(NerdMatrix<T>(size_, size_));
    EXPECT_EQ(m.rows(), size_);
    EXPECT_EQ(m.cols(), size_);
    EXPECT_NE(m.data(), nullptr);

    // check if data is copy
    for (int i = 0; i < size_ * size_; ++i) {
      EXPECT_EQ(m[i], T{0});
    }
  }

  void TestGetRow() {
    NerdMatrix<T> m(size_, size_);
    RandomInitialize(m);

    for (int i = 0; i < m.rows(); ++i) {
      NerdMatrix<T> r = m.GetRow(i);
      EXPECT_EQ(r.rows(), 1);
      EXPECT_EQ(r.cols(), m.cols());
      for (int j = 0; j < m.cols(); ++j) {
        EXPECT_EQ(r[j], m[i * m.cols() + j]);
      }
    }

    EXPECT_THROW(m.GetRow(-1), nerd::except::NerdMatrixException);
    EXPECT_THROW(m.GetRow(size_), nerd::except::NerdMatrixException);
  }

  void TestGetCol() {
    NerdMatrix<T> m(size_, size_);
    RandomInitialize(m);

    for (int j = 0; j < m.cols(); ++j) {
      NerdMatrix<T> c = m.GetCol(j);
      EXPECT_EQ(c.rows(), m.rows());
      EXPECT_EQ(c.cols(), 1);
      for (int i = 0; i < m.rows(); ++i) {
        EXPECT_EQ(c[i], m[i * m.cols() + j]);
      }
    }

    EXPECT_THROW(m.GetCol(-1), nerd::except::NerdMatrixException);
    EXPECT_THROW(m.GetCol(size_), nerd::except::NerdMatrixException);
  }

  void TestGetRows() {
    NerdMatrix<T> m(size_, size_);
    RandomInitialize(m);

    int rowsize = m.rows() / 2;
    for (int i = 0; i < m.rows() / 2; ++i) {
      NerdMatrix<T> r = m.GetRows(i, i + rowsize);
      EXPECT_EQ(r.rows(), rowsize);
      EXPECT_EQ(r.cols(), m.cols());
      for (int k = 0; k < rowsize; ++k) {
        for (int j = 0; j < m.cols(); ++j) {
          EXPECT_EQ(r[k * r.cols() + j], m[(k + i) * m.cols() + j]);
        }
      }
    }

    EXPECT_THROW(m.GetRows(-1, 0), nerd::except::NerdMatrixException);
    EXPECT_THROW(m.GetRows(0, 0), nerd::except::NerdMatrixException);
    EXPECT_THROW(m.GetRows(size_, size_), nerd::except::NerdMatrixException);
    EXPECT_THROW(m.GetRows(1, 1), nerd::except::NerdMatrixException);
    EXPECT_THROW(m.GetRows(size_, 1), nerd::except::NerdMatrixException);
    EXPECT_THROW(m.GetRows(0, size_ + 1), nerd::except::NerdMatrixException);
  }

  void TestGetCols() {
    NerdMatrix<T> m(size_, size_);
    RandomInitialize(m);

    int colsize = m.cols() / 2;
    for (int j = 0; j < m.cols() / 2; ++j) {
      NerdMatrix<T> r = m.GetCols(j, j + colsize);
      EXPECT_EQ(r.rows(), m.rows());
      EXPECT_EQ(r.cols(), colsize);
      for (int i = 0; i < m.rows(); ++i) {
        for (int k = 0; k < colsize; ++k) {
          EXPECT_EQ(r[i * r.cols() + k], m[i * m.cols() + j + k]);
        }
      }
    }

    EXPECT_THROW(m.GetCols(-1, 0), nerd::except::NerdMatrixException);
    EXPECT_THROW(m.GetCols(0, 0), nerd::except::NerdMatrixException);
    EXPECT_THROW(m.GetCols(size_, size_), nerd::except::NerdMatrixException);
    EXPECT_THROW(m.GetCols(1, 1), nerd::except::NerdMatrixException);
    EXPECT_THROW(m.GetCols(size_, 1), nerd::except::NerdMatrixException);
    EXPECT_THROW(m.GetCols(0, size_ + 1), nerd::except::NerdMatrixException);
  }

  void TestSubMatrix() {
    NerdMatrix<T> m(size_, size_);
    RandomInitialize(m);

    int rowsize = m.rows() / 2;
    int colsize = m.cols() / 2;
    for (int i = 0; i < m.rows() / 2; ++i) {
      for (int j = 0; j < m.cols() / 2; ++j) {
        NerdMatrix<T> sub = m.SubMatrix(i, i + rowsize, j, j + colsize);
        for (int k = 0; k < rowsize; ++k) {
          for (int l = 0; l < colsize; ++l) {
            EXPECT_EQ(sub[k * sub.cols() + l], m[(i + k) * m.cols() + j + l]);
          }
        }
      }
    }

    EXPECT_THROW(m.SubMatrix(-1, 0, 0, 0), nerd::except::NerdMatrixException);
    EXPECT_THROW(m.SubMatrix(0, 0, 0, 0), nerd::except::NerdMatrixException);
    EXPECT_THROW(m.SubMatrix(size_, size_, size_, size_),
                 nerd::except::NerdMatrixException);
  }
};

typedef ::testing::Types<float, double, complex<float>, complex<double>>
    TestTypes;
TYPED_TEST_CASE(TestNerdMatrix, TestTypes);

TYPED_TEST(TestNerdMatrix, DefaultConstructor) {
  this->TestDefaultConstructor();
}

TYPED_TEST(TestNerdMatrix, Constructor) {
  this->size_ = 1024;
  this->TestConstructor();
}

TYPED_TEST(TestNerdMatrix, CopyConstructor) {
  this->size_ = 1024;
  this->TestCopyConstructor();
}

TYPED_TEST(TestNerdMatrix, MoveConstructor) {
  this->size_ = 1024;
  this->TestMoveConstructor();
}

TYPED_TEST(TestNerdMatrix, CopyOperator) {
  this->size_ = 1024;
  this->TestCopyOperator();
}

TYPED_TEST(TestNerdMatrix, MoveOperator) {
  this->size_ = 1024;
  this->TestMoveOperator();
}

TYPED_TEST(TestNerdMatrix, GetRow) {
  this->size_ = 1024;
  this->TestGetRow();
}

TYPED_TEST(TestNerdMatrix, GetCol) {
  this->size_ = 1024;
  this->TestGetCol();
}

TYPED_TEST(TestNerdMatrix, GetRows) {
  this->size_ = 256;
  this->TestGetRows();
}

TYPED_TEST(TestNerdMatrix, GetCols) {
  this->size_ = 256;
  this->TestGetCols();
}

TYPED_TEST(TestNerdMatrix, SubMatrix) {
  this->size_ = 128;
  this->TestSubMatrix();
}
