#include <gtest/gtest.h>
#include <chrono>
#include <complex>
#include <random>
#include "NerdMatrix/except/nerdmatrix_exception.h"
#include "NerdMatrix/matrix.h"

using std::complex;
using nerd::matrix::NerdMatrix;

template <typename T>
class TestNerdMatrix : public ::testing::Test {
 public:
  void TestDefaultConstructor() {
    NerdMatrix<T> m;

    EXPECT_EQ(m.rows(), 0);
    EXPECT_EQ(m.cols(), 0);
    EXPECT_EQ(m.data(), nullptr);
  }

  void TestConstructor() {
    NerdMatrix<float> small(1, 1);
    EXPECT_EQ(small.rows(), 1);
    EXPECT_EQ(small.cols(), 1);
    EXPECT_NE(small.data(), nullptr);

    NerdMatrix<float> medium(1024, 1024);
    EXPECT_EQ(medium.rows(), 1024);
    EXPECT_EQ(medium.cols(), 1024);
    EXPECT_NE(medium.data(), nullptr);

    for (int i = 0; i < 1024 * 1024; ++i) {
      EXPECT_EQ(medium[i], 0);
    }

    // exception cases
    EXPECT_THROW(NerdMatrix<float> empty1(0, 1),
                 nerd::except::NerdMatrixException);
    EXPECT_THROW(NerdMatrix<float> empty2(1, 0),
                 nerd::except::NerdMatrixException);
    EXPECT_THROW(NerdMatrix<float> empty3(0, 0),
                 nerd::except::NerdMatrixException);
  }

  void TestCopyConstructor() {
    NerdMatrix<float> m1(1024, 1024);
    // randomize data
    std::mt19937 gen(
        std::chrono::system_clock::now().time_since_epoch().count());
    std::normal_distribution<float> dist(0, 1.0f);
    for (int i = 0; i < 1024 * 1024; ++i) {
      m1[i] = dist(gen);
    }

    NerdMatrix<float> m2 = m1;

    EXPECT_EQ(m1.rows(), m2.rows());
    EXPECT_EQ(m1.cols(), m2.cols());
    EXPECT_NE(m1.data(), m2.data());

    // check if data is copy
    for (int i = 0; i < 1024 * 1024; ++i) {
      EXPECT_EQ(m1[i], m2[i]);
    }
  }

  void TestMoveConstructor() {
    NerdMatrix<float> m(std::move(NerdMatrix<float>(1024, 1024)));

    EXPECT_EQ(m.rows(), 1024);
    EXPECT_EQ(m.cols(), 1024);
    EXPECT_NE(m.data(), nullptr);

    EXPECT_THROW(NerdMatrix<float> empty(std::move(NerdMatrix<float>(0, 0))),
                 nerd::except::NerdMatrixException);
  }

  void TestCopyOperator() {
    NerdMatrix<float> m1(1024, 1024);
    // randomize data
    std::mt19937 gen(
        std::chrono::system_clock::now().time_since_epoch().count());
    std::normal_distribution<float> dist(0, 1.0f);
    for (int i = 0; i < 1024 * 1024; ++i) {
      m1[i] = dist(gen);
    }

    NerdMatrix<float> m2(1024, 1024);
    m2 = m1;

    EXPECT_EQ(m1.rows(), m2.rows());
    EXPECT_EQ(m1.cols(), m2.cols());
    EXPECT_NE(m1.data(), m2.data());

    // check if data is copy
    for (int i = 0; i < 1024 * 1024; ++i) {
      EXPECT_EQ(m1[i], m2[i]);
    }
  }

  void TestMoveOperator() {
    NerdMatrix<float> m(1024, 1024);
    // randomize data
    std::mt19937 gen(
        std::chrono::system_clock::now().time_since_epoch().count());
    std::normal_distribution<float> dist(0, 1.0f);
    for (int i = 0; i < 1024 * 1024; ++i) {
      m[i] = dist(gen);
    }

    EXPECT_EQ(m.rows(), 1024);
    EXPECT_EQ(m.cols(), 1024);
    EXPECT_NE(m.data(), nullptr);

    m = std::move(NerdMatrix<float>(1024, 1024));

    // check if data is copy
    for (int i = 0; i < 1024 * 1024; ++i) {
      EXPECT_EQ(m[i], 0);
    }
  }
};

typedef ::testing::Types<float, double, complex<float>, complex<double>>
    TestTypes;
TYPED_TEST_CASE(TestNerdMatrix, TestTypes);

TYPED_TEST(TestNerdMatrix, DefaultConstructor) {
  this->TestDefaultConstructor();
}

TYPED_TEST(TestNerdMatrix, Constructor) { this->TestConstructor(); }

TYPED_TEST(TestNerdMatrix, CopyConstructor) { this->TestCopyConstructor(); }

TYPED_TEST(TestNerdMatrix, MoveConstructor) { this->TestMoveConstructor(); }

TYPED_TEST(TestNerdMatrix, CopyOperator) { this->TestCopyOperator(); }

TYPED_TEST(TestNerdMatrix, MoveOperator) { this->TestMoveOperator(); }
