#include <gtest/gtest.h>
#include <chrono>
#include <random>
#include "NerdMatrix/except/nerdmatrix_exception.h"
#include "NerdMatrix/matrix.h"
#include "NerdMatrix/ops/add.h"
#include "NerdMatrix/ops/sub.h"
#include "NerdMatrix/ops/mul.h"
#include "NerdMatrix/ops/div.h"

using nerd::matrix::NerdMatrix;

template <typename T>
void RandomInit(NerdMatrix<T>& m) {
  std::mt19937 gen(std::chrono::system_clock::now().time_since_epoch().count());
  std::normal_distribution<float> dist(0, 1.0f);
  for (int i = 0; i < m.rows() * m.cols(); ++i) {
    m[i] = dist(gen);
  }
}

TEST(TestNerdMatrixElementwiseOps, Add) {
  using namespace nerd::matrix::ops;

  int size = 1024;
  NerdMatrix<float> a(size, size);
  NerdMatrix<float> b(size, size);
  RandomInit(a);
  RandomInit(b);

  NerdMatrix<float> c = a + b;

  for (int i = 0 ; i < size ; ++i) {
    EXPECT_EQ(c[i], a[i] + b[i]);
  }

  NerdMatrix<float> d(size, size -1);
  EXPECT_THROW(c = c + d, nerd::except::IllegalOpsException);
}

TEST(TestNerdMatrixElementwiseOps, Sub) {
  using namespace nerd::matrix::ops;

  int size = 1024;
  NerdMatrix<float> a(size, size);
  NerdMatrix<float> b(size, size);
  RandomInit(a);
  RandomInit(b);

  NerdMatrix<float> c = a - b;

  for (int i = 0 ; i < size ; ++i) {
    EXPECT_EQ(c[i], a[i] - b[i]);
  }

  NerdMatrix<float> d(size, size -1);
  EXPECT_THROW(c = c - d, nerd::except::IllegalOpsException);
}

TEST(TestNerdMatrixElementwiseOps, Mul) {
  using namespace nerd::matrix::ops;

  int size = 1024;
  NerdMatrix<float> a(size, size);
  NerdMatrix<float> b(size, size);
  RandomInit(a);
  RandomInit(b);

  NerdMatrix<float> c = a * b;

  for (int i = 0 ; i < size ; ++i) {
    EXPECT_EQ(c[i], a[i] * b[i]);
  }

  NerdMatrix<float> d(size, size -1);
  EXPECT_THROW(c = c * d, nerd::except::IllegalOpsException);
}

TEST(TestNerdMatrixElementwiseOps, Div) {
  using namespace nerd::matrix::ops;

  int size = 1024;
  NerdMatrix<float> a(size, size);
  NerdMatrix<float> b(size, size);
  RandomInit(a);
  RandomInit(b);

  NerdMatrix<float> c = a / b;

  for (int i = 0 ; i < size ; ++i) {
    EXPECT_EQ(c[i], a[i] / b[i]);
  }

  NerdMatrix<float> d(size, size -1);
  EXPECT_THROW(c = c / d, nerd::except::IllegalOpsException);
}
