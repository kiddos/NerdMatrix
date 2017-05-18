#include <gtest/gtest.h>
#include "NerdMatrix/except/nerdmatrix_illegal_ops.h"
#include "NerdMatrix/matrix.h"
#include "NerdMatrix/factory.h"
#include "NerdMatrix/ops/concat.h"

using std::complex;
using nerd::matrix::NerdMatrix;
using namespace nerd::matrix::factory;
using namespace nerd::matrix::ops;

template <typename T>
class TestNerdMatrixConcat : public ::testing::Test {
 public:
  int size_;
  void TestConcat() {
    NerdMatrix<T> a = Randn<T>(size_ + 1, size_ * 2, 0, 1.0);
    NerdMatrix<T> b = Randn<T>(size_ + 2, size_ * 2, 0, 1.0);
    NerdMatrix<T> c;

    Concat(a, b, c, 0);
    EXPECT_EQ(c.rows(), a.rows() + b.rows());
    EXPECT_EQ(c.cols(), a.cols());
    EXPECT_EQ(c.cols(), b.cols());

    for (int i = 0 ; i < a.rows() ; ++i) {
      for (int j = 0 ; j < c.cols() ; ++j) {
        EXPECT_EQ(c[i * c.cols() + j], a[i * a.cols() + j]);
      }
      int offset = a.rows();
      for (int j = 0 ; j < c.cols() ; ++j) {
        EXPECT_EQ(c[(i + offset) * c.cols() + j], b[i * b.cols() + j]);
      }
    }

    a = Randn<T>(size_ * 2, size_ + 1, 0, 1.0);
    b = Randn<T>(size_ * 2, size_ + 2, 0, 1.0);
    Concat(a, b, c, 1);
    EXPECT_EQ(c.rows(), a.rows());
    EXPECT_EQ(c.rows(), b.rows());
    EXPECT_EQ(c.cols(), a.cols() + b.cols());

    for (int i = 0 ; i < c.rows() ; ++i) {
      for (int j = 0 ; j < a.cols() ; ++j) {
        EXPECT_EQ(c[i * c.cols() + j], a[i * a.cols() + j]);
      }
      int offset = a.cols();
      for (int j = 0 ; j < b.cols() ; ++j) {
        EXPECT_EQ(c[i * c.cols() + j + offset], b[i * b.cols() + j]);
      }
    }

    a = Randn<T>(size_ + 1, size_ * 2 + 1, 0, 1.0);
    b = Randn<T>(size_ + 2, size_ * 2, 0, 1.0);
    EXPECT_THROW(Concat(a, b, c, 0), nerd::except::IllegalOpsException);

    a = Randn<T>(size_ * 2 + 1, size_ + 1, 0, 1.0);
    b = Randn<T>(size_ * 2, size_ + 2, 0, 1.0);
    EXPECT_THROW(Concat(a, b, c, 1), nerd::except::IllegalOpsException);
  }
};

typedef ::testing::Types<float, double, complex<float>, complex<double>>
  TestTypes;
TYPED_TEST_CASE(TestNerdMatrixConcat, TestTypes);

TYPED_TEST(TestNerdMatrixConcat, Concat) {
  this->size_ = 256;
  this->TestConcat();
}
