/*
 *  Test for single row/col operations
 */

#include <iostream>
#include <time.h>
#include <stdlib.h>
#include "../src/mat.h"
#include "../src/mat_io.h"

using core::mat;
using std::cout;
using std::endl;

int main(void) {
  srand(time(NULL));

  double r = static_cast<double>(rand() % 1000) / 1000.0;
  double r2 = static_cast<double>(rand() % 1000) / 1000.0;
  const size_t size = 1000;

  mat<double> m1(1, size);
  m1.addrow(0, r);
  m1.subtractrow(0, r2);
  m1.multiplyrow(0, r);
  m1.dividerow(0, r);
  for (uint32_t i = 0 ; i < m1.ncols ; ++i) {
    if (m1.data[i] != (r-r2)*r/r) {
      return -1;
    }
  }

  mat<double> m2(size, 1);
  m2.addcol(0, r);
  m2.subtractcol(0, r2);
  m2.multiplycol(0, r);
  m2.dividecol(0, r);
  for (uint32_t i = 0 ; i < m2.ncols ; ++i) {
    if (m2.data[i] != (r-r2)*r/r) {
      return -1;
    }
  }

  return 0;
}
