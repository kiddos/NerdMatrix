/*
 *  Test for transpose and inverse function
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

  const int size = 3;
  double data1[] = {-24, 18, 5, 20, -15, -4, -5, 4, 1};
  double data2[] = {1, 2, 3, 0, 1, 4, 5, 6, 0};

  mat<double> m1(size, size, data1);
  cout << m1;
  cout << m1.inv().inv().inv();
  mat<double> m2(size, size, data2);
  cout << m2;

  return 0;
}
