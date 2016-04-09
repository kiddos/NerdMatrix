#include <iostream>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

#include "common.h"
#include "../src/mat.h"

#include <armadillo>

using std::cout;
using std::endl;
using core::mat;

bool dequal(double a, double b) {
  const double diff = 0.00000001;
  if ((a > b && (a - b) <= diff) || (a < b && (b - a) <= diff) || a == b)
    return true;
  return false;
}

int main(int argc, char *argv[]) {
  srand(time(NULL));
  const uint32_t maxRows = 1000, maxCols = 1000;
  const uint32_t nrows = rand() % maxRows + 1;
  const uint32_t ncols = rand() % maxCols + 1;

  arma::mat A = arma::randu<arma::mat>(nrows, ncols);
  arma::mat B = arma::randu<arma::mat>(nrows, ncols);
  arma::mat C = arma::randu<arma::mat>(ncols, nrows);

  double *dataA = new double[nrows * ncols];
  double *dataB = new double[nrows * ncols];
  double *dataC = new double[ncols * nrows];
  for (uint64_t i = 0 ; i < nrows ; i ++) {
    for (uint32_t j = 0 ; j < ncols ; j ++) {
      dataA[i * ncols + j] = A(i, j);
      dataB[i * ncols + j] = B(i, j);
    }
  }
  for (uint64_t i = 0 ; i < ncols ; i ++) {
    for (uint32_t j = 0 ; j < nrows ; j ++) {
      dataC[i * nrows + j] = C(i, j);
    }
  }

  cout << "matrix dimension: (" << nrows << ", " << ncols << ")" << endl;
  mat<double> testA(nrows, ncols, dataA);
  mat<double> testB(nrows, ncols, dataB);
  mat<double> testC(ncols, nrows, dataC);

  double start = 0;

  start = omp_get_wtime();
  arma::mat addResult = A + B;
  arma::mat subResult = A - B;
  arma::mat mulResult = A * C;
  cout << "arma done, time use: " << (omp_get_wtime() - start) << endl;

  start = omp_get_wtime();
  mat<double> addTest = testA.add(testB);
  mat<double> subTest = testA.subtract(testB);
  mat<double> mulTest = testA.multiply(testC);
  cout << "test done, time use: " << (omp_get_wtime() - start)  << endl;

  for (uint64_t i = 0 ; i < nrows ; i ++) {
    for (uint64_t j = 0 ; j < ncols ; j ++) {
      if (!dequal(addTest(i, j), addResult(i, j))) {
        cout << "wrong addition result" << endl;
        return -1;
      }
      if (!dequal(subTest(i, j), subResult(i, j))) {
        cout << "wrong subtraction result" << endl;
        return -1;
      }
    }
  }
  for (uint64_t i = 0 ; i < mulTest.nrows ; i ++) {
    for (uint64_t j = 0 ; j < mulTest.ncols ; j ++) {
      if (!dequal(mulTest(i, j), mulResult(i, j))) {
        cout << "wrong multiplication result" << endl;
        return -1;
      }
    }
  }

  return 0;
}
