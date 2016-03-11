#include <iostream>
#include <string>
#include <stdint.h>
#include "../src/mat.h"
#include <time.h>
#include <stdlib.h>

using std::cout;
using std::cin;
using std::endl;
using std::string;

typedef unsigned long long int LONG;

int main(int argc, char *argv[]) {
  const uint32_t maxRows = 1000, maxCols = 1000, maxVal = 1000;
  srand(time(NULL));

  const uint32_t nrows = rand() % maxRows + 1;
  const uint32_t ncols = rand() % maxCols + 1;

  int *intData = new int[nrows * ncols];
  for (uint64_t i = 0 ; i < nrows * ncols ; i ++) {
    intData[i] = rand() % maxVal;
  }

  mat<int> intTest(nrows, ncols, intData);

  float *floatData = new float[nrows * ncols];
  for (uint64_t i = 0 ; i < nrows * ncols ; i ++) {
    floatData[i] = rand() % maxVal * 1.0f;
  }
  mat<float> floatTest(nrows, ncols, floatData);

  double *doubleData = new double[nrows * ncols];
  for (uint64_t i = 0 ; i < nrows * ncols ; i ++) {
    doubleData[i] = rand() % maxVal * 1.0;
  }
  mat<double> doubleTest(nrows, ncols, doubleData);

  LONG *longData = new LONG[nrows * ncols];
  for (uint64_t i = 0 ; i < nrows * ncols ; i ++) {
    longData[i] = rand() % maxVal;
  }
  mat<LONG> longTest(nrows, ncols, longData);

  bool intTestSuccess = true, floatTestSuccess = true;
  bool doubleTestSuccess = true, longTestSuccess = true;
  for (uint32_t i = 0 ; i < nrows ; i ++) {
    for (uint32_t j = 0 ; j < ncols ; j ++) {
      intTestSuccess = intTest(i, j) == intData[i * ncols + j];
      floatTestSuccess = floatTest(i, j) == floatData[i * ncols + j];
      doubleTestSuccess = doubleTest(i, j) == doubleData[i * ncols + j];
      longTestSuccess = longTest(i, j) == longData[i * ncols + j];
    }
  }

  int returnCode = 0;
  if (intTestSuccess) {
    cout << "mat<int> test success" << endl;
  } else {
    cout << "mat<int> test failed" << endl;
    returnCode = -1;
  }
  if (floatTestSuccess) {
    cout << "mat<float> test success" << endl;
  } else {
    cout << "mat<float> test failed" << endl;
    returnCode = -1;
  }
  if (doubleTestSuccess) {
    cout << "mat<double> test success" << endl;
  } else {
    cout << "mat<double> test failed" << endl;
    returnCode = -1;
  }
  if (longTestSuccess) {
    cout << "mat<LONG> test success" << endl;
  } else {
    cout << "mat<LONG> test failed" << endl;
    returnCode = -1;
  }

  return returnCode;
}

