#ifndef COMMON_H
#define COMMON_H

#include <iostream>
#include <stdint.h>
#include "../src/mat.h"

using core::mat;

void print_mat(mat<int> m);
void print_mat(mat<double> m);
void print_mat(mat<float> m);
void print_mat(mat<uint64_t> m);

#endif /* end of include guard: COMMON_H */
