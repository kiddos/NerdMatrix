#include <benchmark/benchmark.h>
#include "NerdMatrix/matrix.h"
#include "NerdMatrix/ops/matmul.h"

using std::complex;
using nerd::matrix::NerdMatrix;

template <typename T>
void BenchmarkMatMul(benchmark::State& state) {
  int size = state.range();
  NerdMatrix<T> a(size, size);
  NerdMatrix<T> b(size, size);
  NerdMatrix<T> c(size, size);
  while (state.KeepRunning()) {
    nerd::matrix::ops::MatMul(a, b, c);
  }
  state.SetBytesProcessed(int64_t(state.iterations()) *
                          int64_t(state.range(0)));
}

BENCHMARK_TEMPLATE(BenchmarkMatMul, float)
    ->Unit(benchmark::kMicrosecond)
    ->RangeMultiplier(2)
    ->Range(1 << 0, 1 << 10);

BENCHMARK_TEMPLATE(BenchmarkMatMul, double)
    ->Unit(benchmark::kMicrosecond)
    ->RangeMultiplier(2)
    ->Range(1 << 0, 1 << 10);

BENCHMARK_TEMPLATE(BenchmarkMatMul, complex<float>)
    ->Unit(benchmark::kMicrosecond)
    ->RangeMultiplier(2)
    ->Range(1 << 0, 1 << 10);

BENCHMARK_TEMPLATE(BenchmarkMatMul, complex<double>)
    ->Unit(benchmark::kMicrosecond)
    ->RangeMultiplier(2)
    ->Range(1 << 0, 1 << 10);

BENCHMARK_MAIN();
