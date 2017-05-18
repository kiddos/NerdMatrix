#include <benchmark/benchmark.h>
#include "NerdMatrix/factory.h"
#include "NerdMatrix/matrix.h"
#include "NerdMatrix/ops/matmul.h"

using std::complex;
using nerd::matrix::NerdMatrix;

template <typename T>
void BenchmarkNerdMatrix(benchmark::State& state) {
  int size = state.range();
  NerdMatrix<T> a;
  while (state.KeepRunning()) {
    a = NerdMatrix<T>(size, size);
  }
  state.SetBytesProcessed(int64_t(state.iterations()) *
                          int64_t(state.range(0)));
}

BENCHMARK_TEMPLATE(BenchmarkNerdMatrix, float)
    ->Unit(benchmark::kNanosecond)
    ->RangeMultiplier(4)
    ->Range(1 << 0, 2 << 10);

BENCHMARK_TEMPLATE(BenchmarkNerdMatrix, double)
    ->Unit(benchmark::kNanosecond)
    ->RangeMultiplier(4)
    ->Range(1 << 0, 2 << 10);

BENCHMARK_TEMPLATE(BenchmarkNerdMatrix, complex<float>)
    ->Unit(benchmark::kNanosecond)
    ->RangeMultiplier(4)
    ->Range(1 << 0, 2 << 10);

BENCHMARK_TEMPLATE(BenchmarkNerdMatrix, complex<double>)
    ->Unit(benchmark::kNanosecond)
    ->RangeMultiplier(4)
    ->Range(1 << 0, 2 << 10);

template <typename T>
void BenchmarkNerdMatrixSubMatrix(benchmark::State& state) {
  int size = state.range();
  NerdMatrix<T> a = nerd::matrix::factory::Randn<T>(size * 2, size * 2, 0, 1.0);
  while (state.KeepRunning()) {
    NerdMatrix<T> b =
        a.SubMatrix(size / 4, size * 5 / 4, size / 4, size * 5 / 4);
  }
  state.SetBytesProcessed(int64_t(state.iterations()) *
                          int64_t(state.range(0)));
}

BENCHMARK_TEMPLATE(BenchmarkNerdMatrixSubMatrix, float)
    ->Unit(benchmark::kMicrosecond)
    ->RangeMultiplier(4)
    ->Range(4 << 0, 1 << 10);

BENCHMARK_TEMPLATE(BenchmarkNerdMatrixSubMatrix, double)
    ->Unit(benchmark::kMicrosecond)
    ->RangeMultiplier(4)
    ->Range(4 << 0, 1 << 10);

BENCHMARK_TEMPLATE(BenchmarkNerdMatrixSubMatrix, complex<float>)
    ->Unit(benchmark::kMicrosecond)
    ->RangeMultiplier(4)
    ->Range(4 << 0, 1 << 10);

BENCHMARK_TEMPLATE(BenchmarkNerdMatrixSubMatrix, complex<double>)
    ->Unit(benchmark::kMicrosecond)
    ->RangeMultiplier(4)
    ->Range(4 << 0, 1 << 10);

BENCHMARK_MAIN();
