#include <benchmark/benchmark.h>
#include "NerdMatrix/matrix.h"
#include "NerdMatrix/ops/add.h"
#include "NerdMatrix/ops/sub.h"
#include "NerdMatrix/ops/mul.h"
#include "NerdMatrix/ops/div.h"

using nerd::matrix::NerdMatrix;

static void BenchmarkNerdMatrixAdd(benchmark::State& state) {
  using namespace nerd::matrix::ops;

  int size = state.range();
  NerdMatrix<float> a(size, size);
  NerdMatrix<float> b(size, size);
  NerdMatrix<float> c;
  while (state.KeepRunning()) {
    c = a + b;
  }
  state.SetBytesProcessed(int64_t(state.iterations()) *
              int64_t(state.range(0)));
}

BENCHMARK(BenchmarkNerdMatrixAdd)
  ->Unit(benchmark::kMicrosecond)
  ->RangeMultiplier(2)
  ->Range(8, 8 << 10);

static void BenchmarkNerdMatrixSub(benchmark::State& state) {
  using namespace nerd::matrix::ops;

  int size = state.range();
  NerdMatrix<float> a(size, size);
  NerdMatrix<float> b(size, size);
  NerdMatrix<float> c;
  while (state.KeepRunning()) {
    c = a - b;
  }
  state.SetBytesProcessed(int64_t(state.iterations()) *
              int64_t(state.range(0)));
}

BENCHMARK(BenchmarkNerdMatrixSub)
  ->Unit(benchmark::kMicrosecond)
  ->RangeMultiplier(2)
  ->Range(8, 8 << 10);

static void BenchmarkNerdMatrixMul(benchmark::State& state) {
  using namespace nerd::matrix::ops;

  int size = state.range();
  NerdMatrix<float> a(size, size);
  NerdMatrix<float> b(size, size);
  NerdMatrix<float> c;
  while (state.KeepRunning()) {
    c = a * b;
  }
  state.SetBytesProcessed(int64_t(state.iterations()) *
              int64_t(state.range(0)));
}

BENCHMARK(BenchmarkNerdMatrixMul)
  ->Unit(benchmark::kMicrosecond)
  ->RangeMultiplier(2)
  ->Range(8, 8 << 10);

static void BenchmarkNerdMatrixDiv(benchmark::State& state) {
  using namespace nerd::matrix::ops;

  int size = state.range();
  NerdMatrix<float> a(size, size);
  NerdMatrix<float> b(size, size);
  NerdMatrix<float> c;
  while (state.KeepRunning()) {
    c = a / b;
  }
  state.SetBytesProcessed(int64_t(state.iterations()) *
              int64_t(state.range(0)));
}

BENCHMARK(BenchmarkNerdMatrixDiv)
  ->Unit(benchmark::kMicrosecond)
  ->RangeMultiplier(2)
  ->Range(8, 8 << 10);

BENCHMARK_MAIN();
