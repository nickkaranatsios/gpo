// Copyright 2007 Google Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// Author: Russell L. Smith

/*

This matrix:

          <--------->      <-- B_ blocks across, not counting edge blocks

        [ A##         # ]  <-- Each letter is a block.
        [ BC##        # ]      The ABC...STU blocks are N_*N_.
        [ DEF##       # ]      The VW...XY block are M_*N_.
        [  GHI##      . ]      The Z block is M_*M_.
        [     ...     . ]
        [      JKL##  . ]
        [       MNO## # ]
        [        PQR# # ]
        [         STU # ]
        [               ]
        [ VW.......XY Z ]

is stored as arrays of blocks like this:

        diag_data_ = [ ACFI ... LORU ]
        main_data_ = [ B DE GH ... JK MN PQ ST ]
        edge_data_ = [ VW ... XY ]
        last_data_ = [ Z ]

Only the lower triangle and diagonal of the diagonal blocks is stored, as the
upper triangle is symmetric. This saves a significant chunk of space.

The upper triangle '#' entries are not stored as they are symmetric to the
stored entries.

The number of scalars stored for B_,N_,M_ is

  scalars = diag_count  + main_count  + edge_count + last_count
          = B*N*(N+1)/2 + (2*B-3)*N*N + B*M*N      + M*(M+1)/2

For N=12 and M=6 (a common case) this is

  scalars = B*438 - 411

This means that 1 Gb of memory can hold B=306432, or a matrix of size 3,677,190
(double precision). At a rate of 100 Hz this is 51.1 minutes of data, where one
sample corresponds to one block. The rule of thumb is 300k blocks per Gb, or at
100 Hz about 50 minutes per Gb (this is a 3.7M x 3.7M matrix).

For N=9 and M=3 (another common case),

  scalars = B*234 - 237

Which means 1 Gb of memory holds B=573581, or 1:35 hours at 100 Hz sampling
rate.

*/

// Define this to use floating point exception handling in this code. Uncomment
// it to propagate infinities and NaNs without exceptions.
// NOTE: This is currently turned off as it seems to not work properly
// sometimes, e.g. when the wheel velocity sensors are used sometimes we get a
// divide by zero trap at perfectly innocent floating point load instructions.
//
//#define USE_FP_EXCEPTIONS

#include "gpo/common.h"
#include "gpo/hessian.h"
#include "gpo/fast_block.h"

using namespace GPO;

//***************************************************************************
// Utility

namespace {

// Return the smallest power of two that is larger than n.
int round_up_power_of_2(size_t n) {
  n |= n >> 1;
  n |= n >> 2;
  n |= n >> 4;
  n |= n >> 8;
  n |= n >> 16;
  // n |= n >> 32;         // If size_t is 64 bits
  return n+1;
}

}

//***************************************************************************
// Floating point error handling

#ifdef USE_FP_EXCEPTIONS

namespace {

volatile bool insolve = false;
jmp_buf setjmp_environment;
const char *fp_error_message = 0;

void fpe_signal_handler(int signum, siginfo_t *info, void *context) {
  DCHECK(insolve);
  switch(info->si_code) {
    case FPE_INTDIV:
      fp_error_message = "FP error: Integer divide by zero";
      break;
    case FPE_INTOVF:
      fp_error_message = "FP error: Integer overflow";
      break;
    case FPE_FLTDIV:
      fp_error_message = "FP error: Floating point divide by zero";
      break;
    case FPE_FLTOVF:
      fp_error_message = "FP error: Floating point overflow";
      break;
    case FPE_FLTUND:
      fp_error_message = "FP error: Floating point underflow";
      break;
    case FPE_FLTRES:
      fp_error_message = "FP error: Floating point inexact result";
      break;
    case FPE_FLTINV:
      fp_error_message = "FP error: Floating point invalid operation";
      break;
    case FPE_FLTSUB:
      fp_error_message = "FP error: Subscript out of range";
      break;
    default:
      fp_error_message = "FP error: unknown";
  }
  longjmp(setjmp_environment, 1);
}

void setup_fpe_signal_handling() {
  fp_error_message = 0;
  struct sigaction fpe_sigaction;
  memset(&fpe_sigaction, 0, sizeof(fpe_sigaction));
  fpe_sigaction.sa_sigaction = &fpe_signal_handler;
  fpe_sigaction.sa_flags = SA_RESETHAND | SA_SIGINFO;
  sigaction(SIGFPE, &fpe_sigaction, 0);
}

}

#endif

//***************************************************************************
// HessianMatrix reference implementation - slow but sure!

HessianMatrix::HessianMatrix() {
  Reset();
}

HessianMatrix::~HessianMatrix() {
#ifdef __GPO_HESSIAN_CACHE__
  // For debugging, print cache statistics:
  // printf ("Hessian diag cache size = %u, group size = %u, miss rate = %f%%\n",
  //         diag_cache_.GetCacheSize(), diag_cache_.GetGroupSize(), 100.0 * diag_cache_.GetMissRate());
  // printf ("Hessian main cache size = %u, group size = %u, miss rate = %f%%\n",
  //         main_cache_.GetCacheSize(), main_cache_.GetGroupSize(), 100.0 * main_cache_.GetMissRate());
  // printf ("Hessian edge cache size = %u, group size = %u, miss rate = %f%%\n",
  //         edge_cache_.GetCacheSize(), edge_cache_.GetGroupSize(), 100.0 * edge_cache_.GetMissRate());
#endif

  Free();
}

bool HessianMatrix::Initialize(int B, int N, int M, size_t max_memory) {
  DCHECK(B >= 2);              // This is assumed in many placed below
  Free();
  Reset();

  B_ = B;
  N_ = N;
  M_ = M;
  NN_ = N_*N_;
  NM_ = N_*M_;
  diag_size_ = (N_*(N_+1))/2;

#ifdef __GPO_HESSIAN_CACHE__
  if (max_memory == 0) max_memory = (size_t)-1;

  // Control the cache sizes so that we don't use more than 'max_memory' bytes.
  size_t blocksize_diag = diag_size_*sizeof(double);
  size_t blocksize_main = N_*N_*sizeof(double);
  size_t blocksize_edge = M_*N_*sizeof(double);

  size_t csize_diag = round_up_power_of_2(B_);
  size_t csize_main = round_up_power_of_2(2*B_-3);
  size_t csize_edge = round_up_power_of_2(B_);

  int reduction = 0;
  for (;;) {
    uint64 memory =
      (uint64)blocksize_diag * (uint64)csize_diag +
      (uint64)blocksize_main * (uint64)csize_main +
      (uint64)blocksize_edge * (uint64)csize_edge;
    if (memory <= max_memory) break;
    csize_diag /= 2;
    csize_main /= 2;
    csize_edge /= 2;
    reduction++;
  }
  if (csize_diag < 256) csize_diag = 256;
  if (csize_main < 256) csize_main = 256;
  if (csize_edge < 256) csize_edge = 256;
  diag_cache_.Initialize(blocksize_diag, csize_diag, "DIAG");
  main_cache_.Initialize(blocksize_main, csize_main, "MAIN");
  edge_cache_.Initialize(blocksize_edge, csize_edge, "EDGE");

  // For debugging:
  // printf ("Diag cache size = %d blocks (%d Mb)\n",
  //   csize_diag, (csize_diag*blocksize_diag)/(1024*1024));
  // printf ("Main cache size = %d blocks (%d Mb)\n",
  //   csize_main, (csize_main*blocksize_main)/(1024*1024));
  // printf ("Edge cache size = %d blocks (%d Mb)\n",
  //   csize_edge, (csize_edge*blocksize_edge)/(1024*1024));
  // printf ("Cache size reduction = 2^%d\n", reduction);
#else
  diag_count_ = B_*(N_*(N_+1))/2;
  main_count_ = (2*B_-3)*N_*N_;
  edge_count_ = B_*M_*N_;
  diag_data_ = new double[diag_count_];
  main_data_ = new double[main_count_];
  edge_data_ = new double[edge_count_];

  if (!diag_data_ || !main_data_ || !edge_data_) {
    Free();
    Reset();
    return false;
  }
#endif

  last_count_ = (M_*(M_+1))/2;
  last_data_ = new double[last_count_];
  if (!last_data_) {
    Free();
    Reset();
    return false;
  }

  SetZero();
  return true;
}

void HessianMatrix::SetZero() {
#ifdef __GPO_HESSIAN_CACHE__
  diag_cache_.SetZero();
  main_cache_.SetZero();
  edge_cache_.SetZero();
#else
  memset(diag_data_, 0, diag_count_*sizeof(double));
  memset(main_data_, 0, main_count_*sizeof(double));
  memset(edge_data_, 0, edge_count_*sizeof(double));
#endif
  memset(last_data_, 0, last_count_*sizeof(double));
}

void HessianMatrix::Add(int i, int j, double x) {
  double *element = Access(i, j);
  DCHECK(element);
  *element += x;
}

double * HessianMatrix::Access(int i, int j) {
  //printf ("%d,%d (B=%d N=%d M=%d)\n",i,j,B_,N_,M_);   // For debug

  // Ensure j <= i (swap if necessary)
  if (j > i) {
    int tmp = i;
    i = j;
    j = tmp;
  }

  int sz = B_*N_ + M_;
  if (i < 0 || i >= sz || j < 0 || j > i) return 0;

  if (i >= B_*N_) {
    // In the edge or the last block
    if (j >= B_*N_) {
      // In the last block
      i -= B_*N_;                       // Index within block
      j -= B_*N_;                       // Index within block
      int offset = (i*(i+1))/2 + j;
      DCHECK(offset >= 0 && offset < (M_*(M_+1))/2);
      return last_data_ + offset;
    } else {
      // In an edge block
      int bj = j / N_;                  // Edge block number
      i -= B_*N_;                       // Index within edge block
      j -= bj*N_;                       // Index within edge block
      #ifdef __GPO_HESSIAN_CACHE__
        double *data = (double*) edge_cache_.Access(bj);
        return data + i*N_ + j;
      #else
        int offset = bj*N_*M_ + i*N_ + j;
        DCHECK(offset >= 0 && offset < B_*M_*N_);
        return edge_data_ + offset;
      #endif
    }
  } else {
    // Not in the edge
    int bi = i / N_;                    // Row block
    int bj = j / N_;                    // Column block
    i -= bi*N_;                         // Index within block
    j -= bj*N_;                         // Index within block
    if (bi == bj) {
      // In a diagonal block
      #ifdef __GPO_HESSIAN_CACHE__
        double *data = (double*) diag_cache_.Access(bi);
        return data + (i*(i+1))/2 + j;
      #else
        int offset = bi*(N_*(N_+1))/2 + (i*(i+1))/2 + j;
        DCHECK(offset >= 0 && offset < B_*(N_*(N_+1))/2);
        return diag_data_ + offset;
      #endif
    } else {
      // In an off-diagonal block
      if (bj < bi-2) return 0;          // Sparsity constraint
      #ifdef __GPO_HESSIAN_CACHE__
        double *data = (double*) main_cache_.Access((2*bi-3) + (bj-bi+2));
        return data + i*N_ + j;
      #else
        int offset = ((2*bi-3) + (bj-bi+2)) * (N_*N_) + i*N_ + j;
        DCHECK(offset >= 0 && offset < (2*B_-3)*N_*N_);
        return main_data_ + offset;
      #endif
    }
  }
}

void HessianMatrix::Multiply(double *x, double *y) {
  // This is slow but simple, it should only be used for testing!
  int sz = B_*N_ + M_;
  for (int i = 0; i < sz; i++) {
    double sum = 0;
    for (int j = 0; j < sz; j++) {
      double *element = Access(i, j);
      if (element) sum += (*element) * x[j];
    }
    y[i] = sum;
  }
}

const char *HessianMatrix::Solve(double *vector) {
#ifdef USE_FP_EXCEPTIONS
  // Install a signal handler for floating point error that returns an error
  // message from this function.
  insolve = true;
  setup_fpe_signal_handling();
  if (setjmp(setjmp_environment)) {
    insolve = false;
    return fp_error_message;
  }
  //@@@TODO: disable this after we're done?
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif

  // Readability-improving macros to access the matrix blocks
#ifdef __GPO_HESSIAN_CACHE__
  #define DIAG(i) ((double*) diag_cache_.Access(i))
  #define MAIN(i) ((double*) main_cache_.Access(i))
  #define EDGE(i) ((double*) edge_cache_.Access(i))
#else
  #define DIAG(i) (diag_data_ + diag_size_*(i) )
  #define MAIN(i) (main_data_ + NN_*(i) )
  #define EDGE(i) (edge_data_ + NM_*(i) )
#endif
  #define LAST (last_data_)
  #define VECTOR(i) (vector + N_*(i))

  // Simultaneously factorize the interior, factorize the border, and perform
  // forward substitution on the right hand side vector.

  // Factor the first two block rows
  factor(DIAG(0), N_);
  solve(DIAG(0), EDGE(0), N_, M_);
  solve(DIAG(0), MAIN(0), N_, N_);
  multiply2_lt(DIAG(1), -1, MAIN(0), MAIN(0), N_, N_);
  factor(DIAG(1), N_);
  multiply2(EDGE(1), -1, EDGE(0), MAIN(0), M_, N_, N_);
  solve(DIAG(1), EDGE(1), N_, M_);
  multiply2_lt(LAST, -1, EDGE(0), EDGE(0), M_, N_);
  multiply2_lt(LAST, -1, EDGE(1), EDGE(1), M_, N_);

  // Forward substitution on the first two block rows
  solve(DIAG(0), VECTOR(0), N_, 1);
  multiply2(VECTOR(1), -1, MAIN(0), VECTOR(0), N_, N_, 1);
  solve(DIAG(1), VECTOR(1), N_, 1);
  multiply2(VECTOR(B_), -1, EDGE(0), VECTOR(0), M_, N_, 1);
  multiply2(VECTOR(B_), -1, EDGE(1), VECTOR(1), M_, N_, 1);

  for (int i = 2; i < B_; i++) {
    // Factor a block row
    solve(DIAG(i-2), MAIN(2*i-3), N_, N_);
    multiply2(MAIN(2*i-2), -1, MAIN(2*i-3), MAIN(2*i-4), N_, N_, N_);
    solve(DIAG(i-1), MAIN(2*i-2), N_, N_);
    multiply2_lt(DIAG(i), -1, MAIN(2*i-3), MAIN(2*i-3), N_, N_);
    multiply2_lt(DIAG(i), -1, MAIN(2*i-2), MAIN(2*i-2), N_, N_);
    factor(DIAG(i), N_);
    multiply2(EDGE(i), -1, EDGE(i-2), MAIN(2*i-3), M_, N_, N_);
    multiply2(EDGE(i), -1, EDGE(i-1), MAIN(2*i-2), M_, N_, N_);
    solve(DIAG(i), EDGE(i), N_, M_);
    multiply2_lt(LAST, -1, EDGE(i), EDGE(i), M_, N_);

    // Forward substitution on a block row
    multiply2(VECTOR(i), -1, MAIN(2*i-3), VECTOR(i-2), N_, N_, 1);
    multiply2(VECTOR(i), -1, MAIN(2*i-2), VECTOR(i-1), N_, N_, 1);
    solve(DIAG(i), VECTOR(i), N_, 1);
    multiply2(VECTOR(B_), -1, EDGE(i), VECTOR(i), M_, N_, 1);

#ifndef USE_FP_EXCEPTIONS
    // If we have FP exceptions turned off, try to detect problems by looking
    // at the factored diagonal element and solution vector, as infinities and
    // NaNs will tend to get propagated there.
    for (int j = 0; j < N_; j++) {
      double value1 = DIAG(i)[(j*(j+3))/2];
      double value2 = VECTOR(i)[j];
      double value3 = VECTOR(B_)[j];
      if (isinf(value1)) return "FP error: Infinity encountered (diag)";
      if (isinf(value2)) return "FP error: Infinity encountered (rhs vector)";
      if (isinf(value3)) return "FP error: Infinity encountered (rhs vector end)";
      if (isnan(value1)) return "FP error: NaN encountered (diag)";
      if (isnan(value2)) return "FP error: NaN encountered (rhs vector)";
      if (isnan(value3)) return "FP error: NaN encountered (rhs vector end)";
    }
#endif
  }

  // Factor the last block
  factor(LAST, M_);

  // Forward substitution on the last block
  solve(LAST, VECTOR(B_), M_, 1);

  // Backward substitution on all block rows
  solveT(LAST, VECTOR(B_), M_);
  multiply1(VECTOR(B_-1), -1, EDGE(B_-1), VECTOR(B_), N_, M_, 1);
  solveT(DIAG(B_-1), VECTOR(B_-1), N_);
  multiply1(VECTOR(B_-2), -1, MAIN(2*B_-4), VECTOR(B_-1), N_, N_, 1);
  multiply1(VECTOR(B_-2), -1, EDGE(B_-2), VECTOR(B_), N_, M_, 1);
  solveT(DIAG(B_-2), VECTOR(B_-2), N_);
  for (int i = B_-3; i >= 0; i--) {
    multiply1(VECTOR(i), -1, MAIN(2*i), VECTOR(i+1), N_, N_, 1);
    multiply1(VECTOR(i), -1, MAIN(2*i+1), VECTOR(i+2), N_, N_, 1);
    multiply1(VECTOR(i), -1, EDGE(i), VECTOR(B_), N_, M_, 1);
    solveT(DIAG(i), VECTOR(i), N_);

#ifndef USE_FP_EXCEPTIONS
    // If we have FP exceptions turned off, try to detect problems by looking
    // at the solution vector.
    for (int j = 0; j < N_; j++) {
      double value = VECTOR(i)[j];
      if (isinf(value)) return "FP error: Infinity encountered (backsolve rhs vector)";
      if (isnan(value)) return "FP error: NaN encountered (backsolve rhs vector)";
    }
#endif
  }
#ifdef USE_FP_EXCEPTIONS
  insolve = false;
#endif
  return 0;
}

double HessianMatrix::Get(int i, int j) {
  double *element = Access(i, j);
  if (element) return *element; else return 0;
}

void HessianMatrix::Print(int n, const char *format, const char *zero_format,
                          FILE *fout) {
  int sz = B_*N_+M_;
  int first = 0;
  int last = sz;
  if (n > 0) last = n; else if (n < 0) first = sz + n;
  for (int i = first; i < last; i++) {
    for (int j = first; j < last; j++) {
      double *element = Access(i, j);
      if (element && *element != 0) {
        fprintf(fout, format, *element);
      } else {
        fprintf(fout, zero_format, 0);
      }
    }
    fprintf(fout, "\n");
  }
  fprintf(fout, "\n");
}

void HessianMatrix::Save(const char *filename) {
  FILE *fout = fopen(filename, "wb");
  CHECK(fout);
  Print(0, "%.14g ", "0 ", fout);
  fclose(fout);
}

void HessianMatrix::Reset() {
  B_ = 0;
  N_ = 0;
  M_ = 0;
  NN_ = 0;
  NM_ = 0;
  diag_size_ = 0;
#if !defined(__GPO_HESSIAN_CACHE__)
  diag_data_ = 0;
  main_data_ = 0;
  edge_data_ = 0;
  diag_count_ = 0;
  main_count_ = 0;
  edge_count_ = 0;
#endif
  last_data_ = 0;
  last_count_ = 0;
}

void HessianMatrix::Free() {
#if !defined(__GPO_HESSIAN_CACHE__)
  delete[] diag_data_;
  delete[] main_data_;
  delete[] edge_data_;
#endif
  delete[] last_data_;
}

//***************************************************************************
// Hessian speed test

#if 0
#define B 50000

using namespace GPO;

int main() {
  HessianMatrix H;
  H.Initialize(B, 12, 6);
  for (int i = 0; i < B*12+6; i++) H.Add(i, i, 1);
  double vector[B*12+6];
  memset(vector, 0, sizeof(vector));

  H.Solve(vector);
}

#endif
