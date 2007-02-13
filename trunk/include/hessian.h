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

// Large sparse square symmetric positive definite 'Hessian' matrix.
// The sparsity pattern is block diagonal with 'edges':
//
//      [ XXX            Y ]
//      [ XXXX           Y ]
//      [ XXXXX          Y ]
//      [  XXXXX         Y ]
//      [   XXXXX        Y ]
//      [    XXXXX       Y ]
//      [      ...     ... ]
//      [       XXXXX    Y ]
//      [        XXXXX   Y ]
//      [         XXXXX  Y ]
//      [          XXXXX Y ]
//      [           XXXX Y ]
//      [            XXX Y ]
//      [                  ]
//      [ YYYYY...YYYYYY Z ]
//
// Where typically X is 12x12, Y is 6x12 and Z is 6x6.
// This class is designed to store and factor these matrices efficiently for
// dimensions up to 10-20 million (i.e. as much as will fit into machine
// memory). Actually with cached (out of core) storage we can factor much
// larger matrices, with a corresponding hit in speed.

#ifndef __GPO_HESSIAN_H__
#define __GPO_HESSIAN_H__

#include "gpo/common.h"
#include "gpo/cache.h"

// Define the following macro to use a block caching scheme for Hessian
// matrices that are too large to fit into memory (or even too large for the
// 4 Gb virtual address space).
//#define __GPO_HESSIAN_CACHE__


namespace GPO {

class HessianMatrix {
 public:
  HessianMatrix();
  ~HessianMatrix();

  // Set the size of the matrix, and zero it. This returns false if we would run
  // out of memory (in which case the object state is such that only
  // Initialize() can be subsequently called). The arguments are:
  //   - Number of matrix blocks is B
  //   - Matrix block size is N*N
  //   - Matrix edge size is N*M
  //   - Maximum memory usage is 'max_memory' bytes, 0 means "unlimited".
  // The total size is B*N+M.
  bool Initialize(int B, int N, int M, size_t max_memory = 0);

  // Set all elements to zero
  void SetZero();

  // Add 'x' to row i, column j. If i != j then 'x' is also added to the
  // symmetric position at row j column i. Indexes are zero-based.
  void Add(int i, int j, double x);

  // Return a pointer to element i,j, or return 0 if the element is not in
  // the sparsity pattern.
  double *Access(int i, int j);

  // Functions to access the stored blocks directly. These functions are faster
  // than Access() but you have to be aware of the sparsity pattern.
  double *AccessDiagonalBlock(int bi, int i, int j) {
    DCHECK(bi >= 0 && bi < B_);
    DCHECK(i >= 0 && i < N_);
    DCHECK(j >= 0 && j < N_);
    #ifdef __GPO_HESSIAN_CACHE__
      double *data = (double*) diag_cache_.Access(bi);
      return data + (i*(i+1))/2 + j;
    #else
      return diag_data_ + bi*diag_size_ + (i*(i+1))/2 + j;
    #endif
  }
  double *AccessOffDiagonalBlock(int bi, int bj, int i, int j) {
    DCHECK(bi >= 0 && bi < B_);
    DCHECK(bj >= 0 && bj >= bi-2 && bj < bi);
    DCHECK(i >= 0 && i < N_);
    DCHECK(j >= 0 && j < N_);
    #ifdef __GPO_HESSIAN_CACHE__
      double *data = (double*) main_cache_.Access((2*bi-3) + (bj-bi+2));
      return data + i*N_ + j;
    #else
      return main_data_ + ((2*bi-3) + (bj-bi+2)) * NN_ + i*N_ + j;
    #endif
  }
  double *AccessEdgeBlock(int bj, int i, int j) {
    DCHECK(bj >= 0 && bj < B_);
    DCHECK(i >= 0 && i < M_);
    DCHECK(j >= 0 && j < N_);
    #ifdef __GPO_HESSIAN_CACHE__
      double *data = (double*) edge_cache_.Access(bj);
      return data + i*N_ + j;
    #else
      return edge_data_ + bj*NM_ + i*N_ + j;
    #endif
  }
  double *AccessLastBlock(int i, int j) {
    DCHECK(i >= 0 && i < M_);
    DCHECK(j >= 0 && j < M_);
    return last_data_ + (i*(i+1))/2 + j;
  }

  // Multiply this matrix by 'x' and put the result in 'y'. Both vectors must
  // have size B*N+M. This is used for testing only, so it is slow-but-sure.
  void Multiply(double *x, double *y);

  // Solve the linear system A*x=b, where A is this matrix, 'x' is read from the
  // given vector and 'b' is written back to the given vector. The vector must
  // have the size B*N+M. An in-place Cholesky factorization is performed on
  // the matrix, so on exit the matrix data is changed. Return 0 if the
  // factorization was successful, or a human readable error message if some
  // error was encountered (e.g. overflow, divide by zero, or square root of
  // negative due to a badly conditioned matrix).
  const char *Solve(double *vector);

  // Get the element at row i, column j. This is used for testing only and
  // does not need to be fast.
  double Get(int i, int j);

  // Print the matrix to the given file (for testing). Since this matrix can
  // be extremely large, setting n > 0 prints only the top left n*n block, and
  // setting n < 0 prints only the bottom right (-n)*(-n) block.
  void Print(int n = 0, const char *format = "%10.3e ",
    const char *zero_format = " 0         ", FILE *fout = stdout);

  // Save the matrix to the given filename (for testing).
  void Save(const char *filename);

private:
  // Data is stored by block, and within blocks by row or partial row.
  int B_;                       // Number of blocks
  int N_;                       // Size of each block
  int M_;                       // Edge size
  int NN_;                      // = N_*N_
  int NM_;                      // = N_*M_
  int diag_size_;               // = (N_*(N_+1))/2
#ifdef __GPO_HESSIAN_CACHE__
  Cache diag_cache_;            // Storage of diagonal blocks
  Cache main_cache_;            // Storage of off-diagonal blocks
  Cache edge_cache_;            // Storage of edge blocks
#else
  double *diag_data_;           // Storage of diagonal blocks
  double *main_data_;           // Storage of off-diagonal blocks
  double *edge_data_;           // Storage of edge blocks
  int diag_count_;              // Array size of diag_data_
  int main_count_;              // Array size of main_data_
  int edge_count_;              // Array size of edge_data_
#endif
  double *last_data_;           // Storage of the bottom right block
  int last_count_;              // Array size of last_data_

  void Reset();                 // Reset all fields to zero
  void Free();                  // Free all matrix data
};

}; // namespace GPO

#endif
