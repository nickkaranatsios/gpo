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

// Fast functions for block factor, solve and multiply. All matrices are stored
// in either row-wise format or the following lower triangular compressed
// format:
//      [ 0       ]
//      [ 1 2     ]
//      [ 3 4 5   ]
//      [ ...     ]


#ifndef __GPO_FAST_BLOCK_H__
#define __GPO_FAST_BLOCK_H__

namespace GPO {

// Compute the Cholesky factorization L*L' = A, where A and L are n*n
// matrices stored in lower triangular compressed format.
void factor(double *A, int n);

// Solve by forward substitution L*X' = B', L is n*n, B is m*n.
// L is stored in lower triangular compressed format.
void solve(const double *L, double *B, int n, int m);

// Solve backward substitution L'*x = b, L is n*n, B is n*1.
// L is stored in lower triangular compressed format.
void solveT(const double *L, double *B, int n);

// Compute A += scale*B'*C, where A is p*r, B is q*p, C is q*r.
void multiply1(double *A, double scale, const double *B, const double *C,
               int p, int q, int r);

// Compute A += scale*B*C', where A is p*r, B is p*q, C is r*q.
void multiply2(double *A, double scale, const double *B, const double *C,
               int p, int q, int r);

// Compute A += scale*B*C', where A is p*p, B is p*q, C is p*q.
// A is stored in lower triangular compressed format.
void multiply2_lt(double *A, double scale, const double *B, const double *C,
                  int p, int q);

// Reference functions, for testing
void solve_ref(const double *L, double *B, int n, int m);
void multiply2_ref(double *A, double scale, const double *B, const double *C,
                   int p, int q, int r);
void multiply2_lt_ref(double *A, double scale, const double *B, const double *C,
                      int p, int q);

}; // namespace GPO

#endif
