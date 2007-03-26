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

// Matrix class for unit testing only. The API and implementation does
// everything in the simplest possible way, there is no optimization at all.
// Correctness is *much* more importance here.

#ifndef __GPO_TEST_MATRIX_H__
#define __GPO_TEST_MATRIX_H__

#include "common.h"

namespace GPO {

class TestMatrix {
 public:
  // Make default 0x0 matrix.
  TestMatrix();

  // Construct zero matrix of the given size
  TestMatrix(int rows, int cols);

  // Construct a copy of the given matrix
  TestMatrix(const TestMatrix &);

  // Create a copy of the given data - element (i,j) is data[i*rowskip+j*colskip]
  TestMatrix(int rows, int cols, double *data, int rowskip, int colskip);

  // Destructor
  ~TestMatrix();

  // Accessors
  int Rows() const { return n_; }
  int Cols() const { return m_; }
  double *Data() { return data_; }

  // Set the size of the matrix, clearing all elements to 0.
  void SetSize(int rows, int cols);

  // Reference an element (element 0,0 is at the top left of the matrix).
  double & operator() (int i, int j);

  // Copy one matrix to another.
  void operator = (const TestMatrix &);

  // Set all entries in the matrix to the given scalar.
  void operator = (double);

  // Return the transposed matrix.
  TestMatrix Transpose();

  // Return a permuted submatrix of this matrix, made up of the rows in p
  // and the columns in q. p has np elements, q has nq elements.
  TestMatrix Select(int np, int *p, int nq, int *q);

  // Matrix operators
  TestMatrix operator + (const TestMatrix &);
  TestMatrix operator - (const TestMatrix &);
  TestMatrix operator - ();
  TestMatrix operator * (const TestMatrix &);
  void operator += (const TestMatrix &);
  void operator -= (const TestMatrix &);

  // Clear the upper and lower triangle of the matrix.
  void ClearUpperTriangle();
  void ClearLowerTriangle();

  // Set all elements of the matrix to a random number in [-range,range].
  void MakeRandom(double range);

  // Print the matrix data to the given file with the given per-element format.
  void Print(char *fmt = "%10.4f ", FILE *f = stdout);

  // Print the nonzeros of the matrix to the given file.
  void PrintNonzeros(double tolerance = 1e-8, int blocksize = 0,
                     FILE *f = stdout);

  // Return the maximum difference between this matrix and M.
  double MaxDifference(const TestMatrix &M);

  // Return sqrt(sum(sum((this_matrix - M).^2)))
  double Distance(const TestMatrix &M);

 private:
  int n_, m_;           // Matrix dimension (n rows, m columns)
  double *data_;        // n*m elements allocated on the heap
};

}; // namespace GPO

#endif
