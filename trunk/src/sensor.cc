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

#include "gpo/sensor.h"
#include "gpo/fast_block.h"
#include "gpo/hessian.h"
#include "gpo/auto_diff.h"


using namespace GPO;

Sensor::Sensor() {
  next_ = 0;
  toffset_ = 0;
  goffset_ = 0;
}

Sensor::~Sensor() {
}



double Sensor::ComputeDifferential(StateInfo &state, double *gradient,
  HessianMatrix *H, int block, int total_blocks)
{
  SensorInfo info;
  memset(&info, 0, sizeof(info));
  GetInfo(&info);

  int num_deps = 3*state.tsize + state.gsize;
  Number prev_step[state.tsize];
  Number this_step[state.tsize];
  Number next_step[state.tsize];
  Number globals[state.gsize];
  Number error[info.error_size];
  double JT[info.error_size * num_deps];        // J^T stored by row
  double error_d[info.error_size];
  IndexSet Jcolset[num_deps];                   // Sets of nonzero indexes in each column of J
  double E = 0;

  // Compute the range of measurements to process
  int block_size = ((info.num_measurements-1) / total_blocks) + 1;
  int start_index = block * block_size;
  if (start_index >= info.num_measurements) {
    // The total blocks was probably larger than the number of measurements.
    return 0;
  }
  int end_index = start_index + block_size;
  if (end_index >= info.num_measurements) end_index = info.num_measurements;

  for (int i = start_index; i < end_index; i++) {
    int index = MeasurementIndex(i);
    int j = index * state.tsize;

    // Set dependent variable values and dependencies.
    for (int k = 0; k < state.tsize; k++) {
      prev_step[k].SetDependent(k, state.x[j - state.tsize + k]);
      this_step[k].SetDependent(k + state.tsize, state.x[j + k]);
      next_step[k].SetDependent(k + 2*state.tsize, state.x[j + state.tsize + k]);
    }
    for (int k = 0; k < state.gsize; k++) {
      globals[k].SetDependent(k + 3*state.tsize, state.x[state.size - state.gsize + k]);
    }

    // Compute the error vector and all its partial derivatives.
    ComputeErrorVector(i, prev_step, this_step, next_step, globals, error);

    // Extract the sets of nonzero indexes in each column of J.
    // In most cases this could be done just ONCE, but how can we be sure of
    // that? This is not too expensive anyway.
    memset(Jcolset, 0, sizeof(Jcolset));
    for (int k = 0; k < info.error_size; k++) {
      IndexSet set = error[k].GetNonzeroDerivativeSet();
      FOR_ALL_INDEXES(l, set)
        Jcolset[l] |= (1 << k);
      }
    }

    // Extract the Jacobian matrix and create a column index set for all
    // columns that are nonzero in all Jacobian rows.
    memset(JT, 0, info.error_size * num_deps * sizeof(double));  //@@@TODO: Really necessary?
    IndexSet Jset = 0;
    for (int k = 0; k < info.error_size; k++) {
      error[k].GetNonzeroDerivatives(num_deps, JT+k, info.error_size);
      Jset |= error[k].GetNonzeroDerivativeSet();
    }

    // Copy the error vector to a double array and accumulate scalar error
    for (int k = 0; k < info.error_size; k++) {
      error_d[k] = to_double(error[k]);
      E += error_d[k] * error_d[k];
    }

    // Compute the gradient vector = J^T * inv(C) * (f(x) - z)
    int kofs = j - state.tsize;        // Add this to k to get index into gradient
    FOR_ALL_INDEXES(k, Jset)
      if (k >= 3*state.tsize) kofs = state.size - state.gsize - 3*state.tsize;
      double sum = 0;
      double *pk = JT + k * info.error_size;
      double *pe = error_d;
      int m = info.error_size;
      while (m >= 3) {          //@@@TODO: This optimization hardly makes a difference
        sum += pk[0]*pe[0] + pk[1]*pe[1] + pk[2]*pe[2];
        pk += 3;
        pe += 3;
        m -= 3;
      }
      while (m > 0) {
        sum += (*pk++) * (*pe++);
        m--;
      }
      gradient[k + kofs] += sum;
    }

    // Compute the Hessian matrix = J^T * inv(C) * J (lower triangle only).
    // If J is sparse then H will be sparse and we shouldn't compute all
    // its values. This is especially true because sensors will not use
    // all the globals.
    FOR_ALL_INDEXES(k, Jset)
      int kblock = k/state.tsize;                       // Block number in k
      int ldofs = kblock*state.tsize;                   // Where the diagonal will start in index l
      FOR_ALL_INDEXES(l, Jset & index_range(k+1))       // Lower triangle
        int lblock = l/state.tsize;                     // Block number in l
        IndexSet nonzeros = Jcolset[k] & Jcolset[l];
        if (nonzeros) {
          double sum = 0;
          double *pk = JT + k * info.error_size;
          double *pl = JT + l * info.error_size;

          /*
          // An alternative method that may be quicker under some circumstances.
          FOR_ALL_INDEXES(m,nonzeros)
            sum += pk[m] * pl[m];
          }
          if (fabs(sum) > 1e-12) {              // Arbitrary tolerance
            H->Add (k + kofs, l + lofs, sum);
          }
          */

          while (nonzeros) {
            sum += (*pk++) * (*pl++);
            nonzeros >>= 1;
          }

          /*
          // The simple (but slow) way to add the sum to the Hessian matrix:
          // Outer loop:
          kofs = j - state.tsize;       // Add this to k to get index into H
          if(k >= 3*state.tsize) kofs = gstart - 3*state.tsize;
          // Inner loop:
          if(l >= 3*state.tsize) {
            H->Add(k + kofs, l + gstart - 3*state.tsize, sum);
          } else {
            H->Add(k + kofs, l + j - state.tsize, sum);
          }
          */

          // The faster way to add the sum to the Hessian matrix:
          if (k >= 3*state.tsize) {
            // Edge block or last block
            if (l >= 3*state.tsize) {
              // Last block
              H->AccessLastBlock(k-ldofs, l-ldofs)[0] += sum;
            } else {
              // Edge block
              H->AccessEdgeBlock(index + lblock-1, k-ldofs, l - lblock*state.tsize)[0] += sum;
            }
          } else {
            // Main block or diagonal block
            if (l >= ldofs) {
              // Diagonal block
              H->AccessDiagonalBlock(index + kblock-1, k-ldofs, l-ldofs)[0] += sum;
            } else {
              // Main block
              H->AccessOffDiagonalBlock(index + kblock-1, index + lblock-1, k-ldofs, l - lblock*state.tsize)[0] += sum;
            }
          }
        }
      }
    }
  }
  return E;
}
