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

// Base class for sensor models.

#ifndef __GPO_SENSOR_H__
#define __GPO_SENSOR_H__

#include "gpo/test_matrix.h"

namespace GPO {

class Number;
class HessianMatrix;

enum {
  // Standard offsets into state vector blocks.
  // We guarantee that the first six elements are always X,Y,Z,E1,E2,E3.
  STATE_X  = 0,
  STATE_Y,
  STATE_Z,
  STATE_E1,                     // 3D euler angles
  STATE_E2,
  STATE_E3,
  POSE_BLOCK_SIZE,              // Number of variables in the 'pose block'
  // Sensor-specific state variables come after this
};


// Information about a sensor that is passed from the sensor to the optimizer.

struct SensorInfo {
  // The number of extra state variables introduced per time step.
  int num_states;

  // The number of extra global state variables introduced (not per time step).
  int num_global_states;

  // The number of discrete measurements made by this sensor. Each measurement is
  // associated with one time step index.
  int num_measurements;

  // The size of the error vector for each measurement.
  int error_size;
};


// Information about the state vector that is passed to from the optimizer
// to the sensor Compute*() functions.

struct StateInfo {
  // The state vector itself.
  const double *x;

  // The size of the state vector
  int size;

  // The size of each time step's block in the state vector.
  int tsize;

  // The size of the global area of the state vector
  int gsize;
};


// A sensor object maps real-world measured data to constraints on the state
// vector.

class Sensor {
 public:
  Sensor();
  virtual ~Sensor();

  // Return information about this sensor in the SensorInfo structure. On entry
  // all fields of info are zeroed.
  virtual void GetInfo(SensorInfo *info) = 0;

  // For the i'th measurement, return the timestep index that it measures.
  // This must be greater than zero and less than the largest timestep - this
  // restriction allows us to avoid the complexity of handling edge effects in
  // the code. The largest timestep is one plus the maximum measurement index
  // of all sensors. The measurement indexes are assumed to be monotonically
  // increasing.
  virtual int MeasurementIndex(int i) = 0;

  // Compute the sensor error vector for measurement i, using automatic
  // differentiation. The error vector = C^(-0.5) * (f(x)-z)
  virtual void ComputeErrorVector(int i, const Number *prev_step,
    const Number *this_step, const Number *next_step, const Number *globals,
    Number *error) = 0;

  // Compute differential quantities for this sensor using automatic
  // differentiation: the gradient vector and the Hessian matrix.
  //   - The gradient vector is: J^T * inv(C) * (f(x) - z)
  //   - The Hessian matrix is: J^T * inv(C) * J
  // The computed values are added (not written) to H and gradient.
  // Also return the error: (f(x)-z)' * inv(C) * (f(x)-z).
  // For the default arguments block=0 and total_blocks=1, all measurements are
  // used in the computation. However a subset of all measurements can be used
  // by specifying total_blocks > 1, so the measurements in block number 'block'
  // of 'total_blocks' total blocks are used. This allows the caller to get more
  // locality of reference in writes to the big Hessian matrix from many
  // sensors.
  double ComputeDifferential(StateInfo &state, double *gradient,
    HessianMatrix *H, int block = 0, int total_blocks = 1);

  // State vector information functions. Note that these can ONLY be called from
  // the Compute*() functions, as they only return valid data during
  // optimization.
  inline int TOffset() const { return toffset_; }
  inline int GOffset() const { return goffset_; }

  // Set the state vector information. These function should only be used for
  // testing.
  inline void SetTOffset(int ofs) { toffset_ = ofs; }
  inline void SetGOffset(int ofs) { goffset_ = ofs; }

private:
  friend class PoseOptimizer;
  Sensor *next_;        // Linked list of sensors
  int toffset_;         // Offset of sensor's variables in time step block.
  int goffset_;         // Offset of sensor's global variables in global variable block.
};

}; // namespace GPO

#endif
