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

The accelerometer measurement vector contains:
     3*(n-2) accelerometer measurements
     3*(n-1) bias measurements

The accelerometer bias constraint says that the bias drift
velocity is zero (plus noise): (bias(i+1) - bias(i))/h = 0

Note that the first element of the measurements array is ignored.

*/

#include "gpo/common.h"
#include "gpo/accelerometer.h"
#include "gpo/pmath.h"
#include "gpo/rotation.h"
#include "gpo/auto_diff.h"

using namespace GPO;

// Define this macro to model the accelerometer bias, comment it out to ignore
// the bias.
#define MODEL_BIAS


AccelerometerSensor::AccelerometerSensor(int num_samples, Data *samples,
                                         double stepsize, double D) {
  CHECK(num_samples > 0 && samples && stepsize > 0 && D > 0);
  num_samples_ = num_samples;
  samples_ = samples;
  stepsize_ = stepsize;
  gravity_[0] = 0;
  gravity_[1] = 0;
  gravity_[2] = 0;
  D_ = D;
}

AccelerometerSensor::~AccelerometerSensor() {
  delete[] samples_;
}

void AccelerometerSensor::SetGravity(double x, double y, double z) {
  gravity_[0] = x;
  gravity_[1] = y;
  gravity_[2] = z;
}

void AccelerometerSensor::GetInfo(SensorInfo *info) {
#ifdef MODEL_BIAS
  info->num_states = 3;                         // Bias states
  info->error_size = 6;
#else
  info->error_size = 3;
#endif
  info->num_measurements = num_samples_-1;
}

int AccelerometerSensor::MeasurementIndex(int i) {
  return i+1;
}

void AccelerometerSensor::ComputeErrorVector(int i, const Number *prev_step,
  const Number *this_step, const Number *next_step, const Number *globals,
  Number *error)
{
  Number dot_p[3], R[9];
#ifdef MODEL_BIAS
  const Number *bias = this_step + TOffset();
#else
  Number bias[3];
  bias[0] = 0;
  bias[1] = 0;
  bias[2] = 0;
#endif
  double h2 = 1.0 / (stepsize_*stepsize_);
  for (int j = 0; j < 3; j++) {
    dot_p[j] = (prev_step[j] - this_step[j]*2.0 + next_step[j]) * h2 + gravity_[j];
  }
  eulerR(this_step + STATE_E1, R);
  FAST_MULTIPLY1_331(error, = , R, dot_p);
  error[0] -= samples_[i+1].ax;
  error[1] -= samples_[i+1].ay;
  error[2] -= samples_[i+1].az;
  for (int j = 0; j < 3; j++) {
    error[j] -= bias[j];
    error[j] /= samples_[i+1].std_dev;
  }

#ifdef MODEL_BIAS
  // Bias. Note that to keep things simple we don't compute the difference of
  // the last two pose points in the state vector, but this doesn't affect the
  // final solution much.
  const Number *prev_bias  = prev_step + TOffset();
  for (int j = 0; j < 3; j++) {
    error[3+j] = bias[j] - prev_bias[j];
    error[3+j] /= D_;
  }
#endif
}
