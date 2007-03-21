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

The rate gyro measurement vector contains:
     3*(n-1) rate gyro measurements
     3*(n-1) bias measurements

The rate gyro bias constraint says that the bias drift
velocity is zero (plus noise): (bias(i+1) - bias(i))/h = 0

Note that the first entry of the measurements array is ignored.

*/

#include "gpo/common.h"
#include "gpo/rate_gyro.h"
#include "gpo/pmath.h"
#include "gpo/rotation.h"
#include "gpo/auto_diff.h"

using namespace GPO;

// Define this macro to model the rate gyro bias, comment it out to ignore
// the bias.
//#define MODEL_BIAS   //@@@TODO: Uncomment this line?

// Define this to model an integrating gyro such as the Honeywell HG1700,
// undefine it to model a generic gyro.
#define INTEGRATING_GYRO


RateGyroSensor::RateGyroSensor(int num_samples, Data *samples, double stepsize,
                               double D) {
  CHECK(num_samples > 0 && samples && stepsize > 0 && D > 0);
  num_samples_ = num_samples;
  samples_ = samples;
  stepsize_ = stepsize;
  spin_[0] = 0;
  spin_[1] = 0;
  spin_[2] = 0;
  D_ = D;
}

RateGyroSensor::~RateGyroSensor() {
}

void RateGyroSensor::SetSpin(double x, double y, double z) {
  spin_[0] = x;
  spin_[1] = y;
  spin_[2] = z;
#ifdef INTEGRATING_GYRO
  spin_[0] *= stepsize_;
  spin_[1] *= stepsize_;
  spin_[2] *= stepsize_;
#endif
}

void RateGyroSensor::GetInfo(SensorInfo *info) {
#ifdef MODEL_BIAS
  info->num_states = 3;                         // Bias states
  info->error_size = 6;
#else
  info->error_size = 3;
#endif
  info->num_measurements = num_samples_ - 1;
}

int RateGyroSensor::MeasurementIndex(int i) {
  return i+1;
}

void RateGyroSensor::ComputeErrorVector(int i, const Number *prev_step,
  const Number *this_step, const Number *next_step, const Number *globals,
  Number *error)
{
  // ********** Integrating gyro
#ifdef INTEGRATING_GYRO
  //@@@TODO: Add bias state to the integrating gyro model.
  //@@@TODO: Should we use rotation matrix half way between R1 and R2 in
  //         some places?

  // The transformation corresponding to this measurement is between the
  // previous time step and this time step.
  //
  // For roll, pitch and yaw defined by deltas in the euler angles, we have:
  //   * roll rotation vector is the middle column of R1 (or R2).
  //   * pitch rotation vector is (cos(yaw),-sin(yaw),0) (global).
  //   * yaw rotation vector is (0,0,-1) (global).

  // Compute R1, the orientation at the previous timestep.
  Number R1[9];
  eulerR(prev_step + STATE_E1, R1);

  // Compute the change in pitch[0], roll[1] and yaw[2] between the previous
  // and current time steps. Adjust for discontinuities.
  Number delta[3];
  delta[0] = this_step[STATE_E1+0] - prev_step[STATE_E1+0];
  delta[1] = this_step[STATE_E1+1] - prev_step[STATE_E1+1];
  delta[2] = this_step[STATE_E1+2] - prev_step[STATE_E1+2];
  for (int j = 0; j < 3; j++) {
    double d = to_double(delta[j]);
    if (d > 3.0) {
      delta[j] -= 2.0 * M_PI;
    } else if (d < -3.0) {
      delta[j] += 2.0 * M_PI;
    }
  }
  
  // Compute synthetic rate gyro measurements.
  Number average_yaw = (prev_step[STATE_E1+2] + this_step[STATE_E1+2]) * 0.5;
  Number pitch_x = cos(average_yaw);    // X component of pitch vector
  Number pitch_y = -sin(average_yaw);   // Y component of pitch vector
  Number px = pitch_x * delta[0];       // (px,py,0) is rotation vector for
  Number py = pitch_y * delta[0];       //   change in pitch.
  // Map pitch rotation vector to the three gyros as R' * [px;py;0].
  // Map roll rotation vector to the three gyros as R' * R(:,2)*delta[1].
  // Map yaw rotation vector to the three gyros as R' * [0;0;-1]*delta[2].
  error[0] = px*R1[0] + py*R1[3] - R1[6]*delta[2];
  error[1] = px*R1[1] + py*R1[4] - R1[7]*delta[2] + delta[1];
  error[2] = px*R1[2] + py*R1[5] - R1[8]*delta[2];
  FAST_MULTIPLY1_331(error, += , R1, spin_);

  // Compute difference between current estimated and actual rates, and adjust
  // for standard deviations.
  error[0] -= samples_[i].rx;
  error[1] -= samples_[i].ry;
  error[2] -= samples_[i].rz;
  for (int j = 0; j < 3; j++) {
    error[j] /= samples_[i].std_dev;
  }

#else // INTEGRATING_GYRO
  // ********** Generic gyro

  Number e[3], dot_e[3], R[9];
  const Number *step1 = prev_step;
  const Number *step2 = this_step;
#ifdef MODEL_BIAS
  Number bias[3];
  const Number *bias1 = step1 + TOffset();
  const Number *bias2 = step2 + TOffset();
#endif
  // The leapfrog euler angles are the averages
  for (int j = 0; j < 3; j++) {
    e[j] = (step1[STATE_E1 + j] + step2[STATE_E1 + j])*0.5;
  }
#ifdef MODEL_BIAS
  // The leapfrog biases are the averages
  for (int j = 0; j < 3; j++) {
    bias[j] = (bias1[j] + bias2[j])*0.5;
  }
#endif
  // The leapfrog euler angle discrete derivatives
  for (int j = 0; j < 3; j++) {
    dot_e[j] = (step2[STATE_E1 + j] - step1[STATE_E1 + j]) / stepsize_;
  }
  // Interpolate adjacent samples
  double interp_samples[3];
  interp_samples[0] = 0.5*(samples_[i-1].rx + samples_[i].rx);
  interp_samples[1] = 0.5*(samples_[i-1].ry + samples_[i].ry);
  interp_samples[2] = 0.5*(samples_[i-1].rz + samples_[i].rz);
  // Compute the rate error
  eulerR(e, R);
  eulerRTEv(e, dot_e, error);
  FAST_MULTIPLY1_331(error, += , R, spin_);
  for (int j = 0; j < 3; j++) {
#ifdef MODEL_BIAS
    error[j] -= bias[j];
#endif
    error[j] -= interp_samples[j];
    error[j] /= samples_[i].stddev;
  }

#ifdef MODEL_BIAS
  // Bias. Note that to keep things simple we don't compute the difference of
  // the last two pose points in the state vector, but this doesn't affect the
  // final solution much.
  for (int j = 0; j < 3; j++) {
    error[3+j] = bias2[j] - bias1[j];
    error[3+j] /= D_;
  }
#endif
#endif // INTEGRATING_GYRO
}
