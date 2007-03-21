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

#include "gpo/common.h"
#include "gpo/gps_velocity.h"
#include "gpo/gps.h"
#include "gpo/pmath.h"
#include "gpo/rotation.h"
#include "gpo/fast_block.h"
#include "gpo/auto_diff.h"

using namespace GPO;


GPSVelocitySensor::GPSVelocitySensor(int num_samples, Data *samples,
                                     double stepsize, GPSSensor *gps_position) {
  CHECK(num_samples > 0 && samples);

  num_samples_ = num_samples;
  samples_ = samples;
  stepsize_ = stepsize;
  gps_position_ = gps_position;
}

GPSVelocitySensor::~GPSVelocitySensor() {
}

void GPSVelocitySensor::GetInfo(SensorInfo *info) {
  // We need a global vector for the center of GPS relative to the IMU.
  if (gps_position_) {
    // Share the center vector between the GPS position and velocity sensors.
    info->num_global_states = 0;
  } else {
    // Allocate our own vector.
    info->num_global_states = 3;
  }
  info->num_measurements = num_samples_;
  info->error_size = 3;
}

int GPSVelocitySensor::MeasurementIndex(int i) {
  return samples_[i].index;
}

void GPSVelocitySensor::ComputeErrorVector(int i, const Number *prev_step,
  const Number *this_step, const Number *next_step, const Number *globals,
  Number *error)
{
  Number R[9], w[3], dot_e[3], c[3];
  int goffset = gps_position_ ? gps_position_->GOffset() : GOffset();
  const Number *center_gps = globals + goffset;
  const Number *e = this_step + STATE_E1;
  // Compute a central difference approximation to de/dt
  for (int j = 0; j < 3; j++) {
    dot_e[j] = (next_step[STATE_E1+j] - prev_step[STATE_E1+j]) / (2.0 * stepsize_);
  }
  // Compute an approximation to the angular velocity.
  eulerEv(e, dot_e, w);
  // Compute the center_gps vector in the world frame.
  eulerR(e, R);
  FAST_MULTIPLY0_331(c, = , R, center_gps);
  // Compute the total velocity
  FAST_CROSS(error, = , w, c);
  for (int j = 0; j < 3; j++) {
    error[j] += (next_step[j] - prev_step[j]) / (2.0 * stepsize_);
  }
  // Compute the velocity error
  error[0] -= samples_[i].vx;
  error[1] -= samples_[i].vy;
  error[2] -= samples_[i].vz;
  error[0] /= samples_[i].sx;
  error[1] /= samples_[i].sy;
  error[2] /= samples_[i].sz;
}
