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

#include "gpo/axel.h"
#include "gpo/pmath.h"
#include "gpo/rotation.h"
#include "gpo/auto_diff.h"

using namespace GPO;

AxelConstraint::AxelConstraint(int num_samples, WheelEncoder::Data *samples, double stepsize, double noise_scale) {
  num_samples_ = num_samples;
  samples_ = samples;
  stepsize_ = stepsize;
  noise_scale_ = noise_scale;
  min_noise_ = 0.03;
  axel_[0] = 0;
  axel_[1] = 0;
  axel_[2] = 0;
  c_cross_l_[0] = 0;
  c_cross_l_[1] = 0;
  c_cross_l_[2] = 0;
}

AxelConstraint::~AxelConstraint() {
}

void AxelConstraint::SetVectors(Vector3_d &axel, Vector3_d &wheel) {
  memcpy(axel_, axel, sizeof(axel));
  // Normalize the axel vector
  double length = sqrt(sqr(axel_[0]) + sqr(axel_[1]) + sqr(axel_[2]));
  for (int i = 0; i < 3; i++) axel_[i] /= length;
  // Compute c_cross_l
  FAST_CROSS(c_cross_l_, = , wheel, axel_);
}

void AxelConstraint::GetInfo(SensorInfo *info) {
  info->num_measurements = num_samples_;
  info->error_size = 1;
}

int AxelConstraint::MeasurementIndex(int i) {
  return samples_[i].index;
}

void AxelConstraint::ComputeErrorVector(int i, const Number *prev_step,
  const Number *this_step, const Number *next_step, const Number *globals,
  Number *error)
{
  Number dot_p[3], dot_e[3], w[3], R[9], axel_world[3], c_cross_l_world[3];
  double noise = noise_scale_ * samples_[i].wheel_velocity;
  if (noise < min_noise_) noise = min_noise_;
  // Compute a central difference approximation to dp/dt
  for (int j = 0; j < 3; j++) {
    dot_p[j] = (next_step[j] - prev_step[j]) / (2.0 * stepsize_);
  }
  // Compute a central difference approximation to de/dt
  for (int j = 0; j < 3; j++) {
    dot_e[j] = (next_step[STATE_E1 + j] - prev_step[STATE_E1 + j]) / (2.0 * stepsize_);
  }
  // Compute an approximation to the angular velocity.
  eulerEv(this_step + STATE_E1, dot_e, w);
  // Compute the axel and c_cross_l vectors in the world frame
  eulerR(this_step + STATE_E1, R);
  FAST_MULTIPLY0_331(axel_world, = , R, axel_);
  FAST_MULTIPLY0_331(c_cross_l_world, = , R, c_cross_l_);
  // Compute axel slip
  error[0] = (FAST_DOT(axel_world, dot_p) + FAST_DOT(c_cross_l_world, w)) / noise;
}
