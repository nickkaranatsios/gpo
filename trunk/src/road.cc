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

//@@@TODO: Combine this with the axel constraint?

#include "gpo/road.h"
#include "gpo/pmath.h"
#include "gpo/rotation.h"
#include "gpo/auto_diff.h"

using namespace GPO;

RoadConstraint::RoadConstraint(int num_samples, WheelEncoder::Data *samples, double stepsize, double noise_scale) {
  num_samples_ = num_samples;
  samples_ = samples;
  stepsize_ = stepsize;
  noise_scale_ = noise_scale;
  min_noise_ = 0.03;
  road_[0] = 0;
  road_[1] = 0;
  road_[2] = 0;
}

RoadConstraint::~RoadConstraint() {
}

void RoadConstraint::SetRoadVector(Vector3_d &road) {
  memcpy(road_, road, sizeof(road));
  // Normalize the road vector
  double length = sqrt(sqr(road_[0]) + sqr(road_[1]) + sqr(road_[2]));
  for (int i = 0; i < 3; i++) road_[i] /= length;
}

void RoadConstraint::GetInfo(SensorInfo *info) {
  info->num_measurements = num_samples_;
  info->error_size = 1;
}

int RoadConstraint::MeasurementIndex(int i) {
  return samples_[i].index;
}

void RoadConstraint::ComputeErrorVector(int i, const Number *prev_step,
  const Number *this_step, const Number *next_step, const Number *globals,
  Number *error)
{
  Number R[9], dot_p[3], road_world[3];
  double noise = noise_scale_ * samples_[i].wheel_velocity;
  if (noise < min_noise_) noise = min_noise_;
  // Compute a central difference approximation to dp/dt
  for (int j = 0; j < 3; j++) {
    dot_p[j] = (next_step[j] - prev_step[j]) / (2.0 * stepsize_);
  }
  // Compute the road vector in the world frame
  eulerR(this_step + STATE_E1, R);
  FAST_MULTIPLY0_331(road_world, = , R, road_);
  // Compute road error
  error[0] = FAST_DOT(road_world, dot_p);
  error[0] /= noise;
}
