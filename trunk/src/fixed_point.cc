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

#include "gpo/fixed_point.h"
#include "gpo/pmath.h"
#include "gpo/rotation.h"
#include "gpo/auto_diff.h"

using namespace GPO;


FixedPoint::FixedPoint(int index, double pos[3], double euler[3],
                       double std_pos, double std_euler) {
  index_ = index;
  memcpy(pos_, pos, 3*sizeof(double));
  memcpy(euler_, euler, 3*sizeof(double));
  std_pos_ = std_pos;
  std_euler_ = std_euler;
}

FixedPoint::~FixedPoint() {
}

void FixedPoint::GetInfo(SensorInfo *info) {
  info->num_measurements = 1;
  info->error_size = 6;
}

int FixedPoint::MeasurementIndex(int i) {
  return index_;
}

void FixedPoint::ComputeErrorVector(int i, const Number *prev_step,
  const Number *this_step, const Number *next_step, const Number *globals,
  Number *error)
{
  error[0] = (this_step[STATE_X] - pos_[0]) / std_pos_;
  error[1] = (this_step[STATE_Y] - pos_[1]) / std_pos_;
  error[2] = (this_step[STATE_Z] - pos_[2]) / std_pos_;
  
  // Ensure that the Euler angle deltas are in the range -pi..pi.
  double e[3];
  for (int i=0; i < 3; i++) {
    double x = to_double(this_step[STATE_E1 + i]);
    double y = euler_[i];
    while (y-x > M_PI) y -= 2.0 * M_PI;
    while (y-x < -M_PI) y += 2.0 * M_PI;
    e[i] = y;
  }
  
  error[3] = (this_step[STATE_E1] - e[0]) / std_euler_;
  error[4] = (this_step[STATE_E2] - e[1]) / std_euler_;
  error[5] = (this_step[STATE_E3] - e[2]) / std_euler_;
}
