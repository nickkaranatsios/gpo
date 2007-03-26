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
  memcpy(pos_, pos, sizeof(pos));
  memcpy(euler_, euler, sizeof(euler));
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
  error[3] = (this_step[STATE_E1] - euler_[0]) / std_euler_;
  error[4] = (this_step[STATE_E2] - euler_[1]) / std_euler_;
  error[5] = (this_step[STATE_E3] - euler_[2]) / std_euler_;
}
