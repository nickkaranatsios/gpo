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

#include "gpo/gps.h"
#include "gpo/pmath.h"
#include "gpo/rotation.h"
#include "gpo/auto_diff.h"

using namespace GPO;


GPSSensor::GPSSensor(int num_samples, Data *samples) {
  CHECK(num_samples > 0 && samples);
  num_samples_ = num_samples;
  samples_ = samples;
  center_[0] = 0;
  center_[1] = 0;
  center_[2] = 0;
}

GPSSensor::~GPSSensor() {
}

void GPSSensor::GetInfo(SensorInfo *info) {
  info->num_measurements = num_samples_;
  info->error_size = 3;
}

int GPSSensor::MeasurementIndex(int i) {
  return samples_[i].index;
}

void GPSSensor::ComputeErrorVector(int i, const Number *prev_step,
  const Number *this_step, const Number *next_step, const Number *globals,
  Number *error)
{
  Number R[9];
  const Number *e = this_step + STATE_E1;
  eulerR(e, R);
  FAST_MULTIPLY0_331(error, = , R, center_)
  error[0] = (error[0] + this_step[STATE_X] - samples_[i].x) / samples_[i].sx;
  error[1] = (error[1] + this_step[STATE_Y] - samples_[i].y) / samples_[i].sy;
  error[2] = (error[2] + this_step[STATE_Z] - samples_[i].z) / samples_[i].sz;
}

void GPSSensor::SetAntennaCenterOffset(double x, double y, double z) {
  center_[0] = x;
  center_[1] = y;
  center_[2] = z;
}
