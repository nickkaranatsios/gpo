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

// Road constraint model, which constrains one degree of freedom in
// linear velocity.

#ifndef __GPO_ROAD_H__
#define __GPO_ROAD_H__

#include "gpo/common.h"
#include "gpo/sensor.h"
#include "gpo/wheel.h"

namespace GPO {

class RoadConstraint : public Sensor {
 public:
  // Create the sensor with the given wheel encoder measurements array, which
  // must be allocated on the heap. This object DOES NOT take ownership of the
  // array, i.e. this object will not delete it.
  // The sensor noise standard deviation is noise_scale times the absolute
  // value of the wheel velocity measurement at each timestep.
  RoadConstraint(int num_samples, WheelEncoder::Data *samples, double stepsize, double noise_scale);

  ~RoadConstraint();

  // Set the road vector relative to the pose frame of reference.
  void SetRoadVector(Vector3_d &road);

  // Set/get the minumum noise that the sensor will produce for each sample.
  void SetMinimumNoise(double n) { min_noise_ = n; }
  double GetMinimumNoise() const { return min_noise_; }

  // Sensor virtual functions
  void GetInfo(SensorInfo *info);
  int MeasurementIndex(int i);
  void ComputeErrorVector(int i, const Number *prev_step,
    const Number *this_step, const Number *next_step, const Number *globals,
    Number *error);

private:
  int num_samples_;
  WheelEncoder::Data *samples_;
  double stepsize_;
  double noise_scale_;          // Noise scale
  double min_noise_;            // Minimum noise (1-sigma)
  Vector3_d road_;              // Road vector, normalized to unit length.
};

}; // namespace GPO

#endif
