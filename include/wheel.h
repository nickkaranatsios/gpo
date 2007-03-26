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

// Wheel encoder model.

#ifndef __GPO_WHEEL_H__
#define __GPO_WHEEL_H__

#include "gpo/common.h"
#include "gpo/sensor.h"

namespace GPO {

class WheelEncoder : public Sensor {
 public:
  // The user must supply an array of the following measurements.
  struct Data {
    // The time step for which wheel encoder data is available (must be > 0).
    int index;
    // Wheel velocity, in radians/second
    double wheel_velocity;
  };

  // Create the sensor with the given measurements array, which must be valid
  // for the lifetime of this object. The sensor noise standard deviation is
  // noise_scale times the absolute value of the wheel velocity measurement at
  // each timestep.
  WheelEncoder(int num_samples, Data *samples, double stepsize,
    double noise_scale);

  ~WheelEncoder();

  // Set the axel vector, the road normal vector, and the wheel connection
  // vector relative to the pose frame of reference. The axel and normal vectors
  // must be perpendicular to each other. The wheel vector can be any point
  // along the axel.
  void SetVectors(Vector3_d &axel, Vector3_d &road_normal, Vector3_d &wheel);

  // Set/get the minumum noise that the sensor will produce for each sample.
  void SetMinimumNoise(double n) { min_noise_ = n; }
  double GetMinimumNoise() const { return min_noise_; }

  // Sensor virtual functions
  void GetInfo(SensorInfo *info);
  int MeasurementIndex(int i);
  void ComputeErrorVector(int i, const Number *prev_step,
    const Number *this_step, const Number *next_step, const Number *globals,
    Number *error);

  // Return the position of the wheel smallness (1/radius) in the state vector.
  int GetWheelSmallnessOffset() const;

private:
  int num_samples_;
  Data *samples_;
  double stepsize_;
  double noise_scale_;          // Noise scale
  double min_noise_;            // Minimum noise (1-sigma)
  Vector3_d pvec_, wvec_;       // sensor = pvec * dp/dt + wvec * w
};

}; // namespace GPO

#endif
