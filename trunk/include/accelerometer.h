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

// Accelerometer sensor model.

#ifndef __GPO_ACCELEROMETER_H__
#define __GPO_ACCELEROMETER_H__

#include "gpo/common.h"
#include "gpo/sensor.h"

namespace GPO {

// Model of a 3-axis accelerometer that produces one reading per time step.
// This sensor also models the accelerometer bias drift.

class AccelerometerSensor : public Sensor {
 public:
  // The user must supply an array of the following measurements.
  struct Data {
    double ax, ay, az;          // In m/s^2
    double std_dev;             // Noise standard deviation for ax,ay,az
  };

  // Create the sensor with the given measurements array, which must be valid
  // for the lifetime of this object. The samples[i] entry corresponds to time
  // step i. The first entry of this array is ignored.
  // D is the drifting bias *per-step* noise standard deviation.
  AccelerometerSensor(int num_samples, Data *samples, double stepsize,
                      double D);

  ~AccelerometerSensor();

  // Set the external gravity vector
  void SetGravity(double x, double y, double z);

  // Sensor virtual functions
  void GetInfo(SensorInfo *info);
  int MeasurementIndex(int i);
  void ComputeErrorVector(int i, const Number *prev_step,
    const Number *this_step, const Number *next_step, const Number *globals,
    Number *error);

private:
  int num_samples_;
  Data *samples_;
  double stepsize_;
  Vector3_d gravity_;   // Gravity vector
  double D_;            // Drifting bias *per-step* noise standard deviation
};

}; // namespace GPO

#endif
