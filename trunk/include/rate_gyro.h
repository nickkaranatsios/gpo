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

// Rate gyro sensor model. The gyro model produces either three angular velocity
// readings or three delta angle readings per time step (an "integrating gyro").
// The gyro bias drift is also modeled.


#ifndef __GPO_RATE_GYRO_H__
#define __GPO_RATE_GYRO_H__

#include "gpo/sensor.h"

namespace GPO {

class RateGyroSensor : public Sensor {
 public:
  // The user must supply an array of the following measurements.
  struct Data {
    double rx, ry, rz;  // Rate vector, in radians/s (or delta angle in radians)
    double std_dev;     // Noise standard deviation for rx,ry,rz
  };

  // Create the sensor with the given measurements array, which must be
  // allocated on the heap. The samples[i] entry corresponds to time step i.
  // The first entry of this array is ignored.
  // This object takes ownership of the array and will delete it when destroyed.
  // C is the sensor noise standard deviation.
  // D is the drifting bias *per-step* noise standard deviation.
  RateGyroSensor(int num_samples, Data *samples, double stepsize,
                 double D);

  ~RateGyroSensor();

  void SetSpin(double x, double y, double z);

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
  Vector3_d spin_;      // Earth spin vector
  double D_;            // Drifting bias *per-step* noise standard deviation
};

}; // namespace GPO

#endif
