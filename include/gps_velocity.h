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

// GPS doppler velocity sensor model.

#ifndef __GPO_GPS_VELOCITY_H__
#define __GPO_GPS_VELOCITY_H__

#include "gpo/common.h"
#include "gpo/sensor.h"

namespace GPO {

class GPSSensor;

// A GPS doppler velocity sensor that measures the velocity at a subset of
// the poses in the state vector.

class GPSVelocitySensor : public Sensor {
 public:
  // The user must supply an array of the following measurements.
  struct Data {
    // The time step for which GPS velocity data is available. This must be > 0.
    int index;
    // Cartesian velocity m/s.
    double vx, vy, vz;
    // Standard deviations for the x, y and z variables, in meters.
    double sx, sy, sz;
  };

  // Create the sensor with the given measurements array, which must be valid
  // for the lifetime of this object
  GPSVelocitySensor(int num_samples, Data *samples, double stepsize);

  ~GPSVelocitySensor();

  // Sensor virtual functions
  void GetInfo(SensorInfo *info);
  int MeasurementIndex(int i);
  void ComputeErrorVector(int i, const Number *prev_step,
    const Number *this_step, const Number *next_step, const Number *globals,
    Number *error);

  // Set the position of the GPS antenna relative to the IMU measurement frame.
  void SetAntennaCenterOffset(double x, double y, double z);

private:
  int num_samples_;
  Data *samples_;
  double stepsize_;
  double center_[3];
};

}; // namespace GPO

#endif
