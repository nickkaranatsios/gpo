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

// GPS sensor model.

#ifndef __GPO_GPS_H__
#define __GPO_GPS_H__

#include "gpo/sensor.h"

namespace GPO {

// A GPS sensor that measures some (not all) of the poses.

class GPSSensor : public Sensor {
 public:
  // The user must supply an array of the following measurements.
  struct Data {
    // The time step for which GPS data is available (must be > 0).
    int index;
    // Cartesian coordinates in meters (NOTE: not lattitide, longitude etc).
    double x, y, z;
    // Standard deviations for the x, y and z variables, in meters.
    double sx, sy, sz;
  };

  // Create the sensor with the given measurements array, which must be
  // allocated on the heap.
  GPSSensor(int num_samples, Data *samples);

  ~GPSSensor();

  // Sensor virtual functions
  void GetInfo(SensorInfo *info);
  int MeasurementIndex(int i);
  void ComputeErrorVector(int i, const Number *prev_step,
    const Number *this_step, const Number *next_step, const Number *globals,
    Number *error);

  // Return the position of the GPS antenna center vector in the state vector.
  int GetAntennaCenterOffset() const;

private:
  int num_samples_;
  Data *samples_;
};

}; // namespace GPO

#endif
