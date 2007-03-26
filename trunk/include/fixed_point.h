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

#ifndef __GPO_FIXED_POINT_H__
#define __GPO_FIXED_POINT_H__

#include "gpo/sensor.h"

namespace GPO {

class FixedPoint : public Sensor {
 public:
  // Create a fixed point constraint at measurement index 'i' (>1), with the
  // given position, orientation and standard deviations.
  FixedPoint(int index, double pos[3], double euler[3], double std_pos,
             double std_euler);
  ~FixedPoint();

  // Sensor virtual functions
  void GetInfo(SensorInfo *info);
  int MeasurementIndex(int i);
  void ComputeErrorVector(int i, const Number *prev_step,
    const Number *this_step, const Number *next_step, const Number *globals,
    Number *error);

private:
  int index_;                   // Measurement index that is fixed
  double pos_[3];               // Fixed position
  double euler_[3];             // Fixed orientation (Euler angles)
  double std_pos_;              // Position standard deviation (m)
  double std_euler_;            // Euler angle standard deviation (radians)
};

}; // namespace Poser

#endif
