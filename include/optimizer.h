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

// The least-squares optimizer core.

#ifndef __GPO_OPTIMIZER_H__
#define __GPO_OPTIMIZER_H__

#include "gpo/test_matrix.h"

namespace GPO {

class Sensor;
class SensorInfo;
class HessianMatrix;


class PoseOptimizer {
 public:
  PoseOptimizer();
  ~PoseOptimizer();

  // Set the verbosity level, which controls the messages that are printed
  // during optimization. 0 prints nothing, 1 prints only errors, and 2 prints
  // informational messages as well.
  void SetVerbosity(int verbosity) { verbosity_ = verbosity; }

  // Set the maximum memory usage in bytes.
  void SetMaxMemory(int max_memory) { max_memory_ = max_memory; }

  // Add a new sensor. This object takes ownership of the sensor, all sensors
  // are deleted automatically when this object is destroyed.
  void AddSensor(Sensor *sensor);

  // Run the optimization. An approximation to the solution is given in the
  // state vector, or if no approximation is known then this vector can be 0.
  // This might take a while. Go grab a cup of coffee or something.
  // If 0 is returned the optimization was successful, otherwise a human
  // readable error message is returned.
  const char *Optimize(double *state);

  // Get information about the size of the state vector for the currently
  // attached sensors.
  void GetStateVectorInfo(int *num_timesteps, int *block_size, int *global_size) const;

private:
  Sensor *first_sensor_;        // Linked list of sensors
  int verbosity_;               // Verbosity level
  size_t max_memory_;           // Maximum memory usage (bytes, 0 = infinity)

  // Assign toffset and goffset to each sensor.
  void AssignSensorOffsets() const;

  // Validate sensor measurement indexes - they should all be greater than zero
  // and less than the largest time step index.
  void ValidateMeasurementIndexes(int num_timesteps) const;
};

}; // namespace GPO

#endif
