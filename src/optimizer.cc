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

#include "gpo/optimizer.h"
#include "gpo/sensor.h"
#include "gpo/hessian.h"
#include "gpo/rotation.h"

//@@@TODO: Note: this code *could* use already computed values of H when steps
// are rejected and new steps must be taken from the old step. This is faster
// for small matrices, but our H matrix is so large that we can't have two copies.

using namespace GPO;

//***************************************************************************
// PoseOptimizer

PoseOptimizer::PoseOptimizer() {
  first_sensor_ = 0;
  verbosity_ = 0;
  max_memory_ = 0;
}

PoseOptimizer::~PoseOptimizer() {
  // Delete all attached sensors
  Sensor *snext;
  for (Sensor *s = first_sensor_; s; s = snext) {
    snext = s->next_;
    delete s;
  }
}

void PoseOptimizer::AddSensor(Sensor *sensor) {
  sensor->next_ = first_sensor_;
  first_sensor_ = sensor;
  AssignSensorOffsets();
}

const char *PoseOptimizer::Optimize(double *state) {
  // Constants that control the behavior of the optimizer
  const int kMaxIter = 300;             // Maximum number of iterations to do
  const double kLambdaStart = 1e-5;     // Starting value of lambda
  const double kLambdaScale = 10.0;     // How to scale lambda
                                        // Note: 10^0.2=1.58489319246111
  const double kLambdaBigScale = 10.0;  // A bigger scale for lambda
  const double kMaxLambda = 1e50;       // Maximum value of lambda
  const double kMinLambda = 1e-20;      // Minimum value of lambda
  const double kScaleTol = 1e-6;        // Convergence tolerance (scale)
  const double kAbsTol = 1e-9;          // Convergence tolerance (absolute)
  const int kMeasurementBlocks = 1000;  // @@@TODO: Arbitrary

  // Misc variables
  double *original_state = state;       // Because we swap state vectors
  const char *return_value = 0;         // The return value of this function

  // Validate the sensors and collect state vector size information
  int num_timesteps, block_size, num_global;
  GetStateVectorInfo(&num_timesteps, &block_size, &num_global);

  // Allocate various matrices and vectors
  int state_size = num_timesteps*block_size + num_global;
  HessianMatrix H;
  H.Initialize(num_timesteps, block_size, num_global, max_memory_);
  double *gradient = new double[state_size];
  double *trial_state = new double[state_size];

  // Print stats
  GPO_LOG("Starting pose optimizer, %d timesteps, %d variables\n", num_timesteps, state_size);

  // Setup the state info that will be passed to the sensors
  StateInfo state_info;
  state_info.x = state;
  state_info.size = state_size;
  state_info.tsize = block_size;
  state_info.gsize = num_global;

  // Assign toffset and goffset to each sensor, then do some checks.
  AssignSensorOffsets();
  ValidateMeasurementIndexes(num_timesteps);

  // Do Levenberg-Marquardt optimization
  double lambda = kLambdaStart;
  double error = 0;             // Error at 'state'
  bool evaled_at_state = false; // True if Hessian, gradient currently evaluated at 'state'
  int bad_counter = 0;          // Measures steps since the last bad factorization
  for (int iteration = 0; iteration < kMaxIter; iteration++) {
    // 'state' is the best known state, search for a better one.

    if (!evaled_at_state) {
      // Compute the Hessian, gradient and total error for all sensors.
      error = 0;
      H.SetZero();
      state_info.x = state;
      memset(gradient, 0, state_size * sizeof(double));
      for (int i = 0; i < kMeasurementBlocks; i++) {
        for (Sensor *s = first_sensor_; s; s = s->next_) {
          error += s->ComputeDifferential(state_info, gradient, &H, i, kMeasurementBlocks);
        }
      }
      GPO_LOG("GPO step %d: Computed H and gradient at 'state', error = %.6g\n", iteration, error);
    }

    // Add lambda times the identity to the Hessian matrix
    for (int i = 0; i < state_size; i++) {
      H.Add(i, i, lambda);
    }

    // Compute the LM update by solving trial_state = state - inv(H)*gradient
    const char *solve_error = H.Solve(gradient);
    if (solve_error) {
      // If we couldn't do it, reject this step
      if (verbosity_ >= 1) {
        GPO_LOG("GPO step %d: lambda=%.2e, factor and solve failed (%s)\n",
                  iteration, lambda, solve_error);
      }
      lambda *= kLambdaBigScale;
      evaled_at_state = false;
      bad_counter = 4;
      continue;
    }
    bad_counter--;
    GPO_LOG("GPO step %d: factored H+lambda*I successfully\n", iteration);

    // Compute the trial state
    for (int i = 0; i < state_size; i++) trial_state[i] = state[i] - gradient[i];

    // Compute the Hessian, gradient and total error at the trial state.
    double trial_error = 0;
    H.SetZero();
    state_info.x = trial_state;
    memset(gradient, 0, state_size * sizeof(double));
    for (int i = 0; i < kMeasurementBlocks; i++) {
      for (Sensor *s = first_sensor_; s; s = s->next_) {
        trial_error += s->ComputeDifferential(state_info, gradient, &H, i, kMeasurementBlocks);
      }
    }
    GPO_LOG("GPO step %d: Computed H and gradient at 'trial_state', error = %.6g\n", iteration, trial_error);

    // If the error is smaller then accept this step, otherwise reject it.
    // Make sure that infinite/NaN errors are not accepted.
    if (isfinite(trial_error) && trial_error <= error) {
      if (verbosity_ >= 2) {
        GPO_LOG("GPO step %d: lambda=%.2e, old error=%.6g, new error=%.6g (accepted) - "
                  "Fractional difference = %e\n",
                  iteration, lambda, error, trial_error, (error - trial_error)/trial_error);
      }

      // Update to the new state. Don't actually copy the vectors,
      // just swap the pointers.
      double *tmp = state;
      state = trial_state;
      trial_state = tmp;

      // If the fractional change in the error is small enough then we have
      // converged (hopefully).
      // @@@TODO: maybe check that x_new-x is small too?
      if ((error - trial_error)/trial_error < kScaleTol || trial_error < kAbsTol) {
        goto done;
      }
      evaled_at_state = true;
      error = trial_error;

      // Scale lambda for next time, but only if we have not recently factored
      // a bad matrix.
      if (bad_counter <= 0) {
        lambda /= kLambdaScale;
      }

      // Don't make lambda *too* small, to protect against rank deficiency of H.
      // @@@TODO: Do we really need this??? - FP error catching should help here.
      if (lambda < kMinLambda) lambda = kMinLambda;
    } else {
      // The error actually increased, so don't accept this step
      if (verbosity_ >= 2) {
        GPO_LOG("GPO step %d: lambda=%.2e, new error=%.6g (rejected, > %.6g)\n",
                  iteration, lambda, trial_error, error);
      }
      lambda *= kLambdaBigScale;
      evaled_at_state = false;

      // If lambda goes too large then inv(A+lambda*I) will be effectively zero,
      // so no more progress can be made. We regard this as failure.
      if (lambda > kMaxLambda) {
        return_value = "The optimizer is not making progress, "
                       "lambda is too large";
        goto done;
      }
    }
  }

  return_value = "The optimizer did not converge";

done:
  if (state != original_state) {
    memcpy(original_state, state, state_size * sizeof(double));
    delete[] state;
  } else {
    delete[] trial_state;
  }
  delete[] gradient;
  return return_value;
}

void PoseOptimizer::GetStateVectorInfo(int *num_timesteps, int *block_size,
                                       int *global_size) const {
  CHECK(first_sensor_);                 // Make sure some sensors are attached
  *block_size = POSE_BLOCK_SIZE;        // State variables per timestep
  *global_size = 0;                     // Global params at end of state vector
  *num_timesteps = 0;
  int last_step = 0;
  for (Sensor *s = first_sensor_; s; s = s->next_) {
    SensorInfo sensor_info;
    memset(&sensor_info, 0, sizeof(SensorInfo));
    s->GetInfo(&sensor_info);
    *block_size += sensor_info.num_states;
    *global_size += sensor_info.num_global_states;
    int last = s->MeasurementIndex(sensor_info.num_measurements - 1) + 1;
    if (last > last_step) last_step = last;
  }
  (*num_timesteps) = last_step + 1;
}

void PoseOptimizer::AssignSensorOffsets() const {
  // Collect state vector size information
  int num_timesteps, block_size, num_global;
  GetStateVectorInfo(&num_timesteps, &block_size, &num_global);

  int toffset = POSE_BLOCK_SIZE;
  int goffset = 0;
  for (Sensor *s = first_sensor_; s; s = s->next_) {
    s->toffset_ = toffset;
    s->goffset_ = goffset;
    SensorInfo sensor_info;
    memset(&sensor_info, 0, sizeof(SensorInfo));
    s->GetInfo(&sensor_info);
    toffset += sensor_info.num_states;
    goffset += sensor_info.num_global_states;
  }
}

void PoseOptimizer::ValidateMeasurementIndexes(int num_timesteps) const {
  for (Sensor *s = first_sensor_; s; s = s->next_) {
    SensorInfo sensor_info;
    memset(&sensor_info, 0, sizeof(SensorInfo));
    s->GetInfo(&sensor_info);
    for (int i = 0; i < sensor_info.num_measurements; i++) {
      int mi = s->MeasurementIndex(i);
      CHECK(mi > 0);
      CHECK(mi < num_timesteps-1);
    }
  }
}
