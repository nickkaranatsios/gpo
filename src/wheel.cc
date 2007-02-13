// Copyright 2007 Google Inc.
//
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

#include "gpo/wheel.h"
#include "gpo/pmath.h"
#include "gpo/rotation.h"
#include "gpo/auto_diff.h"

using namespace GPO;

WheelEncoder::WheelEncoder(int num_samples, Data *samples, double stepsize, double noise_scale) {
  num_samples_ = num_samples;
  samples_ = samples;
  stepsize_ = stepsize;
  noise_scale_ = noise_scale;
  min_noise_ = 0.03;
  pvec_[0] = 0;
  pvec_[1] = 0;
  pvec_[2] = 0;
  wvec_[0] = 0;
  wvec_[1] = 0;
  wvec_[2] = 0;
}

WheelEncoder::~WheelEncoder() {
}

void WheelEncoder::SetVectors(Vector3_d &axel, Vector3_d &road_normal, Vector3_d &wheel) {
  Vector3_d ell, n;
  // Normalize the axel and normal vectors
  double length = sqrt(sqr(axel[0]) + sqr(axel[1]) + sqr(axel[2]));
  for (int i = 0; i < 3; i++) ell[i] = axel[i] / length;
  length = sqrt(sqr(road_normal[0]) + sqr(road_normal[1]) + sqr(road_normal[2]));
  for (int i = 0; i < 3; i++) n[i] = road_normal[i] / length;
  // Compute the travel vector
  FAST_CROSS(pvec_, = , n, ell);
  // Compute the angular velocity dependent vector
  double k = FAST_DOT(wheel, ell);
  for (int i = 0; i < 3; i++) wvec_[i] = k * n[i];
}

void WheelEncoder::GetInfo(SensorInfo *info) {
  info->num_global_states = 1;          // Wheel 'smallness' (1/radius)
  info->num_measurements = num_samples_;
  info->error_size = 1;
}

int WheelEncoder::MeasurementIndex(int i) {
  return samples_[i].index;
}

void WheelEncoder::ComputeErrorVector(int i, const Number *prev_step,
  const Number *this_step, const Number *next_step, const Number *globals,
  Number *error)
{
  Number dot_p[3], dot_e[3], w[3], R[9], pvec_world[3], wvec_world[3];
  double noise = noise_scale_ * samples_[i].wheel_velocity;
  if (noise < min_noise_) noise = min_noise_;
  // Compute a central difference approximation to dp/dt
  for (int j = 0; j < 3; j++) {
    dot_p[j] = (next_step[j] - prev_step[j]) / (2.0 * stepsize_);
  }
  // Compute a central difference approximation to de/dt
  for (int j = 0; j < 3; j++) {
    dot_e[j] = (next_step[STATE_E1 + j] - prev_step[STATE_E1 + j]) / (2.0 * stepsize_);
  }
  // Compute an approximation to the angular velocity.
  eulerEv(this_step + STATE_E1, dot_e, w);
  // Compute the pvec and wvec vectors in the world frame
  eulerR(this_step + STATE_E1, R);
  FAST_MULTIPLY0_331(pvec_world, = , R, pvec_);
  FAST_MULTIPLY0_331(wvec_world, = , R, wvec_);
  // Compute the turnrate error
  const Number *smallness = globals + GOffset();
  error[0] = (*smallness) * (FAST_DOT(pvec_world, dot_p) + FAST_DOT(wvec_world, w));
  error[0] -= samples_[i].wheel_velocity;
  error[0] /= noise;
}

int WheelEncoder::GetWheelSmallnessOffset() const {
  return GOffset();
}
