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

// Utility functions for rotations (Euler angles, quaternions and exponential
// maps).

#ifndef __GPO_ROTATION_H__
#define __GPO_ROTATION_H__

#include "gpo/common.h"

namespace GPO {

class Number;


// ********** Quaternions

// Convert a quaternion to a rotation matrix.
void RtoQ(const Matrix3x3_d &R, Quaternion_d &q);

// Convert a rotation matrix into a quaternion.
void QtoR(const Quaternion_d &q, Matrix3x3_d &R);

// ********** Euler angles

// @@@TODO: Name these functions consistently. The QtoR() etc functions should
// probably have long names like ConvertQuaternionToEuler().

// Convert a quaternion to euler angles
void ConvertQuaternionToEuler(Quaternion_d &q, Vector3_d &e);

// Convert euler angles to a quaternion
void ConvertEulerToQuaternion(Vector3_d &e, Quaternion_d &q);

// Convert euler angles in 'e' to a rotation matrix 'R'.
void eulerR(const Vector3_d &e, Matrix3x3_d &R);
void eulerR(const Number *e, Number *R);

  // Compute the Euler derivative matrix from euler angle parameters
void eulerE(const Vector3_d &e, Matrix3x3_d &E);
void eulerE(const Number *e, Number *E);
void eulerEv(const Number *e, Number *v, Number *ret);

// Compute R'*E from the euler angles e
void eulerRTE(const Vector3_d &e, Matrix3x3_d &RTE);
void eulerRTE(const Number *e, Number *RTE);
void eulerRTEv(const Number *e, Number *v, Number *ret);

// ********** Exponential maps

// Convert the exponential map coordinates w into a quaternion q.
void WtoQ(const Vector3_d &w, Quaternion_d &q);

// Convert the quaternion q into exponential map coordinates w.
void QtoW(const Quaternion_d &q, Vector3_d &w);

// Compute R = expm(hat(w)), i.e. convert the exponential map
// coordinates w into a rotation matrix.
void WtoR(const Vector3_d &w, Matrix3x3_d &R);

// Compute w = nohat(logm(R)), i.e. convert the rotation matrix
// into exponential map coordinates w. Note that there are many
// exponential map vectors that map to a single R, this returns
// the one inside the sphere of radius pi.
void RtoW(const Matrix3x3_d &R, Vector3_d &w);

}; // namespace GPO

#endif
