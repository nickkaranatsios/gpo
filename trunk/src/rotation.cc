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

#include "gpo/rotation.h"
#include "gpo/pmath.h"
#include "gpo/auto_diff.h"

namespace GPO {

//***************************************************************************
// Utility

namespace {

// Return sin(x)/x. This has a singularity at 0 so special handling is needed
// for small arguments.
inline double sinc(double x)
{
  // If |x| < 1e-4 then use a taylor series expansion. This two term expansion
  // is actually accurate to one LS bit within this range if double precision
  // is being used - so don't worry!
  if (fabs(x) < 1.0e-4) {
    return 1.0 - x*x*0.166666666666666666667;
  } else {
    return sin(x)/x;
  }
}

// Return asin(x)/x. This has a singularity at 0 so special handling is needed
// for small arguments.
inline double asinc(double x) {
  if (fabs(x) < 1e-3) {
    double x2 = x*x;
    return 1 + x2*(1.0/6.0 + (3.0/40.0)*x2);
  } else {
    return asin(x)/x;
  }
}

}

//***************************************************************************
// Quaternions

void QtoR(const Quaternion_d &q, Matrix3x3_d &R) {
  double qq1 = 2*q[1]*q[1];
  double qq2 = 2*q[2]*q[2];
  double qq3 = 2*q[3]*q[3];
  R[0] = 1 - qq2 - qq3;
  R[1] = 2*(q[1]*q[2] - q[0]*q[3]);
  R[2] = 2*(q[1]*q[3] + q[0]*q[2]);
  R[3] = 2*(q[1]*q[2] + q[0]*q[3]);
  R[4] = 1 - qq1 - qq3;
  R[5] = 2*(q[2]*q[3] - q[0]*q[1]);
  R[6] = 2*(q[1]*q[3] - q[0]*q[2]);
  R[7] = 2*(q[2]*q[3] + q[0]*q[1]);
  R[8] = 1 - qq1 - qq2;
}

void RtoQ(const Matrix3x3_d &R, Quaternion_d &q) {
  double tr, s;
  tr = R[0] + R[4] + R[8];
  if (tr >= 0) {
    s = sqrt(tr + 1);
    q[0] = 0.5 * s;
    s = 0.5 / s;
    q[1] = (R[7] - R[5]) * s;
    q[2] = (R[2] - R[6]) * s;
    q[3] = (R[3] - R[1]) * s;
  } else {
    // Find the largest diagonal element and jump to the appropriate case
    if (R[4] > R[0]) {
      if (R[8] > R[4]) goto case_2;
      goto case_1;
    }
    if (R[8] > R[0]) goto case_2;
    goto case_0;

  case_0:
    s = sqrt((R[0] - (R[4] + R[8])) + 1);
    q[1] = 0.5 * s;
    s = 0.5 / s;
    q[2] = (R[1] + R[3]) * s;
    q[3] = (R[6] + R[2]) * s;
    q[0] = (R[7] - R[5]) * s;
    return;

  case_1:
    s = sqrt((R[4] - (R[8] + R[0])) + 1);
    q[2] = 0.5 * s;
    s = 0.5 / s;
    q[3] = (R[5] + R[7]) * s;
    q[1] = (R[1] + R[3]) * s;
    q[0] = (R[2] - R[6]) * s;
    return;

  case_2:
    s = sqrt((R[8] - (R[0] + R[4])) + 1);
    q[3] = 0.5 * s;
    s = 0.5 / s;
    q[1] = (R[6] + R[2]) * s;
    q[2] = (R[5] + R[7]) * s;
    q[0] = (R[3] - R[1]) * s;
    return;
  }
}

//***************************************************************************
// Euler angles

void ConvertQuaternionToEuler(Quaternion_d &q, Vector3_d &e) {
  // Convert the quaternion to a rotation matrix then extract the euler angles.
  Matrix3x3_d R;
  QtoR(q, R);
  e[0] = asin(R[7]);           // Pitch (+/- 90 degrees)
  e[1] = atan2(-R[6], R[8]);
  e[2] = atan2(R[1], R[4]);
}

void ConvertEulerToQuaternion(Vector3_d &e, Quaternion_d &q) {
  // Convert the euler angles to a rotation matrix, then convert the rotation
  // matrix to a quaternion.
  Matrix3x3_d R;
  eulerR(e, R);
  RtoQ(R, q);
}

void eulerR(const Vector3_d &e, Matrix3x3_d &R) {
  double s1 = sin(e[0]);
  double s2 = sin(e[1]);
  double s3 = sin(e[2]);
  double c1 = cos(e[0]);
  double c2 = cos(e[1]);
  double c3 = cos(e[2]);
  FAST_SET_MATRIX3(R,
    c3*c2+s3*s1*s2, s3*c1, c3*s2-s3*s1*c2,
    -s3*c2+c3*s1*s2, c3*c1, -s3*s2-c3*s1*c2,
    -c1*s2, s1, c1*c2);
}


void eulerR(const Number *e, Number *R) {
  //@@@TODO: Combine this with eulerR above (templated?)
  Number s1 = sin(e[0]);
  Number s2 = sin(e[1]);
  Number s3 = sin(e[2]);
  Number c1 = cos(e[0]);
  Number c2 = cos(e[1]);
  Number c3 = cos(e[2]);
  Number s3c2 = s3*c2;
  Number s1c3 = s1*c3;
  Number s2s3 = s2*s3;
  FAST_SET_MATRIX3(R,
    c3*c2+s1*s2s3, s3*c1, c3*s2-s3c2*s1,
    -s3c2+s1c3*s2, c3*c1, -s2s3-s1c3*c2,
    -c1*s2, s1, c1*c2);
}


void eulerE(const Vector3_d &e, Matrix3x3_d &E) {
  double s1 = sin(e[0]);
  double s3 = sin(e[2]);
  double c1 = cos(e[0]);
  double c3 = cos(e[2]);
  FAST_SET_MATRIX3(E, c3, s3*c1, 0, -s3, c3*c1, 0, 0, s1, -1);
}


void eulerE(const Number *e, Number *E) {
  //@@@TODO: Combine this with eulerE above.
  Number s1 = sin(e[0]);
  Number s3 = sin(e[2]);
  Number c1 = cos(e[0]);
  Number c3 = cos(e[2]);
  FAST_SET_MATRIX3(E, c3, s3*c1, 0, -s3, c3*c1, 0, 0, s1, -1);
}


void eulerEv(const Number *e, Number *v, Number *ret) {
  Number s1 = sin(e[0]);
  Number s3 = sin(e[2]);
  Number c1 = cos(e[0]);
  Number c3 = cos(e[2]);
  ret[0] = v[0]*c3 + v[1]*s3*c1;
  ret[1] = -v[0]*s3 + v[1]*c3*c1;
  ret[2] = v[1]*s1 - v[2];
}


void eulerRTE(const Vector3_d &e, Matrix3x3_d &RTE) {
  double s1 = sin(e[0]);
  double s2 = sin(e[1]);
  double c1 = cos(e[0]);
  double c2 = cos(e[1]);
  FAST_SET_MATRIX3(RTE, c2, 0, c1*s2, 0, 1, -s1, s2, 0, -c1*c2);
}


void eulerRTE(const Number *e, Number *RTE) {
  Number s1 = sin(e[0]);
  Number s2 = sin(e[1]);
  Number c1 = cos(e[0]);
  Number c2 = cos(e[1]);
  FAST_SET_MATRIX3(RTE, c2, 0, c1*s2, 0, 1, -s1, s2, 0, -c1*c2);
}


// R^T*E*vector
void eulerRTEv(const Number *e, Number *v, Number *ret) {
  Number s1 = sin(e[0]);
  Number s2 = sin(e[1]);
  Number c1 = cos(e[0]);
  Number c2 = cos(e[1]);
  ret[0] = v[0]*c2 + v[2]*c1*s2;
  ret[1] = v[1] - v[2]*s1;
  ret[2] = v[0]*s2 - v[2]*c1*c2;
}

//***************************************************************************
// Exponential maps

void WtoQ(const Vector3_d &w, Quaternion_d &q) {
  double a = 0.5 * sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
  double k = 0.5 * sinc(a);
  q[0] = cos(a);
  q[1] = k * w[0];
  q[2] = k * w[1];
  q[3] = k * w[2];
}

void QtoW(const Quaternion_d &q, Vector3_d &w) {
  double k;
  double len = sqrt(q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
  k = 2.0 * asinc(len);

  // Quaternions have double cover on the space of rotations (q and -q give
  // the same rotation), so make sure that we use the q or -q that gives
  // rotation angles in the range -pi..pi.
  if (q[0] < 0) k = -k;

  w[0] = q[1]*k;
  w[1] = q[2]*k;
  w[2] = q[3]*k;
}

void WtoR(const Vector3_d &w, Matrix3x3_d &R) {
  // Convert w into a quaternion and then convert that into a rotation
  // matrix. An alternative way is to compute:
  //    R = I + sinc(a)*W + 0.5*(sinc(a/2)*W)^2
  // where W = crossmat(w). This is about as much computation and is a
  // little bit more complicated.
  Quaternion_d q;
  WtoQ(w, q);
  QtoR(q, R);
}

void RtoW(const Matrix3x3_d &R, Vector3_d &w) {
  // The following approach works almost everywhere, but fails near the
  // singularities at angle = N*pi :
  //
  //    double trace = R[0] + R[4] + R[8];
  //    double angle = acos((trace-1.0)*0.5);
  //    double k = 0.5/sinc(angle);
  //    w[0] = (R[7]-R[5])*k;
  //    w[1] = (R[2]-R[6])*k;syms
  //    w[2] = (R[3]-R[1])*k;

  // Instead we use a more robust approach: convert R into a quaternion
  // and that convert that into exponential map coordinates.
  Quaternion_d q;
  RtoQ(R, q);
  QtoW(q, w);
}

}; // namespace GPO
