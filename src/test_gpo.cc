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

// Unit test for everything

#include "gpo/gpo.h"
#include "gpo/common.h"
#include "gpo/rotation.h"
#include "gpo/hessian.h"
#include "gpo/fast_block.h"
#include "gpo/pmath.h"
#include "gpo/auto_diff.h"

using namespace GPO;

// Anonymous namespace for globals defined below
namespace {

//***************************************************************************
// Utility

double rand_real() {
  return double(random()) / double(RAND_MAX);
}

void print_error(const char *name, double error) {
  printf("\t%-12s error = %.4e (%s)\n", name, error,
    ((error > 1e-6) || isinf(error) || isnan(error)) ? "*** BAD ***" : "ok");
}

double get_time() {
  struct timeval tv;
  gettimeofday(&tv, 0);
  return tv.tv_sec + tv.tv_usec / 1e6;
}

//***************************************************************************
// Test Euler angle conversion functions

void test_euler_angle_conversion() {
  printf("Running %s()\n", __FUNCTION__);

  // Convert euler angles to quaternions and back
  double max_error = 0;
  for (int i = 0; i < 100000; i++) {
    Vector3_d e, e2;
    Quaternion_d q;
    e[0] = (M_PI * 0.98) * (rand_real() - 0.5);   // Pitch - avoid singularity
    e[1] = 2.0 * M_PI * (rand_real() - 0.5);
    e[2] = 2.0 * M_PI * (rand_real() - 0.5);

    ConvertEulerToQuaternion(e, q);
    ConvertQuaternionToEuler(q, e2);

    for (int j = 0; j < 3; j++) {
      double E = fabs(e[j] - e2[j]);
      if (E > max_error) max_error = E;
    }
  }
  print_error("E -> Q -> E", max_error);

  // Convert quaternions to euler angles and back
  max_error = 0;
  for (int i = 0; i < 100000; i++) {
    double len = 0;
    Vector3_d e;
    Quaternion_d q, q2;
    for (int j = 0; j < 4; j++) {
      q[j] = rand_real() - 0.5;
      len += q[j]*q[j];
    }
    for (int j = 0; j < 4; j++) q[j] /= sqrt(len);

    ConvertQuaternionToEuler(q, e);
    ConvertEulerToQuaternion(e, q2);

    // Make sure q and q2 have the same sign, since there are two quaternion
    // representations for one rotation.
    if (q[0]*q2[0] < 0) {
      for (int j = 0; j < 4; j++) q2[j] *= -1;
    }

    for (int j = 0; j < 4; j++) {
      double E = fabs(q[j] - q2[j]);
      if (E > max_error) max_error = E;
    }
  }
  print_error("Q -> E -> Q", max_error);

}

//***************************************************************************
// Test exponential map functions

// Check that a rotation matrix has an invertible representations as
// exponential map coordinates.

double test_R_to_W_to_R(Matrix3x3_d &R) {
  // Convert R to an exponential map and back
  Vector3_d w;
  RtoW(R, w);
  Matrix3x3_d R2;
  WtoR(w, R2);

  // Check the conversion
  double E = 0;
  for (int i = 0; i < 9; i++) E += fabs(R[i] - R2[i]);
  return E;
}

// Run the conversion tests

void test_exponential_map_conversions() {
  printf("Running %s()\n", __FUNCTION__);

  // Check that WtoR() makes real rotation matrices
  double E = 0;
  for (int iteration = 0; iteration < 100000; iteration++) {
    // Generate random exponential map coordinates
    Vector3_d w;
    for (int i = 0; i < 3; i++) w[i] = (rand_real()-0.5) * 10.0;

    // Convert it to a rotation matrix
    Matrix3x3_d R;
    WtoR(w, R);

    // Check that the rotation matrix is orthonormal: R*R'=I
    Matrix3x3_d I;
    FAST_MULTIPLY1_333(I, = , R, R);
    I[0] -= 1.0;
    I[4] -= 1.0;
    I[8] -= 1.0;
    for (int i = 0; i < 9; i++) E += fabs(I[i]);
  }
  print_error("R*R' = I?", E);

  // Check that RtoW() is the inverse of WtoR()
  double max_error = 0;
  for (int iteration = 0; iteration < 100000; iteration++) {
    // Generate random exponential map coordinates in the sphere of radius pi
    Vector3_d w;
    double length;
    do {
      for (int i = 0; i < 3; i++) w[i] = (rand_real()-0.5) * 2.0 * M_PI;
      length = sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
    } while (length >= M_PI-0.01);

    // Convert it to a rotation matrix and back
    Matrix3x3_d R;
    WtoR(w, R);
    Vector3_d w2;
    RtoW(R, w2);

    // Check the conversion
    for (int i = 0; i < 3; i++) {
      double E = fabs(w[i] - w2[i]);
      if (E > max_error) max_error = E;
    }
  }
  print_error("w -> R -> w", max_error);

  // Check that some random rotation matrices have invertible representations
  // as exponential maps.
  E = 0;
  for (int iteration = 0; iteration < 100000; iteration++) {
    // Generate random rotation matrices from random quaternions.
    Quaternion_d q;
    for (int i = 0; i < 4; i++) q[i] = (rand_real()-0.5);
    double scale = 1.0 / sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
    for (int i = 0; i < 4; i++) q[i] *= scale;
    Matrix3x3_d R;
    QtoR(q, R);

    // Convert it to an exponential map and back.
    E += test_R_to_W_to_R(R);
  }
  print_error("R -> w -> R", E);

  // Check that some troublesome 180-degree rotation matrices have
  // invertible representations as exponential maps.
  E = 0;
  Matrix3x3_d R;
  FAST_SET_MATRIX3(R, -1, 0, 0, 0, -1, 0, 0, 0, 1);
  E += test_R_to_W_to_R(R);
  FAST_SET_MATRIX3(R, 1, 0, 0, 0, -1, 0, 0, 0, -1);
  E += test_R_to_W_to_R(R);
  FAST_SET_MATRIX3(R, -1, 0, 0, 0, 1, 0, 0, 0, -1);
  E += test_R_to_W_to_R(R);
  print_error("Hard R->w->R", E);
}

//***************************************************************************
// Test fast_block functions

void test_solve() {
  printf("Running %s()\n", __FUNCTION__);
  double max_error = 0;
  for (int n = 1; n <= 30; n++) {
    for (int m = 1; m <= 30; m++) {
      // Create random matrices
      int sz = (n*(n+1))/2;
      double L[sz], Bref[n*m], Bopt[n*m];
      for (int i = 0; i < sz; i++) L[i] = rand_real()*2.0 - 1.0;
      for (int i = 0; i < n; i++) L[((i+1)*(i+2))/2-1] = rand_real() + 1.0;
      for (int i = 0; i < n*m; i++) Bref[i] = Bopt[i] = rand_real()*2.0 - 1.0;

      // Solve with the reference and optimized code
      solve_ref(L, Bref, n, m);
      solve(L, Bopt, n, m);

      // Find the maximum per-element error
      for (int i = 0; i < n*m; i++) {
        double E = fabs(Bref[i] - Bopt[i]);
        if (E > max_error) max_error = E;
      }
    }
  }

  print_error("solve", max_error);
}

void test_multiply2() {
  printf("Running %s()\n", __FUNCTION__);
  double error = 0;
  for (int p = 1; p <= 15; p++) {
    for (int q = 1; q <= 15; q++) {
      for (int r = 1; r <= 15; r++) {
        // Create random matrices
        double B[p*q], C[q*r];
        for (int i = 0; i < p*q; i++) B[i] = rand_real()*2.0 - 1.0;
        for (int i = 0; i < q*r; i++) C[i] = rand_real()*2.0 - 1.0;

        // Multiply with the reference and optimized code
        double Aref[p*r], Aopt[p*r];
        for (int i = 0; i < p*r; i++) Aref[i] = Aopt[i] = rand_real()*2.0 - 1.0;
        multiply2_ref(Aref, 1.2, B, C, p, q, r);
        multiply2(Aopt, 1.2, B, C, p, q, r);

        // Accumulate the error
        for (int i = 0; i < p*r; i++) error += fabs(Aref[i] - Aopt[i]);
      }
    }
  }

  print_error("multiply2", error);
}

void test_multiply2_lt() {
  printf("Running %s()\n", __FUNCTION__);
  double error = 0;
  for (int p = 1; p <= 15; p++) {
    for (int q = 1; q <= 15; q++) {
      // Create random matrices
      double B[p*q], C[q*p];
      for (int i = 0; i < p*q; i++) B[i] = rand_real()*2.0 - 1.0;
      for (int i = 0; i < q*p; i++) C[i] = rand_real()*2.0 - 1.0;

      // Multiply with the reference and optimized code
      int sz = (p*(p+1))/2;
      double Aref[sz], Aopt[sz];
      for (int i = 0; i < sz; i++) Aref[i] = Aopt[i] = rand_real()*2.0 - 1.0;
      multiply2_lt_ref(Aref, 1.2, B, C, p, q);
      multiply2_lt(Aopt, 1.2, B, C, p, q);

      // Accumulate the error
      for (int i = 0; i < sz; i++) error += fabs(Aref[i] - Aopt[i]);
    }
  }

  print_error("multiply2_lt", error);
}

//***************************************************************************
// Test Hessian factorization

// Test the accuracy of the factorization for a given B,N,M.

double test_hessian_helper(int B, int N, int M) {
  int sz = B*N + M;                     // Total matrix size

  // Setup vectors
  double x[sz], b[sz], btest[sz];
  for (int i = 0; i < sz; i++) b[i] = rand_real()*2.0-1.0;

  // Setup a random hessian (diagonally dominant)
  HessianMatrix H, Horig;
  H.Initialize(B, N, M);
  Horig.Initialize(B, N, M);
  for (int i = 0; i < sz; i++) {
    double val = rand_real()+20.0;
    H.Add(i, i, val);
    Horig.Add(i, i, val);
  }
  for (int i = 0; i < sz; i++) {
    for (int j = 0; j < i; j++) {
      double *element = H.Access(i, j);
      if (element) {
        double val = rand_real()*2.0-1.0;
        *element = val;
        element = Horig.Access(i, j);
        *element = val;
      }
    }
  }

  // Factor and solve
  memcpy(x, b, sz*sizeof(double));
  const char *error_message = H.Solve(x);
  if (error_message) {
    printf ("\tSolve error message: %s (*** BAD ***)\n", error_message);
  }

  // Test solution
  Horig.Multiply(x, btest);
  double error = 0;
  for (int i = 0; i < sz; i++) error += fabs(b[i] - btest[i]);
  //printf ("B=%d N=%d M=%d, error = %e\n",B,N,M,error);
  return error;
}

// Test the accuracy of the factorization for a range of B,N,M.

void test_hessian() {
  printf("Running %s()\n", __FUNCTION__);
  double max_error = 0;
  for (int B = 2; B <= 15; B++) {
    for (int N = 1; N <= 15; N++) {
      for (int M = 1; M <= 15; M++) {
        double error = test_hessian_helper(B, N, M);
        if (error > max_error) max_error = error;
      }
    }
  }
  print_error("factor", max_error);
}

//***************************************************************************
// Test sensor computation of the measurement vectors, Jacobians and Hessians

const int kNumSamples = 50;

// Test the sensor's automatic differentiation by inspecting the Jacobian for
// each measurement. Return the full Jacobian matrix and error vector.

void test_automatic_differentiation(Sensor &sensor, SensorInfo &info,
  TestMatrix &x, TestMatrix &Jret, TestMatrix &error_ret)
{
  double max_err = 0;

  // Stepsize for numerical derivatives
  const double h = 1.0e-6;

  // Compute sensor information
  int steps = sensor.MeasurementIndex(info.num_measurements-1) + 2;
  int tsize = POSE_BLOCK_SIZE + info.num_states;
  int gsize = info.num_global_states;
  Jret.SetSize(info.error_size * info.num_measurements, steps*tsize + gsize);
  error_ret.SetSize(info.error_size * info.num_measurements, 1);

  for (int i = 0; i < info.num_measurements; i++) {
    int index = sensor.MeasurementIndex(i);

    // Set up the Numbers
    int num_deps = 3*tsize + gsize;
    Number prev_step[tsize];
    Number this_step[tsize];
    Number next_step[tsize];
    Number globals[gsize];
    Number error[info.error_size];

    // Jacobian matrix and its numerical-differentiation estimate
    TestMatrix J(info.error_size, num_deps);
    TestMatrix Jest(info.error_size, num_deps);

    // Set dependent variable values and dependencies.
    for (int k = 0; k < tsize; k++) {
      prev_step[k].SetDependent(k, x((index-1)*tsize + k, 0));
      this_step[k].SetDependent(k + tsize, x((index)*tsize + k, 0));
      next_step[k].SetDependent(k + 2*tsize, x((index+1)*tsize + k, 0));
    }
    for (int k = 0; k < gsize; k++) {
      globals[k].SetDependent(k + 3*tsize, x(steps*tsize + k, 0));
    }

    // Compute the error vector and all its partial derivatives.
    sensor.ComputeErrorVector(i, prev_step, this_step, next_step, globals, error);

    // Extract the Jacobian matrix.
    for (int k = 0; k < info.error_size; k++) {
      error[k].GetNonzeroDerivatives(num_deps, &J(k, 0));
    }

    // Put J in the right slots in Jret
    for (int k = 0; k < info.error_size; k++) {
      for (int m = 0; m < 3*tsize; m++) {
        Jret(k + i*info.error_size, m + (index-1)*tsize) = J(k, m);
      }
      for (int m = 0; m < gsize; m++) {
        Jret(k + i*info.error_size, m + steps*tsize) = J(k, m + 3*tsize);
      }
    }

    // Put error in the right slots in error_ret
    for (int k = 0; k < info.error_size; k++) {
      error_ret(k + i*info.error_size, 0) = to_double(error[k]);
    }

    // Now compute an estimate of J using numerical differentiation.
    for (int m = 0; m < tsize; m++) {
      Number error1[info.error_size], error2[info.error_size];

      Number save = prev_step[m];
      prev_step[m] -= h;
      sensor.ComputeErrorVector(i, prev_step, this_step, next_step, globals, error1);
      prev_step[m] += 2*h;
      sensor.ComputeErrorVector(i, prev_step, this_step, next_step, globals, error2);
      for (int k = 0; k < info.error_size; k++) {
        Jest(k, m) = to_double( (error2[k] - error1[k]) / (2*h) );
      }
      prev_step[m] = save;

      save = this_step[m];
      this_step[m] -= h;
      sensor.ComputeErrorVector(i, prev_step, this_step, next_step, globals, error1);
      this_step[m] += 2*h;
      sensor.ComputeErrorVector(i, prev_step, this_step, next_step, globals, error2);
      for (int k = 0; k < info.error_size; k++) {
        Jest(k, m + tsize) = to_double( (error2[k] - error1[k]) / (2*h) );
      }
      this_step[m] = save;

      save = next_step[m];
      next_step[m] -= h;
      sensor.ComputeErrorVector(i, prev_step, this_step, next_step, globals, error1);
      next_step[m] += 2*h;
      sensor.ComputeErrorVector(i, prev_step, this_step, next_step, globals, error2);
      for (int k = 0; k < info.error_size; k++) {
        Jest(k, m + 2*tsize) = to_double( (error2[k] - error1[k]) / (2*h) );
      }
      next_step[m] = save;
    }
    for (int m = 0; m < gsize; m++) {
      Number error1[info.error_size], error2[info.error_size];

      Number save = globals[m];
      globals[m] -= h;
      sensor.ComputeErrorVector(i, prev_step, this_step, next_step, globals, error1);
      globals[m] += 2*h;
      sensor.ComputeErrorVector(i, prev_step, this_step, next_step, globals, error2);
      for (int k = 0; k < info.error_size; k++) {
        Jest(k, m + 3*tsize) = to_double( (error2[k] - error1[k]) / (2*h) );
      }
      globals[m] = save;
    }

    double err = J.Distance(Jest);
    if (err > max_err) max_err = err;
  }
  print_error("Jacobian", max_err);
}

void test_sensor(Sensor &sensor) {
  // Get sensor information
  SensorInfo sensor_info;
  memset(&sensor_info, 0, sizeof(sensor_info));
  sensor.GetInfo(&sensor_info);
  int steps = sensor.MeasurementIndex(sensor_info.num_measurements-1) + 2;

  // Make a random state vector
  TestMatrix x(steps*(POSE_BLOCK_SIZE + sensor_info.num_states) +
    sensor_info.num_global_states, 1);
  x.MakeRandom(1.0);

  // Setup the StateInfo data to pass to the sensor
  StateInfo state_info;
  state_info.x = x.Data();
  state_info.size = x.Rows();
  state_info.tsize = POSE_BLOCK_SIZE + sensor_info.num_states;
  state_info.gsize = sensor_info.num_global_states;

  // Setup the sensor
  sensor.SetTOffset(POSE_BLOCK_SIZE);
  sensor.SetGOffset(0);

  // Test automatic differentiation
  TestMatrix J, error;
  test_automatic_differentiation(sensor, sensor_info, x, J, error);

  // Compute the sensor's gradient and Hessian matrix
  HessianMatrix H2;
  H2.Initialize(steps, state_info.tsize, sensor_info.num_global_states);
  TestMatrix gradient(x.Rows(), 1);
  double scalar_error = sensor.ComputeDifferential(state_info, gradient.Data(), &H2);

  // Copy the hessian matrix into TestMatrix format
  TestMatrix H(x.Rows(), x.Rows());
  for (int i = 0; i < x.Rows(); i++) {
    for (int j = 0; j < x.Rows(); j++) {
      H(i, j) = H2.Get(i, j);
    }
  }

  // Check the error vector
  TestMatrix error_good = error.Transpose() * error;
  TestMatrix scalar_error2(1, 1);
  scalar_error2(0, 0) = scalar_error;
  print_error("Scalar error", error_good.Distance(scalar_error2));

  // Check the gradient vector
  TestMatrix gradient_good = J.Transpose() * error;
  print_error("Gradient", gradient_good.Distance(gradient));

  // Check the Hessian matrix
  TestMatrix H_good = J.Transpose() * J;
  print_error("Hessian", H_good.Distance(H));
}

void test_gps_sensor() {
  printf("Running %s()\n", __FUNCTION__);

  // Initialize sensor, zero all the measurement data but randomize the
  // covariance.
  GPSSensor::Data *samples = new GPSSensor::Data[kNumSamples-2];
  memset(samples, 0, (kNumSamples-2) * sizeof(GPSSensor::Data));
  for (int i = 0; i < kNumSamples-2; i++) {
    samples[i].index = i+1;
    samples[i].sx = rand_real() + 10.0;
    samples[i].sy = rand_real() + 10.0;
    samples[i].sz = rand_real() + 10.0;
  }
  GPSSensor sensor(kNumSamples-2, samples);

  test_sensor(sensor);
}

void test_gps_velocity_sensor() {
  printf("Running %s()\n", __FUNCTION__);

  // Initialize sensor, zero all the measurement data but randomize the
  // covariance.
  GPSVelocitySensor::Data *samples = new GPSVelocitySensor::Data[kNumSamples-2];
  memset(samples, 0, (kNumSamples-2) * sizeof(GPSVelocitySensor::Data));
  for (int i = 0; i < kNumSamples-2; i++) {
    samples[i].index = i+1;
    samples[i].sx = rand_real() + 10.0;
    samples[i].sy = rand_real() + 10.0;
    samples[i].sz = rand_real() + 10.0;
  }
  GPSVelocitySensor sensor(kNumSamples-2, samples, 0.1*rand_real() + 0.1, 0);

  test_sensor(sensor);
}

void test_accelerometer_sensor() {
  printf("Running %s()\n", __FUNCTION__);

  // Initialize sensor, zero all the measurement data but randomize the
  // covariances.
  AccelerometerSensor::Data *samples = new AccelerometerSensor::Data[kNumSamples-1];
  memset(samples, 0, (kNumSamples-1) * sizeof(AccelerometerSensor::Data));
  for (int i = 0; i < kNumSamples-1; i++) {
    samples[i].std_dev = rand_real() + 1.0;
  }
  AccelerometerSensor sensor(kNumSamples-1, samples, 0.1*rand_real() + 0.1,
    rand_real() + 1.0);
  sensor.SetGravity(rand_real()-0.5, rand_real()-0.5, rand_real()-0.5);

  test_sensor(sensor);
}

void test_rate_gyro_sensor() {
  printf("Running %s()\n", __FUNCTION__);

  // Initialize sensor, zero all the measurement data but randomize the
  // covariances.
  RateGyroSensor::Data *samples = new RateGyroSensor::Data[kNumSamples-1];
  memset(samples, 0, (kNumSamples-1) * sizeof(RateGyroSensor::Data));
  for (int i = 0; i < kNumSamples-1; i++) {
    samples[i].std_dev = rand_real() + 1.0;
  }
  RateGyroSensor sensor(kNumSamples-1, samples, 0.1*rand_real() + 0.1,
    rand_real() + 1.0);
  sensor.SetSpin(rand_real()-0.5, rand_real()-0.5, rand_real()-0.5);

  test_sensor(sensor);
}

void test_axel_sensor() {
  printf("Running %s()\n", __FUNCTION__);

  WheelEncoder::Data *samples = new WheelEncoder::Data[kNumSamples-2];
  for (int i = 0; i < kNumSamples-2; i++) {
    samples[i].index = i+1;
    samples[i].wheel_velocity = rand_real() - 0.5;
  }

  Vector3_d axel, wheel;
  for (int i = 0; i < 3; i++) axel[i] = rand_real() - 0.5;
  for (int i = 0; i < 3; i++) wheel[i] = rand_real() - 0.5;
  AxelConstraint sensor(kNumSamples-2, samples, 0.1*rand_real() + 0.1, rand_real() + 1.0);
  sensor.SetVectors(axel, wheel);
  test_sensor(sensor);
}

void test_wheel_sensor() {
  printf("Running %s()\n", __FUNCTION__);

  WheelEncoder::Data *samples = new WheelEncoder::Data[kNumSamples-2];
  for (int i = 0; i < kNumSamples-2; i++) {
    samples[i].index = i+1;
    samples[i].wheel_velocity = rand_real() - 0.5;
  }

  Vector3_d axel, normal, wheel;
  for (int i = 0; i < 3; i++) axel[i] = rand_real() - 0.5;
  for (int i = 0; i < 3; i++) normal[i] = rand_real() - 0.5;
  for (int i = 0; i < 3; i++) wheel[i] = rand_real() - 0.5;
  WheelEncoder sensor(kNumSamples-2, samples, 0.1*rand_real() + 0.1, rand_real() + 1.0);
  sensor.SetVectors(axel, normal, wheel);
  test_sensor(sensor);
}

void test_road_sensor() {
  printf("Running %s()\n", __FUNCTION__);

  WheelEncoder::Data *samples = new WheelEncoder::Data[kNumSamples-2];
  for (int i = 0; i < kNumSamples-2; i++) {
    samples[i].index = i+1;
    samples[i].wheel_velocity = rand_real() - 0.5;
  }

  Vector3_d road;
  for (int i = 0; i < 3; i++) road[i] = rand_real() - 0.5;
  RoadConstraint sensor(kNumSamples-2, samples, 0.1*rand_real() + 0.1, rand_real() + 1.0);
  sensor.SetRoadVector(road);
  test_sensor(sensor);
}

//***************************************************************************

}; // anonymous namespace


int main(int argc, char **argv) {
  srandom(time(0));

  test_euler_angle_conversion();
  test_exponential_map_conversions();

  test_solve();
  test_multiply2();
  test_multiply2_lt();

  test_gps_sensor();
  test_gps_velocity_sensor();
  test_accelerometer_sensor();
  test_rate_gyro_sensor();
  test_axel_sensor();
  test_wheel_sensor();
  test_road_sensor();

  test_hessian();

  return 0;
}
