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

#include "gpo/common.h"
#include "gpo/test_matrix.h"

using namespace GPO;


TestMatrix::TestMatrix() : n_(0), m_(0), data_(0) {
}


TestMatrix::TestMatrix(int rows, int cols) {
  n_ = 0;
  m_ = 0;
  data_ = 0;
  SetSize(rows, cols);
}


TestMatrix::TestMatrix(const TestMatrix &a) {
  n_ = a.n_;
  m_ = a.m_;
  data_ = new double[n_*m_];
  memcpy(data_, a.data_, n_*m_*sizeof(double));
}


TestMatrix::TestMatrix(int rows, int cols,
                        double *data, int rowskip, int colskip) {
  CHECK(rows >= 1 && cols >= 1);
  n_ = rows;
  m_ = cols;
  data_ = new double[n_*m_];
  for (int i = 0; i < n_; i++) {
    for (int j = 0; j < m_; j++) data_[i*m_+j] = data[i*rowskip + j*colskip];
  }
}


TestMatrix::~TestMatrix() {
  if (data_) delete[] data_;
}


void TestMatrix::SetSize(int rows, int cols) {
  if (data_) delete[] data_;
  CHECK(rows >= 1 && cols >= 1);
  n_ = rows;
  m_ = cols;
  data_ = new double[n_*m_];
  for (int i = 0; i < n_*m_; i++) data_[i] = 0;
}


double & TestMatrix::operator() (int i, int j) {
  CHECK(i >= 0 && i < n_ && j >= 0 && j < m_);
  return data_[i*m_+j];
}


void TestMatrix::operator = (const TestMatrix &a) {
  if (data_) delete[] data_;
  n_ = a.n_;
  m_ = a.m_;
  data_ = new double[n_*m_];
  memcpy(data_, a.data_, n_*m_*sizeof(double));
}


void TestMatrix::operator = (double a) {
  for (int i = 0; i < n_*m_; i++) data_[i] = a;
}


TestMatrix TestMatrix::Transpose() {
  TestMatrix r(m_, n_);
  for (int i = 0; i < n_; i++) {
    for (int j = 0; j < m_; j++) r.data_[j*n_+i] = data_[i*m_+j];
  }
  return r;
}


TestMatrix TestMatrix::Select(int np, int *p, int nq, int *q) {
  CHECK(np >= 1 && nq >= 1);
  TestMatrix r(np, nq);
  for (int i = 0; i < np; i++) {
    for (int j = 0; j < nq; j++) {
      CHECK(p[i] >= 0 && p[i] < n_ && q[i] >= 0 && q[i] < m_);
      r.data_[i*nq+j] = data_[p[i]*m_+q[j]];
    }
  }
  return r;
}


TestMatrix TestMatrix::operator + (const TestMatrix &a) {
  CHECK(n_ == a.n_ && m_ == a.m_);
  TestMatrix r(n_, m_);
  for (int i = 0; i < n_*m_; i++) r.data_[i] = data_[i] + a.data_[i];
  return r;
}


TestMatrix TestMatrix::operator - (const TestMatrix &a) {
  CHECK(n_ == a.n_ && m_ == a.m_);
  TestMatrix r(n_, m_);
  for (int i = 0; i < n_*m_; i++) r.data_[i] = data_[i] - a.data_[i];
  return r;
}


TestMatrix TestMatrix::operator - () {
  TestMatrix r(n_, m_);
  for (int i = 0; i < n_*m_; i++) r.data_[i] = -data_[i];
  return r;
}


TestMatrix TestMatrix::operator * (const TestMatrix &a) {
  CHECK(m_ == a.n_);
  TestMatrix r(n_, a.m_);
  for (int i = 0; i < n_; i++) {
    for (int j = 0; j < a.m_; j++) {
      double sum = 0;
      for (int k = 0; k < m_; k++) sum += data_[i*m_+k] * a.data_[k*a.m_+j];
      r.data_[i*a.m_+j] = sum;
    }
  }
  return r;
}


void TestMatrix::operator += (const TestMatrix &a) {
  CHECK(n_ == a.n_ && m_ == a.m_);
  for (int i = 0; i < n_*m_; i++) data_[i] += a.data_[i];
}


void TestMatrix::operator -= (const TestMatrix &a) {
  CHECK(n_ == a.n_ && m_ == a.m_);
  for (int i = 0; i < n_*m_; i++) data_[i] -= a.data_[i];
}


void TestMatrix::ClearUpperTriangle() {
  CHECK(n_ == m_);
  for (int i = 0; i < n_; i++) {
    for (int j = i+1; j < m_; j++) data_[i*m_+j] = 0;
  }
}


void TestMatrix::ClearLowerTriangle() {
  CHECK(n_ == m_);
  for (int i = 0; i < n_; i++) {
    for (int j = 0; j < i; j++) data_[i*m_+j] = 0;
  }
}


void TestMatrix::MakeRandom(double range) {
  for (int i = 0; i < n_; i++) {
    for (int j = 0; j < m_; j++)
      data_[i*m_+j] = (double(random())/double(RAND_MAX)*2.0 - 1.0)*range;
  }
}


void TestMatrix::Print(char *fmt, FILE *f) {
  for (int i = 0; i < n_; i++) {
    for (int j = 0; j < m_; j++) fprintf(f, fmt, data_[i*m_+j]);
    fprintf(f, "\n");
  }
}


void TestMatrix::PrintNonzeros(double tolerance, int blocksize, FILE *f) {
  for (int i = 0; i < n_; i++) {
    for (int j = 0; j < m_; j++) {
      fputc((data_[i*m_+j] > tolerance) ? 'X' : '.', f);
      if (blocksize > 0 && ((j+1) % blocksize) == 0) fputc(' ', f);
    }
    fputc('\n', f);
    if (blocksize > 0 && ((i+1) % blocksize) == 0) fputc('\n', f);
  }
}


double TestMatrix::MaxDifference(const TestMatrix &M) {
  CHECK(n_ == M.n_ && m_ == M.m_);
  double max = 0;
  for (int i = 0; i < n_; i++) {
    for (int j = 0; j < m_; j++) {
      double diff = fabs(data_[i*m_+j] - M.data_[i*m_+j]);
      if (diff > max) max = diff;
    }
  }
  return max;
}

double TestMatrix::Distance(const TestMatrix &M) {
  CHECK(n_ == M.n_ && m_ == M.m_);
  double sum = 0;
  for (int i = 0; i < n_; i++) {
    for (int j = 0; j < m_; j++) {
      double diff = data_[i*m_+j] - M.data_[i*m_+j];
      sum += diff*diff;
    }
  }
  return sqrt(sum);
}
