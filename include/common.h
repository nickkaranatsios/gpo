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

// Common things.

#ifndef __GPO_COMMON_H__
#define __GPO_COMMON_H__

#include "gpo/config.h"

namespace GPO {

// Doubles should be aligned to 8 byte boundaries to be fast on the x86, but
// the gcc ABI only aligns them to 4. Enforce the alignment by using the
// following type.
#ifdef __GNUC__
  typedef double __attribute__((aligned(8))) AlignedDouble;
#else
  typedef double AlignedDouble;
#endif

// If gcc __PRETTY_FUNCTION__ not available, make a static alternative.
#ifndef __GNUC__
  #define __PRETTY_FUNCTION__ "<function>"
#endif

// An error handling function - print the error message and die.
void gpo_panic(const char *msg, ...)
  #ifdef __GNUC__
    __attribute__((format(printf, 1, 2))) __attribute__((noreturn))
  #endif
  ;

// CHECK() is like assert() except that it is always compiled in even
// when building optimized. DCHECK() is like assert().
#define CHECK(cond) { \
  if (!(cond)) { \
    GPO::gpo_panic("In %s, assertion '%s' failed (at %s:%d)", \
           __PRETTY_FUNCTION__, #cond, __FILE__, __LINE__); \
  } \
}
#define DCHECK(cond) CHECK(cond)

// A logging macro, used in various places to report status
#define GPO_LOG(msg, args...) fprintf(stderr, msg, ##args);

// A macro to report a fatal error and die
#define GPO_PANIC(msg, args...) gpo_panic(msg, ##args);

// Square the number.
extern inline double sqr(double x) {
  return x*x;
}

typedef double Vector3_d[3];
typedef double Quaternion_d[4];
typedef double Matrix3x3_d[9];

}; // namespace GPO

#endif
