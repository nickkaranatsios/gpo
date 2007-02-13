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

// Automatic differentiation.
// There are a number of code configuration options controlled by #defines.

//@@@TODO: This code is x86-32 specific. Port it to (a) x86-64 (b) other architectures.

#ifndef __GPO_AUTO_DIFF_H__
#define __GPO_AUTO_DIFF_H__

#include "gpo/common.h"

namespace GPO {

// Configuration:

// Define this to have up to 64 dependent variables. Undefine it to have up to
// 32 dependent variables (this is faster).
//#define INDEX_SET_64

// If this is defined each Number keeps a pointer to an external differential
// vector. If this is not defined that vector is stored within Number.
// NOTE: if this is undefined, USE_REFCOUNTING must be undefined too.
//#define EXTERNAL_VECTOR

// If EXTERNAL_VECTOR is defined, this does reference counting and copy-on-write
// of the differential vectors.
#ifdef EXTERNAL_VECTOR
//#define USE_REFCOUNTING
#endif

// Define this to use DVector memcpy instead of a sparse index set copy.
//#define USE_MEMCPY


#ifdef INDEX_SET_64
  // Support for index set size = 64

  const int IndexSetSize = 64;
  typedef uint64 IndexSet;

  // Return nonzero if index 'i' is in the set
  inline uint64 index_in_set(IndexSet set, int i) {
    return set & (1ULL << i);
  }

  // Return a mask with indexes 0..i-1 set
  inline IndexSet index_range(int i) {
    if (i >= 64) return IndexSet(-1);
    return(1ULL << i) -1LL;
  }

  // Iterate over all set indexes
  #define _FOR_ALL_INDEXES(i, mask, todo) { \
    for (uint32 __m = mask; __m; __m = __m & (__m-1)) { \
      uint32 i = lowest_bit_index_32(__m); \
      todo \
    } \
    for (uint32 __m = (mask & 0xffffffff00000000ULL) >> 32; __m; __m = __m & (__m-1)) { \
      uint32 i = 32+lowest_bit_index_32(__m); \
      todo \
    } \
  }

#else
  // Support for index set size = 32

  const int IndexSetSize = 32;
  typedef uint32 IndexSet;

  // Return nonzero if index 'i' is in the set
  inline uint32 index_in_set(IndexSet set, int i) {
    return set & (1 << i);
  }

  // Return a mask with indexes 0..i-1 set
  inline IndexSet index_range(int i) {
    if (i >= 32) return IndexSet(-1);
    return(1 << i) -1;
  }

  // Iterate over all set indexes
  #define _FOR_ALL_INDEXES(i, mask, todo) \
    for (uint32 __m = mask; __m; __m = __m & (__m-1)) { \
      uint32 i = lowest_bit_index_32(__m); \
      todo \
    }

#endif


// Two ways to get the index of the lowest set bit in a 32 bit integer.
// The BSF instruction is (on the author's P4) a bit faster, but is
// x86-specific. The table-based way is more portable.

#if GPO_PENTIUM_AVAILABLE

inline int lowest_bit_index_32(uint32 n) {
  int result;
  asm("bsf %[n],%[result]" : [result] "=r" (result) : [n] "rm" (n) : "cc" );
  return result;
}

#else

inline int lowest_bit_index_32(uint32 n) {
  static int table[32] = {
    0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8,
    31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9};
  return table[((n & -n) * 0x077CB531UL) >> 27];
}

#endif


// Cast a double pointer to the encapsulating DVector.
#define TO_DVECTOR(d) ((DVector*) ((char*)(d) - GPO_OFFSETOF(DVector, data)))


// Exportable iterator.
// @@@TODO: Does not handle 64 bit indexes.
// @@@TODO: The types of i and __m should be int,IndexSet.
// @@@TODO: Don't propagate this macro to files that don't need it.
#define FOR_ALL_INDEXES(i, mask) \
  for (IndexSet __m = (mask); __m; __m = __m & (__m-1)) { \
    int i = lowest_bit_index_32(__m);


// Number class that has the same semantics as 'double', except that it also
// remembers partial derivatives with respect to dependent variables.

class Number {
 public:
  //---------- Initialization

  Number() {
    // Doubles declared without initialization are undefined, so there's no
    // reason this class should be different (unnecessary initialization would
    // waste time).
    //   x_ = 0;
    //   mask_ = 0;
#ifdef EXTERNAL_VECTOR
    diff_ = alloc_vector();             // @@@TODO: If first operation is assignment, this alloc is unnecessary,
#endif
  }

  // Conversion from double to Number. Implicit conversion can hide errors,
  // so make this explicit.
  explicit Number(double x) {
    x_ = x;
    mask_ = 0;
#ifdef EXTERNAL_VECTOR
    diff_ = alloc_vector();
#endif
  }

  // Copy constructor
  Number(const Number &b) {
    x_ = b.x_;
    mask_ = b.mask_;
#ifdef USE_REFCOUNTING
    diff_ = b.diff_;
    DVector *v = TO_DVECTOR(diff_);
    v->refs++;
#else
  #ifdef EXTERNAL_VECTOR
    diff_ = alloc_vector();
  #endif

  #ifdef USE_MEMCPY
    memcpy(diff_, b.diff_, IndexSetSize * sizeof(double));
  #else
    _FOR_ALL_INDEXES(i, mask_, diff_[i] = b.diff_[i];)
  #endif
#endif
  }

  // Destructor
  ~Number() {
#ifdef EXTERNAL_VECTOR
    release_vector();
#endif
  }

  // Assign a scalar, all derivatives are set to 0.
  void operator = (double x) {
    x_ = x;
    mask_ = 0;
  }

  // Assign to another Number.
  void operator = (const Number &b) {
    // Make sure self-assignment is handled properly
    x_ = b.x_;
    mask_ = b.mask_;
#ifdef USE_REFCOUNTING
    DVector *v = TO_DVECTOR(b.diff_);
    v->refs++;
    release_vector();
    diff_ = b.diff_;
#else
  #ifdef USE_MEMCPY
    memcpy(diff_, b.diff_, IndexSetSize * sizeof(double));
  #else
    _FOR_ALL_INDEXES(i, mask_, diff_[i] = b.diff_[i];)
  #endif
#endif
  }

  // We *could* have a cast-to-double operator, but in practice this hides
  // errors. Instead, use the to_double() functions below.
  //
  //    operator double() const { return x_; }

  // Make this Number into a dependent variable of index 'i', and assign it
  // a value.
  void SetDependent(int i, double value) {
    CHECK(i < IndexSetSize);
#ifdef USE_REFCOUNTING
    DVector *v = TO_DVECTOR(diff_);
    if (v->refs > 1) {
      //@@@TODO: This part not tested yet
      release_vector();
      diff_ = alloc_vector();
    }
#endif
    x_ = value;
    mask_ = 1 << i;
    diff_[i] = 1;
  }

  //---------- Accessors

  // Return the derivative d/dk[i].
  double Derivative(int i) const {
    return index_in_set(mask_, i) ? diff_[i] : 0.0;
  }

  // For the first n derivatives d/dk[0]..d/dk[n-1], write the nonzero ones
  // into the deriv[0],deriv[skip],deriv[2*skip],... array.
  void GetNonzeroDerivatives(int n, double *deriv, int skip = 1) const {
    _FOR_ALL_INDEXES(i, mask_, deriv[i*skip] = diff_[i];)
  }

  // Get the index set of nonzero derivatives.
  IndexSet GetNonzeroDerivativeSet() {
    return mask_;
  }

  //---------- Plus operator

  Number operator+ (const Number &y) const {
    Number ret;
    ret.x_ = (x_ + y.x_);
    ret.mask_ = mask_ | y.mask_;
    _FOR_ALL_INDEXES(i, mask_ & (~y.mask_), ret.diff_[i] = diff_[i];)
    _FOR_ALL_INDEXES(i, mask_ & y.mask_, ret.diff_[i] = diff_[i] + y.diff_[i];)
    _FOR_ALL_INDEXES(i, y.mask_ & (~mask_), ret.diff_[i] = y.diff_[i];)
    return ret;
  }

  Number operator+ (double y) const {
    Number ret;
    ret.x_ = x_ + y;
    ret.mask_ = mask_;
    _FOR_ALL_INDEXES(i, mask_, ret.diff_[i] = diff_[i];)
    return ret;
  }

  void operator += (const Number &y) {
    x_ += y.x_;
    _FOR_ALL_INDEXES(i, mask_ & y.mask_, diff_[i] += y.diff_[i];)
    _FOR_ALL_INDEXES(i, y.mask_ & (~mask_), diff_[i] = y.diff_[i];)
    mask_ |= y.mask_;
  }

  void operator += (double y) {
    x_ += y;
  }

  //---------- Minus operator

  Number operator- (const Number &y) const {
    Number ret;
    ret.x_ = x_ - y.x_;
    ret.mask_ = mask_ | y.mask_;
    _FOR_ALL_INDEXES(i, mask_ & (~y.mask_), ret.diff_[i] = diff_[i];)
    _FOR_ALL_INDEXES(i, mask_ & y.mask_, ret.diff_[i] = diff_[i] - y.diff_[i];)
    _FOR_ALL_INDEXES(i, y.mask_ & (~mask_), ret.diff_[i] = -y.diff_[i];)
    return ret;
  }

  Number operator- (double y) const {
    Number ret;
    ret.x_ = x_ - y;
    ret.mask_ = mask_;
    _FOR_ALL_INDEXES(i, mask_, ret.diff_[i] = diff_[i];)
    return ret;
  }

  void operator -= (const Number &y) {
    x_ -= y.x_;
    _FOR_ALL_INDEXES(i, mask_ & y.mask_, diff_[i] -= y.diff_[i];)
    _FOR_ALL_INDEXES(i, y.mask_ & (~mask_), diff_[i] = -y.diff_[i];)
    mask_ |= y.mask_;
  }

  void operator -= (double y) {
    x_ -= y;
  }

  //---------- Unary minus operator

  Number operator- () const {
    Number ret;
    ret.x_ = -x_;
    ret.mask_ = mask_;
    _FOR_ALL_INDEXES(i, ret.mask_, ret.diff_[i] = -diff_[i];)
    return ret;
  }

  //---------- Multiply operator

  Number operator* (const Number &y) const {
    Number ret;
    ret.x_ = x_ * y.x_;
    ret.mask_ = mask_ | y.mask_;
    _FOR_ALL_INDEXES(i, mask_ & (~y.mask_), ret.diff_[i] = y.x_*diff_[i];)
    _FOR_ALL_INDEXES(i, mask_ & y.mask_, ret.diff_[i] = y.x_*diff_[i] + x_*y.diff_[i];)
    _FOR_ALL_INDEXES(i, y.mask_ & (~mask_), ret.diff_[i] = x_*y.diff_[i];)
    return ret;
  }

  Number operator* (double y) const {
    Number ret;
    ret.x_ = x_ * y;
    ret.mask_ = mask_;
    _FOR_ALL_INDEXES(i, mask_, ret.diff_[i] = diff_[i] * y;)
    return ret;
  }

  void operator *= (const Number &y) {
    //@@@TODO: Has not been tested yet.
    x_ *= y.x_;
    _FOR_ALL_INDEXES(i, mask_ & y.mask_, diff_[i] = y.x_*diff_[i] + x_*y.diff_[i];)
    _FOR_ALL_INDEXES(i, y.mask_ & (~mask_), diff_[i] *= y.x_;)
    mask_ |= y.mask_;
  }

  void operator *= (double y) {
    x_ *= y;
    _FOR_ALL_INDEXES(i, mask_, diff_[i] *= y;)
  }

  //---------- Divide operator

  Number operator/ (double y) const {
    Number ret;
    ret.x_ = x_ / y;
    ret.mask_ = mask_;
    _FOR_ALL_INDEXES(i, mask_, ret.diff_[i] = diff_[i] / y;)
    return ret;
  }

  void operator /= (double y) {
    x_ /= y;
    _FOR_ALL_INDEXES(i, mask_, diff_[i] /= y;)
  }

  //---------- Debugging

  void PrintDerivatives(const char *name) const {
    printf("%s=%g ", name, x_);
    for (int i = 0; i < IndexSetSize; i++) {
      if (Derivative(i) != 0) printf("d/dk%d=%g ", i, Derivative(i));
    }
    printf("\n");
  }

  //---------- Data
private:
  double x_;                    // The value of the number
  IndexSet mask_;               // Bitmask indicating nonzero differentials
#ifdef EXTERNAL_VECTOR
  double *diff_;                // Vector of partial differentials
#else
  double diff_[IndexSetSize];   // Vector of partial differentials
#endif

  // Friends
  friend double to_double(const Number &y);
  friend Number sin(const Number &y);
  friend Number cos(const Number &y);

  //---------- Support for externally stored vectors

  struct DVector {              // Vector of partial differentials
    DVector *next;              // Next item in allocated or free list
    int refs;                   // References to this vector, if allocated
    AlignedDouble data[IndexSetSize];
  };

  static DVector *first_free;   // Head pointer for list of free vectors

  static double *alloc_vector() __attribute__((malloc)) {
    if (!first_free) AllocateMoreDVectors();
    DVector *v = first_free;
#ifdef USE_REFCOUNTING
    v->refs = 1;
#endif
    first_free = first_free->next;
    return v->data;
  }

  void release_vector() {
    DVector *v = TO_DVECTOR(diff_);
#ifdef USE_REFCOUNTING
    v->refs--;
    if (v->refs == 0) {
#endif
      v->next = first_free;
      first_free = v;
#ifdef USE_REFCOUNTING
    }
#endif
  }

  static void AllocateMoreDVectors();
};


// Convert a Number into a double.
inline double to_double(const Number &y) {
  return y.x_;
}

inline Number sin(const Number &y) {
  Number ret;
  ret.x_ = ::sin(y.x_);
  double k = ::cos(y.x_);
  ret.mask_ = y.mask_;
  _FOR_ALL_INDEXES(i, ret.mask_, ret.diff_[i] = k * y.diff_[i];)
  return ret;
}

inline double sin(double x) {
  return ::sin(x);
}

inline Number cos(const Number &y) {
  Number ret;
  ret.x_ = ::cos(y.x_);
  double k = -::sin(y.x_);
  ret.mask_ = y.mask_;
  _FOR_ALL_INDEXES(i, ret.mask_, ret.diff_[i] = k * y.diff_[i];)
  return ret;
}

inline double cos(double x) {
  return ::cos(x);
}


// Don't pollute the macro namespace.

#undef TO_DVECTOR
#undef INDEX_SET_64
#undef EXTERNAL_VECTOR
#undef USE_REFCOUNTING
#undef USE_MEMCPY
#undef _FOR_ALL_INDEXES

}; // namespace GPO

#endif
