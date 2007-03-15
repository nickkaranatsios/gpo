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

/*

INTERNAL macros to manipulate small matrices and vectors. In particular we care
about the following operations (the notation here is that 'a' is a 3x1 vector
and 'A' is a 3x3 matrix).

        k = a dot b                -- Dot product
        k = a cross b              -- Cross product
        a = B * c
        a += B * c
        a -= B * c
        a = scale * B * c
        a += scale * B * c
        a -= scale * B * c
        a = B' * c
        a += B' * c
        a -= B' * c
        a = scale * B' * c
        a += scale * B' * c
        a -= scale * B' * c
        A = B * C
        A += B * C
        A -= B * C
        A = scale * B * C
        A += scale * B * C
        A -= scale * B * C
        A = B' * C
        A += B' * C
        A -= B' * C
        A = scale * B' * C
        A += scale * B' * C
        A -= scale * B' * C
        A = B * C'
        A += B * C'
        A -= B * C'
        A = scale * B * C'
        A += scale * B * C'
        A -= scale * B * C'

A small number of macros is used to represent this variety of functions by
supplying the assignment operation as a macro argument 'op' that is inserted
into the correct place in the text expansion, like this:

     A = B * C               op <-- =
     A += B * C              op <-- +=
     A -= B * C              op <-- -=
     A = scale * B * C       op <-- =scale*
     A += scale * B * C      op <-- +=scale*
     A -= scale * B * C      op <-- -=scale*

The convention for naming matrix multiplication macros is MULTIPLYx_pqr
where x indicates the matrix that is transposed (if any) and p,q,r are the
sizes of the three matrices:
     x=0:   A = B  * C   (sizes: A:p*r B:p*q C:q*r)
     x=1:   A = B' * C   (sizes: A:p*r B:q*p C:q*r)
     x=2:   A = B  * C'  (sizes: A:p*r B:p*q C:r*q)
Cases x=1 and x=2 are equivalent to saying that the operation is A=B*C but
B or C are stored in column-wise format.

All matrices are stored row-wise, e.g. a 3x3 matrix is stored as follows:
     [0] [1] [2]
     [3] [4] [5]
     [6] [7] [8]

Operations where A, B and C have a mixture of element types are handled
transparently (since the macros are not type sensitive).

----------

C++ purists will complain that these macros should actually be replaced by
classes and inline functions. However, that results in either more complex
source, or code that takes longer to run. Let's examine why.

The first option is to simply implement each of the desired operations with a
separate inline function. This results in a combinatorial explosion of such
functions because of the variety of assignment operators, and therefore the
implementation is more error prone and difficult to maintain.

The second option is to use operator overloading on a matrix class such as
Matrix3x3_d. Unfortunately C++ compilers will generate matrix temporaries for
expressions such as 'A+=B*C', pessimizing the generated code.

The third option is to collapse many of the assigment operations together by
implementing functions of the form:

        A = k1*A + k2*B*C

where k1 and k2 and set to {-1,0,1} to match the desired operation. This
generates rather sub-optimal code on gcc-4.1 as that compiler is unable to do
constant propagation on floating point non-scalar (array, pointer) variables in
expressions. This is probably due to IEE754 rules about propagation of NaNs.
For example, the following code:

        void test3 (double a[1], double b[1]) {
                a[0] = 0.0*a[0] + 0.0*b[0];
        }

is compiled by gcc-4.1 into (x86 assembly):

      _Z5test3PdS_:
        xorpd %xmm0, %xmm0
        movlpd (%rdi), %xmm1
        mulsd %xmm0, %xmm1            <-- multiply by zero operation!
        mulsd (%rsi), %xmm0           <-- multiply by zero operation!
        addsd %xmm0, %xmm1
        movsd %xmm1, (%rdi)
        ret

You can see that the two multiplies-by-zero are explicitly represented in the
assembly code. The consequence is that in the common case of k1=0 and k2=1
there would be two redundant loads and multiplies per element, which pessimizes
the code somewhat. This may be improved by using the -ffast-math option, but
that can change the behavior of some code in unexpected ways.

The macro approach can be easily used on vectors that are embedded in the middle
of larger vectors, e.g.:

        double *foo = new double [1000];
        MULTIPLY1_331 (target,+=,matrix,foo+500);       // Use foo[500..502]

Another issue with the class/function based approach is that it requires
templatization to handle cases where the matrix/vector variables have different
floating point or integer types. This adds an extra layer of complexity and
introduces further impediments to compiler optimization.

----------

The following is done to prevent accidental misuse of these macros:
  - They're all prefixed with FAST_, to help prevent namespace collisions.
  - All these macros take arguments, so the compiler will catch cases where
    the macro names are used as variables (a syntax error will result).
  - Cases where the macro names are used as functions will be caught if the
    call site uses a different number of arguments. If the same number of
    arguments is used then the definition site will result in a syntax error,
    e.g. void FAST_DOT(int i, int j) expands to:
        void ((int i)[0]*(int j)[0] + (int i)[1]*(int j)[1] + ...)

*/

#ifndef __GPO_PMATH_H__
#define __GPO_PMATH_H__

// Dot product.
#define FAST_DOT(a, b) ((a)[0]*(b)[0] + (a)[1]*(b)[1] + (a)[2]*(b)[2])

// Cross product
#define FAST_CROSS(a, op, b, c) { \
  (a)[0] op((b)[1]*(c)[2] - (b)[2]*(c)[1]); \
  (a)[1] op((b)[2]*(c)[0] - (b)[0]*(c)[2]); \
  (a)[2] op((b)[0]*(c)[1] - (b)[1]*(c)[0]); \
}

// Set the elements in a 3x3 matrix
#define FAST_SET_MATRIX3(M, a0, a1, a2, a3, a4, a5, a6, a7, a8) { \
  (M)[0] = (a0); (M)[1] = (a1); (M)[2] = (a2); \
  (M)[3] = (a3); (M)[4] = (a4); (M)[5] = (a5); \
  (M)[6] = (a6); (M)[7] = (a7); (M)[8] = (a8); \
}

// Set the elements in a 3x3 matrix to hat(v)
#define FAST_SET_CROSSMAT(M, v) { \
  (M)[0] =     0; (M)[1] = -v[2]; (M)[2] =  v[1]; \
  (M)[3] =  v[2]; (M)[4] =     0; (M)[5] = -v[0]; \
  (M)[6] = -v[1]; (M)[7] =  v[0]; (M)[8] =     0; \
}

// a = B*c
#define FAST_MULTIPLY0_331(a, op, B, c) FAST_MULTIPLY0_331s(a, op, B, c, 1)
#define FAST_MULTIPLY0_331s(a, op, B, c, askip) { \
  (a)[0]         op((B)[0]*(c)[0] + (B)[1]*(c)[1] + (B)[2]*(c)[2]); \
  (a)[askip]     op((B)[3]*(c)[0] + (B)[4]*(c)[1] + (B)[5]*(c)[2]); \
  (a)[2*(askip)] op((B)[6]*(c)[0] + (B)[7]*(c)[1] + (B)[8]*(c)[2]); \
}

// a = B'*c
#define FAST_MULTIPLY1_331(a, op, B, c) FAST_MULTIPLY1_331s(a, op, B, c, 1)
#define FAST_MULTIPLY1_331s(a, op, B, c, askip) { \
  (a)[0]         op((B)[0]*(c)[0] + (B)[3]*(c)[1] + (B)[6]*(c)[2]); \
  (a)[askip]     op((B)[1]*(c)[0] + (B)[4]*(c)[1] + (B)[7]*(c)[2]); \
  (a)[2*(askip)] op((B)[2]*(c)[0] + (B)[5]*(c)[1] + (B)[8]*(c)[2]); \
}

// A = B*C
#define FAST_MULTIPLY0_333(A, op, B, C) { \
  (A)[0] op((B)[0]*(C)[0] + (B)[1]*(C)[3] + (B)[2]*(C)[6]); \
  (A)[1] op((B)[0]*(C)[1] + (B)[1]*(C)[4] + (B)[2]*(C)[7]); \
  (A)[2] op((B)[0]*(C)[2] + (B)[1]*(C)[5] + (B)[2]*(C)[8]); \
  (A)[3] op((B)[3]*(C)[0] + (B)[4]*(C)[3] + (B)[5]*(C)[6]); \
  (A)[4] op((B)[3]*(C)[1] + (B)[4]*(C)[4] + (B)[5]*(C)[7]); \
  (A)[5] op((B)[3]*(C)[2] + (B)[4]*(C)[5] + (B)[5]*(C)[8]); \
  (A)[6] op((B)[6]*(C)[0] + (B)[7]*(C)[3] + (B)[8]*(C)[6]); \
  (A)[7] op((B)[6]*(C)[1] + (B)[7]*(C)[4] + (B)[8]*(C)[7]); \
  (A)[8] op((B)[6]*(C)[2] + (B)[7]*(C)[5] + (B)[8]*(C)[8]); \
}

// A = B'*C
#define FAST_MULTIPLY1_333(A, op, B, C) { \
  (A)[0] op((B)[0]*(C)[0] + (B)[3]*(C)[3] + (B)[6]*(C)[6]); \
  (A)[1] op((B)[0]*(C)[1] + (B)[3]*(C)[4] + (B)[6]*(C)[7]); \
  (A)[2] op((B)[0]*(C)[2] + (B)[3]*(C)[5] + (B)[6]*(C)[8]); \
  (A)[3] op((B)[1]*(C)[0] + (B)[4]*(C)[3] + (B)[7]*(C)[6]); \
  (A)[4] op((B)[1]*(C)[1] + (B)[4]*(C)[4] + (B)[7]*(C)[7]); \
  (A)[5] op((B)[1]*(C)[2] + (B)[4]*(C)[5] + (B)[7]*(C)[8]); \
  (A)[6] op((B)[2]*(C)[0] + (B)[5]*(C)[3] + (B)[8]*(C)[6]); \
  (A)[7] op((B)[2]*(C)[1] + (B)[5]*(C)[4] + (B)[8]*(C)[7]); \
  (A)[8] op((B)[2]*(C)[2] + (B)[5]*(C)[5] + (B)[8]*(C)[8]); \
}

#endif
