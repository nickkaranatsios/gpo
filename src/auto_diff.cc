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

// Automatic differentiation, support for externally stored vectors.

#include "gpo/auto_diff.h"

using namespace GPO;


Number::DVector *Number::first_free = 0;


void Number::AllocateMoreDVectors() {
  CHECK(first_free == 0);
  const int num_to_alloc = 256;
  //printf ("Alloc %d\n",num_to_alloc);
  DVector *arena = new DVector[num_to_alloc];
  for (int i = 0; i < num_to_alloc-1; i++) {
    arena[i].next = arena + i + 1;
  }
  arena[num_to_alloc-1].next = 0;
  first_free = arena;
}
