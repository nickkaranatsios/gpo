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
#include "gpo/gravimetry.h"


double GPO::ComputeGravity(double lat, double lng, double alt) {
  // Use the International Gravity Formula 1967 with a height correction.
  //@@@TODO: Improve this approximation. The IGF67 is for terrain at sea level,
  //         and the height correction is for free air only. Here 'alt' is
  //         taken to be height above sea level, but we really want it to mean
  //         height above the WGS84 ellipsoid. Thus this formula is only okay
  //         for coastal terrain.
  double phi = lat * (M_PI / 180.0);
  double sin_phi = sin(phi);
  double sin_squared_phi = sin_phi * sin_phi;
  double sin_2phi = sin(2.0 * phi);
  double sin_squared_2phi = sin_2phi * sin_2phi;
  return 9.780318 * (1.0 + 0.0053024*sin_squared_phi
         - 0.0000058*sin_squared_2phi) - 3.086e-6*alt;
}
