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

// Gravimetric utility functions.

#ifndef __GPO_GRAVIMETRY_H__
#define __GPO_GRAVIMETRY_H__

namespace GPO {

// Given a position on the Earth's surface, return the magnitude of the
// acceleration due to gravity vector at that point. The lattitude and longitude
// are given in degrees, the altitude is given in meters above the WGS-84
// ellipsoid.
double ComputeGravity(double lat, double lng, double alt);

}; // namespace GPO

#endif
