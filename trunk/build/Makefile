# Copyright 2007 Google Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Author: Russell L. Smith

#############################################################################
# Filenames

LIB=libgpo.a
LIB_SRC=optimizer.cc hessian.cc fast_block.cc gps.cc gps_velocity.cc \
	accelerometer.cc rate_gyro.cc axel.cc wheel.cc road.cc \
	test_matrix.cc rotation.cc common.cc gravimetry.cc sensor.cc \
	auto_diff.cc cache.cc
APP_SRC=test_gpo.cc
CONFIG_H=../include/config.h

#############################################################################
# Directories and derived files

VPATH=../src
APP=$(APP_SRC:%.cc=%.exe)
LIB_OBJ=$(LIB_SRC:%.cc=%.o)
APP_OBJ=$(APP_SRC:%.cc=%.o)

#############################################################################
# Compiler configuration

CPP=g++ -g -Wall -O2 -fomit-frame-pointer
# For profiling: -fno-inline -fno-default-inline -pg

all: $(LIB) $(APP)

$(LIB): $(CONFIG_H) $(LIB_OBJ)
	ar rcv $@ $(LIB_OBJ)

$(APP): $(LIB) $(APP_OBJ)
	$(CPP) -o $@ $(APP_OBJ) $(LIB)

%.o: %.cc
	$(CPP) -c -I.. $< -MD

$(CONFIG_H):
	gcc -o configure.exe configurator.c
	./configure.exe $(CONFIG_H) "g++ -o " ./

clean:
	-rm -f *.o *.d *.exe *.a $(CONFIG_H)

-include *.d
