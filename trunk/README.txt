
Google Pose Optimizer (GPO)
---------------------------

Author:
	Russell L. Smith (russ.ls@gmail.com)

Home page and documentation:
	http://code.google.com/p/gpo/

Short description:
	The Google pose optimizer allows reconstruction of the pose of a
	sensor platform (i.e. its position and orientation) based on
	information from sensors such as GPS, accelerometers and rate
	gyroscopes. GPO does not provide real-time localization in the
	way that a Kalman filter would, instead it generates the pose as
	a result of a large off-line optimization. This produces better
	results.

License:
	This software is licensed under the Apache License, Version 2.0,
	see the file COPYING.txt for details.

Installation:
	GPO currently compiles on environments that have the GNU c++ compiler
	and GNU make. This includes Linux and Cygwin. Other environments such
	as Microsoft Visual C++ will be added if there is demand.

	To build the library and a test program, just type 'make'. This will
	generate build/libgpo.a and build/test_gpo.exe. To clean everything
	up type 'make clean'. The build uses some configuration detection
	logic that is hopefully reasonably general and robust. If it fails in
	your case, please file a bug.
