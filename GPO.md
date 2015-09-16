# Google Pose Optimizer #

The Google pose optimizer (GPO) is a C++ library that allows reconstruction of the pose of a sensor platform (i.e. its position and orientation over time) based on information from sensors such as GPS, accelerometers and rate gyroscopes.

The traditional way to solve this problem is to use a Kalman filter, which will give you the  pose in real time. GPO uses an off-line optimization approach pioneered by James Diebel at Stanford ([His software is here](http://ai.stanford.edu/~diebel/smoother.html)). The advantage of off-line optimization is that smoother and more accurate results are produced. The optimization problem is to take noisy sensor data over time and compute the pose at each time step that is most likely to have produced those measurements.

GPO has a sensor plugin architecture that allows new sensor models to be easily incorporated in to the optimization. The sensor models currently supported are:
  * GPS position and doppler velocity.
  * Accelerometers.
  * Rate gyroscopes.
  * Wheel encoders.
  * Simple road constraint models.

Note that GPO does not handle file I/O or any other application-specific messy details, it is purely a computational engine (numbers go in, numbers come out).

## Documentation ##

Yes, some documentation would be nice! For now here are some RandomNotes, eventually a user's manual will be written.