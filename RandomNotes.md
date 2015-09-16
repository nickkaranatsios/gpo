# Random Notes #

These will eventually get moved to a real user's manual.

  * Coordinate Systems: The pose optimization occurs in Cartesian coordinates that are not explicity related in any way to the Earth. Meters and radians are the standard units. It is the job of the sensor modules to convert from Earth-related coordinates to Cartesian coordinates if necessary.

  * The State Vector: The state vector is made up of the pose at each time step, extra sensor-specific parameters for each time step, and a number of global parameters that are the same for all time steps. The pose position is a cartesian x,y,z, the pose rotation is three Euler angles. The state vector looks like:
```
     x,y,z,e1,e2,e3, x1,x2..,         <-- For time step 0
     ...
     x,y,z,e1,e2,e3, x1,x2..,         <-- For time step N-1
     p1,p2,p3...                      <-- Global parameters
```

