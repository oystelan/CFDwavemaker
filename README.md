# CFDwavemaker v2
wave kinematics package for CFD initialization of CFD domains or boundaries.
The code supports OpenMP and has some special features for fast initialization of large domains built in.

Supported wave theories: 
- linear wave theory (irregular & short crested waves)
- second order wave theory (Sharma & Dean, with taylor expansion above z=0) support for irregular & short crested waves 
- Stokes 5th regular waves (long crested)
- wave maker theory (piston type wave maker), reading paddle motion from input file (long crested)

outputs: 
- surface elevation at any position (x,y)
- velocity components

Version2 will replace the old version as soon as tests have been carried out, but is not ready for use, yet....
New to Version 2:
- The entire code have been restructured to make use of classes
- previous stokes c code is now converted to a c++ class.
To come:
- a new wave spectrum class for generating waves directly from various wave spectra will soon be included.
- a spreading function class
- deans stream function wave
- hinghed wave maker theory

