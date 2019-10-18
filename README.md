# CFDwavemaker
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

The master branch is the old version of CFDwavemaker (v1)
A new version is under development under branch "version2".
