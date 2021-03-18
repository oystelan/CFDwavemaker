# CFDwavemaker v2.1.6
wave kinematics package for CFD initialization of CFD domains or boundaries.
The code supports OpenMP and has some special features for fast initialization of large domains built in.

User manual for CFDwavemaker available at
http://www.hydrodynamics.no/CFDwavemaker

Supported wave theories: 
- linear wave theory (irregular & short crested waves)
- second order wave theory (Sharma & Dean, with taylor expansion above z=0) support for irregular & short crested waves 
- Stokes 5th regular waves (long crested)
- wave maker theory (piston type wave maker), reading paddle motion from input file (long crested)
- Spectral-wave-data library, which allows for extension of 

outputs: 
- surface elevation at any position (x,y)
- velocity components (ux, uy, uz)

Special:
- Grid interpolation of wave kinematics for faster initialization when using second order kinematics


## Build instructions

1. navigate into the folder src/
2. run make and the serial version of CFDwavemaker will be built. Alternatively provide an alternative build option:
    - make openmp
    - make openmp_swd
    - make openmp_vtk

## Version log (a coarse one)

Version 2.1.6 release:

Version 2.1.5 release:
- added probes as output such that x,y,z coordinates may be specified to dump kinematics to file.
- Fixed CFDwavemaker for windows. now compiles without SWD.
- some cleanup of the code.

Version 2.1.4 release:
- Improved performance with openmp and irregular waves
- minor bug fixes
- added parameter [vtk output]/timelabel and [lsgrid]/init_only.

Version 2.1.3 highlights (since version 1)
- complete restructuring of input file format
- implementation of a new Lagrangian stretched grid interpolation scheme
- implemented support for spectral wave data files (.swd) through the external library Spectral-Wave-Data (https://github.com/SpectralWaveData/spectral_wave_data). This facilitates the use of higher order spectral methods for CFD initialization and at boundaries(HOSM)
- General cleanup of code and removal of old leftovers...
- The entire code have been restructured to make use of classes
- previous stokes c code is now converted to a c++ class.
- fixed a large number of holes and bugs in the original code