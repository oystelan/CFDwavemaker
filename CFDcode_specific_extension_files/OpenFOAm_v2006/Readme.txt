CFDwavemaker is now been implemented as an extension to the official waveModels implementation in OpenFOAM-v2006.
CFD wavemaker may now be applied as

1. on boundaries
2. as euler-overlay (adapted implementation of Tormod Landet, which was originally made for openfoam.org version of openfoam).


Installation instructions
=========================

1. install openfoam-v2006 from source. Download and unpack to a suitable directory.
2. copy all files and folders in the extension_src folder to the OpenFOAM-v2006/src directory.
3. copy a precompiled fresh version of CFDwavemaker (libCFDwavemaker_openmp.a) into the folder: OpenFOAM-v2006/src/waveModels/waveGenerationModels/derived/CFDwavemaker/ (replace the old version already located there)
4. locate OpenFOAM-v2006/src/waveModels/Make/files and add the following two lines to the list of files.

waveGenerationModels/derived/CFDwavemaker/CFDwavemakerWaveModel.C
waveDampToIncident/WaveDampToIncident.C

4. Now, compile openfoam. Step 1-3 may be done after openfoam has been compiled, just remember to recompile the wave library. This is most easily done by running OpenFOAM-v2006/src/Allrun.
If you compile openfoam for the first time, be patient. It takes a while! (typically half a day).


