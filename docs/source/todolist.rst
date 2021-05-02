Todo-list
=========

Here is a list of planned upgrades,ideas and improvements to the library:

- Update manual (always ongoing)
- Implement sloping bottom (variable water depth)
- upgrade second order velocity profile model to a consistent second order wheeler type stretching , Ref. :cite:`smit2017nonlinear` 
- At some point it is probably a good idea to add a Cmake file for this project.
- some variables in the input file should be calculated automatically based on given input values, i.e. default values should be auto.
- implement wind turbulence models for two-phase problems
- add more illustrations in the manual (figures)
- implement time-trace dump for specified list of points
- Implement dopple shift to take into account a constant current. Redefine x og y for a given current (U, V)
x' = x – Ut
y' = y - Vt
t' = t

Task completed
--------------

- Implement Higher order spectral method connection through the swd library `spectral_wave_data`_ 
- Implement grid interpolation scheme for boundaries
- update grid interpolation structure to vertical lagrangian type, which follows the wave surface, to make more efficient and accurate.

Ideas rejected (or posponed)
----------------------------

- Add commonly used wave spectra and directional spreading functions for easy generation without specifying manually frequency amplitudes, phase etc. (work started)
- preprogrammed constrained new wave. kind of a nice to have.

.. _`spectral_wave_data`: https://github.com/SpectralWaveData/spectral_wave_data

