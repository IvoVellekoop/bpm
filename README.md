# bpm
supplementary scripts for arXiv:2002.05279. Note that this repository may not be actively maintained!

This repository contains a rudimentary beam propagation method, using the angular spectral method.
The scripts were used for the model-based wavefront shaping work, reported in arXiv:2002.05279

Field.m = Class representing a 2-D electric field. The field can be propagated through structure with a 3-D refractive index distribution n using the propagate() method.

Unit.m = Helper class for working with units

SizedArray.m = Helper class to have a MATLAB array with an attached unit. Allows plotting (automatically putting the correct units on the axes) and ffts (automatically converting units from real space to k-space correctly)
