#      _GLaSS_                    #



# Background #

**statement**: This is the Generalised Lensing and Shear Spectra (_Glass_) code. It is a module for Cosmosis. See https://bitbucket.org/joezuntz/cosmosis/wiki/Home for more detials. 
Please contact Peter Taylor (**peterllewlyntaylor at gmail dot com**) 
if there are problems.
If you use _GLaSS_ please remember to cite the
papers listed below.


**name**: _GLaSS_

**version**: 1

**purpose**: Computes spherical Bessel, tomographic and generalised 
lensing spectra. 

**attribution**: Peter Taylor, Mullard Space Science Laboratory,
University College London, 2018

**cite**: 

- __Preparing for the Cosmic Shear Data Flood: Optimal Data Extraction and
Simulation Requirements for Stage IV Dark Energy Experiments__: Generalised Shear Spectra are introduced.
- __Testing the Cosmic Shear Spatially-Flat Universe Approximation with GLaSS__: _GLaSS_, is introduced.



**python dependencies**

- numpy (for fast computation build against a linear algebra package like BLAS or Apple Acccelerate) 
- os
- ctypesGsl
- math
- itertools


# Running _GLaSS_ #

 _GLaSS_ is called by Cosmosis. Run-mode options are specified in a Cosmosis ini file. _GLaSS_ can either be run like any other Cosmosis module or interactively in Python. See the jupyter notebooks to see how to run interactively.


# Program Structure #
_GLaSS_ is run internally in Cosmosis. Cosmosis first calls `setup(options)` in `GLaSS.py`. This reads in and returns all the run-mode options that
are specified in the Cosmosis ini file. Bessel function data is also computed inside this function as needed.`setup(options)`
is only run once, each time Cosmosis is called.
Then `exec(block, config)` is called. This function reads in the power spectrum data either from the Cosmosis pipeline or
an external source.	It then computes the lensing power spectra and returns this to the Cosmosis data block. The calculation
is performed by the `main` function in `lensing_calculation.py`. All the functions that are actually used internally to perform the
calculation of the lensing spectra are defined in `lensing_calculation.py`.

# Run-mode Options #
To see which run-mode options to use have a look at the ini files and jupyter lab notebook provided with _GLaSS_. It is also
worth looking at the .py files called in the jupyter notebook. A complete list of run-mode options
can be found by having a look at the code comments in the `setup(options)` function in `GLaSS.py`.

# Output #
_GLaSS returns the l samples, z samples, r(z) samples and k samples that were taken. Note that k is reutrned in units of Mpc ^ -1,
not h Mpc ^ -1, which is the conevention. The shot noise and the lensing spectra are returned in a 2D array with dimensions: (resolution, resolution * number_of_l_modes)
for spherical-bessel and arbitrary weighted lensing and (number_of_tomographic_bins, number_of_tomographic_bins * number_of_l_modes) for 
tomographic lensing spectra. The jupyter notebook shows how to exctract the relevant data.

 







   
   
   
   
   
