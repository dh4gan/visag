# visag - VIscous Semi-Analytic self-Gravitating discs
======================================================

This repository computes the viscous evolution of protostellar discs in 1D.  It either assumes a fixed Shakura-Sunyaev alpha parameter, or computes alpha as produced by self-gravitating disc turbulence in local thermodynamic equilibrium.

Features
----------------

* Viscous evolution using Shakura-Sunyaev alpha prescription
* Models run with fixed alpha or alpha due to gravitational instabilities
* Mass loss due to photoevaporative winds (Owen et al 2011, Alexander and Armitage+++)
* Tidal torques due to planetary bodies (Type I and Type II migration)


Future Features
---------------

* Gravitational Interactions between embedded planets
* Layered accretion (GI/MRI mixtures)
* Accretion from external envelope/other sources
* (integration into `dh4gan/grapus`)

If you plan to use this code for your own research, or if you would like to contribute to this repository then please get in touch with a brief description of what you would like to do.  I will seek co-authorship for any subsequent publications.


Compiling and Running
---------------------

The code is written in FORTRAN 90 throughout. The supplied Makefile compiles the 
program using gfortran through the command

`> make`

And run using the command

`> ./visag`

The input parameters are specified in `visag.params` - an example file is given in `params/`

Plotting
--------

The output files can be plotted using Python scripts found in the `plot/` directory

The scripts were developed in Python 2.7, and depend on numpy and matplotlib

License
-------

This code is licensed under the MIT license - see LICENSE.md
