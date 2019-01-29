# visag - VIscous Semi-Analytic self-Gravitating discs

This repository computes the viscous evolution of protostellar discs in 1D.  It either assumes a fixed Shakura-Sunyaev alpha parameter, or computes alpha as produced by self-gravitating disc turbulence in local thermodynamic equilibrium.

Features
----------------

* Viscous evolution using Shakura-Sunyaev alpha prescription
* Models run with fixed alpha or alpha due to gravitational instability (assuming local thermodynamic equilibrium)
* Mass loss due to X-Ray photoevaporative winds (Owen et al 2011, Alexander and Armitage 2007,2009)
* Tidal torques due to planetary bodies, modelling both Type I and Type II migration (Nayakshin 2015)

Features in development
------------------------

* Gravitational interactions between embedded planets (see `nbody` branch)

Future Features
---------------

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

The input parameters are specified in `visag.params` - an example file is given in `params/`.  The code also reads `myeos.dat` from the run directory (this can also be found in `params/`).

If planets are to be added, the user must specify a file to be read in `visag.params` that contains the planet data.  The planet file must have the following format:

number of planets
mass semimajor_axis
mass semimajor_axis

Masses must be in Jupiter masses, semimajor axis in AU.

Plotting
--------

The output files can be plotted using Python scripts found in the `plot/` directory

The scripts were developed in Python 2.7, and depend on numpy and matplotlib

License
-------

This code is licensed under the MIT license - see LICENSE.md
