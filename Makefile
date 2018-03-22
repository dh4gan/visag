#####################################################
###                                               ###
###  Makefile for semi analytic disc routine   	  ###
###  (updated throughout to Fortran 90)           ###
###                                               ###
###         Duncan H. Forgan 22/04/2010           ###
###       				          ###
###                                               ###
#####################################################

# Compiler variables for WKMR/DHF:
FC     = gfortran


# For real*8 files
FFLAGS = -O3 -fdefault-real-8 -fbounds-check

# For 
#FFLAGS = -O3 -frecord-marker=4  -Wall

# Create object files:
%.o: %.f
	$(FC) $(FFLAGS) -c $<
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

SOURCESAF90 = grav_module.f90 mri_module.f90 planet_module.f90 \
	      unit_module.f90 wind_module.f90 main.f90 \
	      compute_planet_torques.f90 eos.f90 eosread.f90 eos_T.f90 eos_cs.f90 disc_properties.f90 \
              evolve.f90 evolve_layer.f90 layer_properties.f90 midplane_properties.f90 \
	      midplane_properties_grav.f90 midplane_properties_fixedalpha.f90 \
	      midplane_properties_grav_fixedQ.f90  migrate_planets.f90 \
              luminosity.f90 setup.f90 set_accrete.f90 write_dump.f90 	      

OBJECTSA    = $(SOURCESAF90:.f90=.o)

# Create executable files:
build: semi_analytic_disc

semi_analytic_disc:  $(OBJECTSA)
	$(FC) $(FFLAGS) -o $@ $(OBJECTSA)

# Clean statements:
clean: 
	\rm *.o *.mod semi_analytic_disc

# End Makefile
