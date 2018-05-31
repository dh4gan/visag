#####################################################
###                                               ###
###             Makefile for visag   	  	  ###
###                                               ###
###         Duncan H. Forgan 27/04/2018           ###
###       				          ###
###                                               ###
#####################################################

# Compiler variables:
FC     = gfortran
VPATH = src/main/ src/eos/ src/io/ src/layer/ src/planet/ src/setup/ src/wind/

# For real*8 files
FFLAGS = -O3 -fdefault-real-8 -fbounds-check

# Create object files:
%.o: %.f
	$(FC) $(FFLAGS) -c $<
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

SOURCESAF90 = grav_module.f90 mri_module.f90 planet_module.f90 \
	      unit_module.f90 wind_module.f90 main.f90 \
	      calc_typeI_migration.f90 compute_planet_torques.f90 \
	      compute_wind.f90 disc_properties.f90 \
	      eos.f90 eosread.f90 eos_T.f90 eos_cs.f90 \
              evolve.f90 evolve_layer.f90 \
	      layer_properties.f90 midplane_properties.f90 \
	      midplane_properties_grav.f90 midplane_properties_fixedalpha.f90 \
	      midplane_properties_grav_fixedQ.f90  migrate_planets.f90 \
              luminosity.f90 setup.f90 \
	      setup_accrete.f90 setup_planets.f90 setup_wind.f90 timestep.f90 \
	      wind_profiles.f90 xraywind.f90 write_dump.f90 	      

OBJECTSA    = $(SOURCESAF90:.f90=.o)

# Create executable files:
build: visag

visag:  $(OBJECTSA)
	$(FC) $(FFLAGS) -o $@ $(OBJECTSA)

# Clean statements:
clean: 
	\rm *.o *.mod visag

# End Makefile
