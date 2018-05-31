#!/usr/local/bin/python

# Python script makes a movie of a disc variable

import io_disc as io
import filefinder as ff

prefix = ff.get_file_prefix('*.log')

# Input parameters
planets = raw_input('Add planets? (y/n) ')

add_planets = False
if('y' in planets or 'Y' in planets):
    add_planets = True

io.plot_profile_multifiles_variable(prefix,add_planets=add_planets)
