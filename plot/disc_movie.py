#!/usr/local/bin/python

# Python script makes a movie of a disc variable

import io_disc as io


# Input parameters
prefix = raw_input('What is the file prefix? ')

io.plot_profile_multifiles_variable(prefix)
