#!/usr/local/bin/python

# Python script plots profile data from semi_analytic_disc runs

#import matplotlib.pyplot as plt
import numpy as np
import io_disc as io
import filefinder as ff


prefix = raw_input("Please enter the file prefix: ")

logfile = ff.find_sorted_local_input_files(prefix+"*.log")

###################################
# Plot log data
####################################

logdata = io.plot_log_data(logfile)

print 'Plotting complete for file ',logfile