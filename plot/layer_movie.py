#!/usr/local/bin/python

# Python script makes a movie of the two disc layers

import matplotlib.pyplot as plt
import numpy as np

# Input parameters
prefix = raw_input('What is the file prefix? ')
initial = input('Starting filenumber? ')
nfiles = input('Final filenumber? ')

nprofcol=10
nlayercol=11
nlogcol=10

# Read logfile to get times

logfile = prefix+'.log'

# Open log file, to find time of dump

logdata = np.genfromtxt(logfile)

time = logdata[:,0]

# Go interactive

plt.ion()

fig1 = plt.figure()
ax = fig1.add_subplot(111)


# Begin loop
for i in range(initial,nfiles):

    ax.set_xlabel('r (AU)')
    ax.set_ylabel(r'$\Sigma$ (g cm$^{-2}$)')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(1e1,1e6)
    layerfile=prefix+'layer.'+str(i+1)

    # Read data
    print 'Reading file ',layerfile
    layer = np.genfromtxt(layerfile,skip_header=1)
    nrow = layer.size/nlayercol
    layer.reshape(nrow,nlayercol)

    sigma = layer[:,1]

    line1 = ax.plot(layer[:,0],layer[:,3], 'g')
    line2 = ax.plot(layer[:,0],layer[:,1], 'b')
    line3=ax.plot(layer[:,0],layer[:,2], 'r')
    
    ax.text(0.9, 0.9,'t = '+str(time[i])+' yr', bbox=dict(edgecolor='black',facecolor='none'), horizontalalignment='center', verticalalignment='center',transform = ax.transAxes)

    
    plt.draw()

    outputfile=prefix+'sigmas_'+str(i+1)+'.png'

    ax.clear()
#  



