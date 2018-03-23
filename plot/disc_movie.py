#!/usr/local/bin/python

# Python script makes a movie of a disc variable

import matplotlib.pyplot as plt
import numpy as np
import os

outputstring=['r','sigma','cs','kappa','gamma','mu', 'T', 'tau', 'alpha', 'Q']
logplot = [False,True,True,True,False,False,True,True,True,False]

ylabels = [r'r (AU)',r'$\Sigma$ (g cm $^{-2}$)',r' $c_s$ (cm s$^{-1}$)',r'$\kappa$ (cm$^{2}$ g$^{-1}$)',
           r'$\gamma$',r'$\mu$',r'$T_c$ (K)',r'$ \tau $',r'$\alpha_{g}$',r'$ Q $']
xlabel = ylabels[0]

ymin = []
ymax = []

for i in range(len(outputstring)):
    ymin.append(0.0)
    ymax.append(0.0)

# Set canonical minimum, maximum y limits

ymin[1]=1.0e0
ymax[1]= 1.0e6

ymin[6] = 1.0e0
ymax[6] = 2.0e4

ymin[8] = 1.0e-8
ymax[8] = 1.0e0

ymin[9]=0.0
ymax[9] = 10.0

# Input parameters
prefix = raw_input('What is the file prefix? ')
print 'Now select variable to plot: here are the choices'

for i in range(len(outputstring)):
    print str(i+1)+': '+outputstring[i]

var = input('Which variable to plot?')

var = var-1

initial = input('Starting filenumber? ')
nfiles = input('Final filenumber? ')

nprofcol=10
nlayercol=11
nlogcol=10

print 'You have selected ',outputstring[var]

# Read logfile to get times

logfile = prefix+'.log'

# Open log file, to find time of dump

logdata = np.genfromtxt(logfile)

time = logdata[:,0]

# Go interactive

plt.ion()

fig1 = plt.figure()
ax = fig1.add_subplot(111)
ax.set_xlabel('r (AU)')
ax.set_ylabel(outputstring[var])


# Begin loop
for i in range(initial,nfiles):

    proffile=prefix+'profile.'+str(i+1)
    if(logplot[var]):
        ax.set_xscale('log')
        ax.set_yscale('log')
    if(ymin[var]!=ymax[var]):
        ax.set_ylim(ymin[var],ymax[var])
    # Read data
    print 'Reading file ',proffile
    data = np.genfromtxt(proffile,skip_header=1)
    nrow = data.size/nprofcol
    data.reshape(nrow,nprofcol)

    sigma = data[:,1]

    line1 = ax.plot(data[:,0],data[:,var])    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabels[var])
    #ax.text(0.9, 0.9,'t = '+str(time[i])+' yr', bbox=dict(edgecolor='black',facecolor='none'), horizontalalignment='center', verticalalignment='center',transform = ax.transAxes)
    plt.savefig(outputstring[var]+str(i+1)+'.png', format='png')
    plt.draw()

    ax.clear()

# Now make an animated gif

giffile = outputstring[var]+'.gif'

print 'Converting output files to animated gif'
print 'If this fails, check local commands for ImageMagick'
os.system('/usr/bin/convert ' +outputstring[var]+ '*.png '+giffile )

print 'File ', giffile, ' written'




