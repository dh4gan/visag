#!/usr/local/bin/python

# Python script plots data from semi analytic disc runs

import matplotlib.pyplot as plt
import numpy as np

nlayercol=11
nlogcol = 10

prefix = raw_input('What is the file prefix? ')
fileno = input('Which dump is to be plotted? ')

# Input parameters

fileno = int(fileno)

layer = prefix+'layer.'+str(fileno)
logfile = prefix+'.log'

# Open log file, to find time of dump

logdata = np.genfromtxt(logfile)

time = logdata[fileno-1,0]

print 'Time is ',time ,' years'

# Open profile file and read

layerdata = np.genfromtxt(layer,skip_header=1)
layerdata.reshape(layerdata.size/nlayercol,nlayercol)

# Set up plot data for loop

# Filenames
outputstring=['r','sigma','sigma_m', 'sigma_tot', 'cs_m','kappa_m','gamma_m','mu_m', 'tau_m', 'nu_m', 'alpha_m']

# Plot labels
ylabels = [r'r (AU)',r'$\Sigma_g$ (g cm $^{-2}$)',r'$\Sigma_m$ (g cm $^{-2}$)',r'$\Sigma_{tot}$ (g cm $^{-2}$)',
           r' $c_{s,m}$ (cm s$^{-1}$)',r'$\kappa_m$ (cm$^{2}$ g$^{-1}$)',r'$\gamma_m$',r'$\mu_m$',r'$ \tau_m $',r'$\nu_{m}$',r'$ \alpha_m $']

# Log the y axis? True/False
setylog=[False, True, True,True,True,True,False,False,True,True,True]

# y limits - set defaults first
ymin=[]
ymax=[]
for i in range(nlayercol):
    ymin.append(0.0)
    ymax.append(0.0)


# Now define any non-default limits


# Now generate layer plots

for i in range(1,nlayercol):
    print 'Plotting ',outputstring[i]
    fig1 = plt.figure()
    ax = fig1.add_subplot(111)
    ax.set_xlabel(ylabels[0])
    ax.set_xscale('log')
    
    ax.set_ylabel(ylabels[i])
    
    if setylog:
        ax.set_yscale('log')
    
    if ymin[i]!=ymax[i]:
        print 'Setting limits ',ymin[i],ymax[i]
        ax.set_ylim(ymin[i],ymax[i])
    
    line = ax.plot(layerdata[:,0],layerdata[:,i], label='Time = '+str(time)+'yr')
    ax.legend()
    
    outputfile = prefix+outputstring[i]+'.ps'    
    plt.savefig(outputfile, format='ps')

# Finish layer plots

print 'Plotting complete for file ',layer

# Now plot log data

outputstring = ['t', 'mgrav','mmag','mdisc','grav_max','mag_max','sig_max','tot_lum','mdot_grav', 'mdot_mag']

# Plot labels
ylabels = [r't (yr)',r'$M_{grav}$ ($M_{\odot}$)',r'$M_{MRI}$ ($M_{\odot}$)',r'$M_{disc}$ ($M_{\odot}$)',
           r' $\Sigma_{grav,max}$ (g cm$^{-2}$)',r' $\Sigma_{MRI,max}$ (g cm$^{-2}$)',r' $\Sigma_{max}$ (g cm$^{-2}$)',
           r'$L_{tot}$ ($L_{\odot}$)', r'$\dot{M}_{grav}$ ($M_{\odot} yr^{-1}$)', r'$\dot{M}_{mag}$ $(M_{\odot} yr^{-1}$)']

# Log the y axis? True/False
setylog=[False, False, True,True,]

# y limits - set defaults first
ymin=[]
ymax=[]
for i in range(nlogcol):
    ymin.append(0.0)
    ymax.append(0.0)

for i in range(1,nlogcol):
    print 'Plotting ',outputstring[i]
    fig1 = plt.figure()
    ax = fig1.add_subplot(111)
    ax.set_xlabel(ylabels[0])
    ax.set_xscale('log')
    
    ax.set_ylabel(ylabels[i])
    
    if setylog:
        ax.set_yscale('log')
    
    if ymin[i]!=ymax[i]:
        print 'Setting limits ',ymin[i],ymax[i]
        ax.set_ylim(ymin[i],ymax[i])
    
    line = ax.plot(logdata[:,0],logdata[:,i])
    
    outputfile = prefix+outputstring[i]+'time.ps'    
    plt.savefig(outputfile, format='ps')

# Extra plots outside of the loop: plot mdot vs sigmax for grav and MRI layers

print 'Plotting mdot vs sigma for grav layer'

fig1 = plt.figure()
ax = fig1.add_subplot(111)
ax.set_xlabel(ylabels[4])
ax.set_ylabel(ylabels[8])
ax.set_xscale('log')
ax.set_yscale('log')
            
line = ax.scatter(logdata[:,4],logdata[:,8])
    
outputfile = 'mdot_vs_sigma_grav.ps'    
plt.savefig(outputfile, format='ps')

print 'Plotting mdot vs sigma for MRI layer'

figmag = plt.figure()
ax = figmag.add_subplot(111)
ax.set_xlabel(ylabels[5])
ax.set_ylabel(ylabels[9])
ax.set_xscale('log')
ax.set_yscale('log')
            
line = ax.scatter(logdata[:,5],logdata[:,9])
    
outputfile = 'mdot_vs_sigma_mri.ps'
plt.savefig(outputfile, format='ps')
plt.show()
# Finish layer plots

print 'Plotting complete for file ',logfile