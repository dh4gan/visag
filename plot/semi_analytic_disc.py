#!/usr/local/bin/python

# Python script plots data from semi analytic disc runs

#import matplotlib.pyplot as plt
import numpy as np
from multigraph import multigraph, multigraph_legend,multigraph_legend_points

nprofcol=11
nlayercol=11
nlogcol = 11

# Input file data

prefix = raw_input('What is the file prefix? ')
fileno = input('Which dump is to be plotted? ')
layerchoice = raw_input('Plot layer data? (y/n) ')
planetchoice = raw_input('Plot planet data? (y/n) ')

fileno = int(fileno)

use_layer = True
use_planet = True

if 'y' not in layerchoice: use_layer = False
if 'y' not in planetchoice: use_planet = False

# Create filenames

profile = prefix+'profile.'+str(fileno)
spectrum = prefix+'spectrum.'+str(fileno)
layer = prefix+'layer.'+str(fileno)
planets = prefix+'planets.'+str(fileno)

logfile = prefix+'.log'

# Open log file, to find time of dump

logdata = np.genfromtxt(logfile)

time = logdata[fileno-1,0]

print 'Time is ',time ,' years'


# If using planet data, load this into file

if use_planet:
    f = open(planets, 'r')
    line = f.readline()
    arr = np.fromstring(line.strip(), dtype=int, sep=" ")

    nplanet = arr[0]
    nactive = arr[1]

    print 'Number of planets: ',nplanet
    print 'Those of which are active: ',nactive
    
    active = np.zeros(nplanet)
    mp = np.zeros(nplanet)
    ap = np.zeros(nplanet)

    for i in range(nplanet):
        line = f.readline()
        arr = np.fromstring(line.strip(), sep=" ")
        
        active[i] = arr[0]
        mp[i] = arr[1]
        ap[i] = arr[2]

        print 'Planet ',i+1,', mass, semimajor axis: ', mp[i],ap[i]
##################################
# Generate profile data plots
##################################

# Open profile file and read

profdata = np.genfromtxt(profile,skip_header=1)
profdata.reshape(profdata.size/nprofcol,nprofcol)

# Set up plot data for multigraph function

# Filenames
outputstring=['r','sigma','cs','kappa','gamma','mu', 'T', 'tau', 'nu','alpha', 'Q']

for i in range(len(outputstring)):
    outputstring[i] = prefix+outputstring[i]

# Plot labels
ylabels = [r'r (AU)',r'$\Sigma$ (g cm $^{-2}$)',r' $c_s$ (cm s$^{-1}$)',r'$\kappa$ (cm$^{2}$ g$^{-1}$)',
           r'$\gamma$',r'$\mu$',r'$T_c$ (K)',r'$ \tau $',r'$\nu_g$',r'$\alpha_{g}$',r'$ Q $']
xlabel = ylabels[0]

legendstring=[]
# Legend Label
for i in range(nprofcol):
    legendstring.append('Time = '+str(time)+' yr')

# Log the y axis? True/False
setylog=[True, True, True,True,False,False,True,True,True,True,True]

# y limits - set defaults first
ymin=[]
ymax=[]
for i in range(nprofcol):
    ymin.append(0.0)
    ymax.append(0.0)

# Now define any non-default limits

# sigma
ymin[1] = 1.0e1
ymax[1] = 1.0e6

# alpha
ymin[9] = 1.0e-5
ymax[9] = 1.0e0

# Q
ymin[10] = 1.0e0
ymax[10] = 1.0e2

# Now generate plots

# If no planets, use the simple legend function: otherwise use function with points



if use_planet:
    xpoints = np.zeros(nprofcol)
    ypoints = np.zeros(nprofcol)

    for i in range(nprofcol):        
        if ymin[i]!=ymax[i]: 
            ypoints[i] = 2.0*ymin[i]
        else:
            ypoints[i] = 2.0*np.min(profdata[:,i])
    
    multigraph_legend_points(profdata,nprofcol,xlabel,ylabels,setylog,ymin,ymax,outputstring,legendstring,ap,ypoints,mp)
else:
    multigraph_legend(profdata,nprofcol,xlabel,ylabels,setylog,ymin,ymax,outputstring,legendstring)


if use_layer:
    ###############################
    # Generate layer plots
    ###############################

    # Open layer file and read

    layerdata = np.genfromtxt(layer,skip_header=1)
    layerdata.reshape(layerdata.size/nlayercol,nlayercol)

    # Set up plot data for multigraph function

    # Filenames
    outputstring=['r','sigma','sigma_m', 'sigma_tot', 'cs_m','kappa_m','gamma_m','mu_m', 'tau_m', 'nu_m', 'alpha_m']

    for i in range(len(outputstring)):
        outputstring[i] = prefix+outputstring[i]

    # Plot labels
    ylabels = [r'r (AU)',r'$\Sigma_g$ (g cm $^{-2}$)',r'$\Sigma_m$ (g cm $^{-2}$)',r'$\Sigma_{tot}$ (g cm $^{-2}$)',
           r' $c_{s,m}$ (cm s$^{-1}$)',r'$\kappa_m$ (cm$^{2}$ g$^{-1}$)',r'$\gamma_m$',r'$\mu_m$',r'$ \tau_m $',r'$\nu_{m}$',r'$ \alpha_m $']
    xlabel= ylabels[0]


    # Legends

    legendstring=[]
    for i in range(nlayercol):
        legendstring.append('Time = '+str(time)+' yr')

    # Log the y axis? True/False
    setylog=[True, True, True,True,True,True,False,False,True,True,True]

    # y limits - set defaults first
    ymin=[]
    ymax=[]
    for i in range(nlayercol):
        ymin.append(0.0)
        ymax.append(0.0)


    # Now define any non-default limits

    # Now generate layer plots

    multigraph_legend(layerdata,nlayercol,xlabel,ylabels,setylog,ymin,ymax,outputstring,legendstring)

###################################
# Plot log data
####################################


outputstring = ['t', 'dt','mgrav','mmag','mdisc','grav_max','mag_max','sig_max','tot_lum','mdot_grav', 'mdot_mag']

for i in range(len(outputstring)):
    outputstring[i] = prefix+outputstring[i]

# Plot labels
ylabels = [r't (yr)',r'dt (yr)',r'$M_{grav}$ ($M_{\odot}$)',r'$M_{MRI}$ ($M_{\odot}$)',r'$M_{disc}$ ($M_{\odot}$)',
           r' $\Sigma_{grav,max}$ (g cm$^{-2}$)',r' $\Sigma_{MRI,max}$ (g cm$^{-2}$)',r' $\Sigma_{max}$ (g cm$^{-2}$)',
           r'$L_{tot}$ ($L_{\odot}$)', r'$\dot{M}_{grav}$ ($M_{\odot} yr^{-1}$)', r'$\dot{M}_{mag}$ $(M_{\odot} yr^{-1}$)']

xlabel = ylabels[0]

# Log the y axis? True/False
setylog=[False, False, False,True,True,True,True,True,True,True,True]

# y limits - set defaults first
ymin=[]
ymax=[]
for i in range(nlogcol):
    ymin.append(0.0)
    ymax.append(0.0)

# Set any non-default limits
    
multigraph(logdata,nlogcol,xlabel,ylabels,setylog,ymin,ymax,outputstring)   



print 'Plotting complete for file ',logfile