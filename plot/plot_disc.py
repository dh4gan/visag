#!/usr/local/bin/python

# Python script plots profile data from semi_analytic_disc runs

#import matplotlib.pyplot as plt
import numpy as np
import io_disc as io
import filefinder as ff
from multigraph import multigraph, multigraph_legend,multigraph_legend_points

nlogcol = 11
# Input file data

prefix = raw_input('What is the file prefix? ')

profilefile = ff.find_sorted_local_input_files(prefix+'*profile*')

# Get filenumber for plotting
fileno = profilefile.rsplit('.',1)[1]

layerchoice = raw_input('Plot layer data? (y/n) ')
planetchoice = raw_input('Plot planet data? (y/n) ')

use_layer = True
use_planet = True

if 'y' not in layerchoice: use_layer = False
if 'y' not in planetchoice: use_planet = False

# Create filenames

spectrumfile = prefix+'_spectrum.'+fileno
layerfile = prefix+'_layer.'+fileno
planetfile = prefix+'_planets.'+fileno
logfile = prefix+'.log'

ifile = int(fileno)

##################################
# Generate profile data plots
##################################

if use_planet:
    time,profdata,nplanet,nactive,active,mp,ap = io.plot_profile_data_planets(profilefile,planetfile)
else:
    time,profdata = io.plot_profile_data(profilefile)




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