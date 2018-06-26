#
# io_disc.py
# Contains helpful dictionaries and other global variables for plotting 
# Also contains functions to read output files and plot data
#

import numpy as np
import matplotlib.pyplot as plt
import filefinder as ff
from multigraph import multigraph, multigraph_legend,multigraph_legend_points

nprofcol = 11
nlayercol = 11
nlogcol = 11


# Variables for profile file

profilekeys = ['r','sigma','cs','kappa','gamma','mu', 'T', 'tau', 'nu','alpha', 'Q']

profilelabels = [r'r (AU)',r'$\Sigma$ (g cm $^{-2}$)',r' $c_s$ (cm s$^{-1}$)',r'$\kappa$ (cm$^{2}$ g$^{-1}$)',
           r'$\gamma$',r'$\mu$',r'$T_c$ (K)',r'$ \tau $',r'$\nu_g$',r'$\alpha_{g}$',r'$ Q $']
profilexlabel = profilelabels[0]

# Log the y axis? True/False
profileylog=[True, True, True,True,False,False,True,True,True,True,True]

# y limits - set defaults first
profileymin=[]
profileymax=[]

for i in range(nprofcol):
    profileymin.append(0.0)
    profileymax.append(0.0)

# Now define any non-default limits

# sigma
profileymin[1] = 1.0e1
profileymax[1] = 1.0e6

#cs
profileymin[2] = 1.0e1
profileymax[2] = 1.0e7

#kappa
profileymin[3] = 1.0e-10
profileymax[3] = 1.0e5

# alpha
profileymin[9] = 1.0e-5
profileymax[9] = 1.0e0

# Q
profileymin[10] = 1.0e0
profileymax[10] = 1.0e2


# Variables for log data

logkeys = ['t', 'dt','mdisc', 'tot_lumin', 'sig_max', 'mgrav', 'mmag', 'grav_max', 'mag_max', 'mdot_grav', 'mdot_mag','mdot_wind' ]


loglabels = [r't (yr)', r'dt (yr)', r'$M_{disc}$ ($M_{\odot}$)',
             r'$L_{tot}$ ($L_{\odot}$)',r' $\Sigma_{grav,max}$ (g cm$^{-2}$)',
             r'$M_{grav}$ ($M_{\odot}$)',r'$M_{MRI}$ ($M_{\odot}$)',
             r' $\Sigma_{grav,max}$ (g cm$^{-2}$)', r' $\Sigma_{MRI,max}$ (g cm$^{-2}$)',
             r'$\dot{M}_{grav}$ ($M_{\odot} yr^{-1}$)', r'$\dot{M}_{mag}$ $(M_{\odot} yr^{-1}$)',r'$\dot{M}_{wind}$ $(M_{\odot} yr^{-1}$)']

logxlabel = loglabels[0]

# Log the y axis? True/False
logylog=[False, False, False,True,True,True,True,True,True,True,True,True]

# y limits - set defaults first
logymin=[]
logymax=[]
for i in range(nlogcol):
    logymin.append(0.0)
    logymax.append(0.0)


# Some global variables for planet plotting
mearth = 0.00314 # Earth mass in Jupiter masses

# Minimum and maximum planet sizes (in pixels)
minplanetsize = 50
maxplanetsize = 200

# Colours for different classes
earthcolour = '#0099ff'
neptunecolour = '#60b28c'
jupitercolour = 'red'
BDcolour = '#663300'

################################
# END OF VARIABLE DEFINITIONS
################################


def read_profile(profilefile):
    '''Reads profile data from file'''

    f = open(profilefile,'r')
    line = f.readline()
    arr = np.fromstring(line.strip(), dtype=float, sep=" ")
    
    f.close()
    time = arr[0]
    ngrid = int(arr[1])    
    
    print 'Reading file ',profilefile
    print 'Time: '+str(time)+ ' yr'
    
    profdata = np.genfromtxt(profilefile,skip_header=1)
    profdata.reshape(profdata.size/nprofcol,nprofcol)

    return time, profdata


def read_planets(planetfile,verbose=True):
    '''Reads planetary data from file'''

    f = open(planetfile, 'r')
    line = f.readline()
    arr = np.fromstring(line.strip(), dtype=float, sep=" ")
    
    time = arr[0]
    nplanet = int(arr[1])
    nactive = int(arr[2])
    
    if(verbose):
        print 'Number of planets: ',nplanet
        print 'Those of which are active: ',nactive
    
    active = np.zeros(nplanet)
    mp = np.zeros(nplanet)
    ap = np.zeros(nplanet)
    tmig = np.zeros(nplanet)
    
    for i in range(nplanet):
        line = f.readline()
        arr = np.fromstring(line.strip(), sep=" ")
        
        active[i] = arr[0]
        mp[i] = arr[1]
        ap[i] = arr[2]
        tmig[i] = arr[3]
    
        if(verbose):
            print active[i],mp[i],ap[i],tmig[i]

    return time,nplanet,nactive, active,mp,ap,tmig


def get_planet_size_and_colour(nplanet,mp):
    '''Given a planet mass, returns a colour and size'''

    planetcolours = []

    planetsizes = 100*mp[:]

    for i in range(nplanet):
        if(planetsizes[i]<minplanetsize):
            planetsizes[i]=minplanetsize
        if(planetsizes[i]>maxplanetsize):
            planetsizes[i]=maxplanetsize

    for i in range(nplanet):
        if(mp[i]<1.5*mearth):
            planetcolours.append(earthcolour)
        elif(mp[i]>1.5*mearth and mp[i]<5*mearth):
            planetcolours.append(neptunecolour)
        elif(mp[i]>5*mearth and mp[i]<13.0):
            planetcolours.append(jupitercolour)
        else:
            planetcolours.append(BDcolour)

    return planetsizes,planetcolours

def read_log(logfile):
    '''Reads the .log file'''
    return np.genfromtxt(logfile)


def plot_profile_multifiles_variable(prefix, add_planets=False):
    '''Reads multiple profile files and plots a specific variable'''

    filenames = ff.find_sorted_local_input_fileset(prefix+"*profile*")
    nfiles = len(filenames)

    initial = input('Starting filenumber? ')
    final = input('Final filenumber? ')

    print 'Now select variable to plot: here are the choices'

    for i in range(len(profilekeys)):
        print str(i+1)+': '+profilekeys[i]

    var = input('Which variable (1-'+str(len(profilekeys))+')? ')

    var = var-1

    if(final>nfiles):
        print "Limiting count to available files"
        final = nfiles

    fig1 = plt.figure()
    ax = fig1.add_subplot(111)

    if(add_planets):
        planetfiles = ff.find_sorted_local_input_fileset(prefix+"*planets*")

    initial = initial-1
    final = final - 1
    for i in range(initial,final):
        
        time,profdata = read_profile(filenames[i])
        
        if(add_planets):
            t,nplanet,nactive, active,mp,ap,tmig = read_planets(planetfiles[i],verbose=False)

            # Setup planet points for plotting
            xpoints = np.zeros(nplanet)
            ypoints = np.zeros(nplanet)

        
            if profileymin[var]!=profileymax[var]:
                ypoints[:] = 2.0*profileymin[var]
            else:
                ypoints[:] = 2.0*np.min(profdata[:,var])

         
            planetsizes,planetcolours = get_planet_size_and_colour(nplanet,mp)                

        if(profileylog[var]):
            ax.set_xscale('log')
            ax.set_yscale('log')
        if(profileymin[var]!=profileymax[var]):
            ax.set_ylim(profileymin[var],profileymax[var])

        line1 = ax.plot(profdata[:,0],profdata[:,var])
        ax.set_xlabel(profilexlabel)
        ax.set_ylabel(profilelabels[var])
        ax.text(0.9, 0.9,'t = '+str(np.round(time,2))+' yr',
                bbox=dict(edgecolor='black',facecolor='none'), horizontalalignment='center',
                verticalalignment='center',transform = ax.transAxes)

        if(add_planets):
            ax.scatter(ap,ypoints,s=planetsizes,facecolor=planetcolours)

        outputfile =profilekeys[var]+'_'+filenames[i]+'.png'

        print 'Saving to ',outputfile
        plt.savefig(outputfile, format='png')
        ax.clear()


def plot_profile_data(profilefile):
    '''Reads a given profile file and plots all variables at that snapshot'''

    time,profdata = read_profile(profilefile)
    
    # Set up plot data for multigraph function

    # Filenames
    profileoutputstring=[]

    for i in range(nprofcol):
        profileoutputstring.append(profilekeys[i]+'_'+profilefile)

    legendstring=[]

    # Legend Label
    for i in range(nprofcol):
        legendstring.append('t = '+str(np.round(time,2))+' yr',)


    multigraph_legend(profdata,nprofcol,profilexlabel,profilelabels,profileylog,profileymin,profileymax,profileoutputstring,legendstring)

    return time,profdata

def plot_profile_data_planets(profilefile,planetfile):
    '''Reads a given profile file and plots all variables (and planets) at that snapshot'''

    time,profdata = read_profile(profilefile)
    
    t,nplanet,nactive,active,mp,ap,tmig = read_planets(planetfile)

    if(np.abs(time-t)>1.0e-30):
        print "Warning: times of profile/planet files don't match"
    
    # Set up plot data for multigraph function
    
    # Filenames
    profileoutputstring=[]

    for i in range(nprofcol):
        profileoutputstring.append(profilekeys[i]+'_'+profilefile)
    
    legendstring=[]
    
    # Legend Label
    for i in range(nprofcol):
        legendstring.append('t = '+str(np.round(time,2))+' yr',)

    # Setup planet points for plotting
    xpoints = np.zeros((nplanet,nprofcol))
    ypoints = np.zeros((nplanet,nprofcol))
   
        
    for i in range(nprofcol):
        if profileymin[i]!=profileymax[i]:
            ypoints[:,i] = 2.0*profileymin[i]
        else:
            ypoints[:,i] = 2.0*np.min(profdata[:,i])

    planetsizes,planetcolours = get_planet_size_and_colour(nplanet,mp)

    print ap, ypoints, planetcolours
    multigraph_legend_points(profdata,nprofcol,profilexlabel,profilelabels,profileylog,profileymin,profileymax,profileoutputstring,legendstring,ap,ypoints,planetsizes,planetcolours)

    return time,profdata,nplanet,nactive,active,mp,ap


def plot_log_data(logfile):
    ''' Plots log file data'''

    logdata = read_log(logfile)
    
    logoutputstring = []
    for i in range(nprofcol):
        logoutputstring.append(logkeys[i]+'_'+logfile)

    multigraph(logdata,nlogcol,logxlabel,loglabels,logylog,logymin,logymax,logoutputstring)

    return logdata


def obtain_planet_tracks(prefix):
    '''Reads all planetary data and creates tracks for each planet'''
    planetfiles = ff.find_sorted_local_input_fileset(prefix+"*planets*")

    time,nplanet,nactive,active,mp,ap,tmig = read_planets(planetfiles[1])

    time_all = np.zeros(len(planetfiles)-1)
    active_all = np.zeros((nplanet,len(planetfiles)-1))
    ap_all = np.zeros((nplanet,len(planetfiles)-1))
    mp_all = np.zeros((nplanet,len(planetfiles)-1))
    tmig_all= np.zeros((nplanet,len(planetfiles)-1))
    
    for i in range(len(planetfiles)-1):
        
        time_all[i],nplanet,nactive,active, mp,ap,tmig = read_planets(planetfiles[i+1],verbose=False)
        active_all[:,i] = active[:]
        mp_all[:,i] = mp[:]
        ap_all[:,i] = ap[:]
        tmig_all[:,i] = tmig[:]

    return time_all, nplanet,active_all,ap_all, mp_all,tmig_all

        
