# Module contains dictionaries etc to help with plotting routines

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

# alpha
profileymin[9] = 1.0e-5
profileymax[9] = 1.0e0

# Q
profileymin[10] = 1.0e0
profileymax[10] = 1.0e2


# Variables for log data

logkeys = ['t', 'dt','mdisc', 'tot_lumin', 'sig_max', 'mgrav', 'mmag', 'grav_max', 'mag_max', 'mdot_grav', 'mdot_mag' ]


loglabels = [r't (yr)', r'dt (yr)', r'$M_{disc}$ ($M_{\odot}$)',
             r'$L_{tot}$ ($L_{\odot}$)',r' $\Sigma_{grav,max}$ (g cm$^{-2}$)',
             r'$M_{grav}$ ($M_{\odot}$)',r'$M_{MRI}$ ($M_{\odot}$)',
             r' $\Sigma_{grav,max}$ (g cm$^{-2}$)', r' $\Sigma_{MRI,max}$ (g cm$^{-2}$)',
             r'$\dot{M}_{grav}$ ($M_{\odot} yr^{-1}$)', r'$\dot{M}_{mag}$ $(M_{\odot} yr^{-1}$)']

logxlabel = loglabels[0]

# Log the y axis? True/False
logylog=[False, False, False,True,True,True,True,True,True,True,True]

# y limits - set defaults first
logymin=[]
logymax=[]
for i in range(nlogcol):
    logymin.append(0.0)
    logymax.append(0.0)



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
    
    for i in range(nplanet):
        line = f.readline()
        arr = np.fromstring(line.strip(), sep=" ")
        
        active[i] = arr[0]
        mp[i] = arr[1]
        ap[i] = arr[2]
    
        if(verbose):
            print active[i],mp[i],ap[i]

    return time,nplanet,nactive, active,mp,ap


def read_log(logfile):
    return np.genfromtxt(logfile)


def plot_profile_multifiles_variable(prefix, add_planets=False):

    filenames = ff.find_sorted_local_input_fileset(prefix+"*profile*")
    nfiles = len(filenames)

    print 'Now select variable to plot: here are the choices'

    for i in range(len(profilekeys)):
        print str(i+1)+': '+profilekeys[i]

    var = input('Which variable to plot?')

    var = var-1

    initial = input('Starting filenumber? ')
    final = input('Final filenumber? ')

    if(final>nfiles):
        print "Limiting count to available files"
        final = nfiles

    fig1 = plt.figure()
    ax = fig1.add_subplot(111)

    if(add_planets):
        planetfiles = ff.find_sorted_local_input_fileset(prefix+"*planets*")

    for i in range(initial,final):
        
        time,profdata = read_profile(filenames[i])
        
        if(add_planets):
            t,nplanet,nactive, active,mp,ap = read_planets(planetfiles[i],verbose=False)

            # Setup planet points for plotting
            xpoints = np.zeros(nplanet)
            ypoints = np.zeros(nplanet)
        
            if profileymin[var]!=profileymax[var]:
                ypoints[:] = 2.0*profileymin[var]
            else:
                ypoints[:] = 2.0*np.min(profdata[:,var])

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
            ax.scatter(ap,ypoints,s=10*mp,facecolor='red')

        outputfile =profilekeys[var]+'_'+filenames[i]+'.png'

        print 'Saving to ',outputfile
        plt.savefig(outputfile, format='png')
        ax.clear()



def plot_profile_data(profilefile):

    time,profdata = read_profile(profilefile)
    
    # Set up plot data for multigraph function

    # Filenames
    profileoutputstring=[]

    for i in range(nprofcol):
        profileoutputstring.append(profilekeys[i]+'_'+profilefile)

    legendstring=[]

    # Legend Label
    for i in range(nprofcol):
        legendstring.append('Time = '+str(time)+' yr')


    multigraph_legend(profdata,nprofcol,profilexlabel,profilelabels,profileylog,profileymin,profileymax,profileoutputstring,legendstring)

    return time,profdata

def plot_profile_data_planets(profilefile,planetfile):
    
    time,profdata = read_profile(profilefile)
    
    t,nplanet,nactive,active,mp,ap = read_planets(planetfile)

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
        legendstring.append('Time = '+str(time)+' yr')

    # Setup planet points for plotting
    xpoints = np.zeros((nplanet,nprofcol))
    ypoints = np.zeros((nplanet,nprofcol))
        
    for i in range(nprofcol):
        if profileymin[i]!=profileymax[i]:
            ypoints[:,i] = 2.0*profileymin[i]
        else:
            ypoints[:,i] = 2.0*np.min(profdata[:,i])


    print ap, ypoints
    multigraph_legend_points(profdata,nprofcol,profilexlabel,profilelabels,profileylog,profileymin,profileymax,profileoutputstring,legendstring,ap,ypoints,mp)

    return time,profdata,nplanet,nactive,active,mp,ap


def plot_log_data(logfile):

    logdata = read_log(logfile)
    
    logoutputstring = []
    for i in range(nprofcol):
        logoutputstring.append(logkeys[i]+'_'+logfile)

    multigraph(logdata,nlogcol,logxlabel,loglabels,logylog,logymin,logymax,logoutputstring)

    return logdata


def obtain_planet_tracks(prefix):

    planetfiles = ff.find_sorted_local_input_fileset(prefix+"*planets*")

    time,nplanet,nactive,active,mp,ap = read_planets(planetfiles[0])

    time_all = np.zeros(len(planetfiles))
    active_all = np.zeros((nplanet,len(planetfiles)))
    ap_all = np.zeros((nplanet,len(planetfiles)))
    mp_all = np.zeros((nplanet,len(planetfiles)))
    
    for i in range(len(planetfiles)):
        
        time_all[i],nplanet,nactive,active, mp,ap = read_planets(planetfiles[i],verbose=False)
        active_all[:,i] = active[:]
        mp_all[:,i] = mp[:]
        ap_all[:,i] = ap[:]

    return time_all, active_all,ap_all, mp_all

        
