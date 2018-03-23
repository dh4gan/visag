# Module contains dictionaries etc to help with plotting routines

import numpy as np
from multigraph import multigraph, multigraph_legend,multigraph_legend_points

nprofcol = 11
nlayercol = 11
nlogcol = 11

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

# Now generate plots


def read_profile(profilefile):
    '''Reads profile data from file'''

    f = open(profilefile,'r')
    line = f.readline()
    arr = np.fromstring(line.strip(), dtype=float, sep=" ")
    
    f.close()
    time = arr[0]
    ngrid = int(arr[1])
    
    print 'Reading file ',profilefile
    print 'Time: 'str(time)+ ' yr'
    
    profdata = np.genfromtxt(profilefile,skip_header=1)
    profdata.reshape(profdata.size/nprofcol,nprofcol)

    return time, profdata


def read_planets(planetfile):
    '''Reads planetary data from file'''

    f = open(planetsfile, 'r')
    line = f.readline()
    arr = np.fromstring(line.strip(), dtype=int, sep=" ")
    
    time = arr[0]
    nplanet = arr[1]
    nactive = arr[2]
    
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

    return nplanet,nactive, active,mp,ap


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
    
    nplanet,nactive,active,mp,ap = read_planets(planetfile)
    
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
    xpoints = np.zeros(nprofcol)
    ypoints = np.zeros(nprofcol)
        
    for i in range(nprofcol):
        if ymin[i]!=ymax[i]:
            ypoints[i] = 2.0*ymin[i]
        else:
            ypoints[i] = 2.0*np.min(profdata[:,i])

    multigraph_legend_points(profdata,nprofcol,profilexlabel,profilelabels,profileylog,profileymin,profileymax,profileoutputstring,legendstring,ap,ypoints,mp)

    return time,profdata,nplanet,nactive,active,mp,ap



