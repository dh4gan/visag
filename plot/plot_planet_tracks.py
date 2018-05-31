# Written 28/3/18 by dh4gan
# Reads all planets files and computes semimajor axis evolution

import io_disc as io
import filefinder as ff
import matplotlib.pyplot as plt
from numpy import amax,amin

prefix = ff.get_file_prefix('*.log')

time, nplanet,active,ap, mp,tmig = io.obtain_planet_tracks(prefix)

# Semimajor axis
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)

ax1.set_ylim(amin(ap[:,:])*0.9,amax(ap[:,:])*1.1)

for i in range(nplanet):
    ax1.plot(time,ap[i,:], label='Planet '+str(i+1))

ax1.set_ylabel('Semimajor axis (AU)',fontsize=22)
ax1.set_xlabel('Time (yr)',fontsize=22)
ax1.legend()



# Migration Timescale
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

ax2.set_ylim(amin(tmig[:,:])*0.9,amax(tmig[:,:])*1.1)
ax2.set_yscale('log')

for i in range(nplanet):
    ax2.plot(time,tmig[i,:], label='Planet '+str(i+1))

ax2.set_ylabel('Migration Timescale (yr)',fontsize=22)
ax2.set_xlabel('Time (yr)',fontsize=22)
ax2.legend()


plt.show()
