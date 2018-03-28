# Written 28/3/18 by dh4gan
# Reads all planets files and computes semimajor axis evolution

import io_disc as io
import matplotlib.pyplot as plt

prefix = raw_input('What is the file prefix? ')

time, active,ap, mp = io.obtain_planet_tracks(prefix)

print time,ap

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(time,ap[0,:])
ax1.set_ylabel('Semimajor axis (AU)',fontsize=22)
ax1.set_xlabel('Time (yr)',fontsize=22)
plt.show()
