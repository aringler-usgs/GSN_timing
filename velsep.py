#!/usr/bin/env python
import matplotlib.pyplot as plt
import pickle
import numpy as np


import matplotlib as mpl
mpl.rc('font', family = 'serif')
mpl.rc('font', serif = 'Times')
mpl.rc('text', usetex = True)
mpl.rc('font', size = 20)


vels = np.arange(3.5,6,1/100.)

seps = np.arange(0,0.3, 1./1000.)


x, y = np.meshgrid(vels, seps)

Z = (y/x)*1000.


fig =plt.figure(1, figsize=(9,9))
plt.pcolormesh(x, y, Z)
plt.ylabel('Sensor Separation (km)')
plt.xlabel('Apparent Seismic Velocity (km/s)')
plt.colorbar(orientation='horizontal', label='Error in timing (ms)')
plt.savefig('sep.png', format='PNG')
plt.savefig('sep.pdf', format='PDF')