#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from obspy.signal.cross_correlation import correlate, xcorr_max
from obspy.core import UTCDateTime, Trace
import scipy

import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times')
mpl.rc('text', usetex=True)
mpl.rc('font',size=16)

shifts, snrs, vals, vals2, shifts2 = [],[],[],[],[]

fig = plt.figure(1, figsize=(8,12))
for n in range(5,40):
    snrs.append(n)
    t = np.arange(0, 60, 1./40.)

    # 1 Hz sine wave
    sig1 = np.cos(2*np.pi*t*1.)
    snr = 10.**(n/10.)
    var = sig1.var()

    sig1 += np.sqrt(var/snr)*np.random.standard_normal(len(sig1))


    sig2 = np.cos(2*np.pi*(t+3./1000.)*1.)
    sig2 += np.sqrt(var/snr)*np.random.standard_normal(len(sig2))




    stats1 = {'network': 'XX', 'station': 'SIG1', 'location': '',
             'channel': 'BHZ', 'npts': len(sig1), 'sampling_rate': 40.,
             'mseed': {'dataquality': 'D'}}

    stats1['starttime'] = UTCDateTime()
    stats2 = {'network': 'XX', 'station': 'SIG2', 'location': '',
             'channel': 'BHZ', 'npts': len(sig2), 'sampling_rate': 40.,
             'mseed': {'dataquality': 'D'}}

    stats2['starttime'] = UTCDateTime()
    tr2 = Trace(data=sig2, header=stats2)
    tr1 = Trace(data=sig1, header=stats1)


    cc1 = correlate(sig1, sig2, 50)
    shift, val = xcorr_max(cc1)

    shifts.append(shift)
    vals.append(val)

    tr1.resample(1000.)
    tr2.resample(1000.)

    cc2 = correlate(tr1.data, tr2.data, 500)
    shift, val = xcorr_max(cc2)
    print(shift)
    #if shift > 1./1000.:
    #    shift = 1./1000.
    #if shift < -1./1000.:
    #    shift = -1./1000.
    print(val)
    shifts2.append(shift)
    vals2.append(val)
vals = np.array(vals)
vals2 = np.array(vals2)
shifts2 = np.array(shifts2)
shifts = np.array(shifts)



# fig = plt.figure(1)
# plt.plot(tr1.data)
# plt.plot(tr2.data)
# plt.show()

# import sys
# sys.exit()



fig = plt.figure(1, figsize=(12,8))
plt.subplot(2,1,1)
plt.plot(snrs, 1000*shifts/40., '.',label='1 Hz Sine Wave at 40 sps')
plt.plot(snrs, 1000*shifts2/1000., '.', label='1 Hz Sine Wave at 1000 sps')
plt.xlabel('Signal-to-Noise Ratio (dB)')
plt.ylabel('Maximum Lag (ms)')
plt.ylim(-0.01*1000, 0.01*1000)
plt.text(-2, 10, '(a)')
plt.legend(loc=4)
plt.subplot(2,1,2)
plt.plot(snrs, vals, '.',label='1 Hz Sine Wave at 40 sps')
plt.plot(snrs, vals2, '.', label='1 Hz Sine Wave at 1000 sps')
plt.xlabel('Signal-to-Noise Ratio (dB)')
plt.ylabel('Correlation (Normalized)')
plt.legend(loc=4)
plt.text(-2, 1, '(b)')
plt.ylim(0.8, 1.01)
plt.savefig('TEST.png', format='PNG', dpi=400)
plt.savefig('TEST.pdf', format='PDF', dpi=400)


