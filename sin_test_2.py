#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from obspy.signal.cross_correlation import correlate, xcorr_max
from obspy.core import UTCDateTime, Trace

import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times')
mpl.rc('text', usetex=True)
mpl.rc('font',size=14)

fig = plt.figure(1, figsize=(8,12))

pers = [4, 5, 6, 7, 8]
for per in pers:
    phaseangs = np.arange(-5,5,.1)
    vals = []
    for phaserr in phaseangs:

        t = np.arange(0, 60, 1./40.)

        # 1 Hz sine wave
        sig1 = np.cos(2*np.pi*t*(1/per))

        sig2 = np.cos(2*np.pi*(t)*(1/per) + phaserr*np.pi/180.)

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
        print(shift)
        print(val)

        tr1.resample(1000.)
        tr2.resample(1000.)

        cc2 = correlate(tr1.data, tr2.data, 500)
        shift, val = xcorr_max(cc2)
        print(shift)
        print(val)
        vals.append(shift)


    plt.plot(phaseangs, vals, '.', label=str(per) + ' s ', alpha=0.7)
plt.legend()
plt.xlabel('Phase Response Error (degrees)')
plt.ylabel('Apparent Timing Error (ms)')
plt.xlim((-5,5))
plt.savefig('Error.png')


# fig = plt.figure(1, figsize=(8,12))
# plt.subplot(3,1,1)
# plt.plot(t, sig1, label='Initial Signal')
# plt.plot(t, sig2, '.', alpha=0.4, label='2 ms Timing Error')
# plt.plot(t, sig1-sig2, label='Difference')
# plt.xlabel('Time (s)')
# plt.ylabel('Amplitude (Normalized)')
# plt.xlim((min(t), max(t)))
# plt.text(-7, 1., '(a)')
# plt.legend(loc=2)
# plt.subplot(3,1,2)
# plt.plot(tr1.times(), tr1.data, label='Initial Signal at 1000 Hz')
# plt.plot(tr1.times(), tr2.data, '.',ms =2., alpha=0.2, label='2 ms Timing Error at 1000 Hz')
# plt.plot(tr1.times(), tr1.data-tr2.data, label='Difference at 1000 Hz')
# plt.xlabel('Time (s)')
# plt.xlim((min(t), max(t)))
# plt.text(-7, 1., '(b)')
# plt.legend(loc=2)
# plt.ylabel('Amplitude (Normalized)')
# plt.subplot(3,1,3)
# plt.plot((np.arange(len(cc2))-500)/1000., cc2, label='Cross-Correlation 1000 Hz')
# plt.plot([3./1000., 3./1000.], [-1.,1.], label='Maximum Lag')
# plt.plot([0./1000., 0./1000.], [-1.,1.], label='Zero Lag')
# plt.xlim((-0.2,0.2))
# plt.ylim((-1.,1.))
# plt.text(-.245, 1., '(c)')
# #plt.plot(cc1-cc2, label='Difference')
# plt.xlabel('Lag (s)')
# plt.ylabel('Correlation (Normalized)')
# plt.legend(loc=2)
# plt.savefig('Response_Phase.pdf', format='PDF')
# plt.savefig('Response_Phase.png', format='PNG')