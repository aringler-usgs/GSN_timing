#!/usr/bin/env python
from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime, Trace
from scipy import signal
from obspy.signal.cross_correlation import correlate, xcorr_max
import numpy as np
import matplotlib.pyplot as plt
import pickle
import os
client = Client()

import matplotlib as mpl
mpl.rc('font', family = 'serif')
mpl.rc('font', serif = 'Times')
mpl.rc('text', usetex = True)
mpl.rc('font', size = 20)

# def calc_year_sta(year, sta, stime):
#     net, sta = sta.split('.')
#     stime = UTCDateTime(str(year) + '-001T00:00:00')
#     etime = UTCDateTime(str(year+1) + '-001T00:00')
#     ctime = stime
#     times, shifts, vals = [],[], []
#     while ctime < etime:
#         print(ctime)
#         #try:
#         if True:
#             rand = np.random.randint(0,60*60*24)
#             cnt = 0
#             while cnt <= 4:
#                 try:
#                 #if True:
#                     st = client.get_waveforms(net, sta,'*', 'BHZ', ctime + rand, ctime + 60 + rand, attach_response=True)
#                     st.remove_response()
#                     break
#                 except:
#                     cnt += 1
#             fig = plt.figure(1, figsize=(10,14))
#             plt.subplot(3,1,1)
#             plt.plot(st[0].times(), st[0].data*10**6)
#             plt.ylabel('Velocity ($\mu m/s$)')
#             plt.text(-9, np.max(st[0].data*10**6), '(a)')
#             plt.xlim(min(st[0].times()), max(st[0].times()))
#             plt.xlabel('Time (s)')
#             plt.subplot(3,1,2)
#             plt.plot(st[1].times(), st[1].data*10**6)
#             plt.xlim(min(st[1].times()), max(st[1].times()))
#             plt.ylabel('Velocity ($\mu m/s$)')
#             plt.xlabel('Time (s)')
#             plt.text(-9, np.max(st[0].data*10**6), '(b)')
#             st.filter('bandpass', freqmax=1/4., freqmin=1./8.)
#             st.merge(fill_value=0)
#             st.resample(1000)
#             st.sort()
#             print(st)
#             tr1 = st.select(location='00')[0]
#             tr2 = st.select(location='10')[0]
#             cc = correlate(tr1.data, tr2.data, 500)
#             plt.subplot(3,1,3)
#             plt.plot((np.arange(len(cc))-500)/1000., cc)
#             plt.xlim((min((np.arange(len(cc))-500)/1000.), max((np.arange(len(cc))-500)/1000.) ))
#             plt.ylabel('Correlation')
#             plt.xlabel('Lag (s)')
#             plt.text(-.65, 1., '(c)')
#             plt.savefig('example.png', format='PNG')
#             import sys
#             sys.exit()
#             shift, val = xcorr_max(cc)
#             shifts.append(shift)

#             vals.append(val)
#             times.append(str(ctime.year) + ', ' + str(ctime.julday) + ', ' + str(ctime.hour) + ', ' + str(ctime.minute) + str((ctime+rand).second))
#         #except:
#         #    pass
#         ctime += 24*60*60

#     with open(net + '_' + sta + '_' + str(year) + '.pickle2', 'wb') as f:
#         pickle.dump([shifts, vals, times], f)
#     return


fig = plt.figure(1, figsize=(12,16))
net, sta = 'IU', 'SJG'
stime1 = UTCDateTime('2017-150T00:00:00')
stime2 = UTCDateTime('2018-200T00:00:00')
inv1 = client.get_stations(network=net, station=sta,
                                starttime=stime1,
                                endtime=stime1 + 3*24*60*60, level='response', channel='BHZ')

inv2 = client.get_stations(network=net, station=sta,
                                starttime=stime2,
                                endtime=stime2 + 3*24*60*60, level='response', channel='BHZ')


resp1_00 = inv1.get_response(net + '.' + sta + '.00.BHZ', stime1)
resp1_10 = inv1.get_response(net + '.' + sta + '.10.BHZ', stime1)

resp2_00 = inv2.get_response(net + '.' + sta + '.00.BHZ', stime2)
resp2_10 = inv2.get_response(net + '.' + sta + '.10.BHZ', stime2)


tf1_00, f1 = resp1_00.get_evalresp_response(1./20,2**12, start_stage=1, end_stage=1)
tf1_10, f2 = resp1_10.get_evalresp_response(1./40,2**12,start_stage=1, end_stage=1 )

tf2_00, f1 = resp2_00.get_evalresp_response(1./20,2**12,start_stage=1, end_stage=1 )
tf2_10, f2 = resp2_10.get_evalresp_response(1./40,2**12,start_stage=1, end_stage=1 )


plt.subplot(3,1,1)
plt.semilogx(f1, np.angle(tf1_00)*180/np.pi, label=sta + ' 00 BHZ ' +str(stime1.year) + '-' + str(stime1.julday).zfill(3))
plt.semilogx(f2, np.angle(tf1_10)*180/np.pi, label=sta + ' 10 BHZ ' +str(stime1.year) + '-' + str(stime1.julday).zfill(3))
plt.semilogx(f1, np.angle(tf2_00)*180/np.pi, label=sta + ' 00 BHZ ' +str(stime2.year) + '-' + str(stime2.julday).zfill(3))
plt.semilogx(f2, np.angle(tf2_10)*180/np.pi, label=sta + ' 10 BHZ ' +str(stime2.year) + '-' + str(stime2.julday).zfill(3))
plt.text(0.004,88, '(a)')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Phase (degrees)')
plt.xlim((0.01,10.))
plt.legend(ncol=2, loc=9)
plt.subplot(3,1,2)
plt.semilogx(f1, np.angle(tf1_00)*180/np.pi-np.angle(tf2_00)*180/np.pi, label=sta + ' 00 BHZ')
plt.semilogx(f2, np.angle(tf1_10)*180/np.pi-np.angle(tf2_10)*180/np.pi, label=sta + ' 10 BHZ')
#plt.semilogx(f1, np.angle(tf2_00)*180/np.pi, label=net + ' ' + sta + ' 00 BHZ ' +str(stime2.year) + '-' + str(stime2.julday).zfill(3))
#plt.semilogx(f2, np.angle(tf2_10)*180/np.pi, label=net + ' ' + sta + ' 10 BHZ ' +str(stime2.year) + '-' + str(stime2.julday).zfill(3))
plt.xlabel('Frequency (Hz)')
plt.ylabel('Phase Difference (degrees)')
plt.xlim((0.01,10.))
plt.ylim((-10,10))
plt.legend(ncol=2,loc=9)
plt.text(0.004,11, '(b)')

plt.subplot(3,1,3)
pers = [.1, 1, 4,  6,  8, 30]
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
        print('IMPORTANT')
        print(val)
        vals.append(shift)


    plt.plot(phaseangs, vals, '.', label=str(per) + ' s ', alpha=0.7)
plt.legend(loc=9, ncol=5)
plt.xlabel('Phase Response Error (degrees)')
plt.ylabel('Apparent Timing Error (ms)')
plt.xlim((-5,5))
plt.text(-6.3, 60, '(c)')









st1 = client.get_waveforms(net, sta,'*', 'BHZ', stime1, stime1 + 60)
st1.remove_response(inventory=inv1)
st1.filter('bandpass', freqmax=1/4., freqmin=1./8.)
st1.merge(fill_value=0)
st1.resample(1000)
st1.sort()
tr1 = st1.select(location='00')[0]
tr2 = st1.select(location='10')[0]
print(st1)
cc = correlate(tr1.data, tr2.data, 500)
shift, val = xcorr_max(cc)
print(shift)
st2 = client.get_waveforms(net, sta,'*', 'BHZ', stime2, stime2 + 60, attach_response=True)
st2.remove_response(inventory=inv2)
st2.filter('bandpass', freqmax=1/4., freqmin=1./8.)
st2.merge(fill_value=0)
st2.resample(1000)
st2.sort()
tr1 = st2.select(location='00')[0]
tr2 = st2.select(location='10')[0]
print(st2)
cc = correlate(tr1.data, tr2.data, 500)
shift, val = xcorr_max(cc)
print(shift)
plt.savefig('Plot_' + sta + 'TEST.png')
plt.savefig('figure9.pdf', format='PDF', dpi=400)









