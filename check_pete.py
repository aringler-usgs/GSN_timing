#!/usr/bin/env python
from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime
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

def calc_year_sta(year, sta):
    net, sta = sta.split('.')
    stime = UTCDateTime(str(year) + '-001T00:00:00')
    etime = UTCDateTime(str(year+1) + '-001T00:00')
    ctime = stime
    times, shifts, vals = [],[], []
    while ctime < etime:
        print(ctime)
        #try:
        if True:
            rand = np.random.randint(0,60*60*24)
            cnt = 0
            while cnt <= 4:
                try:
                #if True:
                    st = client.get_waveforms(net, sta,'*', 'BHZ', ctime + rand, ctime + 60 + rand, attach_response=True)
                    st.remove_response()
                    break
                except:
                    cnt += 1
            fig = plt.figure(1, figsize=(10,14))
            plt.subplot(3,1,1)
            plt.plot(st[0].times(), st[0].data*10**6)
            plt.ylabel('Velocity ($\mu m/s$)')
            plt.text(-9, np.max(st[0].data*10**6), '(a)')
            plt.xlim(min(st[0].times()), max(st[0].times()))
            plt.xlabel('Time (s)')
            plt.subplot(3,1,2)
            plt.plot(st[1].times(), st[1].data*10**6)
            plt.xlim(min(st[1].times()), max(st[1].times()))
            plt.ylabel('Velocity ($\mu m/s$)')
            plt.xlabel('Time (s)')
            plt.text(-9, np.max(st[0].data*10**6), '(b)')
            st.filter('bandpass', freqmax=1/4., freqmin=1./8.)
            st.merge(fill_value=0)
            st.resample(1000)
            st.sort()
            print(st)
            tr1 = st.select(location='00')[0]
            tr2 = st.select(location='10')[0]
            cc = correlate(tr1.data, tr2.data, 500)
            plt.subplot(3,1,3)
            plt.plot((np.arange(len(cc))-500)/1000., cc)
            plt.xlim((min((np.arange(len(cc))-500)/1000.), max((np.arange(len(cc))-500)/1000.) ))
            plt.ylabel('Correlation')
            plt.xlabel('Lag (s)')
            plt.text(-.65, 1., '(c)')
            #plt.savefig('example.png', format='PNG')
            plt.savefig('example.pdf', format='PDF', dpi=400)
            import sys
            sys.exit()
            shift, val = xcorr_max(cc)
            shifts.append(shift)

            vals.append(val)
            times.append(str(ctime.year) + ', ' + str(ctime.julday) + ', ' + str(ctime.hour) + ', ' + str(ctime.minute) + str((ctime+rand).second))
        #except:
        #    pass
        ctime += 24*60*60

    with open(net + '_' + sta + '_' + str(year) + '.pickle2', 'wb') as f:
        pickle.dump([shifts, vals, times], f)
    return

inv = client.get_stations(network="IU", station="TUC",
                                starttime=UTCDateTime('2000-001T00:00:00'),
                                endtime=UTCDateTime('2020-001T00:00:00'))
stas = []
for net in inv:
    for sta in net:
        stas.append(net.code + '.' + sta.code)

stas = list(set(stas))

def calc_sta(sta):
    for year in range(2014,2020):
        if os.path.exists(sta.replace('.','_') + '_' + str(year) + '.pickle2'):
            if os.stat(sta.replace('.','_') + '_' + str(year) + '.pickle2').st_size < 500:
                pass
            else:
                continue
        
        calc_year_sta(year, sta)
        print('Finished: ' + sta + ' ' + str(year))
    return



calc_sta('IU.TUC')

#from multiprocessing import Pool
#pool = Pool(1)
#pool.map(calc_sta, stas)






