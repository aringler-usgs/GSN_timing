#!/usr/bin/env python
import glob
import matplotlib.pyplot as plt
import pickle
import numpy as np

import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times')
mpl.rc('text', usetex=True)
mpl.rc('font',size=20)


stas = glob.glob('*SJG*.pickle')
stas = list(set([sta.split('_')[1] for sta in stas]))

for sta in stas:
    fig = plt.figure(1, figsize=(12,12))
    files = glob.glob('*' + sta + '*.pickle')
    shiftsall, valsall, timesall = [], [], []
    for cfile in files:
        print(cfile)
        
        with open(cfile, 'rb') as f:
            shifts, vals, times = pickle.load(f)
            shiftsall += shifts
            valsall += vals
            timesall += times
    valsall = np.array(valsall)
    shiftsall = np.array(shiftsall)
    timesall = np.array(timesall)

    shiftsall = shiftsall[(valsall >= 0.98)]
    timesall = timesall[(valsall >=0.98)]
    valsall = valsall[(valsall>=0.98)]
    shiftsall[shiftsall >= 100] = 100.
    shiftsall[shiftsall <= -100] = -100.
    goodtimes =[]
    for time in timesall:
        time = time.replace(' ','')
        time = time.split(',')
        goodtimes.append(float(time[0]) + float(time[1])/365.25)
    goodtimes = np.array(goodtimes)
    plt.subplot(2,1,1)
    plt.title(sta + ' 00 versus 10 BHZ')
    plt.plot(goodtimes, shiftsall , '.')
    plt.text(2009, 100, '(a)' )
    plt.ylabel('Lag (ms)')
    plt.ylim((-105,105))
    plt.xticks(np.arange(1998,2022,2))
    plt.xlim((min(goodtimes)-0.01, max(goodtimes)+0.1))
    plt.subplot(2,1,2)
    plt.text(2009, 1., '(b)')
    plt.plot(goodtimes, valsall, '.')
    plt.xticks(np.arange(1998,2022,2))
    plt.ylim((0.98, 1.))
    plt.xlim((min(goodtimes)-0.01, max(goodtimes)+0.1))
    plt.xlabel('Time (Year)')
    plt.ylabel('Correlation')
    plt.savefig(sta + '.png', format='PNG')
    plt.savefig(sta + '.pdf', format='PDF')
    plt.clf()
    plt.close()
    




