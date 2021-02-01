#!/usr/bin/env python
import glob
import matplotlib.pyplot as plt
import pickle
import numpy as np
from scipy.stats import sem, t

import matplotlib as mpl
mpl.rc('font', family = 'serif')
mpl.rc('font', serif = 'Times')
mpl.rc('text', usetex = True)
mpl.rc('font', size = 20)


files = glob.glob('*.pickle')
shiftsall, valsall, timesall = [], [], []
for cfile in files:
    with open(cfile, 'rb') as f:
        shifts, vals, times = pickle.load(f)
        shiftsall += shifts
        valsall += vals
        timesall += times

valsall = np.array(valsall)
shiftsall = np.array(shiftsall)
timesall = np.array(timesall)
print(len(shiftsall))
shiftsall = shiftsall[(valsall >= 0.98)]
shiftsall[shiftsall >= 100] = 100.
shiftsall[shiftsall <= -100] = -100.
print(len(shiftsall))


fig = plt.figure(1, figsize=(12,12))
plt.subplot(2,1,1)
plt.text(-140,0.09, '(a)')
plt.hist(shiftsall, bins=200, range=(-100,100), density=True)
plt.xlabel('Timing Error (ms)')
plt.ylabel('Hits (Normalized)')
print(np.std(shiftsall))
print(np.std(shiftsall)*2.5)

files = glob.glob('*2019*.pickle')
shiftsall, valsall, timesall = [], [], []
for cfile in files:
    with open(cfile, 'rb') as f:
        shifts, vals, times = pickle.load(f)
        shiftsall += shifts
        valsall += vals
        timesall += times

valsall = np.array(valsall)
shiftsall = np.array(shiftsall)
timesall = np.array(timesall)
print(len(shiftsall))
shiftsall = shiftsall[(valsall >= 0.98)]
shiftsall[shiftsall >= 100] = 100.
shiftsall[shiftsall <= -100] = -100.
print(len(shiftsall))



plt.subplot(2,1,2)
plt.hist(shiftsall, bins=200, range=(-100,100), density=True)
plt.text(-140,0.135, '(b)')
plt.xlabel('Timing Error (ms)')
plt.ylabel('Hits (Normalized)')

plt.savefig('hist.png', format='PNG')


print(np.std(shiftsall))
print(np.std(shiftsall)*2.5)


#plt.savefig('histo2019.png', format='PNG')


# std_err = sem(shiftsall)
# h = std_err *t.ppf((1+0.99)/2., len(shiftsall)-1)

# print(np.mean(shiftsall) + h)
# print(np.mean(shiftsall) - h)
# print(np.mean(shiftsall))
# print(h)