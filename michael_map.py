#!/usr/bin/env python
import glob
import sys
import matplotlib.pyplot as plt
import numpy as np
from obspy.core import UTCDateTime
from obspy.clients.fdsn import Client
import matplotlib as mpl
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import cartopy.feature as cfeature
import matplotlib.font_manager
import cartopy as cart

# Importing and applying font
mpl.rc('font', family = 'serif')
mpl.rc('font', serif = 'Times') 
mpl.rc('text', usetex = True)
mpl.rc('font', size=16)

client = Client()
stime = UTCDateTime('2000-001T00:00:00')
etime = UTCDateTime('2020-001T00:00:00')




def setupmap(handle):


    handle.add_feature(cfeature.LAND)
    handle.add_feature(cfeature.OCEAN)
    handle.add_feature(cfeature.COASTLINE)
    handle.add_feature(cfeature.BORDERS, linestyle=':')
    handle.add_feature(cfeature.LAKES, alpha=0.5)
    handle.add_feature(cfeature.RIVERS)
    handle.add_feature(cfeature.STATES, edgecolor='gray')
    return handle
    

fig = plt.figure(1, figsize=(12,10)) 
ax = fig.add_subplot(1,1,1, projection = ccrs.Robinson())
ax.coastlines()
ax.add_feature(cart.feature.LAND, zorder=2, edgecolor='k')
ax.set_global()
ax.coastlines()

for idx, nets in enumerate(['IU', 'IC', 'II', 'CU']):
    nd = True
    inv = client.get_stations(starttime=stime, endtime=etime, network=nets,
        channel="LHZ", level='response')
    for net in inv:
        for sta in net:
            print(sta)
            try:
            #if True:
                sncl = net.code + '.' + sta.code + '.' + sta[0].location_code + '.' + sta[0].code
                coors = inv.get_coordinates(sncl, etime)
                if nd: 
                    plt.plot(coors['longitude'], coors['latitude'], marker='^', color='C' + str(idx), transform=ccrs.Geodetic(), zorder=10, label=net.code)
                    nd= False
                else:
                    plt.plot(coors['longitude'], coors['latitude'], marker='^', color='C' + str(idx), transform=ccrs.Geodetic(), zorder=10)
            except:
                pass


plt.legend(ncol=4, bbox_to_anchor=(0.5, -0.2), loc='lower center')


#fig.legend(ncol=4, bbox_to_anchor=(0.95,0.13), )  
#cb_ax = fig.add_axes([0.2, 0.05, 0.6, 0.03])
#cbar = fig.colorbar(im, cax=cb_ax, orientation='horizontal')

#cbar.set_label('Day of Clock Quality Below 80%') 

plt.savefig('Map_of_clock.png', format='PNG', dpi=400)
plt.show()