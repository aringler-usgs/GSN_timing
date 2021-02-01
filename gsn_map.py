#!/usr/bin/env python
from obspy.core import UTCDateTime
from obspy.clients.fdsn import Client
import matplotlib.pyplot as plt
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
mpl.rc('font', size=18)
client = Client('IRIS')
# Initial networks



fig = plt.figure(1, figsize=(10,8))

stime = UTCDateTime('2020-001T00:00:00')
etime = stime + 6*24*60*60
networks = ['II', 'IU', 'CU', 'IC']
for idx, netname in enumerate(networks):
    inv = client.get_stations(starttime=stime, endtime=etime, network=netname,
            station = "*")
    lats, lons = [], []
    for net in inv:
        for sta in net:
            lats.append(sta.latitude)
            lons.append(sta.longitude)
    if idx == 0:
        ax = fig.add_subplot(1,1,1, projection = ccrs.Robinson())           
        ax.coastlines()
        ax.set_title('Global Seismographic Network 2020')
        ax.add_feature(cart.feature.LAND, zorder=2, edgecolor='k')
        ax.set_global()
        ax.coastlines()
    if idx > 1 :
        im = ax.plot(lons, lats, '.', color='C' + str(idx+3), transform=ccrs.Geodetic(), zorder=4, markersize=14, alpha=.5, label=netname )
    else:
        im = ax.plot(lons, lats, '.', color='C' + str(idx+1), transform=ccrs.Geodetic(), zorder=4, markersize=14, alpha=.5, label=netname )
plt.legend(ncol=4, bbox_to_anchor=(0.5, -0.2), loc='lower center')






plt.savefig('GSNmap.pdf', format='PDF', dpi=400)
plt.savefig('GSNmap.png', format='PNG', dpi=400)