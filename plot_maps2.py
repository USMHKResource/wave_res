import wave_res as wr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
import matplotlib.pyplot as plt
plt.ion()

"""
This script is used to show the clockwise vs. counterclockwise sense of the contours. Change which region you show by changing the next line.

########### NOTE ###########
The most important results here: Hawaii is CCW, while all other regions are CW.

This means that we need to switch the 'offshore' and 'onshore' flux calculations.
############################
"""

rinf = wr.RegionInfo('ak')

inds = rinf.con_defs['EEZ'][0]

plt.figure(1)
plt.clf()
plt.plot(rinf.gridlonlat[0][inds], rinf.gridlonlat[1][inds])
plt.scatter(rinf.gridlonlat[0][inds], rinf.gridlonlat[1][inds], c=np.arange(len(inds)))
plt.colorbar()
