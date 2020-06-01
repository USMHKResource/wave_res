"""
Monthly Averaged Data
"""

from __future__ import division,print_function

import os,glob
from os.path import abspath as asbp
import numpy as np
import pylab as pl
import netCDF4, netcdftime
import itertools
import matplotlib as mpl
import scipy as spi
import scipy.interpolate
import datetime
import h5py

import cmocean

import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter,LatitudeFormatter
from cartopy.io import shapereader
import cartopy as crpy

import pynmd.physics.waves as gwaves
import pynmd.data.angles as gangles
import pynmd.data.signal as gsignal
import pynmd.tools.gtime as gtime
import pynmd.plotting as gplot

# ==============================================================================
# Load WW3 Output
# ==============================================================================

region = 'wc'

# Load Wind --------------------------------------------------------------------
wndFld = ('/pic/projects/fvwetland/gabriel/waveEnergyResource/' + 
          'assessment/hindcast/resource/' +
          'python/wave_res/ww3Averages/')
if region == 'wc':
    wndFile = wndFld + 'ww3.bulk.wc_10m.monthly_avg_wnd.nc'

nc = netCDF4.Dataset(wndFile,'r')
utime = netcdftime.utime(nc.variables['time'].units)
otWnd = utime.num2date(nc.variables['time'][:])
uwnd = nc.variables['uwnd'][:].data
vwnd = nc.variables['vwnd'][:].data
lon = nc.variables['longitude'][:].data
lat = nc.variables['latitude'][:].data
nc.close()
mes = np.array([aa.month for aa in otWnd])
[lon,lat] = np.meshgrid(lon,lat)

# Load local and potential resource --------------------------------------------
srcFld = ('/pic/projects/fvwetland/gabriel/waveEnergyResource/assessment/' + 
          'hindcast/resource/python/wave_res/results/')

# Load natural resource (freshly computed)
ff = h5py.File(srcFld + '/baseline/' + region + '.local-area.h5','r')
locNat = {'mes':np.array([np.int(aa[-2:]) for aa in ff['time'][:]]),
          'verts':ff['verts'][:],
          'resource':ff['sin'][:] + ff['sds'][:] + ff['snl'][:],
          'xy':ff['xy'][:]}
ff.close()

# Load potential resource
ff = h5py.File(srcFld + '/extraction/' + region + '.local-area.h5','r')
locPot = {'mes':np.array([np.int(aa[-2:]) for aa in ff['time'][:]]),
          'verts':ff['verts'][:],
          'resource':ff['sin'][:] + ff['sds'][:] + ff['snl'][:],
          'xy':ff['xy'][:]}
ff.close()


# Load the EEZ poligons --------------------------------------------------------
eezFld = ('/pic/projects/fvwetland/gabriel/waveEnergyResource/' + 
          'assessment/hindcast/resource/python/wave_res/py/resourceMaps/')
eezPoli = np.genfromtxt(eezFld + '/' + region + '_eez.points',delimiter=',')

# Close the polygon
if region == 'wc':
    # Close the polygon
    eezPoli = np.r_[eezPoli,np.array([[260,eezPoli[-1,1]]]),
                    np.array([[260,eezPoli[0,1]]])]
    # wcPoliExt[:,0] = gangles.wrapto360(wcPoliExt[:,0])

# Mask inside the Polygon
tmpxy = np.array([lon.flatten(),lat.flatten()]).T
flag = gplot.points_inside.inside_poly(tmpxy,eezPoli)   # True if inside
flag = np.logical_not(flag)                             # Flip this
flag = np.reshape(flag,lon.shape)                       # Reshape 

# Final mask
uwnd[:,flag] = np.NAN
vwnd[:,flag] = np.NAN

# ==============================================================================
# Make figures
# ==============================================================================
fs = 8
pl.rcParams['font.family'] = 'sans-serif'
pl.rcParams['font.sans-serif'] = 'Helvetica'
outFld = ('/mnt/fvwetland/gabriel/waveEnergyResource/' + 
          'assessment/energyExtraction/basin/94-figures/paper/')
prop_cycle = pl.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

# West Coast Domain ------------------------------------------------------------

# Create projection
plotProj = ccrs.PlateCarree(central_longitude=180)
dataProj = ccrs.PlateCarree()

# States --------------------------------------------------
borderDataset = ('/pic/projects/birthright/garc525/' + 
                'data/shoreline/gshhg/WDBII_shp/l/')
shp = shapereader.Reader(asbp(borderDataset + '/WDBII_border_l_L2.shp'))
states = crpy.feature.ShapelyFeature(shp.geometries(),
                                        dataProj,
                                        facecolor='none',
                                        edgecolor='darkgray')
# Main Coastline --------------------------------------------
coastDataset = ('/pic/projects/birthright/garc525/' + 
                'data/shoreline/gshhg/GSHHS_shp/l/')
shp = shapereader.Reader(asbp(coastDataset + '/GSHHS_l_L1.shp'))
coast = crpy.feature.ShapelyFeature(shp.geometries(),
                                    dataProj,
                                    facecolor='lightgray',
                                    edgecolor='darkgray',
                                    alpha=0.5)
# Countries --------------------------------------------------
borderDataset = ('/pic/projects/birthright/garc525/' + 
                'data/shoreline/gshhg/WDBII_shp/l/')
shp = shapereader.Reader(asbp(borderDataset + '/WDBII_border_l_L1.shp'))
countries = crpy.feature.ShapelyFeature(shp.geometries(),
                                        dataProj,
                                        facecolor='none',
                                        edgecolor='darkgray')

# Make the figure 
hf1 = pl.figure(figsize=(6.5,7.0))
season = {'Winter':[12,1,2],'Spring':[3,4,5],
          'Summer':[6,7,8],'Autumn':[9,10,11],
          'Yearly':np.arange(1,13)}
if region == 'wc':
    extents = [-130,-117,30.0,50.0]

for ii,aa in enumerate(['Winter','Spring','Summer','Autumn','Yearly']):
    
    # Wind Speed -------------------------------------------------
    if aa != 'Yearly':
        aveInd = np.logical_or(mes == season[aa][0],
                               np.logical_or(mes == season[aa][1],
                                             mes == season[aa][2]))
    else:
        aveInd = np.ones_like(mes,dtype=np.bool)

    ax3 = pl.subplot2grid((3,5),(0,ii),projection=plotProj)
    #ax1 = pl.subplot(1,1,1,projection=plotProj)
    ax3.set_extent(extents,crs=dataProj)

    ax3.add_feature(coast)
    #ax3.add_feature(states)
    #ax3.add_feature(countries)    

    wndMag = np.mean((uwnd[aveInd,...]**2 + vwnd[aveInd,...]**2)**0.5,axis=0)
    hcp3 = ax3.pcolormesh(lon,lat,wndMag,
                          cmap=pl.cm.gist_ncar,transform=dataProj)
    hcp3.set_clim([0,10])   
    q3 = ax3.quiver(lon[::10,::10],lat[::10,::10],
                    np.mean(uwnd[aveInd,::10,::10],axis=0)/wndMag[::10,::10],
                    np.mean(vwnd[aveInd,::10,::10],axis=0)/wndMag[::10,::10],
                    units='x',transform=dataProj,scale=2.0)
    
    ax3.set_title(aa,fontsize=fs)

    if ii == 0:
        ax3.text(-0.07, 0.55, 'Wind', va='bottom', ha='center',
                 rotation='vertical', rotation_mode='anchor',
                 transform=ax3.transAxes,fontsize=fs)        

    # Local resource --------------------------------------------
    ax0 = pl.subplot2grid((3,5),(1,ii),projection=plotProj)
    ax0.set_extent(extents,crs=dataProj)

    #ax0.add_feature(coast)
    #ax0.add_feature(states)
    #ax0.add_feature(countries)

    hcp0 = ax0.tripcolor(locNat['lon'],locNat['lat'],locNat['verts'],
                         np.mean(locNat['resource'][aveInd,...],axis=0),
                         cmap=pl.cm.gist_ncar,transform=dataProj,
                         shading='gouraud')
    #pl.colorbar(hcp0,ax=ax0)
    hcp0.set_clim([0,0.2])
    if ii == 0:
        ax0.text(-0.07, 0.55, 'Local Resource', va='bottom', ha='center',
                 rotation='vertical', rotation_mode='anchor',
                 transform=ax0.transAxes,fontsize=fs)        
        
        #ax0.text(-122,46.5,'WA',fontsize=fs,transform=dataProj)
        #ax0.text(-122,43.5,'OR',fontsize=fs,transform=dataProj)
        #ax0.text(-121.5,37,'CA',fontsize=fs,transform=dataProj)

    # Potential resource --------------------------------------------
    ax1 = pl.subplot2grid((4,6),(1,ii),projection=plotProj)
    ax1.set_extent(extents,crs=dataProj)

    ax1.add_feature(coast)
    ax1.add_feature(states)
    ax1.add_feature(countries)

    hcp1 = ax1.tripcolor(locPot['lon'],locPot['lat'],locPot['verts'],
                         np.mean(locPot['resource'][aveInd,...],axis=0),
                         cmap=pl.cm.gist_ncar,transform=dataProj,
                         shading='gouraud')
    #pl.colorbar(hcp1,ax=ax1)
    hcp1.set_clim([0,0.2])
    if ii == 0:
        ax1.text(-0.07, 0.55, 'Potential Resource', va='bottom', ha='center',
                 rotation='vertical', rotation_mode='anchor',
                 transform=ax1.transAxes,fontsize=fs)        



# Finishing touches
hf1.subplots_adjust(bottom=0.01,top=0.99,left=0.03,right=0.91,
                    wspace=0.00,hspace=0.05)
pl.pause(2)

# Manually Add Colorbars
ax0Pos = ax0.get_position()
cax0 = hf1.add_axes([0.92,ax0Pos.y0,0.015,ax0Pos.height])
hcb0 = pl.colorbar(hcp0,cax=cax0)
hcb0.ax.tick_params(labelsize=fs)
hcb0.ax.set_title(r'[W/m$^{2}$]',fontsize=fs)

ax1Pos = ax1.get_position()
cax1 = hf1.add_axes([0.92,ax1Pos.y0,0.015,ax1Pos.height])
hcb1 = pl.colorbar(hcp1,cax=cax1)
hcb1.ax.tick_params(labelsize=fs)
hcb1.ax.set_title(r'[W/m$^{2}$]',fontsize=fs)

ax2Pos = ax2.get_position()
cax2 = hf1.add_axes([0.92,ax2Pos.y0,0.02,ax2Pos.height])
hcb2 = pl.colorbar(hcp2,cax=cax2)
hcb2.ax.tick_params(labelsize=fs)
hcb2.ax.set_title('[kW/m]',fontsize=fs)

ax3Pos = ax3.get_position()
cax3 = hf1.add_axes([0.92,ax3Pos.y0,0.02,ax3Pos.height])
hcb3 = pl.colorbar(hcp3,cax=cax3)
hcb3.ax.tick_params(labelsize=fs)
hcb3.ax.set_title('[m/s]',fontsize=fs)

# Save the figure
hf1.savefig(outFld + '/westCoast_owp_monthly.png',dpi=500)

# ==============================================================================
# Compute average wind speeds over the domain
# ==============================================================================

wmag = []
for aa in range(len(uwnd)):
    wmag.append((uwnd[aa]**2 + vwnd[aa]**2)**0.5)
wmagMean = np.array([np.nanmean(aa) for aa in wmag])