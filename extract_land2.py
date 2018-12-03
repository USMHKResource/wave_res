import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
plt.ion()

dat = Dataset('data/ak_grid.nc')

lat = dat.variables['latitude'][:].data
lon = dat.variables['longitude'][:].data
depth = dat.variables['dpt'][0].data
bad = dat.variables['dpt'][0].mask
depth[bad] = np.NaN

fig = plt.figure(11)
fig.clf()
ax = fig.add_axes([.1, .1, .8, .8])
pcol = ax.pcolormesh(lon, lat, depth, vmin=0, vmax=500, cmap=plt.get_cmap('viridis_r'))
plt.colorbar(pcol)

