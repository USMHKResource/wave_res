import paths as p
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation

d0 = Dataset(p.tmpdir / 'ww3.ak.20090101.1d.nc', 'r')

lon, lat = d0.variables['lon'][:], d0.variables['lat'][:]
lon[lon > 0] -= 360

tri = Triangulation(lon, lat)

# tags = d0.variables['tag'][:]
# con = []
# last_start = None
# c = 0
# for tag in tags:
#     start = tag[:5]
#     if start != last_start:
#         con.append(start)
#         print('{} has {} points.'.format(last_start, c))
#         c = 0
#     c += 1
#     last_start = start

fignum = 10
fig = plt.figure(10)
fig.clf()
fig, ax = plt.subplots(1, 1, num=fignum)

# Plot the triangulation.
ax.set_aspect('equal')
ax.plot(lon, lat, 'b.', ms=1)
ax.set_title('triplot of Delaunay triangulation')
