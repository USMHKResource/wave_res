import paths as p
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from py_sphere_Voronoi import voronoi_utility as sv
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

d0 = Dataset(p.srcdir / '2009/src/ww3.ak.20090101.1d.nc', 'r')

lon, lat = d0.variables['lon'][:].data, d0.variables['lat'][:].data
lon[lon > 0] -= 360

tri = Triangulation(lon, lat)

# earth radius in meters
re = 6.371e6
re = 1.0

theta = np.deg2rad((90 - lat))
rtp = np.stack([re * np.ones_like(lon), np.deg2rad(lon), theta]).T

xyz = sv.convert_spherical_array_to_cartesian_array(rtp)

vor = sv.Voronoi_Sphere_Surface(xyz, re)
d_verts = vor.voronoi_region_vertices_spherical_surface()

#plot the Voronoi diagram
fig = plt.figure()
fig.set_size_inches(2,2)
ax = fig.add_subplot(111, projection='3d')
for generator_index, voronoi_region in d_verts.iteritems():
   random_color = colors.rgb2hex(np.random.rand(3))
   #fill in the Voronoi region (polygon) that contains the generator:
   polygon = Poly3DCollection([voronoi_region],alpha=1.0)
   polygon.set_color(random_color)
   ax.add_collection3d(polygon)
ax.set_xlim(-1,1);ax.set_ylim(-1,1);ax.set_zlim(-1,1);
(-1, 1)
(-1, 1)
(-1, 1)
ax.set_xticks([-1,1]);ax.set_yticks([-1,1]);ax.set_zticks([-1,1]); 
plt.tick_params(axis='both', which='major', labelsize=6)


a = vor.voronoi_region_surface_areas_spherical_surface()

error
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
