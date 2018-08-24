import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
import scipy as sp
from py_sphere_Voronoi import voronoi_utility
plt.ion()

#pin down the pseudo random number generator (prng) object to avoid certain pathological generator sets
prng = np.random.RandomState(117) #otherwise, would need to filter the random data to ensure Voronoi diagram is possible
#produce 1000 random points on the unit sphere using the above seed
random_coordinate_array = voronoi_utility.generate_random_array_spherical_generators(1000,1.0,prng)
#produce the Voronoi diagram data
voronoi_instance = voronoi_utility.Voronoi_Sphere_Surface(random_coordinate_array,1.0)
dictionary_voronoi_polygon_vertices = voronoi_instance.voronoi_region_vertices_spherical_surface()
#plot the Voronoi diagram
fig = plt.figure()
fig.set_size_inches(2,2)
ax = fig.add_subplot(111, projection='3d')
for generator_index, voronoi_region in dictionary_voronoi_polygon_vertices.iteritems():
   random_color = colors.rgb2hex(sp.rand(3))
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
