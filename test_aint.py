from wave_res import area_integral as aint
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Delaunay
plt.ion()

xy = np.random.uniform(size=100000).reshape(2, -1)

tri = Delaunay(xy.T)

fig = plt.figure(2);fig.clf()
ax = fig.subplots(1, 1)
ax.plot(xy[0], xy[1], 'k.')

z = xy[0] + xy[1]

zsum, area = aint.sumtriangles(xy.T, z, tri.simplices)

print(zsum)
