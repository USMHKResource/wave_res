from cartopy import crs
import numpy as np
import gis
import matplotlib.pyplot as plt

pc = crs.PlateCarree()
sphere = crs.Globe(ellipse='sphere')

# This is the radius used by the 'sphere' Proj4 projection (cartopy
# under the hood), as seen here:
# https://github.com/OSGeo/proj.4/blob/master/src/pj_ellps.c
Re = 6370997.0

proj0 = crs.AlbersEqualArea(-149.76, 61,
                            standard_parallels=[40, 63],
                            globe=sphere)
proj = crs.AlbersEqualArea(-149.76, 61,
                            standard_parallels=[40, 63])

triang = np.array([[-150.88, 58.6],
                   [-148.45, 57.73],
                   [-151.69, 56.09]]).T

#triang[0] -= 0
#triang[1] += -15.


def calc_excess(triang):
    _tr = np.deg2rad(triang[:, [0, 1, 2, 0]])
    b0 = np.exp(1j * gis.bearing(_tr, 'radians'))
    b1 = np.exp(1j * gis.bearing(_tr[:, ::-1], 'radians'))[[0, 2, 1]]
    return np.abs((np.angle(b0 / b1)).sum()) - np.pi


def calc_tri_area(triang):
    x = triang[0]
    y = triang[1]
    return np.abs(x[0] * (y[1] - y[2]) +
                  x[1] * (y[2] - y[0]) +
                  x[2] * (y[0] - y[1])) / 2


_tri = triang[:, [0, 1, 2, 0]]
fignum = 10
fig = plt.figure(10)
fig.clf()
ax = plt.axes(projection=pc)
ax.coastlines()
ax.plot(_tri[0], _tri[1], 'r-')
ax.plot(_tri[0, 0], _tri[1, 0], 'b.')
ax.plot(_tri[0, :2], _tri[1, :2], 'b-')
params = proj0.proj4_params
ax.plot(params['lon_0'], params['lat_0'], 'r+', ms=20, mew=4)
ax.axhline(params['lat_1'], linestyle='--', color='k')
ax.axhline(params['lat_2'], linestyle='--', color='k')
ax.set_ylim([50, 80])
ax.set_xlim([-180, -130])

A = calc_excess(triang) * Re ** 2

xy0 = proj0.transform_points(pc, triang[0], triang[1]).T[:2]
Ap0 = calc_tri_area(xy0)

xy = proj.transform_points(pc, triang[0], triang[1]).T[:2]
Ap = calc_tri_area(xy)

print("\n"
      "The spherical area is:           {: 8.3f} km^2\n"
      "The projected spherical area is: {: 8.3f} km^2 ({:.3f} %)\n"
      "The projected WGS84 area is:     {: 8.3f} km^2 ({:.3f} %)\n"
      .format(A / 1e6,
              Ap0 / 1e6, ((A - Ap0) / A) * 100,
              Ap / 1e6, ((A - Ap) / A) * 100,
      )
)
