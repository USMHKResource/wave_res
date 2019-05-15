import wave_res as wr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

rinf = wr.RegionInfo('hi')

fig = plt.figure(10)
fig.clf()
ax = plt.axes(projection=rinf.proj)

inds_all = rinf.con_defs['eez']
for inds in inds_all:
    x, y = rinf.gridlonlat[:, inds]
    ax.scatter(x, y, color='b')

# inds_all = rinf.con_defs['borders']
# for inds in inds_all:
#     x, y = rinf.gridlonlat[:, inds]
#     ax.scatter(x, y, 5, color='r', )

inds_all = rinf.con_defs['200']
for inds in inds_all:
    x, y = rinf.gridlonlat[:, inds]
    ax.scatter(x, y, color='m')
