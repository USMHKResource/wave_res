import matplotlib.pyplot as plt
import proj
import cartopy.feature as cfeature
import base
plt.ion()
reload(proj)

scale = '50m'
land = cfeature.NaturalEarthFeature('physical', 'land', scale,
                                    edgecolor='face', facecolor=cfeature.COLORS['land'])

region = 'at'
rinf = base.RegionInfo(region)
prj = proj.proj[region]
fignum = 10
fig = plt.figure(fignum)
fig.clf()
ax = plt.axes(projection=prj)
ax.add_feature(land)
ax.add_feature(cfeature.STATES.with_scale(scale))
ax.set_extent((prj.lonlim + prj.latlim), crs=proj.pc)
for cid in rinf.con_defs:
    ll = rinf.get_lonlat(cid)
    color = '0.5'
    if cid == 'EEZ':
        color = 'r'
    ax.plot(ll[0], ll[1], transform=proj.pc, color=color)

#dnow = 
