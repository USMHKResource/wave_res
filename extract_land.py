import matplotlib.pyplot as plt
plt.ion()
import proj
import cartopy.feature as cfeature
import base
import numpy as np

scale = '50m'
land = cfeature.NaturalEarthFeature('physical', 'land', scale,
                                    edgecolor='face', facecolor=cfeature.COLORS['land'])


def overlap(r0, r1):
    return (within(r0[0], r1) or within(r0[1], r1) or
            within(r1[0], r0) or within(r1[1], r0))


def within(val, r):
    return (r[0] < val) & (val < r[1])

   
def plot_lines(rinf, ax, **kwargs):
    eez_kw = kwargs.pop('eez_kw', {})
    kwargs.update(transform=proj.pc)
    for cid in rinf.con_defs:
        tmp = kwargs.copy()
        if cid == 'eez':
            tmp.update(eez_kw)
        for ll in rinf.get_lonlat(cid):
            ax.plot(ll[0], ll[1], **tmp)

region = 'wc'

prj = proj.proj[region]
fignum = 31
fig = plt.figure(fignum)
fig.clf()
ax = plt.axes(projection=prj)
feat = ax.add_feature(land)
ax.add_feature(cfeature.STATES.with_scale(scale))
ax.set_extent((prj.lonlim + prj.latlim), crs=proj.pc)

rinf = base.RegionInfo(region, use_old_con_defs=False)
prj = proj.proj[region]
ax.set_extent((prj.lonlim + prj.latlim), crs=proj.pc)

plot_lines(rinf, ax, color='m', 
           eez_kw=dict(lw=3))

cons = []

for idx, d in enumerate(land.geometries()):
    b = d.bounds
    lonb = [b[0], b[2]]
    latb = [b[1], b[3]]
    if overlap(lonb, prj.lonlim) and overlap(latb, prj.latlim):
        #print("Hello!")
        print('{}: {}'.format(idx, d.area))
        for id2, p in enumerate(d):
            lldat = np.stack(p.boundary.xy)
            cons.append(lldat)

            ax.plot(lldat[0], lldat[1], '.', transform=proj.pc)

    
