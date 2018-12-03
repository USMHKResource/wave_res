import matplotlib.pyplot as plt
plt.ion()
import proj
import cartopy.feature as cfeature
import base
reload(proj)
reload(base)

scale = '50m'
land = cfeature.NaturalEarthFeature('physical', 'land', scale,
                                    edgecolor='face', facecolor=cfeature.COLORS['land'])

   
def plot_lines(rinf, ax, **kwargs):
    eez_kw = kwargs.pop('eez_kw', {})
    kwargs.update(transform=proj.pc)
    for cid in rinf.con_defs:
        tmp = kwargs.copy()
        if cid == 'eez':
            tmp.update(eez_kw)
        for ll in rinf.get_lonlat(cid):
            ax.plot(ll[0], ll[1], **tmp)


def show_atlantic():
    region = 'at'
    rinf = base.RegionInfo(region)
    rinfec = base.RegionInfo('ec')
    rinfgm = base.RegionInfo('gm')
    prj = proj.proj[region]
    fignum = 10
    fig = plt.figure(fignum)
    fig.clf()
    ax = plt.axes(projection=prj)
    ax.add_feature(land)
    ax.add_feature(cfeature.STATES.with_scale(scale))
    ax.set_extent((prj.lonlim + prj.latlim), crs=proj.pc)

    plot_lines(rinfec, ax, color='m',
               eez_kw=dict(lw=3))
    plot_lines(rinfgm, ax, color='c',
               eez_kw=dict(lw=3))


def show_prusvi():

    region = 'prusvi'
    rinf = base.RegionInfo(region)
    prj = proj.proj[region]
    fignum = 11
    fig = plt.figure(fignum)
    fig.clf()
    ax = plt.axes(projection=prj)

    ax.set_extent((prj.lonlim + prj.latlim), crs=proj.pc)

    land2 = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                         edgecolor='face', facecolor=cfeature.COLORS['land'])
    land2._kwargs['facecolor'] = [0.5, 0.5, 0.5]
    #ax.add_feature(land2)
    # ax.add_feature(cfeature.STATES.with_scale(scale))

    plot_lines(rinf, ax, color='m',
               eez_kw=dict(lw=3))


def show_ak(old_con=False):
    region = 'ak'
    rinf = base.RegionInfo(region, use_old_con_defs=old_con)
    prj = proj.proj[region]
    fignum = 12
    fig = plt.figure(fignum)
    fig.clf()
    ax = plt.axes(projection=prj)
    ax.add_feature(land)
    ax.add_feature(cfeature.STATES.with_scale(scale))
    ax.set_extent((prj.lonlim + prj.latlim), crs=proj.pc)

    plot_lines(rinf, ax, color='m', 
               eez_kw=dict(lw=3))



def show_hi(old_con=False):
    region = 'hi'
    rinf = base.RegionInfo(region, use_old_con_defs=old_con)
    prj = proj.proj[region]
    fignum = 13
    fig = plt.figure(fignum)
    fig.clf()
    ax = plt.axes(projection=prj)
    # ax.add_feature(land)
    # ax.add_feature(cfeature.STATES.with_scale(scale))
    ax.set_extent((prj.lonlim + prj.latlim), crs=proj.pc)

    plot_lines(rinf, ax, color='m', 
               eez_kw=dict(lw=3))
