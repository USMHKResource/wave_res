import matplotlib.pyplot as plt
plt.ion()
from base import RegionInfo
import numpy as np
import cartopy.feature as cfeature


scale = '50m'
land = cfeature.NaturalEarthFeature('physical', 'land', scale,
                                    edgecolor='face',
                                    facecolor=cfeature.COLORS['land'])


def setup_figure(rinf, fignum, **kwargs):
    # First setup the figure
    prj = rinf.proj
    fig = plt.figure(fignum)
    fig.clf()
    ax = plt.axes(projection=prj)
    ax.add_feature(land)
    ax.add_feature(cfeature.STATES.with_scale(scale))
    ax.set_extent((prj.lonlim + prj.latlim), crs=rinf._proj_pc)
    return fig, ax


def mindist(pt, line):
    d2 = ((pt[:, None] - line) ** 2).sum(0)
    idx = np.argmin(d2)
    return idx, np.sqrt(d2[idx])


def polyinds_wc():

    out = {}

    rinf = RegionInfo('wc')
    fig, ax = setup_figure(rinf, 1000 + 30)

    for rng in np.arange(10, 200, 10):
        ky = '{:03d}'.format(rng)

        xy = rinf.get_contour(ky, xy=True)[0]

        poly_inds = []

        border_ind = 0
        bxy = rinf.get_contour('borders', xy=True)[border_ind]
        binds = rinf.con_defs['borders'][border_ind]
        id0 = mindist(xy[:, 0], bxy)[0]
        poly_inds.append(binds[:id0 + 1])
        poly_inds.append(rinf.con_defs[ky][0])

        border_ind = 1
        bxy = rinf.get_contour('borders', xy=True)[border_ind]
        binds = rinf.con_defs['borders'][border_ind]
        ide = mindist(xy[:, -1], bxy)[0]
        poly_inds.append(binds[ide:])

        pinds = np.hstack(poly_inds).tolist()
        boundary = np.hstack([
            rinf.gridxy[:, pinds],
            rinf.transform(rinf.mainland)])
        ax.plot(boundary[0], boundary[1], 'r-')

        out[ky] = pinds
    return out


def polyinds_ak():
    out = {}
    rinf = RegionInfo('ak')
    fig, ax = setup_figure(rinf, 1000 + 31)

    for rng in np.arange(10, 200, 10):
        ky = '{:03d}'.format(rng)

        seg_i = 0
        xy = rinf.get_contour(ky, xy=True)[seg_i]

        poly_inds = []

        # link_seg_with_border(rinf, ky, 0, 0, 2)
        border_ind = 0
        bxy = rinf.get_contour('borders', xy=True)[border_ind]
        binds = rinf.con_defs['borders'][border_ind]
        id0, d = mindist(xy[:, 0], bxy)
        poly_inds.append(binds[:id0 + 1])
        poly_inds.append(rinf.con_defs[ky][seg_i])

        if rng > 81:
            border_ind = 1
            binds = rinf.con_defs['borders'][border_ind]
            bxy = rinf.get_contour('borders', xy=True)[border_ind]

            ide, d = mindist(xy[:, -1], bxy)

            seg_i += 1
            xy = rinf.get_contour(ky, xy=True)[seg_i]
            id0, d = mindist(xy[:, 0], bxy)
            poly_inds.append(binds[ide:id0 + 1])
            poly_inds.append(rinf.con_defs[ky][seg_i])

        border_ind = 2
        bxy = rinf.get_contour('borders', xy=True)[border_ind]
        binds = rinf.con_defs['borders'][border_ind]
        ide, d = mindist(xy[:, -1], bxy)
        poly_inds.append(binds[ide:])

        pinds = np.hstack(poly_inds).tolist()
        boundary = np.hstack([
            rinf.gridxy[:, pinds],
            rinf.transform(rinf.mainland)])
        ax.plot(boundary[0], boundary[1], 'r-')

        other = rinf.con_defs[ky][seg_i + 1:]
        for itmp, slc in enumerate(other):
            slc = other[itmp] = slc + [slc[0]]
            ax.plot(rinf.gridxy[0, slc], rinf.gridxy[1, slc], 'b-')
        out[ky] = [pinds, ] + other
    out['EEZ'] = out['eez'] = [rinf.con_defs['EEZ'][0], ]
    boundary = rinf.gridxy[:, out['EEZ'][0]]
    ax.plot(boundary[0], boundary[1], '-', color='y', alpha=0.5, lw=7)

    return out


def polyinds_at():
    out = {}
    rinf = RegionInfo('at')
    fig, ax = setup_figure(rinf, 1000 + 39)
    rinf.con_defs

    for inds in rinf.con_defs['borders']:
        x, y = rinf.gridxy[:, inds]
        ax.plot(x, y, 'm-', lw=3, zorder=10)
    inds = rinf.con_defs['EEZ'][0]
    x, y = rinf.gridxy[:, inds]
    ax.plot(x, y, 'k-', lw=1, zorder=12)

    for rng in np.arange(10, 200, 10):
        ky = '{:03d}'.format(rng)

        seg_i = 0
        xy = rinf.get_contour(ky, xy=True)[seg_i]

        poly_inds = []

        # link_seg_with_border(rinf, ky, 0, 0, 2)
        border_ind = 0
        bxy = rinf.get_contour('borders', xy=True)[border_ind]
        binds = rinf.con_defs['borders'][border_ind]
        id0, d = mindist(xy[:, 0], bxy)
        poly_inds.append(binds[:id0 + 1])
        poly_inds.append(rinf.con_defs[ky][seg_i])

        if rng > 21:
            border_ind = 1
            binds = rinf.con_defs['borders'][border_ind]
            bxy = rinf.get_contour('borders', xy=True)[border_ind]

            ide, d = mindist(xy[:, -1], bxy)

            seg_i += 1
            xy = rinf.get_contour(ky, xy=True)[seg_i]
            id0, d = mindist(xy[:, 0], bxy)
            poly_inds.append(binds[ide:id0 + 1])
            poly_inds.append(rinf.con_defs[ky][seg_i])

        if rng == 190:
            poly_inds[2] = range(133, 259)

        border_ind = 2
        bxy = rinf.get_contour('borders', xy=True)[border_ind]
        binds = rinf.con_defs['borders'][border_ind]
        ide, d = mindist(xy[:, -1], bxy)
        poly_inds.append(binds[ide:])

        pinds = np.hstack(poly_inds).tolist()
        boundary = np.hstack([
            rinf.gridxy[:, pinds],
            rinf.transform(rinf.mainland)])
        ax.plot(boundary[0], boundary[1], 'r-')

        other = rinf.con_defs[ky][seg_i + 1:]
        for itmp, slc in enumerate(other):
            slc = other[itmp] = slc + [slc[0]]
            ax.plot(rinf.gridxy[0, slc], rinf.gridxy[1, slc], 'b-')
        out[ky] = [pinds, ] + other
    out['EEZ'] = out['eez'] = [rinf.con_defs['EEZ'][0], ]
    boundary = rinf.gridxy[:, out['EEZ'][0]]
    ax.plot(boundary[0], boundary[1], '-', color='y', alpha=0.5, lw=7)

    return out


def polyinds_prusvi():
    out = {}
    rinf = RegionInfo('prusvi')
    fig, ax = setup_figure(rinf, 1000 + 32)

    cmap = plt.get_cmap('tab20')

    for rng in range(10, 200, 10):
        ky = '{:03d}'.format(rng)
        seg_i = 0
        xy = rinf.get_contour(ky, xy=True)[seg_i]
        ax.plot(xy[0, 0], xy[1, 0], 'r+')

        poly_inds = []
        # link_seg_with_border(rinf, ky, 0, 0, 2)
        border_ind = 0
        bxy = rinf.get_contour('borders', xy=True)[border_ind]
        ax.plot(bxy[0, 0], bxy[1, 0], 'b+')
        binds = rinf.con_defs['borders'][border_ind]
        id0, d = mindist(xy[:, -1], bxy)
        poly_inds.append(rinf.con_defs[ky][seg_i])

        if rng < 171:
            seg_i += 1
            xy = rinf.get_contour(ky, xy=True)[seg_i]
            ax.plot(xy[0, 0], xy[1, 0], 'rx')
            ide, d = mindist(xy[:, 0], bxy)
            poly_inds.append(binds[id0:ide])
            poly_inds.append(rinf.con_defs[ky][seg_i])

        if rng != 10:
            id0, d = mindist(xy[:, -1], bxy)
            seg_i = 0
            xy = rinf.get_contour(ky, xy=True)[seg_i]
            ide, d = mindist(xy[:, 0], bxy)
            poly_inds.append(binds[id0:ide])
            seg_i = 1

        poly_inds.append(rinf.con_defs[ky][0][0])

        pinds = np.hstack(poly_inds).tolist()
        boundary = rinf.gridxy[:, pinds]
        ax.plot(boundary[0], boundary[1], '-', color=cmap((rng - 5) / 200.),
                alpha=0.5)

        other = rinf.con_defs[ky][seg_i + 1:]
        for itmp, slc in enumerate(other):
            slc = other[itmp] = slc + [slc[0]]
            ax.plot(rinf.gridxy[0, slc], rinf.gridxy[1, slc], 'b-')
        out[ky] = [pinds, ] + other
    rng = 200

    out['EEZ'] = out['eez'] = [rinf.con_defs['EEZ'][0] +
                               [rinf.con_defs['EEZ'][0][0], ], ]
    boundary = rinf.gridxy[:, out['EEZ'][0]]
    ax.plot(boundary[0], boundary[1], '-',
            color=cmap((rng - 5) / 200.), alpha=0.5)

    return out


#def polyinds_hi():
if False:
    out = {}
    rinf = RegionInfo('hi')
    fig, ax = setup_figure(rinf, 1000 + 32)

    cmap = plt.get_cmap('tab20')

    for rng in range(10, 200, 10):
        ky = '{:03d}'.format(rng)
        seg_i = 0
        xy = rinf.get_contour(ky, xy=True)[seg_i]
        ax.plot(xy[0, 0], xy[1, 0], 'r+')

        poly_inds = []
        # link_seg_with_border(rinf, ky, 0, 0, 2)
        border_ind = 0
        bxy = rinf.get_contour('borders', xy=True)[border_ind]
        ax.plot(bxy[0, 0], bxy[1, 0], 'b+')
        binds = rinf.con_defs['borders'][border_ind]
        id0, d = mindist(xy[:, -1], bxy)
        poly_inds.append(rinf.con_defs[ky][seg_i])

        if rng < 51:
            pass
        elif rng == 60:
            seg_i += 1
            poly_inds[0] = poly_inds[0][:-7]
            poly_inds.append(rinf.con_defs[ky][seg_i][13:])
        elif rng == 70:
            seg_i += 1
            poly_inds[0] = poly_inds[0][:-7]
            poly_inds.append(rinf.con_defs[ky][seg_i][15:])
        elif rng == 80:
            poly_inds[0] = poly_inds[0][15:-5]
        elif rng == 90:
            poly_inds[0] = poly_inds[0][17:-3]

        poly_inds.append(poly_inds[0])

        pinds = np.hstack(poly_inds).tolist()
        boundary = rinf.gridxy[:, pinds]
        ax.plot(boundary[0], boundary[1], '-', color=cmap((rng - 5) / 200.),
                alpha=0.5)

        other = rinf.con_defs[ky][seg_i + 1:]
        for itmp, slc in enumerate(other):
            slc = other[itmp] = slc + [slc[0]]
            ax.plot(rinf.gridxy[0, slc], rinf.gridxy[1, slc], 'b-')
        out[ky] = [pinds, ] + other
    rng = 200
    out['EEZ'] = out['eez'] = [rinf.con_defs['EEZ'][0] +
                               [rinf.con_defs['EEZ'][0][0], ], ]
    boundary = rinf.gridxy[:, out['EEZ'][0]]
    ax.plot(boundary[0], boundary[1], '-',
            color=cmap((rng - 5) / 200.), alpha=0.5)

    #return out


def run_all():

    poly_defs = {}

    poly_defs['wc'] = polyinds_wc()
    poly_defs['ak'] = polyinds_ak()
    poly_defs['prusvi'] = polyinds_prusvi()
    poly_defs['at'] = polyinds_at()
    ##### IMPORTANT ####
    # HI is not yet included in the poly defs b/c it's boundaries are wonky.

    return poly_defs
