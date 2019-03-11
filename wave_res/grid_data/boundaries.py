"""
This module calculates the 'clipping paths' (boundaries) for each
distance (range path) from the coast. These boundaries can be combined
with the land boundaries (extracted by the ``extract_land.py`` script)
to create complete clipping paths. These clipping paths are used to
clip triangles in the area integration.

Each function ``polyinds_<region>`` defined below returns a list of
lists. The first list (in the outer list) is always the 'main
boundary', and subsequent ones are additional (island) boundaries.
The island boundaries are always loops around islands that do not
intersect the 'main boundary'.

The 'main boundary' for each range path is created by finding the
points on the 'borders' paths (defined in ``base.py``) that are
closest to the endpoints of each range path. The indices for the
segment of the 'border' between land and that point is then joined
with the indices for that range path, and these are returned to define
the main boundary. If the boundary is a loop, the first index in this
list is the same as the last. If the boundary intersects land (the
mainland), then the first and last indices are different. The mainland
is intended to be joined to this data in later processing.

Each function is essentially hard-coded to join range paths with the
borders in the correct order. If you want to add a new grid, or if the
definitions of the underlying grids change you'll need to carefully
create the boundaries for that grid. The minimal automation involved
is in finding the point along the border that is nearest the endpoint
of each range path.

"""


import matplotlib.pyplot as plt
import numpy as np
import setpath
from base import RegionInfo
import proj

plt.ion()


def setup_figure(rinf, fignum, **kwargs):
    # First setup the figure
    prj = rinf.proj
    fig = plt.figure(fignum)
    fig.clf()
    ax = plt.axes(projection=prj)
    ax.add_feature(proj.land)
    ax.add_feature(proj.states)
    ax.set_extent((prj.lonlim + prj.latlim), crs=rinf._proj_pc)
    return fig, ax


def mindist(pt, line):
    d2 = ((pt[:, None] - line) ** 2).sum(0)
    idx = np.argmin(d2)
    return idx, np.sqrt(d2[idx])


def polyinds_wc():
    """The west coast boundary definitions are the simplest because
    there are no islands, and exactly two borders (at each end) that
    intersect all range paths.
    """
    out = {}

    rinf = RegionInfo('wc')
    fig, ax = setup_figure(rinf, 1000 + 30)

    # loop over range paths
    for rng in np.arange(10, 200, 10):
        ky = '{:03d}'.format(rng)

        xy = rinf.get_contour(ky, xy=True)[0]

        # Initialize the list of indices for this range.
        poly_inds = []

        border_ind = 0
        bxy = rinf.get_contour('borders', xy=True)[border_ind]
        binds = rinf.con_defs['borders'][border_ind]
        # Find the border (bxy) point nearest the start of the range
        # path (xy)
        id0 = mindist(xy[:, 0], bxy)[0]
        # Now add the indices of the indices of the border
        poly_inds.append(binds[:id0 + 1])
        # Then the range path indices
        poly_inds.append(rinf.con_defs[ky][0])

        border_ind = 1
        bxy = rinf.get_contour('borders', xy=True)[border_ind]
        binds = rinf.con_defs['borders'][border_ind]
        # Find this border (bxy) point nearest the end of the range
        # path (xy)
        ide = mindist(xy[:, -1], bxy)[0]
        # append the border indices
        poly_inds.append(binds[ide:])

        # Stack the all together.
        pinds = np.hstack(poly_inds).tolist()
        boundary = np.hstack([
            rinf.gridxy[:, pinds],
            rinf.transform(rinf.mainland)])
        ax.plot(boundary[0], boundary[1], 'r-')

        # The output is a dict (range) of lists of lists for each range.
        # But, the west coast is simple, so there are no island boundaries.
        out[ky] = [pinds, ]
    # The 'EEZ' range path IS the boundary for 'eez' and 'EEZ'.
    out['EEZ'] = out['eez'] = [rinf.con_defs['EEZ'][0], ]
    return out


def polyinds_ak():
    """The Alaska boundary definitions have many 'island' boundaries
    and three 'borders', so there is some hardcoding to treat
    different ranges differently.
    """
    out = {}
    rinf = RegionInfo('ak')
    fig, ax = setup_figure(rinf, 1000 + 31)

    for rng in np.arange(10, 200, 10):
        ky = '{:03d}'.format(rng)

        seg_i = 0
        xy = rinf.get_contour(ky, xy=True)[seg_i]

        poly_inds = []

        border_ind = 0
        bxy = rinf.get_contour('borders', xy=True)[border_ind]
        binds = rinf.con_defs['borders'][border_ind]
        id0, d = mindist(xy[:, 0], bxy)
        poly_inds.append(binds[:id0 + 1])
        poly_inds.append(rinf.con_defs[ky][seg_i])

        # ranges >= 80 start+end on borders 0 and 2.
        if rng > 81:
            # ranges > 80 intersect the middle (1) border.
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

        # Collect island boundaries.
        other = rinf.con_defs[ky][seg_i + 1:]
        for itmp, slc in enumerate(other):
            # Add the first index to the end
            slc = other[itmp] = slc + [slc[0]]
            # Plot islands in blue
            ax.plot(rinf.gridxy[0, slc], rinf.gridxy[1, slc], 'b-')
        # Add the islands to main boundary
        out[ky] = [pinds, ] + other
    out['EEZ'] = out['eez'] = [rinf.con_defs['EEZ'][0], ]
    boundary = rinf.gridxy[:, out['EEZ'][0]]
    ax.plot(boundary[0], boundary[1], '-', color='y', alpha=0.5, lw=7)

    return out


def polyinds_at():
    """The Atlantic boundary definitions have some tricky borders.
    """
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
            # hardcode this boundary because the 'border' definitions
            # are a bit strange in the gulf.
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


def polyinds_gm():
    """The Gulf boundary definitions have some tricky borders.
    """
    out = {}
    rinf = RegionInfo('gm')
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

        border_ind = 0
        bxy = rinf.get_contour('borders', xy=True)[border_ind]
        binds = rinf.con_defs['borders'][border_ind]
        id0, d = mindist(xy[:, 0], bxy)
        poly_inds.append(binds[:id0 + 1])
        poly_inds.append(rinf.con_defs[ky][seg_i])

        if rng == 190:
            # hardcode this boundary because the 'border' definitions
            # are a bit strange in the gulf.
            poly_inds[0] += range(228, 259)

        border_ind = 1
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


def polyinds_ec():
    """The east coast borders are as simple as the west coast.
    """
    out = {}
    rinf = RegionInfo('ec')
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

        border_ind = 0
        bxy = rinf.get_contour('borders', xy=True)[border_ind]
        binds = rinf.con_defs['borders'][border_ind]
        id0, d = mindist(xy[:, 0], bxy)
        poly_inds.append(binds[:id0 + 1])
        poly_inds.append(rinf.con_defs[ky][seg_i])

        border_ind = 1
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
    """The Puerto Rico US Virgin Islands boundary definitions do not
    inersect land (all are loops).
    many 'island' boundaries and three 'borders', so there is some
    hardcoding to treat different ranges differently.
    """
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
        border_ind = 0
        bxy = rinf.get_contour('borders', xy=True)[border_ind]
        ax.plot(bxy[0, 0], bxy[1, 0], 'b+')
        binds = rinf.con_defs['borders'][border_ind]
        id0, d = mindist(xy[:, -1], bxy)
        poly_inds.append(rinf.con_defs[ky][seg_i])

        if rng < 171:
            # This joins the next segment
            seg_i += 1
            xy = rinf.get_contour(ky, xy=True)[seg_i]
            ax.plot(xy[0, 0], xy[1, 0], 'rx')
            ide, d = mindist(xy[:, 0], bxy)
            poly_inds.append(binds[id0:ide])
            poly_inds.append(rinf.con_defs[ky][seg_i])

        if rng != 10:
            # The first range loops all the way around.
            # The rest need fall into this `if` and are joined with the
            # first border.
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


def polyinds_hi():
    """This function is mostly a pass-through because the contours
    around the islands are mostly defined in the setup_data.fix_gaps
    function. There is no mainland to loop in here.
    """
    out = {}
    rinf = RegionInfo('hi')
    fig, ax = setup_figure(rinf, 1000 + 32)

    cmap = plt.get_cmap('tab20')

    for rng in range(10, 200, 10):
        ky = '{:03d}'.format(rng)
        seg_i = 0
        xy = rinf.get_contour(ky, xy=True)[seg_i]
        ax.plot(xy[0, 0], xy[1, 0], 'r+')

        poly_inds = rinf.con_defs[ky][seg_i]
        border_ind = 0
        bxy = rinf.get_contour('borders', xy=True)[border_ind]
        ax.plot(bxy[0, 0], bxy[1, 0], 'b+')

        pinds = np.hstack(poly_inds).tolist()
        boundary = rinf.gridxy[:, pinds]
        ax.plot(boundary[0], boundary[1], '-', color=cmap((rng - 5) / 200.),
                alpha=1, zorder=500-rng)

        other = rinf.con_defs[ky][seg_i + 1:]
        for itmp, slc in enumerate(other):
            slc = other[itmp] = slc + [slc[0]]
            ax.plot(rinf.gridxy[0, slc], rinf.gridxy[1, slc], 'b-')
        out[ky] = [pinds, ] + other

    rng = 200
    seg_i = 0
    ky = 'EEZ'
    xy = rinf.get_contour(ky, xy=True)[seg_i]
    ax.plot(xy[0, 0], xy[1, 0], 'r+')
    out['EEZ'] = out['eez'] = rinf.con_defs['EEZ']
    boundary = rinf.gridxy[:, out['EEZ'][0]]
    ax.plot(boundary[0], boundary[1], '-',
            color=cmap((rng - 5) / 200.), alpha=0.5)

    return out


def run_all():
    plt.ioff()
    poly_defs = {}

    poly_defs['wc'] = polyinds_wc()
    poly_defs['ak'] = polyinds_ak()
    poly_defs['at'] = polyinds_at()
    poly_defs['gm'] = polyinds_gm()
    poly_defs['ec'] = polyinds_ec()
    poly_defs['prusvi'] = polyinds_prusvi()
    poly_defs['hi'] = polyinds_hi()
    plt.close('all')
    plt.ion()
    return poly_defs
