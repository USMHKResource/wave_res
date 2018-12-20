import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
import shapely.geometry as sg
import numpy as np
from copy import copy
import setpath
from base import RegionInfo
import proj


def rinf2clip(rinf, con):
    isla = []
    for isl in rinf.islands:
        isla.append(sg.Polygon(rinf.transform(isl).T))
    polygons = []
    for id0, inds in enumerate(rinf.bounds[con]):
        xy = rinf.gridxy[:, inds]
        if id0 == 0 and inds[0] != inds[-1]:
            xy = np.hstack((
                xy,
                rinf.transform(rinf.mainland)))
        p = sg.Polygon(shell=xy.T)
        itmp = []
        for isl in isla:
            if isl.within(p):
                itmp.append(np.array(isl.boundary.xy).T)
            p = sg.Polygon(shell=xy.T, holes=itmp)
        polygons.append(p)
    return sg.MultiPolygon(polygons)


def clip_tri(clip, tri, inplace=False):

    if not inplace:
        tri = copy(tri)
    tri_index = np.ones(tri.vertices.shape[0], 'bool')
    for id0, v in enumerate(tri.vertices):
        tnow = sg.Polygon(tri.points[v])
        if not clip.contains(tnow):
            tri_index[id0] = False
    tri.vertices = tri.vertices[tri_index]
    tri.simplices = tri.simplices[tri_index]
    return tri


def plot_clip(clip, ax=None, **kwargs):
    if ax is None:
        ax = plt.gca()
    if isinstance(clip, sg.MultiPolygon):
        for c in clip:
            plot_clip(c, ax=ax, **kwargs)
        return
    x, y = clip.exterior.coords.xy
    ax.plot(x, y, **kwargs)
    for i in clip.interiors:
        x, y = i.coords.xy
        ax.plot(x, y, **kwargs)


def plot_tri(tri, ax=None, **kwargs):
    if ax is None:
        ax = plt.gca()
    x, y = tri.points[tri.vertices, :].T
    x = x[[0, 1, 2, 0], :]
    y = y[[0, 1, 2, 0], :]
    ax.plot(x, y, **kwargs)


def tri_diff(tri1, tri2):
    if isinstance(tri1, Delaunay):
        out = copy(tri1)
        tri1 = tri1.vertices
    else:
        out = False
    if isinstance(tri2, Delaunay):
        tri2 = tri2.vertices
    v1 = {tuple(v) for v in tri1}
    v2 = {tuple(v) for v in tri2}
    vo = np.array(list(v1 - v2))
    if out:
        out.vertices = vo
        out.simplices = vo
    else:
        out = vo
    return out


def calc_triangles(region):
    """Calculate the vertices (a.k.a. simplices) indices for `region`.
    """

    print("Calculating triangles for '{}' region...".format(region.upper()))
    # Get the relavent data for this region.
    rinf = RegionInfo(region)
    tri = Delaunay(rinf.allxy.T)

    # Return a dictionary of 'simplices' indexes to the 'allxy' points
    triout = {}

    for rng in range(10, 201, 10):
        print("  {} nm".format(rng))
        # The 'range key'
        rky = '{:03d}'.format(rng)
        if rng == 200:
            rky = 'EEZ'

        # get the clipping path
        clip = rinf2clip(rinf, rky)
        # clip the triangles
        tri_clip = clip_tri(clip, tri)
        triout[rky] = tri_clip.simplices
        trilast = copy(tri_clip)
    print('Done.')
    triout['200'] = triout['eez'] = triout['EEZ']
    return triout


def run_all(plot=False):
    tri_defs = {}
    for idr, region in enumerate(['wc', 'ec', 'gm', 'at', 'ak', 'prusvi']):
        if plot:
            fignum = False
        else:
            fignum = 1200 + idr
        tri_defs[region] = calc_triangles(region, fignum)
    return tri_defs


def run_diff_tri_dict(tri_dict):

    # tri_dict_diff
    tdd = {}
    for rgn in tri_dict:
        tdd[rgn] = {}
        tdd[rgn]['010'] = lastrgn = tri_dict[rgn]['010']
        for r in range(20, 201, 10):
            rky = '{:03d}'.format(r)
            if r == 200:
                rky = 'EEZ'
            td_now = tri_dict[rgn][rky]
            tdd[rgn][rky] = tri_diff(td_now, lastrgn)
            lastrgn = td_now
        # Set '200' and 'eez' defs to match 'EEZ'
        tdd[rgn]['eez'] = tdd[rgn]['200'] = tdd[rgn]['EEZ']
    return tdd


def show_grid(region, fignum):

    rinf = RegionInfo(region)
    prj = rinf.proj
    tri = Delaunay(rinf.allxy.T)

    # Setup the figure.
    fig = plt.figure(fignum)
    fig.clf()
    ax = plt.axes(projection=prj)
    ax.add_feature(proj.land)
    ax.add_feature(proj.states)
    ax.set_extent((prj.lonlim + prj.latlim), crs=proj.pc)

    plot_tri(tri,
             color='0.7', alpha=0.2,
             ax=ax)

    clrs = plt.get_cmap('tab20')

    for irng, rng in enumerate(range(10, 201, 10)):
        rky = '{:03d}'.format(rng)
        if rky == '200':
            rky = 'EEZ'
        tmp = copy(tri)
        inds = rinf.tri_inds[rky]
        tmp.vertices = inds
        tmp.simplices = inds
        plot_tri(tmp,
                 color=clrs((rng - 5) / 200.),
                 ax=ax)
    return fig, ax


if __name__ == '__main__':

    for idreg, region in enumerate(['wc', 'ec', 'gm', 'at', 'ak', 'prusvi']):
        fig, ax = show_grid(region, 2000 + idreg)
        fig.savefig('fig/RegionDef-{}-01.png'.format(region), dpi=300)
