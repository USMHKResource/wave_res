import setpath
from base import RegionInfo
import numpy as np
from tools import group
import shapely.geometry as sg
from shapely.ops import nearest_points
import create_triangles as tri
import gzip
import cPickle as pkl
from alaska_boem import boem, LL2XY

import proj as _proj
proj = _proj.proj['ak']

import matplotlib.pyplot as plt
plt.ion()


do_these = ['GOA', 'KOD', 'SHU', 'ALA', # Gulf of AK
#            'COK', # Cook Inlet N/A because no data.
            'BOW', 'GEO', 'NAL', 'ALB', 'NAV', 'MAT', 'NOR'] # Bering Sea
    

ll2xy = LL2XY(proj)
    
def setup_figure(rinf, fignum, **kwargs):
    # First setup the figure
    prj = rinf.proj
    fig = plt.figure(fignum, **kwargs)
    fig.clf()
    ax = plt.axes(projection=prj)
    ax.add_feature(_proj.land)
    ax.add_feature(_proj.states)
    ax.set_extent((prj.lonlim + prj.latlim), crs=rinf._proj_pc)
    return fig, ax


rinf = RegionInfo('ak')



if True:

    fig, ax = setup_figure(rinf, 102, figsize=[10, 8])

    inside = {}
    inds = {}

    for ky in boem:
        ak_reg = boem[ky]
        pargs = dict(facecolor='none', zorder=10)
        if ky in do_these:
            pargs['edgecolor'] = 'b'
        else:
            pargs['edgecolor'] = '0.5'
        ax.add_geometries([ak_reg.sg], crs=proj, **pargs)
        #ax.plot(*ak_reg.sg.exterior.xy)
        center = ak_reg.sg.centroid.xy
        ax.text(center[0][0], center[1][0], ky, ha='center', va='center')
        inside['ak.' + ky] = np.zeros(rinf.gridxy[0].shape, dtype='bool')
        for idx, (x, y) in enumerate(rinf.gridxy.T):
            inside['ak.' + ky][idx] = ak_reg.sg.contains(sg.Point(x, y))

        inds['ak.' + ky] = set(np.nonzero(inside['ak.' + ky])[0])

    for cd in rinf.con_defs:
        for id, idxs in enumerate(rinf.con_defs[cd]):
            if id == 0:
                color = '0.75'
            else:
                color = '0.85'
            ax.plot(*rinf.gridxy[:, idxs], color=color, zorder=3)

 
def create_cons(inds):
    if len(inds) == 0:
        return None
    dtmp = np.diff(inds)
    out = [[inds[0]], ]
    il = 0
    for idx, d in enumerate(dtmp):
        id1 = idx + 1
        if d == 1:
            out[il] += [inds[id1], ]
        elif d == 2:
            out[il] += [inds[id1] - 1, inds[id1]]
        # elif d == 3:
        #     out[il] += [inds[id1] - 2, inds[id1] - 1, inds[id1]]
        elif d < 0 and idx == (len(dtmp) - 1):
            # This is closing an island loop
            out[il] += [inds[id1]]
        else:
            #out[il] += [inds[idx] + 1]
            il += 1
            out.append([inds[id1]])
    return out
 
def clip_contours(cd_in, INDS):
    kys = cd_in.keys()
    kys.sort()
    cd_out = {}
    for ky in kys:
        if ky.lower() in ['borders']:
            continue
        cd_now = cd_in[ky]
        cd_out[ky] = []
        for idx in range(len(cd_now)):
            tmp = create_cons([val for val in cd_now[idx] if val in INDS])
            if tmp is not None:
                cd_out[ky] += tmp
            # if idx > 0:
            #     if (cd_now[idx][0] in INDS) and (cd_now[idx][-1] in INDS):
            #         import pdb
            #         pdb.set_trace()
    return cd_out

def show_cons_region(region='ak.ALA'):

    dnow = CON_DEFS_OUT[region]

    clr = ['r', 'b', 'g', 'y', 'm', 'c']

    lines = []
    
    for idx, ky in enumerate(sorted(dnow)):
        inow = dnow[ky]
        for i in inow:
            lines += ax.plot(*rinf.gridxy[:, i], color=clr[idx % 6], zorder=3)
    return lines


def hide_lines(lines, unhide=False):
    for l in lines:
        l.set_visible(unhide)


plot_lines = {}
CON_DEFS_OUT = {}
for ky in inds.keys():
    CON_DEFS_OUT[ky] = clip_contours(rinf.con_defs, inds[ky])
    plot_lines[ky] = show_cons_region(ky)

for reg in plot_lines:
    if reg in ['ak.BOW', 'ak.SHU', 'ak.GOA', 'ak.GEO', 'ak.NOR']:
        #Only plot these ones because we want to see that the boundaries don't overlap
        hide_lines(plot_lines[reg])
    else:
        hide_lines(plot_lines[reg], unhide=True)

# with open('Boundaries_ec-subregions.pkl', 'w') as fl:
#     pkl.dump(BOUNDS, fl)

with open('Contour_Ranges_ak-boem-planning-areas.pkl', 'w') as fl:
    pkl.dump(CON_DEFS_OUT, fl)



#########
# CREATE CLIPPING PATHS for clipping delaunay triangles.
 
def make_clip(rng='010'):
    pgon0 = []
    isls = []
    for inds in rinf.con_defs[rng]:
        if inds[0] != inds[-1]:
            pgon0 += inds
        else:
            isls.append(inds)
    xy = rinf.gridxy[:, pgon0]
    extra_xy = ll2xy.transform_points(np.array([[-168.98, 65.64, ],
                                                [-136.91, 64.26],
                                                [-130.95, 54.85]]).T)[:, :2].T
    pgon0 = sg.Polygon(np.hstack((xy, extra_xy)).T)
    pg_isls = []
    for isl in isls:
        pg_isls.append(sg.Polygon(rinf.gridxy[:, isl].T))
    return sg.MultiPolygon([pgon0] + pg_isls)

def do_make_clip():
    """This creates the clipping paths/masks that are used to identify each 
    """
    
    
    CLIP_MASK_OUT = {}

    range_mask = {}
    for rng in range(10, 201, 10):
        ky = '{:03d}'.format(rng)
        range_mask[ky] = make_clip(ky).buffer(1)

    ## For some reason these Multi-Polygons have some bad polygons at the end...
    #range_mask['050'] = range_mask['050'][:-1]
    #
    #range_mask['020'] = range_mask['020'].buffer(1)

    for reg in do_these:
        msk = CLIP_MASK_OUT['ak.' + reg] = {}
        for rng in range(10, 201, 10):
            ky = '{:03d}'.format(rng)
            clip = range_mask[ky].intersection(boem[reg].sg)
            msk[ky] = clip


    with gzip.open('ClippingPaths_ak-boem-planning-areas.pkl.gz', 'w') as fl:
        pkl.dump(CLIP_MASK_OUT, fl)

    return CLIP_MASK_OUT

#### UNCOMMENT THIS LINE TO RECREATE CLIPPING PATHS:
CLIP_MASK_OUT = do_make_clip()


# You will need to run this block twice if you change any files above,
# b/c these are calculated from rinf directly.

if False:

    TRI_DEFS = {}
    for ky in do_these:
        TRI_DEFS['ak.' + ky] = tri.calc_triangles('ak.' + ky)
    TRI_DEFS_DIFF = tri.run_diff_tri_dict(TRI_DEFS)


    with gzip.open(str('DiffTriangles_ak-boem-planning-areas.pkl.gz'), 'w') as fl:
        pkl.dump(TRI_DEFS_DIFF, fl)

if __name__ == '__main__':

        
    for idx, ky in enumerate(do_these):
        fig, ax = tri.show_grid('ak.' + ky, 2000 + idx)
