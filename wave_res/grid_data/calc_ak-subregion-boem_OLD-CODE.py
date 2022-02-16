"""
This script has old code for locating 'borders' between boem planning areas.
"""

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


##########
# THIS CODE is for finding the 'borders' between regions, but it's not complete, and it doesn't work very well, and I think it's unnecessary because I've decided to use the 'clipping paths' instead.

def find_nearest_lat(con_def, gridlonlat, lat):
    return con_def[0][np.argmin(np.abs(gridlonlat[1, con_def[0]] - lat))]

def find_nearest_line(con_def, grid, line):
    tmp = sg.LineString(grid[:, con_def[0]].T)
    point = tmp.intersection(line)
    tmp = np.array(tmp.xy)
    delta = (tmp[0] - point.xy[0][0]) ** 2 + (tmp[1] - point.xy[1][0]) ** 2
    return con_def[0][np.argmin(delta)]

def find_nearest_line2(con_def, grid, line):
    _tmp = sg.LineString(grid[:, con_def[0]].T)
    point = _tmp.intersection(line)
    if point.is_empty:
        return None
    tmp = np.array(_tmp.xy)
    if isinstance(point, sg.MultiPoint):
        pts = point.geoms
    else:
        pts = [point, ]
    out = []
    for p in pts:
        delta = (tmp[0] - p.xy[0][0]) ** 2 + (tmp[1] - p.xy[1][0]) ** 2
        out.append(con_def[0][np.argmin(delta)])
    return out

border = {}
for id0, ky0 in enumerate(do_these):
    sgeom = boem[ky0].sg
    for id1, ky1 in enumerate(do_these):
        if ky0 == ky1 or (ky1, ky0) in border:
            continue
        if isinstance(sgeom, sg.multipolygon.MultiPolygon):
            mline = []
            for gm in sgeom:
                line = gm.buffer(100).intersection(boem[ky1].sg.boundary)
                if line.length > 0:
                    mline.append(line)
            if len(mline) == 0:
                continue
            elif len(mline) > 1:
                mline = sg.shape(mline)
            else:
                mline = mline[0]
        else:
            mline = sgeom.buffer(100).intersection(boem[ky1].sg.boundary)
        bnow = border[(ky0, ky1)] = []
        for cdky in sorted(rinf.con_defs.keys()):
            try:
                int(cdky)
            except:
                continue
            val = find_nearest_line2(rinf.con_defs[cdky], rinf.gridxy, mline)
            if val is not None:
                bnow += [val[0], ]
            ax.plot(*rinf.gridxy[:, bnow], color='r')
    #     if id1 > 2:
    #         break
    # if id0 > 3:
    #     break
