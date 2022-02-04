import setpath
from base import RegionInfo
import numpy as np
from tools import group
import shapely.geometry as sg
from shapely.ops import nearest_points
import create_triangles as tri
import gzip
import cPickle as pkl
from alaska_boem import boem

import proj as _proj
proj = _proj.proj['ak']

import matplotlib.pyplot as plt
plt.ion()

do_these = ['GOA', 'KOD', 'SHU', 'ALA', # Gulf of AK
#            'COK', # Cook Inlet N/A because no data.
            'BOW', 'GEO', 'NAL', 'ALB', 'NAV', 'MAT', 'NOR'] # Bering Sea
            

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
            tmp = [val for val in cd_now[idx] if val in INDS]
            if len(tmp) > 0:
                if (np.diff(tmp) > 1).any():
                    import pdb
                    pdb.set_trace()
                cd_out[ky].append(tmp)
            
    return cd_out

CON_DEFS_OUT = {}
for ky in inds.keys():
    CON_DEFS_OUT[ky] = clip_contours(rinf.con_defs, inds[ky])


dnow = CON_DEFS_OUT['ak.ALA']

clr = ['r', 'b', 'g', 'y', 'm', 'c']

for idx, ky in enumerate(sorted(dnow)):
    inow = dnow[ky]
    for i in inow:
        ax.plot(*rinf.gridxy[:, i], color=clr[idx % 6], zorder=3)
    #break

error

# Add the borders
# CA
CON_DEFS_OUT['ec.ne']['borders'] = [
    rinf.con_defs['borders'][0],
    border['ne_ma'][::-1]
]
# OR
CON_DEFS_OUT['ec.ma']['borders'] = [
    border['ne_ma'], 
    border['ma_se'][::-1]
]
# WA
CON_DEFS_OUT['ec.se']['borders'] = [
    border['ma_se'],
    rinf.con_defs['borders'][1],
]


LAND_DATA = {}
mainland_xy = ll2xy(*rinf.mainland)
NE_LAND_INDS = np.zeros_like(mainland_xy[:, 0], dtype='bool')
for idx, xy in enumerate(mainland_xy[:, :2]):
    NE_LAND_INDS[idx] = poly_NE.contains(sg.Point(xy))
    
LAND_DATA['ec.ne'] = (
    rinf.mainland[:, NE_LAND_INDS],
    []
)

LAND_DATA['ec.ma'] = (
    rinf.mainland[:, ~NE_LAND_INDS & (rinf.mainland[1] > lat_VA_NC)],
    [] # Initialize empty list, filled below
)
LAND_DATA['ec.se'] = (
    rinf.mainland[:, rinf.mainland[1] < lat_VA_NC],
    [] # Initialize empty list, filled below
)

# Now sort islands
for idx, isl in enumerate(rinf.islands):
    isl_ = sg.Polygon(ll2xy(*isl))
    if poly_NE.intersection(isl_):
        # Double check that it is fully contained in New England
        assert poly_NE.contains(isl_)
        LAND_DATA['ec.ne'][1].append(isl)
        continue
    elif (isl[1] > lat_VA_NC).all():
        LAND_DATA['ec.ma'][1].append(isl)
    elif (isl[1] < lat_VA_NC).all():
        LAND_DATA['ec.se'][1].append(isl)
    else:
        print("One island is on the border?", idx)


BOUNDS = {}
region = 'ec.ne'
bnow = BOUNDS[region] = {}
cdnow = CON_DEFS_OUT[region]

bin = rinf.bounds
for idx, ind in enumerate(range(10, 201, 10)):
    ky = '{:03d}'.format(ind)
    if ind == 200:
        ky = 'eez'
    maxval = cdnow['borders'][1][-(idx + 1)]
    maxind = np.nonzero(~(np.array(bin[ky][0])<maxval))[0][0]
    bnow[ky] = [bin[ky][0][:maxind] + cdnow['borders'][1][-(idx + 1):]]
bnow['200'] = bnow['EEZ'] = bnow['eez']

region = 'ec.se'
bnow = BOUNDS[region] = {}
cdnow = CON_DEFS_OUT[region]

for idx, ind in enumerate(range(10, 201, 10)):
    ky = '{:03d}'.format(ind)
    if ind == 200:
        ky = 'eez'
    minval = cdnow['borders'][0][idx]
    minind = np.nonzero((np.array(bin[ky][0]) >= minval))[0][0]
    bnow[ky] = [cdnow['borders'][0][:idx] + bin[ky][0][minind:]]
bnow['200'] = bnow['EEZ'] = bnow['eez']

region = 'ec.ma'
bnow = BOUNDS[region] = {}
cdnow = CON_DEFS_OUT[region]
n_border = CON_DEFS_OUT[region]['borders'][0]
s_border = CON_DEFS_OUT[region]['borders'][1]
for idx, ind in enumerate(range(10, 201, 10)):
    ky = '{:03d}'.format(ind)
    if ind == 200:
        ky = 'eez'
    
    bnow[ky] = [n_border[:(idx)] + cdnow[ky][0] + s_border[-(idx + 1):][1:]]

bnow['200'] = bnow['EEZ'] = bnow['eez']

# with open('Boundaries_ec-subregions.pkl', 'w') as fl:
#     pkl.dump(BOUNDS, fl)

# with open('Contour_Ranges_ec-subregions.pkl', 'w') as fl:
#     pkl.dump(CON_DEFS_OUT, fl)

# with open('Land_Data_ec-subregions.pkl', 'w') as fl:
#     pkl.dump(LAND_DATA, fl)


# You will need to run the script twice if you change any files above, b/c these are calculated from rinf directly.
TRI_DEFS = {}
for ky in subregions:
    TRI_DEFS[ky] = tri.calc_triangles(ky)
TRI_DEFS_DIFF = tri.run_diff_tri_dict(TRI_DEFS)

    
# with gzip.open(str('DiffTriangles_ec-subregions.pkl.gz'), 'w') as fl:
#     pkl.dump(TRI_DEFS_DIFF, fl)
