import setpath
from base import RegionInfo
import numpy as np
import shapely.geometry as sg
from shapely.ops import nearest_points
import proj as _proj
import matplotlib.pyplot as plt
proj = _proj.proj['ec']
import create_triangles as tri
import gzip
import cPickle as pkl

def ll2xy(lon, lat):
    return proj.transform_points(_proj.pc, lon, lat)

rinf = RegionInfo('ec')

subregions = ['ec.ne', 'ec.ma', 'ec.se']

lonlat_NE_MA = np.array(
    [[-67.69955, -72.02103],
     [36.66689, 41.41920]])
xy_NE_MA = ll2xy(*lonlat_NE_MA)

line_NE_MA = sg.LineString(xy_NE_MA)

poly_NE = sg.Polygon(ll2xy(
    *np.hstack((lonlat_NE_MA,
                   np.array(
                       [[-71.43, -57.79],
                        [49.14, 41.57]])
                )
               ))
)

plt.plot(*poly_NE.exterior.xy)

lat_VA_NC = 36.5505

inside_ne = np.zeros(rinf.gridxy[0].shape, dtype='bool')
for idx, (x, y) in enumerate(rinf.gridxy.T):
    inside_ne[idx] = poly_NE.contains(sg.Point(x, y))

###########
# Find the points in each region
inds = {}
inds['ec.ne'] = set(np.nonzero(inside_ne)[0])
inds['ec.se'] = set(np.nonzero(lat_VA_NC > rinf.gridlonlat[1])[0])
inds['ec.ma'] = set(np.nonzero(
    ~inside_ne & (rinf.gridlonlat[1] > lat_VA_NC)
)[0])

border = {}
border['ne_ma'] = []
border['ma_se'] = []
def find_nearest_lat(con_def, gridlonlat, lat):
    return con_def[0][np.argmin(np.abs(gridlonlat[1, con_def[0]] - lat))]

def find_nearest_line(con_def, grid, line):
    tmp = sg.LineString(grid[:, con_def[0]].T)
    point = tmp.intersection(line)
    tmp = np.array(tmp.xy)
    delta = (tmp[0] - point.xy[0][0]) ** 2 + (tmp[1] - point.xy[1][0]) ** 2
    return con_def[0][np.argmin(delta)]

kys = rinf.con_defs.keys()
kys.sort()
for ky in kys:
    try:
        int(ky)
    except:
        continue
    else:
        val = find_nearest_lat(rinf.con_defs[ky], rinf.gridlonlat, lat_VA_NC)
        border['ma_se'].append(val)
        val = find_nearest_line(rinf.con_defs[ky], rinf.gridxy, line_NE_MA)
        border['ne_ma'].append(val)

inds['ec.ne'].update(border['ne_ma'])
inds['ec.ma'].update(border['ne_ma'])
inds['ec.ma'].update(border['ma_se'])
inds['ec.se'].update(border['ma_se'])

CON_DEFS_OUT = {}
CON_DEFS_OUT['ec.ne'] = {}
CON_DEFS_OUT['ec.ma'] = {}
CON_DEFS_OUT['ec.se'] = {}

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
                cd_out[ky].append(tmp)
    return cd_out

for ky in CON_DEFS_OUT.keys():
    CON_DEFS_OUT[ky] = clip_contours(rinf.con_defs, inds[ky])


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

with open('Boundaries_ec-subregions.pkl', 'w') as fl:
    pkl.dump(BOUNDS, fl)

with open('Contour_Ranges_ec-subregions.pkl', 'w') as fl:
    pkl.dump(CON_DEFS_OUT, fl)

with open('Land_Data_ec-subregions.pkl', 'w') as fl:
    pkl.dump(LAND_DATA, fl)


# You will need to run the script twice if you change any files above, b/c these are calculated from rinf directly.
TRI_DEFS = {}
for ky in subregions:
    TRI_DEFS[ky] = tri.calc_triangles(ky)
TRI_DEFS_DIFF = tri.run_diff_tri_dict(TRI_DEFS)

    
with gzip.open(str('DiffTriangles_ec-subregions.pkl.gz'), 'w') as fl:
    pkl.dump(TRI_DEFS_DIFF, fl)
