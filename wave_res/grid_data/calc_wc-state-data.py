import cPickle as pkl
import setpath
from base import RegionInfo
import numpy as np
import create_triangles as tri
import gzip

rinf = RegionInfo('wc')

lat_WA_OR = 46.2614
lat_OR_CA = 42.0000

#error

###########
# Find the points in each region
inds = {}
inds['wc.wa'] = set(np.nonzero(rinf.gridlonlat[1] > lat_WA_OR)[0])
inds['wc.or'] = set(np.nonzero((lat_WA_OR > rinf.gridlonlat[1]) &
                               (rinf.gridlonlat[1] > lat_OR_CA))[0])
inds['wc.ca'] = set(np.nonzero(rinf.gridlonlat[1] < lat_OR_CA)[0])

###########
# Find the points on the borders
border = {}
border['wa_or'] = []
border['or_ca'] = []
def find_nearest_lat(con_def, gridlonlat, lat):
    return con_def[0][np.argmin(np.abs(gridlonlat[1, con_def[0]] - lat))]

kys = rinf.con_defs.keys()
kys.sort()
for ky in kys:
    try:
        int(ky)
    except:
        continue
    else:
        val = find_nearest_lat(rinf.con_defs[ky], rinf.gridlonlat, lat_WA_OR)
        border['wa_or'].append(val)
        val = find_nearest_lat(rinf.con_defs[ky], rinf.gridlonlat, lat_OR_CA)
        border['or_ca'].append(val)

# Join borders to inds
inds['wc.or'].update(border['wa_or'])
inds['wc.or'].update(border['or_ca'])
inds['wc.wa'].update(border['wa_or'])
inds['wc.ca'].update(border['or_ca'])


CON_DEFS_OUT = {}
CON_DEFS_OUT['wc.wa'] = {}
CON_DEFS_OUT['wc.or'] = {}
CON_DEFS_OUT['wc.ca'] = {}


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

for ky in ['wc.wa', 'wc.or', 'wc.ca']:
    CON_DEFS_OUT[ky] = clip_contours(rinf.con_defs, inds[ky])

# Add the borders
# CA
CON_DEFS_OUT['wc.ca']['borders'] = [
    rinf.con_defs['borders'][0],
    border['or_ca'][::-1]
]
# OR
CON_DEFS_OUT['wc.or']['borders'] = [
    border['or_ca'],
    border['wa_or'][::-1]
]
# WA
CON_DEFS_OUT['wc.wa']['borders'] = [
    border['wa_or'],
    rinf.con_defs['borders'][1],
]


LAND_DATA = {}
LAND_DATA['wc.wa'] = (
    rinf.mainland[:, rinf.mainland[1] > lat_WA_OR],
    [] # No Islands in WA
)
LAND_DATA['wc.or'] = [
    rinf.mainland[:, (lat_WA_OR > rinf.mainland[1]) & (rinf.mainland[1] > lat_OR_CA)],
    [] # No Islands in OR
]
# Drop the first 11 points, they go up the CR, and mess up clipping later.
LAND_DATA['wc.or'][0] = LAND_DATA['wc.or'][0][:, 11:]
LAND_DATA['wc.or'] = tuple(LAND_DATA['wc.or'])

LAND_DATA['wc.ca'] = (
    rinf.mainland[:, lat_OR_CA > rinf.mainland[1]],
    rinf.islands
)


BOUNDS = {}
bnow = BOUNDS['wc.ca'] = {}

bin = rinf.bounds
kys = bin.keys()
kys.sort()
for idx, ind in enumerate(range(10, 201, 10)):
    ky = '{:03d}'.format(ind)
    if ind == 200:
        ky = 'eez'
    maxval = CON_DEFS_OUT['wc.ca']['borders'][1][-(idx + 1)]
    maxind = np.nonzero(~(np.array(bin[ky][0])<maxval))[0][0]
    bnow[ky] = [bin[ky][0][:maxind] + CON_DEFS_OUT['wc.ca']['borders'][1][-(idx + 1):]]

bnow['200'] = bnow['EEZ'] = bnow['eez']


bnow = BOUNDS['wc.wa'] = {}

bin = rinf.bounds
kys = bin.keys()
kys.sort()
for idx, ind in enumerate(range(10, 201, 10)):
    ky = '{:03d}'.format(ind)
    if ind == 200:
        ky = 'eez'
    minval = CON_DEFS_OUT['wc.wa']['borders'][0][idx]
    minind = np.nonzero((np.array(bin[ky][0]) >= minval))[0][0]
    bnow[ky] = [CON_DEFS_OUT['wc.wa']['borders'][0][:idx] + bin[ky][0][minind:]]
bnow['200'] = bnow['EEZ'] = bnow['eez']

bnow = BOUNDS['wc.or'] = {}
con = CON_DEFS_OUT['wc.or']
s_border = CON_DEFS_OUT['wc.or']['borders'][0]
n_border = CON_DEFS_OUT['wc.or']['borders'][1]
for idx, ind in enumerate(range(10, 201, 10)):
    ky = '{:03d}'.format(ind)
    if ind == 200:
        ky = 'eez'
    
    bnow[ky] = [s_border[:(idx)] + con[ky][0] + n_border[-(idx + 1):][1:]]

bnow['200'] = bnow['EEZ'] = bnow['eez']

with open('Boundaries_wc-states.pkl', 'w') as fl:
    pkl.dump(BOUNDS, fl)

with open('Contour_Ranges_wc-states.pkl', 'w') as fl:
    pkl.dump(CON_DEFS_OUT, fl)

with open('Land_Data_wc-states.pkl', 'w') as fl:
    pkl.dump(LAND_DATA, fl)

TRI_DEFS = {}
for ky in ['wc.wa', 'wc.or', 'wc.ca']:
    TRI_DEFS[ky] = tri.calc_triangles(ky)
TRI_DEFS_DIFF = tri.run_diff_tri_dict(TRI_DEFS)

    
with gzip.open(str('DiffTriangles_wc-states.pkl.gz'), 'w') as fl:
    pkl.dump(TRI_DEFS_DIFF, fl)
