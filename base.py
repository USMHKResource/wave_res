import cPickle as _pkl
import paths as p
import numpy as np
import proj

regions = {'wc': 'wc',
           'ec': 'at', 'gm': 'at', 'at': 'at',
           'hi': 'hi',
           'ak': 'ak',
           'prusvi': 'prusvi'}

# Use set() to get unique values.
source_regions = list(set([sr for sr in regions.itervalues()]))

conids = ['EEZ'] + ['{:03d}'.format(n) for n in range(10, 200, 10)]


with open(str(p.projdir / 'data/Contour_Ranges2.pkl'), 'r') as fl:
    con_defs = _pkl.load(fl)
with open(str(p.projdir / 'data/Contour_Ranges.pkl'), 'r') as fl:
    con_defs_old = _pkl.load(fl)
with open(str(p.projdir / 'data/GridLonLat.pkl'), 'r') as fl:
    gridlonlat = _pkl.load(fl)
with open(str(p.projdir / 'data/FreqBins.pkl'), 'r') as fl:
    freqbins = _pkl.load(fl)
with open(str(p.projdir / 'data/LandData.pkl'), 'r') as fl:
    land_data = _pkl.load(fl)
        
# This is the outer-boundary of the EEZ (not including the Canada
# + Mexico borders)
con_defs['wc']['eez'] = slice(34, 148)

con_defs['ec']['eez'] = slice(32, 147)
con_defs['gm']['eez'] = slice(227, 282)
con_defs['prusvi']['eez'] = slice(0, 7)
con_defs['ak']['eez'] = [slice(48, 392), slice(439, 541)]

# These are from GGM's def of HI EEZ (it's a bit wonky b/c the HI EEZ
# includes the other islands in the HI chain (e.g., Midway)).
con_defs['hi']['eez'] = [np.r_[np.arange(0,67,1),
                               np.array([3081, 2907, 2736, 2567, 2400, 2236,
                                         2075, 1607, 1456, 1221, 1077, 1078,
                                         1079, 1080, 1081, 1082, 1083, 1084,
                                         1085, 1086, 1087, 1088, 1090, 1091,
                                         950,  949,  948,  947,  946, 1068,
                                         1069, 1070, 1071, 1072, 1073, 1074,
                                         1075, 1076, 1220, 1455, 1760, 1916,
                                         2074, 2566, 2735, 3080]),
                               np.arange(340,453,1),
                               np.array([0])],
                         ]


class RegionInfo(object):
    """A class for region info.
    """

    mainland = None
    islands = None
    
    _proj_pc = proj.pc

    def __init__(self, region, use_old_con_defs=False):
        # 'region' must be in regions.
        if region not in regions:
            raise Exception("Invalid region '{}', please choose from ".format(region) +
                            ("{}, " * len(regions)).format(*regions)[:-2] + '.')
        self.source_region = regions[region]
        self.region = region
        self.gridlonlat = gridlonlat[self.source_region]
        if use_old_con_defs:
            self.con_defs = con_defs_old[region]
        else:
            self.con_defs = con_defs[region]
        self.freqbins = freqbins[self.source_region]
        if region in proj.proj:
            self.proj = proj.proj[region]
            self.gridxy = self.transform(self.gridlonlat)
        if region in land_data:
            self.mainland, self.islands = land_data[region]

    def get_contour(self, conid, xy=False):
        """conid must be in '010', '020', ... '200', or 'EEZ' """
        if xy:
            data = self.gridxy
        else:
            data = self.gridlonlat
        if isinstance(self.con_defs[conid], list):
            out = []
            for slc in self.con_defs[conid]:
                if isinstance(slc, slice) and slc.step == 1j:
                    slc = range(slc.start, slc.stop) + [slc.start]
                out.append(data[:, slc])
        else:
            out = [data[:, self.con_defs[conid]], ]
        return out

    def transform(self, points, inverse=False):
        """Transform points from lonlat to xy coordinates, 
        or the reverse (if inverse==True).
        """
        if inverse:
            pout = proj.pc
            pin = self.proj
        else:
            pout = self.proj
            pin = proj.pc
        return pout.transform_points(pin, points[0], points[1]).T[:2]
