import numpy as np
import proj
import grid_data as gdat

regions = {'wc': 'wc',
           'wc.or': 'wc', 'wc.wa': 'wc', 'wc.ca': 'wc',
           'ec': 'at', 'gm': 'at', 'at': 'at',
           'ec.ne': 'at', 'ec.ma': 'at', 'ec.se': 'at',
           'hi': 'hi',
           'ak': 'ak',
           'prusvi': 'prusvi'}

# Use set() to get unique values.
source_regions = list(set([sr for sr in regions.itervalues()]))

conids = ['EEZ'] + ['{:03d}'.format(n) for n in range(10, 200, 10)]

con_defs = gdat.con_defs
# This is the outer-boundary of the EEZ (not including the Canada
# + Mexico borders)
con_defs['wc']['eez'] = [range(34, 148), ]
con_defs['wc']['borders'] = [range(0, 35), range(147, 176)]

con_defs['at']['eez'] = [range(32, 147), range(227, 282)]
con_defs['at']['borders'] = [range(0, 33), range(133, 228), range(282, 303)]

con_defs['ec']['eez'] = [range(32, 147), ]
con_defs['ec']['borders'] = [range(0, 33), range(133, 191) + [849, 499]]
con_defs['gm']['eez'] = [range(227, 282), ]
con_defs['gm']['borders'] = [[499, 849] + range(190, 228), range(282, 303)]
con_defs['prusvi']['eez'] = [range(0, 7), ]
con_defs['prusvi']['borders'] = [range(7, 117) + [0], ]
con_defs['ak']['eez'] = [range(48, 392), range(439, 541)]
con_defs['ak']['borders'] = [range(0, 49), range(391, 440), range(540, 616)]

# These are from GGM's def of HI EEZ (it's a bit wonky b/c the HI EEZ
# includes the other islands in the HI chain (e.g., Midway)).
con_defs['hi']['eez'] = [np.r_[range(0,67,1),
                               np.array([3081, 2907, 2736, 2567, 2400, 2236,
                                         2075, 1607, 1456, 1221, 1077, 1078,
                                         1079, 1080, 1081, 1082, 1083, 1084,
                                         1085, 1086, 1087, 1088, 1090, 1091,
                                         950,  949,  948,  947,  946, 1068,
                                         1069, 1070, 1071, 1072, 1073, 1074,
                                         1075, 1076, 1220, 1455, 1760, 1916,
                                         2074, 2566, 2735, 3080]),
                               range(340,453,1),
                               np.array([0])].tolist(),
                         ]
con_defs['hi']['borders'] = [np.array([3081, 2907, 2736, 2567, 2400, 2236,
                                       2075, 1607, 1456, 1221, 1077, 1078,
                                       1079, 1080, 1081, 1082, 1083, 1084,
                                       1085, 1086, 1087, 1088, 1090, 1091,
                                       950,  949,  948,  947,  946, 1068,
                                       1069, 1070, 1071, 1072, 1073, 1074,
                                       1075, 1076, 1220, 1455, 1760, 1916,
                                       2074, 2566, 2735, 3080]),]


class RegionInfo(object):
    """A class for region info.
    """

    mainland = None
    islands = None

    _proj_pc = proj.pc

    def __init__(self, region, use_old_con_defs=False):
        # 'region' must be in regions.
        if region not in regions:
            raise Exception("Invalid region '{}', please choose from "
                            .format(region) +
                            ("{}, " * len(regions)).format(*regions)[:-2] + '.')
        self.source_region = regions[region]
        self.region = region
        self.gridlonlat = gdat.gridlonlat[self.source_region]
        if use_old_con_defs:
            self.con_defs = gdat.con_defs_old[region]
        else:
            self.con_defs = con_defs[region]
        if '200' not in self.con_defs and 'eez' in self.con_defs:
            self.con_defs['200'] = self.con_defs['eez']
        self.freqbins = gdat.freqbins[self.source_region]
        if region in proj.proj:
            self.proj = proj.proj[region]
            self.gridxy = self.transform(self.gridlonlat)
        if region in gdat.land_data:
            self.mainland, self.islands = gdat.land_data[region]
        if region in gdat.bounds:
            self.bounds = gdat.bounds[region]
        if region in gdat.tri_defs:
            self.tri_inds = gdat.tri_defs[region]

    @property
    def allxy(self, ):
        """This is the spatial grid including land (in projection coordinates).
        """
        allxy = [self.gridxy, ]
        if self.mainland is not None:
            allxy.append(self.transform(self.mainland))
        for isl in self.islands:
            allxy.append(self.transform(isl))
        return np.hstack(allxy)

    def get_contour(self, conid, xy=False):
        """conid must be in '010', '020', ... '200', or 'EEZ' """
        if xy:
            data = self.gridxy
        else:
            data = self.gridlonlat
        if isinstance(self.con_defs[conid], list):
            out = []
            for slc in self.con_defs[conid]:
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
