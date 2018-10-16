import cPickle as _pkl
import paths as p

regions = {'wc': 'wc',
           'ec': 'at', 'gm': 'at', 'at': 'at',
           'hi': 'hi',
           'ak': 'ak',
           'prusvi': 'prusvi'}

# Use set() to get unique values.
source_regions = list(set([sr for sr in regions.iteritems()]))

conids = ['EEZ'] + ['{:03d}'.format(n) for n in range(10, 200, 10)]


with open(str(p.projdir / 'data/Contour_Ranges.pkl'), 'r') as fl:
    con_defs = _pkl.load(fl)
with open(str(p.projdir / 'data/GridLonLat.pkl'), 'r') as fl:
    gridlonlat = _pkl.load(fl)
with open(str(p.projdir / 'data/FreqBins.pkl'), 'r') as fl:
    freqbins = _pkl.load(fl)
        
# This is the outer-boundary of the EEZ (not including the Canada
# + Mexico borders)
con_defs['wc']['eez'] = slice(34, 148)

con_defs['ec']['eez'] = slice(32, 147)
con_defs['gm']['eez'] = slice(227, 282)
con_defs['prusvi']['eez'] = slice(0, 7)
con_defs['ak']['eez'] = [slice(48, 392), slice(439, 541)]

class RegionInfo(object):
    """A class for region info.
    """
    def __init__(self, region):
        # 'region' must be in regions.
        if region not in regions:
            raise Exception("Invalid region '{}', please choose from ".format(region) +
                            ("{}, " * len(regions)).format(*regions)[:-2] + '.')
        self.source_region = regions[region]
        self.region = region
        self.gridlonlat = gridlonlat[self.source_region]
        self.con_defs = con_defs[region]
        self.freqbins = freqbins[self.source_region]

    def get_lonlat(self, conid):
        """conid must be in '010', '020', ... '200', or 'EEZ' """
        if isinstance(self.con_defs[conid], list):
            out = []
            for slc in self.con_defs[conid]:
                out.append(self.gridlonlat[:, slc])
        else:
            out = [self.gridlonlat[:, self.con_defs[conid]], ]
        return out
