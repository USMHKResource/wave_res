import cPickle as _pkl
import paths as p

regions = ['wc', 'at', 'hi', 'ak', 'prusvi']
conids = ['EEZ'] + ['{:03d}'.format(n) for n in range(10, 200, 10)]

with open(str(p.projdir / 'data/Contour_Ranges.pkl'), 'r') as fl:
    con_defs = _pkl.load(fl)
with open(str(p.projdir / 'data/GridLonLat.pkl'), 'r') as fl:
    gridlonlat = _pkl.load(fl)
with open(str(p.projdir / 'data/FreqBins.pkl'), 'r') as fl:
    freqbins = _pkl.load(fl)


class RegionInfo(object):
    """A class for region info.
    """
    def __init__(self, region):
        # 'region' must be in regions.
        if region not in regions:
            raise Exception("Invalid region '{}', please choose from ".format(region) +
                            ("{}, " * len(regions)).format(*regions)[:-2] + '.')
        self.region = region
        self.gridlonlat = gridlonlat[region]
        self.con_defs = con_defs[region]
        self.freqbins = freqbins[region]

    def get_lonlat(self, conid):
        """conid must be in '010', '020', ... '200', or 'EEZ' """
        return self.gridlonlat[:, self.con_defs[conid]]
        
# This is the outer-boundary of the EEZ (not including the Canada
# + Mexico borders)
con_defs['wc']['EEZ'] = slice(34, 148)
