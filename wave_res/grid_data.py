import cPickle as _pkl
import gzip as gzip
from paths import pkgdir

GDDir = pkgdir / 'grid_data/'

with open(str(GDDir / 'Contour_Ranges2.pkl'), 'r') as fl:
    con_defs = _pkl.load(fl)
with open(str(GDDir / 'Contour_Ranges.pkl'), 'r') as fl:
    con_defs_old = _pkl.load(fl)
with open(str(GDDir / 'GridLonLat.pkl'), 'r') as fl:
    gridlonlat = _pkl.load(fl)
with open(str(GDDir / 'FreqBins.pkl'), 'r') as fl:
    freqbins = _pkl.load(fl)
with open(str(GDDir / 'LandData.pkl'), 'r') as fl:
    land_data = _pkl.load(fl)
with open(str(GDDir / 'Boundaries.pkl'), 'r') as fl:
    bounds = _pkl.load(fl)
with gzip.open(str(GDDir / 'DiffTriangles.pkl.gz'), 'r') as fl:
    tri_defs = _pkl.load(fl)


