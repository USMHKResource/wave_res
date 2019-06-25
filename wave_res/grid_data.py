import os
import cPickle as _pkl
import gzip as gzip
from paths import pkgdir

Dir = os.getcwd()+'\\grid_data\\'

with open(Dir+'Contour_Ranges2.pkl', 'r') as fl:
    con_defs = _pkl.load(fl)
with open(Dir+'Contour_Ranges.pkl', 'r') as fl:
    con_defs_old = _pkl.load(fl)
with open(Dir+'GridLonLat.pkl', 'r') as fl:
    gridlonlat = _pkl.load(fl)
with open(Dir+'FreqBins.pkl', 'r') as fl:
    freqbins = _pkl.load(fl)
with open(Dir+'LandData.pkl', 'r') as fl:
    land_data = _pkl.load(fl)
with open(Dir+'Boundaries.pkl', 'r') as fl:
    bounds = _pkl.load(fl)
with gzip.open(Dir+'DiffTriangles.pkl.gz', 'r') as fl:
    tri_defs = _pkl.load(fl)


