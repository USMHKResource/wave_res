import cPickle as _pkl
import gzip as gzip
from paths import pkgdir

GDDir = pkgdir / 'grid_data/'

with open(str(GDDir / 'Contour_Ranges2.pkl'), 'r') as fl:
    con_defs = _pkl.load(fl)
with open(str(GDDir / 'Contour_Ranges_wc-states.pkl'), 'r') as fl:
    _con_defs = _pkl.load(fl)
con_defs.update(_con_defs)
with open(str(GDDir / 'Contour_Ranges_ec-subregions.pkl'), 'r') as fl:
    _con_defs = _pkl.load(fl)
con_defs.update(_con_defs)
with open(str(GDDir / 'Contour_Ranges_ak-boem-planning-areas.pkl'), 'r') as fl:
    _con_defs = _pkl.load(fl)
con_defs.update(_con_defs)

with open(str(GDDir / 'Contour_Ranges.pkl'), 'r') as fl:
    con_defs_old = _pkl.load(fl)

with open(str(GDDir / 'GridLonLat.pkl'), 'r') as fl:
    gridlonlat = _pkl.load(fl)
with open(str(GDDir / 'FreqBins.pkl'), 'r') as fl:
    freqbins = _pkl.load(fl)

with open(str(GDDir / 'LandData.pkl'), 'r') as fl:
    land_data = _pkl.load(fl)
with open(str(GDDir / 'Land_Data_wc-states.pkl'), 'r') as fl:
    _land_data = _pkl.load(fl)
land_data.update(_land_data)
with open(str(GDDir / 'Land_Data_ec-subregions.pkl'), 'r') as fl:
    _land_data = _pkl.load(fl)
land_data.update(_land_data)

with open(str(GDDir / 'Boundaries.pkl'), 'r') as fl:
    bounds = _pkl.load(fl)
with open(str(GDDir / 'Boundaries_wc-states.pkl'), 'r') as fl:
    _bounds = _pkl.load(fl)
bounds.update(_bounds)
with open(str(GDDir / 'Boundaries_ec-subregions.pkl'), 'r') as fl:
    _bounds = _pkl.load(fl)
bounds.update(_bounds)

# Thesea are used in-place of bounds, when they are present for a region.
with gzip.open(str(GDDir / 'ClippingPaths_ak-boem-planning-areas.pkl.gz'), 'r') as fl:
    clips = _pkl.load(fl)

with gzip.open(str(GDDir / 'DiffTriangles.pkl.gz'), 'r') as fl:
    tri_defs = _pkl.load(fl)
with gzip.open(str(GDDir / 'DiffTriangles_wc-states.pkl.gz'), 'r') as fl:
    _tri_defs = _pkl.load(fl)
tri_defs.update(_tri_defs)
with gzip.open(str(GDDir / 'DiffTriangles_ec-subregions.pkl.gz'), 'r') as fl:
    _tri_defs = _pkl.load(fl)
tri_defs.update(_tri_defs)
with gzip.open(str(GDDir / 'DiffTriangles_ak-boem-planning-areas.pkl.gz'), 'r') as fl:
    _tri_defs = _pkl.load(fl)
tri_defs.update(_tri_defs)
