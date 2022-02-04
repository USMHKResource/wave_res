import setpath
import shapely.geometry as sg
from functools import partial
import pyproj
import pycrs
import shapefile
import numpy as np
import proj as _proj
proj = _proj.proj['ak']
from shapely.ops import transform


#project = partial(pyproj.transform, )

#converpyproj.Proj(proj.proj4_init)

sfile = shapefile.Reader('boem_plan_area/AK_PLAN.shp')

geoj = sfile.__geo_interface__


def ll2xy(lon, lat):
    return proj.transform_point(lon, lat, _proj.pc)

class ShapeBoundary2(object):

    def __init__(self, name, shape, proj=proj):
        self.name = name
        self.proj = proj
        self.sg = transform(ll2xy, shape)

boem = {}

for feat in geoj['features']:
    nm = feat['properties']['MMS_PLAN_A']
    boem[nm] = ShapeBoundary2(nm, sg.shape(feat['geometry']))
    boem[nm].full_name = feat['properties']['TEXT_LABEL']
    
