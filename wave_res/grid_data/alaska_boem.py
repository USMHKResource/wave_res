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


class LL2XY(object):

    def __init__(self, proj_out, proj_in=_proj.pc):
        self._prj_in = proj_in
        self._prj_out = proj_out
        self._ifunc = lambda x, y: self._prj_in.transform_point(x, y, self._prj_out)
        self._func = lambda x, y: self._prj_out.transform_point(x, y, self._prj_in)
        
    def transform_shape(self, shape, inverse=False):
        from shapely.ops import transform as _tr_
        if inverse:
            func = self._ifunc
        else:
            func = self._func
        return _tr_(func, shape)

    def transform_points(self, points, inverse=False):
        if inverse:
            return self._prj_in.transform_points(self._prj_out, *points)
        else:
            return self._prj_out.transform_points(self._prj_in, *points)


class ShapeBoundary2(object):

    def __init__(self, name, shape, proj=proj):
        self._ll2xy = LL2XY(proj)
        self.name = name
        self.proj = proj
        self.sg = self._ll2xy.transform_shape(shape)

boem = {}

for feat in geoj['features']:
    nm = feat['properties']['MMS_PLAN_A']
    boem[nm] = ShapeBoundary2(nm, sg.shape(feat['geometry']))
    boem[nm].full_name = feat['properties']['TEXT_LABEL']
    
