import setpath
import shapely.geometry as sg
import pyproj
import pycrs
import shapefile
import numpy as np
import proj as _proj
proj = _proj.proj['ak']

abbrev_defs = {'Koniag Inc': 'Koniag',
               'The Aleut Corp': 'Aleut', 
               'Chugach Alaska Corp': 'Chugach',
               'Ahtna Inc': 'Ahtna', 
               'Sealaska Corp': 'Sealaska',
               'Arctic Slope Regional Corp': 'Arctic',
               'Calista Corp': 'Calista',
               'Doyon Ltd': 'Doyon',
               'Bristol Bay Native Corp': 'BBNC',
               'Cook Inlet Region Inc': 'CIRI',
               'N.A.N.A. Regional Corp': 'NANA',
               'Bering Straits Native Corp': 'BSNC',
               }
abbrevs = abbrev_defs.values()


# This converts the data in the shapefile to lat/lon
crs = pycrs.load.from_file('ansca_corp/ansca_corp.prj')
proj0 = pyproj.Proj(crs.to_proj4())

sfile = shapefile.Reader('ansca_corp/ansca_corp.shp')

ancsa = {}

class ShapeBoundary(object):

    def __init__(self, name, lonlat, proj=proj):
        self.name = name
        self.lonlat = lonlat
        self.proj = proj
        self.xy = np.array(proj.transform_points(_proj.pc, lonlat[:, 0], lonlat[:, 1]))
        self.sg = sg.Polygon(self.xy)

for sr in sfile.shapeRecords():
    nm = abbrev_defs[sr.record[1]]
    xy = np.array(sr.shape.points).copy()
    lonlat = np.array(proj0(xy[:, 0],xy[:, 1], inverse=True)).T
    ancsa[nm] = ShapeBoundary(nm, lonlat)
    ancsa[nm].full_name = sr.record[1]
