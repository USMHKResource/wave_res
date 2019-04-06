The average resource along US WC for 2009, based on data from GH, is:

    The total '1way' resource is:   37.0 GW
    The total 'bidir' resource is:   39.0 GW
    The total 'trad' resource is:   34.9 GW
    The total 'unit' resource is:   55.7 GW


Computing dA
======

After spending some time digging into spherical geometry (i.e., learning about [spherical polygons](http://mathworld.wolfram.com/SphericalPolygon.html) and trying to use `scipy.spatial.SphericalVoronoi`), I did a check to compare the area of a spherical triangle to the area of the same triangle using an Albers Equal Area projection (based on the same spherical globe). The difference in area is <1% so long as you're not spanning more than 40 degrees of latitude or so (and you choose your standard parallels to span the data). In fact, the difference in area between a WGS84 ellipsoid and a sphere (both projected) is generally larger than the difference between true spherical area and a spherical projection. This makes me think that a WGS84 ellipsoid (the gold standard) projection is more accurate than doing this in spherical coordinates. This is all demonstrated [here](https://github.com/USMHKResource/wave_res/tree/test-projection-distortion/check_proj_distortion.py).

Next was the task of how *to do* the area integral itself. There are several methods, but I ended up deciding to go with a 'triangle integral' approach, that uses triangles as defined by a Delaunay triangulation, for example. A standard Delaunay triangulation, however, creates a 'convex' mesh of the grid points, which will span over land and empty (of data) areas of water that we do not actually want to integrate over. Therefore, it is necessary to 'clip' the Delaunay triangulation to only 

Notes for GGM
======

Need to clear tmpdir, to rebuild those files.
