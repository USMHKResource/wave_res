# Clipping the area integration

A critical -- and tedious -- piece of computing an area integral of the source terms is to estimate the differential area of each data point (i.e., `dA`). The approach taken ends up uses a 'triangle' integration method, where the triangles are defined by a Delaunay-triangulation. The 'convex' Delaunay triangulation is then clipped to only include the triangles that are within the study area. A significant amount of work has gone into properly clipping these triangles for each region so that we can estimate area integrals for each range. This was accomplished, for each region, following these steps:

1. `proj`.py - Establish a projections that is used throughout all stages of the analysis for each region. 
2. `extract_land.py` gets the coastline data (from the projection) for each region. There is some 'hardcoding' here to get the right pieces of coastline for each region. The data from this step is stored in `data/LandData.pkl`, for quick use/access.
3. `boundaries.py` defines the 'boundaries' of each range. The data from this step is stored in `data/Boundaries.pkl`, for quick retrieval.
   - For regions that intersect the 'mainland', this basically pieces a range contour together with the correct borders, and then follows the coastline path back to the other border, and then out to the start of that range again.
   - For regions that are islands (PRUSVI, HI), this means looping the boundaries with the borders as necessary.
   - There is *a lot* of hardcoding in this script.
4. `create_triangles.py` this clips the triangles from the full Delaunay triangulation to generate the sets of triangles that are within each 'range'. The data from this step is stored in `data/DiffTriangles.pkl.gz`. 

All of these steps can be repeated by running the `setup_data.py` script, but it is not necessary because the data files store the relavent information, and are loaded/used appropriately by the other tools in the toolbox.

