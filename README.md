This repo is for calculating the US wave resource.

Getting started
==========

## Install dependencies

Build a python environment for this project from the `environment.yml` file included:

    conda env create environment.yml

## Set path of source data

In order for these tools to work, you must set the `srcdir` variable in the `paths.py` file to point to the root directory where the model data is located. You should also set the `tmpdir` variable to point to a location for storing processed/intermediate files.


Terms and Acronyms
=======

## Cases
These are the WW3 model runs that implement different amounts of energy extraction, and other sensitivity tests.

## Regions
These are the regions over which we are computing the wave resource. These are:
    AK - Alaska
    AT - Atlantic (which is further divided into 'EC' and 'GM')
    EC - East Coast
    GM - Gulf of Mexico
    HI - Hawaii
    PRUSVI - Puerto Rico and U.S. Virgin Islands
    WC - West Coast

## Ranges
The 'distance from shore' contours (in nautical miles, or nmi) along which we compute the remote resource, and along which the data is stored. These are (10 nmi, 20 nmi, 30 nmi, ... 200 nmi). The 200 nmi range is also known as the 'eez'. The range contours do not include the borders. The 'EEZ' contour (different from 'eez' range, i.e. case-sensitive) does include the borders.

## Borders
The edges of the U.S. EEZ that are adjacent to other national boundaries (vs. the 'law of the sea' open ocean that has no political designation).

## EEZ
Exclusive Economic Zone

## Clipping Boundaries
The area between the coastline and a range contour. These are used to determine which Delaunay triangules are within each range.

Calculating Remote Resource
=========

`calc_remote.py` is a command-line tool for calculating the remote resource for each region. Take a look at the docstring for that script for more info. The `run_calc_remote.py` is a script for using the functions in `calc_remote.py` in a Python script (rather than as a command-line tool).

These tools calculates the remote resource based on the data in the `{year}/eez/ww3.{region}.{year}{month}_spec.nc` data files. They create *temporary* files (in `paths.tmpdir`) that dramatically improve the speed of doing this calculation after the first time (for each region).

Calculating Local Resource
===========

## Data Setup

### Process `.src` files into `.nc`

Prior to calculating the local resource, the source term output files (text files ending with `.src`), must be converted to netCDF. Use the `read_source.py` script to accomplish this task. This gives a ~8x reduction in the size of these files (from 500GB to 70GB), and dramatically increases the speed of accessing that data.

### Clipping the area integration

A critical -- and tedious -- piece of computing an area integral of the source terms is to estimate the differential area of each data point (i.e., `dA`). The approach taken ends up uses a 'triangle' integration method, where the triangles are defined by a Delaunay-triangulation. The 'convex' Delaunay triangulation is then clipped to only include the triangles that are within the study area. A significant amount of work has gone into properly clipping these triangles for each region so that we can estimate area integrals for each range. This was accomplished, for each region, following these steps:

1. `proj`.py - Establish a projections that is used throughout all stages of the analysis for each region. 
2. `extract_land.py` gets the coastline data (from the projection) for each region. There is some 'hardcoding' here to get the right pieces of coastline for each region. The data from this step is stored in `data/LandData.pkl`, for quick use/access.
3. `boundaries.py` defines the 'boundaries' of each range. The data from this step is stored in `data/Boundaries.pkl`, for quick retrieval.
   - For regions that intersect the 'mainland', this basically pieces a range contour together with the correct borders, and then follows the coastline path back to the other border, and then out to the start of that range again.
   - For regions that are islands (PRUSVI, HI), this means looping the boundaries with the borders as necessary.
   - There is *a lot* of hardcoding in this script.
4. `create_triangles.py` this clips the triangles from the full Delaunay triangulation to generate the sets of triangles that are within each 'range'. The data from this step is stored in `data/DiffTriangles.pkl.gz`. 

All of these steps can be repeated by running the `setup_data.py` script, but it is not necessary because the data files store the relavent information, and are loaded/used appropriately by the other tools in the toolbox.

## Computing local resource

The `calc_local.py` script currently runs 
