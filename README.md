This repo is for calculating the US wave resource.

Getting started
==========

Build a python environment for this project from the `environment.yml` file included:

    conda env create environment.yml

In order for most of these scripts to work, you will also need to set the `srcdir` variable in the `paths.py` file to point to the root directory where the model data is located. You should also set the `tmpdir` variable to point to a location for storing processed/intermediate files.

Data Processing
========

The `read_source.py` file reads 1-D source-term files (text files ending with `.src`), and writes the data to netCDF. This gives a ~8x reduction in the size of these files (from 500GB to 70GB), and dramatically increases the speed of accessing that data.

Calculating Remote Resource
=========

`calc_remote.py` is a command-line tool for calculating the remote resource for each region. Take a look at the docstring for that script for more info. The `run_calc_remote.py` is a script for using the functions in `calc_remote.py` in a Python script (rather than as a command-line tool).

These tools calculates the remote resource based on the data in the `{year}/eez/ww3.{region}.{year}{month}_spec.nc` data files. They create *temporary* files (in `paths.tmpdir`) that dramatically improve the speed of doing this calculation after the first time (for each region).

Calculating Local Resource
===========

After spending some time digging into spherical geometry (i.e., learning about [spherical polygons](http://mathworld.wolfram.com/SphericalPolygon.html) and trying to use `scipy.spatial.SphericalVoronoi`), I did a check to compare the area of a spherical triangle to the area of the same triangle using an Albers Equal Area projection (based on the same spherical globe). The difference in area is <1% so long as you're not spanning more than 40 degrees of latitude or so (and you choose your standard parallels to span the data). In fact, the difference in area between a WGS84 ellipsoid and a sphere (both projected) is generally larger than the difference between true spherical area and a spherical projection. This makes me think that a WGS84 ellipsoid (the gold standard) projection is more accurate than doing this in spherical coordinates. This is all demonstrated [here](https://github.com/USMHKResource/wave_res/tree/test-projection-distortion/check_proj_distortion.py).
