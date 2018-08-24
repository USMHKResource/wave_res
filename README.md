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

The `calc_remote.py` file calculates the remote resource from the `{year}/eez/ww3.{region}.{year}{month}_spec.nc` data files.
