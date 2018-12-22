This repo is for calculating the US wave resource.

Getting started
==========

## Install dependencies

Build a python environment for this project from the `environment.yml` file included:

    conda env create environment.yml

## Set path of source data

In order for these tools to work, you must set the `srcdir` variable in the `paths.py` file to point to the root directory where the model data is located. You should also set the `tmpdir` variable to point to a location for storing processed/intermediate files (or just use the local one provided).

## Run

Execute the `run_all.py` script to calculate results for the scenarios and regions defined in the `run_these` variable. Output from this script are written to the `results/` folder as `hdf5` files, which can be loaded using `pyDictH5.load`.


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
