Next
====

Organize the repo
- Make the tools a package?
- Move 'data processing' files into a folder (a subpackage?)
- Test that data-processing files still work.
- Create top-level callable scripts, or a single `run_all.py` script (see below)?

For discussion with GGM
====

## File structure

Hopefully GGM has a consistent file structure? e.g., something like:

Directional data spectrum files:

    <project_root_data_folder>/<scenario>/<year>/eez/ww3.<region>.<year><month>_spec.nc

Source term data files:

    <project_root_data_folder>/<scenario>/<year>/src/ww3.<region>.<year><month><day>.1d.nc


## Cleanup `calc_remote.py`
- Continue to use temporary files? If so, we'll need to capture the 'case' in the output location.
- Would it be useful to create a *command-line callable script*, something like:

        calc_remote.py --scenario extract100 --region wc --start 1980-01-01 --end 2010-12-31 --output output/remote/wc.x100.alltime.h5

This would compute the wave energy flux for all ranges (contours), and store it in `output/wef/wc.x100.alltime.h5` for later analysis.

## Cleanup `calc_local.py`
- Allow for computing Delaunay at runtime when the grid/region (clipping, etc.) isn't defined in `base.py`? (i.e., for GGM's idealized sensitivity studies clipping may not be important on 'square' domains.)
- Would it be useful to create a *command-line callable script*, something like:

        calc_local.py --scenario extract100 --region wc --start 1980-01-01 --end 2010-12-31 --output output/local/wc.x100.alltime.h5

This would compute the source terms for all ranges, and store it in
`output/local/wc.x100.alltime.h5` for later analysis. Or, instead of
creating a callable script, should we just have a Python script that
runs the processes the data for all of the regions+scenarios that we
want to look at?


Misc. Ideas (for later?)
======

## Integrating Source Terms

- Check against different integration methods.
  - Voronoi
    - Calc Voronoi and clip it
    - Calculate areas
    - Integrate
  - Other ideas?
