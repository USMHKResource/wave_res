Next
====

- Update source term loading to use correct netcdf format (from PNNL servers)
- Create a top-level `run_all.py` script that processes all of the data files.
  - calc_remote still uses 'temporary' files (but need to make them specific to each scenario)
  - Also write results files
  - Test that this script works.

## File structure

Directional data spectrum files:

    <project_root_data_folder>/<scenario>/<year>/eez/ww3.<region>.<year><month>_spec.nc

Source term data files:

    <project_root_data_folder>/<scenario>/<year>/src/ww3.<region>.<year><month><day>.1d.nc


This would compute the wave energy flux for all ranges (contours), and store it in `output/wef/wc.x100.alltime.h5` for later analysis.



Misc. Ideas (for later?)
======

## Integrating Source Terms

- Check against different integration methods.
  - Voronoi
    - Calc Voronoi and clip it
    - Calculate areas
    - Integrate
  - Other ideas?
