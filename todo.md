IN PROGRESS
====
- Fix contour data
  - break contours where distance is greater than 40 km
    - 'Wrap' circles where distance is < 50(?) km (see `show_ak()` plot in `show_grid.py`)

Integrate Source Terms
========

- Check against different integration methods.
  - Voronoi
    - Calc Voronoi and clip it
    - Calculate areas
    - Integrate
  - Other ideas?
