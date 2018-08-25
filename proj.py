from cartopy import crs

pc = crs.PlateCarree()

proj = {}
proj['wc'] = crs.AlbersEqualArea(
    #central_longitude=-124.74, central_latitude=48.38,  # ~Cape Flattery
    central_longitude= -124.41, central_latitude=40.44,  # ~Cape Mendocino
    standard_parallels=[32, 49]
)
