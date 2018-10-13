from cartopy import crs

pc = crs.PlateCarree()

proj = {}
proj['wc'] = crs.AlbersEqualArea(
    #central_longitude=-124.74, central_latitude=48.38,  # ~Cape Flattery
    central_longitude=-124.41, central_latitude=40.44,  # ~Cape Mendocino
    standard_parallels=[32, 49]
)
proj['wc'].lonlim = [-133, -116]
proj['wc'].latlim = [31.8, 49.5]

proj['at'] = crs.AlbersEqualArea(
    central_longitude=-75.52, central_latitude=35.22,  # ~Cape Hatteras
    standard_parallels=[25, 44]
)
proj['at'].lonlim = [-98, -65.]
proj['at'].latlim = [23, 45.]
