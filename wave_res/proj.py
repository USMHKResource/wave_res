from cartopy import crs
import cartopy.feature as cfeature

scale = '50m'
land = cfeature.NaturalEarthFeature('physical', 'land', scale,
                                    edgecolor='face',
                                    facecolor=cfeature.COLORS['land'])
states = cfeature.STATES.with_scale(scale)

pc = crs.PlateCarree()


proj = {}
proj['wc'] = crs.AlbersEqualArea(
    #central_longitude=-124.74, central_latitude=48.38,  # ~Cape Flattery
    central_longitude=-124.41, central_latitude=40.44,  # ~Cape Mendocino
    standard_parallels=[32, 49]
)
proj['wc'].lonlim = [-130, -116]
proj['wc'].latlim = [30, 49.5]

proj['wc.wa'] = proj['wc.or'] = proj['wc.ca'] = proj['wc']

proj['at'] = crs.AlbersEqualArea(
    central_longitude=-75.52, central_latitude=35.22,  # ~Cape Hatteras
    standard_parallels=[25, 44]
)
proj['at'].lonlim = [-98, -65.]
proj['at'].latlim = [23, 45.]

proj['ec'] = crs.AlbersEqualArea(
    central_longitude=-75.52, central_latitude=35.22,  # ~Cape Hatteras
    standard_parallels=[25, 44]
)
proj['ec'].lonlim = [-84, -65]
proj['ec'].latlim = [23, 45]

proj['ec.ne'] = proj['ec.ma'] = proj['ec.se'] = proj['ec']

proj['gm'] = crs.AlbersEqualArea(
    central_longitude=-97.15, central_latitude=25.96,  # ~Gulf end of US-Mexico border
    standard_parallels=[25, 44]
)
proj['gm'].lonlim = [-98, -78]
proj['gm'].latlim = [23, 32]

proj['prusvi'] = crs.AlbersEqualArea(
    central_longitude=-66.125, central_latitude=18.472,  # Punta Del Morro
    standard_parallels=[18, 20]
)
proj['prusvi'].lonlim = [-70, -63.]
proj['prusvi'].latlim = [14, 22.]


proj['ak'] = crs.AlbersEqualArea(
    central_longitude=-150.966, central_latitude=59.196,  # Gore Point
    standard_parallels=[56, 60]
)
proj['ak'].lonlim = [-190, -135]
proj['ak'].latlim = [48, 69]

for ky in ['ALA', 'ALB', 'BOW', 'COK', 'GEO', 'GOA', 'KOD', 'MAT', 'NAL', 'NAV', 'NOR', 'SHU']:
    proj['ak.' + ky] = proj['ak']


proj['hi'] = crs.AlbersEqualArea(
    central_longitude=-157.80, central_latitude=21.25,  # Diamond head (Oahu)
    standard_parallels=[20, 22]
)
proj['hi'].lonlim = [-182, -150.]
proj['hi'].latlim = [15, 33.]
