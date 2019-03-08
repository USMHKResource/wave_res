from netCDF4 import Dataset
import numpy as np
import paths as p
import area_integral as aint
from wavespeed.wave_dispersion import rho, gravity
import pyDictH5
import base

# Wave Flux Units
# ########
#   s: wave-height var spec (m2 s)
# J = rho g cg s df = (kg/m3) (m/s2) (m/s) (m2 s) (1/s) = kg m /s3
# The line integral of this gives power
# P = \int J dl = kg m2 /s3

# Wave Source Term Units
# ########
#   S: source spec (m2)
# Q = rho g S df = (kg/m3) (m/s2) (m2) (1/s) = kg /s3
# The area integral of this gives power
# P = \int Q dA

source_terms = ['sin', 'snl', 'sds', 'sbt', 'sice', 'stot']

class LocalTotals(pyDictH5.data):
    pass


class LocalResults(pyDictH5.data):

    def hourly_average(self):
        out = LocalTotals()
        out['area'] = self['area']
        out['range'] = self['range']
        for ky in source_terms:
            out[ky] = np.average(self[ky],
                                 weights=self['Nhour'],
                                 axis=0)
        return out


def load_source(scenario, region, dt):
    if isinstance(dt, basestring):
        dt = np.datetime64(dt).astype("O")
    elif isinstance(dt, np.datetime64):
        dt = dt.astype('O')
    dat = Dataset(
        p.srcdir +
        '/{scenario}/usa/{year}/src/ww3.{region}.{year}{month:02d}_tab.nc'
        .format(scenario=scenario, region=region,
                year=dt.year, month=dt.month))
    return dat


def calc_local(scenario, region, dates,
               terms=source_terms):
    """Calculate the local resource for `scenario` in `region` for
    `months`.
    
    Parameters
    ==========
    scenario : string
         the 'scenario' (case) to evaluate.
    region : string {'wc', 'ec', 'at', 'gm', 'prusvi', 'hi'}
         the region
    months : list/array of np.datetime64.
         The months to calculate.
    source_terms : list of strings {'sin', '}

    Returns
    =======
    results : a dictionary-like object (based on pyDictH5.data)
       The results object contains the following keys:
         'time': (n_months) the month
         'Nhour': (n_months) The number of hours in each month
         'range': (n_ranges) The range of each contour [nautical miles]
         'area': (n_ranges) The area between each contour [meters^2]
         'sin': (n_months, n_ranges) Monthly averaged 'wind input'
                source term. [watts]
         'sds': (n_months, n_ranges) Monthly averaged 'dissipation'
                source term. [watts]
         'snl': (n_months, n_ranges) Monthly averaged 'non-linear'
                source term. [watts]
         'sbt': (n_months, n_ranges) Monthly averaged 'bottom
                friction' source term. [watts]
         'sbt': (n_months, n_ranges) Monthly averaged 'ice'
                source term. [watts]
         'stot': (n_months, n_ranges) Monthly averaged 'total'
                source term. [watts]
    Notes
    =====
    The source terms are totals for the area between that contour
    range, and the range inshore of it. So, to get the source term for
    the entire EEZ, you need to sum them.  For example:

          sin_total = local['sin'].sum(-1)

    """

    rinf = base.RegionInfo(region)
    xy = rinf.allxy.T
    df = np.diff(rinf.freqbins)
    src = np.zeros(xy.shape[0], dtype=np.float32)
    n_grid = rinf.gridxy.shape[1]

    # # Load the first dataset to compute the grid
    # dat = load_source(region, dates[0])
    # xy = prj.transform_points(
    #     proj.pc,
    #     dat.variables['lon'][:],
    #     dat.variables['lat'][:])[:, :2]
    # verts = Delaunay(xy).vertices
    # dat.close()

    # Initialize the output
    out = LocalResults()
    out['time'] = dates
    out['Nhour'] = np.zeros(len(out['time']), dtype=np.uint16)
    for ky in terms:
        out[ky] = np.zeros((len(dates), 20),
                           dtype=np.float32)

    out['range'] = np.arange(10, 201, 10)

    out['area'] = np.zeros(len(out['range']), dtype=np.float32)
    for irng, rng in enumerate(range(10, 201, 10)):
        rky = '{:03d}'.format(rng)
        verts = rinf.tri_inds[rky]
        # Now integrate
        _, area = aint.sumtriangles(xy, np.zeros_like(xy), verts)
        out['area'][irng] = area

    # Load data and perform integral
    for idt, dt in enumerate(dates):
        print("      Integrating {}...".format(dt))
        dat = load_source(scenario, region, dt)
        out['Nhour'][idt] = len(dat.variables['time'])
        for ky in terms:
            # average time + integrate frequency
            src[:n_grid] = (dat.variables[ky][:].mean(0) *
                            df[None, :]).sum(-1)
            for irng, rng in enumerate(range(10, 201, 10)):
                rky = '{:03d}'.format(rng)
                verts = rinf.tri_inds[rky]
                # Now integrate
                zsum, _ = aint.sumtriangles(xy, src, verts)
                out[ky][idt, irng] = zsum
        dat.close()
    # Return units of Watts
    for ky in terms:
        out[ky] *= rho * gravity
    return out