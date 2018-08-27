from netCDF4 import Dataset
import numpy as np
import paths as p
import proj
import area_integral as aint
from scipy.spatial import Delaunay
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


def load_source(region, dt):
    if isinstance(dt, basestring):
        dt = np.datetime64(dt).astype("O")
    elif isinstance(dt, np.datetime64):
        dt = dt.astype('O')
    dat = Dataset(
        p.srcdir /
        '{year}/src/ww3.{region}.{year}{month:02d}{day:02d}.1d.nc'
        .format(year=dt.year, month=dt.month, day=dt.day, region=region))
    return dat


def integrate_source(region, dates,
                     terms=['Sin', 'Snl', 'Sds', 'Sbt',
                            'Sice', 'Stot']):

    prj = proj.proj[region]
    df = np.diff(base.freqbins[region])
    # Load the first dataset to compute the grid
    dat = load_source(region, dates[0])
    xy = prj.transform_points(
        proj.pc,
        dat.variables['lon'][:],
        dat.variables['lat'][:])[:, :2]
    verts = Delaunay(xy).vertices
    dat.close()

    # Initialize the output
    out = pyDictH5.data()
    out['time'] = dates
    for ky in terms:
        out[ky] = np.zeros(len(dates),
                           dtype=np.float32)

    # Load data and perform integral
    for idt, dt in enumerate(dates):
        print("Integrating {}...".format(dt))
        dat = load_source(region, dt)
        for ky in terms:
            # average time + integrate frequency
            src = (dat.variables[ky][:].mean(0) *
                   df[None, :]).sum(-1)
            # Now integrate area
            zsum, _ = aint.sumtriangles(xy, src, verts)
            out[ky][idt] = zsum
        dat.close()

    # Return units of Watts
    for ky in terms:
        out[ky] *= rho * gravity
    return out


if __name__ == '__main__':

    region = 'wc'
    terms = ['Sin', 'Snl', 'Sds', 'Sbt', 'Sice', 'Stot']

    # A set of test dates.
    dates = np.arange(np.datetime64('2009-01-01'),
                      np.datetime64('2009-01-10'))
    tag = '2009test'

    # The full 2009 dataset.
    dates = np.arange(np.datetime64('2009-01-01'),
                      np.datetime64('2010-01-01'))
    tag = '2009'

    Psrc = integrate_source(region, dates, terms=terms)

    Psrc.to_hdf5(str(p.tmpdir /
                     'integrated_source.{region}.{tag}.h5'
                     .format(region=region, tag=tag)))

    print("The area-integrated source terms for "
          "the '{region}' are:"
          .format(region=region))
    print(" averaged from {} to {}".format(dates[0], dates[-1]))
    print("------------------------")
    for ky in terms:
        print("   {} : {: 8.2f} GW".format(ky, Psrc[ky].mean() / 1e9))
