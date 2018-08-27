import paths as p
from netCDF4 import Dataset
import wavespeed.wave_dispersion as wavespd
import gis
import numpy as np
import pyDictH5


def _concatenate_id(array):
    out = []
    for row in array:
        out.append(''.join(row))
    return np.array(out)


def load(region, month):
    if isinstance(month, basestring):
        month = np.datetime64(month).astype("O")
    elif isinstance(month, np.datetime64):
        month = month.astype('O')
    dat = Dataset(
        p.srcdir /
        '{year:04d}/eez/ww3.{region}.{year:04d}{month:02d}_spec.nc'
        .format(year=month.year, region=region, month=month.month),
        'r')
    return dat


def calc_wef(indat):
    """Load data and calculate the average wave energy flux for a
    given region and month.

    Parameters
    ==========
    indat : a netCDF4.Dataset
       The dataset containing the data to be processed.

    Returns
    =======
    data : `pyDictH5.data` object
       A dictionary-like object containing the wave energy flux
       data. This includes:
         spec (NxFxM) [m^2/Hz/rad]: The wave energy directional spectrum.
         lon (N) (deg): longitude of data.
         lat (N) (deg): latitude of data.
         f (F) [Hz]: frequency of spectrum
         fbins (F+1) [Hz]: the bounds of each spectral bin
         direction M [deg]: The "direction-to" of wave propagation.
         Nhour () [hours]: Number of hours in average.
         wef (NxM) [W/rad/m]: directional wave energy flux
         conid (N) [3-char string]: contour identifier

    """

    v = indat.variables

    data = pyDictH5.data()
    data['spec'] = v['efth'][:, con_inds].data.mean(0)
    # For some reason lat/lon are also fn's of time, but there is no
    # additional info there.
    data['lon'] = v['longitude'][0, con_inds].data
    data['lat'] = v['latitude'][0, con_inds].data
    data['f'] = v['frequency'][:].data
    data['fbins'] = np.hstack((v['frequency1'][:].data,
                               [v['frequency2'][-1].data]))
    data['direction'] = v['direction'][:].data
    data['depth'] = v['dpt'][0, con_inds].data
    # These should all be hours, but this makes sure...
    dt = np.float32((v['time'][1] - v['time'][0]) * 24)
    data['Nhour'] = len(v['time']) * dt
    cg = wavespd.c_group0(data['f'][None, :], data['depth'][:, None])
    df = np.diff(data['fbins'])
    data['wef'] = (wavespd.rho * wavespd.gravity *
                   (data['spec'] * cg[:, :, None] * df[None, :, None]).sum(1))
    data['conid'] = _concatenate_id(v['station_name'][:, 2:5].data)
    return data


def integrate_wef(lon, lat, wef, direction):
    """Integrate the wave energy flux along a line.

    Parameters
    ==========
    lon : len(N) array (degrees)
        The longitude of the integration contour.

    lat : len(N) array (degrees)
        The latitude of the integration contour.

    wef : NxM array (W / m / rad)
        The directional wave energy flux.

    direction : len(M) array (degrees True)
        The direction the waves are propagating toward. Note: this is
        'heading direction' (0 degrees is True north)

    Returns
    =======
    tot : dict
        A dictionary of totals, containing the total for each of the
        following integration modes:
         '1way' : 'One-way mode' integrates fluxes that cross the
                  contour from left to right (if you are moving from
                  the beginning toward the end of the contour). Fluxes
                  from right to left are ignored.
         '-1way' : Is the sum of fluxes crossing the contour from left
                   to right.
         'bdir' : 'Bidirectional' mode adds fluxes from both directions.
         'trad' : 'Traditional' dot-product mode adds fluxes from one
                  direction, and subtracts them in the other.
         'unit' : 'Unit-circle' mode does not perform a
                  dot-product. This method was used in the Phase 1
                  resource assessment.

    Notes
    =====
    The right-to-left '1way' flux is tot['bdir'] - tot['1way'], or
    tot['1way'] - tot['trad'].
    """
    # The direction is the 'True heading' the waves are propagating
    # towards, so we convert that to radian polar-angles (positive
    # clockwise from east)
    theta = 90 - direction
    theta[theta < 0] += 360
    theta *= np.pi / 180
    d, midp = gis.diffll(np.stack((lon, lat)))
    norm = d * np.exp(1j * -np.pi / 2)
    dang = np.abs(np.median(np.diff(theta)))
    # Average the wave energy flux between two points, and give it a
    # complex direction
    _wef = (wef[1:] + wef[:-1]) * np.exp(1j * theta) / 2

    # Integrate.
    # Notes
    #  - `norm` has length equal to the distance between points
    #  - conjugate subtracts the normal angle
    #  - Take the `.real` component to get contour-normal fluxes
    flux = (_wef * np.conj(norm)[:, None]) * dang
    # Sum for all of the different methods.
    # Note: this is a 'double sum' (direction and line-integral)
    out = {}
    out['trad'] = flux.real.sum()
    out['bdir'] = np.abs(flux.real).sum()
    out['unit'] = np.abs(flux).sum()
    _flux = flux.real.copy()
    _flux[_flux < 0] = 0
    out['1way'] = _flux.sum()
    return out


if __name__ == '__main__':

    months = np.arange(np.datetime64('2009-01'), np.datetime64('2010-01'))
    region = 'wc'
    # This is the outer-boundary of the EEZ (not including the Canada
    # + Mexico borders)
    con_inds = slice(34, 148)

    tot = {}
    dat = {}
    for m in months:
        m_ = m.astype('O')
        tempname = (p.tmpdir / 'ww3.{region}.{year}{month:02d}_wef.nc'
                    .format(region=region, year=m_.year, month=m_.month))
        print('Processing file {}'.format(tempname.name))
        if tempname.is_file():
            dat[m] = dnow = pyDictH5.load(str(tempname))
        else:
            ncdat = load_source('wc', m)
            dat[m] = dnow = calc_wef(ncdat, 'ALL')
            dnow.to_hdf5(str(tempname))
        tot[m] = calc = integrate_wef(
            dnow['lon'][con_inds], dnow['lat'][con_inds],
            dnow['wef'][con_inds], dnow['direction'])

    final = {}
    Nh = []
    for k in tot[tot.keys()[0]]:
        final[k] = []
    for ky in tot:
        d = dat[ky]
        t = tot[ky]
        Nh.append(d['Nhour'])
        for k in t:
            final[k].append(t[k])

    Nh = np.array(Nh)
    final_sum = {}
    for k in final:
        final[k] = np.array(final[k])
        final_sum[k] = (final[k] * Nh).sum() / Nh.sum()
        print("The '{}' resource is: {: 5.1f} GW"
              .format(k, final_sum[k] / 1e9))
