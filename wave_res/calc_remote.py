import paths as p
from netCDF4 import Dataset
import wavespeed.wave_dispersion as wavespd
import gis
import numpy as np
import pyDictH5
import warnings
import base
warnings.filterwarnings("ignore", category=RuntimeWarning)

wef_int_modes = ['trad', '1way', 'bdir', 'unit']


class RemoteTotals(pyDictH5.data):
    pass


class RemoteResults(pyDictH5.data):

    def calc_totals(self):
        """Compute hourly average, and integrate in frequency.
        """
        out = RemoteTotals()
        out['range'] = self['range']
        out['length'] = self['length']
        df = np.diff(self['fbins'])
        for ky in wef_int_modes:
            out[ky] = (np.average(self[ky],
                                  weights=self['Nhour'],
                                  axis=0) * df[:, None]).sum(0)
        return out


def _concatenate_id(array):
    out = []
    for row in array:
        out.append(''.join(row))
    return np.array(out)


def load(scenario, region, month):
    if isinstance(month, basestring):
        month = np.datetime64(month).astype("O")
    elif isinstance(month, np.datetime64):
        month = month.astype('O')
    dat = Dataset(
        p.srcdir /
        '{scenario}/usa/{year:04d}/eez/ww3.{region}.{year:04d}{month:02d}_spec.nc'
        .format(scenario=scenario,
                year=month.year, region=region, month=month.month),
        'r')
    return dat


def load_processed(scenario, region, month, overwrite=False):
    """Process large data files, compute wave energy flux, and store
    small temporary data files."""
    # create directory <tmpdir>/{scenario}/ (if it doesn't already exist)
    p.mkdir(str(p.tmpdir / '{}'.format(scenario)))
    m_ = month.astype('O')
    tempname = (p.tmpdir / '{scenario}/ww3.{region}.{year}{month:02d}_wef.nc'
                .format(scenario=scenario, region=region,
                        year=m_.year, month=m_.month))
    if tempname.is_file() and not overwrite:
        dat = pyDictH5.load(str(tempname))
    else:
        print('Processing file {}'.format(tempname.name))
        ncdat = load(scenario, region, month)
        dat = calc_wef(ncdat)
        dat.to_hdf5(str(tempname))
    return dat


def calc_wef(indat):
    """Reduce the source-data to the variables of interest.

    Parameters
    ==========
    indat : a netCDF4.Dataset
       The dataset containing the data to be processed.

    Returns
    =======
    data : `pyDictH5.data` object
       A dictionary-like object containing the wave energy flux
       data. This includes:
         lon (N) (deg): longitude of data
         lat (N) (deg): latitude of data
         f (F) [Hz]: frequency of spectrum
         fbins (F+1) [Hz]: the bounds of each spectral bin
         direction M [deg]: The "direction-to" of wave propagation
         Nhour () [hours]: Number of hours in average
         depth (N) [m]: the water depth
         cg (NxF) [m/s]: the group velocity
         wef (NxFxM) [W/rad/m/Hz]: directional wave energy flux vs. freq
         conid (N) [3-char string]: contour identifier

    """

    v = indat.variables

    data = pyDictH5.data()
    # efth has dims: (hourly time, station/location, freq, direction)
    # So, this averages away time hourly-time
    _spec = v['efth'][:].data.mean(0)
    # For some reason lat/lon are also fn's of time, but there is no
    # additional info there.
    data['lon'] = v['longitude'][0].data
    data['lat'] = v['latitude'][0].data
    data['f'] = v['frequency'][:].data
    data['fbins'] = np.hstack((v['frequency1'][:].data,
                               [v['frequency2'][-1].data]))
    data['direction'] = v['direction'][:].data
    data['depth'] = v['dpt'][0].data
    # These should all be hours, but this makes sure...
    dt = np.float32((v['time'][1] - v['time'][0]) * 24)
    data['Nhour'] = len(v['time']) * dt
    data['cg'] = wavespd.c_group0(data['f'][None, :], data['depth'][:, None])
    data['wef'] = (wavespd.rho * wavespd.gravity *
                   (_spec * data['cg'][:, :, None]))
    data['conid'] = _concatenate_id(v['station_name'][:, 2:5].data)
    return data


def integrate_wef(lon, lat, wef, direction, sum_axes=(0, -1), int_order=1):
    """Integrate the wave energy flux along a line.

    Parameters
    ==========
    lon : len(N) array (degrees)
        The longitude of the integration contour.

    lat : len(N) array (degrees)
        The latitude of the integration contour.

    wef : NxFxM array (W / m / rad / Hz)
        The directional wave energy flux.

    direction : len(M) array (degrees True)
        The direction the waves are propagating toward. Note: this is
        'heading direction' (0 degrees is True north)

    sum_axes : tuple of axes to sum-along
        The default is (0, -1), which for most data is probably
        position and direction.

    int_order : the integration order (1 or 2)
        First (1) order integration is done at the midpoints between
        the gridpoints. Second order integration is done at the
        gridpoints.

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
    # _wef has shape: N_points, N_freq, N_direction (AB: CHECKTHIS)
    if int_order == 1:
        _norm = norm
        _wef = (wef[1:] + wef[:-1]) * np.exp(1j * theta) / 2
    elif int_order == 2:
        _norm = np.zeros_like(wef)
        _norm[0], _norm[-1] = norm[0], norm[-1]
        _norm[1:-1] = (norm[1:] + norm[:-1]) / 2
        _wef = wef
    else:
        raise ValueError("Invalid value for int_order: this must be 1 or 2")

    # Integrate.
    # Notes
    #  - `norm` has length equal to the distance between points
    #  - conjugate subtracts the normal angle
    #  - Take the `.real` component to get contour-normal fluxes
    flux = (_wef * np.conj(_norm)[:, None, None]) * dang
    # Sum for all of the different methods.
    # Note: this is a 'double sum' (line-integral, direction)
    #       but we don't sum in the frequency direction.
    trad = flux.real.sum(sum_axes)
    bdir = np.abs(flux.real).sum(sum_axes)
    unit = np.abs(flux).sum(sum_axes)
    _flux = flux.real.copy()
    _flux[_flux < 0] = 0
    oneway = _flux.sum(sum_axes)
    return trad, oneway, bdir, unit

    # offshore_flux = bdir - oneway

def con_length(lon, lat):
    d, midp = gis.diffll(np.stack((lon, lat)))
    return np.abs(d).sum()


def calc_remote(scenario, region, months):
    """Calculate the remote resource for `scenario` in `region` for
    `months`.
    
    Parameters
    ==========
    scenario : string
         the 'scenario' (case) to evaluate.
    region : string {'wc', 'ec', 'at', 'gm', 'prusvi', 'hi'}
         the region
    months : list/array of np.datetime64.
         The months to calculate.

    Returns
    =======
    This returns a dictionary-like object (based on pyDictH5.data)
    containing:
        'time': (n_months) the month
        'Nhour': (n_months) The number of hours in each month
        'range': (n_ranges) The range of each contour [nautical miles]
        'length': (n_ranges) The length of each contour [meters]
        '1way': (n_months, n_ranges) Monthly averaged wave energy
                flux using the '1way' method [watts]
        'trad': (n_months, n_ranges) Monthly averaged wave energy
                flux using a traditional dot-product [watts]
        'bdir': (n_months, n_ranges) Monthly averaged wave energy
                flux using a bi-directional dot-product [watts]
        'unit': (n_months, n_ranges) Monthly averaged wave energy
                flux using the unit-circle method [watts]
    """
    rinf = base.RegionInfo(region)

    ranges = np.arange(10, 201, 10)

    N_f = len(rinf.freqbins) - 1

    # initialize the output
    out = RemoteResults()
    for ky in wef_int_modes:
        out[ky] = np.empty((len(months), N_f, len(ranges)), dtype=np.float32)
    out['time'] = months
    out['range'] = ranges
    out['Nhour'] = np.empty((len(months)), dtype=np.uint16)
    out['length'] = np.zeros(ranges.shape, dtype=np.float32)
    out['fbins'] = rinf.freqbins

    for imo, mo in enumerate(months):
        # This creates and loads temporary files as needed.
        dnow = load_processed(scenario, rinf.source_region, mo)
        out['Nhour'][imo] = dnow['Nhour']
        for irng, rng in enumerate(ranges):
            tmp = np.zeros((4, N_f), dtype=np.float32)
            rky = '{:03d}'.format(rng)
            con_inds = rinf.con_defs[rky]
            for ci in con_inds:
                #plt.plot(dnow['lon'][ci], dnow['lat'][ci])
                tmp += integrate_wef(
                    dnow['lon'][ci], dnow['lat'][ci],
                    dnow['wef'][ci], dnow['direction'])
                if imo == 0:
                    out['length'][irng] += con_length(
                        dnow['lon'][ci], dnow['lat'][ci])
            for iky, ky in enumerate(wef_int_modes):
                out[ky][imo, :, irng] = tmp[iky]
    return out
