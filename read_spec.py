import paths as p
from netCDF4 import Dataset
import wavespeed.wave_dispersion as wavespd
import gis
import numpy as np
import matplotlib.pyplot as plt
import pyDictH5
plt.ion()

months = np.arange(np.datetime64('2009-01'), np.datetime64('2010-01'))

region = 'wc'
month = '2009-04'
contour = 'ALL'

def _concatenate_id(array):
    out = []
    for row in array:
        out.append(''.join(row))
    return np.array(out)


def calc_wef(region, month, contour='ALL'):
    """Calculate the average wave energy flux for a given zone and month.

    Parameters
    ==========
    region : string {ak, at, hi, prusvi, wc}

    month : numpy.datetime64, string (e.g., '2009-05'), or datetime object

    contour : string or integer {'ALL' (default), 'EEZ', 10, 20, 30, ... 190}

    Returns
    =======
    data : dict
           A dictionary containing the wave energy flux data.
    """

    if not isinstance(contour, basestring):
        # This is to handle integers
        contour = '{:03d}'.format(contour)
    contour = contour.upper()
    if isinstance(month, basestring):
        month = np.datetime64(month).astype("O")
    elif isinstance(month, np.datetime64):
        month = month.astype('O')
    indat = Dataset(
        p.srcdir /
        '{year:04d}/eez/ww3.{region}.{year:04d}{month:02d}_spec.nc'
        .format(year=month.year, region=region, month=month.month),
        'r')
    v = indat.variables
    conid = _concatenate_id(v['station_name'][:, 2:5].data)
    if contour == 'ALL':
        con_inds = slice(None)
    else:
        tmp = np.nonzero(conid == contour)[0][[0, -1]]
        con_inds = slice(tmp[0], tmp[1] + 1)

    data = pyDictH5.data()
    data['spec'] = v['efth'][:, con_inds].data.mean(0)

    # For some reason lat/lon are also fn's of time, but there is no
    # additional info there.
    data['lon'] = v['longitude'][0, con_inds].data
    data['lat'] = v['latitude'][0, con_inds].data
    data['f'] = v['frequency'][:].data
    data['fbins'] = np.hstack((v['frequency1'][:].data, [v['frequency2'][-1].data]))
    data['direction'] = v['direction'][:].data
    data['depth'] = v['dpt'][0, con_inds].data
    # These should all be hours, but this makes sure...
    dt = np.float32((v['time'][1] - v['time'][0]) * 24)
    data['Nhour'] = len(v['time']) * dt
    cg = wavespd.c_group0(data['f'][None, :], data['depth'][:, None])
    df = np.diff(data['fbins'])
    data['wef'] = (wavespd.rho * wavespd.gravity *
                   (data['spec'] * cg[:, :, None] * df[None, :, None]).sum(1))
    data['conid'] = conid
    return data

def integrate_wef(lonlat, wef, direction):
    # The direction is the 'True heading' the waves are propagating
    # towards, so we convert that to radian polar-angles (positive
    # clockwise from east)
    theta = 90 - direction
    theta[theta < 0] += 360
    theta *= np.pi / 180
    d, midp = gis.diffll(lonlat)
    norm = d * np.exp(1j * -np.pi / 2)
    dang = np.abs(np.median(np.diff(theta)))
    _wef = (wef[1:] + wef[:-1]) * np.exp(1j * theta) / 2
    flux = (_wef * np.conj(norm)[:, None]).real * dang
    out = {}
    out['trad'] = flux.sum()
    out['bdir'] = np.abs(flux).sum()
    flux[flux < 0] = 0
    out['1way'] = flux.sum()
    return out

if __name__ == '__main__':

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
            dat[m] = dnow = calc_wef('wc', m, 'ALL')
            dnow.to_hdf5(str(tempname))
        tot[m] = calc = integrate_wef(np.stack((dnow['lon'], dnow['lat'])), dnow['wef'], dnow['direction'])


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
        print("The '{}' resource is: {: 5.1f} GW".format(k, final_sum[k] / 1e9))
