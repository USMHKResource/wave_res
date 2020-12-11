import pyDictH5 as pdh5
import wave_res as wr
import wave_res.paths as p
import wave_res.calc_remote as cr
import wave_res.gis as gis
import numpy as np
import pandas as pd
from base import _factor

threshold = 0

datdir = p.tmpdir / 'baseline_spatial'

units = 'TWh/yr'

def calc_tech(region, contour, trange, threshold=0):

    rinf = wr.RegionInfo(region)
    s_region = rinf.source_region
    data = pdh5.load(str(datdir / 'ww3.{}.wef_spatial.h5'.format(s_region)))
    tinds = (trange[0] <= data['time']) & (data['time'] < trange[1])

    rcon = '{:03d}'.format(contour)
    out = 0
    
    for cid in rinf.con_defs[rcon]:
        lon, lat = data['lon'][cid], data['lat'][cid]
        norm, midp = gis.diffll(np.stack((lon, lat)))
        L = np.abs(norm)
        trad, oneway, bdir, unit = cr.integrate_wef(lon, lat, data['wef'][cid][:, tinds].mean(1), data['direction'], sum_axes=(-1, ))
        if region == 'hi' and rcon == '200':
            oneway = bdir - oneway
        use_this = bdir
        #use_this = oneway
        use_this -= threshold * L
        use_this[use_this < 0] = 0
        out += use_this.sum()
    return out

def checkval(region, contour, trange):
    print("Checking region: {}; contour: {}; time-range: {}".format(region, contour, trange))

    checkdat = pdh5.load('results/freq.fcut/baseline/{}.remote-totals.h5'.format(region))

    tinds = (trange[0] <= checkdat['time']) & (checkdat['time'] < trange[1])
    con_ind = np.nonzero(contour == g['range'])[0][0]

    checkval = (checkdat['1way'][tinds][..., con_ind] * np.diff(checkdat['fbins'])[None, :]).sum(1).mean(0)

    techval = calc_tech(region, contour, trange, 0)
    
    if (checkval - techval) / checkval < 1e-3:
        print("... Good!")
    else:
        print("!!! FAIL !!!")

if __name__ == '__main__':


    trange = [np.datetime64('1979-01'), np.datetime64('2011-01')]
    cases = [(10, 0), (10, 1000), (10, 5000), (10, 8000),
             (200, 0), (200, 1000), (200, 5000), (200, 8000)]

    regions = ['wc', 'wc.wa', 'wc.or', 'wc.ca',
               'ec', 'ec.ne', 'ec.ma', 'ec.se',
               'hi', 'ak', 'gm', 'prusvi']
    results = pd.DataFrame(index=regions, columns=cases)
    
    for region in regions:
        print("Running region: {}".format(region))
        for con, threshold in cases:
            print("...Running case: {}".format((con, threshold)))
            results.loc[region, (con, threshold)] = calc_tech(region, con, trange, threshold) * _factor
        
    results.to_csv('results/Technical_Results-bidir.csv')
