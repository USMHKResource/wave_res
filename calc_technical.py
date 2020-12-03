import pyDictH5 as pdh5
from copy import deepcopy
import numpy as np
import wave_res as wr
import wave_res.calc_remote as cr
from netCDF4 import Dataset
import wave_res.base as wrb



tmp_ = cr.load_processed('baseline', 'wc', np.datetime64('2009-01'))

wef = (tmp_['wef']*np.diff(tmp_['fbins'])[None,:,None]).sum(1)


def calc_normal_fluxes(scenario, region, months, ranges=[10], threshold_W_m=0):
    rinf = base.RegionInfo(region)
    out = ConNormFluxGroup()
    out['time'] = months
    out['Nhour'] = np.empty((len(months)), dtype=np.uint16)

    for rky in rinf.con_defs.keys():
        out[rky] = ContourNormalFluxes()
        N_x = 0
        for ci in rinf.con_defs[rky]:
            N_x += len(ci)
        for int_ky in cr.wef_int_modes:
            out[rky][int_ky] = np.zeros(len(months), N_x)

    for imo, mo in enumerate(months):
        dnow = load_processed(scenario, rinf.source_region, mo)
        # Integrate in frequency (direction persists)
        wef = (dnow['wef'] * np.diff(rinf.freqbins)[None, :, None]).sum(1)
        out['Nhour'][imo] = dnow['Nhour']
        for rky in rinf.con_defs:
            dtmp = out[rky]
            con_inds = rinf.con_defs[rky]
            for ci in con_inds:
                (dtmp['trad'], dnow['1way'],
                 dnow['bidir'], dnow['unit']) += cr.integrate_wef(
                     dnow['lon'][ci], dnow['lat'][ci],
                     dnow['wef'][ci], dnow['direction'],
                     sum_axes=(-1, ) )
                con_lengt
                

regions = {'ak', 'ec', 'wc', 'hi', 'gm'}
remote_terms = ['trad', '1way', 'bdir', 'unit']

for ireg, region in enumerate(regions):

    tmp = pdh5.load('results/freq.fcut/{}/{}.remote-totals.h5'.format('baseline', region))

    dat = int_freq(tmp, remote_terms)

    # [0] for the 10 nm data
    flux = dat['1way'][:, 0] / dat['length'][0]
    break
