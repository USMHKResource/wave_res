import numpy as np
import wave_res.calc_remote as cr
from wave_res import paths
import pyDictH5 as pdh5

run_these = ['ak', 'hi', 'wc', 'at', 'prusvi']
#run_these = ['ak']
#run_these = ['hi']
#run_these = ['wc']
#run_these = ['at']
#run_these = ['prusvi']

months = np.arange(np.datetime64('1979-01'), np.datetime64('2011-01'))
###months = np.arange(np.datetime64('2009-01'), np.datetime64('2009-03'))

outdir = paths.tmpdir / 'baseline_spatial'

paths.mkdir(str(outdir))

print("This job is process these regions:")
print(run_these)
print('------------- Begin Loop ----------------')

for region in run_these:

    print("Processing data for region: {}".format(region))

    _tmp = cr.load_processed('baseline', region, np.datetime64('2009-01'))
    N_ang = len(_tmp['direction'])
    N_x = len(_tmp['lon'])
    N_t = len(months)
    
    out = _tmp.copy()
    out.pop('f')
    out.pop('fbins')
    out.pop('wef')
    out.pop('cg')
    out['wef'] = np.zeros((N_x, N_t, N_ang), dtype=np.float32)
    out['Nhour'] = np.zeros(N_t, dtype=np.float32)
    out['time'] = months

    for imo, month in enumerate(months):
        print(month)
        dnow = cr.load_processed('baseline', region, month)
        out['wef'][:, imo] = (dnow['wef'] * np.diff(dnow['fbins'])[None, :, None]).sum(1)

    out.to_hdf5(str(outdir / 'ww3.{}.wef_spatial.h5'.format(region)))
