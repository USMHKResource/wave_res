import pyDictH5 as pdh5
import numpy as np
import matplotlib.pyplot as plt

regions = {'ak', 'at', 'prusvi', 'wc', 'hi'}
scenario = 'baseline'
scenario = 'extraction'

remote = {}
local = {}
rtot = {}
ltot = {}

unit = 'GW'
unit = 'TWh/yr'
_factordict = {
    'GW': 1e9,
    'TWh/yr': 365 * 24 * 1e-12,
}
factor = _factordict[unit]

for ireg, region in enumerate(regions):
    rd = remote[region] = pdh5.load('results/{}/{}.remote-totals.h5'.format(scenario, region))
    ld = local[region] = pdh5.load('results/{}/{}.local-totals.h5'.format(scenario, region))

    rd['oneway'] = rd['1way']

    rtot[region] = {m: (np.average(rd[m][:, -1], weights=rd['Nhour']) * factor)
                    for m in ['oneway', 'unit', 'bdir', 'trad']}
    ltot[region] = {m: (np.average(ld[m][:, 2:-2].sum(-1), weights=ld['Nhour']) * factor)
                    for m in ['sbt', 'sds', 'snl', 'stot', 'sin', 'sice']}

print("")

print("Remote Totals ({})".format(unit))
print("============")
print(("{:10s}|" + "{:>8s}" * 4).format("", 'one-way', 'trad', 'unit', 'bidir'))
print("-" * 44)
for ireg, region in enumerate(regions):
    rt = rtot[region]
    print("{:10s}: {oneway: 7.4g} {trad: 7.4g} {unit: 7.4g} {bdir: 7.4g}"
          .format(region, **rt))


print("Local Totals ({})".format(unit))
print("============")
print(("{:10s}|" + "{:>8s}" * 6).format("", 'stot', 'sin', 'sds', 'snl', 'sice', 'sbt'))
print("-" * 44)
for ireg, region in enumerate(regions):
    lt = ltot[region]
    print("{:10s}: {stot: 7.4g} {sin: 7.4g} {sds: 7.4g} {snl: 7.4g} {sice: 7.4g} {sbt: 7.4g}"
          .format(region, **lt))
