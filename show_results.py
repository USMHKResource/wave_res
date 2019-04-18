import pyDictH5 as pdh5
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy


if __name__ == "__main__":

    regions = {'ak', 'at', 'prusvi', 'wc', 'hi'}

    # Sum over these regions to get the total.
    totregions = ['ak', 'at', 'prusvi', 'wc', 'hi']

    # Initialize dictionaries for 'baseline' totals
    remote0 = {}
    local0 = {}
    # These are time-weighted averages
    rtot0 = {}
    ltot0 = {}

    # Initialize dictionaries for 'extraction' totals
    remoteX = {}
    localX = {}
    # These are time-weighted averages
    rtotX = {}
    ltotX = {}

    unit = 'GW'
    unit = 'TWh/yr'
    _factordict = {
        'GW': 1e-9,
        'TWh/yr': 365 * 24 * 1e-12,
    }
    factor = _factordict[unit]

    irange = -1

    for ireg, region in enumerate(regions):
        rd0 = remote0[region] = pdh5.load('results/{}/{}.remote-totals.h5'
                                        .format('baseline', region))
        ld0 = local0[region] = pdh5.load('results/{}/{}.local-totals.h5'
                                        .format('baseline', region))
        rdX = remoteX[region] = pdh5.load('results/{}/{}.remote-totals.h5'
                                        .format('extraction', region))
        ldX = localX[region] = pdh5.load('results/{}/{}.local-totals.h5'
                                        .format('extraction', region))
        rd0['oneway'] = rd0['1way']
        rdX['oneway'] = rdX['1way']

        rtot0[region] = {m: (np.average(rd0[m][:, irange],
                                        weights=rd0['Nhour']) * factor)
                        for m in ['oneway', 'unit', 'bdir', 'trad']}
        ltot0[region] = {m: (np.average(ld0[m][:, :irange].sum(-1),
                                        weights=ld0['Nhour']) * factor)
                        for m in ['sbt', 'sds', 'snl', 'stot', 'sin', 'sice']}
        rtotX[region] = {m: (np.average(rdX[m][:, irange],
                                        weights=rdX['Nhour']) * factor)
                        for m in ['oneway', 'unit', 'bdir', 'trad']}
        ltotX[region] = {m: (np.average(ldX[m][:, :irange].sum(-1),
                                        weights=ldX['Nhour']) * factor)
                        for m in ['sbt', 'sds', 'snl', 'stot', 'sin', 'sice']}


    for ireg, region in enumerate(totregions):
        # Calculate the total across all regions
        if 'total' not in rtot0:
            rtot0['total'] = deepcopy(rtot0[region])
            ltot0['total'] = deepcopy(ltot0[region])
            rtotX['total'] = deepcopy(rtotX[region])
            ltotX['total'] = deepcopy(ltotX[region])
        else:
            for m in ['oneway', 'unit', 'bdir', 'trad']:
                rtot0['total'][m] += rtot0[region][m]
                rtotX['total'][m] += rtotX[region][m]
            for m in ['sbt', 'sds', 'snl', 'stot', 'sin', 'sice']:
                ltot0['total'][m] += ltot0[region][m]
                ltotX['total'][m] += ltotX[region][m]

    print("")

    print("Remote Totals ({})".format(unit))
    print("..................")
    print(("{:10s}|" + "{:>8s}" * 4).format("", 'one-way', 'trad', 'unit', 'bidir'))
    print("-" * 44)
    for ireg, region in enumerate(regions):
        rt = rtot0[region]
        print("{:10s}: {oneway: 10.4g} {trad: 10.4g} {unit: 10.4g} {bdir: 10.4g}"
            .format(region, **rt))
    print("=======================================================")
    print("{:10s}: {oneway: 10.4g} {trad: 10.4g} {unit: 10.4g} {bdir: 10.4g}"
        .format('TOTAL', **rtot0['total']))

    print("")
    print("Local Baseline Totals ({})".format(unit))
    print("...............................")
    print(("{:10s}|" + "{:>11s}" * 6).format("", 'stot', 'sin', 'sds', 'snl', 'sice', 'sbt'))
    print("-" * 44)
    for ireg, region in enumerate(regions):
        lt = ltot0[region]
        print("{:10s}: {stot: 10.4g} {sin: 10.4g} {sds: 10.4g} {snl: 10.4g} {sice: 10.4g} {sbt: 10.4g}"
            .format(region, **lt))
    print("=============================================================================")
    print("{:10s}: {stot: 10.4g} {sin: 10.4g} {sds: 10.4g} {snl: 10.4g} {sice: 10.4g} {sbt: 10.4g}"
        .format('TOTAL', **ltot0['total']))

    print("")
    print("Local Potential Totals ({})".format(unit))
    print("...............................")
    print(("{:10s}|" + "{:>11s}" * 6).format("", 'stot', 'sin', 'sds', 'snl', 'sice', 'sbt'))
    print("-" * 44)
    for ireg, region in enumerate(regions):
        lt = ltotX[region]
        print("{:10s}: {stot: 10.4g} {sin: 10.4g} {sds: 10.4g} {snl: 10.4g} {sice: 10.4g} {sbt: 10.4g}"
            .format(region, **lt))
    print("=============================================================================")
    print("{:10s}: {stot: 10.4g} {sin: 10.4g} {sds: 10.4g} {snl: 10.4g} {sice: 10.4g} {sbt: 10.4g}"
        .format('TOTAL', **ltotX['total']))

