import pyDictH5 as pdh5
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
import json
import matplotlib.pyplot as plt

regions = {'ak', 'at', 'prusvi', 'wc', 'hi'}

# Sum over these regions to get the total.
totregions = ['ak', 'at', 'prusvi', 'wc', 'hi']

# Initialize dictionaries for 'baseline' totals
remote0 = {}
local0 = {}
# These are time-weighted averages
rtot0 = {}
ltot0 = {}

with open('EPRI_totals.json') as fl:
    epri = json.load(fl)

source_terms = ['sbt', 'sds', 'snl', 'stot', 'sin', 'sice']
remote_terms = ['trad', 'oneway', 'bdir', 'unit']
    
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

if unit != 'TWh/yr':
    # The EPRI data is in TWh/yr, so we adjust units here
    for ky in epri:
        epri[ky] *= factor / _factordict['TWh/yr']        

irange = -1


def zero_pad(arr, n):
    return np.pad(arr, n, mode='constant')


def int_freq(dat, terms):
    for ky in terms:
        dat[ky] = (dat[ky] * np.diff(dat['fbins'])[None, :, None]).sum(1)


def int_freq2(dat, terms):
    df = np.diff(dat['fbins'])[None, :, None]
    f = dat['fbins'][:-1][None, :, None] + df / 2
    for ky in terms:
        dat[ky] = (dat[ky] * df).sum(1)


for ireg, region in enumerate(regions):
    rd0 = remote0[region] = pdh5.load('frequencyResults/{}/{}.remote-totals.h5'
                                      .format('baseline', region))
    ld0 = local0[region] = pdh5.load('frequencyResults/{}/{}.local-totals.h5'
                                     .format('baseline', region))
    rdX = remoteX[region] = pdh5.load('frequencyResults/{}/{}.remote-totals.h5'
                                      .format('extraction', region))
    ldX = localX[region] = pdh5.load('frequencyResults/{}/{}.local-totals.h5'
                                     .format('extraction', region))
    rd0['oneway'] = rd0['1way']
    rdX['oneway'] = rdX['1way']

    int_freq(rd0, remote_terms)
    int_freq(rdX, remote_terms)
    int_freq2(ld0, source_terms)
    int_freq2(ldX, source_terms)
    
    
    rtot0[region] = {m: (np.average(rd0[m][:, irange],
                                    weights=rd0['Nhour']) * factor)
                     for m in remote_terms}
    ltot0[region] = {m: (np.average(ld0[m][:, :irange].sum(-1),
                                    weights=ld0['Nhour']) * factor)
                     for m in source_terms}
    rtotX[region] = {m: (np.average(rdX[m][:, irange],
                                    weights=rdX['Nhour']) * factor)
                     for m in remote_terms}
    ltotX[region] = {m: (np.average(ldX[m][:, :irange].sum(-1),
                                    weights=ldX['Nhour']) * factor)
                     for m in source_terms}


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


if False:
    region = 'wc'
    rem = remote0[region]
    remX = remoteX[region]
    lc0 = local0[region]
    lcX = localX[region]
    fig = plt.figure(10);fig.clf()
    ax = fig.subplots(1, 1)
    dist = rem['range']
    dat = np.average(rem['oneway'], weights=rem['Nhour'], axis=0) * factor
    ax.plot(dist, dat, 'b-')
    dat = np.average(remX['oneway'], weights=rem['Nhour'], axis=0) * factor
    ax.plot(dist, dat, 'b-.')
    dat = np.average(remX['bdir'] - remX['oneway'], weights=rem['Nhour'], axis=0) * factor
    ax.plot(dist, dat, 'b:')
    #ax.plot(dist[:-1], np.diff(dat))
    dat = np.average(lc0['stot'], weights=rem['Nhour'], axis=0) * factor
    ax.plot(dist, np.cumsum(dat))
    dat = np.average(lcX['stot'], weights=rem['Nhour'], axis=0) * factor
    ax.plot(dist, np.cumsum(dat))

    fig = plt.figure(11);fig.clf()
    ax = fig.subplots(1, 1)
    dist = rem['range']
    dat_one = np.average(rem['oneway'], weights=rem['Nhour'], axis=0) * factor
    dat_bdr = np.average(rem['bdir'], weights=rem['Nhour'], axis=0) * factor
    ax.plot(dist[:-1], np.diff(dat_one), 'b-')
    ax.plot(dist[:-1], np.diff(dat_bdr - dat_one), 'r-')
    ax.plot(dist[:-1], -(np.diff(dat_one) - np.diff(dat_bdr - dat_one)), 'r--')
    
    dat = np.average(lc0['stot'], weights=rem['Nhour'], axis=0) * factor
    ax.plot(dist, dat, 'k-')
    dat = np.average(lcX['stot'], weights=rem['Nhour'], axis=0) * factor
    ax.plot(dist, dat, 'k--')
    
if True:
    region = 'wc'
    rem = remote0[region]
    lc0 = local0[region]
    lcX = localX[region]
    fig = plt.figure(10);fig.clf()
    ax = fig.subplots(1, 1)
    dist = rem['range']
    dat = np.average(rem['oneway'], weights=rem['Nhour'], axis=0) * factor
    ax.plot(dist, dat, 'b-')
    dat = np.average(rem['bdir'] - rem['oneway'], weights=rem['Nhour'], axis=0) * factor
    ax.plot(dist, dat, 'b:')
    #ax.plot(dist[:-1], np.diff(dat))
    dat = np.average(lc0['stot'], weights=rem['Nhour'], axis=0) * factor
    ax.plot(dist, np.cumsum(dat))
    dat = np.average(lcX['stot'], weights=rem['Nhour'], axis=0) * factor
    ax.plot(dist, np.cumsum(dat))

    fig = plt.figure(11);fig.clf()
    ax = fig.subplots(1, 1)
    dat_one = np.average(rem['oneway'], weights=rem['Nhour'], axis=0) * factor
    dat_bdr = np.average(rem['bdir'], weights=rem['Nhour'], axis=0) * factor
    ax.plot(dist, -np.diff(zero_pad(dat_one, (1, 0))), 'b-')
    ax.plot(dist, np.diff(zero_pad(dat_bdr - dat_one, (1, 0))), 'r-')
    ax.plot(dist, np.diff(zero_pad(dat_bdr - dat_one, (1, 0)))-np.diff(zero_pad(dat_one, (1, 0))), 'm-')
    dat = np.average(lc0['stot'], weights=rem['Nhour'], axis=0) * factor
    ax.plot(dist, dat, 'k-')
    dat = np.average(lcX['stot'], weights=rem['Nhour'], axis=0) * factor
    ax.plot(dist, dat, 'k--')
    ax.set_ylim([-20, 20])
    ax.axhline(0,color='k',linestyle=':')
    fig.savefig('fig/Flux2Sourceterms.png')


def print_results():

    print("")

    print("Remote Totals ({})".format(unit))
    print("#" * 55)
    print(("{:10s}|" + "{:>11s}" * 4)
          .format("", 'one-way', 'trad', 'unit', 'bidir'))
    print("-" * 55)
    for ireg, region in enumerate(regions):
        rt = rtot0[region]
        print("{:10s}: {oneway: 10.4g} {trad: 10.4g} {unit: 10.4g} {bdir: 10.4g}"
              .format(region, **rt))
    print("=" * 55)
    print("{:10s}: {oneway: 10.4g} {trad: 10.4g} {unit: 10.4g} {bdir: 10.4g}"
          .format('TOTAL', **rtot0['total']))

    print("")
    print("Local Baseline Totals ({})".format(unit))
    print("#" * 77)
    print(("{:10s}|" + "{:>11s}" * 6)
          .format("", 'stot', 'sin', 'sds', 'snl', 'sice', 'sbt'))
    print("-" * 77)
    for ireg, region in enumerate(regions):
        lt = ltot0[region]
        print("{:10s}: {stot: 10.4g} {sin: 10.4g} {sds: 10.4g} {snl: 10.4g} {sice: 10.4g} {sbt: 10.4g}"
              .format(region, **lt))
    print("=" * 77)
    print("{:10s}: {stot: 10.4g} {sin: 10.4g} {sds: 10.4g} {snl: 10.4g} {sice: 10.4g} {sbt: 10.4g}"
          .format('TOTAL', **ltot0['total']))

    print("")
    print("Local Potential Totals ({})".format(unit))
    print("#" * 77)
    print(("{:10s}|" + "{:>11s}" * 6)
          .format("", 'stot', 'sin', 'sds', 'snl', 'sice', 'sbt'))
    print("-" * 77)
    for ireg, region in enumerate(regions):
        lt = ltotX[region]
        print("{:10s}: {stot: 10.4g} {sin: 10.4g} {sds: 10.4g} {snl: 10.4g} {sice: 10.4g} {sbt: 10.4g}"
              .format(region, **lt))
    print("=" * 77)
    print("{:10s}: {stot: 10.4g} {sin: 10.4g} {sds: 10.4g} {snl: 10.4g} {sice: 10.4g} {sbt: 10.4g}"
          .format('TOTAL', **ltotX['total']))

print_results()
