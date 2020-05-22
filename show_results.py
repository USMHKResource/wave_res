import pyDictH5 as pdh5
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
import json
import matplotlib.pyplot as plt
import pandas as pd

regions = {'ak', 'ec', 'gm', 'prusvi', 'wc', 'hi'}
#regions = {'wc', 'ak'}

# Sum over these regions to get the total.
totregions = ['ak', 'ec', 'gm', 'prusvi', 'wc', 'hi']
#totregions = ['wc', 'ak', ]
# Initialize dictionaries for 'baseline' totals
remote0 = {}
local0 = {}
# These are time-weighted averages
rtot0Int = {}
ltot0Int = {}

rtot0Bin = {}
ltot0Bin = {}
        

with open('EPRI_totals.json') as fl:
    epri = json.load(fl)

source_terms = ['sbt', 'sds', 'snl', 'stot', 'sin', 'sice']
remote_terms = ['trad', 'oneway', 'bdir', 'unit']
    
# Initialize dictionaries for 'extraction' totals
remoteX = {}
localX = {}
# These are time-weighted averages
rtotXInt = {}
ltotXInt = {}

rtotXBin = {}
ltotXBin = {}


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

irange = 20;range_tag = 'Total'
irange = 1;range_tag = 'Nearshore'
irange = 2;range_tag = 'Nearshore-2'


def zero_pad(arr, n):
    return np.pad(arr, n, mode='constant')


def int_freq(dat, terms):
    rtnDat = deepcopy(dat)
    for ky in terms:
        rtnDat[ky] = (dat[ky] * np.diff(dat['fbins'])[None, :, None]).sum(1)
    return rtnDat

def freq_bin(dat, terms):
    rtnDat = deepcopy(dat)
    for ky in terms:
        rtnDat[ky] = (dat[ky] * np.diff(dat['fbins'])[None, :, None])
    return rtnDat

def int_freq2(dat, terms):
    df = np.diff(dat['fbins'])[None, :, None]
    f = dat['fbins'][:-1][None, :, None] + df / 2
    rtnDat = deepcopy(dat)
    for ky in terms:
        rtnDat[ky] = (dat[ky] * df).sum(1)
    return rtnDat

def freq_bin2(dat, terms):
    df = np.diff(dat['fbins'])[None, :, None]
    f = dat['fbins'][:-1][None, :, None] + df / 2
    rtnDat = deepcopy(dat)
    for ky in terms:
        rtnDat[ky] = (dat[ky] * df)
    return rtnDat


if __name__ == '__main__':

    for ireg, region in enumerate(regions):

        rd0 = remote0[region] = pdh5.load('results/freq.fcut/{}/{}.remote-totals.h5'
                                        .format('baseline', region))
        ld0 = local0[region] = pdh5.load('results/freq.fcut/{}/{}.local-totals.h5'
                                        .format('baseline', region))
        rdX = remoteX[region] = pdh5.load('results/freq.fcut/{}/{}.remote-totals.h5'
                                        .format('extraction', region))
        ldX = localX[region] = pdh5.load('results/freq.fcut/{}/{}.local-totals.h5'
                                        .format('extraction', region))
        if region == 'hi' and irange == 20:
            # The nearshore flux doesn't need this, but the EEZ does.
            # What about other ranges?
            rd0['1way'] = rd0['bdir'] - rd0['1way']
            rdX['1way'] = rdX['bdir'] - rdX['1way']
            rd0['trad'] *= -1
            rdX['trad'] *= -1

        rd0['oneway'] = rd0['1way']
        rdX['oneway'] = rdX['1way']

        rd0Int = int_freq(rd0, remote_terms)
        rdXInt = int_freq(rdX, remote_terms)


        rd0Bin = freq_bin(rd0, remote_terms)
        rdXBin = freq_bin(rdX, remote_terms)

        ld0Int = int_freq2(ld0, source_terms)
        ldXInt = int_freq2(ldX, source_terms)

        ld0Bin = freq_bin2(ld0, source_terms)
        ldXBin = freq_bin2(ldX, source_terms)

        # Integral Averages
        rtot0Int[region] = {m: (np.average(rd0Int[m][:, irange - 1],
                                        weights=rd0Int['Nhour']) * factor)
                        for m in remote_terms}
        ltot0Int[region] = {m: (np.average(ld0Int[m][:, :irange].sum(-1),
                                        weights=ld0Int['Nhour']) * factor)
                        for m in source_terms}
        rtotXInt[region] = {m: (np.average(rdXInt[m][:, irange - 1],
                                        weights=rdXInt['Nhour']) * factor)
                        for m in remote_terms}
        ltotXInt[region] = {m: (np.average(ldXInt[m][:, :irange].sum(-1),
                                        weights=ldXInt['Nhour']) * factor)
                        for m in source_terms}
        rtot0Int[region]['length'] = rd0['length'][irange - 1]
        rtotXInt[region]['length'] = rdX['length'][irange - 1]
        ltot0Int[region]['area'] = ld0['area'][:irange].sum()
        ltotXInt[region]['area'] = ldX['area'][:irange].sum()
        
        # Binned Averages
        rtot0Bin[region] = {m: (np.average(rd0Bin[m][:,:, irange - 1],
                                        weights=rd0Bin['Nhour'],axis=0) * factor)
                        for m in remote_terms}
        ltot0Bin[region] = {m: (np.average(ld0Bin[m][:,:, :irange].sum(-1),
                                        weights=ld0Bin['Nhour'],axis=0) * factor)
                        for m in source_terms}
        rtotXBin[region] = {m: (np.average(rdXBin[m][:,:, irange - 1],
                                        weights=rdXBin['Nhour'],axis=0) * factor)
                        for m in remote_terms}
        ltotXBin[region] = {m: (np.average(ldXBin[m][:,:, :irange].sum(-1),
                                        weights=ldXBin['Nhour'],axis=0) * factor)
                        for m in source_terms}

        if False:
            region = 'wc'
            binnedData,total = rtot0Bin, rtot0Int
            terms = [remote_terms,source_terms]
            fig = plt.figure()
            ax = plt.axes()
            mf = [np.mean([freq[i],freq[i+1]]) for i in range(len(freq)-1)]
            #for term in terms[0]:
            Norm = binnedData[region]['bdir']-binnedData[region]['oneway']#/total[region][term]
            ax.plot(mf,Norm,label='diff')
            plt.title(str(irange))
            plt.legend()

    for ireg, region in enumerate(totregions):
        # Calculate the offshore flux.
        rtot0Int[region]['off'] = rtot0Int[region]['bdir'] - rtot0Int[region]['oneway']
        rtotXInt[region]['off'] = rtotXInt[region]['bdir'] - rtotXInt[region]['oneway']
        # Calculate the total across all regions
        if 'total' not in rtot0Int:
            rtot0Int['total'] = deepcopy(rtot0Int[region])
            ltot0Int['total'] = deepcopy(ltot0Int[region])
            rtotXInt['total'] = deepcopy(rtotXInt[region])
            ltotXInt['total'] = deepcopy(ltotXInt[region])
        else:
            for m in ['oneway', 'off', 'unit', 'bdir', 'trad', 'length']:
                rtot0Int['total'][m] += rtot0Int[region][m]
                rtotXInt['total'][m] += rtotXInt[region][m]
            for m in ['sbt', 'sds', 'snl', 'stot', 'sin', 'sice', 'area']:
                ltot0Int['total'][m] += ltot0Int[region][m]
                ltotXInt['total'][m] += ltotXInt[region][m]

    if False:
        region = 'wc'
        binnedData,total = ltotXBin, ltotXInt
        terms = [remote_terms,source_terms]
        fig = plt.figure()
        ax = plt.axes()
        mf = [np.mean([freq[i],freq[i+1]]) for i in range(len(freq)-1)]
        for term in terms[1]:
            Norm = binnedData[region][term]/total[region][term]
            ax.plot(mf,Norm,label=term)
        plt.legend()
        plt.show()


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

    if False:
        ifreq = 5
        cut = lambda x: x[:,ifreq,:]
        region = 'wc'
        rem = remoteX[region]
        lc0 = local0[region]
        lcX = localX[region]
        fig = plt.figure(10);fig.clf()
        ax = fig.subplots(1, 1)
        dist = rem['range']
        dat = np.average(cut(rem['oneway']), weights=rem['Nhour'], axis=0) * factor
        ax.plot(dist, dat, 'b-',label='rmX oneway')
        dat = np.average(cut(rem['bdir']) - cut(rem['oneway']), weights=rem['Nhour'], axis=0) * factor
        ax.plot(dist, dat, 'b:',label='rmX bdir-oneway')
        #ax.plot(dist[:-1], np.diff(dat))
        dat = np.average(cut(lc0['stot']), weights=rem['Nhour'], axis=0) * factor
        ax.plot(dist, np.cumsum(dat),label='lc0 stot')
        dat = np.average(cut(lcX['stot']), weights=rem['Nhour'], axis=0) * factor
        ax.plot(dist, np.cumsum(dat),label='lcX stot')
        plt.title('irange '+str(irange)+', for '+region)

        fig = plt.figure(11);fig.clf()
        ax = fig.subplots(1, 1)


        dat_one = np.average(cut(rem['oneway']), weights=rem['Nhour'], axis=0) * factor
        dat_bdr = np.average(cut(rem['bdir']), weights=rem['Nhour'], axis=0) * factor
        ax.plot(dist, -np.diff(zero_pad(dat_one, (1, 0))), 'b-',label='-df(oneway)')
        ax.plot(dist, np.diff(zero_pad(dat_bdr - dat_one, (1, 0))), 'r-',label='df(bdr-oneway)')
        ax.plot(dist, np.diff(zero_pad(dat_bdr - dat_one, (1, 0)))-np.diff(zero_pad(dat_one, (1, 0))), 'm-',
                                                label='df(bdr-oneway)-df(oneway)')
        dat = np.average(cut(lc0['stot']), weights=rem['Nhour'], axis=0) * factor
        ax.plot(dist, dat, 'k-',label='lc0 stot')
        dat = np.average(cut(lcX['stot']), weights=rem['Nhour'], axis=0) * factor
        ax.plot(dist, dat, 'k--',label='lcX stot')
        ax.set_ylim([-20, 20])
        ax.axhline(0,color='k',linestyle=':')
        plt.title('Flux source terms at frequency '+str(np.mean([rem['fbins'][ifreq],rem['fbins'][ifreq+1]])))
        plt.legend()
        #fig.savefig('fig/Flux2Sourceterms.png')
        plt.show()

    if False:
        cut = lambda x: x[:,ifreq,:]
        region = 'wc'
        rem = remote0[region]
        lc0 = local0[region]
        lcX = localX[region]
        dist = rem['range']

        ifreqs = range(5,25)
        for ifreq in ifreqs:
            fig = plt.figure(11);fig.clf()
            ax = fig.subplots(1, 1)
            if True:
                bw = np.diff(rem['fbins'])[ifreq]

            dat_one = np.average(cut(rem['oneway']), weights=rem['Nhour'], axis=0) * factor * bw
            dat_bdr = np.average(cut(rem['bdir']), weights=rem['Nhour'], axis=0) * factor * bw


            off = dat_bdr - dat_one
            osf = dat_one
            print(off,osf)

            #off = off[:-1]
            #osf = osf[:-1]

            print(off.shape,osf.shape)

            # Way one
            #out_of_box = off[1::] + osf[:-1]
            #into_box = osf[1::] + off[:-1]
            # Way 2
            into_box = off[1::] + osf[:-1]
            out_of_box = osf[1::] + off[:-1]

            net_flux = into_box - out_of_box



            print('Net Flux {}, for freq: {}'.format(net_flux, ifreq))

            #ax.plot(dist, -np.diff(zero_pad(dat_one, (1, 0))), 'b-',label='-df(oneway)')
            #ax.plot(dist, np.diff(zero_pad(dat_bdr - dat_one, (1, 0))), 'r-',label='df(bdr-oneway)')
            ax.plot(dist, zero_pad(net_flux,(1,0)), 'r-', label='net flux')

            ax.plot(dist, np.diff(zero_pad(dat_bdr - dat_one, (1, 0)))-np.diff(zero_pad(dat_one, (1, ))), 'm-',
                                                    label='df(bdr-oneway)-df(oneway)')

            dat = np.average(cut(lc0['stot']), weights=rem['Nhour'], axis=0) * factor * bw 
            ax.plot(dist, dat, 'k-',label='lc0 stot')

            '''
            dat = np.average(cut(lc0['sin']), weights=rem['Nhour'], axis=0) * factor * bw
            ax.plot(dist, dat, 'r-*',label='lc0 sin')
            dat = np.average(cut(lc0['sds']), weights=rem['Nhour'], axis=0) * factor * bw
            ax.plot(dist, dat, 'r--',label='lc0 sds')

            dat = np.average(cut(lc0['snl']), weights=rem['Nhour'], axis=0) * factor * bw
            ax.plot(dist, dat, 'b--',label='lc0 snl')
            dat = np.average(cut(lc0['sbt']), weights=rem['Nhour'], axis=0) * factor * bw
            ax.plot(dist, dat, 'g--',label='lc0 sbt')
            '''

            dat = np.average(cut(lcX['stot']), weights=rem['Nhour'], axis=0) * factor * bw
            ax.plot(dist, dat, 'k--',label='lcX stot')
            ax.set_ylim([-0.5, 2.5])
            ax.axhline(0,color='k',linestyle=':')
            plt.title('Flux source terms at frequency '+str(np.mean([rem['fbins'][ifreq],rem['fbins'][ifreq+1]])))
            plt.legend()
            #fig.savefig('fig/Flux2Sourceterms_f'+str(np.mean([rem['fbins'][ifreq],rem['fbins'][ifreq+1]]))+'.png')
            plt.show()

    if True:
        pass


    def print_results():

        print("")

        print("Remote Totals ({})".format(unit))
        print("#" * 66)
        print(("{:10s}|" + "{:>11s}" * 5)
              .format("", 'one-way', 'off', 'trad', 'unit', 'bidir'))
        print("-" * 66)
        for ireg, region in enumerate(regions):
            rt = rtot0Int[region]
            print("{:10s}: {oneway: 10.4g} {off: 10.4g} {trad: 10.4g} {unit: 10.4g} {bdir: 10.4g}"
                  .format(region, **rt))
        print("=" * 66)
        print("{:10s}: {oneway: 10.4g} {off: 10.4g} {trad: 10.4g} {unit: 10.4g} {bdir: 10.4g}"
              .format('TOTAL', **rtot0Int['total']))

        # print('')
        # print("Remote ***Potential*** Totals ({})".format(unit))
        # print("#" * 66)
        # print(("{:10s}|" + "{:>11s}" * 5)
        #       .format("", 'one-way', 'off', 'trad', 'unit', 'bidir'))
        # print("-" * 66)
        # for ireg, region in enumerate(regions):
        #     rt = rtotXInt[region]
        #     print("{:10s}: {oneway: 10.4g} {off: 10.4g} {trad: 10.4g} {unit: 10.4g} {bdir: 10.4g}"
        #           .format(region, **rt))
        # print("=" * 66)
        # print("{:10s}: {oneway: 10.4g} {off: 10.4g} {trad: 10.4g} {unit: 10.4g} {bdir: 10.4g}"
        #       .format('TOTAL', **rtotXInt['total']))

        print("")
        print("Local Baseline Totals ({})".format(unit))
        print("#" * 77)
        print(("{:10s}|" + "{:>11s}" * 6)
              .format("", 'stot', 'sin', 'sds', 'snl', 'sice', 'sbt'))
        print("-" * 77)
        for ireg, region in enumerate(regions):
            lt = ltot0Int[region]
            print("{:10s}: {stot: 10.4g} {sin: 10.4g} {sds: 10.4g} {snl: 10.4g} {sice: 10.4g} {sbt: 10.4g}"
                  .format(region, **lt))
        print("=" * 77)
        print("{:10s}: {stot: 10.4g} {sin: 10.4g} {sds: 10.4g} {snl: 10.4g} {sice: 10.4g} {sbt: 10.4g}"
              .format('TOTAL', **ltot0Int['total']))

        print("")
        print("Local Potential Totals ({})".format(unit))
        print("#" * 77)
        print(("{:10s}|" + "{:>11s}" * 6)
              .format("", 'stot', 'sin', 'sds', 'snl', 'sice', 'sbt'))
        print("-" * 77)
        for ireg, region in enumerate(regions):
            lt = ltotXInt[region]
            print("{:10s}: {stot: 10.4g} {sin: 10.4g} {sds: 10.4g} {snl: 10.4g} {sice: 10.4g} {sbt: 10.4g}"
                  .format(region, **lt))
        print("=" * 77)
        print("{:10s}: {stot: 10.4g} {sin: 10.4g} {sds: 10.4g} {snl: 10.4g} {sice: 10.4g} {sbt: 10.4g}"
              .format('TOTAL', **ltotXInt['total']))

    print_results()

    # Write-out the results.
    df0 = pd.DataFrame(rtot0Int).T.join(pd.DataFrame(ltot0Int).T)
    dfX = pd.DataFrame(rtotXInt).T.join(pd.DataFrame(ltotXInt).T)

    df0 = df0.reindex(['wc', 'ec', 'hi', 'ak', 'gm', 'prusvi', 'total'])
    dfX = dfX.reindex(['wc', 'ec', 'hi', 'ak', 'gm', 'prusvi', 'total'])

    col_order = ['length', 'oneway', 'off', 'bdir', 'trad', 'unit',
                 'area', 'sbt', 'sds', 'sice', 'sin', 'snl', 'stot']
    
    df0.to_csv(
        'results/{}_Baseline_Results.csv'.format(range_tag),
        columns=col_order)
    dfX.to_csv(
        'results/{}_Extraction_Results.csv'.format(range_tag),
        columns=col_order)
