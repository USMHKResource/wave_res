import pyDictH5 as pdh5
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
import json
import matplotlib.pyplot as plt


#regions = {'ak', 'at', 'prusvi', 'wc', 'hi'}
regions = {'wc'}

# Sum over these regions to get the total.
# totregions = ['ak', 'at', 'prusvi', 'wc', 'hi']
totregions = ['wc']
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

iranges = range(10,20)


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


for irange in iranges:
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

        rd0Int = int_freq(rd0, remote_terms)
        rdXInt = int_freq(rdX, remote_terms)


        rd0Bin = freq_bin(rd0, remote_terms)
        rdXBin = freq_bin(rdX, remote_terms)

        ld0Int = int_freq2(ld0, source_terms)
        ldXInt = int_freq2(ldX, source_terms)

        ld0Bin = freq_bin2(ld0, source_terms)
        ldXBin = freq_bin2(ldX, source_terms)
        
        # Intergal Averages
        rtot0Int[region] = {m: (np.average(rd0Int[m][:, irange],
                                        weights=rd0Int['Nhour']) * factor)
                        for m in remote_terms}
        ltot0Int[region] = {m: (np.average(ld0Int[m][:, :irange].sum(-1),
                                        weights=ld0Int['Nhour']) * factor)
                        for m in source_terms}
        rtotXInt[region] = {m: (np.average(rdXInt[m][:, irange],
                                        weights=rdXInt['Nhour']) * factor)
                        for m in remote_terms}
        ltotXInt[region] = {m: (np.average(ldXInt[m][:, :irange].sum(-1),
                                        weights=ldXInt['Nhour']) * factor)
                        for m in source_terms}

        # Binned Averages
        rtot0Bin[region] = {m: (np.average(rd0Bin[m][:,:, irange],
                                        weights=rd0Bin['Nhour'],axis=0) * factor)
                        for m in remote_terms}
        ltot0Bin[region] = {m: (np.average(ld0Bin[m][:,:, :irange].sum(-1),
                                        weights=ld0Bin['Nhour'],axis=0) * factor)
                        for m in source_terms}
        rtotXBin[region] = {m: (np.average(rdXBin[m][:,:, irange],
                                        weights=rdXBin['Nhour'],axis=0) * factor)
                        for m in remote_terms}
        ltotXBin[region] = {m: (np.average(ldXBin[m][:,:, :irange].sum(-1),
                                        weights=ldXBin['Nhour'],axis=0) * factor)
                        for m in source_terms}


## Normal comaprisons of flux v source terms and distance from shore
if False:
    cut = lambda x: x[:,ifreq,:]
    region = 'wc'
    rem, lc0, lcX = remoteX[region], local0[region], localX[region]
    dist,fbins,weights = rem['range'],rem['fbins'],rem['Nhour']

    freqs = [np.mean(f) for f in zip(fbins,fbins[1:])]
    for ifreq,freq in enumerate(freqs):
        fig = plt.figure(11);fig.clf()
        ax = fig.subplots(1, 1)
        if True:
            bw = np.diff(fbins)[ifreq]

        dat_one = np.average(cut(rem['oneway']), weights=weights, axis=0) * factor * bw
        dat_bdr = np.average(cut(rem['bdir']), weights=weights, axis=0) * factor * bw
        
        off = dat_bdr - dat_one
        osf = dat_one

        into_box = off[1::] + osf[:-1]
        out_of_box = osf[1::] + off[:-1]

        net_flux = into_box - out_of_box

        print(f'Net Flux {net_flux}, for freq: {freq} Hz')

        #ax.plot(dist, -np.diff(zero_pad(dat_one, (1, 0))), 'b-',label='-df(oneway)')
        #ax.plot(dist, np.diff(zero_pad(dat_bdr - dat_one, (1, 0))), 'r-',label='df(bdr-oneway)')
        ax.plot(dist, zero_pad(net_flux,(1,0)), 'r-', label='net flux')

        ax.plot(dist, np.diff(zero_pad(dat_bdr - dat_one, (1, 0)))-np.diff(zero_pad(dat_one, (1, 0))), 'm-',
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
    '''
    Here I want to look at the sum of the source terms when compared to fluxes based
    on the WW3 cut off frequency. 

    This should have the effect of removing high energy in the source terms at high frequency
    Or at least thats the idea. Cutoff is based on:
        f_c = f_FM * f_m = 2.5 * f_m

    f_m == frequency modelled I believe
    '''
    cut = lambda x: x[:,ifreq,:]
    avg = lambda x,y: np.average(x,weights=weights,axis=0).sum(axis=0)*factor*y

    region = 'wc'
    rem, lc0, lcX = remoteX[region], local0[region], localX[region]
    dist,fbins,weights = rem['range'],rem['fbins'],rem['Nhour']

    freqs = [np.mean(f) for f in zip(fbins,fbins[1:])]
    fluxes, sources0, sourcesX = [],[],[]
    for ifreq,freq in enumerate(freqs):
        bandWidth = np.diff(fbins)[ifreq]
        
        fc = 2.5*freq
        ifc = freqs <= fc

        lc0Sums = avg(lc0['stot'][:,ifc,:],bandWidth)
        lcXSums = avg(lcX['stot'][:,ifc,:],bandWidth)
        dat_one, dat_bdr = (avg(rem['oneway'][:,ifc,:],bandWidth),
                            avg(rem['bdir'][:,ifc,:],bandWidth))

        off = dat_bdr - dat_one
        osf = dat_one
        into_box = off[1::] + osf[:-1]
        out_of_box = osf[1::] + off[:-1]
        net_flux = into_box - out_of_box

        sources0.append(lc0Sums)
        sourcesX.append(lcXSums)
        fluxes.append(net_flux)

        if False:
            fig = plt.figure(11);fig.clf()
            ax = fig.subplots(1, 1)
            ax.plot(dist, lc0Sums, 'r-', label='lc0 stot')
            ax.plot(dist, zero_pad(net_flux,(1,0)), 'b-', label='net flux')
            ax.set_ylim([-0.5, 50])
            ax.axhline(0,color='k',linestyle=':')
            plt.title(f'Flux / Source terms, f = {freq:.5f} Hz, f_c = {fc} Hz')
            plt.legend()
            #fig.savefig('fig/Flux2Sourceterms_f'+str(np.mean([rem['fbins'][ifreq],rem['fbins'][ifreq+1]]))+'.png')
            plt.show()
        
    if True:
        Sum = lambda x: np.sum(np.array(x),axis=0)
        print(np.cumsum(sources0,axis=0))
        sources0, sourcesX, fluxes = Sum(sources0), Sum(sourcesX), Sum(fluxes)

        fig = plt.figure(11);fig.clf()
        ax = fig.subplots(1, 1)
        ax.plot(dist, sources0, 'r-', label='lc0 stot')
        ax.plot(dist, sourcesX, 'g-', label='lcX stot')
        ax.plot(dist, zero_pad(fluxes,(1,0)), 'b-', label='net flux')
        ax.set_ylim([-0.5, 800])
        ax.axhline(0,color='k',linestyle=':')
        plt.title(f'Integral Flux / Source terms')
        plt.legend()
        fig.savefig(f'fig/Flux2Sourceterms_fc.png')
        plt.show()