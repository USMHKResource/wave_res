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

# calculate fluxes
def calc_fluxes(dat_one, dat_bdr):
    off = dat_bdr - dat_one
    osf = dat_one
    into_box = off[1::] + osf[:-1]
    out_of_box = osf[1::] + off[:-1]
    return into_box - out_of_box

Sum = lambda x: np.sum(np.array(x),axis=0)

'''
PLOT ACTIVATING VARIABLES
'''
noFcCut = True # run calculations with no fc cutting
FcCut = True # run calculations with fc cutting
singleFreq = False # calculate data with single frequency, if False freq are summed
freqBinPlot = False # plot fluxes and sources for each frequency bin
sumPlot = False # plot fluxes and sources over summed frequency ranges
absDiff = False # plot the absolute difference between lc0-rem0, lcX-remX in both cases
RMSE = True # calculate root mean squared calculations


if True:
    '''
    Normal comaprisons of flux v source terms and distance from shore
    '''
    cut = lambda x: x[:,ifreq,:]
    region = 'wc'
    rem0, remX , lc0, lcX = (remote0[region],remoteX[region],
                                 local0[region], localX[region])

    fluxes0, fluxesX, sources0, sourcesX = [],[],[],[]
    
    dist,fbins,weights = rem0['range'],rem0['fbins'],rem0['Nhour']

    freqs = [np.mean(f) for f in zip(fbins,fbins[1:])]
    for ifreq,freq in enumerate(freqs):

        if True:
            bw = np.diff(fbins)[ifreq]

        if singleFreq:
            # Chooses a particular frequency to work with
            dat_one0 = np.average(cut(rem0['oneway']), weights=weights, axis=0) * factor * bw
            dat_bdr0 = np.average(cut(rem0['bdir']), weights=weights, axis=0) * factor * bw
            dat_oneX = np.average(cut(remX['oneway']), weights=weights, axis=0) * factor * bw
            dat_bdrX = np.average(cut(remX['bdir']), weights=weights, axis=0) * factor * bw
            dat0 = np.average(cut(lc0['stot']), weights=weights, axis=0) * factor * bw
            datX = np.average(cut(lcX['stot']), weights=weights, axis=0) * factor * bw
            sources0.append(dat0)
            sourcesX.append(datX)
            flux0, fluxX = (calc_fluxes(dat_one0, dat_bdr0),
                            calc_fluxes(dat_oneX, dat_bdrX))
            fluxes0.append(flux0)
            fluxesX.append(fluxX)

        else:
            # total sum over all frequencies
            dat_one0 = np.average(rem0['oneway'], weights=weights, axis=0).sum(axis=0) * factor * bw
            dat_bdr0 = np.average(rem0['bdir'], weights=weights, axis=0).sum(axis=0) * factor * bw
            dat_oneX = np.average(remX['oneway'], weights=weights, axis=0).sum(axis=0) * factor * bw
            dat_bdrX = np.average(remX['bdir'], weights=weights, axis=0).sum(axis=0) * factor * bw
            dat0 = np.average(lc0['stot'], weights=weights, axis=0).sum(axis=0) * factor * bw
            datX = np.average(lcX['stot'], weights=weights, axis=0).sum(axis=0) * factor * bw
        
            sources0.append(dat0)
            sourcesX.append(datX)
            flux0, fluxX = (calc_fluxes(dat_one0, dat_bdr0),
                            calc_fluxes(dat_oneX, dat_bdrX))
            fluxes0.append(flux0)
            fluxesX.append(fluxX)

        if freqBinPlot:
            fig = plt.figure(11);fig.clf()
            ax = fig.subplots(1, 1)
            ax.plot(dist, dat0, 'r-',label='lc0 stot')
            ax.plot(dist, datX, 'r--',label='lcX stot')
            ax.plot(dist, zero_pad(flux0,(1,0)), 'b-', label='net flux 0')
            ax.plot(dist, zero_pad(fluxX,(1,0)), 'b--', label='net flux X')
            if singleFreq:
                ax.set_ylim([-0.5, 2.5])
            else:
                ax.set_ylim([-0.5, 50])
            ax.axhline(0,color='k',linestyle=':')
            plt.title(f'Flux / Source terms, f={freq:.5f} Hz, summed = {freqBinPlot}')
            plt.legend()
            #fig.savefig('fig/Flux2Sourceterms_f'+str(np.mean([rem['fbins'][ifreq],rem['fbins'][ifreq+1]]))+'.png')
            plt.show()

    sources0, sourcesX, fluxes0, fluxesX = (Sum(sources0), Sum(sourcesX), 
                                                Sum(fluxes0),Sum(fluxesX))
    noCut = [sources0, sourcesX, fluxes0, fluxesX]

    if sumPlot: # Plots the total sum over frequencies for fluxes and sources.

        fig = plt.figure(11,figsize=(8,4));fig.clf()
        ax = fig.subplots(1, 1)
        ax.plot(dist, sources0, 'r-', label='lc0 stot')
        ax.plot(dist, sourcesX, 'r--', label='lcX stot')
        ax.plot(dist, zero_pad(fluxes0,(1,0)), 'b-', label='net flux 0')
        ax.plot(dist, zero_pad(fluxesX,(1,0)), 'b--', label='net flux X')
        ax.set_ylim([-100, 600])
        ax.axhline(0,color='k',linestyle=':')
        plt.title(f'Integral Flux / Source terms, no Frequaency cutoff')
        plt.legend(loc=1)
        #fig.savefig(f'fig/Flux2Sourceterms_fc.png')
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
    # Convenience function for the time averaging
    avg = lambda x,y: np.average(x,weights=weights,axis=0).sum(axis=0)*factor*y

    region = 'wc'
    rem0, remX, lc0, lcX = remote0[region], remoteX[region], local0[region], localX[region]
    dist,fbins,weights = rem0['range'],rem0['fbins'],rem0['Nhour']

    freqs = [np.mean(f) for f in zip(fbins,fbins[1:])]
    fluxes0, fluxesX, sources0, sourcesX = [],[],[],[]
    
    # Loop through all frequencies
    for ifreq,freq in enumerate(freqs):
        bandWidth = np.diff(fbins)[ifreq] # bandwidth scaling for each freq bin
        
        # Calculate the cutoff frequency for the given center frequency.
        # Logic: energy produced from waves of this ~frequency are allowed to 
        #        transfer up to fc. So the fluxes calculated may be artificially
        #        cut by WW3. It would seem then that we should exclude the fluxes
        #        AND source terms at these frequencies.
        fc = 2.5*freq
        ifc = freqs <= fc

        # cut the source and flux terms based on the fc
        lc0Sums = avg(lc0['stot'][:,ifc,:],bandWidth)
        lcXSums = avg(lcX['stot'][:,ifc,:],bandWidth)
        dat_one0, dat_bdr0 = (avg(rem0['oneway'][:,ifc,:],bandWidth),
                            avg(rem0['bdir'][:,ifc,:],bandWidth))
        dat_oneX, dat_bdrX = (avg(remX['oneway'][:,ifc,:],bandWidth),
                            avg(remX['bdir'][:,ifc,:],bandWidth))

        # append over frequency
        sources0.append(lc0Sums)
        sourcesX.append(lcXSums)

        flux0, fluxX = (calc_fluxes(dat_one0, dat_bdr0),
                        calc_fluxes(dat_oneX, dat_bdrX))
        fluxes0.append(flux0)
        fluxesX.append(fluxX)

        if freqBinPlot: # Plots fluxes and sources for each frequency and fc
            fig = plt.figure(11);fig.clf()
            ax = fig.subplots(1, 1)
            ax.plot(dist, lc0Sums, 'r-', label='lc0 stot')
            ax.plot(dist, lcXSums, 'r--', label='lcX stot')
            ax.plot(dist, zero_pad(flux0,(1,0)), 'b-', label='net flux 0')
            ax.plot(dist, zero_pad(fluxX,(1,0)), 'b--', label='net flux X')
            ax.set_ylim([-0.5, 50])
            ax.axhline(0,color='k',linestyle=':')
            plt.title(f'Flux / Source terms, f = {freq:.5f} Hz, f_c = {fc:.5f} Hz')
            plt.legend()
            #fig.savefig(f'fig/Flux2Sourceterms_f={ffreq:.3f}Hz.png')
            plt.show()
        
    sources0, sourcesX, fluxes0, fluxesX = (Sum(sources0), Sum(sourcesX), 
                                                Sum(fluxes0),Sum(fluxesX))
    Cut = [sources0, sourcesX, fluxes0, fluxesX]

    if sumPlot: # Plots the total sum over frequencies for fluxes and sources.

        fig = plt.figure(11,figsize=(8,4));fig.clf()
        ax = fig.subplots(1, 1)
        ax.plot(dist, sources0, 'r-', label='lc0 stot')
        ax.plot(dist, sourcesX, 'r--', label='lcX stot')
        ax.plot(dist, zero_pad(fluxes0,(1,0)), 'b-', label='net flux 0')
        ax.plot(dist, zero_pad(fluxesX,(1,0)), 'b--', label='net flux X')
        ax.set_ylim([-100, 600])
        ax.axhline(0,color='k',linestyle=':')
        plt.title(f'Integral Flux / Source terms')
        plt.legend(loc=1)
        #fig.savefig(f'fig/Flux2Sourceterms_fc.png')
        plt.show()

noCut0diff = np.array(noCut[0])-zero_pad(noCut[2],(1,0))
noCutXdiff = np.array(noCut[1])-zero_pad(noCut[3],(1,0))
Cut0diff = np.array(Cut[0])-zero_pad(Cut[2],(1,0))
CutXdiff = np.array(Cut[1])-zero_pad(Cut[3],(1,0))

# cross terms to look at the extracted remote and local params
noCutCrossXdiff = np.array(noCut[0])-zero_pad(noCut[3],(1,0))
cutCrossXdiff = np.array(Cut[0])-zero_pad(Cut[3],(1,0))

if absDiff:
    a = np.abs 

    fig = plt.figure(11,figsize=(8,4));fig.clf()
    ax = fig.subplots(1, 1)
    ax.plot(dist, a(noCut0diff), 'k-', label='no fc Cut (lc0-rem0)')
    ax.plot(dist, a(noCutXdiff), 'k--', label='no fc Cut (lcX-remX)')
    ax.plot(dist, a(Cut0diff), 'g-', label='fc Cut (lc0-rem0)')
    ax.plot(dist, a(CutXdiff), 'g--', label='fc Cut (lcX-remX)')
    ax.plot(dist, a(noCutCrossXdiff), 'm-', label='no fc Cut (lc0-remX)')
    ax.plot(dist, a(cutCrossXdiff), 'm--', label='fc Cut (lc0-remX)')
    ax.set_ylim([-100, 600])
    ax.axhline(0,color='k',linestyle=':')
    plt.title(f'Absolute Difference')
    plt.legend(loc=1)
    #fig.savefig(f'fig/Flux2Sourceterms_fc.png')
    plt.show()

if RMSE:
    rmse = lambda x: np.sqrt(np.mean(x**2))

    print(f'With no freq cut RMSE(lc0-rem0) = {rmse(noCut0diff[:-1])}')
    print(f'With no freq cut RMSE(lcX-remX) = {rmse(noCutXdiff[:-1])}')
    print(f'With freq cut RMSE(lc0-rem0) = {rmse(Cut0diff[:-1])}')
    print(f'With freq cut RMSE(lcX-remX) = {rmse(CutXdiff[:-1])}')
    print(f'With no freq cut the cross term RMSE(lc0-remX) = {rmse(noCutCrossXdiff[:-1])}')
    print(f'With freq cut the cross term RMSE(lc0-remX) = {rmse(cutCrossXdiff[:-1])}')