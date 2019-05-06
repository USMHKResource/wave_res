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
baseComparison = True # run the base comparison
noFcCut = False # run calculations with no fc cutting
FcCut = False # run calculations with fc cutting
singleFreq = True # calculate data with single frequency, if False freq are summed
freqBinPlot = False # plot fluxes and sources for each frequency bin
freqBinPlotFull = True # plot fluxes and sources for each frequency bin
sumPlot = False # plot fluxes and sources over summed frequency ranges
sumPlotFull = False # plot fluxes and sources over summed frequency ranges
absDiff = False # plot the absolute difference between lc0-rem0, lcX-remX in both cases
RMSE = False # calculate root mean squared calculations
cumulativeSums = True
vSumPlot = True


if baseComparison:
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
            bwTot = np.diff(fbins)[:,np.newaxis]

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
            '''
            # no Bandwidth scaling
            dat_one0 = np.average(rem0['oneway'], weights=weights, axis=0).sum(axis=0) * factor# * bw
            dat_bdr0 = np.average(rem0['bdir'], weights=weights, axis=0).sum(axis=0) * factor #* bw
            dat_oneX = np.average(remX['oneway'], weights=weights, axis=0).sum(axis=0) * factor# * bw
            dat_bdrX = np.average(remX['bdir'], weights=weights, axis=0).sum(axis=0) * factor #* bw
            dat0 = np.average(lc0['stot'], weights=weights, axis=0).sum(axis=0) * factor #* bw
            datX = np.average(lcX['stot'], weights=weights, axis=0).sum(axis=0) * factor #* bw
            '''
            # Correct application of bandwidth scaling
            dat_one0 = (np.average(rem0['oneway'], weights=weights, axis=0)/bwTot).sum(axis=0) * factor
            dat_bdr0 = (np.average(rem0['bdir'], weights=weights, axis=0)/bwTot).sum(axis=0) * factor
            dat_oneX = (np.average(remX['oneway'], weights=weights, axis=0)/bwTot).sum(axis=0) * factor
            dat_bdrX = (np.average(remX['bdir'], weights=weights, axis=0)/bwTot).sum(axis=0) * factor 
            dat0 = (np.average(lc0['stot'], weights=weights, axis=0)/bwTot).sum(axis=0) * factor 
            datX = (np.average(lcX['stot'], weights=weights, axis=0)/bwTot).sum(axis=0) * factor
            
            sources0.append(dat0)
                
            sourcesX.append(datX)
            flux0, fluxX = (calc_fluxes(dat_one0, dat_bdr0),
                            calc_fluxes(dat_oneX, dat_bdrX))
            fluxes0.append(flux0)
            fluxesX.append(fluxX)
        

        if freqBinPlotFull:
            fig = plt.figure(11,figsize=[8,3]);fig.clf()
            ax = fig.subplots(1, 1)
            ax.plot(dist, dat0, 'r-',label='Baseline Source Totals')
            ax.plot(dist, datX, 'r--',label='Extraction Source Totals')
            ax.plot(dist, zero_pad(flux0,(1,0)), 'b-', label='Baseline Net Flux')
            ax.plot(dist, zero_pad(fluxX,(1,0)), 'b--', label='Extraction Net Flux')
            if singleFreq:
                ax.set_ylim([-0.5, 2.5])
            else:
                ax.set_ylim([-0.5, 100000])
            ax.set_xlabel('Distance form Shore (miles)')
            ax.set_ylabel('')
            ax.axhline(0,color='k',linestyle=':')
            plt.title(f'Flux / Source terms, f={freq:.5f} Hz')#, summed = {freqBinPlot}')
            plt.legend()
            plt.grid()
            plt.tight_layout()
            #fig.savefig('fig/Flux2Sourceterms_f'+str(np.mean([rem['fbins'][ifreq],rem['fbins'][ifreq+1]]))+'.png')
            plt.show()

    
    sources0, sourcesX, fluxes0, fluxesX = (Sum(sources0), Sum(sourcesX), 
                                                Sum(fluxes0),Sum(fluxesX))
    noCut = [sources0, sourcesX, fluxes0, fluxesX]

    print(f'Numpy Post {fluxes0}')

    if sumPlotFull: # Plots the total sum over frequencies for fluxes and sources.

        fig = plt.figure(11,figsize=(8,4));fig.clf()
        ax = fig.subplots(1, 1)
        ax.plot(dist, sources0, 'r-', label='Baseline Source Totals')
        ax.plot(dist, sourcesX, 'r--', label='Extraction Source Totals')
        ax.plot(dist, zero_pad(fluxes0,(1,0)), 'b-', label='Baseline Net Flux')
        ax.plot(dist, zero_pad(fluxesX,(1,0)), 'b--', label='Extraction Net Flux')
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
    avg = lambda x,y: (np.average(x,weights=weights,axis=0)*y).sum(axis=0)*factor

    region = 'wc'
    rem0, remX, lc0, lcX = remote0[region], remoteX[region], local0[region], localX[region]
    dist,fbins,weights = rem0['range'],rem0['fbins'],rem0['Nhour']

    freqs = [np.mean(f) for f in zip(fbins,fbins[1:])]
    fluxes0, fluxesX, sources0, sourcesX = [],[],[],[]
    
    # Loop through all frequencies
    for ifreq,freq in enumerate(freqs):
        # LOL No: bandWidth = np.diff(fbins)[ifreq] # bandwidth scaling for each freq bin
        
        # Calculate the cutoff frequency for the given center frequency.
        # Logic: energy produced from waves of this ~frequency are allowed to 
        #        transfer up to fc. So the fluxes calculated may be artificially
        #        cut by WW3. It would seem then that we should exclude the fluxes
        #        AND source terms at these frequencies.
        fc = 1*freq
        ifc = freqs <= fc
        bandWidth = np.diff(fbins)[ifc,np.newaxis] 
        print(bandWidth)

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
            fig = plt.figure(11,figsize=[8,3]);fig.clf()
            ax = fig.subplots(1, 1)
            ax.plot(dist, lc0Sums, 'r-', label='Baseline Source Totals')
            ax.plot(dist, lcXSums, 'r--', label='Extraction Source Totals')
            ax.plot(dist, zero_pad(flux0,(1,0)), 'b-', label='Baseline Net Flux')
            ax.plot(dist, zero_pad(fluxX,(1,0)), 'b--', label='Baseline Net Flux')
            ax.set_ylim([-1, 10])
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

# Calculate differences
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

if cumulativeSums:
    '''
    Here I want to look at the convergence (if any) of the cummulative sums of flux and
    source terms. 

    It seems from looking at the individual frequencies that this might occur 
    '''

    '''
    # This is a curiousity test of the speed of parallel numba execution on cpu and gpu
    # vs numpy functionality. Mainly to satisfy an itch to stay in numpy framework but
    # faster. 
    # 
    # It looks like (at least on a nvidia Quadro P600) the gpu is not the greatest for raw calculation,
    # there are probably some optimizing techniques for low memory but massive calculation 
    # requirements.
    # 
    # Also assigning values to a prebuilt array within a njit function is much faster that returning
    # the array. aka types = 'void(..., ...)' 

    from timeit import default_timer as timer
    import numba as nb
    from numba import cuda, njit, prange

    types = 'int32[:,::1](int32[::1],int32[::1],int32[:,::1])'
    gutypes = 'void(int32[::1],int32[::1],int32[:,::1])'
    @njit(types,parallel=True,nogil=True,cache=True)
    def Average(freqs,weights,output):
        lf, lw = len(freqs),len(weights)
        for i in prange(lf):
            for j in prange(lw):
                output[i,j] = freqs[i]*weights[j]
        return output
                
    gutypes = 'void(int32[::1],int32[::1],int32[:,::1],int32,int32)'
    @cuda.jit(gutypes)
    def cudaAverage(freqs,weights,output,lf,lw):
        for i in prange(lf):
            for j in prange(lw):
                output[i,j] = freqs[i]*weights[j]
    
    arraySize = 5000
    freqs,weights,factor,bw = np.arange(arraySize),np.arange(arraySize),np.arange(arraySize),2
    output = np.zeros_like(np.outer(freqs,weights))

    start = timer()
    outnb = Average(freqs,weights,output)
    timing = timer()-start
    print(f'outNB {outnb.shape}, {timing}')
    
    start = timer()
    lf, lw = len(freqs),len(weights)
    cudaAverage(freqs,weights,output,lf,lw)
    timing = timer()-start
    print(f'cuda {output.shape}, {timing}')
    
    start = timer()
    outnp = np.multiply.outer(freqs,weights)
    timing = timer()-start
    print(f'outNP {outnp.shape}, {timing}')

    '''

    
    from timeit import default_timer as timer
    from numba import jit, njit, prange, typeof
    parallel = False # here for numba execution. Most cases are False

    wmtypes = 'void(float32[:,:,::1],uint16[::1],float64,float32[::1],float32[::1],float32[::1],float32[::1])'
    @njit(wmtypes,parallel=parallel,nogil=True,cache=True,fastmath=True)
    def weighted_Average(data,weights,factor,bandWidth,frequencies,frange,output):
        times, freqs, dists = data.shape[0], frequencies, data.shape[-1]
        Min, Max = frange
        weightSum = np.sum(weights)
        for dist in prange(dists):
            for freq in prange(freqs.shape[0]):
                f = freqs[freq]
                if f >= Min and f <= Max:
                    tempfreq = 0
                    for time in prange(times):
                        tempfreq += data[time,freq,dist]*weights[time]
                    output[dist] += (tempfreq/weightSum)*bandWidth[freq]*factor
                else:
                    output[dist] += 0 
            
    

    fbtypes = 'void(float32[::1],float32[::1])'
    # Note here, best preformance comes with parallel=False
    # arrays are too small to benefit
    @njit(fbtypes,parallel=parallel,nogil=True,cache=True,fastmath=True)
    def freqsFromBins(fbins,frequencies):
        for binEdge in prange(fbins.shape[0]-1):
            frequencies[binEdge] = (fbins[binEdge+1]+fbins[binEdge])/2
                        
    cftypes = 'void(float32[::1],float32[::1],float32[::1])'
    @njit(cftypes,parallel=parallel,nogil=True,cache=True,fastmath=True)
    def calc_fluxes(ow, bdr, flux):
        flux[:] = ((bdr[1::]-ow[1::])+ow[:-1])-(ow[1::]+(bdr[:-1]-ow[:-1]))
    
        

    ''' EXECUTION FLAGS '''
    Time = True
    numba = True
    numpy  = False
    cumSum = False
    single = False

    import seaborn
    frange = np.array([0,0.17],dtype='float32')

    region = 'wc'
    rem0, remX, lc0, lcX = remote0[region], remoteX[region], local0[region], localX[region]
    dist,fbins,weights = rem0['range'],rem0['fbins'],rem0['Nhour']
    fluxes0, fluxesX, sources0, sourcesX = [],[],[],[]
    bandwidths = np.diff(fbins)

    if numba:
        start = timer()
        frequencies = np.zeros_like(fbins[:-1])
        freqsFromBins(fbins,frequencies)
        timing = timer()-start
        if Time: print(f'Calc Frequencies numba {frequencies.shape}, {timing}')
        print(frequencies)

    if numpy:
        start = timer()    
        freqs = [np.mean(f) for f in zip(fbins,fbins[1:])]
        timing = timer()-start
        if Time: print(f'Calc Frequencies numpy {np.array(freqs).shape}, {timing}')

    if numba and single:
        wAvg = np.zeros_like(dist,dtype='float32')
        start = timer()
        weighted_Average(lc0['stot'],weights,factor,bandwidths,frequencies,frange,wAvg)
        timing = timer()-start
        if Time: print(f'Weighted Average numba {wAvg.shape}, {timing}')

    if numpy and single:
        start = timer()
        for ifreq,freq in enumerate(freqs):
            lc0Sums = avg(lc0['stot'][:,ifc,:],bandwidths[ifreq])
            sources0.append(lc0Sums)
        sources0 = Sum(sources0)
        timing = timer()-start
        if Time: print(f'Weighted Average numpy {sources0.shape}, {timing}')

    if numba and not single:
        if cumSum:
            pass
        else:
            start = timer()

            z = np.zeros_like(dist,dtype='float32')
            ow0, owX, bdr0, bdrX, sources0, sourcesX = (z.copy() for i in range(6))
            fluxes0, fluxesX = z.copy()[:-1], z.copy()[:-1]
            dataSets = ((lc0['stot'],sources0),(lcX['stot'],sourcesX),
                            (rem0['oneway'],ow0),(remX['oneway'],owX),
                            (rem0['bdir'],bdr0),(remX['bdir'],bdrX))
            for data in dataSets:
                weighted_Average(data[0],weights,factor,bandwidths,frequencies,frange,data[1])
                
            calc_fluxes(ow0,bdr0,fluxes0)
            calc_fluxes(owX,bdrX,fluxesX)

            if Time: print(f'numba terms calculated in {timing} s')

            Max = np.max([np.max(sources0),np.max(sourcesX),
                            np.max(zero_pad(fluxes0,(1,0))),
                            np.max(zero_pad(fluxesX,(1,0)))])

            if vSumPlot: # Plots the total sum over frequencies for fluxes and sources.

                fig = plt.figure(11,figsize=(8,4));fig.clf()
                ax = fig.subplots(1, 1)
                ax.plot(dist, sources0/Max, 'r-', label='BaseLine Source Totals')
                ax.plot(dist, sourcesX/Max, 'r--', label='Extraction Source Totals')
                ax.plot(dist, zero_pad(fluxes0,(1,0))/Max,'b-', label='Baseline Net Flux')
                ax.plot(dist, zero_pad(fluxesX,(1,0))/Max,'b--', label='Extraction Net Flux')
                ax.set_ylim([0, 1.1])
                ax.set_xlabel(f'Distance from Shore (miles)')
                ax.set_ylabel(f'Normalized Sums over Frequency')
                ax.axhline(0,color='k',linestyle=':')
                plt.grid()
                plt.title(f'Integral Flux and Source terms \n'+\
                    f' Frequency Range: {frange[0]:.3f} - {frange[1]:.3f} Hz')
                plt.legend(loc=1)
                plt.tight_layout()
                #fig.savefig(f'fig/Flux2Sourceterms_fc.png')
                plt.show()



    if numpy and not single:
        # Loop through all frequencies
        for ifreq,freq in enumerate(freqs):
            bandWidth = bandwidths[ifreq]

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

        print(fluxes0)