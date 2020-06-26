import pyDictH5 as pdh5
import pathlib2 as pl
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumtrapz
from collections import defaultdict
import base
plt.ion()
import pdb

flag = defaultdict(lambda: False)

flag['show remote'] = True
flag['show local'] = True
#flag['show each'] = True

srcdir = pl.Path('./results/freq.fcut/').resolve()

colors = {'ak': 'k', 'wc': 'b', 'hi': 'g',
          'ec': 'r', 'prusvi': 'purple', 'gm': 'm'}

labels = {'wc': 'West Coast',
          'hi': 'Hawaii',
          'ec': 'East Coast',
          'ak': 'Alaska',
          'gm': 'Gulf of Mexico'}

tag = '02'
plot_regions = ['wc', 'hi', 'ec', 'ak', 'gm']

# tag = '03'
# plot_regions = ['wc', 'hi', 'ec']


def _translate_fbins(fbins):
    Tbins = fbins ** -1
    df = np.diff(fbins)
    f = fbins[:-1] + df / 2
    dT = -np.diff(Tbins)
    T = f ** -1
    return Tbins, T, dT, f, df

    
def calc_cpsd_points(spec, fbins, pvals):
    pvals = 1 - np.array(pvals)
    Tbins, T, dT, f, df = _translate_fbins(fbins)
    cspec = np.cumsum(spec * dT)
    cspec /= cspec[-1]
    Tvals = np.interp(pvals, cspec, Tbins[1:])
    Svals = np.interp(Tvals, T[::-1], spec[::-1])
    return Tvals, Svals
    


def spec_freq2period(spec, fbins, norm=False):
    """Convert the units of a frequency spectrum (W/Hz) to a
    period-spectrum (W/s).
    """
    
    Tbins, T, dT, f, df = _translate_fbins(fbins)
    specT = spec * (f ** 2)
    
    # Check that the integrals match!
    If = np.sum(spec * df)
    IT = np.sum(specT * dT)
    assert (IT - If) / If < 0.01, "The integrals don't match!?"

    if norm:
        specT /= IT
    #specT = np.hstack([[0], specT])
    #pdb.set_trace()
    return specT, T


def int_period(specT, fbins, axis=0):
    _, _, dT, _, _ = _translate_fbins(fbins)
    #specT = specT[1:]
    return (specT * dT).sum(axis)
    

############## Begin Plotting ############

# Initialize figures
figL0 = plt.figure(202)
figL0.clf()
figL0, axL0 = plt.subplots(1, 1, num=figL0.number)

figR0 = plt.figure(102)
figR0.clf()
figR0, axR0 = plt.subplots(1, 1, num=figR0.number)

figC0 = plt.figure(300)
figC0.clf()
figC0, axC0 = plt.subplots(1, 1, num=figC0.number)
axC0.lg0 = []
axC0.lg1 = []

ls = {'natural': '-', 'potential': ':'}
# zorder = {'baseline': 1, 'extraction': 0}
# alpha = {'baseline': 1, 'extraction': 0.5}
# width_factor = {'baseline': 1, 'extraction': 0.7}

alld = {}

for idx, region in enumerate(plot_regions):

    alld[region] = dreg = {}
    
    dreg['remote'] = pdh5.load(str(srcdir / '{case}/{region}.{type}-totals.h5'.format(case='baseline', region=region, type='remote')))
    dreg['natural'] = pdh5.load(str(srcdir / '{case}/{region}.{type}-totals.h5'.format(case='baseline', region=region, type='local')))
    dreg['potential'] = pdh5.load(str(srcdir / '{case}/{region}.{type}-totals.h5'.format(case='extraction', region=region, type='local')))

    fbins = dreg['remote']['fbins']
    
    specT = {}
    # Grab the edge of the EEZ, and take a time average, then convert to T-spec
    specT['remote'], specT['T'] = spec_freq2period(dreg['remote']['1way'][:, :, -1].mean(0),
                                                   fbins)
    # Sum over the EEZ, average in time, then convert to T-spec
    specT['natural'], _ = spec_freq2period(dreg['natural']['stot'].sum(-1).mean(0),
                                           fbins)
    specT['potential'], _ = spec_freq2period(dreg['potential']['stot'].sum(-1).mean(0),
                                             fbins)

    axR0.plot(specT['T'], specT['remote'] / int_period(specT['remote'], fbins),
              color=colors[region],
              label=labels[region],
              lw=3,
    )
    
    # axC0.plot(specT['T'], specT['remote'] / int_period(specT['remote'], fbins),
    #          color=colors[region],
    #          label=labels[region],
    # )
    drem = specT['remote']
    dpot = specT['remote'] + specT['potential']
    dnat = specT['remote'] + specT['natural']
    norm = int_period(dpot, fbins)
    drem /= int_period(drem, fbins)
    dpot /= norm
    axC0.lg0 += axC0.plot(specT['T'],
                          drem,
                          color=colors[region],
                          label=labels[region],
                          linestyle='-',
                          lw=3,
    )

    pvals = [0.05, 0.95]
    
    Tvals, Svals = calc_cpsd_points(drem, fbins, pvals)
    #pdb.set_trace()
    axC0.plot(Tvals, Svals, '.', color='w', ms=3, mew=2, zorder=5)
    axC0.plot(Tvals, Svals, '+', color=colors[region], ms=8, zorder=10, mew=1)
    
    axC0.plot(specT['T'],
              dpot,
              color=colors[region],
              lw=1,
              #linestyle='--',
              alpha=0.5,
    )
    Tvals, Svals = calc_cpsd_points(dpot, fbins, pvals)
    axC0.plot(Tvals, Svals, '.', color='w', ms=2, mew=1, zorder=2)
    axC0.plot(Tvals, Svals, '+', color=colors[region], ms=5, zorder=3, mew=0.7,
              alpha=0.5)
    
    # axC0.fill_between(specT['T'], y1=dnat, y2=dpot,
    #                   color=colors[region],
    #                   alpha=0.2,
    # )
    # axC0.fill_between(specT['T'], y1=drem, y2=dnat,
    #                   color=colors[region],
    #                   alpha=0.4,
    # )

    drem = specT['remote']
    dpot = specT['remote'] + specT['potential']
    dnat = specT['remote'] + specT['natural']
    norm = int_period(dpot, fbins)
    drem /= norm # Normalize by TOTAL POTENTIAL (less confusing!)
    dnat /= norm # Normalize by TOTAL POTENTIAL (less confusing!)
    dpot /= norm

    if flag['show each']:

        figCR = plt.figure(3000 + idx)
        figCR.clf()
        figCR, axCR = plt.subplots(1, 1, num=figCR.number)

        axCR.plot(specT['T'],
                  drem,
                  color='k',
                  label=labels[region],
                  linestyle='-'
        )
        axCR.plot(specT['T'],
                  dnat,
                  color='k',
                  linestyle='--'
        )
        axCR.plot(specT['T'],
                  dpot,
                  color='k',
                  linestyle='-.'
        )
        axCR.fill_between(specT['T'], y1=dnat, y2=dpot,
                          color='k',
                          alpha=0.2,
        )
        axCR.fill_between(specT['T'], y1=drem, y2=dnat,
                          color='k',
                          alpha=0.4,
        )
        axCR.set_title("{}: Total resource frequency breakdown".format(labels[region]))

    for case in ['natural', 'potential']:

        if case == 'natural':
            label = labels[region]
        else:
            label = None

        axL0.plot(specT['T'], specT[case] / int_period(specT[case], fbins),
                 color=colors[region],
                 label=label,
                 linestyle=ls[case]
        )

axR0.set_xlim([0, 30])
axR0.set_xticks(np.arange(0, 31, 5))
axR0.legend()
axR0.set_ylabel('Normalized Energy Distribution [1/s]')
axR0.set_title("Remote Resource")
axR0.set_xlabel("Wave Period [s]")
axR0.set_ylim([0, 0.2])
base.savefig(figR0, 'RemoteResource_Freq{}'.format(tag), dpi=300)

axL0.set_xlim([0, 30])
axL0.set_xticks(np.arange(0, 31, 5))
axL0.xaxis.grid(True)
axL0.axhline(0, linestyle=':', color='k')
axL0.set_title('Local Resource')
axL0.set_ylabel('Normalized Energy Distribution [1/s]')
axL0.set_xlabel('Wave Period [s]')
axL0.legend()
base.savefig(figL0, 'LocalResource_Freq{}'.format(tag), dpi=300)

axC0.lg1 += axC0.plot(np.NaN, np.NaN, linestyle='-', color='k', lw=3, label='Remote')
axC0.lg1 += axC0.plot(np.NaN, np.NaN, linestyle='-', color='k', lw=1, alpha=0.5, label='Remote + Potential')
axC0.set_ylim([0, 0.2])
axC0.set_title('Frequency Dependence of Wave Energy')
axC0.set_ylabel('Normalized Energy Distribution [1/s]')
axC0.legend(axC0.lg0, [h.get_label() for h in axC0.lg0])
figC0.legend(axC0.lg1, [h.get_label() for h in axC0.lg1], loc='upper center', bbox_to_anchor=[0.5, 0.88])
base.savefig(figC0, 'TotalResource_Freq{}'.format(tag))
