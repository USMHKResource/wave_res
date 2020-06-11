import pyDictH5 as pdh5
import pathlib2 as pl
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumtrapz
from collections import defaultdict
plt.ion()

flag = defaultdict(lambda: False)

flag['show remote'] = True
flag['show local'] = True

srcdir = pl.Path('../results/freq.fcut/').resolve()

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

# Number of plotted regions
# [width, spacing, offset]
bar_plot_coefs = {3: [1.2, 1.6666, 0.8],
                  4: [0.8, 1, 0.8],
                  5: [0.7, 0.8, 0.7]}[len(plot_regions)]

if 'ak' in plot_regions:
    local_ylim = np.array([-5, 100])
    remote_ylim = [0, 60]
else:
    remote_ylim = [0, 25]
    local_ylim = [-5, 20]
# plot_regions = ['wc', 'hi', 'ec', 'ak']
# # [width, spacing, offset]


def spec_freq2period(spec, fbins, norm=False):
    """Convert the units of a frequency spectrum (W/Hz) to a
    period-spectrum (W/s).
    """
    
    Tbins = fbins ** -1
    df = np.diff(fbins)
    f = fbins[:-1] + df / 2
    specT = spec * (f ** 2)
    dT = -np.diff(Tbins)
    T = f ** -1
    
    # Check that the integrals match!
    If = np.sum(spec * df)
    IT = np.sum(specT * dT)
    assert (IT - If) / If < 0.01, "The integrals don't match!?"

    if norm:
        specT /= IT
    return specT, T


def int_period(specT, fbins, axis=0):
    Tbins = fbins ** -1
    dT = -np.diff(Tbins)
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
    )

    axC0.plot(specT['T'], specT['remote'] / int_period(specT['remote'], fbins),
             color=colors[region],
             label=labels[region],
    )
    
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
figR0.savefig('../fig/RemoteResource_Freq{}.pdf'.format(tag))
figR0.savefig('../fig/RemoteResource_Freq{}.png'.format(tag), dpi=300)

axL0.set_xlim([0, 30])
axL0.set_xticks(np.arange(0, 31, 5))
axL0.xaxis.grid(True)
axL0.axhline(0, linestyle=':', color='k')
axL0.set_title('Local Resource')
axL0.set_ylabel('Normalized Energy Distribution [1/s]')
axL0.set_xlabel('Wave Period [s]')
axL0.legend()
figL0.savefig('../fig/LocalResource_Freq{}.pdf'.format(tag))
figL0.savefig('../fig/LocalResource_Freq{}.png'.format(tag), dpi=300)
