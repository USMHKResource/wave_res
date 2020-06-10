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


def spec_freq2period(spec, fbins, norm=True):
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
    
if flag['show remote']:

    fig0 = plt.figure(102)
    fig0.clf()
    fig0, ax0 = plt.subplots(1, 1, num=fig0.number)

    for idx, region in enumerate(plot_regions):

        dat = pdh5.load(str(srcdir / '{case}/{region}.{type}-totals.h5'.format(case='baseline', region=region, type='remote')))

        # Grab the edge of the EEZ, and take a time average
        dnow = dat['1way'][:, :, -1].mean(0)
        specT, T_center = spec_freq2period(dnow, dat['fbins'])
                
        ax0.plot(T_center, specT,
                 color=colors[region],
                 label=labels[region])
        
        # ax0.bar(Tbins[1:] - bar_plot_coefs[2] - idx * bar_plot_coefs[1],
        #         np.diff(Int2),
        #         color=colors[region],
        #         width=bar_plot_coefs[0],
        #         label=labels[region])

    ax0.set_xlim([0, 30])
    ax0.set_xticks(np.arange(0, 31, 5))
    #ax0.xaxis.grid(True)
    ax0.legend()
    ax0.set_ylabel('Normalized Energy Distribution [1/s]')
    ax0.set_title("Remote Resource")
    ax0.set_xlabel("Wave Period [s]")
    ax0.set_ylim([0, 0.2])
    fig0.savefig('../fig/RemoteResource_Freq{}.pdf'.format(tag))
    fig0.savefig('../fig/RemoteResource_Freq{}.png'.format(tag), dpi=300)


ls = {'baseline': '-', 'extraction': ':'}
# zorder = {'baseline': 1, 'extraction': 0}
# alpha = {'baseline': 1, 'extraction': 0.5}
# width_factor = {'baseline': 1, 'extraction': 0.7}

if flag['show local']:

    fig0 = plt.figure(202)
    fig0.clf()
    fig0, ax0 = plt.subplots(1, 1, num=fig0.number)

    for idx, region in enumerate(plot_regions):

        for case in ['baseline', 'extraction']:
            dat = pdh5.load(str(srcdir / '{case}/{region}.{type}-totals.h5'.format(case=case, region=region, type='local')))


            # Sum over the EEZ, average in time
            dnow = dat['stot'][:, :, :].sum(-1).mean(0)

            specT, T_center = spec_freq2period(dnow, dat['fbins'])
            
            if case == 'baseline':
                label = labels[region]
            else:
                label = None

            ax0.plot(Tbins_center, specT,
                     color=colors[region],
                     label=label,
                     linestyle=ls[case]
            )

            # axI.plot(1. / f, Int)
            # Int2 = np.interp(Fbins, f, Int)
            # axI.plot(1. / Fbins, Int2, '+')

            # ax0.bar(Tbins[1:] - bar_plot_coefs[2] - idx * bar_plot_coefs[1],
            #         np.diff(Int2),
            #         color=colors[region],
            #         width=bar_plot_coefs[0] * width_factor[case],
            #         label=label,
            #         hatch=hatch[case], zorder=zorder[case],
            #         alpha=alpha[case])

    ax0.set_xlim([0, 30])
    ax0.set_xticks(np.arange(0, 31, 5))
    ax0.xaxis.grid(True)
    ax0.axhline(0, linestyle=':', color='k')
    ax0.set_title('Local Resource')
    ax0.set_ylabel('Normalized Energy Distribution [1/s]')
    ax0.set_xlabel('Wave Period [s]')
    ax0.legend()
    #ax0.set_ylim(local_ylim)
    fig0.savefig('../fig/LocalResource_Freq{}.pdf'.format(tag))
    fig0.savefig('../fig/LocalResource_Freq{}.png'.format(tag), dpi=300)
