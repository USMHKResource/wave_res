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
          'ak': 'Alaska'}

plot_regions = ['wc', 'hi', 'ec']
# [width, spacing, offset]
bar_plot_coefs = [1.2, 1.6666, 0.8]
# plot_regions = ['wc', 'hi', 'ec', 'ak']
# # [width, spacing, offset]
# bar_plot_coefs = [0.8, 1, 0.8]

if flag['show remote']:

    figI = plt.figure(101)
    figI.clf()
    figI, axI = plt.subplots(1, 1, num=figI.number)

    fig0 = plt.figure(102)
    fig0.clf()
    fig0, ax0 = plt.subplots(1, 1, num=fig0.number)

    for idx, region in enumerate(plot_regions):

        dat = pdh5.load(str(srcdir / '{case}/{region}.{type}-totals.h5'.format(case='baseline', region=region, type='remote')))

        df = np.diff(dat['fbins'])
        f = dat['fbins'][:-1] + df / 2

        Tbins = np.arange(0, 31, 5)
        Fbins = 1. / Tbins


        dnow = dat['1way'][:,:,-1].mean(0) * df / 1e9
        #dnow = np.hstack(([0], dnow, [0]))
        #f2 = np.hstack(([1. / 30], f, [0.5]))

        Int = np.cumsum(dnow[::-1])[::-1]

        axI.plot(1. / f, Int)
        Int2 = np.interp(Fbins, f, Int)
        axI.plot(1. / Fbins, Int2, '+')

        ax0.bar(Tbins[1:] - bar_plot_coefs[2] - idx * bar_plot_coefs[1],
                np.diff(Int2),
                color=colors[region],
                width=bar_plot_coefs[0],
                label=labels[region])

    ax0.set_xlim([0, 30])
    ax0.set_xticks(np.arange(0, 31, 5))
    ax0.xaxis.grid(True)
    ax0.legend()
    ax0.set_ylabel('[GW]')
    ax0.set_title("Remote Resource")
    ax0.set_xlabel("Wave Period [s]")
    ax0.set_ylim([0, 25])
    fig0.savefig('../fig/RemoteResource_Freq01.pdf')
    fig0.savefig('../fig/RemoteResource_Freq01.png', dpi=300)


hatch = {'baseline': None, 'extraction': None}
zorder = {'baseline': 1, 'extraction': 0}
alpha = {'baseline': 1, 'extraction': 0.5}
width_factor = {'baseline': 1, 'extraction': 0.7}

if flag['show local']:

    figI = plt.figure(201)
    figI.clf()
    figI, axI = plt.subplots(1, 1, num=figI.number)

    fig0 = plt.figure(202)
    fig0.clf()
    fig0, ax0 = plt.subplots(1, 1, num=fig0.number)

    for idx, region in enumerate(plot_regions):

        for case in ['baseline', 'extraction']:
            dat = pdh5.load(str(srcdir / '{case}/{region}.{type}-totals.h5'.format(case=case, region=region, type='local')))

            df = np.diff(dat['fbins'])
            dnow = (dat.sds + dat.snl + dat.sin).sum(-1)
            dnow *= dat.Nhour[:, None]
            dnow = dnow.sum(0) / dat.Nhour.sum() / 1e9
            dnow *= df

            Int = np.cumsum(dnow[::-1])[::-1]

            axI.plot(1. / f, Int)
            Int2 = np.interp(Fbins, f, Int)
            axI.plot(1. / Fbins, Int2, '+')

            if case == 'baseline':
                label = labels[region]
            else:
                label = None
            ax0.bar(Tbins[1:] - bar_plot_coefs[2] - idx * bar_plot_coefs[1],
                    np.diff(Int2),
                    color=colors[region],
                    width=bar_plot_coefs[0] * width_factor[case],
                    label=label,
                    hatch=hatch[case], zorder=zorder[case],
                    alpha=alpha[case])

    ax0.set_xlim([0, 30])
    ax0.set_xticks(np.arange(0, 31, 5))
    ax0.xaxis.grid(True)
    ax0.axhline(0, linestyle=':', color='k')
    ax0.set_title('Local Resource')
    ax0.set_ylabel('[GW]')
    ax0.set_xlabel('Wave Period [s]')
    ax0.legend()
    ax0.set_ylim([-5, 20])
    fig0.savefig('../fig/LocalResource_Freq01.pdf')
    fig0.savefig('../fig/LocalResource_Freq01.png', dpi=300)
