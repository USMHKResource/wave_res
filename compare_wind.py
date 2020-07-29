import zipfile
import numpy
import pdb
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import base

wm = base.WrapMonths(shift=6)

def load_zip_csv(zipname, csvname, **kwargs):
    """Read a zipped csv.
    
    **kwargs are inputs to pd.read_csv
    """
    with zipfile.ZipFile(zipname) as zf:
        with zf.open(csvname) as fl:
            dat = pd.read_csv(fl, **kwargs)
    return dat


gray = load_zip_csv('results/powerWindGrayland.zip', 'powerWindGrayland.csv', skiprows=2)
mdco = load_zip_csv('results/powerWindMendocino.zip', 'powerWindMendocino.csv', skiprows=2)


gray['U'] = np.abs(gray['u'] + 1j * gray['v'])
mdco['U'] = np.abs(mdco['u'] + 1j * mdco['v'])

g = gray.groupby('Month').mean()
m = mdco.groupby('Month').mean()

if False:

    wcolor = 'b'
    jcolor = 'r'

    fig1 = plt.figure(101)
    fig1.clf()
    fig1, axs = plt.subplots(2, 1, num=fig1.number, sharex=True)
    axWm = axs[0]
    axJm = plt.twinx(axWm)
    axWg = axs[1]
    axJg = plt.twinx(axWg)

    axWm.plot(wm.x, wm(m['U']), color=wcolor)
    axJm.plot(wm.x, wm(m['j']) / 1000, color=jcolor)
    axWg.plot(wm.x, wm(g['U']), lw=2, color=wcolor)
    axJg.plot(wm.x, wm(g['j']) / 1000, lw=2, color=jcolor)

    axWm.set_ylabel('Wind Speed [m/s]')
    axWg.set_ylabel('Wind Speed [m/s]')
    axJm.set_ylabel('J [kW/m]')
    axJg.set_ylabel('J [kW/m]')

    axJm.set_ylim([0, 80])
    axJg.set_ylim([0, 80])
    axWm.set_ylim([4, 8])
    axWg.set_ylim([4, 8])

    for ax in [axJm, axJg]:
        ax.spines['right'].set_color(jcolor)
        ax.tick_params(axis='y', colors=jcolor)
        ax.yaxis.label.set_color(jcolor)
        ax.spines['left'].set_visible(False)

    for ax in [axWm, axWg]:
        ax.spines['left'].set_color(wcolor)
        ax.spines['right'].set_visible(False)
        ax.tick_params(axis='y', colors=wcolor)
        ax.yaxis.label.set_color(wcolor)

    ax.set_xticks(wm.xticks)
    ax.set_xticklabels(wm.labels)
    ax.set_xlim(wm.xlim)


if True:

    wcolor = 'g'
    jcolor = 'b'

    g_args = dict(lw=2, alpha=0.7, ls='--')
    m_args = dict(lw=3, alpha=1, ls='-')
    
    fig2 = plt.figure(102)
    fig2.clf()
    fig2, axJ = plt.subplots(1, 1, num=fig2.number)
    axW = plt.twinx(axW)

    axW.plot(wm.x, wm(m['U']), color=wcolor, **m_args)
    axJ.plot(wm.x, wm(m['j']) / 1000, color=jcolor, **m_args)
    axW.plot(wm.x, wm(g['U']), color=wcolor, **g_args)
    axJ.plot(wm.x, wm(g['j']) / 1000, color=jcolor, **g_args)

    axJ.plot(np.NaN, np.NaN, color='k', label='Cape Mendocino', **m_args)
    axJ.plot(np.NaN, np.NaN, color='k', label='Northern Washington', **g_args)
    
    axW.set_ylabel('Wind Speed [m/s]')
    axW.set_ylabel('Wind Speed [m/s]')
    axJ.set_ylabel('J [kW/m]')
    axJ.set_ylabel('J [kW/m]')

    axJ.set_ylim([0, 80])
    axJ.set_ylim([0, 80])
    axW.set_ylim([4, 8])
    axW.set_ylim([4, 8])

    axJ.spines['left'].set_color(jcolor)
    axJ.spines['right'].set_visible(False)
    axJ.tick_params(axis='y', colors=jcolor)
    axJ.yaxis.label.set_color(jcolor)

    axW.spines['right'].set_color(wcolor)
    axW.spines['left'].set_visible(False)
    axW.tick_params(axis='y', colors=wcolor)
    axW.yaxis.label.set_color(wcolor)

    axJ.set_xticks(wm.xticks)
    axJ.set_xticklabels(wm.labels)
    axW.set_xlim(wm.xlim)
    axJ.set_xlim(wm.xlim)

    ax = axJ
    ax.axvspan(wm.season_edges[1], wm.season_edges[2], facecolor='b', zorder= -6, alpha=0.05)
    ax.axvspan(wm.xlim[0], wm.season_edges[0], facecolor='r', zorder= -6, alpha=0.05)
    ax.axvspan(wm.season_edges[-2], wm.xlim[-1], facecolor='r', zorder= -6, alpha=0.05)

    axJ.legend(loc='center', bbox_to_anchor=[0.5, 0.25], prop=dict(size='large'))

    fig2.savefig('fig/CompareWind01.pdf')
    fig2.savefig('fig/CompareWind01.png')
