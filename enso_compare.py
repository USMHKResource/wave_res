import numpy as np
import base as b
import matplotlib.pyplot as plt
from enso_data import enso
import colormaps
reload(colormaps)
from colormaps import seasons

dnow = b.TotResult('wc')

if False:
    #####
    # Just using a simple annual average
    # Doesn't seem to work
    
    # According to:
    # https://psl.noaa.gov/enso/enso_101.html
    # >How long does ENSO last?
    #   El Niño typically lasts 9–12 months, and La Niña typically lasts
    #   1–3 years. Both tend to develop during March–June, reach peak
    #   intensity during December–April, and then weaken during
    #   May–July. However, prolonged El Niño episodes have lasted 2 years,
    #   and even as long as 3-4 years.
    # So, I'll do my average from May-April?

    dtmp = dnow.remote.int_freq()['1way'][:, -1]

    # May-April
    etmp = enso['oni'][4:-8].reshape((-1, 12))
    dtmp = dtmp[4:-8].reshape((-1, 12))
    t = dnow.time[4:-8]    

    fig = plt.figure(10)
    fig.clf()
    ax = plt.gca()

    ax.plot(t, etmp.flatten())
    ax.plot(t, dtmp.flatten() / 1000)

    fig = plt.figure(11)
    fig.clf()
    ax = plt.gca()

    ax.plot(etmp.mean(1), dtmp.mean(1), '.')

if True:
    # Get the full data
    L = dnow.remote._data['baseline']['length'][-1] * (365 * 24 * 1e-9)
    
    dtmp = dnow.remote.int_freq()['1way'][:, -1] / L
    davg = dnow.remote.avg_annual()['1way'][:, -1] / L

    # Calculate annomalies
    annom = dtmp.reshape((-1, 12)) - davg[None, :]

    # Smooth them a bit
    window = np.hanning(9)
    window = np.ones(5)
    # Normalize the window
    window /= window.sum()
    dplt = np.convolve(annom.flatten(), window, mode='same')

    t = dnow.time
    
    fig = plt.figure(110)
    fig.clf()
    ax = plt.gca()

    nlag = 2

    pdo = np.convolve(enso['pdo'], np.hanning(5), mode='same')

    ax.plot(t, enso['oni'], label='ONI')
    #ax.plot(t, enso['meiv2'], label='mei.v2')
    ax.plot(t, pdo, label='pdo')
    factor = 500
    plot_args = dict(color='0.5',
                     label='{} Resource'.format(dnow.name),
                     lw=2,
                     zorder=-5
    )
    if True:
        ax.plot(t, dplt / factor, **plot_args)
    else:
        # Lag the timeseries.
        if nlag > 0:
            ax.plot(t[:-nlag], dplt[nlag:] / factor, **plot_args)
        elif nlag == 0:
            ax.plot(t, dplt / factor, **plot_args)
        else:
            ax.plot(t[-nlag:], dplt[:nlag] / factor, **plot_args)
    ax.legend()
    #ax.set_title('{}'.format(dnow.name))

    fig = plt.figure(111)
    fig.clf()
    ax = plt.gca()

    eplt = enso['oni']

    if nlag > 0:
        eslc, dslc, = slice(None, -nlag), slice(nlag, None)
    elif nlag == 0:
        eslc, dslc = slice(None), slice(None)
    else:
        eslc, dslc = slice(-nlag, None), slice(None, nlag)

    cmap = seasons
    coff = 11
    
    
    sc = ax.scatter(eplt[eslc], dplt[dslc],
                    c=(dnow.time.month[dslc] - coff) % 12,
                    s=10,
                    cmap=cmap, vmin=0.5, vmax=12.5, )
    cbar = plt.colorbar(sc)
    cbar.set_ticks(np.arange(1, 13))
    wm = b.WrapMonths(coff)
    cbar.set_ticklabels(wm.labels)

    ax.axhline(0, color='0.6', lw=1, zorder=-10, ls='--')
    ax.axvline(0, color='0.6', lw=1, zorder=-10, ls='--')
    
        
    if nlag > 0:
        fit_func = np.array([eplt[:-nlag],
                             #eplt[:-nlag] ** 2,
                             np.exp(eplt[:-nlag]),
                             np.ones_like(eplt[:-nlag])]).T
        fit_vals = dplt[nlag:]
        

    bval, res, rank, s = np.linalg.lstsq(fit_func, fit_vals, rcond=None)

    x = np.arange(eplt.min(), eplt.max(), 0.1)
    #xfit = np.array([x, x ** 2, np.ones_like(x)]).T
    xfit = np.array([x, np.exp(x), np.ones_like(x)]).T
    y = (xfit * bval[None, :]).sum(1)
    
    ax.plot(x, y, 'k-', lw=3, zorder=-2)
    R2 = (1-res/(fit_vals.var()*len(fit_vals)))[0]
    ax.text(0.1, 0.9, '$R^2={:.2f}$'.format(R2),
            transform=ax.transAxes,
            color='k', size='x-large')
    

    ax.set_ylabel('Wave Energy Flux Anomaly (kW/m)')
    ax.set_xlabel(u'Oceanic Ni\xf1o Index')

    b.savefig(fig, 'ENSO-Comparison.{}'.format(dnow.region))
    
    fig = plt.figure(112)
    fig.clf()
    ax = plt.gca()

    nlag = 0

    eplt = pdo
    
    if nlag > 0:
        ax.plot(eplt[:-nlag], dplt[nlag:], '.')
    elif nlag == 0:
        ax.plot(eplt, dplt, '.')
    else:
        ax.plot(eplt[-nlag:], dplt[:nlag], '.')
        

    nlag = 2

    if nlag > 0:
        fit_func = np.array([enso['oni'][:-nlag],
                             #pdo[nlag:],
                             np.exp(enso['oni'][:-nlag]),
                             #enso['oni'][:-nlag] ** 2,
                             #np.exp(1j * np.pi / 6 * np.array(dnow.time.month[:-nlag])), 
                             np.ones_like(pdo[nlag:])]).T
        fit_vals = dplt[nlag:]
        

    bval, res, rank, s = np.linalg.lstsq(fit_func, fit_vals, rcond=None)
    
    fig = plt.figure(113)
    fig.clf()
    ax = plt.gca()

    fit = (fit_func*bval[None,:]).sum(1)

    ax.plot(fit, fit_vals, '.')
    R2 = (1-res/(fit_vals.var()*len(dplt)))[0]
    ax.text(0.1, 0.9, '$R^2={:.3f}$'.format(R2), transform=ax.transAxes)
