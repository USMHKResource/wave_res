import base as b
import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.ion()

wm = b.WrapMonths(shift=6)

if 'dat' not in vars():
    dat = {ky: b.TotResult(ky)
           for ky in
           ['wc', 'hi', 'ak', 'ec', 'gm', 'prusvi']}

######
# Plot the full REMOTE timeseries
fig0 = plt.figure(100)
fig0.clf()
ax0 = plt.gca()

plot_these = ['wc', 'hi', 'ak', 'ec', 'gm', ]

for ky in plot_these:
    dnow = dat[ky]
    tot = dnow.remote.total()['1way'][-1]
    
    ax0.plot(dnow.time, dnow.remote.int_freq()['1way'][:, -1] / tot, color=dnow.color, label=dnow.name)

ax0.legend()


######
# Plot the TOTAL annual cycle

fig1 = plt.figure(101)
fig1.clf()
ax1 = plt.gca()
ax = ax1
#plot_these = ['wc', 'hi', 'ak', 'ec', 'gm', 'prusvi']
plot_these = ['wc', 'hi', 'ak', 'ec', 'gm', ]

for ky in plot_these:
    dnow = dat[ky]
    dtmp = (dnow.remote.avg_annual()['1way'][:, -1] +
           dnow.local.avg_annual()['stot'].sum(1))

    dtmp /= dtmp.mean()
    ax.plot(wm.x, wm(dtmp), label=dnow.name, color=dnow.color)

ax.set_xticks(wm.xticks)
ax.set_xticklabels(wm.labels)
ax.set_xlim(wm.xlim)
ax.legend(loc='lower center', bbox_to_anchor=[0.5, 0.11], facecolor='w', framealpha=1)
#ax.legend(loc='upper right', bbox_to_anchor=[0.99, 0.99])
ax.axhline(1, linestyle='--', color='k', linewidth=2, zorder=-3)
ax.axhline(0.5, linestyle='--', color='0.6', linewidth=1, zorder=-3)
ax.axhline(1.5, linestyle='--', color='0.6', linewidth=1, zorder=-3)

for val in wm.season_edges:
    ax.axvline(val, color='0.6', linestyle=':',
                zorder=-5, linewidth=1)
ax.axvline(12, color='0.3', linestyle='--', zorder= -5, linewidth=1)

for val, lbl in zip(wm.season_centers, wm.season_labels):
    if lbl == 'summer':
        val += 0.3
    ax.text(val, 0.1, lbl, ha='center')

ax.axvspan(wm.season_edges[1], wm.season_edges[2], facecolor='b', zorder= -6, alpha=0.05)
ax.axvspan(wm.xlim[0], wm.season_edges[0], facecolor='r', zorder= -6, alpha=0.05)
ax.axvspan(wm.season_edges[-2], wm.xlim[-1], facecolor='r', zorder= -6, alpha=0.05)
    
ax.set_yticks(np.arange(0.0, 2.5, 0.5))
ax.set_ylim([0, 2])
ax.set_title("Annual cycle of total resource")

fig1.savefig('fig/AnnualCycle01.pdf')
fig1.savefig('fig/AnnualCycle01.png', dpi=200)


######
# Plot the REMOTE annual cycle

fig2 = plt.figure(102)
fig2.clf()
ax2 = plt.gca()
ax = ax2
plot_these = ['wc', 'hi', 'ak', 'ec', 'gm', ]

for ky in plot_these:
    dnow = dat[ky]
    dtmp = dnow.remote.avg_annual()['1way'][:, -1]
    dtmp /= dtmp.mean()

    ax.plot(wm.x, wm(dtmp), label=dnow.name, color=dnow.color)


ax.set_xticks(wm.xticks)
ax.set_xticklabels(wm.labels)
ax.set_xlim(wm.xlim)
ax.legend(loc='lower center', bbox_to_anchor=[0.5, 0.11])
ax.axhline(1, linestyle='--', color='k', linewidth=2, zorder=-3)
ax.axhline(0.5, linestyle='--', color='0.6', linewidth=1, zorder=-3)
ax.axhline(1.5, linestyle='--', color='0.6', linewidth=1, zorder=-3)

for val in wm.season_edges:
    ax.axvline(val, color='0.6', linestyle=':',
                zorder=-5, linewidth=1)
ax.axvline(12, color='0.3', linestyle='--', zorder= -5, linewidth=1)

for val, lbl in zip(wm.season_centers, wm.season_labels):
    if lbl == 'summer':
        val += 0.3
    ax.text(val, 0.1, lbl, ha='center')

ax.set_title("Annual cycle of remote resource")

ax.set_yticks(np.arange(0.0, 2.5, 0.5))
ax.set_ylim([0, 2])

ax.axvspan(wm.season_edges[1], wm.season_edges[2], facecolor='b', zorder= -6, alpha=0.05)
ax.axvspan(wm.xlim[0], wm.season_edges[0], facecolor='r', zorder= -6, alpha=0.05)
ax.axvspan(wm.season_edges[-2], wm.xlim[-1], facecolor='r', zorder= -6, alpha=0.05)

fig2.savefig('fig/AnnualCycle02.pdf')
fig2.savefig('fig/AnnualCycle02.png', dpi=200)


######
# Plot the inter-annual variability
fig3 = plt.figure(103)
fig0.clf()
ax3 = plt.gca()

ax = ax3

plot_these = ['wc', 'hi', 'ak', 'ec', 'gm', ]

t = np.arange(np.datetime64('1979'),np.datetime64('2011'))

for ky in plot_these:
    dnow = dat[ky]

    dtmp = (dnow.remote.avg_yearly()['1way'][:, -1] +
           dnow.local.avg_yearly()['stot'].sum(1))

    
    dtmp /= dtmp.mean()
    ax.plot(t, dtmp, label=dnow.name, color=dnow.color)

ax.set_ylim([0, 2])

