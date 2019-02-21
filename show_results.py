import pyDictH5 as pdh5
import numpy as np
import matplotlib.pyplot as plt

regions = {'ak', 'at', 'prusvi', 'wc'}

remote = {}
local = {}
rtot = {}
ltot = {}

unit = 'GW'
unit = 'TWh/yr'

for ireg, region in enumerate(regions):
    rd = remote[region] = pdh5.load('tmpdata/results/baseline/{}.remote-totals.h5'.format(region))
    ld = local[region] = pdh5.load('tmpdata/results/baseline/{}.local-totals.h5'.format(region))

    if unit == 'GW':
        factor = 1e-9
    elif unit == 'TWh/yr':
        factor = 365 * 24 * 1e-12
    else:
        raise Exception("Invalid unit specification.")
    rd['oneway'] = rd['1way']
        
    rtot[region] = {m: (np.average(rd[m][:, -1], weights=rd['Nhour']) * factor)
                    for m in ['oneway', 'unit', 'bdir', 'trad']}
    ltot[region] = {m: (np.average(ld[m][:, 2:-2].sum(-1), weights=ld['Nhour']) * factor)
                    for m in ['sbt', 'sds', 'snl', 'stot', 'sin', 'sice']}

print("")

print("Remote Totals ({})".format(unit))
print("============")
print(("{:10s}|" + "{:>8s}" * 4).format("", 'one-way', 'trad', 'unit', 'bidir'))
print("-" * 44)
for ireg, region in enumerate(regions):
    rt = rtot[region]
    print("{:10s}: {oneway: 7.4g} {trad: 7.4g} {unit: 7.4g} {bdir: 7.4g}"
          .format(region, **rt))


print("Local Totals ({})".format(unit))
print("============")
print(("{:10s}|" + "{:>8s}" * 6).format("", 'stot', 'sin', 'sds', 'snl', 'sice', 'sbt'))
print("-" * 44)
for ireg, region in enumerate(regions):
    lt = ltot[region]
    print("{:10s}: {stot: 7.4g} {sin: 7.4g} {sds: 7.4g} {snl: 7.4g} {sice: 7.4g} {sbt: 7.4g}"
          .format(region, **lt))


fctr = 1e-9
fignum = 10
fig = plt.figure(fignum)
fig.clf()
fig, ax = plt.subplots(1, 1, num=fignum)
t = local['ak'].time
dnow = local['ak'].stot.sum(-1) * fctr
pos = dnow > 0
ax.plot(t[~pos], -dnow[~pos], '_', color='r', label='-AK')
ax.semilogy(t[pos], dnow[pos], 'g+', label='AK')
dnow = local['wc'].stot.sum(-1) * fctr
ax.semilogy(t, dnow, 'b.', label='WC')
ax.set_ylabel('Local Wave Resoruce (GW)')
ax.legend()
ax.set_yticks(np.logspace(0,11,12))
ax.yaxis.grid(True)
ax.set_xlabel("Time")
fig.savefig('fig/AK-Resource.time.png')

fctr = 1e-9
fignum = 12
fig = plt.figure(fignum)
fig.clf()
fig, ax = plt.subplots(1, 1, num=fignum)
r = local['ak'].range
dnow = local['ak'].stot.mean(0) * fctr
pos = dnow > 0
ax.plot(r[~pos], -dnow[~pos], '_', color='r', label='-AK')
ax.semilogy(r[pos], dnow[pos], 'g+', label='AK')
dnow = local['wc'].stot.mean(0) * fctr
ax.semilogy(r, dnow, 'b.', label='WC')
ax.set_ylabel('Local Wave Resoruce (GW)')
ax.legend()
ax.set_yticks(np.logspace(0,11,12))
ax.yaxis.grid(True)
ax.set_xlabel("Distance from shore")
fig.savefig("fig/AK-Resource.distance.png")
