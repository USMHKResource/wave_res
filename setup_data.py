import calc_remote as cr
import base
import numpy as np
import paths as p
import cPickle as pkl
import gis
from copy import deepcopy


def create_info():
    con_defs = {}
    freqbins = {}
    gridlonlat = {}
    for region in base.source_regions:
        print("Processing data for region '{}'...".format(region))
        dat = cr.load(region, '2009-01')
        v = dat.variables
        con_defs[region] = {}
        cid = cr._concatenate_id(v['station_name'][:, 2:5].data)
        freqbins[region] = np.hstack((v['frequency1'][:].data,
                                      [v['frequency2'][-1].data]))
        # The first dim in lon/lat is time, and it doesn't change.
        ll = gridlonlat[region] = np.vstack((v['longitude'][0, :].data,
                                             v['latitude'][0, :].data))
        dat.close()
        for id in base.conids:
            tmp = np.nonzero(cid == id)[0][[0, -1]]
            con_defs[region][id] = slice(tmp[0], tmp[1] + 1)

        if region == 'ak':
            ll[0][ll[0] > 0] -= 360

        # Split the 'at' grid into 'ec' and 'gm' subgrids
        if region == 'at':
            con_defs['ec'] = {}
            con_defs['gm'] = {}
            gm_inds = ((ll[1] < 29) & (ll[0] < -80.9) |
                       (ll[0] < -82))
            for id in base.conids:
                slc = con_defs['at'][id]
                n = gm_inds[slc].sum()
                sep = slc.stop - n
                delta = 0
                if id in ['010', '020', 'EEZ']:
                    # Subtract 1 so that the two regions connect
                    delta = -1
                con_defs['gm'][id] = slice(sep + delta, slc.stop)
                con_defs['ec'][id] = slice(slc.start, sep)

    with open(str(p.projdir / 'data/Contour_Ranges.pkl'), 'w') as fl:
        pkl.dump(con_defs, fl)
    with open(str(p.projdir / 'data/FreqBins.pkl'), 'w') as fl:
        pkl.dump(freqbins, fl)
    with open(str(p.projdir / 'data/GridLonLat.pkl'), 'w') as fl:
        pkl.dump(gridlonlat, fl)


def fix_gaps():

    r_max = 40000

    new_con_defs = {}
    for reg in base.regions:
        inf = base.RegionInfo(reg, use_old_con_defs=True)
        cdefs = deepcopy(inf.con_defs)
        for id in inf.con_defs.keys():
            if id == 'eez':
                # These are created explicitly, and might have breaks already inserted.
                continue
            for ll in inf.get_lonlat(id):
                # Calculate the distance betweeen points
                r = np.abs(gis.diffll(ll)[0])
                # Find the distances greater than r_max
                inds = np.nonzero(r > r_max)[0]
                if len(inds) == 0 and isinstance(cdefs[id], slice):
                    cdefs[id] = [cdefs[id], ]
                else:
                    tmp = []
                    last = start = cdefs[id].start
                    for i in inds:
                        r_ends = np.abs(gis.diffll(
                            inf.gridlonlat[:, [last, start + i]])[0])
                        if r_ends < r_max:
                            step = 1j
                        else:
                            step = 1
                        slc_now = slice(last, start + i + 1, step)
                        tmp.append(slc_now)
                        last = start + i + 1
                    tmp.append(slice(last, cdefs[id].stop, step))
                    cdefs[id] = tmp
        new_con_defs[reg] = cdefs
    with open(str(p.projdir / 'data/Contour_Ranges2.pkl'), 'w') as fl:
        pkl.dump(new_con_defs, fl)


if __name__ == '__main__':
    create_info()
