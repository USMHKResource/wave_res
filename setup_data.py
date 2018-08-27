import calc_remote as cr
import base
import numpy as np
import paths as p
import cPickle as pkl


def create_info():
    con_defs = {}
    freqbins = {}
    for region in base.regions:

        dat = cr.load(region, '2009-01')
        v = dat.variables
        con_defs[region] = {}
        cid = cr._concatenate_id(v['station_name'][:, 2:5].data)
        freqbins[region] = np.hstack((v['frequency1'][:].data,
                                      [v['frequency2'][-1].data]))
        dat.close()
        for id in base.conids:
            tmp = np.nonzero(cid == id)[0][[0, -1]]
            con_defs[region][id] = slice(tmp[0], tmp[1] + 1)

    with open(str(p.projdir / 'data/Contour_Ranges.pkl'), 'w') as fl:
        pkl.dump(con_defs, fl)
    with open(str(p.projdir / 'data/FreqBins.pkl'), 'w') as fl:
        pkl.dump(freqbins, fl)


if __name__ == '__main__':
    create_info()
