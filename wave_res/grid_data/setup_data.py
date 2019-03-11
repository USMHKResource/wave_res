import setpath
import numpy as np
import cPickle as pkl
from copy import deepcopy
import gzip
# Import wave_res modules
import calc_remote as cr
import gis
import base


def create_info():
    """
    This function extracts the grid data from the directional wave
    spectra data files. It extracts: lat/lon, frequency, and defines
    the contour_range indices.
    """
    con_defs = {}
    freqbins = {}
    gridlonlat = {}
    for region in base.source_regions:
        print("Processing data for region '{}'...".format(region))
        dat = cr.load('baseline', region, '2009-01')
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

        if region in ['ak', 'hi']:
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

    with open(str('Contour_Ranges.pkl'), 'w') as fl:
        pkl.dump(con_defs, fl)
    with open(str('FreqBins.pkl'), 'w') as fl:
        pkl.dump(freqbins, fl)
    with open(str('GridLonLat.pkl'), 'w') as fl:
        pkl.dump(gridlonlat, fl)


def check_loop(inf, s, r_max):
    r_ends = np.abs(gis.diffll(
        inf.gridlonlat[:, [s[0], s[-1]]])[0])
    if r_ends < r_max:
        s = s + [s[0]]
    return s


def fix_gaps():
    """This function fixes the contour definitions (indices) to insert
    breaks where a contour should not continue, and creates loops
    where they should be.

    This is done by looking at the distance between points.
    - When the distance between adjacent points on a contour is
      greater than r_max, than we assume it is a new contour.
    - When the distance between endpoints of a contour are closer
      together than 2*r_max, than we assume it is a loop.
    """
    r_max = 40000

    new_con_defs = {}
    for reg in base.regions:
        rinf = base.RegionInfo(reg, use_old_con_defs=True)
        cdefs = deepcopy(rinf.con_defs)
        taile = []
        for id in np.sort(rinf.con_defs.keys()):
            if id == 'eez':
                # These are created explicitly, and might have breaks already inserted.
                continue
            for ll in rinf.get_contour(id):
                # Calculate the distance betweeen points
                r = np.abs(gis.diffll(ll)[0])
                # Find the distances greater than r_max
                inds = np.nonzero(r > r_max)[0]
                if len(inds) == 0 and isinstance(cdefs[id], slice):
                    cdefs[id] = [np.mgrid[cdefs[id]].tolist(), ]
                else:
                    tmp = []
                    last = start = cdefs[id].start
                    for i in inds:
                        tmp.append(check_loop(rinf,
                                              range(last, start + i + 1),
                                              2 * r_max))
                        last = start + i + 1
                    tmp.append(
                        check_loop(rinf,
                                   range(last, cdefs[id].stop),
                                   2 * r_max))
                    cdefs[id] = tmp
            poly_inds = []
            if reg == 'hi':
                if id == 'EEZ':
                    rng = 200
                else:
                    rng = int(id)
                seg_i = 0
                poly_inds.append(cdefs[id][seg_i])
                if rng < 51:
                    pass
                elif rng == 60:
                    seg_i += 1
                    poly_inds[0] = poly_inds[0][:-7]
                    poly_inds.append(cdefs[id][seg_i][13:])
                elif rng == 70:
                    seg_i += 1
                    poly_inds[0] = poly_inds[0][:-7]
                    poly_inds.append(cdefs[id][seg_i][15:])
                elif rng == 80:
                    poly_inds[0] = poly_inds[0][15:-5]
                elif rng == 90:
                    poly_inds[0] = poly_inds[0][17:-3]
                elif rng == 100:
                    poly_inds[0] = poly_inds[0][17:-2]
                elif rng == 200:
                    poly_inds[0] = poly_inds[0][:50] + taile + poly_inds[0][-113:]
                elif rng >= 110:
                    poly_inds[0] = poly_inds[0][17:]
            
                # if rng > 90:
                #     tail0 += [poly_inds[0][0]]
                #     poly_inds[0] = tail0_last + poly_inds[0]
                if (100 < rng) & (rng < 200):
                    taile += [poly_inds[0][-1], ]
                    poly_inds[0] = poly_inds[0] + taile_last[::-1]
                taile_last = taile

                if rng >= 40:
                    # Patch the start back on
                    poly_inds = poly_inds + [poly_inds[0][0], ]
                    cdefs[id] = [np.hstack(poly_inds).tolist(), ] + cdefs[id][seg_i + 1:]

        new_con_defs[reg] = cdefs

    with open(str('Contour_Ranges2.pkl'), 'w') as fl:
        pkl.dump(new_con_defs, fl)


if __name__ == '__main__':
    # Create Grid data
    create_info()
    # Insert breaks and loops
    fix_gaps()

    # This extracts the land data that is relavent for each data point
    import extract_land
    land_data = extract_land.run_all()
    with open(str('LandData.pkl'), 'w') as fl:
        pkl.dump(land_data, fl)

    # This sets boundaries (polygons) of each contour.
    import boundaries
    bounds = boundaries.run_all()
    with open(str('Boundaries.pkl'), 'w') as fl:
        pkl.dump(bounds, fl)

    # This finds the triangles inside of each border
    import create_triangles as tri
    tri_defs = tri.run_all()
    # # This stores the full Triangle.pkl.gz file.
    # with gzip.open(str('Triangles.pkl.gz'), 'w') as fl:
    #     pkl.dump(tri_defs, fl)
    tri_defs_diff = tri.run_diff_tri_dict(tri_defs)
    with gzip.open(str('DiffTriangles.pkl.gz'), 'w') as fl:
        pkl.dump(tri_defs_diff, fl)
