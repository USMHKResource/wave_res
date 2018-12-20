"""
This script grabs the coastline data from the NaturalEarth dataset that overlap with the model output, and writes it to a data file.
"""
from __future__ import division 
import matplotlib.pyplot as plt
import shapely.geometry as sg
import numpy as np
import setpath
import proj
from base import RegionInfo


land = proj.land

def plot_lines(rinf, ax, **kwargs):
    eez_kw = kwargs.pop('eez_kw', {})
    kwargs.update(transform=proj.pc)
    for cid in rinf.con_defs:
        tmp = kwargs.copy()
        if cid == 'eez':
            tmp.update(eez_kw)
        for ll in rinf.get_contour(cid):
            ax.plot(ll[0], ll[1], **tmp)


def overlap(r0, r1):
    return (within(r0[0], r1) or within(r0[1], r1) or
            within(r1[0], r0) or within(r1[1], r0))


def within(val, r):
    return (r[0] < val) & (val < r[1])


def argmaxd(d):
    return max(d.iterkeys(), key=lambda key: d[key])


def argmind(d):
    return min(d.iterkeys(), key=lambda key: d[key])


def argsortd(d, reverse=False):
    return sorted(d, key=d.get, reverse=reverse)


def argnearest(pt, line):
    d2 = ((pt[:, None] - line) ** 2).sum(0)
    return np.argmin(d2)


def get_coastline_data(lonlim, latlim, max_area=5800, plot_kwargs=None):
    cons = {}
    areas = {}

    for idx, d in enumerate(land.geometries()):
        b = d.bounds
        lonb = np.array([b[0], b[2]])
        latb = [b[1], b[3]]
        if (overlap(lonb, lonlim) or overlap(lonb - 360, lonlim) and
                overlap(latb, latlim)):
            for id2, p in enumerate(d):
                try:
                    lldat = np.stack(p.boundary.xy)
                except:
                    continue
                if d.area > max_area:
                    continue
                areas[idx] = d.area
                cons[idx] = lldat
                if plot_kwargs is not None:
                    ax.plot(lldat[0], lldat[1], **plot_kwargs)

    area_order = argsortd(areas, reverse=True)
    coastline = []
    for k in area_order:
        coastline.append(cons[k])
    return coastline


def get_mainland(lleez, ll_land, prj):
    # Grab the segment of mainland that is closest to the endpoints of the eez.
    xy0 = np.array(prj.transform_point(lleez[0, 0], lleez[1, 0], proj.pc))
    xye = np.array(prj.transform_point(lleez[0, -1], lleez[1, -1], proj.pc))
    xy_land = prj.transform_points(proj.pc, ll_land[0], ll_land[1]).T[:2]
    d0 = np.sqrt((xy0[0] - xy_land[0]) ** 2 + (xy0[1] - xy_land[1]) ** 2) / 1000
    de = np.sqrt((xye[0] - xy_land[0]) ** 2 + (xye[1] - xy_land[1]) ** 2) / 1000
    ie = np.argmin(de)
    i0 = np.argmin(d0)
    if i0 < ie:
        step = -1
    else:
        step = 1
    slc = slice(ie, i0 + step, step)
    return np.array(ll_land[:, slc], dtype=np.float32)


def get_islands(boundary_ll, land_list, prj):
    # eez is the lonlat eez data
    # mainland is the lonlat mainland data
    # land_list is a list of lonlat land data
    # prj is the projection
    boundary_xy = prj.transform_points(proj.pc,
                                       boundary_ll[0],
                                       boundary_ll[1]).T[:2]
    boundary_poly = sg.Polygon(boundary_xy.T)
    # Grab the islands that are inside the EEZ.
    islands = []
    if isinstance(land_list, dict):
        land_list = list(land_list.itervalues())
    for land in land_list:
        xy = prj.transform_points(proj.pc, land[0], land[1])[:, :2]
        pnow = sg.Polygon(xy)
        if boundary_poly.contains(pnow):
            islands.append(land)
    return islands

def setup_figure(rinf, fignum, **kwargs):
    region = rinf.region
    # First setup the figure
    prj = proj.proj[region]
    fig = plt.figure(fignum)
    fig.clf()
    ax = plt.axes(projection=prj)
    feat = ax.add_feature(land)
    ax.add_feature(proj.states)
    ax.set_extent((prj.lonlim + prj.latlim), crs=proj.pc)
    plot_lines(rinf, ax, **kwargs)
    return fig, ax


def get_land_wc():

    rinf = RegionInfo('wc', use_old_con_defs=False)
    prj = proj.proj[rinf.region]

    fig, ax = setup_figure(rinf, 31,
                           color='m',
                           eez_kw=dict(lw=3), alpha=0.5)

    # First grab the EEZ boundary (including borders)
    lleez = rinf.get_contour('EEZ')[0]

    # Now grab the coastline data
    coastline = get_coastline_data(prj.lonlim, prj.latlim, )

    # Now pick-out the segment of mainland that matches it
    mainland = get_mainland(lleez, coastline[0], prj)

    # Construct the boundary
    boundary_ll = np.hstack((lleez, mainland))
    # ... and plot it
    ax.plot(boundary_ll[0], boundary_ll[1], 'r-', transform=proj.pc)

    # Now grab and plot the islands
    islands = get_islands(boundary_ll, coastline, prj)
    for isl in islands:
        ax.plot(isl[0], isl[1], transform=proj.pc, color='b')

    return mainland, islands


def get_land_hi():

    rinf = RegionInfo('hi', use_old_con_defs=False)
    prj = proj.proj[rinf.region]

    fig, ax = setup_figure(rinf, 32,
                           color='m',
                           eez_kw=dict(lw=3), alpha=0.5)

    # First grab the EEZ boundary (including borders)
    lleez = rinf.get_contour('eez')[0]

    # Now grab the coastline data
    coastline = get_coastline_data(prj.lonlim, prj.latlim, )

    # There is no mainland
    mainland = None

    # lleez is the boundary
    ax.plot(lleez[0], lleez[1], 'r-', transform=proj.pc)

    # Now grab and plot the islands
    islands = get_islands(lleez, coastline, prj)
    for isl in islands:
        ax.plot(isl[0], isl[1], transform=proj.pc, color='b')

    return mainland, islands


def get_land_ak():

    rinf = RegionInfo('ak', use_old_con_defs=False)
    prj = proj.proj[rinf.region]

    fig, ax = setup_figure(rinf, 33,
                           color='m',
                           eez_kw=dict(lw=3), alpha=0.5)

    # First grab the EEZ boundary (including borders)
    lleez = rinf.get_contour('EEZ')[0]
    # Now grab the coastline data
    lonlim = np.array(prj.lonlim)
    lonlim[1] += 5
    coastline = get_coastline_data(lonlim, prj.latlim, )

    # Now pick-out the segment of mainland that matches it
    mainland = get_mainland(lleez, coastline[0], prj)
    # Crop out Cook Inlet
    i0 = argnearest(np.array([-153.32,58.88]),mainland)
    i1 = argnearest(np.array([-151.75,59.18]),mainland)
    mainland = np.hstack((mainland[:, :i0], mainland[:, i1:]))

    # Construct the boundary
    boundary_ll = np.hstack((lleez, mainland))
    # ... and plot it
    ax.plot(boundary_ll[0], boundary_ll[1], 'r-', transform=proj.pc)

    # Now grab and plot the islands
    islands = get_islands(boundary_ll, coastline, prj)
    for isl in islands:
        ax.plot(isl[0], isl[1], transform=proj.pc, color='b')

    return mainland, islands


def get_land_at():

    rinf = RegionInfo('at', use_old_con_defs=False)
    prj = proj.proj[rinf.region]

    fig, ax = setup_figure(rinf, 34,
                           color='m',
                           eez_kw=dict(lw=3), alpha=0.5)

    # First grab the EEZ boundary (including borders)
    lleez = rinf.get_contour('EEZ')[0]
    # Now grab the coastline data
    coastline = get_coastline_data(prj.lonlim, prj.latlim, )

    # Now pick-out the segment of mainland that matches it
    mainland = get_mainland(lleez, coastline[0], prj)
    # crop out some large bays
    for ll0, ll1 in [
            [(-89.16, 29.38), (-88.70, 30.33)], # Miss. River Delta (LA)
            [(-75.95, 37.14), (-76.00, 36.92)], # Chesapeake
            [(-75.53, 35.82), (-76.35, 34.94)], # NC coast
    ]:
        i0 = argnearest(np.array(ll0),mainland)
        i1 = argnearest(np.array(ll1),mainland)
        if i0 > i1:
            tmp = i0;i0 = i1;i1 = tmp
        mainland = np.hstack((mainland[:, :i0], mainland[:, i1:]))

    # Construct the boundary
    boundary_ll = np.hstack((lleez, mainland))
    # ... and plot it
    ax.plot(boundary_ll[0], boundary_ll[1], 'r-', transform=proj.pc)

    # Now grab and plot the islands
    islands = get_islands(boundary_ll, coastline, prj)
    for isl in islands:
        ax.plot(isl[0], isl[1], transform=proj.pc, color='b')

    return mainland, islands


def get_land_gm():

    rinf = RegionInfo('gm', use_old_con_defs=False)
    prj = proj.proj[rinf.region]

    fig, ax = setup_figure(rinf, 35,
                           color='m',
                           eez_kw=dict(lw=3), alpha=0.5)

    # First grab the EEZ boundary (including borders)
    lleez = rinf.get_contour('EEZ')[0]
    # Now grab the coastline data
    coastline = get_coastline_data(prj.lonlim, prj.latlim, )

    # Now pick-out the segment of mainland that matches it
    mainland = get_mainland(lleez, coastline[0], prj)
    # crop out some large bays
    for ll0, ll1 in [
            [(-89.16, 29.38), (-88.70, 30.33)], # Miss. River Delta (LA)
    ]:
        i0 = argnearest(np.array(ll0),mainland)
        i1 = argnearest(np.array(ll1),mainland)
        if i0 > i1:
            tmp = i0;i0 = i1;i1 = tmp
        mainland = np.hstack((mainland[:, :i0], mainland[:, i1:]))

    # Construct the boundary
    boundary_ll = np.hstack((lleez, mainland))
    # ... and plot it
    ax.plot(boundary_ll[0], boundary_ll[1], 'r-', transform=proj.pc)

    # Now grab and plot the islands
    islands = get_islands(boundary_ll, coastline, prj)
    for isl in islands:
        ax.plot(isl[0], isl[1], transform=proj.pc, color='b')

    return mainland, islands


def get_land_ec():

    rinf = RegionInfo('ec', use_old_con_defs=False)
    prj = proj.proj[rinf.region]

    fig, ax = setup_figure(rinf, 36,
                           color='m',
                           eez_kw=dict(lw=3), alpha=0.5)

    # First grab the EEZ boundary (including borders)
    lleez = rinf.get_contour('EEZ')[0]
    # Now grab the coastline data
    coastline = get_coastline_data(prj.lonlim, prj.latlim, )

    # Now pick-out the segment of mainland that matches it
    mainland = get_mainland(lleez, coastline[0], prj)
    # crop out some large bays
    for ll0, ll1 in [
            [(-75.95, 37.14), (-76.00, 36.92)], # Chesapeake
            [(-75.53, 35.82), (-76.35, 34.94)], # NC coast
    ]:
        i0 = argnearest(np.array(ll0),mainland)
        i1 = argnearest(np.array(ll1),mainland)
        if i0 > i1:
            tmp = i0;i0 = i1;i1 = tmp
        mainland = np.hstack((mainland[:, :i0], mainland[:, i1:]))

    # Construct the boundary
    boundary_ll = np.hstack((lleez, mainland))
    # ... and plot it
    ax.plot(boundary_ll[0], boundary_ll[1], 'r-', transform=proj.pc)

    # Now grab and plot the islands
    islands = get_islands(boundary_ll, coastline, prj)
    for isl in islands:
        ax.plot(isl[0], isl[1], transform=proj.pc, color='b')

    return mainland, islands


def get_land_prusvi():

    rinf = RegionInfo('prusvi', use_old_con_defs=False)
    prj = proj.proj[rinf.region]

    fig, ax = setup_figure(rinf, 37,
                           color='m',
                           eez_kw=dict(lw=3), alpha=0.5)

    # First grab the EEZ boundary (including borders)
    lleez = rinf.get_contour('EEZ')[0]

    # Now grab the coastline data
    coastline = get_coastline_data(prj.lonlim, prj.latlim, )

    # There is no mainland
    mainland = None

    # lleez is the boundary
    ax.plot(lleez[0], lleez[1], 'r-', transform=proj.pc)

    # Now grab and plot the islands
    islands = get_islands(lleez, coastline, prj)
    for isl in islands:
        ax.plot(isl[0], isl[1], transform=proj.pc, color='b')

    return mainland, islands


def run_all():

    land_data = {}

    land_data['wc'] = get_land_wc()
    land_data['hi'] = get_land_hi()
    land_data['ak'] = get_land_ak()
    land_data['at'] = get_land_at()
    land_data['gm'] = get_land_gm()
    land_data['ec'] = get_land_ec()
    land_data['prusvi'] = get_land_prusvi()

    return land_data
