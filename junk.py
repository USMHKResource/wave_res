# This code was to try to find the intersection point between the
# distance-from-shore contours, and the 'borders'. However, I don't
# think it is necessary because I don't want to add points to the
# dataset (I don't have source data there.)

def argnearest(pt, line):
    d2 = ((pt[:, None] - line) ** 2).sum(0)
    return np.argmin(d2)


def line_intersection(seg, mline):
    """Find the point on `mline` where an infinite extension of `seg` intersects it.
    If there are multiple intersections, choose the one closest to the first point in seg.

    seg : (2, 2)

    mline : (2, N)
    """

    def line(p):
        B, A= np.diff(p, 1, 1)
        C = p[0, :-1] * p[1, 1:] - p[0, 1:] * p[1, :-1]
        return -A, B, -C

    def intersection(L1, L2):
        D  = L1[0] * L2[1] - L1[1] * L2[0]
        Dx = L1[2] * L2[1] - L1[1] * L2[2]
        Dy = L1[0] * L2[2] - L1[2] * L2[0]
        out = np.empty((2, len(L2[2])))
        out[0] = Dx / D
        out[1] = Dy / D
        return out

    Lseg = line(seg)
    Lm = line(mline)

    pts = intersection(Lseg, Lm)
    gds = (((mline[0, :-1] < pts[0]) & (pts[0] < mline[0, 1:])) |
           ((mline[0, 1:] < pts[0]) & (pts[0] < mline[0, :-1])))
    out = pts[:, gds]
    if out.shape[1] > 1:
        ind = argnearest(seg[:, 0], out)
        out = out[:, ind]
    else:
        out = out[:, 0]
    return out

ax.plot(bounds[0], bounds[1], 'r-', transform=proj.pc)

ax.plot(pgy.boundary.xy[0], pgy.boundary.xy[1], 'b-')

brd = rinf.get_contour('borders')[0]
brd = rinf.proj.transform_points(proj.pc, brd[0], brd[1]).T[:2]
con = rinf.get_contour('190')[0]
con = rinf.proj.transform_points(proj.pc, con[0], con[1]).T[:2]

tmp = line_intersection(con[:2, :2], brd)

ipt = argnearest(con[:, 0], brd)
# ax.plot(lleez[0, 0], lleez[1, 0], 'r+', transform=proj.pc)
# ax.plot(lleez[0, -1], lleez[1, -1], 'r+', transform=proj.pc)
