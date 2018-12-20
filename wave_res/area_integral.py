import numpy as np


def sumtriangles(xy, z, triangles):
    """ integrate scattered data, given a triangulation
    zsum, areasum = sumtriangles( xy, z, triangles )
    In:
        xy: npt, dim data points in 2d, 3d ...
        z: npt data values at the points, scalars or vectors
        triangles: ntri, dim+1 indices of triangles or simplexes, as from
http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.Delaunay.html
    Out:
        zsum: sum over all triangles of (area * z at midpoint).
            Thus z at a point where 5 triangles meet
            enters the sum 5 times, each weighted by that triangle's area / 3.
        areasum: the area or volume of the convex hull of the data points.
            For points over the unit square, zsum outside the hull is 0,
            so zsum / areasum would compensate for that.
            Or, make sure that the corners of the square or cube are in xy.
    """
    # I found this function here (what a lifesaver!)
    # https://stackoverflow.com/questions/5941113/looking-for-python-package-for-numerical-integration-over-a-tessellated-domain

    # z concave or convex => under or overestimates
    npt, dim = xy.shape
    ntri, dim1 = triangles.shape
    assert npt == len(z), "shape mismatch: xy %s z %s" % (xy.shape, z.shape)
    assert dim1 == dim + 1, "triangles ? %s" % triangles.shape
    zsum = np.zeros(z[0].shape)
    areasum = 0
    dimfac = np.prod(np.arange(1, dim + 1))
    for tri in triangles:
        corners = xy[tri]
        t = corners[1:] - corners[0]
        if dim == 2:
            area = abs(t[0, 0] * t[1, 1] - t[0, 1] * t[1, 0]) / 2
        else:
            area = abs(np.linalg.det(t)) / dimfac  # v slow
        zsum += area * z[tri].mean(axis=0)
        areasum += area
    return (zsum, areasum)


if __name__ == "__main__":
    from time import time
    from scipy.spatial import Delaunay

    npt = 500
    dim = 2
    seed = 1

    np.set_printoptions(2, threshold=100, edgeitems=5, suppress=True)
    np.random.seed(seed)

    points = np.random.uniform(size=(npt, dim))
    z = points  # vec; zsum should be ~ constant
    # z = points[:,0]
    t0 = time()
    tessellation = Delaunay(points)
    t1 = time()
    triangles = tessellation.vertices  # ntri, dim+1
    zsum, areasum = sumtriangles(points, z, triangles)
    t2 = time()

    print "%s: %.0f msec Delaunay, %.0f msec sum %d triangles:  zsum %s  areasum %.3g" % (
        points.shape, (t1 - t0) * 1000, (t2 - t1) * 1000,
        len(triangles), zsum, areasum )

    xy = np.array([[0., 0], [1, 0], [1, 1],
                   [3, 1], [3, 0], [0, 1]])
    z = np.array([1., 1, 1, 3, 2, 1])
    tri = Delaunay(xy)
    zsum, areasum = sumtriangles(xy, z, tri.vertices)
    print(zsum)
    
    import matplotlib.pyplot as plt
    plt.clf()
    plt.plot(xy[:,0],xy[:,1],'.')
