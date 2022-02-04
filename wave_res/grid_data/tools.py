import numpy as np

def group(bl, min_length=0):
    """Find continuous segments in a boolean array.

    Parameters
    ----------
    bl : |np.ndarray| (dtype='bool')
      The input boolean array.
    min_length : int (optional)
      Specifies the minimum number of continuos points to consider a
      `group` (i.e. that will be returned).

    Returns
    -------
    out : np.ndarray(slices,)
      a vector of slice objects, which indicate the continuous
      sections where `bl` is True.

    Notes
    -----
    This function has funny behavior for single points.  It will
    return the same two indices for the beginning and end.

    """
    if not any(bl):
        return np.empty(0)
    vl = np.diff(bl.astype('int'))
    ups = np.nonzero(vl == 1)[0] + 1
    dns = np.nonzero(vl == -1)[0] + 1
    if bl[0]:
        if len(ups) == 0:
            ups = np.array([0])
        else:
            ups = np.concatenate((np.arange([0]), [len(ups)]))
    if bl[-1]:
        if len(dns) == 0:
            dns = np.array([len(bl)])
        else:
            dns = np.concatenate((dns, [len(bl)]))
    out = np.empty(len(dns), dtype='O')
    idx = 0
    for u, d in zip(ups, dns):
        if d - u < min_length:
            continue
        out[idx] = slice(u, d)
        idx += 1
    return out[:idx]
