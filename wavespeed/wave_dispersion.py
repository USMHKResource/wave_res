import numpy as np
import scipy.interpolate as interp
from scipy.optimize import brentq
import pickle as pkl
import gzip
import os

"""The surface gravity wave dispersion relation.

All frequencies, f, are in Hz, and wavenumbers, k are in rad/m. Depth is in meters.

"""

gravity = 9.81

thisdir = os.path.dirname(os.path.abspath(__file__)) + '/'

try:
    with gzip.open(thisdir + 'dispersion_fit.pkl.gz', 'rb',
                   compresslevel=3) as gz:
        k_func = pkl.load(gz)
except IOError:
    if __name__ != '__main__':
        raise Exception("Run this module as a script to generate the "
                        "inverse-dispersion lookup-table.")


def disp(k, d):
    """Surface gravity wave dispersion relation.

    Parameters
    ----------
    k : wavenumber, rad/m

    d : depth, m

    Returns
    -------
    f : frequency, Hz
    """
    return 1. / (2 * np.pi) * np.sqrt(gravity * k * np.tanh(k * d))


def inv_disp(f, d):
    """The *inverse* surface gravity wave dispersion relation.

    Parameters
    ----------
    f : frequency, Hz

    d : depth, m

    Returns
    -------
    k : wavenumber, rad/m

    Notes
    -----
    This function works by using a lookup table, and using a linear
    fit between points. The lookup table is created when this module
    is run as a script.
    """
    return k_func((d, f))


def c_group0(f, d):
    k = inv_disp(f, d)
    return np.pi * f / k * (1 + (2 * k * d) / np.sinh(2 * k * d))


def c_group(k, d):
    f = disp(k, d)
    return np.pi * f / k * (1 + (2 * k * d) / np.sinh(2 * k * d))


def compare_disp():
    # This looks at dispersion relation errors
    dvals = np.arange(1, 50)[:, None] * 20
    l = np.arange(1, 51) * 5.
    kvals = 2 * np.pi / l
    fvals = disp(kvals, dvals)

    k2 = inv_disp(fvals, dvals)
    f2 = disp(k2, dvals)

    f_err = (fvals - f2) / fvals
    f_err_rms = np.sqrt((f_err ** 2).mean())

    print("")
    print("The RMS freq error is: {:06.4f}%".format(f_err_rms * 100))
    print("The max |freq| error is: {:06.4f}% (ind: {})"
          .format(np.abs(f_err).max() * 100, f_err.argmax()))

    cvals = c_group(kvals, dvals)

    c2 = c_group0(f2, dvals)
    cg_err = (cvals - c2) / cvals
    cg_err_rms = np.sqrt((cg_err ** 2).mean())

    print("")
    print("The RMS c_g error is: {:06.4f}%".format(cg_err_rms * 100))
    print("The max |c_g| error is: {:06.4f}% (ind: {})"
          .format(np.abs(cg_err).max() * 100, cg_err.argmax()))


if __name__ == '__main__':

    # ######
    # This was an old method that was too slow
    # d = np.logspace(-1, 4.5, 300)
    # k = np.logspace(-5.5, 2, 300)
    # k, d = np.meshgrid(k, d)
    # k_func = interp.LinearNDInterpolator((disp(k, d).flatten(),
    #                                       d.flatten()),
    #                                      k.flatten())
    # # k_func = interp.NearestNDInterpolator((disp(k, d).flatten(),
    # #                                        d.flatten()),
    # #                                       k.flatten())

    # def inv_disp(f, d):
    #     return k_func(f, d)

    def rt_fn(k, f, d):
        return disp(k, d) - f

    f = np.logspace(-2, 1.5, 500)
    d = np.logspace(-1, 4.5, 500)
    karr = np.empty((len(d), len(f)))
    for idf, freq in enumerate(f):
        for idd, depth in enumerate(d):
            karr[idd, idf] = brentq(rt_fn, 1e-9, 1e9, (freq, depth))

    k_func = interp.RegularGridInterpolator((d, f), karr, method='linear',
                                            # Allow extrapolation:
                                            bounds_error=False, fill_value=None)

    compare_disp()

    with gzip.open('dispersion_fit.pkl.gz', 'wb', compresslevel=3) as gz:
        pkl.dump(k_func, gz)
