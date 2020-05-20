import numpy as np
import os

try:
    _thisdir = os.path.dirname(os.path.realpath(__file__))
except NameError:
    _thisdir = './'

enso = {}
enso['time'] = np.arange(np.datetime64('1979-01'),
                      np.datetime64('2011-01'))

####
# Load the ONI data
_dtmp = np.loadtxt(_thisdir + '/oni.data', skiprows=1, max_rows=63)

_range = slice(29, 61)
# This grabs the 1979-2010 segment of data
enso['oni'] = _dtmp[_range, 1:].flatten()

# Confirm we're getting the right data
assert (_dtmp[_range, 0] == np.arange(1979, 2011)).all()

####
# Load the meiv2 data
_dtmp = np.loadtxt(_thisdir + '/meiv2.data', skiprows=1, max_rows=32)

enso['meiv2'] = _dtmp[:, 1:].flatten()

# Confirm we're getting the right data
assert (_dtmp[:, 0] == np.arange(1979, 2011)).all()

####
# Load the PDO data (not technically enso, but we're using it for similar purposes)
_dtmp = np.loadtxt(_thisdir + '/pdo.csv', skiprows=2, delimiter=',')

_range = slice(1500, 1884)

enso['pdo'] = _dtmp[_range, 1]

assert (_dtmp[_range.start, 0] == 197901 and _dtmp[_range.stop, 0] == 201101)
