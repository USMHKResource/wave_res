import paths as p
from netCDF4 import Dataset
import matplotlib.pyplot as plt

d0 = Dataset(p.srcdir + '2009/eez/ww3.ak.200901_spec.nc', 'r')

fignum = 10
fig = plt.figure(10)
fig.clf()
fig, ax = plt.figure(1, 1)


# for idx in 
