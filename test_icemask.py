from wave_res.calc_local import calc_local
import numpy as np

months = np.arange(np.datetime64('1979-01'),
                   np.datetime64('1979-07'))

local = calc_local('baseline', 'ak', months)

local.to_hdf5('ak-icemask-test.h5')
