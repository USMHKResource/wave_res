import pyDictH5 as pdh5
import numpy as np


source_terms = ['sbt', 'sds', 'snl', 'stot', 'sin', 'sice']
remote_terms = ['trad', '1way', 'bdir', 'unit']
iranges = range(20)

# Choose units here.
unit = 'GW'
unit = 'TWh/yr'
_factordict = {
    'GW': 1e-9,
    'TWh/yr': 365 * 24 * 1e-12,
}
_factor = _factordict[unit]

regions = {'ak', 'wc', 'ec', 'hi', 'gm', 'prusvi', }

colors = {'ak': 'k',
          'wc': 'b',
          'ec': 'r',
          'hi': 'g',
          'gm': 'm',
          'prusvi': 'lightcoral',
}

names = {'ak': 'Alaska',
         'ec': 'East Coast',
         'wc': 'West Coast',
         'prusvi': 'Puerto Rico + U.S.V.I.',
         'gm': 'Gulf of Mexico',
         'hi': 'Hawaii',
}


class _Result(object):
    """Abstract base class for results objects"""

    @property
    def color(self, ):
        return colors[self.region]

    @property
    def name(self, ):
        return names[self.region]


class TotResult(_Result):
    """Container for results objects

    Mostly one uses the ``.remote`` and ``.local`` attributes here.
    """
    
    def __init__(self, region, data_source='results/freq.fcut/'):
        self.region = region
        self.remote = RemoteResult(region, data_source=data_source)
        self.local = LocalResult(region, data_source=data_source)
        self.f = self.remote.f
        self.time = self.remote._data['baseline']['time']
        self.unit = unit

    
        
class Result(_Result):

    @property
    def _Nhour(self, ):
        return self._data['baseline']['Nhour']
    
    def int_freq(self, unit='TWh/yr', source='baseline'):
        """Integrate frequency information."""
        out = {}
        dat = self._data[source]
        for ky in self._terms:
            out[ky] = (dat[ky] * self.df[None, :, None]).sum(1) * _factor
        return out

    def avg_time(self, source='baseline'):
        """A time-average only (frequency info retained).
        """
        out = {}
        dat = self._data[source]
        for ky in self._terms:
            out[ky] = np.average(dat[ky],
                                 weights=self._Nhour,
                                 axis=0) * _factor


    def avg_yearly(self, source='baseline'):
        dat = self.int_freq()
        Nhour = np.tile(self._Nhour.reshape((-1, 12))[:, :, None], (1, 1, 20))
        for ky in self._terms:
            d = dat[ky]
            shp = list(d.shape)
            shp = [-1, 12] + shp[1:]
            dat[ky] = np.average(dat[ky].reshape(shp),
                                 weights=Nhour,
                                 axis=1)
        return dat

    def avg_annual(self, source='baseline'):
        dat = self.int_freq()
        for ky in self._terms:
            d = dat[ky]
            shp = list(d.shape)
            shp = [-1, 12] + shp[1:]
            dat[ky] = np.average(dat[ky].reshape(shp), axis=0)
        return dat

    def total(self, source='baseline'):
        dat = self.int_freq(source=source)
        for ky in self._terms:
            dat[ky] = np.average(dat[ky],
                                 weights=self._data[source]['Nhour'],
                                 axis=0) # no `* _factor` b/c that happens in int_freq
        return dat

class RemoteResult(Result):

    def __init__(self, region, data_source='results/freq.fcut/'):
        self.region = region
        self._data_source = data_source
        self._terms = remote_terms
        self._data = {'baseline': pdh5.load(data_source + '{}/{}.remote-totals.h5'
                                .format('baseline', region)),
                      'extraction': pdh5.load(data_source + '{}/{}.remote-totals.h5'
                                .format('extraction', region)),
                                }
        self.df = np.diff(self._data['baseline']['fbins'])
        self.f = self._data['baseline']['fbins'][:-1] + self.df / 2


class LocalResult(Result):

    def __init__(self, region, data_source='results/freq.fcut/'):
        self.region = region
        self._data_source = data_source
        self._terms = source_terms
        self._data = {'baseline': pdh5.load(data_source + '{}/{}.local-totals.h5'
                                .format('baseline', region)),
                      'extraction': pdh5.load(data_source + '{}/{}.local-totals.h5'
                                .format('extraction', region)),
                                }
        self.df = np.diff(self._data['baseline']['fbins'])
        self.f = self._data['baseline']['fbins'][:-1] + self.df / 2


class WrapMonths(object):

    def __init__(self, shift=0, offset=0.5):
        self._n = n = 12
        self.shift = shift
        self._index = list(np.arange(shift, shift + n + 1) % n)
        self._index_seasons = (np.arange(shift, shift + n + 1, 3) % n) / 3
        self.labels = np.array(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'June', 'July', 'Aug', 'Sept', 'Oct', 'Nov', 'Dec'])[self._index]
        self.x = np.arange(n + 1.) + offset + shift
        self.xticks = np.arange(n + 1.) + shift

    @property
    def xlim(self, ):
        return [self.x[0], self.x[-1]]
        
    def __call__(self, y):
        if not isinstance(y, np.ndarray):
            y = np.array(y)
        return y[self._index]

    @property
    def season_edges(self, ):
        return np.arange(3, 16, 3) + 21. / 30.0 + self.shift - 1

    @property
    def season_labels(self):
        return np.array(['spring', 'summer', 'fall', 'winter', ])[self._index_seasons]

    @property
    def season_centers(self):
        out = self.season_edges + 1.5
        out[out > self.xlim[-1]] -= self._n
        return out
        
