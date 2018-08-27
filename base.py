import cPickle as _pkl
import paths as p

regions = ['wc', 'at', 'hi', 'ak', 'prusvi']
conids = ['EEZ'] + ['{:03d}'.format(n) for n in range(10, 200, 10)]

with open(str(p.projdir / 'data/Contour_Ranges.pkl'), 'r') as fl:
    con_defs = _pkl.load(fl)
with open(str(p.projdir / 'data/FreqBins.pkl'), 'r') as fl:
    freqbins = _pkl.load(fl)

# This is the outer-boundary of the EEZ (not including the Canada
# + Mexico borders)
con_defs['wc']['EEZ'] = slice(34, 148)
