from __future__ import print_function
from multiprocessing import Pool
import numpy as np
from netCDF4 import Dataset
from datetime import datetime as dt
import paths as p
import time

# try:
#     from tqdm import tqdm
# except ImportError:
def tqdm(input, **kwargs):
    return input

def timestamp(dt):
    return time.mktime(dt.timetuple())


def read_1d_source(fname):

    out = {}
    print('Reading file: {}'.format(fname))
    with open(str(fname), 'rb') as fd:
        n_blocks, n_freq, terms, units = _check_file(fd)
        n_terms = len(terms) - 1
        time = out['time'] = np.zeros(n_blocks, dtype=np.uint64)
        tag = out['tag'] = np.zeros(n_blocks, dtype='S16')
        lon = out['lon'] = np.zeros(n_blocks, dtype=np.float32)
        lat = out['lat'] = np.zeros(n_blocks, dtype=np.float32)
        depth = out['depth'] = np.zeros(n_blocks, dtype=np.float32)
        ustar = out['ustar'] = np.zeros(n_blocks, dtype=np.float32)
        U10 = out['U10'] = np.zeros(n_blocks, dtype=np.float32)
        data = np.zeros((n_blocks, n_freq, n_terms), dtype=np.float32)
        for idx in tqdm(range(n_blocks)):
            (time[idx],
             tag[idx], lon[idx], lat[idx],
             depth[idx], ustar[idx], U10[idx]) = _read_header(fd)
            tmp = _read_table(fd, n_freq)
            if idx == 0:
                f = tmp[:, 0]
            data[idx] = tmp[:, 1:]
    for ky in out:
        out[ky] = out[ky].reshape(24, -1)
    data = data.reshape(24, -1, n_freq, n_terms)
    out['_units'] = units
    out['_terms'] = terms
    out['f'] = f
    out['data'] = data
    for ky in ['lon', 'lat', 'depth', 'tag']:
        out[ky] = out[ky][0]
    out['time'] = out['time'][:, 0]
    return out


def _check_file(fd):
    lastline = None
    c = 0
    table_start = None
    while True:
        ln = fd.readline()
        c += 1
        if len(ln) > 2 and ln == (' ' + ('-' * (len(ln) - 2)) + '\n'):
            table_start = c
            units = lastline.split()
        if ln == ' \n' and lastline == ' \n':
            break
        if lastline == '\n' and ln.startswith('    f'):
            terms = ln.split()
        lastline = ln
    b = fd.tell()
    fd.seek(0, 2)  # go to the end
    end = fd.tell()
    n_blocks = end // b
    # print('Table start line: {: 10d}\n'
    #       'Number of terms:  {: 10d}\n'
    #       'Lines per block:  {: 10d}\n'
    #       'Bytes per block:  {: 10d}\n'
    #       'Total bytes:      {: 10d}\n'
    #       'Total blocks:     {: 10d}\n'
    #       .format(table_start, len(terms), c, b, end, n_blocks))
    n_freq = c - table_start - 2
    fd.seek(0, 0)
    return n_blocks, n_freq, terms, units


def _read_header(fd):
    line_num = 0
    # We use a loop here, because we can't mix loops
    # (used in genfromtxt) with 'readline'
    for ln in fd:
        if line_num == 0:
            t = int(timestamp(dt.strptime(ln.split(' : ')[-1],
                                          '%Y/%m/%d %H:%M:%S %Z\n')))
        elif line_num == 1:
            tag, ll = ln.split(' : ')[-1].split('  ', 1)
            ll = ll.split()[1:3]
            lon = float(ll[0])
            lat = float(ll[1])
        elif line_num == 2:
            depth = float(ln.split(' : ')[-1].split()[0])
        elif line_num == 3:
            ustar = float(ln.split(' : ')[-1].split()[0])
        elif line_num == 4:
            U10 = float(ln.split(' : ')[-1].split()[0])
            break
        line_num += 1
    return t, tag, lon, lat, depth, ustar, U10


def _read_table(fd, n_freq):
    dat = np.genfromtxt(fd, dtype=np.float32, skip_header=4, max_rows=n_freq)
    next(fd);next(fd); # Skip the next two lines
    return dat


def write2ncdf(fname, data):
    print("Writing to file: {}...".format(fname), end='')
    fd = Dataset(fname, 'w')
    dims = data['data'].shape
    fd.createDimension('time', dims[0])
    fd.createDimension('spatial', dims[1])
    fd.createDimension('freq', dims[2])
    for name in ['lat', 'lon', 'tag', 'depth']:
        var = fd.createVariable(name, data[name].dtype, ['spatial'], zlib=True)
        var[:] = data[name]
    fd.variables['depth'].units = 'm'
    var = fd.createVariable('time', data['time'].dtype, ['time'])
    var[:] = data['time']
    var.units = 'unixtime'
    var = fd.createVariable('freq', data['f'].dtype, ['freq'])
    var[:] = data['f']
    var.units = 'Hz'
    for name in ['ustar', 'U10', ]:
        var = fd.createVariable(name, data[name].dtype, ['time', 'spatial'])
        var[:] = data[name]
    fd.variables['ustar'].units = 'm/s'
    fd.variables['U10'].units = 'm/s'
    for idx, (term, unit) in enumerate(zip(data['_terms'][1:],
                                           data['_units'][1:])):
        var = fd.createVariable(term, data['data'].dtype,
                                ['time', 'spatial', 'freq'],
                                zlib=True)
        var.units = unit[1:-1]  # skip the parenthesis
        var[:] = data['data'][..., idx]
    fd.close()
    print("  Done.")


def read_and_write(fname):
    bname = fname.name.rstrip('.src')
    out = read_1d_source(fname)
    write2ncdf(str(p.tmpdir.joinpath(bname + '.nc')), out)
    return fname


if __name__ == '__main__':
    """This script reads '1D WW3 source files' (directionality integrated
    away), and writes them to gzip-compressed net-cdf files.
    """

    fnames = list(p.srcdir.glob('2009/src/ww3.*.1d.src'))
    #print(fnames)

    pl = Pool(4)
    res = pl.map(read_and_write, fnames[142:])
    print(res)

    # for fname in fnames[140:]:
    #     read_and_write(fname)

    # #fname = p.srcdir + '/2009/src/ww3.ak.20090519.1d.src'
    # fname = p.srcdir + '/2009/src/ww3.ak.20090520.1d.src'
    # bname = fname.rsplit('/', 1)[-1].rstrip('.src')
    # out = read_1d_source(fname)
    # write2ncdf(p.tmpdir + bname + '.nc', out)
