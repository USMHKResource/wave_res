import calc_remote as cr
import numpy as np
import base


if __name__ == '__main__':

    months = np.arange(np.datetime64('2009-01'), np.datetime64('2010-01'))
    region = 'wc'
    contour = 'EEZ'
    con_inds = base.con_defs[region][contour]

    # months = np.arange(np.datetime64('2009-01'), np.datetime64('2010-01'))
    # region = 'at'
    # con_inds = None

    dat = cr.process_and_load(region, months, overwrite=False)

    if con_inds is not None:
        total = cr.calc_total(dat, con_inds)
        cr.print_total(total, region, contour)
