import wave_res as wr
import numpy as np

all_months = np.arange(np.datetime64('2009-09'),
                       np.datetime64('2010-01'))

all_regions = ['ak', 'wc', 'at', 'prusvi']
ranges = np.arange(10, 201, 10)

# In this dictionary, the keys are the scenarios, and the values are lists of regions that should be run for each scenario.
run_these = {
    #'baseline': ['prusvi'],
    'extraction': ['prusvi'],
}


for scenario, regions in run_these.items():

    for regn in regions:

        rinf = wr.RegionInfo(regn)

        # Generate temporary data files, and/or load their data
        dat = wr.remote.process_and_load(scenario, regn, all_months)
