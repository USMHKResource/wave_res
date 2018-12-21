import wave_res.calc_local as cl
reload(cl)
import wave_res as wr
import numpy as np
from wave_res.paths import mkdir

months = np.arange(np.datetime64('2009-09'),
                   np.datetime64('2010-01'))

all_regions = ['ak', 'wc', 'at', 'prusvi']
ranges = np.arange(10, 201, 10)

# In this dictionary, the keys are the scenarios, and the values are lists of regions that should be run for each scenario.
run_these = {
    #'baseline': ['prusvi'],
    'extraction': ['prusvi'],
}


for scenario, REGIONS in run_these.items():

    mkdir('results/{scenario}'.format(scenario=scenario))
    for region in REGIONS:

        rinf = wr.RegionInfo(region)

        print("#### Calculating totals for '{}' scenario, '{}' region..."
              .format(scenario, region))

        print("   Calculating remote resource..."
              .format(scenario, region))
        # Calculate the remote resource
        res = wr.remote.calc_remote(scenario, region, months)
        # This returns a dictionary-like object (based on pyDictH5.data)
        # containing:
        #  'time': (n_months) the month
        #  'Nhour': (n_months) The number of hours in each month
        #  'ranges': (n_ranges) The range of each contour [nautical miles]
        #  'length': (n_ranges) The length of each contour [meters]
        #  '1way': (n_months, n_ranges) Monthly averaged wave energy
        #          flux using the '1way' method [watts]
        #  'trad': (n_months, n_ranges) Monthly averaged wave energy
        #          flux using a traditional dot-product [watts]
        #  'bdir': (n_months, n_ranges) Monthly averaged wave energy
        #          flux using a bi-directional dot-product [watts]
        #  'unit': (n_months, n_ranges) Monthly averaged wave energy
        #          flux using the unit-circle method [watts]

        # Compute the hour-weighted (rather than month-weighted) averages
        tot = res.hourly_average()
        print(
            "    The remote resource at the eez boundary for the '{region}' region, \n"
            "    scenario '{scenario}', from {start} to {end} is:\n"
            "      1-way:  {oway: 8.3f} GW\n"
            "      trad:   {trad: 8.3f} GW\n"
            "      bi-dir: {bdir: 8.3f} GW\n"
            "      unit:   {unit: 8.3f} GW\n"
            .format(
                region=region, scenario=scenario,
                start=months[0], end=months[-1],
                # The last range is the eez (thus -1)
                oway=tot['1way'][-1] / 1e9, trad=tot['trad'][-1] / 1e9,
                bdir=tot['bdir'][-1] / 1e9, unit=tot['unit'][-1] / 1e9,
            )
        )

        res.to_hdf5('results/{scenario}/{region}.remote-totals.h5'
                    .format(scenario=scenario, region=region))
        # You can load this data with res = pyDictH5.load(<fname>)
