import wave_res as wr
import numpy as np
from wave_res.paths import mkdir

# Range to compute the resource
months = np.arange(np.datetime64('2010-01'),
                   np.datetime64('2011-01'))

all_regions = ['wc', 'at', 'prusvi', 'ak', 'hi','gm','ec']
# In this dictionary, the keys are the scenarios, and the values are
# lists of regions that should be run for each scenario.
run_these = {
#    'baseline': all_regions,
#    'extraction': all_regions,
    'baseline':['prusvi'],
#    'extraction':['gm','ec']
}


for scenario, REGIONS in run_these.items():

    mkdir('results/{scenario}'.format(scenario=scenario))
    for region in REGIONS:

        rinf = wr.RegionInfo(region)

        print("#### Calculating totals for '{}' scenario, '{}' region..."
              .format(scenario, region))

        print("   Calculating remote resource...")
        # Calculate the remote resource
        remote = wr.calc_remote(scenario, region, months)
        # This returns a dictionary-like object (based on pyDictH5.data)
        # containing:
        #  'time': (n_months) the month
        #  'Nhour': (n_months) The number of hours in each month
        #  'range': (n_ranges) The range of each contour [nautical miles]
        #  'length': (n_ranges) The length of each contour [meters]
        #  'fbins': (n_freq + 1) The edges of the frequency bins [Hz]
        #  '1way': (n_months, n_freq, n_ranges) Monthly averaged wave energy
        #          flux using the '1way' method [watts/Hz]
        #  'trad': (n_months, n_freq, n_ranges) Monthly averaged wave energy
        #          flux using a traditional dot-product [watts/Hz]
        #  'bdir': (n_months, n_freq, n_ranges) Monthly averaged wave energy
        #          flux using a bi-directional dot-product [watts/Hz]
        #  'unit': (n_months, n_freq, n_ranges) Monthly averaged wave energy
        #          flux using the unit-circle method [watts/Hz]

        # Compute the hour-weighted (rather than month-weighted) averages
        rtot = remote.calc_totals()

        # Print the remote results
        print(
            "    The average remote resource at the eez boundary for the\n"
            "    '{region}' region, scenario '{scenario}',\n"
            "    from {start} to {end} is:\n"
            "      1-way:  {oway: 8.3f} GW\n"
            "      trad:   {trad: 8.3f} GW\n"
            "      bi-dir: {bdir: 8.3f} GW\n"
            "      unit:   {unit: 8.3f} GW\n"
            .format(
                region=region, scenario=scenario,
                start=months[0], end=months[-1],
                # The last range is the eez (thus -1)
                oway=rtot['1way'][-1] / 1e9, trad=rtot['trad'][-1] / 1e9,
                bdir=rtot['bdir'][-1] / 1e9, unit=rtot['unit'][-1] / 1e9,
            )
        )

        remote.to_hdf5('results/{scenario}/{region}.remote-totals.h5'
                       .format(scenario=scenario, region=region))
        # You can load this data with pyDictH5.load(<fname>)

        print("   Calculating local resource...")

        local = wr.calc_local(scenario, region, months,fc=True)
        # This returns a dictionary-like object (based on pyDictH5.data)
        # containing:
        #  'time': (n_months) the month
        #  'Nhour': (n_months) The number of hours in each month
        #  'range': (n_ranges) The range of each contour [nautical miles]
        #  'area': (n_ranges) The area between each contour [meters^2]
        #  'fbins': (n_freq + 1) The edges of the frequency bins [Hz]
        #  'sin': (n_months, n_freq, n_ranges) Monthly averaged 'wind input'
        #         source term. [watts/Hz]
        #  'sds': (n_months, n_freq, n_ranges) Monthly averaged 'dissipation'
        #         source term. [watts/Hz]
        #  'snl': (n_months, n_freq, n_ranges) Monthly averaged 'non-linear'
        #         source term. [watts/Hz]
        #  'sbt': (n_months, n_freq, n_ranges) Monthly averaged 'bottom
        #         friction' source term. [watts/Hz]
        #  'sbt': (n_months, n_freq, n_ranges) Monthly averaged 'ice'
        #         source term. [watts/Hz]
        #  'stot': (n_months, n_freq, n_ranges) Monthly averaged 'total'
        #         source term. [watts/Hz]
        # NOTE: the source terms are totals for the area between that
        # contour range (i.e., in local['range']), and the range
        # inshore of it. So, to get the source term for the entire
        # EEZ, you need to sum them, e.g.:
        #    sin_total = local['sin'].sum(-1)

        ltot = local.calc_totals()

        print(
            "    The average local resource over the EEZ for the\n"
            "    '{region}' region, scenario '{scenario}',\n"
            "    from {start} to {end} is:\n"
            "      sin:  {sin: 8.3f} GW\n"
            "      sds:  {sds: 8.3f} GW\n"
            "      snl:  {snl: 8.3f} GW\n"
            "      sbt:  {sbt: 8.3f} GW\n"
            "      sice: {sice: 8.3f} GW\n"
            "      stot: {stot: 8.3f} GW\n"
            .format(
                region=region, scenario=scenario,
                start=months[0], end=months[-1],
                # Sum over ranges to get total
                sin=ltot['sin'].sum(-1) / 1e9, sds=ltot['sds'].sum(-1) / 1e9,
                snl=ltot['snl'].sum(-1) / 1e9, sbt=ltot['sbt'].sum(-1) / 1e9,
                sice=ltot['sice'].sum(-1) / 1e9, stot=ltot['stot'].sum(-1) / 1e9,
            )
        )

        local.to_hdf5('results/{scenario}/{region}.local-totals.h5'
                      .format(scenario=scenario, region=region))
