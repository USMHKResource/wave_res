import wave_res as wr
import numpy as np
from wave_res.paths import mkdir

# Range to compute the resource
months = np.arange(np.datetime64('1979-01'),
                   np.datetime64('2011-01'))

all_regions = ['wc', 'at', 'prusvi', 'ak', 'hi']

# If True: Stores results that are not area integrated
storeSpatial = True

# In this dictionary, the keys are the scenarios, and the values are
# lists of regions that should be run for each scenario.
run_these = {
    'baseline': ['wc','at'],
    'extraction': ['wc','at']
    # 'extraction':all_regions,
    # 'baseline':all_regions
}


for scenario, REGIONS in run_these.items():

    mkdir('results/{scenario}'.format(scenario=scenario))
    for region in REGIONS:
        
        rinf = wr.RegionInfo(region)

        print("#### Calculating LOCAL totals for '{}' scenario, '{}' region..."
              .format(scenario, region))

        print("   Calculating local resource...")

        local = wr.calc_local(scenario, region, months, 
                              storeSpatial=storeSpatial)
        # This returns a dictionary-like object (based on pyDictH5.data)
        # containing:
        #  'time': (n_months) the month
        #  'Nhour': (n_months) The number of hours in each month
        #  'range': (n_ranges) The range of each contour [nautical miles]
        #  'area': (n_ranges) The area between each contour [meters^2]
        #  'sin': (n_months, n_ranges) Monthly averaged 'wind input'
        #         source term. [watts]
        #  'sds': (n_months, n_ranges) Monthly averaged 'dissipation'
        #         source term. [watts]
        #  'snl': (n_months, n_ranges) Monthly averaged 'non-linear'
        #         source term. [watts]
        #  'sbt': (n_months, n_ranges) Monthly averaged 'bottom
        #         friction' source term. [watts]
        #  'sbt': (n_months, n_ranges) Monthly averaged 'ice'
        #         source term. [watts]
        #  'stot': (n_months, n_ranges) Monthly averaged 'total'
        #         source term. [watts]
        # NOTE: the source terms are totals for the area between that
        # contour range (i.e., in local['range']), and the range
        # inshore of it. So, to get the source term for the entire
        # EEZ, you need to sum them, e.g.:
        #    sin_total = local['sin'].sum(-1)

        if storeSpatial:
            local.to_hdf5('results/{scenario}/{region}.local-area.h5'
                          .format(scenario=scenario, region=region))
            continue

        ltot = local.hourly_average()

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
