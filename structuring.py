import pyDictH5 as pdh5
import itertools as it 
import numpy as np
import pandas as pd 
import xarray as xr


if __name__ == "__main__":
     unit = 'TWh/yr'

    mfs = 200

    regions = ['ak', 'at', 'prusvi', 'wc', 'hi'] 
    baseExt = ['baseline','extraction']
    remLoc = ['remote','local']

    files = [[[f'results/{be}/{region}.{rl}-totals.h5',region]
            for region in regions] for be in baseExt for rl in remLoc]
    
    loadedData = {'rtot0':[{region:rePack(pdh5.load(f)) 
                            for f,region in files[0]},False],
                    'ltot0':[{region:rePack(pdh5.load(f)) 
                            for f,region in files[1]},True],
                    'rtotX':[{region:rePack(pdh5.load(f)) 
                            for f,region in files[2]},False],
                    'ltotX':[{region:rePack(pdh5.load(f)) 
                            for f,region in files[3]},True]}