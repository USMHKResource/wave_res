#import cPickle as pkl 
import pickle as pkl 

import numpy as np
import pandas as pd
import cartopy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import shapely.geometry as sgeom






if __name__ == "__main__":

    
    from boundaries import *
    region = 'ak'
    data = run_all()
    grids = pkl.load(open('GridLonLat.pkl','r'))
    
    for region in data.keys():
        try:
            idx = data[region]['EEZ']
            lon,lat = grids[region][0][idx[0]],grids[region][1][idx[0]]
            df = pd.DataFrame(np.array([lon,lat]))
            df.to_csv('latlon_EEZ_{}.csv'.format(region))
        except KeyError:
            pass
    '''
    #with open('temp.pkl','rb') as temp:
    #    lon,lat = pkl.load(temp)
    
    #df = pd.read_csv('temp.csv')
    #lon,lat = df.iloc[:,0], df.iloc[:,1]

    #plt.figure()
    #plt.scatter(lon,lat)
    #plt.show()

    #track = sgeom.LineString(zip(lon, lat))

    fig = plt.figure(figsize=(8, 12))
    ax = plt.axes(projection=ccrs.Geostationary(),
                    #central_longitude=np.mean(lon),
                    #central_latitude=np.mean(lat))
                    )
                    
    ax.coastlines()

    plt.show()
    '''