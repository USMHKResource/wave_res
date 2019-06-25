import geopandas as gpd 
import matplotlib.pyplot as plt


if __name__ == "__main__":

    base = 'C:\\Users\\abharath\\Documents\\coastlines\\GSHHS_shp\\h\\'
    eez = 'C:\\Users\\abharath\\Documents\\World_EEZ\\'
    
    fname = 'GSHHS_h_L1.shp'
    eezs = 'eez_boundaries_v10.shp'

    fig, ax = plt.subplots(1, 1)

    a = gpd.read_file(base+fname)
    b = gpd.read_file(eez+eezs)
    a.plot(ax=ax)
    b.plot(ax=ax)
    plt.show()