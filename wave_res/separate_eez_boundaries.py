import numpy as np
import proj
import grid_data as gdat
import geopandas as gpd
import matplotlib.pyplot as plt

if __name__ == "__main__":

    base = 'C:\\Users\\abharath\\Documents\\coastlines\\GSHHS_shp\\c\\'
    eez = 'C:\\Users\\abharath\\Documents\\World_EEZ\\'
    usM = 'USMaritimeLimitsNBoundaries.shp'
    fname = 'GSHHS_c_L1.shp'
    eezs = 'eez_boundaries_v10.shp'

    fig, ax = plt.subplots(1, 1)
    a = gpd.read_file(base+fname)
    b = gpd.read_file(eez+eezs)
    c = gpd.read_file(eez+usM)
    a.plot(ax=ax)
    b.plot(ax=ax)
    c.plot(ax=ax,color='red')
    alpha = 0.9

    """
    con_defs = gdat.con_defs
    latlon = gdat.gridlonlat
    bounds = gdat.bounds
    
    #'''
    ### WC
    con_defs['wc']['eez'] = [range(34, 148), ]
    con_defs['wc']['brdr_MEX'] = [range(0, 35), ] 
    con_defs['wc']['brdr_CAN'] = [range(147, 176), ]

    ax.scatter(latlon['wc'][0,range(34, 148)],latlon['wc'][1,range(34, 148)],c='r',alpha=alpha)
    ax.scatter(latlon['wc'][0,range(0, 35)],latlon['wc'][1,range(0, 35)],c='g',alpha=alpha)
    ax.scatter(latlon['wc'][0,range(147, 176)],latlon['wc'][1,range(147, 176)],c='b',alpha=alpha)
    #'''

    #'''
    ## AT
    con_defs['at']['eez-N'] = [range(32, 147), ]
    con_defs['at']['eez-S'] = [range(227, 282), ]
    con_defs['at']['brdr_CAN'] = [range(0, 33), ]
    con_defs['at']['brdr_BHMS'] = [range(146,193)]
    con_defs['at']['brdr_CUBA'] = [range(193, 228), ]
    con_defs['at']['brdr_MEX'] = [range(281, 303), ]

    ax.scatter(latlon['at'][0,range(32, 147)],latlon['at'][1,range(32, 147)],c='r',alpha=alpha)
    ax.scatter(latlon['at'][0,range(227, 282)],latlon['at'][1,range(227, 282)],c='g',alpha=alpha)
    ax.scatter(latlon['at'][0,range(0, 33)],latlon['at'][1,range(0, 33)],c='b',alpha=alpha)
    ax.scatter(latlon['at'][0,range(193, 228)],latlon['at'][1,range(193, 228)],c='y',alpha=alpha)
    ax.scatter(latlon['at'][0,range(146, 193)],latlon['at'][1,range(146, 193)],c='k',alpha=alpha)
    ax.scatter(latlon['at'][0,range(281, 303)],latlon['at'][1,range(281, 303)],c='m',alpha=alpha)
    #'''
    #'''
    ### AK
    
    con_defs['ak']['eez-N'] = [range(439, 541), ]
    con_defs['ak']['eez-S'] = [range(48, 392), ]
    
    con_defs['ak']['brdr_CAN'] = [range(0, 49), ]
    con_defs['ak']['brdr_RUS-S'] = [range(391, 440), ]
    con_defs['ak']['brdr_RUS-N'] = [range(540, 616), ]
    

    ax.scatter(latlon['ak'][0,range(439, 541)],latlon['ak'][1,range(439, 541)],c='k',alpha=alpha)
    ax.scatter(latlon['ak'][0,range(48, 392)],latlon['ak'][1,range(48, 392)],c='g',alpha=alpha)
    ax.scatter(latlon['ak'][0,range(0, 49)],latlon['ak'][1,range(0, 49)],c='b',alpha=alpha)
    ax.scatter(latlon['ak'][0,range(391, 440)],latlon['ak'][1,range(391, 440)],c='r',alpha=alpha)
    ax.scatter(latlon['ak'][0,range(540, 616)],latlon['ak'][1,range(540, 616)],c='m',alpha=alpha)
    #'''
    #'''
    ### PRUSVI
    con_defs['prusvi']['eez'] = [range(0, 7), ]
    con_defs['prusvi']['brdr_DR'] = [range(74, 116), ]
    con_defs['prusvi']['brdr_VNZ'] = [range(43, 74), ]
    con_defs['prusvi']['brdr_BVI'] = [range(7, 43), ]
    
    ax.scatter(latlon['prusvi'][0,range(7, 43)],latlon['prusvi'][1,range(7, 43)],c='k',alpha=alpha)
    ax.scatter(latlon['prusvi'][0,range(0, 7)],latlon['prusvi'][1,range(0, 7)],c='g',alpha=alpha)
    ax.scatter(latlon['prusvi'][0,range(74,116)],latlon['prusvi'][1,range(74,116)],c='r',alpha=alpha)
    ax.scatter(latlon['prusvi'][0,range(43, 74)],latlon['prusvi'][1,range(43, 74)],c='b',alpha=alpha)
    #'''
    ### HI
    con_defs['hi']['eez'] = [np.r_[range(0,67,1),
                               np.array([3081, 2907, 2736, 2567, 2400, 2236,
                                         2075, 1607, 1456, 1221, 1077, 1078,
                                         1079, 1080, 1081, 1082, 1083, 1084,
                                         1085, 1086, 1087, 1088, 1090, 1091,
                                         950,  949,  948,  947,  946, 1068,
                                         1069, 1070, 1071, 1072, 1073, 1074,
                                         1075, 1076, 1220, 1455, 1760, 1916,
                                         2074, 2566, 2735, 3080]),
                               range(340,453,1),
                               np.array([0])].tolist(),
                         ]
    con_defs['hi']['brdr_MW'] = [np.array([3081, 2907, 2736, 2567, 2400, 2236,
                                       2075, 1607, 1456, 1221, 1077, 1078,
                                       1079, 1080, 1081, 1082, 1083, 1084,
                                       1085, 1086, 1087, 1088, 1090, 1091,
                                       950,  949,  948,  947,  946, 1068,
                                       1069, 1070, 1071, 1072, 1073, 1074,
                                       1075, 1076, 1220, 1455, 1760, 1916,
                                       2074, 2566, 2735, 3080]),]

    ax.scatter(latlon['hi'][0,con_defs['hi']['eez']],latlon['hi'][1,con_defs['hi']['eez']],c='k',alpha=alpha)
    ax.scatter(latlon['hi'][0,con_defs['hi']['brdr_MW']],latlon['hi'][1,con_defs['hi']['brdr_MW']],c='g',alpha=alpha)
    """
    plt.show()