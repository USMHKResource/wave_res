import pyDictH5 as pdh5
import itertools as it 
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd 


'''
Calculating Functions
'''
def rePack(x):
    time = x['time']
    rRange = x['range']
 
    def Td(i,j,time):
        df = pd.DataFrame(j,index=time,columns=rRange)
        df.name = i
        return df 
    def Od(i,j):
        if len(j) == len(time):
            index = time
        else:
            index = None
        return pd.Series(j,name=i,index=index)
 
    dataDict = {i:Td(i,j,time) if j.ndim is 2 else Od(i,j) 
                 for i,j in x.items()}

    if '1way' in dataDict.keys():
        dataDict['oneway'] = dataDict.pop('1way',None)

    return dataDict
         
def integ(x,mfs=200,local=False):
    if local:
        integ = {i:j.loc[:,:mfs].sum(1) for i,j in x.items() if j.ndim is 2}
    else:
        integ = {i:j.loc[:,mfs] for i,j in x.items() if j.ndim is 2}
    
    integ['Nhour'] = x['Nhour']
    
    return integ


def wAverage(x, unit):
    _factordict = {
        'GW': 1e-9,
        'TWh/yr': 365 * 24 * 1e-12,
    }
    factor = _factordict[unit]
    weights = x['Nhour']

    return {i:np.average(j ,weights=weights) * factor 
                for i,j in x.items() if i != 'Nhour'}

def wAverage_grouping(x,unit,group='month'):
    _factordict = {
        'GW': 1e-9,
        'TWh/yr': 365 * 24 * 1e-12,
    }
    factor = _factordict[unit]
    gw = x['Nhour']

    wMean = lambda x,i: np.sum(x['Nhour']*x[i])/np.sum(x['Nhour'])
    
    result = []
    for i,j in x.items():
        if i != 'Nhour':
            j.name = i
            grp = pd.concat([j,gw],axis=1)
            gx = grp.groupby(getattr(grp.index,group)) 
            tempDf = gx.apply(wMean,(i))*factor
            tempDf.name = i
            result.append(tempDf)
    
    return pd.concat(result,axis=1)

'''
Plotting Functions
'''
def singleSource_allLocals_singleGroup(groupTotals,totals,method,
                                        unit = 'TWh/yr',
                                        groupMethod = 'month'):
    
    plt.figure(figsize = [8,3])
    for location,group in groupTotals[totals].items():
        plt.plot(group.index,group[method],label=location)
    plt.ylabel(f'Selected Total ({unit})')
    plt.xlabel(groupMethod)
    plt.title(f'{method} from {totals}')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.show()
    plt.cla()

def allSource_singleLocals_singleGroup(groupTotals,totals,local,
                                        unit = 'TWh/yr',
                                        groupMethod = 'month'):

    groupTotals[totals][local].plot(figsize = [8,3])
    plt.ylabel('Selected Total ('+unit+')')
    plt.xlabel(groupMethod)
    plt.title(f'All Source Terms from {local} from {totals}')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.show()
    plt.cla()

def singleSource_singleLocals_allGroup(groupTotals,sourceTerm,local,
                                        unit = 'TWh/yr',
                                        groupMethod = 'month'):

    plt.figure(figsize = [8,3])
    for group in groupTotals:
        plt.plot(group[local].index,group[local][sourceTerm])



if __name__ == "__main__":

    #unit = 'GW'
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
    
    partitionedData = {name:{region: integ(data[0][region],mfs=mfs,local=data[1]) 
                                for region in regions}
                                for name,data in loadedData.items()}

    totals = {name:pd.DataFrame({region: wAverage(data[region],unit) 
                    for region in regions}).transpose() 
                    for name,data in partitionedData.items()}
    
'''
print("")
print('Remote Totals ({})'.format(unit))
print(totals['rtotX'])
print("")
print('Local Baseline Totals ({})'.format(unit))
print(totals['ltot0'])
print("")
print('Local Potential Totals ({})'.format(unit))
print(totals['ltotX'])


# To save
rtot0.to_csv('/filename.csv')
ltot0.to_csv('/filename.csv')
rtotX.to_csv('/filename.csv')
ltotX.to_csv('/filename.csv')
'''

# Monthly Groupings using datetime indexing through Pandas
# little more complex but not much

#groupMethod = 'month'
#groupTotals = {name:{region:wAverage_grouping(data[region],unit,group=groupMethod) 
                    for region in regions}
                    for name,data in partitionedData.items()}

# A simple plotting routine

#totals, sourceTerm, local = 'ltotX', 'snl', 'ak'

#singleSource_allLocals_singleGroup(groupTotals,totals,sourceTerm)
#allSource_singleLocals_singleGroup(groupTotals,totals,local)
#singleSource_singleLocals_allGroup(groupTotals,sourceTerm,local)






