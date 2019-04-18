import pyDictH5 as pdh5
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd 



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
                 for i,j in x.iteritems()}

    if '1way' in dataDict.keys():
        dataDict['oneway'] = dataDict.pop('1way',None)

    return dataDict
         
def integ(x,mfs=200,local=False):
    if local:
        integ = {i:j.loc[:,:mfs].sum(1) for i,j in x.iteritems() if j.ndim is 2}
    else:
        integ = {i:j.loc[:,mfs] for i,j in x.iteritems() if j.ndim is 2}
    
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
                for i,j in x.iteritems() if i != 'Nhour'}

def wAverage_grouping(x,unit,group='month'):
    _factordict = {
        'GW': 1e-9,
        'TWh/yr': 365 * 24 * 1e-12,
    }
    factor = _factordict[unit]
    gw = x['Nhour']

    wMean = lambda x,i: np.sum(x['Nhour']*x[i])/np.sum(x['Nhour'])
    
    result = []
    for i,j in x.iteritems():
        if i != 'Nhour':
            j.name = i
            grp = pd.concat([j,gw],axis=1)
            gx = grp.groupby(getattr(grp.index,group)) 
            tempDf = gx.apply(wMean,(i))*factor
            tempDf.name = i
            result.append(tempDf)
    
    result = pd.concat(result,axis=1)
    
    return result
            


if __name__ == "__main__":

    regions = ['ak', 'at', 'prusvi', 'wc', 'hi'] 

    #unit = 'GW'
    unit = 'TWh/yr'

    mfs = 200
    
    remote0 = {region:rePack(pdh5.load('results/{}/{}.remote-totals.h5'
                                        .format('baseline', region)))
                                        for region in regions}
    local0 = {region:rePack(pdh5.load('results/{}/{}.local-totals.h5'
                                        .format('baseline', region)))
                                        for region in regions}
    remoteX = {region:rePack(pdh5.load('results/{}/{}.remote-totals.h5'
                                        .format('extraction', region)))
                                        for region in regions}                                    
    localX = {region:rePack(pdh5.load('results/{}/{}.local-totals.h5'
                                        .format('extraction', region)))
                                        for region in regions}
    
    loadedData = {'rtot0':[remote0,False],'ltot0':[local0,True],
                    'rtotX':[remoteX,False],'ltotX':[localX,True]}
    
    partitionedData = {name:{region: integ(data[0][region],mfs=mfs,local=data[1]) 
                                for region in regions}
                                for name,data in loadedData.iteritems()}

    totals = {name:pd.DataFrame({region: wAverage(data[region],unit) 
                    for region in regions}).transpose() 
                    for name,data in partitionedData.iteritems()}
    

print("")
print('Remote Totals ({})'.format(unit))
print(totals['rtotX'])
print("")
print('Local Baseline Totals ({})'.format(unit))
print(totals['ltot0'])
print("")
print('Local Potential Totals ({})'.format(unit))
print(totals['ltotX'])

'''
# To save
rtot0.to_csv('/filename.csv')
ltot0.to_csv('/filename.csv')
rtotX.to_csv('/filename.csv')
ltotX.to_csv('/filename.csv')
'''

# Monthly Groupings using datetime indexing through Pandas
# little more complex but not much

groupMethod = 'month'
groupTotals = {name:{region:wAverage_grouping(data[region],unit,group=groupMethod) 
                    for region in regions}
                    for name,data in partitionedData.iteritems()}

# A simple plotting routine

totals, method = 'ltotX', 'sds'
plotGroup = groupTotals[totals]

fig = plt.figure(figsize = [8,3])
for location,group in plotGroup.iteritems():
    plt.plot(group.index,group[method],label=location)

plt.ylabel('Selected Total ('+unit+')')
plt.xlabel(groupMethod)
plt.title(method+' from '+totals)
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()







