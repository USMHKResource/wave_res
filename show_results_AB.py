import pyDictH5 as pdh5
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd 
from collections import defaultdict

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
        'TWh/yr': 365*24*1e-12,
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
    
    return pd.concat(result,axis=1)


if __name__ == "__main__":

    regions = ['ak', 'at', 'prusvi', 'wc', 'hi'] 

    saveDir = 'preferedSaveDir/'
    #unit = 'GW'
    unit = 'TWh/yr'
    
    # Initial Loading
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
    

    # pick a distance and calculate 32 yr avg
    mfs = 200 #miles from shore
    dataAtMFS = {name:{region: integ(data[0][region],mfs=mfs,local=data[1]) 
                                for region in regions}
                                for name,data in loadedData.iteritems()}

    totals = {name:pd.DataFrame({region: wAverage(data[region],unit) 
                    for region in regions}).transpose() 
                    for name,data in dataAtMFS.iteritems()}
    

    # as a function of distance
    mfss = range(0,201)[::10][1:]
    all_mfs = {}
    for mfs in mfss:
        dataAtMFS = {name:{region: integ(data[0][region],mfs=mfs,local=data[1]) 
                                for region in regions}
                                for name,data in loadedData.iteritems()}

        all_mfs[mfs] = {name:pd.DataFrame({region: wAverage(data[region],unit) 
                    for region in regions}).transpose() 
                    for name,data in dataAtMFS.iteritems()}
    
    # Group as a function of distance
    group, source = 'ltot0', 'stot'
    regVal = defaultdict(list)
    for region in regions:    
        for mfs in mfss:
            regVal[region].append(all_mfs[mfs][group].loc[region,source])
    
    # plot the Grouped data
    fig = plt.figure(figsize = [8,3])
    for location,values in regVal.iteritems():
        plt.plot(mfss,values,label=location)

    plt.ylabel('Selected Total ('+unit+')')
    plt.xlabel('Distance from Shore (miles)')
    plt.title('Full Year Average '+source+' from '+group)
    plt.legend()
    plt.grid()
    plt.tight_layout()
    #plt.savefig(saveDir+'figname.png')
    plt.show()

    # Monthly Groupings using datetime indexing through Pandas
    # little more complex but not much

    mfs = 200 #miles from shore
    dataAtMFS = {name:{region: integ(data[0][region],mfs=mfs,local=data[1]) 
                                for region in regions}
                                for name,data in loadedData.iteritems()}

    groupMethod = 'month'
    groupTotals = {name:{region:wAverage_grouping(data[region],unit,group=groupMethod) 
                        for region in regions}
                        for name,data in dataAtMFS.iteritems()}

    # A simple plotting routine

    totals, source = 'ltotX', 'stot'
    plotGroup = groupTotals[totals]

    fig = plt.figure(figsize = [8,3])
    for location,group in plotGroup.iteritems():
        plt.plot(group.index,group[source],label=location)

    plt.ylabel('Selected Total ('+unit+')')
    plt.xlabel(groupMethod)
    plt.title('Monthly Averaged '+source+' from '+totals+' '+str(mfs)+' miles Offshore')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    #plt.savefig(saveDir+'figname.png')
    plt.show()

    # source term monthly averages as a function of distance from shore

    groupMethod = 'month'
    mfss = range(0,201)[::10][1:]
    all_mfs = {}
    for mfs in mfss:
        dataAtMFS = {name:{region: integ(data[0][region],mfs=mfs,local=data[1]) 
                                for region in regions}
                                for name,data in loadedData.iteritems()}


        all_mfs[mfs] = {name:{region:wAverage_grouping(data[region],unit,group=groupMethod) 
                        for region in regions}
                        for name,data in dataAtMFS.iteritems()}

    # Group as a function of distance based on the month
    group, source, month  = 'ltot0', 'stot', 2
    regVal = defaultdict(list)
    for region in regions:    
        for mfs in mfss:
            regVal[region].append(
                all_mfs[mfs][group][region].loc[month,source])

    # Plot source vs distance from shore for chosen month

    fig = plt.figure(figsize = [8,3])
    for location,values in regVal.iteritems():
        plt.plot(mfss,values,label=location)

    plt.ylabel('Selected Total ('+unit+')')
    plt.xlabel('Distance from Shore (miles)')
    plt.title(source+' from '+group+' in Month '+str(month))
    plt.legend()
    plt.grid()
    plt.tight_layout()
    #plt.savefig(saveDir+'figname.png')
    plt.show()

    # Surface plot 
    group, source, region  = 'ltot0', 'stot', 'hi'
    regVal = []
    for mfs in mfss:
        data = all_mfs[mfs][group][region][source]
        regVal.append(data.values)

    sourceSurface = np.vstack(regVal)
    fig, ax = plt.subplots(constrained_layout=True)
    cont = ax.contourf(range(1,13),mfss,sourceSurface,50)
    cbar = fig.colorbar(cont)
    cbar.ax.set_ylabel('Selected Source Total ('+unit+')')
    ax.set_ylabel('Distance from Shore (miles)')
    ax.set_xlabel('Month of the Year')
    ax.set_title('Monthly Average of '+source+' from '+group+' at '+region)
    ax.grid()
    #plt.savefig(saveDir+'figname.png')
    plt.show()



