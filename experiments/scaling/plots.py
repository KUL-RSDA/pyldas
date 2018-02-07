

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib import gridspec

from pyldas.readers import LDAS_io

from myprojects.timeseries import calc_clim, calc_pentadal_mean

def plot_Tb_clims():

    ts = LDAS_io('ObsFcstAna', exp='US_M36_SMOS_noDA_unscaled').timeseries
    data = pd.DataFrame(index=ts.time)
    data['obs'] = pd.to_numeric(ts['obs_obs'][8,45,30].values,errors='coerce')
    data['fcst'] = pd.to_numeric(ts['obs_fcst'][8,45,30].values,errors='coerce')

    data_plt = data.copy().resample('1d').first()
    data_plt.interpolate(inplace=True)

    plt.figure(figsize=(18, 9))
    gs = gridspec.GridSpec(2, 2, width_ratios=[1, 3])

    ax = plt.subplot(gs[1])
    arr = data_plt['obs']
    yrs = np.unique(arr.index.year)
    df = pd.DataFrame(columns=yrs,index=np.arange(366)+1,dtype='float')
    for yr in yrs:
        df[yr].loc[arr[arr.index.year == yr].index.dayofyear] = arr[arr.index.year == yr].values
        df[yr] = df[yr]
    df.drop(366,axis='index',inplace=True)
    df.plot(ax=ax,linestyle='--',linewidth=0.5)
    arr = data['obs'].copy()
    clim = pd.DataFrame(index=np.arange(365)+1)
    for n in [0,1,2,3,10,]:
        clim[n] = calc_clim(arr,n=n)
    clim['pentad_mean'] = calc_pentadal_mean(arr)[0]
    clim.plot(linewidth=2,ax=ax)
    ax.set_xlim((0,365))
    ax.set_ylim((180,285))
    ax.legend(loc=2)

    ax = plt.subplot(gs[0])
    doys = arr.index.dayofyear.values
    modes = clim.columns.values
    rmse = pd.Series(index=modes)
    rmse[:] = 0.
    for mod in modes:
        n = 0
        for doy in np.arange(365)+1:
            n += len(arr[doys==doy])
            rmse[mod] += ((arr[doys==doy] - clim.loc[doy,mod])**2).sum()
        rmse[mod] = np.sqrt(rmse[mod]/n)
    ax = rmse.plot(ax=ax)
    ax.set_xticks(np.arange(len(modes)+1))
    ax.set_xticklabels(modes)

    ax = plt.subplot(gs[3])
    arr = data_plt['fcst']
    yrs = np.unique(arr.index.year)
    df = pd.DataFrame(columns=yrs,index=np.arange(366)+1,dtype='float')
    for yr in yrs:
        df[yr].loc[arr[arr.index.year == yr].index.dayofyear] = arr[arr.index.year == yr].values
        df[yr] = df[yr]
    df.drop(366,axis='index',inplace=True)
    df.plot(ax=ax,linestyle='--',linewidth=0.5)
    arr = data['fcst'].copy()
    clim = pd.DataFrame(index=np.arange(365)+1)
    for n in [0,1,2,3,10,]:
        clim[n] = calc_clim(arr,n=n)
    clim['pentad_mean'] = calc_pentadal_mean(arr)[0]
    clim.plot(linewidth=2,ax=ax)
    ax.set_xlim((0,365))
    ax.set_ylim((190,250))
    ax.legend(loc=2)

    ax = plt.subplot(gs[2])
    doys = arr.index.dayofyear.values
    modes = clim.columns.values
    rmse = pd.Series(index=modes)
    rmse[:] = 0.
    for mod in modes:
        n = 0
        for doy in np.arange(365)+1:
            n += len(arr[doys==doy])
            rmse[mod] += ((arr[doys==doy] - clim.loc[doy,mod])**2).sum()
        rmse[mod] = np.sqrt(rmse[mod]/n)
    ax = rmse.plot(ax=ax)
    ax.set_xticks(np.arange(len(modes)+1))
    ax.set_xticklabels(modes)

    plt.tight_layout()
    plt.show()

if __name__=='__main__':
    plot_Tb_clims()