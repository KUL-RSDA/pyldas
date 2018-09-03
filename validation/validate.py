
import os

import numpy as np
import pandas as pd

from collections import OrderedDict

from scipy.stats import pearsonr

from pyldas.interface import LDAS_io
from myprojects.readers.insitu import ISMN_io

from myprojects.timeseries import calc_anomaly

from netCDF4 import Dataset

def ncfile_init(fname, lats, lons, runs, species, tags):

    ds = Dataset(fname, 'w', 'NETCDF4')

    dims = ['lat','lon','run','species']
    dimvals = [lats, lons, runs, species]
    chunksizes = [len(lats), len(lons), 1, 1]
    dtypes = ['float32','float32','uint8', 'uint8']

    for dim, dimval, chunksize, dtype in zip(dims, dimvals, chunksizes, dtypes):
        ds.createDimension(dim, len(dimval))
        ds.createVariable(dim, dtype, dimensions=[dim], chunksizes=[chunksize], zlib=True)
        ds.variables[dim][:] = dimval

    for tag in tags:
        if tag.find('innov') != -1:
            ds.createVariable(tag, 'float32', dimensions=dims, chunksizes=chunksizes, fill_value=-9999., zlib=True)
        else:
            ds.createVariable(tag, 'float32', dimensions=dims[0:-1], chunksizes=chunksizes[0:-1], fill_value=-9999., zlib=True)

    return ds


def filter_diagnostics_evaluation():

    result_file = r'D:\work\LDAS\2018-06_rmse_uncertainty\filter_diagnostics.nc'

    cal_DA_clim_innov = LDAS_io('ObsFcstAna','US_M36_SMOS40_DA_cal_scaled')
    cal_DA_seas_innov = LDAS_io('ObsFcstAna','US_M36_SMOS40_DA_cal_scl_errfile')

    cal_DA_clim_incr = LDAS_io('incr','US_M36_SMOS40_DA_cal_scaled')
    cal_DA_seas_incr = LDAS_io('incr','US_M36_SMOS40_DA_cal_scl_errfile')

    runs = OrderedDict([(1,[cal_DA_clim_innov.timeseries, cal_DA_clim_incr.timeseries]),
                        (2,[cal_DA_seas_innov.timeseries, cal_DA_seas_incr.timeseries])])

    tags = ['innov_mean','innov_var',
            'norm_innov_mean','norm_innov_var',
            'n_valid_innov',
            'incr_catdef_mean','incr_catdef_var',
            'incr_rzexc_mean','incr_rzexc_var',
            'incr_srfexc_mean','incr_srfexc_var']

    lons = np.unique(cal_DA_clim_innov.tilecoord['com_lon'].values)
    lats = np.unique(cal_DA_clim_innov.tilecoord['com_lat'].values)[::-1]

    species = cal_DA_clim_innov.timeseries['species'].values

    ds = ncfile_init(result_file, lats, lons, runs.keys(), species, tags)

    for i_run,run in enumerate(runs):
        for i_spc,spc in enumerate(species):

            print 'run %i, species %i' % (i_run,i_spc)

            ds['innov_mean'][:,:,i_run,i_spc] = (runs[run][0]['obs_obs'][i_spc] - runs[run][0]['obs_fcst'][i_spc]).mean(dim='time').values
            ds['innov_var'][:,:,i_run,i_spc] = (runs[run][0]['obs_obs'][i_spc] - runs[run][0]['obs_fcst'][i_spc]).var(dim='time').values
            ds['norm_innov_mean'][:,:,i_run,i_spc] = ((runs[run][0]['obs_obs'][i_spc] - runs[run][0]['obs_fcst'][i_spc]) /
                                                      np.sqrt(runs[run][0]['obs_obsvar'][i_spc] + runs[run][0]['obs_fcstvar'][i_spc])).mean(dim='time').values
            ds['norm_innov_var'][:,:,i_run,i_spc] = ((runs[run][0]['obs_obs'][i_spc] - runs[run][0]['obs_fcst'][i_spc]) /
                                                     np.sqrt(runs[run][0]['obs_obsvar'][i_spc] + runs[run][0]['obs_fcstvar'][i_spc])).var(dim='time').values

            tmp = runs[run][0]['obs_obs'][i_spc].values
            np.place(tmp, ~np.isnan(tmp), 1.)
            np.place(tmp, np.isnan(tmp), 0.)
            ds['n_valid_innov'][:, :, i_run, i_spc] = tmp.sum(axis=2)

        if len(runs[run]) == 2:
            np.place(runs[run][1]['catdef'].values, runs[run][1]['catdef'].values == 0, np.nan)
            np.place(runs[run][1]['rzexc'].values, runs[run][1]['rzexc'].values == 0, np.nan)
            np.place(runs[run][1]['srfexc'].values, runs[run][1]['srfexc'].values == 0, np.nan)
            ds['incr_catdef_mean'][:, :, i_run] = runs[run][1]['catdef'].mean(dim='time').values
            ds['incr_catdef_var'][:, :, i_run] = runs[run][1]['catdef'].var(dim='time').values
            ds['incr_rzexc_mean'][:, :, i_run] = runs[run][1]['rzexc'].mean(dim='time').values
            ds['incr_rzexc_var'][:, :, i_run] = runs[run][1]['rzexc'].var(dim='time').values
            ds['incr_srfexc_mean'][:, :, i_run] = runs[run][1]['srfexc'].mean(dim='time').values
            ds['incr_srfexc_var'][:, :, i_run] = runs[run][1]['srfexc'].var(dim='time').values

    ds.close()


def insitu_evaluation():

    result_file = r'D:\work\LDAS\2018-06_rmse_uncertainty\insitu_evaluation\validation.csv'

    noDA = LDAS_io('xhourly', 'US_M36_SMOS40_noDA_cal_scaled')
    DA_const_err = LDAS_io('xhourly', 'US_M36_SMOS40_DA_cal_scaled')
    DA_varia_err = LDAS_io('xhourly', 'US_M36_SMOS40_DA_cal_scl_errfile')

    ismn = ISMN_io(col_offs=noDA.tilegrids.loc['domain','i_offg'],
                   row_offs=noDA.tilegrids.loc['domain','j_offg'])

    runs = ['noDA', 'DA_const_err','DA_varia_err']
    tss = [noDA.timeseries, DA_const_err.timeseries, DA_varia_err.timeseries]

    variables = ['sm_surface','sm_rootzone','sm_profile']
    modes = ['absolute','longterm','shortterm']

    # ismn.list = ismn.list.iloc[101::]

    i = 0
    for meta, ts_insitu in ismn.iter_stations():
        i += 1
        print '%i/%i' % (i, len(ismn.list))

        res = pd.DataFrame(meta.copy()).transpose()
        col = meta.ease_col
        row = meta.ease_row

        for var in variables:
            for mode in modes:

                if mode == 'absolute':
                    ts_ref = ts_insitu[var].dropna()
                elif mode == 'mean':
                    ts_ref = calc_anomaly(ts_insitu[var], mode).dropna()
                else:
                    ts_ref = calc_anomaly(ts_insitu[var], method='moving_average', longterm=(mode=='longterm')).dropna()

                res['len_' + mode + '_' + var] = len(ts_ref)

                for run,ts_model in zip(runs,tss):

                    ind = (ts_model['snow_mass'][row, col].values == 0)&(ts_model['soil_temp_layer1'][row, col].values > 277.15)
                    ts_mod = ts_model[var][row, col].to_series().loc[ind]
                    ts_mod.index += pd.to_timedelta('2 hours')
                    # TODO: Make sure that time of netcdf file is correct!!

                    if mode == 'absolute':
                        ts_mod = ts_mod.dropna()
                    else:
                        ts_mod = calc_anomaly(ts_mod, method='moving_average', longterm=mode=='longterm').dropna()

                    tmp = pd.DataFrame({1: ts_ref, 2: ts_mod}).dropna()
                    r,p = pearsonr(tmp[1],tmp[2])

                    res['corr_' + run +'_' + mode + '_' + var] = r if (r > 0) & (p < 0.01) else np.nan
                    res['ubrmsd_' + run +'_' + mode + '_' + var] = np.sqrt(((ts_ref-ts_mod)**2).mean())


        if (os.path.isfile(result_file) == False):
            res.to_csv(result_file, float_format='%0.4f')
        else:
            res.to_csv(result_file, float_format='%0.4f', mode='a', header=False)

if __name__=='__main__':
    # filter_diagnostics_evaluation()
    insitu_evaluation()