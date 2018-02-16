
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

    result_file = r'D:\work\LDAS\2018-02_scaling\diagnostics\filter_diagnostics.nc'

    no_DA_innov = LDAS_io('ObsFcstAna','US_M36_SMOS_noDA_unscaled')
    cal_pent_DA_innov = LDAS_io('ObsFcstAna','US_M36_SMOS_DA_calibrated_scaled')
    cal_harm_DA_innov = LDAS_io('ObsFcstAna','US_M36_SMOS_DA_calibrated_harmonic')
    uncal_pent_DA_innov = LDAS_io('ObsFcstAna','US_M36_SMOS_DA_nocal_scaled_pentadal')
    uncal_harm_DA_innov = LDAS_io('ObsFcstAna','US_M36_SMOS_DA_nocal_scaled_harmonic')

    cal_pent_DA_incr = LDAS_io('incr','US_M36_SMOS_DA_calibrated_scaled')
    cal_harm_DA_incr = LDAS_io('incr','US_M36_SMOS_DA_calibrated_harmonic')
    uncal_pent_DA_incr = LDAS_io('incr','US_M36_SMOS_DA_nocal_scaled_pentadal')
    uncal_harm_DA_incr = LDAS_io('incr','US_M36_SMOS_DA_nocal_scaled_harmonic')

    runs = OrderedDict([(1,[no_DA_innov.timeseries]),
                        (2,[cal_pent_DA_innov.timeseries, cal_pent_DA_incr.timeseries]),
                        (3,[cal_harm_DA_innov.timeseries, cal_harm_DA_incr.timeseries]),
                        (4,[uncal_pent_DA_innov.timeseries, uncal_pent_DA_incr.timeseries]),
                        (5,[uncal_harm_DA_innov.timeseries, uncal_harm_DA_incr.timeseries])])

    tags = ['innov_mean','innov_var',
            'norm_innov_mean','norm_innov_var',
            'n_valid_innov',
            'incr_catdef_mean','incr_catdef_var',
            'incr_rzexc_mean','incr_rzexc_var',
            'incr_srfexc_mean','incr_srfexc_var']

    lons = np.unique(no_DA_innov.tilecoord['com_lon'].values)
    lats = np.unique(no_DA_innov.tilecoord['com_lat'].values)[::-1]

    species = no_DA_innov.timeseries['species'].values

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

    result_file = r'D:\work\LDAS\2018-02_scaling\_new\ismn_eval\validation.csv'

    no_DA = LDAS_io('xhourly','US_M36_SMOS_noDA_unscaled')
    cal_pent_DA = LDAS_io('xhourly','US_M36_SMOS_DA_calibrated_scaled')
    cal_harm_DA = LDAS_io('xhourly','US_M36_SMOS_DA_calibrated_harmonic')
    uncal_pent_DA = LDAS_io('xhourly','US_M36_SMOS_DA_nocal_scaled_pentadal')
    uncal_harm_DA = LDAS_io('xhourly','US_M36_SMOS_DA_nocal_scaled_harmonic')

    ismn = ISMN_io(col_offs=no_DA.tilegrids.loc['domain','i_offg'],
                   row_offs=no_DA.tilegrids.loc['domain','j_offg'])

    runs = ['OL', 'DA_cal_pent','DA_cal_harm','DA_uncal_pent', 'DA_uncal_harm']
    tss = [no_DA.timeseries, cal_pent_DA.timeseries, cal_harm_DA.timeseries, uncal_pent_DA.timeseries, uncal_harm_DA.timeseries]

    variables = ['sm_surface','sm_rootzone','sm_profile']
    modes = ['mean','ma','harmonic']

    # ismn.list = ismn.list.iloc[280::]
    i = 0
    for meta, ts_insitu in ismn.iter_stations():

        i += 1
        print '%i/%i' % (i, len(ismn.list))

        res = pd.DataFrame(meta.copy()).transpose()
        col = meta.ease_col
        row = meta.ease_row

        for var in variables:
            for mode in modes:

                ts_ref = calc_anomaly(ts_insitu[var],mode).dropna()
                res['len_' + mode + '_' + var] = len(ts_ref)

                for run,ts_model in zip(runs,tss):

                    ind = (ts_model['snow_mass'][row, col].values == 0)&(ts_model['soil_temp_layer1'][row, col].values > 277.15)
                    ts_mod = ts_model[var][row, col].to_series().loc[ind]
                    ts_mod.index += pd.to_timedelta('2 hours')
                    # TODO: Make sure that time of netcdf file is correct!!

                    ts_mod = calc_anomaly(ts_mod,mode).dropna()

                    tmp = pd.DataFrame({1: ts_ref, 2: ts_mod}).dropna()
                    r,p = pearsonr(tmp[1],tmp[2])

                    res['corr_' + run +'_' + mode + '_' + var] = r if (r > 0) & (p < 0.01) else np.nan
                    res['ubrmsd_' + run +'_' + mode + '_' + var] = np.sqrt(((ts_ref-ts_mod)**2).mean())


        if (os.path.isfile(result_file) == False):
            res.to_csv(result_file, float_format='%0.4f')
        else:
            res.to_csv(result_file, float_format='%0.4f', mode='a', header=False)

if __name__=='__main__':
    insitu_evaluation()