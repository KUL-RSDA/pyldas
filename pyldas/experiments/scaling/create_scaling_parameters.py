
import os

import timeit

import numpy as np
import pandas as pd
import xarray as xr

from pyldas.grids import EASE2
from pyldas.readers import LDAS_io
from pyldas.templates import template_scaling

from myprojects.timeseries import calc_clim, calc_pentadal_mean

def calc_clim_p(ts, n):

    clim = calc_clim(ts, n=n)

    pentads = np.floor((clim.index.values - 1) / 5.)
    clim = clim.groupby(pentads,axis=0).mean()
    clim.index = np.arange(73)+1

    return clim

def calc_vals_per_pent(ts):
    xts = ts.copy().dropna()
    doys = xts.copy().dropna().index.dayofyear.values
    doys[xts.index.is_leap_year & (doys>59)] -= 1
    pentads = np.floor((doys - 1) / 5.)

    res = pd.Series(index=np.arange(73)+1)
    for pent in np.arange(73)+1:
        res.loc[pent] = len(xts.loc[pentads==pent])

    return res.astype('int')

def run():
    froot = r"C:\Users\u0116961\Documents\VSC\vsc_data_copies\scratch_TEST_RUNS\US_M36_SMOS_noDA_unscaled\scaling_files"
    fbase = '7Thv_TbSM_001_SMOS_zscore_stats_2010_p37_2015_p36_hscale_0.00_W_9p_Nmin_20_'

    io = LDAS_io('ObsFcstAna', exp='US_M36_SMOS_noDA_unscaled')
    grid = EASE2()
    dtype = template_scaling()[0]
    n = 3

    tiles = io.tilecoord['tile_id'].values.astype('int32')
    pentads = np.arange(73)+1
    angles = np.array([30., 35., 40., 45., 50., 55., 60.])
    pols = ['V','H']
    orbits = ['A', 'D']

    template = pd.DataFrame(columns=dtype.names, index=tiles).astype('float32')
    template['lon'] = io.tilecoord['com_lon'].values.astype('float32')
    template['lat'] = io.tilecoord['com_lat'].values.astype('float32')
    template['tile_id'] = tiles.astype('int32')

    dummy = np.full([len(tiles),len(pentads),len(angles),len(pols),len(orbits)],-9999)

    coords = {'tile_id': tiles,
              'pentad': pentads,
              'angle': angles,
              'pol': pols,
              'orbit': orbits}

    darr = xr.DataArray(dummy, coords=coords, dims=['tile_id','pentad','angle','pol','orbit'])

    data = xr.Dataset({'m_obs_p':darr.astype('float32'),
                       'm_mod_p':darr.astype('float32'),
                       'm_obs_h':darr.astype('float32'),
                       'm_mod_h':darr.astype('float32'),
                       'N_data':darr.astype('int32')})

    # ----- calculate mean and reshuffle -----
    for i,til in enumerate(tiles):
        print '%i/%i' % (i, len(tiles))
        for pol in pols:
            for ang in angles:
                for orb in orbits:

                    spc = io.get_species(pol=pol, ang=ang, orbit=orb)
                    col,row = grid.tileid2colrow(til)

                    obs = io.timeseries['obs_obs'][spc-1,col,row].to_series()
                    mod = io.timeseries['obs_fcst'][spc-1,col,row].to_series()

                    data['m_obs_p'].sel(tile_id=til, pol=pol, angle=ang, orbit=orb)[:] = calc_pentadal_mean(obs).values
                    data['m_mod_p'].sel(tile_id=til, pol=pol, angle=ang, orbit=orb)[:] = calc_pentadal_mean(mod).values
                    data['m_obs_h'].sel(tile_id=til, pol=pol, angle=ang, orbit=orb)[:] = calc_clim_p(obs, n=n).values
                    data['m_mod_h'].sel(tile_id=til, pol=pol, angle=ang, orbit=orb)[:] = calc_clim_p(mod, n=n).values

                    s = timeit.default_timer()
                    data['N_data'].sel(tile_id=til, pol=pol, angle=ang, orbit=orb)[:] = calc_vals_per_pent(obs).values

    modes = np.array([1, 0])
    sdate = np.array([2010, 1, 1, 0, 0])
    edate = np.array([2015, 12, 31, 0, 0])
    lengths = np.array([len(tiles), len(angles), 1])  # tiles, incidence angles, whatever

    # ----- write output files -----
    for pent in pentads:
        for orb in orbits:

            # pentadal means
            res = template.copy()

            for ang in angles:
                for pol in pols:
                    tmp = data['m_obs_p'].sel(pol=pol, angle=ang, orbit=orb, pentad=pent).to_series()
                    res.loc[tmp.index, 'm_obs_' + pol + '_%i' % ang] = tmp
                    tmp = data['m_mod_p'].sel(pol=pol, angle=ang, orbit=orb, pentad=pent).to_series()
                    res.loc[tmp.index, 'm_mod_' + pol + '_%i' % ang] = tmp
                    tmp = data['N_data'].sel(pol=pol, angle=ang, orbit=orb, pentad=pent).to_series()
                    res.loc[tmp.index, 'N_data_' + pol + '_%i' % ang] = tmp

            res.replace(np.nan,-9999,inplace=True)

            fname = os.path.join(froot, 'pentadal_mean', fbase + orb + '_p%02i' % pent + '.bin')
            fid = open(fname, 'wb')
            io.write_fortran_block(fid, modes)
            io.write_fortran_block(fid, sdate)
            io.write_fortran_block(fid, edate)
            io.write_fortran_block(fid, lengths)
            io.write_fortran_block(fid, angles)
            for f in res.columns.values:
                io.write_fortran_block(fid, res[f].values)
            fid.close()


            # harmonic means
            res = template.copy()

            for ang in angles:
                for pol in pols:
                    tmp = data['m_obs_h'].sel(pol=pol, angle=ang, orbit=orb, pentad=pent).to_series()
                    res.loc[tmp.index, 'm_obs_' + pol + '_%i' % ang] = tmp
                    tmp = data['m_mod_h'].sel(pol=pol, angle=ang, orbit=orb, pentad=pent).to_series()
                    res.loc[tmp.index, 'm_mod_' + pol + '_%i' % ang] = tmp
                    tmp = data['N_data'].sel(pol=pol, angle=ang, orbit=orb, pentad=pent).to_series()
                    res.loc[tmp.index, 'N_data_' + pol + '_%i' % ang] = tmp

            res.replace(np.nan, -9999, inplace=True)

            fname = os.path.join(froot, 'harmonic_mean', fbase + orb + '_p%02i' % pent + '.bin')
            fid = open(fname, 'wb')
            io.write_fortran_block(fid, modes)
            io.write_fortran_block(fid, sdate)
            io.write_fortran_block(fid, edate)
            io.write_fortran_block(fid, lengths)
            io.write_fortran_block(fid, angles)
            for f in res.columns.values:
                io.write_fortran_block(fid, res[f].values)
            fid.close()

