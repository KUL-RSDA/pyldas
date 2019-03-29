
import os
import logging

import numpy as np
import pandas as pd
import xarray as xr

from pyldas.grids import EASE2
from pyldas.interface import LDAS_io
from myprojects.functions import find_files
from pyldas.templates import template_scaling

from myprojects.timeseries import calc_clim_harmonic, calc_pentadal_mean

logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

def calc_clim_p(ts, mode='pentadal', n=3):

    if mode == 'pentadal':
        clim = calc_pentadal_mean(ts)
    else:
        clim = calc_clim_harmonic(ts, n=n)
        pentads = np.floor((clim.index.values - 1) / 5.)
        clim = clim.groupby(pentads,axis=0).mean()
        clim.index = np.arange(73)+1

    return clim

def run(mode='annual'): # annual / longterm


    froot = r"D:\data_sets\LDAS_runs\US_M36_SMOS_noDA_cal_unscaled\obs_scaling_pentadal_annual"
    fbase = '7Thv_TbSM_001_SMOS_zscore_stats_2010_p37_2015_p36_hscale_0.00_W_9p_Nmin_20_'

    io = LDAS_io('ObsFcstAna', exp='US_M36_SMOS_noDA_cal_unscaled')
    dtype = template_scaling()[0]

    tiles = io.grid.tilecoord['tile_id'].values.astype('int32')
    angles = np.array([30., 35., 40., 45., 50., 55., 60.])
    pols = ['V','H']
    orbits = ['A', 'D']

    template = pd.DataFrame(columns=dtype.names, index=tiles).astype('float32')
    template['lon'] = io.grid.tilecoord['com_lon'].values.astype('float32')
    template['lat'] = io.grid.tilecoord['com_lat'].values.astype('float32')
    template['tile_id'] = tiles.astype('int32')

    pentads = np.arange(73)+1

    if mode == 'longterm':
        dummy = np.full([len(tiles),len(pentads),len(angles),len(pols),len(orbits)],-9999)
        coords = {'tile_id': tiles,
                  'pentad': pentads,
                  'angle': angles,
                  'pol': pols,
                  'orbit': orbits}
        darr = xr.DataArray(dummy, coords=coords, dims=['tile_id','pentad','angle','pol','orbit'])
    else:
        years = np.arange(2010, 2017)
        dummy = np.full([len(tiles), len(pentads), len(years), len(angles), len(pols), len(orbits)], -9999)
        coords = {'tile_id': tiles,
                  'pentad': pentads,
                  'year': years,
                  'angle': angles,
                  'pol': pols,
                  'orbit': orbits}
        darr = xr.DataArray(dummy, coords=coords, dims=['tile_id', 'pentad', 'year', 'angle', 'pol', 'orbit'])

    data = xr.Dataset({'m_obs':darr.astype('float32'),
                       'm_mod':darr.astype('float32'),
                       'N_data':darr.astype('int32')})

    # ----- calculate mean and reshuffle -----
    for i,til in enumerate(tiles):
        logging.info('%i/%i' % (i, len(tiles)))
        for pol in pols:
            for ang in angles:
                for orb in orbits:

                    spc = io.get_species(pol=pol, ang=ang, orbit=orb)
                    col,row = io.grid.tileid2colrow(til)

                    obs = io.timeseries['obs_obs'][spc-1,row,col].to_series()
                    mod = io.timeseries['obs_fcst'][spc-1,row,col].to_series()

                    if mode == 'longterm':
                        data['m_obs'].sel(tile_id=til, pol=pol, angle=ang, orbit=orb)[:] = calc_clim_p(obs).values
                        data['m_mod'].sel(tile_id=til, pol=pol, angle=ang, orbit=orb)[:] = calc_clim_p(mod).values
                        data['N_data'].sel(tile_id=til, pol=pol, angle=ang, orbit=orb)[:] = len(obs.dropna())
                    else:
                        for yr in years:
                            data['m_obs'].sel(tile_id=til, pol=pol, angle=ang, orbit=orb, year=yr)[:] = calc_clim_p(obs[obs.index.year==yr]).values
                            data['m_mod'].sel(tile_id=til, pol=pol, angle=ang, orbit=orb, year=yr)[:] = calc_clim_p(mod[obs.index.year==yr]).values
                            data['N_data'].sel(tile_id=til, pol=pol, angle=ang, orbit=orb, year=yr)[:] = len(obs[obs.index.year==yr].dropna())

    modes = np.array([0, 0])
    sdate = np.array([2010, 1, 1, 0, 0])
    edate = np.array([2016, 12, 31, 0, 0])
    lengths = np.array([len(tiles), len(angles), 1])  # tiles, incidence angles, whatever

    # ----- write output files -----
    for pent in pentads:
        for orb in orbits:
            # !!! inconsistent with the definition in the obs_paramfile (species) !!!
            modes[0] = 1 if orb == 'A' else 0

            if mode == 'longterm':
                res = template.copy()
                for ang in angles:
                    for pol in pols:
                        tmp = data['m_obs'].sel(pol=pol, angle=ang, orbit=orb, pentad=pent).to_series()
                        res.loc[tmp.index, 'm_obs_' + pol + '_%i' % ang] = tmp
                        res.loc[tmp.index, 's_obs_' + pol + '_%i' % ang] = tmp
                        tmp = data['m_mod'].sel(pol=pol, angle=ang, orbit=orb, pentad=pent).to_series()
                        res.loc[tmp.index, 'm_mod_' + pol + '_%i' % ang] = tmp
                        res.loc[tmp.index, 's_mod_' + pol + '_%i' % ang] = tmp
                        tmp = data['N_data'].sel(pol=pol, angle=ang, orbit=orb, pentad=pent).to_series()
                        res.loc[tmp.index, 'N_data_' + pol + '_%i' % ang] = tmp

                res.replace(np.nan, -9999, inplace=True)

                fname = os.path.join(froot, fbase + orb + '_p%02i' % pent + '.bin')
                fid = open(fname, 'wb')
                io.write_fortran_block(fid, modes)
                io.write_fortran_block(fid, sdate)
                io.write_fortran_block(fid, edate)
                io.write_fortran_block(fid, lengths)
                io.write_fortran_block(fid, angles)
                for f in res.columns.values:
                    io.write_fortran_block(fid, res[f].values)
                fid.close()
            else:
                for yr in years:
                    res = template.copy()
                    for ang in angles:
                        for pol in pols:
                            tmp = data['m_obs'].sel(pol=pol, angle=ang, orbit=orb, pentad=pent, year=yr).to_series()
                            res.loc[tmp.index, 'm_obs_' + pol + '_%i' % ang] = tmp
                            res.loc[tmp.index, 's_obs_' + pol + '_%i' % ang] = tmp
                            tmp = data['m_mod'].sel(pol=pol, angle=ang, orbit=orb, pentad=pent, year=yr).to_series()
                            res.loc[tmp.index, 'm_mod_' + pol + '_%i' % ang] = tmp
                            res.loc[tmp.index, 's_mod_' + pol + '_%i' % ang] = tmp
                            tmp = data['N_data'].sel(pol=pol, angle=ang, orbit=orb, pentad=pent,
                                                     year=yr).to_series()
                            res.loc[tmp.index, 'N_data_' + pol + '_%i' % ang] = tmp

                    res.replace(np.nan, -9999, inplace=True)

                    fname = os.path.join(froot, fbase + orb + '_p%02i' % pent + '_y%04i' % yr + '.bin')
                    fid = open(fname, 'wb')
                    io.write_fortran_block(fid, modes)
                    io.write_fortran_block(fid, sdate)
                    io.write_fortran_block(fid, edate)
                    io.write_fortran_block(fid, lengths)
                    io.write_fortran_block(fid, angles)
                    for f in res.columns.values:
                        io.write_fortran_block(fid, res[f].values)
                    fid.close()


def replace_orbit_field():

    root = r'C:\Users\u0116961\Documents\VSC\vsc_data_copies\scratch_TEST_RUNS\US_M36_SMOS_noDA_unscaled\obs_scaling'

    for f in find_files(root,'_D_p'):
        data = np.fromfile(f,'>i4')
        data[1] = 0
        data.tofile(f)


if __name__ == '__main__':
    run()