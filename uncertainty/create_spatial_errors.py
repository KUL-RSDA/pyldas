
import os

import numpy as np
import pandas as pd

from pyldas.interface import LDAS_io
from pyldas.templates import template_error_Tb40

def run():

    exp = 'US_M36_SMOS40_noDA_cal_scaled'

    io = LDAS_io('ObsFcstAna', exp)

    froot = r"D:\data_sets\LDAS_runs" + "\\" +  exp + "\\obs_err"
    fbase = 'SMOS_fit_Tb_'

    dtype = template_error_Tb40()[0]

    angles = np.array([40.,])
    orbits = ['A', 'D']

    tiles = io.tilecoord['tile_id'].values.astype('int32')
    ind_lat = io.tilecoord.loc[:, 'j_indg'].values - io.tilegrids.loc['domain', 'j_offg']
    ind_lon = io.tilecoord.loc[:, 'i_indg'].values - io.tilegrids.loc['domain', 'i_offg']

    template = pd.DataFrame(columns=dtype.names, index=tiles).astype('float32')
    template['lon'] = io.tilecoord['com_lon'].values.astype('float32')
    template['lat'] = io.tilecoord['com_lat'].values.astype('float32')

    modes = np.array([0, 0])
    sdate = np.array([2010, 1, 1, 0, 0])
    edate = np.array([2016, 12, 31, 0, 0])
    lengths = np.array([len(tiles), len(angles)])  # tiles, incidence angles, whatever

    # ----- write output files -----
    for orb in orbits:
        # !!! inconsistent with the definition in the obs_paramfile (species) !!!
        modes[0] = 1 if orb == 'A' else 0

        res = template.copy()

        spc = 1 if orb == 'A' else 2 # H polarization
        fcst_errvar = io.timeseries.sel(species=spc).obs_fcstvar.mean('time')
        obs_errstd = np.sqrt(((io.timeseries.sel(species=spc).obs_fcst - io.timeseries.sel(species=spc).obs_obs) ** 2).mean('time') - fcst_errvar)
        np.place(obs_errstd.values, obs_errstd.values < 0, 0)
        np.place(obs_errstd.values, obs_errstd.values > 20, 20)
        res['err_Tbh'] = obs_errstd.values[ind_lat,ind_lon]

        spc = 3 if orb == 'A' else 4 # V polarization
        fcst_errvar = io.timeseries.sel(species=spc).obs_fcstvar.mean('time')
        obs_errstd = np.sqrt(((io.timeseries.sel(species=spc).obs_fcst - io.timeseries.sel(species=spc).obs_obs) ** 2).mean('time') - fcst_errvar)
        np.place(obs_errstd.values, obs_errstd.values < 0, 0)
        np.place(obs_errstd.values, obs_errstd.values > 20, 20)
        res['err_Tbv'] = obs_errstd.values[ind_lat,ind_lon]

        res.replace(np.nan, 6, inplace=True)

        fname = os.path.join(froot, fbase + orb + '.bin')

        fid = open(fname, 'wb')
        io.write_fortran_block(fid, modes)
        io.write_fortran_block(fid, sdate)
        io.write_fortran_block(fid, edate)
        io.write_fortran_block(fid, lengths)
        io.write_fortran_block(fid, angles)

        for f in res.columns.values:
            io.write_fortran_block(fid, res[f].values)
        fid.close()


if __name__ == '__main__':
    run()