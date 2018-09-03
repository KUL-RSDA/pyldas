
import os

import numpy as np
import pandas as pd

from pyldas.interface import LDAS_io
from pyldas.templates import template_error_Tb40

from myprojects.timeseries import calc_anomaly

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

    dims = io.timeseries['obs_obs'].shape

    obs_errstd = np.full(dims[0:-1], 4.)

    # ----- Calculate anomalies -----
    cnt = 0
    for spc in np.arange(dims[0]):
        for lat in np.arange(dims[1]):
            for lon in np.arange(dims[2]):
                cnt += 1
                print '%i / %i' % (cnt, np.prod(dims[0:-1]))

                try:
                    obs = calc_anomaly(io.timeseries['obs_obs'][spc, lat, lon, :].to_dataframe()['obs_obs'],
                                       method='moving_average', longterm=True)
                    fcst = calc_anomaly(io.timeseries['obs_fcst'][spc, lat, lon, :].to_dataframe()['obs_fcst'],
                                        method='moving_average', longterm=True)
                    fcst_errvar = np.nanmean(io.timeseries['obs_fcstvar'][spc, lat, lon, :].values)

                    tmp_obs_errstd = (((obs - fcst) ** 2).mean() - fcst_errvar) ** 0.5
                    if not np.isnan(tmp_obs_errstd):
                        obs_errstd[spc, lat, lon] = tmp_obs_errstd

                except:
                    pass

    np.place(obs_errstd, obs_errstd < 0, 0)
    np.place(obs_errstd, obs_errstd > 20, 20)

    # ----- write output files -----
    for orb in orbits:
        # !!! inconsistent with the definition in the obs_paramfile (species) !!!
        modes[0] = 1 if orb == 'A' else 0

        res = template.copy()

        spc = 0 if orb == 'A' else 1
        res['err_Tbh'] = obs_errstd[spc,ind_lat,ind_lon]

        spc = 2 if orb == 'A' else 3
        res['err_Tbv'] = obs_errstd[spc,ind_lat,ind_lon]

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