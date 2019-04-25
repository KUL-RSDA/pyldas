
import os

import numpy as np
import pandas as pd
import xarray as xr

import platform
if platform.system() in ['Linux', 'Darwin']:
    import matplotlib
    matplotlib.use("TkAgg")

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

from pyldas.interface import LDAS_io
from myprojects.readers.insitu import ISMN_io
from myprojects.timeseries import calc_anomaly

from scipy.stats import pearsonr, spearmanr

def scatterplot_RTMparam_incr_innov_diff():

    outpath = r'D:\work\LDAS\2018-02_scaling\_new\diagnostics'

    diag = xr.open_dataset(r"D:\work\LDAS\2018-02_scaling\_new\diagnostics\filter_diagnostics.nc")

    params_cal = LDAS_io(exp='US_M36_SMOS_DA_calibrated_scaled').read_params('RTMparam')
    params_uncal = LDAS_io(exp='US_M36_SMOS_DA_nocal_scaled_harmonic').read_params('RTMparam')

    tc = LDAS_io().grid.tilecoord
    tg = LDAS_io().grid.tilegrids

    tc.i_indg -= tg.loc['domain', 'i_offg']  # col / lon
    tc.j_indg -= tg.loc['domain', 'j_offg']  # row / lat

    ind_lon,ind_lat = tc.i_indg.values, tc.j_indg.values

    ind_cal = 1
    ind_uncal = 3

    fontsize = 20

    incin = 'incr_'
    mv = 'mean'

    modes = ['catdef_','srfexc_','rzexc_']
    params = ['bh','bv','omega','rgh_hmin','rgh_hmax']


    plt.figure(figsize=(18,9))

    i=0
    for r,mode in enumerate(modes):
        for c,param in enumerate(params):
            i += 1
            ax = plt.subplot(3,5,i)

            # xdiff = (params_cal[param] - params_uncal[param]).values
            # ydiff = diag[incin+mode+mv][:,:,ind_cal].values[ind_lat,ind_lon] - diag[incin+mode+mv][:,:,ind_uncal].values[ind_lat,ind_lon]

            xdiff = params_uncal[param].values
            ydiff = diag[incin+mode+mv][:,:,ind_uncal].values[ind_lat,ind_lon]

            ind_valid = (~np.isnan(xdiff))&(~np.isnan(ydiff))

            xdiff = xdiff[ind_valid]
            ydiff = ydiff[ind_valid]

            s = np.argsort(xdiff)
            xdiff = xdiff[s]
            ydiff = ydiff[s]

            ax.plot(xdiff,ydiff,'o',markersize=3,markerfacecolor='k',markeredgecolor='k')

            fit = np.polyfit(xdiff, ydiff, deg=2)

            ax.plot(xdiff, fit[0] * xdiff**2 + fit[1] * xdiff + fit[2], color='red')

            if param == 'bh':
                ax.set_ylabel(mode[0:-1])
            if mode == 'rzexc_':
                ax.set_xlabel(param)

            corr = pearsonr(xdiff,ydiff)[0]
            rho = spearmanr(xdiff,ydiff)[0]

            ax.set_title('R = %.2f, $\\rho$ = %.2f' % (corr,rho))


    plt.tight_layout()
    plt.show()


def plot_cat_timeseries():

    outpath = r'D:\work\LDAS\2018-02_scaling\_new\ismn_eval\timeseries'

    fname = r"D:\work\LDAS\2018-02_scaling\_new\ismn_eval\validation.csv"
    res = pd.read_csv(fname)

    diff_srf = res['corr_DA_cal_pent_ma_sm_surface'] - res['corr_DA_uncal_pent_ma_sm_surface']
    diff_rz = res['corr_DA_cal_pent_ma_sm_rootzone'] - res['corr_DA_uncal_pent_ma_sm_rootzone']
    diff_prof = res['corr_DA_cal_pent_ma_sm_profile'] - res['corr_DA_uncal_pent_ma_sm_profile']
    ind = (diff_srf > 0.2) | (diff_rz > 0.2) | (diff_prof > 0.2)
    res = res.loc[ind,['network','station','lat', 'lon']]

    ismn = ISMN_io()
    cal = LDAS_io('xhourly','US_M36_SMOS_DA_calibrated_scaled')
    uncal = LDAS_io('xhourly','US_M36_SMOS_DA_nocal_scaled_pentadal')

    variables = ['sm_surface','sm_rootzone','sm_profile']

    for idx,stat in res.iterrows():

        fname = os.path.join(outpath,stat.network+'_'+stat.station+'.png')

        ts_ismn = ismn.read(stat.network,stat.station)
        lat = stat.lat
        lon = stat.lon

        plt.figure(figsize=(17,9))

        for i,var in enumerate(variables):

            ax = plt.subplot(3,1,i+1)

            ts_cal = calc_anomaly(cal.read_ts(var, lon, lat),method='ma')
            ts_cal.index += pd.to_timedelta('2 hours')
            ts_uncal = calc_anomaly(uncal.read_ts(var, lon, lat),method='ma')
            ts_uncal.index += pd.to_timedelta('2 hours')

            df = pd.DataFrame({'cal' : ts_cal,'uncal' : ts_uncal, 'insitu': calc_anomaly(ts_ismn[var],method='ma')}).dropna()
            if len(df) > 0:
                df.plot(ax=ax)
            else:
                continue

            title = 'R(ismn - cal) = %.2f , R(ismn - uncal) = %.2f' % (df.corr().loc['insitu','cal'], df.corr().loc['insitu','uncal'])

            ax.set_title(title,fontsize=12)
            ax.set_xlim('2010-01-01','2016-01-01')
            ax.set_ylim(-0.3,0.3)
            ax.set_xlabel('')

        plt.tight_layout()

        plt.savefig(fname, dpi=150)
        plt.close()


def plot_P_R_estimates():

    exp = 'US_M36_SMOS_DA_nocal_scaled_pentadal'

    io = LDAS_io('ObsFcstAna',exp)

    outpath = r"D:\work\LDAS\2018-02_scaling\uncertainty_estimates"

    lons = io.timeseries.lon.values
    lats = io.timeseries.lat.values

    lons, lats = np.meshgrid(lons, lats)
    llcrnrlat = 24
    urcrnrlat = 51
    llcrnrlon = -128
    urcrnrlon = -64

    cmap = 'jet'

    fontsize = 16

    cbrange = (0, 60)

    for spc in np.arange(len(io.timeseries.species)):

        R_est = ((io.timeseries['obs_obs'][spc, :, :, :] - io.timeseries['obs_fcst'][spc, :, :, :]) *
                 (io.timeseries['obs_obs'][spc, :, :, :] - io.timeseries['obs_ana'][spc, :, :, :])).mean(dim='time')

        P_est = ((io.timeseries['obs_obs'][spc, :, :, :] - io.timeseries['obs_fcst'][spc, :, :, :]) * (io.timeseries['obs_ana'][spc, :, :, :] - io.timeseries['obs_fcst'][spc, :, :, :])).mean(dim='time')

        f = plt.figure(figsize=(16, 5))

        ax = plt.subplot(121)
        m = Basemap(projection='mill', llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat, llcrnrlon=llcrnrlon,
                    urcrnrlon=urcrnrlon, resolution='c')
        m.drawcoastlines(); m.drawcountries(); m.drawstates()
        plt_img = np.ma.masked_invalid(R_est)
        im = m.pcolormesh(lons, lats, plt_img, cmap=cmap, latlon=True)
        im.set_clim(vmin=cbrange[0], vmax=cbrange[1])
        cb = m.colorbar(im, "bottom", size="7%", pad="8%")
        for t in cb.ax.get_xticklabels():
            t.set_fontsize(fontsize)
        for t in cb.ax.get_yticklabels():
            t.set_fontsize(fontsize)
        ax.set_title('R', fontsize=fontsize)

        ax = plt.subplot(122)
        m = Basemap(projection='mill', llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat, llcrnrlon=llcrnrlon,
                    urcrnrlon=urcrnrlon, resolution='c')
        m.drawcoastlines(); m.drawcountries(); m.drawstates()
        plt_img = np.ma.masked_invalid(P_est)
        im = m.pcolormesh(lons, lats, plt_img, cmap=cmap, latlon=True)
        im.set_clim(vmin=cbrange[0], vmax=cbrange[1])
        cb = m.colorbar(im, "bottom", size="7%", pad="8%")
        for t in cb.ax.get_xticklabels():
            t.set_fontsize(fontsize)
        for t in cb.ax.get_yticklabels():
            t.set_fontsize(fontsize)
        ax.set_title('HPH$^T$', fontsize=fontsize)

        plt.tight_layout()

        fname = os.path.join(outpath, 'R_P_est_%i.png' % (spc+1))

        plt.savefig(fname, dpi=f.dpi)
        plt.close()

        # plt.show()


def plot_filter_innov_diff():

    ds = xr.open_dataset(r"D:\work\LDAS\2018-02_scaling\diagnostics\filter_diagnostics.nc")

    # outpath = r"D:\work\LDAS\2018-02_scaling\diagnostics"

    lons = ds.lon.values
    lats = ds.lat.values

    lons, lats = np.meshgrid(lons, lats)

    imgs = [(ds['norm_innov_mean'][:, :, 1, :]-ds['norm_innov_mean'][:, :, 2, :]).mean(dim='species').values,
            (ds['incr_catdef_mean'][:, :, 1] - ds['incr_catdef_mean'][:, :, 2]).values,
            (ds['incr_srfexc_mean'][:, :, 1] - ds['incr_srfexc_mean'][:, :, 2]).values,
            (ds['norm_innov_var'][:, :, 1, :]-ds['norm_innov_var'][:, :, 2, :]).mean(dim='species').values,
            (ds['incr_catdef_var'][:, :, 1] - ds['incr_catdef_var'][:, :, 2]).values,
            (ds['incr_srfexc_var'][:, :, 1] - ds['incr_srfexc_var'][:, :, 2]).values]

    names = ['innov mean (ref-pent)',
             'incr catdef mean (ref-pent)',
             'incr srfexc mean (ref-pent)',
             'innov var (ref-pent)',
             'incr catdef var (ref-pent)',
             'incr srfexc var (ref-pent)']

    llcrnrlat = 24
    urcrnrlat = 51
    llcrnrlon = -128
    urcrnrlon = -64

    cmap = 'RdYlBu'

    fontsize = 16

    f = plt.figure(figsize=(19, 9))

    i = 0
    for img, name in zip(imgs, names):
        i += 1

        perc = np.nanpercentile(img, [1, 99])

        np.place(img, (img < perc[0]) | (img > perc[1]), np.nan)
        img = img.reshape(lons.shape)

        ax = plt.subplot(2, 3, i)

        if name.find('var') != -1:
            if name.find('catdef') != -1:
                cbrange = (-40, 40)
            else:
                cbrange = (-0.8, 0.8)
        else:
            if name.find('catdef') != -1:
                cbrange = (-0.8, 0.8)
            else:
                cbrange = (-0.15, 0.15)

        m = Basemap(projection='mill', llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat, llcrnrlon=llcrnrlon,
                    urcrnrlon=urcrnrlon,
                    resolution='c')

        m.drawcoastlines()
        m.drawcountries()
        m.drawstates()

        plt_img = np.ma.masked_invalid(img)
        im = m.pcolormesh(lons, lats, plt_img, cmap=cmap, latlon=True)

        im.set_clim(vmin=cbrange[0], vmax=cbrange[1])

        cb = m.colorbar(im, "bottom", size="7%", pad="8%")

        for t in cb.ax.get_xticklabels():
            t.set_fontsize(fontsize)
        for t in cb.ax.get_yticklabels():
            t.set_fontsize(fontsize)

        ax.set_title(name, fontsize=fontsize)

    plt.tight_layout()
    plt.show()


def plot_filter_increments(mode='image', var='catdef'):

    ds = xr.open_dataset(r"D:\work\LDAS\2018-06_rmse_uncertainty\filter_diagnostics.nc")

    lons = ds.lon.values
    lats = ds.lat.values

    lons,lats = np.meshgrid(lons,lats)

    imgs = [ds['incr_'+var+'_mean'][:, :, 0].values,
            ds['incr_'+var+'_mean'][:, :, 1].values,
            ds['incr_'+var+'_var'][:, :, 0].values,
            ds['incr_'+var+'_var'][:, :, 1].values]

    names = [var + ' mean (constant error)',
             var + ' mean (variable error)',
             var + ' var (constant error)',
             var + ' var (variable error)']

    f = plt.figure(figsize=(15,9))

    llcrnrlat = 24
    urcrnrlat = 51
    llcrnrlon = -128
    urcrnrlon = -64

    cmap = 'jet'

    fontsize = 16

    i = 0
    for img,name in zip(imgs,names):
        i += 1

        perc = np.nanpercentile(img, [1,99])

        np.place(img,(img<perc[0])|(img>perc[1]),np.nan)

        ax = plt.subplot(2,2,i)

        if mode == 'image':

            cbrange = (-0.8, 0.8) if name.find('mean') != -1 else (0,100)

            m = Basemap(projection='mill', llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat, llcrnrlon=llcrnrlon,
                        urcrnrlon=urcrnrlon,
                        resolution='c')
            m.drawcoastlines()
            m.drawcountries()
            m.drawstates()

            plt_img = np.ma.masked_invalid(img)
            im = m.pcolormesh(lons, lats, plt_img, cmap=cmap, latlon=True)

            im.set_clim(vmin=cbrange[0], vmax=cbrange[1])

            cb = m.colorbar(im, "bottom", size="7%", pad="8%")

            for t in cb.ax.get_xticklabels():
                t.set_fontsize(fontsize)
            for t in cb.ax.get_yticklabels():
                t.set_fontsize(fontsize)

        elif mode == 'hist':

            if name.find('mean') != -1:
                cbrange = (-0.5, 0.5)
                ylim=(0,6)
            else :
                cbrange = (0.2,2)
                ylim=(0,2)

            pd.Series(img.flatten()).dropna().hist(bins=30, range=cbrange, density=True, ax=ax)
            ax.set_ylim(ylim)

            x = cbrange[1] - (cbrange[1] - cbrange[0]) / 4.
            y = ylim[1] - (ylim[1] - ylim[0]) / 10.

            plt.text(x,y,'Mean = %.2f' % pd.Series(img.flatten()).mean(),fontsize=fontsize)

        ax.set_title(name, fontsize=fontsize)


    plt.tight_layout()
    plt.show()




def plot_filter_diagnostics(mode='image',spc_mean=True):

    ds = xr.open_dataset(r"D:\work\LDAS\2018-06_rmse_uncertainty\filter_diagnostics.nc")
    outpath = r"D:\work\LDAS\2018-06_rmse_uncertainty"

    lons = ds.lon.values
    lats = ds.lat.values

    lons,lats = np.meshgrid(lons,lats)


    if spc_mean is True:

        if mode == 'hist':

            imgs = [ds['norm_innov_mean'][:, :, 0, :].values,
                    ds['norm_innov_mean'][:, :, 1, :].values,
                    ds['norm_innov_var'][:, :, 0, :].values,
                    ds['norm_innov_var'][:, :, 1, :].values]

        else:

            imgs = [ds['norm_innov_mean'][:, :, 0, :].mean(dim='species').values,
                    ds['norm_innov_mean'][:, :, 1, :].mean(dim='species').values,
                    ds['norm_innov_var'][:, :, 0, :].mean(dim='species').values,
                    ds['norm_innov_var'][:, :, 1, :].mean(dim='species').values]

    names = ['Norm. innov mean (constant error)',
             'Norm. innov mean (variable error)',
             'Norm. innov var (constant error)',
             'Norm. innov var (variable error)']

    llcrnrlat = 24
    urcrnrlat = 51
    llcrnrlon = -128
    urcrnrlon = -64

    cmap = 'RdYlBu'

    fontsize = 16

    n_spc = 1 if spc_mean is True else len(ds.species)

    for j in np.arange(n_spc):

        if spc_mean is False:
            imgs = [ds['norm_innov_mean'][:, :, 0, j].values,
                    ds['norm_innov_mean'][:, :, 1, j].values,
                    ds['norm_innov_var'][:, :, 0, j].values,
                    ds['norm_innov_var'][:, :, 1, j].values]

        f = plt.figure(figsize=(15,9))

        i = 0
        for img,name in zip(imgs,names):
            i += 1

            # perc = np.nanpercentile(img, [1,99])
            #
            # np.place(img,(img<perc[0])|(img>perc[1]),np.nan)


            ax = plt.subplot(2,2,i)

            if mode == 'image':

                img = img.reshape(lons.shape)

                if name.find('mean') != -1 :
                    cbrange = (-0.6, 0.6)
                else:
                    if i == 3:
                        cbrange = (-8, 10)
                    else:
                        cbrange = (0, 2)

                m = Basemap(projection='mill', llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat, llcrnrlon=llcrnrlon,
                            urcrnrlon=urcrnrlon,
                            resolution='c')
                m.drawcoastlines()
                m.drawcountries()
                m.drawstates()

                plt_img = np.ma.masked_invalid(img)
                im = m.pcolormesh(lons, lats, plt_img, cmap=cmap, latlon=True)

                im.set_clim(vmin=cbrange[0], vmax=cbrange[1])

                cb = m.colorbar(im, "bottom", size="7%", pad="8%")

                for t in cb.ax.get_xticklabels():
                    t.set_fontsize(fontsize)
                for t in cb.ax.get_yticklabels():
                    t.set_fontsize(fontsize)

            elif mode == 'hist':

                if name.find('mean') != -1:
                    cbrange = (-1, 1)
                    ylim=(0,2)
                else:
                    cbrange = (-0.5,2.5)
                    ylim=(0,2)

                pd.Series(img.flatten()).dropna().hist(bins=30, range=cbrange, density=True, ax=ax)
                ax.set_ylim(ylim)

                x = cbrange[1] - (cbrange[1] - cbrange[0]) / 4.
                y = ylim[1] - (ylim[1] - ylim[0]) / 10.

                plt.text(x,y,'Mean = %.2f' % pd.Series(img.flatten()).mean(),fontsize=fontsize)

            ax.set_title(name, fontsize=fontsize)

        plt.tight_layout()

        # if spc_mean is True:
        #     fname = os.path.join(outpath, 'innov_species_mean_'+ mode + '.png')
        # else:
        #     fname = os.path.join(outpath, 'norm_innov_species_mean_'+ mode + '%i.png'%(j+1))
        #
        # plt.savefig(fname, dpi=f.dpi)
        # plt.close()

        plt.show()

def spatial_plot_ismn_stats():
    fname = r"C:\Users\u0116961\Documents\work\ldas_ismn_eval\2018-02_scaling\validation_masked.csv"

    res = pd.read_csv(fname)

    lats = res['lat'].values
    lons = res['lon'].values

    # variables = ['sm_surface','sm_rootzone','sm_profile']
    # runs = ['OL', 'DA_ref', 'DA_pent', 'DA_harm']
    # modes = ['mean', 'ma', 'harmonic']

    layer = 'rootzone'

    ts1 = res['corr_DA_ref_harmonic_sm_' + layer].values - res['corr_OL_harmonic_sm_' + layer].values
    ts2 = res['corr_DA_ref_harmonic_sm_' + layer].values - res['corr_DA_pent_harmonic_sm_' + layer].values
    ts3 = res['corr_DA_ref_harmonic_sm_' + layer].values - res['corr_DA_harm_harmonic_sm_' + layer].values
    ts4 = res['corr_DA_pent_harmonic_sm_' + layer].values - res['corr_DA_harm_harmonic_sm_' + layer].values
    title1 = 'ref - OL'
    title2 = 'ref - pent'
    title3 = 'ref - harm'
    title4 = 'pent - harm'

    # ts1 = res['corr_OL_ma_sm_' + layer].values - res['corr_OL_harmonic_sm_' + layer].values
    # ts2 = res['corr_DA_ref_ma_sm_' + layer].values - res['corr_DA_ref_harmonic_sm_' + layer].values
    # ts3 = res['corr_DA_harm_ma_sm_' + layer].values - res['corr_DA_harm_harmonic_sm_' + layer].values
    # ts4 = res['corr_DA_pent_ma_sm_' + layer].values - res['corr_DA_pent_harmonic_sm_' + layer].values
    # title1 = 'OL (moving-average minus harmonic)'
    # title2 = 'DA ref (moving-average minus harmonic)'
    # title3 = 'DA harmonic (moving-average minus harmonic)'
    # title4 = 'DA pentadal (moving-average minus harmonic)'


    vmin = -0.08
    vmax = 0.08

    marker_size = 25
    cmap='seismic_r'

    llcrnrlat = 24
    urcrnrlat = 51
    llcrnrlon = -128
    urcrnrlon = -64

    figsize = (18, 10)

    plt.figure(num=None, figsize=figsize, dpi=90, facecolor='w', edgecolor='k')

    ax = plt.subplot(221)
    m = Basemap(projection='mill', llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat, llcrnrlon=llcrnrlon, urcrnrlon=urcrnrlon,
                resolution='c')
    m.drawcoastlines()
    m.drawcountries()
    m.drawstates()
    x, y = m(lons, lats)
    sc = ax.scatter(x, y, s=marker_size, c=ts1, marker='o', cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_title(title1)
    plt.colorbar(sc,orientation='horizontal',fraction=0.0475,pad=0.02,aspect=40)

    ax = plt.subplot(222)
    m = Basemap(projection='mill', llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat, llcrnrlon=llcrnrlon, urcrnrlon=urcrnrlon,
                resolution='c')
    m.drawcoastlines()
    m.drawcountries()
    m.drawstates()
    x, y = m(lons, lats)
    sc = ax.scatter(x, y, s=marker_size, c=ts2, marker='o', cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_title(title2)
    plt.colorbar(sc,orientation='horizontal',fraction=0.0475,pad=0.02,aspect=40)

    ax = plt.subplot(223)
    m = Basemap(projection='mill', llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat, llcrnrlon=llcrnrlon, urcrnrlon=urcrnrlon,
                resolution='c')
    m.drawcoastlines()
    m.drawcountries()
    m.drawstates()
    x, y = m(lons, lats)
    sc = ax.scatter(x, y, s=marker_size, c=ts3, marker='o', cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_title(title3)
    plt.colorbar(sc,orientation='horizontal',fraction=0.0475,pad=0.02,aspect=40)

    ax = plt.subplot(224)
    m = Basemap(projection='mill', llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat, llcrnrlon=llcrnrlon, urcrnrlon=urcrnrlon,
                resolution='c')
    m.drawcoastlines()
    m.drawcountries()
    m.drawstates()
    x, y = m(lons, lats)
    sc = ax.scatter(x, y, s=marker_size, c=ts4, marker='o', cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_title(title4)
    plt.colorbar(sc,orientation='horizontal',fraction=0.0475,pad=0.02,aspect=40)

    plt.tight_layout()
    plt.show()

def spatial_plot_TCA_res():

    fname = r"D:\work\LDAS\2018-06_rmse_uncertainty\TCA_evaluation\validation.csv"
    res = pd.read_csv(fname)

    lats = res['lat'].values
    lons = res['lon'].values

    layer = 'surface'

    vmin = 0
    vmax = 0.05

    marker_size = 100
    cmap='hot_r'

    llcrnrlat = 24
    urcrnrlat = 51
    llcrnrlon = -128
    urcrnrlon = -64

    figsize = (18, 10)

    plt.figure(num=None, figsize=figsize, dpi=90, facecolor='w', edgecolor='k')

    fontsize = 20

    # ------------------------------------------------------------------------------------------------------------------
    ax = plt.subplot(111)

    # c = res['RMSE_model_DA_const_err_absolute_sm_surface'].values - res['RMSE_model_DA_varia_err_absolute_sm_surface'].values
    c = res['RMSE_model_DA_varia_err_absolute_sm_surface'].values

    m = Basemap(projection='mill', llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat, llcrnrlon=llcrnrlon, urcrnrlon=urcrnrlon,
                resolution='c')
    m.drawcoastlines()
    m.drawcountries()
    m.drawstates()
    x, y = m(lons, lats)
    sc = ax.scatter(x, y, s=marker_size, c=c, marker='o', cmap=cmap, vmin=vmin, vmax=vmax)
    # ax.set_title('Constant error - variable error', fontsize=fontsize)
    ax.set_title('ubRMSE variable error', fontsize=fontsize)
    cb = plt.colorbar(sc,orientation='horizontal',fraction=0.0475,pad=0.02,aspect=40)
    cb.ax.tick_params(labelsize=fontsize)

    plt.tight_layout()
    plt.show()

def spatial_plot_TCA_beta():

    fname = r"D:\work\LDAS\2018-06_rmse_uncertainty\TCA_evaluation\validation.csv"
    res = pd.read_csv(fname)

    lats = res['lat'].values
    lons = res['lon'].values

    layer = 'surface'

    vmin = 0.4
    vmax = 1.6

    marker_size = 100
    cmap='seismic_r'

    llcrnrlat = 24
    urcrnrlat = 51
    llcrnrlon = -128
    urcrnrlon = -64

    figsize = (18, 10)

    plt.figure(num=None, figsize=figsize, dpi=90, facecolor='w', edgecolor='k')

    fontsize = 20

    # ------------------------------------------------------------------------------------------------------------------
    ax = plt.subplot(111)

    c = res['beta_insitu_DA_varia_err_absolute_sm_surface'].values

    m = Basemap(projection='mill', llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat, llcrnrlon=llcrnrlon, urcrnrlon=urcrnrlon,
                resolution='c')
    m.drawcoastlines()
    m.drawcountries()
    m.drawstates()
    x, y = m(lons, lats)
    sc = ax.scatter(x, y, s=marker_size, c=c, marker='o', cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_title('scaling factor', fontsize=fontsize)
    cb = plt.colorbar(sc,orientation='horizontal',fraction=0.0475,pad=0.02,aspect=40)
    cb.ax.tick_params(labelsize=fontsize)

    plt.tight_layout()
    plt.show()

def plot_ismn_statistics():

    res = pd.read_csv(r"D:\work\LDAS\2018-06_rmse_uncertainty\TCA_evaluation\validation.csv")
    res2 = pd.read_csv(r"D:\work\LDAS\2018-06_rmse_uncertainty\insitu_evaluation\validation.csv")

    # variables = ['sm_surface','sm_rootzone','sm_profile']
    modes = ['absolute',]
    runs = ['noDA', 'DA_const_err', 'DA_varia_err']

    networks = ['SCAN','USCRN']
    res.index = res.network
    res = res.loc[networks,:]
    title = ','.join(networks)
    # title = 'all networks'

    # r_title = 'R (rzsm) '+ title
    ubrmsd_title = 'ubRMSD (ssm) '+ title
    ubrmse_title = 'ubRMSE (ssm) '+ title
    var = 'sm_surface'

    plt.figure(figsize=(10,8))

    offsets = [-0.2,0.0,0.2]
    cols = ['lightblue', 'lightgreen', 'coral']
    # offsets = [-0.3,-0.1,0.1,0.3]
    # cols = ['lightblue', 'lightgreen', 'coral', 'brown']
    fontsize=12

    # ax = plt.subplot(211)
    # plt.grid(color='k', linestyle='--', linewidth=0.25)
    #
    # data = list()
    # ticks = list()
    # pos = list()
    # colors = list()
    #
    # for i,mode in enumerate(modes):
    #     ticks.append(mode)
    #     for col,offs,run in zip(cols,offsets,runs):
    #         tmp_data = res['corr_'+run+'_'+mode+'_'+var].values
    #         tmp_data = tmp_data[~np.isnan(tmp_data)]
    #         data.append(tmp_data)
    #         pos.append(i+1 + offs)
    #         colors.append(col)
    # box = ax.boxplot(data, whis=[5,95], showfliers=False, positions=pos, widths=0.1, patch_artist=True)
    # for patch, color in zip(box['boxes'], colors):
    #     patch.set(color='black', linewidth=2)
    #     patch.set_facecolor(color)
    # for patch in box['medians']:
    #     patch.set(color='black', linewidth=2)
    # for patch in box['whiskers']:
    #     patch.set(color='black', linewidth=1)
    # plt.figlegend((box['boxes'][0:4]),runs,'upper left',fontsize=fontsize)
    # plt.xticks(np.arange(len(modes))+1, ticks,fontsize=fontsize)
    # plt.yticks(fontsize=fontsize)
    # plt.xlim(0.5,len(ticks)+0.5)
    # plt.ylim(0.0,1.0)
    # for i in np.arange(len(modes)):
    #     plt.axvline(i+0.5, linewidth=1, color='k')
    # ax.set_title(r_title ,fontsize=fontsize)

    # ---------------------------------------------------------------------------------------------------------
    ax = plt.subplot(121)
    plt.grid(color='k', linestyle='--', linewidth=0.25)

    data = list()
    ticks = list()
    pos = list()
    colors = list()

    for i, mode in enumerate(modes):
        ticks.append(mode)
        for col, offs, run in zip(cols, offsets, runs):
            tmp_data = res2['ubrmsd_' + run + '_' + mode + '_' + var].values
            tmp_data = tmp_data[~np.isnan(tmp_data)]
            data.append(tmp_data)
            pos.append(i + 1 + offs)
            colors.append(col)
    box = ax.boxplot(data, whis=[5, 95], showfliers=False, positions=pos, widths=0.1, patch_artist=True)
    for patch, color in zip(box['boxes'], colors):
        patch.set(color='black', linewidth=2)
        patch.set_facecolor(color)
    for patch in box['medians']:
        patch.set(color='black', linewidth=2)
    for patch in box['whiskers']:
        patch.set(color='black', linewidth=1)
    plt.figlegend((box['boxes'][0:4]), runs, 'upper left', fontsize=fontsize)
    plt.xticks(np.arange(len(modes)) + 1, ticks, fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.xlim(0.5, len(ticks) + 0.5)
    plt.ylim(0, 0.1)
    for i in np.arange(len(modes)):
        plt.axvline(i + 0.5, linewidth=1, color='k')
    ax.set_title(ubrmsd_title, fontsize=fontsize)

    # ---------------------------------------------------------------------------------------------------------
    ax = plt.subplot(122)
    plt.grid(color='k', linestyle='--', linewidth=0.25)

    data = list()
    ticks = list()
    pos = list()
    colors = list()

    for i,mode in enumerate(modes):
        ticks.append(mode)
        for col,offs,run in zip(cols,offsets,runs):
            tmp_data = res['RMSE_model_'+run+'_'+mode+'_'+var].values
            tmp_data = tmp_data[~np.isnan(tmp_data)]
            data.append(tmp_data)
            pos.append(i+1 + offs)
            colors.append(col)
    box = ax.boxplot(data, whis=[5,95], showfliers=False, positions=pos, widths=0.1, patch_artist=True)
    for patch, color in zip(box['boxes'], colors):
        patch.set(color='black', linewidth=2)
        patch.set_facecolor(color)
    for patch in box['medians']:
        patch.set(color='black', linewidth=2)
    for patch in box['whiskers']:
        patch.set(color='black', linewidth=1)
    plt.figlegend((box['boxes'][0:4]),runs,'upper left',fontsize=fontsize)
    plt.xticks(np.arange(len(modes))+1, ticks,fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.xlim(0.5,len(ticks)+0.5)
    plt.ylim(0,0.1)
    for i in np.arange(len(modes)):
        plt.axvline(i+0.5, linewidth=1, color='k')
    ax.set_title(ubrmse_title,fontsize=fontsize)

    plt.show()

def plot_ensemble_uncertainty_vs_ubrmsd():

    DA_const_err = LDAS_io('ensstd', 'US_M36_SMOS40_DA_cal_scaled')
    DA_varia_err = LDAS_io('ensstd', 'US_M36_SMOS40_DA_cal_scl_errfile')

    t_ana = pd.DatetimeIndex(LDAS_io('ObsFcstAna', 'US_M36_SMOS40_DA_cal_scaled').timeseries.time.values).sort_values()

    res = pd.read_csv(r'D:\work\LDAS\2018-06_rmse_uncertainty\insitu_evaluation\validation.csv',index_col=0)
    res2 = pd.read_csv(r'D:\work\LDAS\2018-06_rmse_uncertainty\TCA_evaluation\validation.csv',index_col=0)

    res['RMSE_model_DA_const_err_absolute_sm_surface'] = res2['RMSE_model_DA_const_err_absolute_sm_surface']
    res['RMSE_model_DA_varia_err_absolute_sm_surface'] = res2['RMSE_model_DA_varia_err_absolute_sm_surface']

    res['ensstd_const_err'] = np.nan
    res['ensstd_varia_err'] = np.nan

    param = 'sm_surface'

    for  idx, vals in res.iterrows():
        print(idx)
        res.loc[idx, 'ensstd_const_err'] = DA_const_err.timeseries[param][vals['ease_row'], vals['ease_col'], :].to_pandas().loc[t_ana - pd.to_timedelta('2 hours')].mean()
        res.loc[idx, 'ensstd_varia_err'] = DA_varia_err.timeseries[param][vals['ease_row'], vals['ease_col'], :].to_pandas().loc[t_ana - pd.to_timedelta('2 hours')].mean()

    xlim = [0,0.12]
    ylim = [0,0.12]

    plt.figure(figsize=(13,6))

    # ---------------------------------------------------------------------------------
    ax = plt.subplot(121)

    xx = res['ensstd_const_err']
    yy = res['ubrmsd_DA_const_err_absolute_sm_surface']
    zz = res['RMSE_model_DA_const_err_absolute_sm_surface']

    a = res[['ubrmsd_DA_const_err_absolute_sm_surface','RMSE_model_DA_const_err_absolute_sm_surface']]
    b = res[['ensstd_const_err',]]
    print(a.apply(lambda col: col.corr(b.ix[:,0], method='spearman'), axis=0))

    ax.plot(xx, yy, 'o', markersize=3, markerfacecolor='k', markeredgecolor='k')
    # (xx - yy).hist(bins=20, range=(-0.2, 0.02))
    # (xx - zz).hist(bins=20, range=(-0.06, 0.06))

    ax.plot(xlim,ylim,'--k')
    ax.set_title('Constant observation error')
    ax.set_xlim(xlim)
    ax.set_ylim(xlim)
    # ax.set_xlabel('ensemble standard deviation minus ubRMSD / TCA RMSE')
    ax.set_xlabel('ensemble standard deviation')
    ax.set_ylabel('ubRMSD')

    # print(np.percentile((xx-yy).dropna(), [5,25,50,75,95]))
    # print(np.percentile(yy.dropna(), [5,25,50,75,95]))


    # ---------------------------------------------------------------------------------

    ax = plt.subplot(122)
    xx = res['ensstd_varia_err']
    yy = res['ubrmsd_DA_varia_err_absolute_sm_surface']
    zz = res['RMSE_model_DA_varia_err_absolute_sm_surface']

    a = res[['ubrmsd_DA_varia_err_absolute_sm_surface', 'RMSE_model_DA_varia_err_absolute_sm_surface']]
    b = res[['ensstd_varia_err', ]]
    print(a.apply(lambda col: col.corr(b.ix[:, 0], method='spearman'), axis=0))

    ax.plot(xx, yy, 'o', markersize=3, markerfacecolor='k', markeredgecolor='k')
    # (xx - yy).hist(bins=20, range=(-0.2, 0.02))
    # (xx - zz).hist(bins=20, range=(-0.06, 0.06))

    ax.plot(xlim, ylim, '--k')
    ax.set_title('Variable observation error')
    ax.set_xlim(xlim)
    ax.set_ylim(xlim)
    # ax.set_xlabel('ensemble standard deviation minus ubRMSD / TCA RMSE')
    ax.set_xlabel('ensemble standard deviation')
    ax.set_ylabel('ubRMSD')

    # print(np.percentile((xx-yy).dropna(), [5,25,50,75,95]))
    # print(np.percentile(yy.dropna(), [5,25,50,75,95]))

    plt.show()

def plot_Tb_uncertainty_vs_ubrmsd():

    res = pd.read_csv(r"D:\work\LDAS\2018-06_rmse_uncertainty\Tb_evaluation\validation.csv", index_col=0)

    xlim = [0,16]
    ylim = [0,4]

    plt.figure(figsize=(13, 6))

    # ---------------------------------------------------------------------------------
    ax = plt.subplot(121)

    xx = res['ubrmsd_const_err']
    yy = res['ensstd_const_err']

    ax.plot(xx, yy, 'o', markersize=3, markerfacecolor='k', markeredgecolor='k')
    # (xx - yy).hist(bins=20, range=(-0.2, 0.02))
    # (xx - zz).hist(bins=20, range=(-0.06, 0.06))

    # ax.plot(xlim,ylim,'--k')
    ax.set_title('Constant observation error')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xlabel('Tb ubRMSD')
    ax.set_ylabel('Tb ensemble standard deviation')
    # ax.set_ylabel('ubRMSD')

    # print(np.percentile((xx - yy).dropna(), [5, 25, 50, 75, 95]))
    # print(np.percentile(yy.dropna(), [5,25,50,75,95]))

    # ---------------------------------------------------------------------------------

    ax = plt.subplot(122)
    xx = res['ubrmsd_varia_err']
    yy = res['ensstd_varia_err']
    ax.plot(xx, yy, 'o', markersize=3, markerfacecolor='k', markeredgecolor='k')
    # (xx - yy).hist(bins=20, range=(-0.2, 0.02))
    # (xx - zz).hist(bins=20, range=(-0.06, 0.06))

    # ax.plot(xlim, ylim, '--k')
    ax.set_title('Variable observation error')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xlabel('Tb ubRMSD')
    ax.set_ylabel('Tb ensemble standard deviation')
    # ax.set_ylabel('ubRMSD')

    print(np.percentile((xx - yy).dropna(), [5, 25, 50, 75, 95]))
    # print(np.percentile(yy.dropna(), [5,25,50,75,95]))

    plt.show()


if __name__=='__main__':
    # plot_filter_diagnostics()
    # plot_ismn_statistics()
    plot_ensemble_uncertainty_vs_ubrmsd()
    # plot_Tb_uncertainty_vs_ubrmsd()
    # spatial_plot_TCA_beta()
    # spatial_plot_TCA_res()
