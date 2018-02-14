
import os

import numpy as np
import pandas as pd
import xarray as xr

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

from pyldas.interface import LDAS_io

def plot_P_R_estimates():

    exp = 'US_M36_SMOS_DA_calibrated_scaled'

    io = LDAS_io('ObsFcstAna',exp)

    outpath = r"C:\Users\u0116961\Documents\work\LDASsa\2018-02_scaling\diagnostics\uncertainty_estimates"

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

        P_est = ((io.timeseries['obs_obs'][spc, :, :, :] - io.timeseries['obs_fcst'][spc, :, :, :]) *
                 (io.timeseries['obs_ana'][spc, :, :, :] - io.timeseries['obs_fcst'][spc, :, :, :])).mean(dim='time')

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

    ds = xr.open_dataset(r"C:\Users\u0116961\Documents\work\LDASsa\2018-02_scaling\diagnostics\filter_diagnostics.nc")

    # outpath = r"C:\Users\u0116961\Documents\work\LDASsa\2018-02_scaling\diagnostics"

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


def plot_filter_increments(mode='image', var='srfexc'):

    ds = xr.open_dataset(r"C:\Users\u0116961\Documents\work\LDASsa\2018-02_scaling\diagnostics\filter_diagnostics.nc")

    lons = ds.lon.values
    lats = ds.lat.values

    lons,lats = np.meshgrid(lons,lats)

    imgs = [ds['incr_'+var+'_mean'][:, :, 1].values,
            ds['incr_'+var+'_mean'][:, :, 2].values,
            ds['incr_'+var+'_mean'][:, :, 3].values,
            ds['incr_'+var+'_var'][:, :, 1].values,
            ds['incr_'+var+'_var'][:, :, 2].values,
            ds['incr_'+var+'_var'][:, :, 3].values]

    names = [var + ' mean (ref)',
             var + ' mean (pentadal)',
             var + ' mean (harmonic)',
             var + ' var (ref)',
             var + ' var (pentadal)',
             var + ' var (harmonic)']

    f = plt.figure(figsize=(19,9))

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

        ax = plt.subplot(2,3,i)

        if mode == 'image':

            cbrange = (-0.8, 0.8) if name.find('mean') != -1 else (0,1)

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

    ds = xr.open_dataset(r"C:\Users\u0116961\Documents\work\LDASsa\2018-02_scaling\diagnostics\filter_diagnostics.nc")
    outpath = r"C:\Users\u0116961\Documents\work\LDASsa\2018-02_scaling\diagnostics"

    lons = ds.lon.values
    lats = ds.lat.values

    lons,lats = np.meshgrid(lons,lats)


    if spc_mean is True:

        if mode == 'hist':

            imgs = [ds['innov_mean'][:, :, 1, :].values,
                    ds['innov_mean'][:, :, 2, :].values,
                    ds['innov_mean'][:, :, 3, :].values,
                    ds['innov_var'][:, :, 1, :].values,
                    ds['innov_var'][:, :, 2, :].values,
                    ds['innov_var'][:, :, 3, :].values]

        else:

            imgs = [ds['innov_mean'][:, :, 1, :].mean(dim='species').values,
                    ds['innov_mean'][:, :, 2, :].mean(dim='species').values,
                    ds['innov_mean'][:, :, 3, :].mean(dim='species').values,
                    ds['innov_var'][:, :, 1, :].mean(dim='species').values,
                    ds['innov_var'][:, :, 2, :].mean(dim='species').values,
                    ds['innov_var'][:, :, 3, :].mean(dim='species').values]

    names = ['innov mean (ref)',
             'innov mean (pentadal)',
             'innov mean (harmonic)',
             'innov var (ref)',
             'innov var (pentadal)',
             'innov var (harmonic)']

    llcrnrlat = 24
    urcrnrlat = 51
    llcrnrlon = -128
    urcrnrlon = -64

    cmap = 'jet'

    fontsize = 16

    n_spc = 1 if spc_mean is True else len(ds.species)

    for j in np.arange(n_spc):

        if spc_mean is False:
            imgs = [ds['norm_innov_mean'][:, :, 1, j].values,
                    ds['norm_innov_mean'][:, :, 2, j].values,
                    ds['norm_innov_mean'][:, :, 3, j].values,
                    ds['norm_innov_var'][:, :, 1, j].values,
                    ds['norm_innov_var'][:, :, 2, j].values,
                    ds['norm_innov_var'][:, :, 3, j].values]

        f = plt.figure(figsize=(19,9))

        i = 0
        for img,name in zip(imgs,names):
            i += 1

            perc = np.nanpercentile(img, [1,99])

            np.place(img,(img<perc[0])|(img>perc[1]),np.nan)
            img = img.reshape(lons.shape)

            ax = plt.subplot(2,3,i)

            if mode == 'image':

                cbrange = (-1.5, 1.5) if name.find('mean') != -1 else (0,100)

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
                    cbrange = (-1.8, 1.8)
                    ylim=(0,1)
                else :
                    cbrange = (10,60)
                    ylim=(0,0.05)

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

    res = pd.DataFrame.from_csv(fname)

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

def plot_ismn_statistics():
    fname = r"C:\Users\u0116961\Documents\work\ldas_ismn_eval\2018-02_scaling\validation_masked.csv"

    res = pd.DataFrame.from_csv(fname)

    # variables = ['sm_surface','sm_rootzone','sm_profile']
    modes = ['mean','ma','harmonic']
    runs = ['OL','DA_ref','DA_pent','DA_harm']

    # networks = np.unique(res.network)
    networks = ['SCAN','USCRN']
    res.index = res.network
    res = res.loc[networks,:]
    title = ','.join(networks)

    # title = 'all networks'
    r_title = 'R (rsm) '+ title
    ubrmsd_title = 'ubRMSD (rsm) '+ title
    var = 'sm_rootzone'

    plt.figure(figsize=(18,9))

    offsets = [-0.3,-0.1,0.1,0.3]
    cols = ['lightblue', 'lightgreen', 'coral', 'brown']
    fontsize=12

    ax = plt.subplot(211)
    plt.grid(color='k', linestyle='--', linewidth=0.25)

    data = list()
    ticks = list()
    pos = list()
    colors = list()

    for i,mode in enumerate(modes):
        ticks.append(mode)
        for col,offs,run in zip(cols,offsets,runs):
            tmp_data = res['corr_'+run+'_'+mode+'_'+var].values
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
    plt.ylim(0.0,1.0)
    for i in np.arange(len(modes)):
        plt.axvline(i+0.5, linewidth=1, color='k')
    ax.set_title(r_title ,fontsize=fontsize)

    # ---------------------------------------------------------------------------------------------------------
    ax = plt.subplot(212)
    plt.grid(color='k', linestyle='--', linewidth=0.25)

    data = list()
    ticks = list()
    pos = list()
    colors = list()

    for i,mode in enumerate(modes):
        ticks.append(mode)
        for col,offs,run in zip(cols,offsets,runs):
            tmp_data = res['ubrmsd_'+run+'_'+mode+'_'+var].values
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
    # plt.figlegend((box['boxes'][0:4]),runs,'upper left',fontsize=fontsize)
    plt.xticks(np.arange(len(modes))+1, ticks,fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.xlim(0.5,len(ticks)+0.5)
    plt.ylim(0,0.1)
    for i in np.arange(len(modes)):
        plt.axvline(i+0.5, linewidth=1, color='k')
    ax.set_title(ubrmsd_title,fontsize=fontsize)


    plt.show()

if __name__=='__main__':
    plot_P_R_estimates()