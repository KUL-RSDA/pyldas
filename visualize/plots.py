
import numpy as np
import pandas as pd

import platform
if platform.system() in ['Linux']:
    import matplotlib
    matplotlib.use("TkAgg")

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

from pyldas.interface import LDAS_io
from pyldas.paths import paths


def plot_catparams():

    outpath = paths().plots / 'catparams'
    if not outpath.exists():
        outpath.mkdir(parents=True)

    exp = 'US_M36_SMOS40_DA_cal_scaled'

    tc = LDAS_io().grid.tilecoord
    tg = LDAS_io().grid.tilegrids

    tc.i_indg -= tg.loc['domain','i_offg'] # col / lon
    tc.j_indg -= tg.loc['domain','j_offg'] # row / lat

    lons = np.unique(tc.com_lon.values)
    lats = np.unique(tc.com_lat.values)[::-1]

    lons, lats = np.meshgrid(lons, lats)

    llcrnrlat = 24
    urcrnrlat = 51
    llcrnrlon = -128
    urcrnrlon = -64
    figsize = (20, 10)
    # cbrange = (-20, 20)
    cmap = 'jet'
    fontsize = 20

    params = LDAS_io(exp=exp).read_params('catparam')

    for param in params:

        fname = outpath / (param + '.png')

        img = np.full(lons.shape, np.nan)
        img[tc.j_indg.values, tc.i_indg.values] = params[param].values
        img_masked = np.ma.masked_invalid(img)

        f = plt.figure(num=None, figsize=figsize, dpi=90, facecolor='w', edgecolor='k')

        m = Basemap(projection='mill',
                    llcrnrlat=llcrnrlat,
                    urcrnrlat=urcrnrlat,
                    llcrnrlon=llcrnrlon,
                    urcrnrlon=urcrnrlon,
                    resolution='c')

        m.drawcoastlines()
        m.drawcountries()
        m.drawstates()

        im = m.pcolormesh(lons, lats, img_masked, cmap=cmap, latlon=True)

        cb = m.colorbar(im, "bottom", size="7%", pad="8%")

        for t in cb.ax.get_xticklabels():
            t.set_fontsize(fontsize)
        for t in cb.ax.get_yticklabels():
            t.set_fontsize(fontsize)

        plt.title(param)

        plt.savefig(fname, dpi=f.dpi)
        plt.close()


def plot_rtm_parameters():

    root = paths().plots / 'RTM_parameters'

    experiments = ['US_M36_SMOS_DA_calibrated_scaled', 'US_M36_SMOS_DA_nocal_scaled_harmonic']

    tc = LDAS_io().grid.tilecoord
    tg = LDAS_io().grid.tilegrids

    tc.i_indg -= tg.loc['domain','i_offg'] # col / lon
    tc.j_indg -= tg.loc['domain','j_offg'] # row / lat

    lons = np.unique(tc.com_lon.values)
    lats = np.unique(tc.com_lat.values)[::-1]

    lons, lats = np.meshgrid(lons, lats)

    llcrnrlat = 24
    urcrnrlat = 51
    llcrnrlon = -128
    urcrnrlon = -64
    figsize = (20, 10)
    # cbrange = (-20, 20)
    cmap = 'jet'
    fontsize = 20

    for exp in experiments:

        outpath = root / exp
        if not outpath.exists:
            outpath.mkdir(exists=True)

        params = LDAS_io(exp=exp).read_params('RTMparam')

        for param in params:

            fname = outpath / (param + '.png')

            img = np.full(lons.shape, np.nan)
            img[tc.j_indg.values, tc.i_indg.values] = params[param].values
            img_masked = np.ma.masked_invalid(img)

            f = plt.figure(num=None, figsize=figsize, dpi=90, facecolor='w', edgecolor='k')

            m = Basemap(projection='mill',
                        llcrnrlat=llcrnrlat,
                        urcrnrlat=urcrnrlat,
                        llcrnrlon=llcrnrlon,
                        urcrnrlon=urcrnrlon,
                        resolution='c')

            m.drawcoastlines()
            m.drawcountries()
            m.drawstates()

            im = m.pcolormesh(lons, lats, img_masked, cmap=cmap, latlon=True)

            # im.set_clim(vmin=cbrange[0], vmax=cbrange[1])

            cb = m.colorbar(im, "bottom", size="7%", pad="8%")

            for t in cb.ax.get_xticklabels():
                t.set_fontsize(fontsize)
            for t in cb.ax.get_yticklabels():
                t.set_fontsize(fontsize)

            plt.title(param)

            plt.savefig(fname, dpi=f.dpi)
            plt.close()


def plot_rtm_parameter_differences():

    outpath = paths().plots / 'RTM_parameters' / 'differences'
    if not outpath.exists:
        outpath.mkdir(exists=True)

    tc = LDAS_io().grid.tilecoord
    tg = LDAS_io().grid.tilegrids

    tc.i_indg -= tg.loc['domain','i_offg'] # col / lon
    tc.j_indg -= tg.loc['domain','j_offg'] # row / lat

    lons = np.unique(tc.com_lon.values)
    lats = np.unique(tc.com_lat.values)[::-1]

    lons, lats = np.meshgrid(lons, lats)

    llcrnrlat = 24
    urcrnrlat = 51
    llcrnrlon = -128
    urcrnrlon = -64
    figsize = (20, 10)
    #
    cmap = 'RdYlBu'
    fontsize = 20

    params_cal = LDAS_io(exp='US_M36_SMOS_DA_calibrated_scaled').read_params('RTMparam')
    params_uncal = LDAS_io(exp='US_M36_SMOS_DA_nocal_scaled_harmonic').read_params('RTMparam')

    for param in params_cal:

        if (param == 'bh')|(param == 'bv'):
            cbrange = (-0.3, 0.3)
        elif (param == 'omega'):
            cbrange = (-0.1, 0.1)
        else:
            cbrange = (-1, 1)

        fname = outpath / (param + '.png')

        img = np.full(lons.shape, np.nan)
        img[tc.j_indg.values, tc.i_indg.values] = params_cal[param].values - params_uncal[param].values
        img_masked = np.ma.masked_invalid(img)

        f = plt.figure(num=None, figsize=figsize, dpi=90, facecolor='w', edgecolor='k')

        m = Basemap(projection='mill',
                    llcrnrlat=llcrnrlat,
                    urcrnrlat=urcrnrlat,
                    llcrnrlon=llcrnrlon,
                    urcrnrlon=urcrnrlon,
                    resolution='c')

        m.drawcoastlines()
        m.drawcountries()
        m.drawstates()

        im = m.pcolormesh(lons, lats, img_masked, cmap=cmap, latlon=True)

        im.set_clim(vmin=cbrange[0], vmax=cbrange[1])

        cb = m.colorbar(im, "bottom", size="7%", pad="8%")

        for t in cb.ax.get_xticklabels():
            t.set_fontsize(fontsize)
        for t in cb.ax.get_yticklabels():
            t.set_fontsize(fontsize)

        plt.title(param)

        plt.savefig(fname, dpi=f.dpi)
        plt.close()


def plot_ease_img(data,tag,
                  llcrnrlat=24,
                  urcrnrlat=51,
                  llcrnrlon=-128,
                  urcrnrlon=-64,
                  figsize=(20,10),
                  cbrange=(-20,20),
                  cmap='jet',
                  title='',
                  fontsize=20):

    io = LDAS_io()

    tc = io.grid.tilecoord

    lons,lats = np.meshgrid(io.grid.ease_lons, io.grid.ease_lats)

    img = np.empty(lons.shape, dtype='float32')
    img.fill(None)

    ind_lat = tc.loc[data.index.values,'j_indg']
    ind_lon = tc.loc[data.index.values,'i_indg']

    img[ind_lat,ind_lon] = data[tag]
    img_masked = np.ma.masked_invalid(img)

    plt.figure(num=None, figsize=figsize, dpi=90, facecolor='w', edgecolor='k')

    m = Basemap(projection='mill',
                llcrnrlat=llcrnrlat,
                urcrnrlat=urcrnrlat,
                llcrnrlon=llcrnrlon,
                urcrnrlon=urcrnrlon,
                resolution='c')

    m.drawcoastlines()
    m.drawcountries()
    m.drawstates()

    im = m.pcolormesh(lons, lats, img_masked, cmap=cmap, latlon=True)

    im.set_clim(vmin=cbrange[0], vmax=cbrange[1])

    cb = m.colorbar(im, "bottom", size="7%", pad="8%")

    for t in cb.ax.get_xticklabels():
        t.set_fontsize(fontsize)
    for t in cb.ax.get_yticklabels():
        t.set_fontsize(fontsize)

    plt.title(title,fontsize=fontsize)

    plt.tight_layout()
    plt.show()


def plot_grid_coord_indices():

    io = LDAS_io('ObsFcstAna', exp='US_M36_SMOS40_noDA_cal_scaled')

    lats = io.images.lat.values
    lons = io.images.lon.values

    plt.figure(figsize=(10, 5))

    llcrnrlat = 24
    urcrnrlat = 51
    llcrnrlon = -128
    urcrnrlon = -64
    m = Basemap(projection='mill', llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat, llcrnrlon=llcrnrlon,
                urcrnrlon=urcrnrlon,
                resolution='c')
    m.drawcoastlines(linewidth=0.5)
    m.drawcountries(linewidth=0.5)
    m.drawstates(linewidth=0.1)

    lats = lats[np.arange(0, len(lats), 15)]
    lons = lons[np.arange(0, len(lons), 15)]
    m.drawparallels(lats, labels=[False, False, False, False], linestyle='--', linewidth=1, color='red')
    m.drawmeridians(lons, labels=[False, False, False, False], linestyle='--', linewidth=1, color='red')

    x = np.zeros(len(lons))
    for i, lon in enumerate(lons):
        x[i], tmp = m(lon, lats[0])

    y = np.zeros(len(lats))
    for i, lat in enumerate(lats):
        tmp, y[i] = m(lons[-1], lat)

    plt.xticks(x[0:-1], np.arange(0, len(lons) - 1) * 15)
    plt.yticks(y[1::], np.arange(1, len(lats)) * 15)

    plt.show()


def plot_ObsFcstAna_image(species=8):

    io = LDAS_io('ObsFcstAna')
    img = io.read_image(2011, 7, 10, 0, 0)

    img = img[img['obs_species']==species]

    tag = 'innov'
    img[tag] = img['obs_obs']-img['obs_fcst']
    img.index = img['obs_tilenum'].values
    plot_ease_img(img, tag)


def plot_model_image():

    io = LDAS_io('xhourly')
    img = io.read_image(2011, 4, 20, 10, 30)

    # tag = 'precipitation_total_surface_flux'
    tag = 'snow_mass'
    # cbrange = (0,0.0001)
    # cbrange = (0,0.6)
    cbrange = (0,100)

    plot_ease_img(img, tag, cbrange=cbrange)

def plot_innov(spc=8, row=35, col=65):

    ts_scl = LDAS_io('ObsFcstAna',exp='US_M36_SMOS_noDA_scaled').timeseries
    ts_usc = LDAS_io('ObsFcstAna',exp='US_M36_SMOS_noDA_unscaled').timeseries

    plt.figure(figsize=(18,11))

    ax1 = plt.subplot(311)
    df = pd.DataFrame(index=ts_scl.time)
    df['obs'] = ts_scl['obs_obs'][spc,row,col].values
    df['fcst'] = ts_scl['obs_fcst'][spc,row,col].values
    df.dropna().plot(ax=ax1)

    ax2 = plt.subplot(312)
    df = pd.DataFrame(index=ts_usc.time)
    df['obs'] = ts_usc['obs_obs'][spc,row,col].values
    df['fcst'] = ts_usc['obs_fcst'][spc,row,col].values
    df.dropna().plot(ax=ax2)

    ax3 = plt.subplot(313)
    df = pd.DataFrame(index=ts_usc.time)
    df['obs_diff'] = ts_scl['obs_obs'][spc,row,col].values - ts_usc['obs_obs'][spc,row,col].values
    df['fcst_diff'] = ts_scl['obs_fcst'][spc,row,col].values - ts_usc['obs_fcst'][spc,row,col].values
    df.dropna().plot(ax=ax3)

    print(len(ts_scl['obs_obs'][spc,row,col].dropna('time')))
    print(len(ts_scl['obs_fcst'][spc,row,col].dropna('time')))
    print(len(ts_usc['obs_obs'][spc,row,col].dropna('time')))
    print(len(ts_usc['obs_fcst'][spc,row,col].dropna('time')))

    plt.tight_layout()
    plt.show()

    ts_scl.close()
    ts_usc.close()

def plot_xarr_img(img,lons,lats,
                  llcrnrlat=24,
                  urcrnrlat=51,
                  llcrnrlon=-128,
                  urcrnrlon=-64,
                  cbrange=(-20,20),
                  cmap='jet',
                  fontsize=16):

    lons, lats = np.meshgrid(lons, lats)
    img_masked = np.ma.masked_invalid(img)
    m = Basemap(projection='mill',
                llcrnrlat=llcrnrlat,
                urcrnrlat=urcrnrlat,
                llcrnrlon=llcrnrlon,
                urcrnrlon=urcrnrlon,
                resolution='c')
    m.drawcoastlines()
    m.drawcountries()
    m.drawstates()
    im = m.pcolormesh(lons, lats, img_masked, cmap=cmap, latlon=True)

    im.set_clim(vmin=cbrange[0], vmax=cbrange[1])
    cb = m.colorbar(im, "bottom", size="7%", pad="8%")
    for t in cb.ax.get_xticklabels():
        t.set_fontsize(fontsize)
    for t in cb.ax.get_yticklabels():
        t.set_fontsize(fontsize)

def plot_obs_uncertainties():

    io = LDAS_io('ObsFcstAna',exp='US_M36_SMOS40_DA_cal_scl_errfile')

    lons = io.images['lon'].values
    lats = io.images['lat'].values

    figsize=(18,9)
    f = plt.figure(num=None, figsize=figsize, dpi=90, facecolor='w', edgecolor='k')

    ax = f.add_subplot(2,2,1)
    obserr = io.images.sel(species=1)['obs_obsvar'].mean('time').values
    plot_xarr_img(obserr,lons,lats,cbrange=[0,250])
    ax.set_title('Asc / H-pol', fontsize=16)

    ax = f.add_subplot(2,2,2)
    obserr = io.images.sel(species=2)['obs_obsvar'].mean('time').values
    plot_xarr_img(obserr,lons,lats,cbrange=[0,250])
    ax.set_title('Dsc / H-pol', fontsize=16)

    ax = f.add_subplot(2,2,3)
    obserr = io.images.sel(species=3)['obs_obsvar'].mean('time').values
    plot_xarr_img(obserr,lons,lats,cbrange=[0,250])
    ax.set_title('Asc / V-pol', fontsize=16)

    ax = f.add_subplot(2,2,4)
    obserr = io.images.sel(species=4)['obs_obsvar'].mean('time').values
    plot_xarr_img(obserr,lons,lats,cbrange=[0,250])
    ax.set_title('Dsc / V-pol', fontsize=16)

    plt.tight_layout()
    plt.show()

def plot_fcst_uncertainties():

    io = LDAS_io('ObsFcstAna',exp='US_M36_SMOS40_noDA_cal_scaled')

    lons = io.images['lon'].values
    lats = io.images['lat'].values

    figsize=(18,9)
    f = plt.figure(num=None, figsize=figsize, dpi=90, facecolor='w', edgecolor='k')

    ax = f.add_subplot(2,2,1)
    obserr = io.images.sel(species=1)['obs_fcstvar'].mean('time').values
    plot_xarr_img(obserr,lons,lats,cbrange=[0,50])
    ax.set_title('Asc / H-pol', fontsize=16)

    ax = f.add_subplot(2,2,2)
    obserr = io.images.sel(species=2)['obs_fcstvar'].mean('time').values
    plot_xarr_img(obserr,lons,lats,cbrange=[0,50])
    ax.set_title('Dsc / H-pol', fontsize=16)

    ax = f.add_subplot(2,2,3)
    obserr = io.images.sel(species=3)['obs_fcstvar'].mean('time').values
    plot_xarr_img(obserr,lons,lats,cbrange=[0,50])
    ax.set_title('Asc / V-pol', fontsize=16)

    ax = f.add_subplot(2,2,4)
    obserr = io.images.sel(species=4)['obs_fcstvar'].mean('time').values
    plot_xarr_img(obserr,lons,lats,cbrange=[0,50])
    ax.set_title('Dsc / V-pol', fontsize=16)

    plt.tight_layout()
    plt.show()



if __name__=='__main__':
    plot_catparams()


# llcrnrlat = -58.,
# urcrnrlat = 78.,
# llcrnrlon = -172.,
# urcrnrlon = 180.,

# clipped global
# llcrnrlat=-58.,
# urcrnrlat=60.,
# llcrnrlon=-132.,
# urcrnrlon=157.,

# USA
# llcrnrlat = 24.,
# urcrnrlat = 51.,
# llcrnrlon = -128.,
# urcrnrlon = -64.,

