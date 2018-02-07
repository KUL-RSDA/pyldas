
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

from pyldas.grids import EASE2

from pyldas.readers import LDAS_io

def plot_ease_img(data,tag,
                  llcrnrlat=24,
                  urcrnrlat=51,
                  llcrnrlon=-128,
                  urcrnrlon=-64,
                  figsize=(20,10),
                  cbrange=(-20,20),
                  cmap='jet',
                  fontsize=20):

    grid = EASE2()

    lons,lats = np.meshgrid(grid.londim,grid.latdim)

    img = np.empty(lons.shape, dtype='float32')
    img.fill(None)

    ind_lat = grid.tilecoord.loc[data.index.values,'j_indg']
    ind_lon = grid.tilecoord.loc[data.index.values,'i_indg']

    img[ind_lat,ind_lon] = data[tag]
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

    tag = 'precipitation_total_surface_flux'
    tag = 'snow_mass'
    cbrange = (0,0.0001)
    cbrange = (0,0.6)
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

    print len(ts_scl['obs_obs'][spc,row,col].dropna('time'))
    print len(ts_scl['obs_fcst'][spc,row,col].dropna('time'))
    print len(ts_usc['obs_obs'][spc,row,col].dropna('time'))
    print len(ts_usc['obs_fcst'][spc,row,col].dropna('time'))

    plt.tight_layout()
    plt.show()

    ts_scl.close()
    ts_usc.close()

def plot_scaling_parameters():

    # fname = r"C:\Users\u0116961\Documents\VSC\vsc_data_copies\scratch_TEST_RUNS\obs_scaling_old\7Thv_TbSM_001_SMOS_zscore_stats_2010_p37_2015_p36_hscale_0.00_W_9p_Nmin_20_A_p38.bin"
    # fname = r"C:\Users\u0116961\Documents\VSC\vsc_data_copies\scratch_TEST_RUNS\US_M36_SMOS_noDA_unscaled\obs_scaling\pentadal_mean\7Thv_TbSM_001_SMOS_zscore_stats_2010_p37_2015_p36_hscale_0.00_W_9p_Nmin_20_A_p38.bin"
    fname = r"C:\Users\u0116961\Documents\VSC\vsc_data_copies\scratch_TEST_RUNS\US_M36_SMOS_noDA_unscaled\obs_scaling\harmonic_mean\7Thv_TbSM_001_SMOS_zscore_stats_2010_p37_2015_p36_hscale_0.00_W_9p_Nmin_20_A_p38.bin"

    io = LDAS_io('scale')

    res = io.read_scaling_parameters(fname=fname)

    angle = 30

    res = res[['lon','lat','m_mod_H_%2i'%angle,'m_mod_V_%2i'%angle,'m_obs_H_%2i'%angle,'m_obs_V_%2i'%angle]]
    res.replace(-9999.,np.nan,inplace=True)

    lats = res['lat'].values
    lons = res['lon'].values

    llcrnrlat = 24
    urcrnrlat = 51
    llcrnrlon = -128
    urcrnrlon = -64

    figsize = (17,8)

    plt.figure(num=None, figsize=figsize, dpi=90, facecolor='w', edgecolor='k')

    ax = plt.subplot(221)
    m = Basemap(projection='mill',llcrnrlat=llcrnrlat,urcrnrlat=urcrnrlat,llcrnrlon=llcrnrlon,urcrnrlon=urcrnrlon,resolution='c')
    m.drawcoastlines()
    m.drawcountries()
    m.drawstates()
    x,y = m(lons,lats)
    ax.scatter(x,y,s=10,c=res['m_obs_H_%2i'%angle].values,marker='o', cmap='jet',vmin=220,vmax=300)

    ax = plt.subplot(222)
    m = Basemap(projection='mill',llcrnrlat=llcrnrlat,urcrnrlat=urcrnrlat,llcrnrlon=llcrnrlon,urcrnrlon=urcrnrlon,resolution='c')
    m.drawcoastlines()
    m.drawcountries()
    m.drawstates()
    x,y = m(lons,lats)
    ax.scatter(x,y,s=10,c=res['m_mod_H_%2i'%angle].values,marker='o', cmap='jet', vmin=220, vmax=300)


    ax = plt.subplot(223)
    m = Basemap(projection='mill', llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat, llcrnrlon=llcrnrlon, urcrnrlon=urcrnrlon,
                resolution='c')
    m.drawcoastlines()
    m.drawcountries()
    m.drawstates()
    x, y = m(lons, lats)
    ax.scatter(x, y, s=10, c=res['m_obs_V_%2i'%angle].values, marker='o', cmap='jet', vmin=220, vmax=300)

    ax = plt.subplot(224)
    m = Basemap(projection='mill', llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat, llcrnrlon=llcrnrlon, urcrnrlon=urcrnrlon,
                resolution='c')
    m.drawcoastlines()
    m.drawcountries()
    m.drawstates()
    x, y = m(lons, lats)
    ax.scatter(x, y, s=10, c=res['m_mod_V_%2i'%angle].values, marker='o', cmap='jet', vmin=220, vmax=300)

    plt.tight_layout()
    plt.show()

if __name__=='__main__':
    plot_scaling_parameters()


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

