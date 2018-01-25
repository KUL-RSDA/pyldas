
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

    ts = LDAS_io('ObsFcstAna').timeseries

    df = pd.DataFrame(index=ts.time)

    df['obs'] = ts['obs_obs'][spc,row,col].values
    df['fcst'] = ts['obs_fcst'][spc,row,col].values

    df.dropna().plot()
    plt.show()

    ts.close()


if __name__=='__main__':
    plot_innov()


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
# llcrnrlat=-58.,
# urcrnrlat=60.,
# llcrnrlon=-132.,
# urcrnrlon=157.,

