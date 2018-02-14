
import numpy as np

from ease_grid.ease2_grid import EASE2_grid

class EASE2(EASE2_grid):
    """
    Class that contains EASE2 grid parameters and LDAS specific
    attributes/functions

    Parameters
    ----------
    res : int
        Gridcell resolution [km]
    map_scale : float
        The exact map scale [km]
    tilecoord : pd.DataFrame
        Metadata information for the tile coordinates of a particular LDAS experiment
    tilegrids : pd.DataFrame
        Metadata information for the tile grids of a particular LDAS experiment

    Attributes
    ----------
    tilecoord : pd.DataFrame
        Metadata information for the tile coordinates of a particular LDAS experiment
    tilegrids : pd.DataFrame
        Metadata information for the tile grids of a particular LDAS experiment
    latdim : np.array
        Array containing the latitudes of grid rows
    londim : np.array
        Array containing the longitudes of grid columns

    """

    def __init__(self,
                 res=36,
                 map_scale=None,
                 tilecoord=None,
                 tilegrids=None):

        res *= 1000

        # Default map scales for 3, 9 and 36 km resolution grids are
        # defined to match the exact SMAP grid resolution
        # (When specifying resolution only, 3 and 9 km would not be exact
        #  sub-fractions of the 36 km grid cells)
        if map_scale is None:
            if res == 36000:
                map_scale = 36032.220840584
            elif res == 9000:
                map_scale = 9008.055210146
            elif res == 3000:
                map_scale = 3002.6850700487
        else:
            map_scale *= 1000

        self.tilecoord = tilecoord
        self.tilegrids = tilegrids

        super(EASE2, self).__init__(res, map_scale=map_scale)

    def colrow2tileid(self, col, row):

        ind_col = self.tilecoord['i_indg']-self.tilegrids.loc['domain','i_offg'] == col
        ind_row = self.tilecoord['j_indg']-self.tilegrids.loc['domain','j_offg'] == row

        return self.tilecoord.loc[ind_col & ind_row,'tile_id'].values[0]

    def tileid2colrow(self, tile_id):

        col = (self.tilecoord.loc[self.tilecoord['tile_id']==tile_id, 'i_indg'] - self.tilegrids.loc['domain','i_offg']).values[0]
        row = (self.tilecoord.loc[self.tilecoord['tile_id']==tile_id, 'j_indg'] - self.tilegrids.loc['domain','j_offg']).values[0]

        return col, row

    def colrow2lonlat(self, col, row, glob=True):
        """ Convert col/row (domain-based) into lon/lat """
        if glob is True:
            return self.londim[col], self.latdim[row]
        else:
            return self.londim[col+self.tilegrids.loc['domain','i_offg']], self.latdim[row+self.tilegrids.loc['domain','j_offg']]

    def lonlat2colrow(self, lon, lat, domain=True):
        """ Find nearest GLOBAL tile (col/row) from any given lon/lat """
        londif = np.abs(self.londim - lon)
        latdif = np.abs(self.latdim - lat)
        col = np.where(np.abs(londif-londif.min())<0.0001)[0][0]
        row = np.where(np.abs(latdif-latdif.min())<0.0001)[0][0]

        if domain is True:
            col -= self.tilegrids.loc['domain','i_offg']
            row -= self.tilegrids.loc['domain','j_offg']

        return col, row

    def lonlat2tilenum(self, lon, lat):
        """" Get the (domain-based) tile number for a given lon/lat """
        col, row = self.lonlat2colrow(lon, lat)
        tilenum = np.where((self.tilecoord['i_indg'] == col)&
                           (self.tilecoord['j_indg'] == row))[0][0]
        return tilenum + 1

