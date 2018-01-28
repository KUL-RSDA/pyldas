
import numpy as np

from ease_grid.ease2_grid import EASE2_grid

from pyldas.functions import find_files
from pyldas.constants import paths
from pyldas.readers import LDAS_io

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
    tileinfo_path : string
        Path for non-default tile grid definitions
        (i.e., for specific experiments/domains)

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

    def __init__(self, res=36, map_scale=None, exp=None, domain=None):

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

        io = LDAS_io(exp=exp, domain=domain)
        self.tilecoord = io.tilecoord
        self.tilegrids = io.tilegrids

        super(EASE2, self).__init__(res, map_scale=map_scale)


    def colrow2lonlat(self, col, row):
        """ Convert col/row (domain-based) into lon/lat """
        return self.londim[col], self.latdim[row]

    def lonlat2colrow(self, lon, lat):
        """ Find nearest tile (col/row) from any given lon/lat """
        londif = np.abs(self.londim - lon)
        latdif = np.abs(self.latdim - lat)
        lon = np.where(np.abs(londif-londif.min())<0.0001)[0][0]
        lat = np.where(np.abs(latdif-latdif.min())<0.0001)[0][0]
        return lon, lat

    def lonlat2tilenum(self, lon, lat):
        """" Get the (domain-based) tile number for a given lon/lat """
        col, row = self.lonlat2colrow(lon, lat)
        tilenum = np.where((self.tilecoord['i_indg'] == col)&
                           (self.tilecoord['j_indg'] == row))[0][0]
        return tilenum

