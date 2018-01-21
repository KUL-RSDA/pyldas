
import numpy as np

from ease_grid.ease2_grid import EASE2_grid

from pyldas.functions import find_files
from pyldas.constants import paths
from pyldas.readers import LDAS_io

class EASE2(EASE2_grid):

    def __init__(self, res=36, map_scale=None, tileinfo_path=None):

        res = res * 1000

        if map_scale is None:
            if res == 36000:
                map_scale = 36032.220840584
            elif res == 9000:
                map_scale = 9008.055210146
            elif res == 3000:
                map_scale = 3002.6850700487

        if tileinfo_path is None:
            tileinfo_path = paths().rc_out

        tilecoord = find_files(tileinfo_path,'tilecoord')
        tilegrids = find_files(tileinfo_path,'tilegrids')

        io = LDAS_io(tilecoord_path=tilecoord, tilegrids_path=tilegrids)
        self.tilecoord = io.tilecoord
        self.tilegrids = io.tilegrids

        super(EASE2, self).__init__(res, map_scale=map_scale)


    def colrow2lonlat(self, col, row):
        return self.londim[col], self.latdim[row]

    def lonlat2colrow(self, lon, lat):
        londif = np.abs(self.londim - lon)
        latdif = np.abs(self.latdim - lat)
        lon = np.where(np.abs(londif-londif.min())<0.0001)[0][0]
        lat = np.where(np.abs(latdif-latdif.min())<0.0001)[0][0]
        return lon, lat

    def lonlat2tilenum(self, lon, lat):
        col, row = self.lonlat2colrow(lon, lat)
        tilenum = np.where((self.tilecoord['i_indg'] == col)&
                           (self.tilecoord['j_indg'] == row))[0][0]
        return tilenum

