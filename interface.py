
import os

import numpy as np
import pandas as pd
import xarray as xr

from netCDF4 import Dataset, date2num
from collections import OrderedDict

from pyldas.grids import EASE2

from pyldas.templates import get_template
from myprojects.functions import find_files, walk_up_folder
from pyldas.paths import paths

def s(line):
    """" Remove quotes from string """
    return line[1:-2]

def b(line):
    """" Turn 'T' and 'F' into True and False, respectively"""
    return True if line[-2] == 'T' else False

class LDAS_io(object):
    """
    Class for reading and writing LDAS specific data
    Default paths are taken from pygldas.constants

    Parameters
    ----------
    param : str
        Name of the parameter for which data should be read
    exp : string
        experiment name (appended to root path)
    domain : string
        domain name (appended to experiment path)

    Attributes
    ----------
    paths : dict
        Dictionary holding the path information for the specified exp/domain
    obsparam : pd.DataFrame
        Metadata about observations
    tilecoord : pd.DataFrame
        Metadata information for the tile coordinates of a particular LDAS experiment
    tilegrids : pd.DataFrame
        Metadata information for the tile grids of a particular LDAS experiment
    param : str
        Name of the parameter for which data is loaded
    files : np.array
        Array containing all names within the specified experiment directory, that match
        the specified parameter
    dates : pd.DatetimeIndex
        Dates corresponding to the files in self.files
    images : xr.Dataset
        netCDF image-chunked image stack (if available)
    timeseries : xr.Dataset
        netCDF timeseries-chunked image stack (if available)

    """

    def __init__(self,
                 param=None,
                 exp=None,
                 domain=None):

        self.paths = paths(exp=exp, domain=domain)

        self.obsparam = self.read_obsparam()
        self.tilecoord = self.read_params('tilecoord')
        self.tilegrids = self.read_params('tilegrids')

        self.grid = EASE2(tilecoord=self.tilecoord, tilegrids=self.tilegrids)

        self.param = param
        if param is not None:

            if param == 'xhourly':
                path = self.paths.__getattribute__('cat')
            else:
                path = self.paths.__getattribute__('exp_root')

            self.files = find_files(path, param)

            if self.files[0].find('images.nc') == -1:
                print 'NetCDF image cube not yet created. Use method "bin2netcdf".'
                self.dates = pd.to_datetime([f[-18:-5] for f in self.files], format='%Y%m%d_%H%M').sort_values()

                # TODO: Currently valid for 3-hourly data only! Times of the END of the 3hr periods are assigned!
                # if self.param == 'xhourly':
                    # self.dates += pd.to_timedelta('2 hours')

                self.dtype, self.hdr, self.length = get_template(self.param)

            else:
                self.images = xr.open_dataset(self.files[0])
                if self.files[1].find('timeseries.nc') == -1:
                    print 'NetCDF time series cube not yet created. Use the NetCDF kitchen sink.'
                else:
                    self.timeseries = xr.open_dataset(self.files[1])


    def read_obsparam(self):
        """ Read the 'obsparam' file. """

        fp = open(find_files(self.paths.rc_out, 'obsparam'))

        lines = fp.readlines()[1::]
        n_lines = len(lines)

        # 30 or 32 fields (before and after two entries for the use of uncertainty maps)
        n_fields = 32 if n_lines == 128 else 30

        n_blocks = n_lines / n_fields

        res = []
        for bl in np.arange(n_blocks) * n_fields:
            if n_fields == 32:
                res.append({'descr': s(lines[bl + 0]),
                            'species': int(lines[bl + 1]),
                            'orbit': int(lines[bl + 2]),
                            'pol': int(lines[bl + 3]),
                            'N_ang': int(lines[bl + 4]),
                            'ang': float(lines[bl + 5]),
                            'freq': float(lines[bl + 6]),
                            'FOV': float(lines[bl + 7]),
                            'FOV_units': s(lines[bl + 8]),
                            'assim': b(lines[bl + 9]),
                            'scale': b(lines[bl + 10]),
                            'getinnov': b(lines[bl + 11]),
                            'RTM_ID': int(lines[bl + 12]),
                            'bias_Npar': int(lines[bl + 13]),
                            'bias_trel': int(lines[bl + 14]),
                            'bias_tcut': int(lines[bl + 15]),
                            'nodata': float(lines[bl + 16]),
                            'varname': s(lines[bl + 17]),
                            'units': s(lines[bl + 18]),
                            'path': s(lines[bl + 19]),
                            'name': s(lines[bl + 20]),
                            'scalepath': s(lines[bl + 21]),
                            'scalename': s(lines[bl + 22]),
                            'errstd': float(lines[bl + 23]),
                            'errstd_file': b(lines[bl + 24]),
                            'path_errstd': s(lines[bl + 25]),
                            'std_normal_max': float(lines[bl + 26]),
                            'zeromean': b(lines[bl + 27]),
                            'coarsen_pert': b(lines[bl + 28]),
                            'xcorr': float(lines[bl + 29]),
                            'ycorr': float(lines[bl + 30]),
                            'adapt': int(lines[bl + 31])})
            else:
                res.append({'descr': s(lines[bl + 0]),
                            'species': int(lines[bl + 1]),
                            'orbit': int(lines[bl + 2]),
                            'pol': int(lines[bl + 3]),
                            'N_ang': int(lines[bl + 4]),
                            'ang': float(lines[bl + 5]),
                            'freq': float(lines[bl + 6]),
                            'FOV': float(lines[bl + 7]),
                            'FOV_units': s(lines[bl + 8]),
                            'assim': b(lines[bl + 9]),
                            'scale': b(lines[bl + 10]),
                            'getinnov': b(lines[bl + 11]),
                            'RTM_ID': int(lines[bl + 12]),
                            'bias_Npar': int(lines[bl + 13]),
                            'bias_trel': int(lines[bl + 14]),
                            'bias_tcut': int(lines[bl + 15]),
                            'nodata': float(lines[bl + 16]),
                            'varname': s(lines[bl + 17]),
                            'units': s(lines[bl + 18]),
                            'path': s(lines[bl + 19]),
                            'name': s(lines[bl + 20]),
                            'scalepath': s(lines[bl + 21]),
                            'scalename': s(lines[bl + 22]),
                            'errstd': float(lines[bl + 23]),
                            'std_normal_max': float(lines[bl + 24]),
                            'zeromean': b(lines[bl + 25]),
                            'coarsen_pert': b(lines[bl + 26]),
                            'xcorr': float(lines[bl + 27]),
                            'ycorr': float(lines[bl + 28]),
                            'adapt': int(lines[bl + 29])})

        return pd.DataFrame(res)

    def get_species(self, pol=None, orbit=None, ang=None):

        pol = 1 if pol == 'H' else 2
        orbit = 1 if orbit == 'A' else 2

        return self.obsparam.loc[(self.obsparam['pol'] == pol) & \
                                 (self.obsparam['orbit'] == orbit) & \
                                 (self.obsparam['ang'] == ang),'species'].values[0]

    def read_fortran_binary(self, fname, dtype,
                            hdr=None,
                            length=None,
                            reg_ftags=True,
                            loc=None,
                            idx=None):
        """
        Class for reading fortran binary files

        Parameters
        ----------
        fname : str
            Name of the file to be read
        dtype : np.dtype
            Template holding parameter names and data types of the file entries
        hdr : int
            Number of (4-byte) header entries to be skipped
        length : int
            If provided together with 'hdr':
                Position of file length information within the header
            else:
                Number of successive data blocks contained in the file
        reg_ftags : Boolean
            If True, a fortran tag (byte) is expected before and after each data field, otherwise
            only before and after each data block (only the case for tilegrids files)
        loc : int
            read only the <loc>-th data element of the file. (only works if reg_ftags is True)
            TODO: currently rather useless for ObsFcstAna, because file-lengths don't make much sense
        idx : str
            If specified, the data of the field named 'idx' will be used as index of the output data frame

        Returns
        -------
        fid : file-id (if return_hdr is True)
            The file ID of the opened file
        hdr : list (if return_hdr is True)
            The file header as list of integers

        data : pd.DataFrame (if return_hdr is False)
            Content of the fortran binary file

        """

        if not os.path.isfile(fname):
            print 'file "', fname, '" not found.'
            return None

        fid = open(fname, 'rb')

        # read header
        if hdr is not None:
            hdr = np.fromfile(fid, dtype='>i4', count=hdr).byteswap().newbyteorder()

            if length is not None:
                length = hdr[length]
            else:
                length = hdr[1]
        else:
            if length is None:
                length = len(self.tilecoord)

        if loc is None:
            data = pd.DataFrame(columns=dtype.names, index=np.arange(length))
        else:
            data = pd.DataFrame(columns=dtype.names, index=(loc,))

        if length==0:
            # print 'Empty file'
            return data

        # read data
        if reg_ftags is True:
            for dt in dtype.names:
                if loc is None:
                    fid.seek(4, 1)  # skip fortran tag
                    data.loc[:, dt] = np.fromfile(fid, dtype=dtype[dt], count=length).byteswap().newbyteorder()
                    fid.seek(4, 1)  # skip fortran tag
                else:
                    fid.seek(4 + 4*loc, 1)
                    data.loc[:, dt] = np.fromfile(fid, dtype=dtype[dt], count=1).byteswap().newbyteorder()
                    fid.seek(4 + 4*length - 4*loc - 4, 1)

        else:
            for i in np.arange(length):
                fid.seek(4, 1)  # skip fortran tag
                for dt in dtype.names:
                    data.loc[i, dt] = np.fromfile(fid, dtype=dtype[dt], count=1)[0]
                fid.seek(4, 1)  # skip fortran tag

        fid.close()

        if idx is not None:
            data.index = data.loc[:, idx].values
            data.drop(idx, axis='columns', inplace=True)

        return data

    def read_params(self, param, fname=None):
        """ Read parameter files (tilegrids, tilecoord, RTMparam, catparam"""

        if fname is None:
            fname = find_files(self.paths.rc_out, param)

        reg_ftags = False if param == 'tilegrids' else True

        dtype, hdr, length = get_template(param)
        data = self.read_fortran_binary(fname, dtype, hdr=hdr, length=length, reg_ftags=reg_ftags)
        data.replace(-9999., np.nan, inplace=True)

        if param == 'tilegrids':
            data.index = ['global', 'domain']
        else:
            # index equals the 'tilenum' which starts at 1!!
            data.index += 1

        return data


    def read_scaling_parameters(self, pentad=1, fname=None, tile_id=None):
        """
        Class for reading scaling files. These hold the observation and model mean and standard deviation, and
        number of observations per per pentade.

        Parameters
        ----------
        pentad : int
            If no tile-id / filename is provided, only one pentad will be read
        fname : str
            The path to a scaling file that should be read
        tile_id : int
            If provided, the scaling files for all pentads will be read for the specified tile_id

        Returns
        -------
        data : pd.DataFrame
            Mean, Std.dev, and N_data for model forecasts and observations.

        """
        # TODO: WRONG TREATMENT OF ORBIT DIRECTION! A/D IS IN THE FILENMAE BEFORE THE PENTADE!

        dtype, hdr, length = get_template('scaling')

        if tile_id is None:
            if fname is not None:
                data = self.read_fortran_binary(fname, dtype, hdr=hdr, length=length)
            else:
                data = self.read_fortran_binary(self.files[pentad], dtype, hdr=hdr, length=length)

        else:
            pentads = np.arange(73)+1

            fields = dtype.names[3:]
            data = pd.DataFrame(columns=fields, index=pentads)

            for pentad in pentads:
                tmp_data = self.read_fortran_binary(self.files[pentad], dtype, hdr=hdr, length=length, idx='tile_id')
                data.loc[pentad,fields] = tmp_data.loc[tile_id, fields].values

        return data

    def read_image(self, yr, mo, da, hr, mi, species=None):
        """"
        Read an image for a given date/time(/species)
        If a netCDF file has been created with self.bin2netcdf, the image will be read from this file,
        otherwise it is read from the fortran binary file

        Returns
        -------
        img : np.ndarray (if netCDF has been created)
            Complete 2D or 3D image (with/without species, depending on the specified parameter)
              pd.DataFrame (if read from fortran binary)
            Only tiles for which data is available

        """

        # If netCDF file has been created/loaded, use xarray indexing functions
        if hasattr(self, 'images'):
            datestr = '%04i-%02i-%02i %02i:%02i' % (yr, mo, da, hr, mi)
            if species is not None:
                img = self.images.sel(species=species, time=datestr).values
            else:
                img = self.images.sel(time=datestr).values

        # Otherwise, read from fortran binary
        else:
            datestr = '%04i%02i%02i_%02i%02i' % (yr, mo, da, hr, mi)
            fname = [f for f in self.files if f.find(datestr) != -1]

            if len(fname) == 0:
                print 'No files found for: "' + datestr + '".'
                return None
            elif len(fname) > 1:
                print 'Multiple files found for: "' + datestr + '".'
            else:
                fname = fname[0]

            img = self.read_fortran_binary(fname, self.dtype, hdr=self.hdr, length=self.length)

        return img

    def read_ts(self, param, col, row, species=None, lonlat=True):

        if lonlat is True:
            col, row = self.grid.lonlat2colrow(col, row, domain=True)

        if species is None:
            ts = self.timeseries[param][row,col,:].to_series()
        else:
            ts = self.timeseries[param].sel(species=species)[row,col,:].to_series()

        return ts

    @staticmethod
    def write_fortran_block(fid, data):
        """" Writes a data block (1D numpy array) into a binary file including fortran tags """

        # force 32 bit if 64 bit (usually for float)
        dtype = data.dtype.newbyteorder('>')
        if dtype.str[2]=='8':
            dtype = np.dtype(dtype.str[0:2]+'4')

        ftag = data.size * dtype.itemsize
        np.array(ftag).astype('>i4').tofile(fid)
        data.astype(dtype).tofile(fid)
        np.array(ftag).astype('>i4').tofile(fid)


    @staticmethod
    def ncfile_init(fname, dimensions, variables):
        """"
        Method to initialize dimensions/variables of a image-chunked netCDF file

        Parameters
        ----------
        fname : str
            Filename of the netCDF file to be created
        dimensions : dict
            Dictionary containing the dimension names and values
        variables : list
            list of variables to be created with the specified dimensions

        Returns
        -------
        ds : fileid
            File ID of the created netCDF file

        """

        ds = Dataset(fname, mode='w')
        timeunit = 'hours since 2000-01-01 00:00'

        # Initialize dimensions
        chunksizes = []
        for key, values in dimensions.iteritems():

            # convert pandas Datetime Index to netCDF-understandable numeric format
            if key == 'time':
                values = date2num(values.to_pydatetime(), timeunit).astype('int32')

            # Files are per default image chunked
            if key in ['lon','lat']:
                chunksize = len(values)
            else:
                chunksize = 1
            chunksizes.append(chunksize)

            dtype = values.dtype
            ds.createDimension(key, len(values))
            ds.createVariable(key,dtype,
                              dimensions=(key,),
                              chunksizes=(chunksize,),
                              zlib=True)
            ds.variables[key][:] = values
        ds.variables['time'].setncattr('units',timeunit)

        # Initialize variables
        for var in variables:
            ds.createVariable(var, 'float32',
                              dimensions=dimensions.keys(),
                              chunksizes=chunksizes,
                              fill_value=-9999.,
                              zlib=True)

        return ds

    def bin2netcdf(self):
        """" Convert fortran binary image into a netCDF data cube """

        out_path = walk_up_folder(self.files[0],3)
        out_file = os.path.join(out_path,self.param + '_images.nc')

        # get variable names from fortran reader template
        variables = get_template(self.param)[0].names

        lons = np.sort(self.tilecoord.groupby('i_indg').first()['com_lon'])
        lats = np.sort(self.tilecoord.groupby('j_indg').first()['com_lat'])[::-1]
        dates = self.dates

        # Innovation file data has an additional 'species' dimension
        if self.param == 'ObsFcstAna':
            # Remove dates which do not contain any data
            nodata = list()
            for i, dt in enumerate(dates):
                print '%d / %d' % (i, len(dates))
                data = self.read_image(dt.year, dt.month, dt.day, dt.hour, dt.minute)
                if len(data) == 0:
                    nodata.append(dt)

            dates = dates.drop(nodata)
            if len(dates) == 0:
                print 'Images do not contain valid data.'
                return

            spc = pd.DataFrame(self.obsparam)['species'].values.astype('uint8')
            dimensions = OrderedDict([('species',spc), ('lat',lats), ('lon',lons), ('time',dates)])
        else:
            dimensions = OrderedDict([('lat',lats), ('lon',lons), ('time',dates)])

        dataset = self.ncfile_init(out_file, dimensions, variables)

        for i,dt in enumerate(dates):
            print '%d / %d' % (i, len(dates))

            data = self.read_image(dt.year,dt.month,dt.day,dt.hour,dt.minute)
            if len(data) == 0:
                continue

            if self.param == 'ObsFcstAna':
                img = np.full((len(spc), len(lats), len(lons)), -9999., dtype='float32')
                ind_lat = self.tilecoord.loc[data['obs_tilenum'].values, 'j_indg'].values - self.tilegrids.loc['domain','j_offg']
                ind_lon = self.tilecoord.loc[data['obs_tilenum'].values, 'i_indg'].values - self.tilegrids.loc['domain','i_offg']
                ind_spc = data['obs_species'].values - 1

            else:
                img = np.full((len(lats),len(lons)), -9999., dtype='float32')
                ind_lat = self.tilecoord.loc[:, 'j_indg'].values - self.tilegrids.loc['domain','j_offg']
                ind_lon = self.tilecoord.loc[:, 'i_indg'].values - self.tilegrids.loc['domain','i_offg']

            for var in variables:
                # replace NaN values with the default -9999. fill Value
                tmp_img = data[var].values.copy()
                np.place(tmp_img, np.isnan(tmp_img), -9999.)

                if self.param == 'ObsFcstAna':
                    img[ind_spc,ind_lat,ind_lon] = tmp_img
                    dataset.variables[var][:,:,:,i] = img
                else:
                    img[ind_lat,ind_lon] = tmp_img
                    dataset.variables[var][:,:,i] = img

        # Save file to disk and loat it as xarray Dataset into the class variable space
        dataset.close()
        self.images = xr.open_dataset(out_file)


if __name__=='__main__':

    io = LDAS_io('ensstd', 'US_M36_SMOS40_DA_cal_scl_errfile_w_std')
    io.bin2netcdf()
    # io = LDAS_io('incr', 'US_M36_SMOS_DA_cal_scaled_yearly')
    # io.bin2netcdf()
    # io = LDAS_io('xhourly', 'US_M36_SMOS_DA_cal_scaled_yearly')
    # io.bin2netcdf()
    # io = LDAS_io('ObsFcstAna', 'US_M36_SMOS_DA_cal_scaled_yearly')
    # io.bin2netcdf()


    # io.read_ts('obs_obs', -113.480529785, 40.691051628, species=1).plot()
