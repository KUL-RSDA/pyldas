
import os

import numpy as np
import pandas as pd
import xarray as xr

from netCDF4 import Dataset, date2num
from collections import OrderedDict

from pyldas.templates import get_template
from pyldas.functions import find_files, walk_up_folder
from pyldas.constants import paths

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
    obsparam_path : str
        Path to the 'obsparam' file
    tilecoord_path : str
        Path to the 'tilecoord' file
    tilecoord_path : str
        Path to the 'tilegrids' file

    Attributes
    ----------
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
                 obsparam_path=None,
                 tilecoord_path=None,
                 tilegrids_path=None):

        self.obsparam = self.read_obsparam(fname=obsparam_path)
        self.tilecoord = self.read_tilecoord(fname=tilecoord_path)
        self.tilegrids = self.read_tilegrids(fname=tilegrids_path)

        self.param = param
        if param is not None:

            if param == 'scale':
                search_dir = 'scalefile_root'
            else:
                search_dir = 'exp_root'

            self.files = find_files(getattr(paths(),search_dir), param)
            if self.files is None:
                print 'No files for parameter: "' + param + '".'
                return

            if self.param != 'scale':
                if self.files[0].find('images.nc') == -1:
                    print 'NetCDF image cube not yet created. Use method "bin2netcdf".'
                    self.dates = pd.to_datetime([f[-18:-5] for f in self.files], format='%Y%m%d_%H%M')
                else:
                    self.images = xr.open_dataset(self.files[0])
                    if self.files[1].find('timeseries.nc') == -1:
                        print 'NetCDF time series cube not yet created. Use the NetCDF kitchen sink.'
                    else:
                        self.timeseries = xr.open_dataset(self.files[1])



    @staticmethod
    def read_obsparam(fname=None):
        """
        Class for reading the 'obsparam' text file holding metadata on the observations.

        Parameters
        ----------
        fname : str
            name of the 'obsparam' file

        Returns
        -------
        res : pd.DataFrame
            Metadata about observations

        """

        if fname is None:
            fname = find_files(paths().rc_out, 'obsparam')

        fp = open(fname)

        lines = fp.readlines()[1::]

        n_fields = 30
        n_blocks = len(lines) / n_fields

        res = []
        for bl in np.arange(n_blocks) * n_fields:
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


    def read_fortran_binary(self, fname, dtype,
                            hdr=None,
                            length=None,
                            reg_ftags=True,
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
        length : str
            If provided together with 'hdr':
                Position of file length information within the header
            else:
                Number of successive data blocks contained in the file
        reg_ftags : Boolean
            If True, a fortran tag (byte) is expected before and after each data field, otherwise
            only before and after each data block (only the case for tilegrids files)
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

        if hdr is not None:
            hdr = np.fromfile(fid, dtype='>i4', count=hdr).byteswap().newbyteorder()

            if length is not None:
                # If hdr & length are specified, length refers to the n-th field in hdr,
                # which contains the number of blocks to be read
                length = hdr[length]
            else:
                # If only hdr is specified, the second field refers to the field
                # which contains the number of blocks to be read
                length = hdr[1]
        else:
            if length is None:
                # If neither length nor hdr are specified, a data block for each (domain) tile is assumed
                length = len(self.tilecoord)

        data = pd.DataFrame(columns=dtype.names, index=np.arange(length))

        if reg_ftags is True:

            for dt in dtype.names:
                fid.seek(4, 1)  # skip fortran tag
                data.loc[:, dt] = np.fromfile(fid, dtype=dtype[dt], count=length).byteswap().newbyteorder()
                fid.seek(4, 1)  # skip fortran tag

        else:
            for i in np.arange(length):
                fid.seek(4, 1)  # skip fortran tag
                for dt in dtype.names:
                    data.loc[i, dt] = np.fromfile(fid, dtype=dtype[dt], count=1)[0].byteswap().newbyteorder()
                fid.seek(4, 1)  # skip fortran tag

        fid.close()

        if idx is not None:
            data.index = data.loc[:, idx].values
            data.drop(idx, axis='columns', inplace=True)

        return data


    def read_tilegrids(self, fname=None):
        """ Read the 'tilegrids' file. """

        if fname is None:
            fname = find_files(paths().rc_out, 'tilegrids')

        dtype, hdr, length = get_template('tilegrids')

        data = self.read_fortran_binary(fname, dtype, length=length, reg_ftags=False)

        # remove leading/trailing blanks from the txt file
        data['gridtype'] = np.char.strip(data['gridtype'].values.astype('str'))
        data.index = ['global', 'domain']

        return data


    def read_tilecoord(self, fname=None):
        """ Read the 'tilecoords' file. """

        if fname is None:
            fname = find_files(paths().rc_out, 'tilecoord')

        dtype, hdr, length = get_template('tilecoord')

        return self.read_fortran_binary(fname, dtype, hdr=hdr)


    def read_scaling_parameters(self, pentad=1, tile_id=None, angles=(40,)):
        """ Read the scaling files. """

        dtype, hdr, length = get_template('scaling')

        if tile_id is None:
            data = self.read_fortran_binary(self.files[pentad], dtype, hdr=hdr, length=length, idx='tile_id')

        else:
            pentads = np.arange(73)+1
            fields = ['m_obs_H_%i' % ang for ang in angles] + \
                     ['m_obs_V_%i' % ang for ang in angles] + \
                     ['m_mod_H_%i' % ang for ang in angles] + \
                     ['m_mod_V_%i' % ang for ang in angles]
                     # ['N_data_H_%i' % ang for ang in angles] + \
                     # ['N_data_V_%i' % ang for ang in angles]

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

            dtype, hdr, length = get_template(self.param)
            img = self.read_fortran_binary(fname, dtype, hdr=hdr, length=length)

        return img

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
                chunksize = 10
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
        out_file = os.path.join(out_path,'.'.join(os.path.basename(self.files[0]).split('.')[:-2]) + '_images.nc')

        # get variable names from fortran reader template
        variables = get_template(self.param)[0].names

        lons = np.sort(self.tilecoord.groupby('i_indg').first()['com_lon'])
        lats = np.sort(self.tilecoord.groupby('j_indg').first()['com_lat'])[::-1]
        dates = self.dates

        # Innovation file data has an additional 'species' dimension
        if self.param == 'ObsFcstAna':
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


# if __name__=='__main__':
#     io = LDAS_io('scale')
#
#     tile_id = 107300
#     data = io.read_scaling_parameters(tile_id=tile_id)
#     print data
#     data.plot()
#     import matplotlib.pyplot as plt
#     plt.show()