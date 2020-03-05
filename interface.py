
import logging

import numpy as np
import pandas as pd
import xarray as xr

from pathlib import Path

from netCDF4 import Dataset, date2num
from collections import OrderedDict

from pyldas.grids import EASE2

from pyldas.templates import get_template
from pyldas.paths import paths

logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

def s(line):
    """" Remove quotes from string """
    return line[1:-2]

def b(line):
    """" Turn 'T' and 'F' into True and False, respectively"""
    return True if line[-2] == 'T' else False

class LDAS_io(object):
    """
    Class for reading and writing LDAS specific data
    Default paths are taken from pyldas.paths

    Parameters
    ----------
    param : str
        Name of the parameter for which data should be read
    exp : string
        experiment name (appended to root path)
    domain : string
        domain name (appended to experiment path)
    root : pathlib.Path
        root path to the experiment directory

    Attributes
    ----------
    paths : dict
        Dictionary holding the path information for the specified exp/domain
    obsparam : pd.DataFrame
        Metadata about observations
    param : str
        Name of the parameter for which data is loaded
    files : np.array
        Array containing all names within the specified experiment directory, that match
        the specified parameter (excluding netcdf files, if already created)
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
                 domain=None,
                 root=None):

        self.paths = paths(exp=exp, domain=domain, root=root)

        try:
            self.obsparam = self.read_obsparam()
        except:
            logging.info('No obsparam file. This is a model-only (No DA) run.')

        tilecoord = self.read_params('tilecoord')
        tilegrids = self.read_params('tilegrids')
        self.grid = EASE2(tilecoord, tilegrids)

        self.param = param
        if param is not None:

            if param == 'xhourly':
                path = self.paths.cat
            else:
                path = self.paths.exp_root

            nc_file = list(path.glob('**/*' + param + '_images.nc'))
            if len(nc_file) == 0:
                logging.warning('NetCDF image cube not yet created. Use method "bin2netcdf".')
            else:
                self.images = xr.open_dataset(nc_file[0])

            nc_file = list(path.glob('**/*' + param + '_timeseries.nc'))
            if len(nc_file) == 0:
                logging.warning('NetCDF time series cube not yet created. Use the NetCDF kitchen sink.')
            else:
                self.timeseries = xr.open_dataset(nc_file[0])

            if param == 'ObsFcstAna':
                self.files = np.sort(list(path.glob('**/*' + param + '.*.bin')))
            else:
                self.files = np.sort(list(path.glob('**/*' + param + '*.bin')))

            if param == 'hscale':
                self.pentads = np.array([f.name[-6:-4] for f in self.files]).astype('int')
                self.orbits = np.array([f.name[-9:-8] for f in self.files])
            else:
                self.dates = pd.to_datetime([f.name[-18:-5] for f in self.files], format='%Y%m%d_%H%M')

            if param == 'ObsFcstAnaEns':
                self.ens_id = np.array([f.name[-42:-38] for f in self.files]).astype('int')

            # TODO: Currently valid for 3-hourly data only! Times of the END of the 3hr periods are assigned!
            # if self.param == 'xhourly':
                # self.dates += pd.to_timedelta('2 hours')

            self.dtype, self.hdr, self.length = get_template(self.param)


    def read_obsparam(self):
        """ Read the 'obsparam' file. """

        fp = open(list(self.paths.rc_out.glob('**/*obsparam*'))[0])

        lines = fp.readlines()[1::]
        n_lines = len(lines)

        # 30 or 32 fields (before and after two entries for the use of uncertainty maps)
        if n_lines == 128 or n_lines == 896:
            n_fields = 32
        else:
            n_fields = 30

        n_blocks = n_lines / n_fields

        # different output scenarios.
        res = []
        for bl in np.arange(n_blocks, dtype='int') * n_fields:
            if n_fields == 32 and n_blocks == 28:
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
                            'flistpath': s(lines[bl + 23]),
                            'flistname': s(lines[bl + 24]),
                            'errstd': float(lines[bl + 25]),
                            'std_normal_max': float(lines[bl + 26]),
                            'zeromean': b(lines[bl + 27]),
                            'coarsen_pert': b(lines[bl + 28]),
                            'xcorr': float(lines[bl + 29]),
                            'ycorr': float(lines[bl + 30]),
                            'adapt': int(lines[bl + 31])})
            elif n_fields == 32 and n_blocks == 4:
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
        data : pd.DataFrame (if return_hdr is False)
            Content of the fortran binary file

        """

        if not Path(fname).exists():
            logging.warning('file "', fname, '" not found.')
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
                length = len(self.grid.tilecoord)

        if loc is None:
            data = pd.DataFrame(columns=dtype.names, index=np.arange(length))
        else:
            data = pd.DataFrame(columns=dtype.names, index=(loc,))

        if length==0:
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
            fname = list(self.paths.rc_out.glob('**/*' + param + '*'))[0]

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


    def read_scaling_parameters(self, pentad=1, fname=None, tile_id=None, sensor='SMAP'):
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
        # TODO: WRONG TREATMENT OF ORBIT DIRECTION! A/D IS IN THE FILENAME BEFORE THE PENTADE!

        dtype, hdr, length = get_template('hscale', sensor=sensor)

        if tile_id is None:
            if fname is not None:
                data = self.read_fortran_binary(fname, dtype, hdr=hdr, length=length)
            else:
                data = self.read_fortran_binary(self.files[pentad], dtype, hdr=hdr, length=length)
            # increase index to match tilecoord indices
            data.index += 1

        else:
            pentads = np.arange(73)+1

            fields = dtype.names[3:]
            data = pd.DataFrame(columns=fields, index=pentads)

            for pentad in pentads:
                tmp_data = self.read_fortran_binary(self.files[pentad], dtype, hdr=hdr, length=length, idx='tile_id')
                data.loc[pentad,fields] = tmp_data.loc[tile_id, fields].values

        return data

    def read_image(self, yr=None, mo=None, da=None, hr=None, mi=None, species=None, ens_id=None, fname=None):
        """"
        Read an image for a given date/time(/species/subregion)
        If a netCDF file has been created with self.bin2netcdf, the image will be read from this file,
        otherwise it is read from the fortran binary file

        Parameters
        ----------
        yr, mo, da, hr, mi : int
            Date for which the image should be read.
        species : int
            If provided, only a specific species will be read
            No effect when reading fortran binaries!
        ens_id : int
            If provided, only a specific ensemble member will be read
            No effect when reading fortran binaries!
        fname: string
            Direct input of the filename (instead of inferring from date/time information)

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
            tmp_img = self.images
            if species is not None:
                tmp_img = tmp_img.sel(species=species)
            if ens_id is not None:
                tmp_img = tmp_img.sel(ens_id=ens_id)
            img = tmp_img.sel(time=datestr).values

        # Otherwise, read from fortran binary
        else:
            if fname is None:
                datestr = '%04i%02i%02i_%02i%02i' % (yr, mo, da, hr, mi)
                fname = [f for f in self.files if f.name.find(datestr) != -1]

                if len(fname) == 0:
                    logging.warning('No files found for: "' + datestr + '".')
                    return None
                elif len(fname) > 1:
                    logging.warning('Multiple files found for: "' + datestr + '".')
                else:
                    fname = fname[0]

            img = self.read_fortran_binary(fname, self.dtype, hdr=self.hdr, length=self.length)

            # set index to match tilecoord indices
            if 'obs_tilenum' in img:
                # ObsFcstAna files
                img.index = img['obs_tilenum'].values
            else:
                # All other files (hopefully)
                img.index += 1

        return img

    def read_ts(self, param, col, row, species=None, ens_id=None, lonlat=True):
        """ Reads a time series from the netCDF time series chunked cube. """

        if lonlat is True:
            col, row = self.grid.lonlat2colrow(col, row, domain=True)

        tmp_ts = self.timeseries[param]
        if species is not None:
            tmp_ts = tmp_ts.sel(species=species)
        if ens_id is not None:
            tmp_ts = tmp_ts.sel(ens_id=ens_id)

        return tmp_ts.isel(lat=row,lon=col).to_pandas()


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
        for dim in dimensions:

            # convert pandas Datetime Index to netCDF-understandable numeric format
            if dim == 'time':
                dimensions[dim] = date2num(dimensions[dim].to_pydatetime(), timeunit).astype('int32')

            # Files are per default image chunked
            if dim in ['lon','lat']:
                chunksize = len(dimensions[dim])
            else:
                chunksize = 1
            chunksizes.append(chunksize)

            dtype = dimensions[dim].dtype
            ds.createDimension(dim, len(dimensions[dim]))
            ds.createVariable(dim,dtype,
                              dimensions=(dim,),
                              chunksizes=(chunksize,),
                              zlib=True)
            ds.variables[dim][:] = dimensions[dim]

        # Coordinate attributes following CF-conventions

        if 'time' in dimensions:
            ds.variables['time'].setncatts({'long_name': 'time',
                                            'units': timeunit})
        ds.variables['lon'].setncatts({'long_name': 'longitude',
                                       'units':'degrees_east'})
        ds.variables['lat'].setncatts({'long_name': 'latitude',
                                        'units':'degrees_north'})

        # Initialize variables
        for var in variables:
            ds.createVariable(var, 'float32',
                              dimensions=list(dimensions.keys()),
                              chunksizes=chunksizes,
                              fill_value=-9999.,
                              zlib=True)

        return ds

    def bin2netcdf(self,
                   overwrite=False,
                   date_from=None,
                   date_to=None,
                   latmin=-90.,
                   latmax=90.,
                   lonmin=-180.,
                   lonmax=180.,
                   out_file=None):

        """"
        Convert fortran binary image into a netCDF data cube.

        Parameters
        ----------
        overwrite : boolean
            If set, an already existing netCDF file will be overwritten
        date_from : string
            Lower time limit for which a netCDF image cube should be generated (string format, e.g., '2010-01-01')
        date_to : string
            Upper time limit for which a netCDF image cube should be generated (string format, e.g., '2010-01-01')
        latmin : float
            Lower latitude limit for which a netCDF image cube should be generated
        latmax : float
            Upper latitude limit for which a netCDF image cube should be generated
        lonmin : float
            Lower longitude limit for which a netCDF image cube should be generated
        lonmax : float
            Upper longitude limit for which a netCDF image cube should be generated
        out_file : string
            Optional alternative path / filename for the created NetCDF image cube

        """

        if out_file is None:
            out_file = self.files[0].parents[3] / (self.param + '_images.nc')

        # remove file if it already exists
        if hasattr(self,'images'):
            if overwrite is False:
                logging.warning('bin2netcdf: NetCDF image file already exists. Use keyword "overwrite" to regenerate.')
                return
            else:
                delattr(self, 'images')
                out_file.unlink()

        # get variable names from fortran reader template
        if self.param == 'hscale':
            variables = get_template(self.param)[0].names[3::]
        else:
            variables = get_template(self.param)[0].names

        # If specified, only generate netCDF file for specific date range
        if self.param == 'hscale':
            pentads = np.unique(self.pentads)
            orbits = np.unique(self.orbits)
        elif self.param == 'ObsFcstAnaEns':
            files = self.files
            dates = self.dates
            ids = self.ens_id
            if date_from is not None:
                files = files[dates >= pd.to_datetime(date_from)]
                ids = ids[dates >= pd.to_datetime(date_from)]
                dates = dates[dates >= pd.to_datetime(date_from)]
            if date_to is not None:
                files = files[dates <= pd.to_datetime(date_to)]
                ids = ids[dates <= pd.to_datetime(date_to)]
                dates = dates[dates <= pd.to_datetime(date_to)]
        else:
            dates = self.dates
            if date_from is not None:
                dates = dates[dates >= pd.to_datetime(date_from)]
            if date_to is not None:
                dates = dates[dates <= pd.to_datetime(date_to)]

        domainlons = self.grid.ease_lons[np.min(self.grid.tilecoord.i_indg):(np.max(self.grid.tilecoord.i_indg)+1)]
        domainlats = self.grid.ease_lats[np.min(self.grid.tilecoord.j_indg):(np.max(self.grid.tilecoord.j_indg)+1)]

        lonmin = domainlons[np.argmin(np.abs(domainlons-lonmin))]
        lonmax = domainlons[np.argmin(np.abs(domainlons-lonmax))]
        latmin = domainlats[np.argmin(np.abs(domainlats-latmin))]
        latmax = domainlats[np.argmin(np.abs(domainlats-latmax))]

        # Use grid lon lat to avoid rounding issues
        tmp_tilecoord = self.grid.tilecoord.copy()
        tmp_tilecoord['com_lon'] = self.grid.ease_lons[self.grid.tilecoord.i_indg]
        tmp_tilecoord['com_lat'] = self.grid.ease_lats[self.grid.tilecoord.j_indg]

        # Clip region based on specified coordinate boundaries
        ind_img = self.grid.tilecoord[(tmp_tilecoord['com_lon']>=lonmin)&(tmp_tilecoord['com_lon']<=lonmax)&
                                 (tmp_tilecoord['com_lat']<=latmax)&(tmp_tilecoord['com_lat']>=latmin)].index
        lons = domainlons[(domainlons >= lonmin) & (domainlons <= lonmax)]
        lats = domainlats[(domainlats >= latmin) & (domainlats <= latmax)]
        i_offg_2 = np.where(domainlons >= lonmin)[0][0]
        j_offg_2 = np.where(domainlats <= latmax)[0][0]

        # Innovation file data has an additional 'species' dimension
        if self.param == 'ObsFcstAna':
            # Remove dates which do not contain any data
            nodata = list()
            for i, dt in enumerate(dates):
                logging.info('%d / %d' % (i, len(dates)))

                data = self.read_image(dt.year, dt.month, dt.day, dt.hour, dt.minute)
                data = data.loc[data.index.intersection(ind_img), :]

                if len(data) == 0:
                    nodata.append(dt)

            dates = dates.drop(nodata)
            if len(dates) == 0:
                logging.warning('Images do not contain valid data.')
                return

            spc = pd.DataFrame(self.obsparam)['species'].values.astype('uint8')
            dimensions = OrderedDict([('time', dates), ('species', spc), ('lat', lats), ('lon', lons)])

        # Innovation ensemble file data has an additional 'species' + 'ens_id' dimension
        elif self.param == 'ObsFcstAnaEns':
            # Remove dates which do not contain any data
            nodata = list()
            for i, fn in enumerate(files):
                logging.info('%d / %d' % (i, len(files)))

                data = self.read_image(fname=fn)
                data = data.loc[data.index.intersection(ind_img), :] # clip subregion

                if len(data) == 0:
                    nodata.append(i)

            files = np.delete(files, nodata)
            dates = np.delete(dates, nodata)
            ids = np.delete(ids, nodata)

            if len(files) == 0:
                logging.warning('Images do not contain valid data.')
                return

            udates = dates.unique()
            uids = np.unique(ids)

            spc = pd.DataFrame(self.obsparam)['species'].values.astype('uint8')
            dimensions = OrderedDict([('time', udates), ('ens_id', uids), ('species', spc), ('lat', lats), ('lon', lons)])
        elif self.param == 'hscale':
            dimensions = OrderedDict([('pentad', pentads), ('orbit', orbits), ('lat', lats), ('lon', lons)])
        else:
            dimensions = OrderedDict([('time', dates), ('lat', lats), ('lon', lons)])

        dataset = self.ncfile_init(out_file, dimensions, variables)

        if self.param == 'hscale':

            for i,(file, pentad, orbit) in enumerate(zip(self.files, self.pentads, self.orbits)):

                logging.info('%d / %d' % (i, len(self.files)))

                data = self.read_scaling_parameters(fname=file)
                data = data.loc[data.index.intersection(ind_img), :]

                if len(data) == 0:
                    continue

                img = np.full((1, 1, len(lats), len(lons)), -9999., dtype='float32')
                ind_pen = np.where(pentads==pentad)[0][0]
                ind_orb = np.where(orbits==orbit)[0][0]
                ind_lat = self.grid.tilecoord.loc[ind_img, 'j_indg'].values - self.grid.tilegrids.loc[
                    'domain', 'j_offg'] - j_offg_2
                ind_lon = self.grid.tilecoord.loc[ind_img, 'i_indg'].values - self.grid.tilegrids.loc[
                    'domain', 'i_offg'] - i_offg_2

                for var in variables:
                    # replace NaN values with the default -9999. fill Value
                    tmp_img = data[var].values.copy()
                    np.place(tmp_img, np.isnan(tmp_img), -9999.)
                    img[0, 0, ind_lat, ind_lon] = tmp_img
                    dataset.variables[var][ind_pen, ind_orb, :, :] = img

        elif self.param == 'ObsFcstAnaEns':

            for i,(fn,dt,ensid) in enumerate(zip(files,dates,ids)):

                logging.info('%d / %d' % (i, len(files)))

                data = self.read_image(fname=fn)
                data = data.loc[data.index.intersection(ind_img), :]

                if len(data) == 0:
                    continue

                ind_dt = np.where(udates==dt)[0][0]
                ind_id = ensid
                ind_spc = data['obs_species'].values - 1
                ind_lat = self.grid.tilecoord.loc[data['obs_tilenum'].values, 'j_indg'].values - self.grid.tilegrids.loc['domain','j_offg'] - j_offg_2
                ind_lon = self.grid.tilecoord.loc[data['obs_tilenum'].values, 'i_indg'].values - self.grid.tilegrids.loc['domain','i_offg'] - i_offg_2

                for j,var in enumerate(variables):
                    # replace NaN values with the default -9999. fill Value
                    img = np.full((len(spc), len(lats), len(lons)), -9999., dtype='float32')
                    tmp_img = data[var].values.copy()
                    np.place(tmp_img, np.isnan(tmp_img), -9999.)
                    img[ind_spc, ind_lat, ind_lon] = tmp_img
                    dataset.variables[var][ind_dt, ind_id, :, :, :] = img

        else:
            for i,dt in enumerate(dates):
                logging.info('%d / %d' % (i, len(dates)))

                data = self.read_image(dt.year, dt.month, dt.day, dt.hour, dt.minute)
                data = data.loc[data.index.intersection(ind_img), :]

                if len(data) == 0:
                    continue

                if self.param == 'ObsFcstAna':
                    img = np.full((len(spc), len(lats), len(lons)), -9999., dtype='float32')
                    ind_lat = self.grid.tilecoord.loc[data['obs_tilenum'].values, 'j_indg'].values - self.grid.tilegrids.loc['domain','j_offg'] - j_offg_2
                    ind_lon = self.grid.tilecoord.loc[data['obs_tilenum'].values, 'i_indg'].values - self.grid.tilegrids.loc['domain','i_offg'] - i_offg_2
                    ind_spc = data['obs_species'].values - 1

                else:
                    img = np.full((len(lats),len(lons)), -9999., dtype='float32')
                    ind_lat = self.grid.tilecoord.loc[ind_img, 'j_indg'].values - self.grid.tilegrids.loc['domain','j_offg'] - j_offg_2
                    ind_lon = self.grid.tilecoord.loc[ind_img, 'i_indg'].values - self.grid.tilegrids.loc['domain','i_offg'] - i_offg_2

                for var in variables:
                    # replace NaN values with the default -9999. fill Value
                    tmp_img = data[var].values.copy()
                    np.place(tmp_img, np.isnan(tmp_img), -9999.)

                    if self.param == 'ObsFcstAna':
                        img[ind_spc,ind_lat,ind_lon] = tmp_img
                        dataset.variables[var][i,:,:,:] = img
                    else:
                        img[ind_lat,ind_lon] = tmp_img
                        dataset.variables[var][i,:,:] = img

        # Save file to disk and loat it as xarray Dataset into the class variable space
        dataset.close()
        self.images = xr.open_dataset(out_file)


if __name__=='__main__':

    date_from = '2016-06-01'
    date_to = '2016-07-01'
    # date_from = None
    # date_to = None

    latmin=35.
    latmax=45.
    lonmin=-90.
    lonmax=-50.
    # latmin=0
    # latmax=90
    # lonmin=-180.
    # lonmax=0.

    param = 'xhourly'
    exp = 'US_M36_SMOS40_DA_cal_scaled'

    io = LDAS_io(param, exp)

    io.bin2netcdf(overwrite=True, date_from=date_from, date_to=date_to, latmin=latmin, latmax=latmax, lonmin=lonmin, lonmax=lonmax)

    # io.read_ts('obs_obs', -113.480529785, 40.691051628, species=1).plot()
