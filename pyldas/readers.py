
import os

import numpy as np
import pandas as pd
import xarray as xr

from struct import unpack
from netCDF4 import Dataset, date2num
from collections import OrderedDict

from pyldas.templates import get_template
from pyldas.functions import find_files, walk_up_folder
from pyldas.constants import paths

def s(line):
    return line[1:-2]

def b(line):
    return True if line[-2] == 'T' else False

class LDAS_io(object):

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
            self.files = find_files(paths().exp_root, param)
            if self.files is None:
                print 'No files for parameter: "' + param + '".'
                return

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
    def read_obsparam(fname):

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

    def read_fortran_binary(self, fname, dtype, hdr=None, length=None, reg_ftags=True, idx=None):

        if not os.path.isfile(fname):
            print 'file "', fname, '" not found.'
            return None

        fid = open(fname, 'rb')

        if hdr is not None:
            hdr = fid.read(4 * hdr)
            if length is None:
                # read header, assumed to be int32 after the fortran tag
                length = unpack('i', hdr[4:8][::-1])[0]
        else:
            if length is None:
                length = len(self.tilecoord)

        data = pd.DataFrame(columns=dtype.names, index=np.arange(length))

        if reg_ftags is True:
            for dt in dtype.names:
                fid.seek(4, 1)  # skip fortran tag
                data.loc[:, dt] = np.fromfile(fid, dtype=dtype[dt], count=length)
                fid.seek(4, 1)  # skip fortran tag

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

    def read_tilegrids(self, fname=None):
        if fname is None:
            fname = find_files(paths().rc_out, 'tilegrids')

        dtype, hdr, length = get_template('tilegrids')

        data = self.read_fortran_binary(fname, dtype, length=length, reg_ftags=False)

        data['gridtype'] = np.char.strip(data['gridtype'].values.astype('str'))
        data.index = ['global', 'domain']

        return data

    def read_tilecoord(self, fname=None):

        if fname is None:
            fname = find_files(paths().rc_out, 'tilecoord')

        dtype, hdr, length = get_template('tilecoord')

        return self.read_fortran_binary(fname, dtype, hdr=hdr)


    def read_image(self, yr, mo, da, hr, mi, species=None):

        if hasattr(self, 'images'):
            datestr = '%04i-%02i-%02i %02i:%02i' % (yr, mo, da, hr, mi)
            if species is not None:
                img = self.images.sel(species=species, time=datestr).values
            else:
                img = self.images.sel(time=datestr).values
            return img

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

            return self.read_fortran_binary(fname, dtype, hdr=hdr, length=length)

    @staticmethod
    def ncfile_init(fname, dimensions, variables):

        ds = Dataset(fname, mode='w')
        timeunit = 'hours since 2000-01-01 00:00'

        # initialize dimensions
        chunksizes = []
        for key, values in dimensions.iteritems():
            if key == 'time':
                values = date2num(values.to_pydatetime(), timeunit).astype('int32')
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

        # initialize variables
        for var in variables:
            ds.createVariable(var, 'float32',
                              dimensions=dimensions.keys(),
                              chunksizes=chunksizes,
                              fill_value=-9999.,
                              zlib=True)

        return ds

    def bin2netcdf(self):

        out_path = walk_up_folder(self.files[0],3)
        out_file = os.path.join(out_path,'.'.join(os.path.basename(self.files[0]).split('.')[:-2]) + '_images.nc')

        variables = get_template(self.param)[0].names

        lons = np.sort(self.tilecoord.groupby('i_indg').first()['com_lon'])
        lats = np.sort(self.tilecoord.groupby('j_indg').first()['com_lat'])[::-1]
        dates = self.dates

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
                tmp_img = data[var].values.copy()
                np.place(tmp_img, np.isnan(tmp_img), -9999.)

                if self.param == 'ObsFcstAna':
                    img[ind_spc,ind_lat,ind_lon] = tmp_img
                    dataset.variables[var][:,:,:,i] = img
                else:
                    img[ind_lat,ind_lon] = tmp_img
                    dataset.variables[var][:,:,i] = img

        dataset.close()


if __name__ == '__main__':

    obj = LDAS_io('xhourly')
    obj.bin2netcdf()

