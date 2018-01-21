
import os

import pandas as pd

from struct import unpack

from pyldas.binary_templates import *
from pyldas.functions import find_files
from pyldas.constants import paths

def s(line):
    return line[1:-2]

def b(line):
    return True if line[-2] == 'T' else False

def read_obsparam(fname):

    if not os.path.isfile(fname):
        print 'file "', fname, '" not found.'
        return None

    fp = open(fname)

    lines = fp.readlines()[1::]

    n_fields = 30
    n_blocks = len(lines)/n_fields

    res = []
    for bl in np.arange(n_blocks)*n_fields:

        res.append({'descr': s([bl+0]),
                    'species': int(lines[bl+1]),
                    'orbit': int(lines[bl+2]),
                    'pol': int(lines[bl+3]),
                    'N_ang': int(lines[bl+4]),
                    'ang': float(lines[bl+5]),
                    'freq': float(lines[bl+6]),
                    'FOV': float(lines[bl+7]),
                    'FOV_units': s(lines[bl+8]),
                    'assim': b(lines[bl+9]),
                    'scale': b(lines[bl+10]),
                    'getinnov': b(lines[bl+11]),
                    'RTM_ID': int(lines[bl+12]),
                    'bias_Npar': int(lines[bl+13]),
                    'bias_trel': int(lines[bl+14]),
                    'bias_tcut': int(lines[bl+15]),
                    'nodata': float(lines[bl+16]),
                    'varname': s(lines[bl+17]),
                    'units': s(lines[bl+18]),
                    'path': s(lines[bl+19]),
                    'name': s(lines[bl+20]),
                    'scalepath': s(lines[bl+21]),
                    'scalename': s(lines[bl+22]),
                    'errstd': float(lines[bl+23]),
                    'std_normal_max': float(lines[bl+24]),
                    'zeromean': b(lines[bl+25]),
                    'coarsen_pert': b(lines[bl+26]),
                    'xcorr': float(lines[bl+27]),
                    'ycorr': float(lines[bl+28]),
                    'adapt': int(lines[bl+29])})

    return res

def read_fortran_binary(fname, dtype, hdr=None, length=None, reg_ftags=True, idx=None):

    if not os.path.isfile(fname):
        print 'file "', fname, '" not found.'
        return None

    fid = open(fname,'rb')

    if hdr is not None:
        hdr = fid.read(4*hdr)
        if length is None:
            # read header, assumed to be int32 after the fortran tag
            length = unpack('i',hdr[4:8][::-1])[0]

    data = pd.DataFrame(columns=dtype.names, index=np.arange(length))

    if reg_ftags is True:
        for dt in dtype.names:
            fid.seek(4,1) # skip fortran tag
            data.loc[:,dt] = np.fromfile(fid, dtype=dtype[dt], count=length)
            fid.seek(4,1) # skip fortran tag

    else:
        for i in np.arange(length):
            fid.seek(4, 1)  # skip fortran tag
            for dt in dtype.names:
                data.loc[i, dt] = np.fromfile(fid, dtype=dtype[dt], count=1)[0]
            fid.seek(4, 1)  # skip fortran tag

    fid.close()

    if idx is not None:
        data.index = data.loc[:,idx].values
        data.drop(idx, axis='columns', inplace=True)

    return data

def read_tilegrids(fname=None):

    if fname is None:
        fname = find_files(paths().rc_out,'tilegrids')

    length = 2
    template = template_tilegrids()

    data = read_fortran_binary(fname, template, length=length, reg_ftags=False)

    data['gridtype'] = np.char.strip(data['gridtype'].values.astype('str'))
    data.index = ['global','domain']

    return data


def read_tilecoord(fname=None):

    if fname is None:
        fname = find_files(paths().rc_out,'tilecoord')

    hdr = 3
    template = template_tilecoord()

    return read_fortran_binary(fname, template, hdr=hdr)


class LDAS_io(object):

    def __init__(self, param):

        self.files = find_files(paths().exp_root, param)
        if self.files is None:
            print 'No files for parameter: "' + param + '".'
            return

    def get_dates(self):
        return pd.to_datetime([file[-18:-5] for file in self.files], format='%Y%m%d_%H%M')

    def read_ObsFcstAna_img(self, yr, mo, da, hr, mi):

        datestr = '%04i%02i%02i_%02i%02i' % (yr, mo, da, hr, mi)
        fname = [f for f in self.files if f.find(datestr) != -1]

        if len(fname) == 0:
            print 'No files found for: "' + datestr + '".'
            return None
        elif len(fname) > 1:
            print 'Multiple files found for: "' + datestr + '".'
        else:
            fname = fname[0]

        hdr = 11
        template = template_ObsFcstAna()
        data = read_fortran_binary(fname, template, hdr=hdr)

        names = ['obs_tilenum','obs_assim','obs_species','obs_obs','obs_obsvar','obs_fcst','obs_fcstvar','obs_ana','obs_anavar']
        data = data[names]
        data.columns = [name[4::] for name in names]

        return data

    def read_ObsFcstAna_ts(self, tilenum, species):

        names = ['assim','obs','obsvar','fcst','fcstvar','ana','anavar']
        dts = self.get_dates()

        data = pd.DataFrame(columns=names, index=dts)
        for dt in dts:
            img = self.read_ObsFcstAna_img(dt.year,dt.month,dt.day,dt.hour,dt.minute)
            tmp = img.loc[(img.tilenum == tilenum) & (img.species == species), names].values
            if len(tmp != 0):
                data.loc[dts==dt,:] = tmp

        return data



