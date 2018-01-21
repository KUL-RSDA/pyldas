
import numpy as np


def get_template(param):

    if param == 'tilegrids':
        dtype, hdr, length = template_tilegrids()

    elif param == 'tilecoord':
        dtype, hdr, length = template_tilecoord()

    elif param == 'ObsFcstAna':
        dtype, hdr, length = template_ObsFcstAna()

    elif param == 'xhourly':
        dtype, hdr, length = template_xhourly()

    else:
        print 'No template found for "' + param + '".'
        dtype, hdr, length = (None, None, None)

    return dtype, hdr, length

def template_tilegrids():

    hdr = None
    length = 2
    dtype = np.dtype([('gridtype', '|S40'),
                      ('intbase', '>i4'),
                      ('i_dir', '>i4'),
                      ('j_dir', '>i4'),
                      ('N_lon', '>i4'),
                      ('N_lat', '>i4'),
                      ('i_offg', '>i4'),
                      ('j_offg', '>i4'),
                      ('ll_lon', '>f4'),
                      ('ll_lat', '>f4'),
                      ('ur_lon', '>f4'),
                      ('ur_lat', '>f4'),
                      ('dlon', '>f4'),
                      ('dlat', '>f4')])

    return dtype, hdr, length

def template_tilecoord():

    hdr = 3
    length = None
    dtype = np.dtype([('tile_id', '>i4'),
                      ('typ', '>i4'),
                      ('pfaf', '>i4'),
                      ('com_lon', '>f4'),
                      ('com_lat', '>f4'),
                      ('min_lon', '>f4'),
                      ('max_lon', '>f4'),
                      ('min_lat', '>f4'),
                      ('max_lat', '>f4'),
                      ('i_indg', '>i4'),
                      ('j_indg', '>i4'),
                      ('frac_cell', '>f4'),
                      ('frac_pfaf', '>f4'),
                      ('area', '>f4'),
                      ('elev', '>f4')])

    return dtype, hdr, length

def template_ObsFcstAna():

    hdr = 11
    length = None
    dtype = np.dtype([('obs_assim', '>i4'),
                      ('obs_species', '>i4'),
                      ('obs_tilenum', '>i4'),
                      ('obs_lon', '>f4'),
                      ('obs_lat', '>f4'),
                      ('obs_obs', '>f4'),
                      ('obs_obsvar', '>f4'),
                      ('obs_fcst', '>f4'),
                      ('obs_fcstvar', '>f4'),
                      ('obs_ana', '>f4'),
                      ('obs_anavar', '>f4')])

    return dtype, hdr, length

def template_xhourly():

    hdr = None
    length = None
    dtype = np.dtype([('sm_surface', '>f4'),
                      ('sm_rootzone', '>f4'),
                      ('sm_profile', '>f4'),
                      ('soil_temp_layer1', '>f4'),
                      ('snow_mass', '>f4'),
                      ('precipitation_total_surface_flux', '>f4')])

    return dtype, hdr, length