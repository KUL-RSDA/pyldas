
import numpy as np

def template_tilegrids():

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

    return dtype

def template_tilecoord():
    
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

    return dtype

def template_ObsFcstAna():

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

    return dtype