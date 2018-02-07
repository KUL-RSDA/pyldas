
import numpy as np


def get_template(param):
    """
    Return templates for reading fortran binary files

    Parameters
    ----------
    param : string
        Name of the parameter for which to return the reading template

    Returns
    -------
    dtype : np.dtype
        Field names and formats of the entries in the binary files
    hdr : int
        Number of (4-byte) header entries to SKIP before the data block
    length : int
        Number of data entries to be read
        If None, this will be inferred from the second header byte

    """

    if param == 'tilegrids':
        dtype, hdr, length = template_tilegrids()

    elif param == 'tilecoord':
        dtype, hdr, length = template_tilecoord()

    elif param == 'ObsFcstAna':
        dtype, hdr, length = template_ObsFcstAna()

    elif param == 'xhourly':
        dtype, hdr, length = template_xhourly()

    elif param == 'scaling':
        dtype, hdr, length = template_scaling()

    else:
        print 'No template found for "' + param + '".'
        dtype, hdr, length = (None, None, None)

    return dtype, hdr, length

def template_scaling(sensor='SMOS'):
    """ Template for reading scaling files. """

    # 23 header fields + 7 incidence angles
    # TODO: allow for a different number of inc. angles when using for SMAP (# angles on hdr pos 20)
    hdr = 32
    length = 19

    if sensor == 'SMOS':
        angles = [30,35,40,45,50,55,60]
    else:
        angles = [40,]

    dtype = np.dtype([('lon', '>f4'),('lat', '>f4'),('tile_id', '>i4')]+
                     [('m_obs_H_%i'%ang, '>f4') for ang in angles] +
                     [('s_obs_H_%i'%ang, '>f4') for ang in angles] +
                     [('m_mod_H_%i'%ang, '>f4') for ang in angles] +
                     [('s_mod_H_%i'%ang, '>f4') for ang in angles] +
                     [('N_data_H_%i'%ang, '>i4') for ang in angles] +
                     [('m_obs_V_%i'%ang, '>f4') for ang in angles] +
                     [('s_obs_V_%i'%ang, '>f4') for ang in angles] +
                     [('m_mod_V_%i'%ang, '>f4') for ang in angles] +
                     [('s_mod_V_%i'%ang, '>f4') for ang in angles] +
                     [('N_data_V_%i'%ang, '>i4') for ang in angles])

    return dtype, hdr, length


def template_tilegrids():
    """" Template for reading the 'tilegrids' binary file """

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
    """" Template for reading the 'tilecoord' binary file """

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
    """" Template for reading innovation files """

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
    """"
    Template for reading xhourly catchment output files
    TODO Include the possibility to specify the Collection ID. Currently: 8 / tavg

    """


    hdr = None
    length = None
    dtype = np.dtype([('sm_surface', '>f4'),
                      ('sm_rootzone', '>f4'),
                      ('sm_profile', '>f4'),
                      ('soil_temp_layer1', '>f4'),
                      ('snow_mass', '>f4'),
                      ('precipitation_total_surface_flux', '>f4')])

    return dtype, hdr, length