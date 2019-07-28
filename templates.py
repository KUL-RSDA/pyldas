
import logging
import numpy as np

logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

def get_template(param, **kwargs):
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

    elif param == 'ObsFcstAnaEns':
        dtype, hdr, length = template_ObsFcstAnaEns()

    elif (param == 'xhourly')|(param=='ensstd'):
        dtype, hdr, length = template_xhourly()

    elif param == 'hscale':
        dtype, hdr, length = template_scaling(**kwargs)

    elif param == 'error':
        dtype, hdr, length = template_error_Tb40()

    elif param == 'RTMparam':
        dtype, hdr, length = template_RTMparam()

    elif param == 'catparam':
        dtype, hdr, length = template_catparam()

    elif (param == 'incr')|(param == 'rstrt'):
        dtype, hdr, length = template_incr_rstrt()

    elif param == 'smosL4SMaup':
        dtype, hdr, length = template_smosL4SMaup()

    else:
        logging.warning('No template found for "' + param + '".')
        dtype, hdr, length = (None, None, None)

    return dtype, hdr, length

def template_catparam():
    """" Template for reading rtm parameter files """

    hdr = None
    length = None
    dtype = np.dtype([('dpth', '>f4'),
                      ('dzsf', '>f4'),
                      ('dzrz', '>f4'),
                      ('dzpr', '>f4'),
                      ('dzgt1', '>f4'),
                      ('dzgt2', '>f4'),
                      ('dzgt3', '>f4'),
                      ('dzgt4', '>f4'),
                      ('dzgt5', '>f4'),
                      ('dzgt6', '>f4'),
                      ('poros', '>f4'),
                      ('cond', '>f4'),
                      ('psis', '>f4'),
                      ('bee', '>f4'),
                      ('wpwet', '>f4'),
                      ('gnu', '>f4'),
                      ('vgwmax', '>f4'),
                      ('vegcls', '>i4'),
                      ('soilcls30', '>i4'),
                      ('soilcls100', '>i4'),
                      ('bf1', '>f4'),
                      ('bf2', '>f4'),
                      ('bf3', '>f4'),
                      ('cdcr1', '>f4'),
                      ('cdcr2', '>f4'),
                      ('ars1', '>f4'),
                      ('ars2', '>f4'),
                      ('ars3', '>f4'),
                      ('ara1', '>f4'),
                      ('ara2', '>f4'),
                      ('ara3', '>f4'),
                      ('ara4', '>f4'),
                      ('arw1', '>f4'),
                      ('arw2', '>f4'),
                      ('arw3', '>f4'),
                      ('arw4', '>f4'),
                      ('tsa1', '>f4'),
                      ('tsa2', '>f4'),
                      ('tsb1', '>f4'),
                      ('tsb2', '>f4'),
                      ('atau', '>f4'),
                      ('btau', '>f4'),
                      ('gravel30', '>f4'),
                      ('orgC30', '>f4'),
                      ('orgC', '>f4'),
                      ('sand30', '>f4'),
                      ('clay30', '>f4'),
                      ('sand', '>f4'),
                      ('clay', '>f4'),
                      ('wpwet30', '>f4'),
                      ('poros30', '>f4')])

    return dtype, hdr, length


def template_RTMparam():
    """" Template for reading rtm parameter files """

    hdr = 3
    length = 1
    dtype = np.dtype([('vegcls', '>i4'),
                      ('soilcls', '>i4'),
                      ('sand', '>f4'),
                      ('clay', '>f4'),
                      ('poros', '>f4'),
                      ('wang_wt', '>f4'),       # transition soil moisture
                      ('wang_wp', '>f4'),       # wilting point
                      ('rgh_hmin', '>f4'),      # soil roughness (at saturation)
                      ('rgh_hmax', '>f4'),      # soil roughness (at/below transition soil moisture)
                      ('rgh_wmin', '>f4'),
                      ('rgh_wmax', '>f4'),
                      ('rgh_Nrh', '>f4'),       # angular dependence of roughness
                      ('rgh_Nrv', '>f4'),
                      ('rgh_polmix', '>f4'),    # polarization mixing ration
                      ('omega', '>f4'),         # scattering albedo
                      ('bh', '>f4'),            # vegetation structure parameter
                      ('bv', '>f4'),
                      ('lewt', '>f4')])         # leave equivalent water thickness

    return dtype, hdr, length


def template_incr_rstrt():
    """" Template for reading increment or restart files """

    hdr = None
    length = None
    dtype = np.dtype([('tc1', '>f4'),
                      ('tc2', '>f4'),
                      ('tc4', '>f4'),
                      ('qa1', '>f4'),
                      ('qa2', '>f4'),
                      ('qa4', '>f4'),
                      ('capac', '>f4'),
                      ('catdef', '>f4'),
                      ('rzexc', '>f4'),
                      ('srfexc', '>f4'),
                      ('ght1', '>f4'),
                      ('ght2', '>f4'),
                      ('ght3', '>f4'),
                      ('ght4', '>f4'),
                      ('ght5', '>f4'),
                      ('ght6', '>f4'),
                      ('wesn1', '>f4'),
                      ('wesn2', '>f4'),
                      ('wesn3', '>f4'),
                      ('htsn1', '>f4'),
                      ('htsn2', '>f4'),
                      ('htsn3', '>f4'),
                      ('sndz1', '>f4'),
                      ('sndz2', '>f4'),
                      ('sndz3', '>f4')])

    return dtype, hdr, length

def template_error_Tb40():
    """ Template for reading uncertainty files. """

    # 23 header fields + 1 incidence angles (normalized to 40 degrees)
    # TODO: allow for a different number of inc. angles when using for SMAP (# angles on hdr pos 20)
    hdr = 25
    length = 19

    dtype = np.dtype([('lon', '>f4'),
                      ('lat', '>f4'),
                      ('err_Tbh', '>f4'),
                      ('err_Tbv', '>f4')])

    return dtype, hdr, length


def template_scaling(sensor='SMAP'):
    """ Template for reading scaling files. """

    # 23 header fields + 7 incidence angles
    # TODO: allow for a different number of inc. angles when using for SMAP (# angles on hdr pos 20)

    if sensor == 'SMOS':
        hdr = 32
        length = 19
        angles = [30,35,40,45,50,55,60]
    else:
        hdr = 26
        length = 19
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

def template_ObsFcstAnaEns():
    """" Template for reading innovation files """

    hdr = 11
    length = None
    dtype = np.dtype([('obs_ensmem', '>i4'),
                      ('obs_species', '>i4'),
                      ('obs_tilenum', '>i4'),
                      ('obs_lon', '>f4'),
                      ('obs_lat', '>f4'),
                      ('obs_obs', '>f4'),
                      ('obs_fcst', '>f4'),
                      ('obs_ana', '>f4')])

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

def template_smosL4SMaup():
    """"
    Template for reading smos L4 soil moisture DA output

    """

    hdr = None
    length = None
    dtype = np.dtype([('srfexc_fcst', '>f4'),
                      ('rzexc_fcst', '>f4'),
                      ('catdef_fcst', '>f4'),
                      ('srfexc_ana', '>f4'),
                      ('rzexc_ana', '>f4'),
                      ('catdef_ana', '>f4')])

    return dtype, hdr, length









