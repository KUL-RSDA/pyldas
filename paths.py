
import os
import platform

class paths(object):
    """
    Class holding the most important LDAS path definitions

    Parameters
    ----------
    root : string
        root path to the experiment directory
    exp : string
        experiment name (appended to root path)
    domain : string
        domain name (appended to experiment path)

    Attributes
    ----------
    exp_root : string
        root path for the experiment (including experiment and domain names)
    scalefile_root : string
        root path to the scaling files (!!currently not correct for the HPC!!)
    ana : string
        Path to DA analysis output
    cat : string
        Path to catchment model output
    rc_out : string
        Path to ancillary output (grid definitions, log files etc.)
    rs : string
        Path to restart files (for continuing processing or spin-up)
    plots : string
        Path to plots

    """

    def __init__(self, root=None, exp=None, domain=None):

        if root is None:
            if platform.system() == 'Windows':
                # default path for local copies on a windows machine
                self.root = os.path.join('D:', 'data_sets', 'LDAS_runs')
            elif platform.system() == 'Darwin':
                # default path for local copies on a Mac
                self.root = os.path.join('/', 'data_sets', 'LDAS_runs')
            else:
                # default path on the HPC
                self.root = os.path.join('/', 'scratch', 'leuven', '320', 'vsc32046', 'output', 'TEST_RUNS')

        # default experiment name
        if exp is None:
            exp = 'US_M36_SMOS_noDA_unscaled'

        # default domain name
        if domain is None:
            domain = 'SMAP_EASEv2_M36_US'

        self.exp_root = os.path.join(self.root,exp,'output',domain)

        self.ana = os.path.join(self.exp_root,'ana')
        self.cat = os.path.join(self.exp_root,'cat')
        self.rc_out = os.path.join(self.exp_root,'rc_out')
        self.rs = os.path.join(self.exp_root,'rs')

        self.plots = os.path.join(self.exp_root,'plots')
