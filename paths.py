
import logging
import getpass
import platform
from pathlib import Path

class paths(object):
    """
    Class holding the most important LDAS path definitions

    Parameters
    ----------
    mode : string
        'LDASsa' or "GEOSldas'
    root : pathlib.Path
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

    def __init__(self, mode, root=None, exp=None, domain=None):

        if root is None:
            uid = getpass.getuser()
            if uid[:3] == 'vsc':
                # default path on the HPC (Tier-1 and Tier-2)
                self.root = Path('/') / 'scratch' / 'leuven' / uid[3:6] / uid / mode
                if not self.root.exists():
                    self.root = Path('/') / 'scratch' / 'leuven' / 'projects' / 'lt1_2020_es_pilot' / 'project_output' / 'rsda' / uid / mode
            else:
                # default path on local machines
                self.root = Path('~').expanduser() / 'data_sets' / f'{mode}_runs'
        else:
            self.root = Path(root)

        # default experiment name
        if exp is None:
            if mode == 'GEOSldas':
                exp = 'NLv4_M36_US_OL_Pcorr'
            else:
                exp = 'US_M36_SMAP_TB_OL_scaled_4K_obserr'

        # default domain name
        if domain is None:
            domains = list(self.root.glob(exp+'/output/*'))
            domains = [d for d in domains if d.name[0] != '.'] # remove hidden folders/files (dot-files)
            if len(domains) != 1:
                logging.error('Domain could not be identified.')
                return
            else:
                domain = domains[0].name

        self.exp_root = self.root / exp / 'output' / domain

        self.ana = self.exp_root / 'ana'
        self.cat = self.exp_root / 'cat'
        self.rc_out = self.exp_root / 'rc_out'
        self.rs = self.exp_root / 'rs'

        # Probably not a good location, but OK for a default option
        self.plots = self.exp_root / 'plots'
