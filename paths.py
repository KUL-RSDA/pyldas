
import logging
import getpass
import platform
from pathlib import Path

class paths(object):
    """
    Class holding the most important LDAS path definitions

    Parameters
    ----------
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

    def __init__(self, root=None, exp=None, domain=None):

        if root is None:
            root = 'D:' if platform.system() == 'Windows' else '/'
            uid = getpass.getuser()
            if uid[:3] == 'vsc':
                # default path on the HPC
                self.root = Path(root) / 'scratch' / 'leuven' / uid[3:6] / uid / 'output' / 'TEST_RUNS'
            else:
                # default path on local machines
                self.root = Path(root) / 'data_sets' / 'LDAS_runs'
        else:
            self.root = Path(root)

        # default experiment name
        if exp is None:
            exp = 'US_M36_SMOS_noDA_unscaled'

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
