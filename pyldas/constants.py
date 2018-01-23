
import os
import platform

class paths(object):

    def __init__(self, root=None, exp=None, domain=None):

        if root is None:
            sys = 'win' if platform.system() == 'Windows' else 'lnx'
            if sys == 'win':
                self.root = r'C:\Users\u0116961\Documents\VSC\vsc_data_copies\scratch_TEST_RUNS'
            else:
                self.root = '/scratch/leuven/320/vsc32046/output/TEST_RUNS'

        if exp is None:
            exp = 'US_M36_EASEv2_atL_DA_7Thv_CalDnew_M'

        if domain is None:
            domain = 'SMAP_EASEv2_M36_US'

        self.exp_root = os.path.join(self.root,exp,'output',domain)

        self.ana = os.path.join(self.exp_root,'ana')
        self.cat = os.path.join(self.exp_root,'cat')
        self.rc_out = os.path.join(self.exp_root,'rc_out')
        self.rs = os.path.join(self.exp_root,'rs')

        self.plots = os.path.join(self.exp_root,'plots')
