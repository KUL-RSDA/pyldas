
import sys
import getpass

# This needs to point to the directory where pyldas is located (intended for the HPC)
uid = getpass.getuser()
sys.path.append('/data/leuven/' + uid[3:6] + '/' + uid + '/python')

from pyldas.interface import LDAS_io

# experiment name and parameters to be converted to netCDF need to be passed when calling this script!
root = sys.argv[1]
exp = sys.argv[2]
domain = sys.argv[3]
for param in sys.argv[4::]:
    io = LDAS_io(param=param, root=root, exp=exp, domain=domain)
    io.bin2netcdf()
