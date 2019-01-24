import sys
sys.path.append(r'/data/leuven/320/vsc32046/python')

from pyldas.interface import LDAS_io

exp='US_M36_SMOS40_DA_cal_scl_prog_std'

io = LDAS_io('ObsFcstAna', exp)
io.bin2netcdf()

io = LDAS_io('xhourly', exp)
io.bin2netcdf()

