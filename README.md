# **PYLDAS**

This package provides tools to convert fortran binary output of LDASsa to more convenient NetCDF output as well as vizualization functionalities for plotting time series and images, in particular catchment model parameters, RTM parameters, etc.

## Installation

Pyldas does not need to be installed as a package, so you can just clone it into your python workspace. There is one convenient batch script (see later) intended for the HPC that assumes pyldas to be at:

`/data/leuven/.../vsc...../python/pyldas/`

It is therefore recommended to clone pyldas into this directory (your VSC number is automatically recognized).

pyldas requires **python 3** and the following python packages, which can be installed using e.g. (mini)conda:

`numpy, pandas, xarray, netCDF4, matplotlib, Basemap, nco`

## Default directories

pyldas assumes a default directory as root directory where your LDASsa experiments are stored. You can, however, specify any other root directory when calling the pyldas interface (see later).

By default, LDAS experiment folders are assumed to be located at:

`/scratch/leuven/.../vsc...../output/TEST_RUNS/` (on the HPC; your VSC number will be automatically recognized)

`D:\data_sets\LDAS_runs` (on a local Windows machine)

`/data_sets/LDAS_runs` (on a local Unix machine)

I recommend to stick to these paths so that you don't always have to pass the root directory.

## Usage

Core of pyldas is the **LDAS_io** class located within pyldas/interface.py. 

### Reading model parameters

When creating an instance of the LDAS_io class, you have to **at least** specify the experiment name (corresponding to the experiment folder), which allows to directly read catchment and RTM model parameters as well as grid specific parameters.

```
from pyldas.interface import LDAS_io

io = LDAS_io(exp='US_M36_SMOS_DA', domain=None, root=None)

print(io.read_params('catparam'))
print(io.read_params('RTMparam'))

print(io.grid.tilecoord)
print(io.grid.tilegrids)
```

It is possible to pass a domain name as well (e.g. 'SMAP_EASEv2_M36_US'), but since it is not common to have more than one domain folders within one experiment folder, it will be recognized automatically if you don't. Also, as mentioned earlier, you can specify an alternative root path to your experiment directories if you don't wan't to stick with the default ones (see above). 

### Dynamic model output 

#### Initialization

To work with dynamic model output, you have to create an LDAS_io instance for each parameter you want to work with. E.g.:

```
from pyldas.interface import LDAS_io

io = LDAS_io(param='xhourly', exp='US_M36_SMOS_DA')
```

Currently supported parameters are `'xhourly', 'ObsFcstAna', 'incr', 'rstrt'`. It is generally easy to extend that list by just defining the field names of the fortran files in pyldas/templates.py. The name of the parameter just has to match the parameter name used in the fortran file names (case sensitive)! Note, however, that it is not yet possible to specify a specific collection ID! This feature should be implemented soon!

!! Note, that LDAS_io requires the existence of the tilegrids, tilecoords, and obsparam files (for DA runs), so if you work on your local machine, you have to at least copy those two/three files sticking with the folder structure on the HPC !!

#### NetCDF file creation

NetCDF files are created with the 'bin2netcdf' function, which converts the fortran image files of the specified parameter into a NetCDF **image** cube:

```
from pyldas.interface import LDAS_io

io = LDAS_io(param='xhourly', exp='US_M36_SMOS_DA')

io.bin2netcdf()
```

It is possible to specify a spatial / temporal subregion, if files are getting too large (see the documentation in pyldas/interface.py).

This routine creates a NetCDF file called 'xhourly_images.nc' in the ens_avg directory where your parameter files are located. This file is stored using **image chunks**, which means that it is very fast for reading individual spatial images, but terribly slow in extracting time series... Therefore, you should use the so-called netcdf kitchen sink to re-chunk the **image cube** into a **time series cube**. This is a command line programm, that is pre-installed on the HPC, which you will need to run there because your local machine will not have enough memory. The general terminal command (obviously to be executed from the directory where the xxx_images.nc is located) is:

```
ncks -4 -L 4 --cnk_dmn time,30000 --cnk_dmn lat,1 --cnk_dmn lon,1 xhourly_images.nc xhourly_timeseries.nc
```

On the HPC, this needs to be executed in the queue. Conveniently, there is a batch file available that automatically calls python, creates the image cubes and then rechunks them into time series cubes for as many parameters as you want. 

This batch files is located at **pyldas/batchscripts/create_timeseries.pbs** and called simply as:

```
qsub create_timeseries.pbs
```

**THIS IS THE ONLY THING YOU NEED TO DO AFTER RUNNING LDAS** which will create all the nice NetCDF files that will make your life so much easier, regardless of whether you want to use python, matlab, or whatever. Just make sure that you specify the correct working directory and all the parameters for which you want to create NetCDF files within the batch script.

There is also a batch script called **pyldas/batchscripts/copy_files.sh** for automatically copying all relevant LDAS parameter and the generated NetCDF files to your local machine, if you prefer to work locally. It is just hard-coded for my file structure, so you'd need to change that.

### Using the time series / image cubes with pyldas

When creating an instance of LDAS_io for a specific parameter, it automatically checks, if the xxx_images.nc and xxx_timeseries.nc files have been created. If so, it directly loads them and allows you to easily read images and time series based on a specify date and location:

```
from pyldas.interface import LDAS_io

io = LDAS_io(param='xhourly', exp='US_M36_SMOS_DA')

ts = read_ts(param, longitude, latitude, species=None)
img = read_image(year, month, day, hour, minute, species=None)
```

where 'param' is the parameter within the NetCDF file that should be read (e.g., catdef); and the species can be optionally specified if available.


### Vizualization

There are some routines available to automatically generate plots of e.g. all catchment parameters, all RTM parameters, etc... For details, check **pyldas/vizualize/plots.py**


