# **PYLDAS**

## Python package for converting and vizualizing LDASsa output.

This package is intended to convert the fortran binary output of LDSAsa to more convenient NetCDF output, and to provide simple vizualization functionalities for plotting time series and images, in particular catchment model parameters, RTM parameters, etc.

### Installation

Currently, pyldas cannot be installed as a package, so it just has to be cloned to the python workspace. There is one convenient batch script (see later) that assumes the python workspace to be at:

`/data/leuven/.../vsc...../python`

It is therefore recommended to clone pyldas into this directory (your VSC number is automatically recognized).

pyldas requires **python 3** and the following standard python packages, which should be installed using (mini)conda:

`numpy, pandas, xarray, netCDF4, matplotlib, Basemap`

### Usage

The basic idea is that pyldas is used to convert LDAS output, i.e. fortran binaries, to more easy-to-work-with NetCDF files. To do that, you don't need a python IDE... There is a batch script to do just that (see below), so if you prefer to work with matlab, cdos, etc., you can just use pyldas from the command line to reformat files... However, pyldas also provides some more vizualization functionality (e.g. plotting time series based on lat/lon, plotting catchment model parameters, etc.), which you can use locally or on the HPC using any python IDE (I recommend pycharm).

#### Default directories

pyldas assumes a few default directories. If you stick to these, pyldas works right away and you don't need to modify anything!

LDAS runs are assumed to be located at:
`/data/leuven/.../vsc...../python`
