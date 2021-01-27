########################################
####     Behzad Asadieh, Ph.D.      ####
####  University of Pennsylvania    ####
####    basadieh@sas.upenn.edu      ####
########################################

### IMPORTANT: The following path needs to be appended so that the neccessary codes can be read from that directory ###
import sys
sys.path.append('/data1/home/basadieh/behzadcodes/behzadlibrary')
from Behzadlib import * ### This is loaded from '/data1/home/basadieh/behzadcodes/behzadlibrary'
from CMIP5lib import *  ### This is loaded from '/data1/home/basadieh/behzadcodes/behzadlibrary'
########################################
import numpy as np
from numpy import zeros, ones, empty, nan, shape
from numpy import isnan, nanmean, nanmax, nanmin
import numpy.ma as ma
from netCDF4 import MFDataset, Dataset, num2date, date2num, date2index
import os
import matplotlib
import matplotlib.mlab as ml
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm, maskoceans
from scipy.interpolate import griddata
import math
import copy
########################################
########################################
direct= '/data2/scratch/cabre/CMIP5/CMIP5_models/ocean_physics/GFDL-ESM2G/historical/mo/'
file_name='thetao_Omon_GFDL-ESM2G_historical_r1i1p1_200101-200512.nc'

file_dire=direct+file_name

dset_t = Dataset(file_dire)
dset_t.variables

variable='thetao'

Lat_orig=np.asarray(dset_t.variables['lat'][:])
Lon_orig=np.asarray(dset_t.variables['lon'][:])

Var_main = np.asarray(dset_t.variables[variable][:,0,:,:]) # Reading variable only at the surface
Var_main_ave = nanmean(Var_main, 0)
Var_main_ave [ Var_main_ave > 1e19 ] = nan

### Simple ploting ###
plt.imshow(np.flipud(Var_main_ave), cmap=plt.cm.jet)
plt.colorbar()
plt.show()

### Basemap Ploting ###
bounds_max=float("{0:.0f}".format(np.nanpercentile(Var_main_ave , 99.99)))
func_plot(Var_main_ave, Lat_orig, Lon_orig, bounds_max, '-', '-', 'cyl', 0)


#### Regridingthe variable ###
lat_n_regrid, lon_n_regrid = 180, 360 # Number of Lat and Lon elements in the regridded data
lon_min_regrid, lon_max_regrid = 0, 360 # Min and Max value of Lon in the regridded data
lat_min_regrid, lat_max_regrid = -90, 90 # Min and Max value of Lat in the regridded data
Lat_regrid_1D, Lon_regrid_1D, Lat_bound_regrid, Lon_bound_regrid = func_latlon_regrid(lat_n_regrid, lon_n_regrid, lat_min_regrid, lat_max_regrid, lon_min_regrid, lon_max_regrid)
Lon_regrid_2D, Lat_regrid_2D = np.meshgrid(Lon_regrid_1D, Lat_regrid_1D)

Var_main_ave_regrid = func_regrid(Var_main_ave, Lat_orig, Lon_orig, Lat_regrid_2D, Lon_regrid_2D)

### Regrided variable ploting ###
bounds_max=float("{0:.0f}".format(np.nanpercentile(Var_main_ave_regrid , 99.99)))
func_plot(Var_main_ave_regrid, Lat_regrid_2D, Lon_regrid_2D, bounds_max, '-', '-', 'cyl', 0)






#####
aaa=[]
for ii in range(1,stream_AMOC_max.shape[0]):

aaa=(stream_AMOC_max-np.nanmean(stream_AMOC_max))/np.std(stream_AMOC_max)


plt.imshow( np.transpose( aaa), cmap=plt.cm.jet)
plt.colorbar()
plt.show()






















