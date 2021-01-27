#################################
####   CMIP5 - CO2 & Phyto   ####
########################################
####     Behzad Asadieh, Ph.D.      ####
####  University of Pennsylvania    ####
####    basadieh@sas.upenn.edu      ####
########################################

### IMPORTANT: The following path needs to be appended so that the neccessary codes can be read from that directory ###
import sys
sys.path.append('/data1/home/basadieh/behzadcodes/behzadlibrary')
from Behzadlib import func_latlon_regrid, func_latlon_regrid_eq, func_oceanlandmask, func_oceanindex, func_regrid, func_plot, func_plot_save, func_plot_bounds_save
from BehzadlibPlot import func_plot2Dcolor_contourf, func_plotline_1var, func_plotline_2var, func_plot2Dcontour_1var, func_plot2Dcontour_2var
from CMIP5lib import netcdf_read
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

dir_pwd = os.getcwd() # Gets the current directory (and in which the code is placed)
dir_data_in1 = ('/data2/scratch/cabre/CMIP5/CMIP5_models/ocean_biogeochemistry/') # Directory to raed raw data from
dir_data_in2 = ('/data2/scratch/cabre/CMIP5/CMIP5_models/ocean_physics/')
dir_data_in4 = ('/data2/scratch/cabre/CMIP5/CMIP5_models/Heat/')
dir_data_in3 = ('/data5/scratch/CMIP5/CMIP5_models/ocean_biogeochemistry/')
dir_figs = (dir_pwd + '/Figures/') # Directory to save processed data
dir_results = (dir_pwd + '/Results/') # Directory to save results

### Regrdridding calculations ###
# creating new coordinate grid, same which was used in interpolation in data processing code
lat_n_regrid, lon_n_regrid = 180, 360 # Number of Lat and Lon elements in the regridded data
lon_min_regrid, lon_max_regrid = 0, 360 # Min and Max value of Lon in the regridded data
lat_min_regrid, lat_max_regrid = -90, 90 # Min and Max value of Lat in the regridded data

Lat_regrid_1D, Lon_regrid_1D, Lat_bound_regrid, Lon_bound_regrid = func_latlon_regrid(lat_n_regrid, lon_n_regrid, lat_min_regrid, lat_max_regrid, lon_min_regrid, lon_max_regrid)
Lon_regrid_2D, Lat_regrid_2D = np.meshgrid(Lon_regrid_1D, Lat_regrid_1D)

Ocean_Land_mask = func_oceanlandmask(Lat_regrid_2D, Lon_regrid_2D) # 1= ocean, 0= land

Ocean_Index = func_oceanindex (Lat_regrid_2D, Lon_regrid_2D) # [0=land] [2=Pacific] [3=Indian Ocean] [6=Atlantic] [10=Arctic] [8=Baffin Bay (west of Greenland)] [9=Norwegian Sea (east of Greenland)] [11=Hudson Bay (Canada)] 
###############################################################################
#### Anna's files
#### Phydiat (P2) reading ###
#dir_data_in1 = ('/data5/scratch/forTiho/') # Directory to raed raw data from
#Var_file_name='diatom' # The variable name to be read from .nc files
#dset1 = Dataset(dir_data_in1+Var_file_name+'-monthly-1990-2010_'+str(GCM)+'.nc') 
##dset1.variables
#Var_name='DIATOM'
#
#Lat_orig=np.asarray(dset1.variables['LAT'][:])
#Lon_orig=np.asarray(dset1.variables['LON'][:])
#
#Phydiat_Z0_Allmonths = np.asarray(dset1.variables[Var_name][:])
#dset1.close() # closing netcdf file
###############################################################################

GCM_Names = ['GFDL-ESM2G']
GCM=GCM_Names[0]

dir_data_in_h=(dir_data_in1+ GCM + '/historical/mo/') # Directory to read Historical data
dir_data_in_r=(dir_data_in1+ GCM + '/rcp85/mo/') # Directory to read RCP8.5 data

### Variables to be edited by user ###
start_date_plt=1990 # Start Date used for plotting
end_date_plt=2010 # End Date used for plotting

start_date_cal_h=1990 # Start Date for Calculations # Define as string, not number
end_date_cal_h=2005 # End Date for Calculations
start_date_cal_r=2006 # Start Date for Calculations # Define as string, not number
end_date_cal_r=2010 # End Date for Calculations

yrs_n=21;

############################
### Phydiat (P2) reading ###
variable='phydiat'
dset1= netcdf_read(dir_data_in_h+variable+'_Omon_'+str(GCM)+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
## dset1.variables  # Shows the variables in the .nc file
## This finds which row in the time variable corresponds to the desired start_date_cal and end_date_cal for this specific calculation
start_date_i,end_date_i = dset1.find_time(dset1.times, start_date_cal_h, end_date_cal_h)
Lat_orig=dset1.y
Lon_orig=dset1.x
dset1.close_ncfile(dset1.fin)

dset2 = MFDataset(dir_data_in_h+variable+'_Omon_'+str(GCM)+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read

Var_dset = dset2.variables[variable] # append phydiat data to variable

Data_allmonths_orig=np.asarray(Var_dset[start_date_i:end_date_i+1,:,:]) # phydiat(time, rlat, rlon) - phydiat of all the desired months at the surface
Data_allmonths_orig[Data_allmonths_orig <= 0] = nan
#### Adding RCP data to the Historical data ####
dset1= netcdf_read(dir_data_in_r+variable+'_Omon_'+str(GCM)+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
start_date_i,end_date_i = dset1.find_time(dset1.times, start_date_cal_r, end_date_cal_r)
dset1.close_ncfile(dset1.fin)
dset2 = MFDataset(dir_data_in_r+variable+'_Omon_'+str(GCM)+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
Var_dset = dset2.variables[variable]
Data_allmonths_orig_r=np.asarray(Var_dset[start_date_i:end_date_i+1,:,:]) # phydiat(time, rlat, rlon) - phydiat of all the desired months at the surface
Data_allmonths_orig_r[Data_allmonths_orig <= 0] = nan
Data_allmonths_orig = np.concatenate((Data_allmonths_orig, Data_allmonths_orig_r), axis=0)
### Regriding data to 180*360 degrees
Phydiat_allmonths = func_regrid(Data_allmonths_orig, Lat_orig, Lon_orig, Lat_regrid_2D, Lon_regrid_2D) # Phydiat at surface - all months 
Phydiat_allmonths[Phydiat_allmonths > 1E19] = nan # Replacing missing values with NaN - Missing values in the original .nc file are filled with 1E+20 values

############################
### Phyc (P1+P2) reading ###
variable='phyc'
dset1= netcdf_read(dir_data_in_h+variable+'_Omon_'+str(GCM)+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read

start_date_i,end_date_i = dset1.find_time(dset1.times, start_date_cal_h, end_date_cal_h) # This finds which row in the time variable corresponds to the desired start_date_cal and end_date_cal for this specific calculation
dset1.close_ncfile(dset1.fin)

dset2 = MFDataset(dir_data_in_h+variable+'_Omon_'+str(GCM)+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
Var_dset = dset2.variables[variable] # append phyc data to variable

Data_allmonths_orig=np.asarray(Var_dset[start_date_i:end_date_i+1,:,:]) # phyc(time, rlat, rlon) - phydiat of all the desired months at the surface
Data_allmonths_orig[Data_allmonths_orig <= 0] = nan
#### Adding RCP data to the Historical data ####
dset1= netcdf_read(dir_data_in_r+variable+'_Omon_'+str(GCM)+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
start_date_i,end_date_i = dset1.find_time(dset1.times, start_date_cal_r, end_date_cal_r)
dset1.close_ncfile(dset1.fin)
dset2 = MFDataset(dir_data_in_r+variable+'_Omon_'+str(GCM)+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
Var_dset = dset2.variables[variable]
Data_allmonths_orig_r=np.asarray(Var_dset[start_date_i:end_date_i+1,:,:]) # phyc(time, rlat, rlon) - phydiat of all the desired months at the surface
Data_allmonths_orig_r[Data_allmonths_orig <= 0] = nan
Data_allmonths_orig = np.concatenate((Data_allmonths_orig, Data_allmonths_orig_r), axis=0)
### Regriding data to 180*360 degrees
Phyc_allmonths = func_regrid(Data_allmonths_orig, Lat_orig, Lon_orig, Lat_regrid_2D, Lon_regrid_2D) # Phyc at surface - all months 
Phyc_allmonths[Phyc_allmonths > 1E19] = nan # Replacing missing values with NaN - Missing values in the original .nc file are filled with 1E+20 values

############################
### Chlorophyll reading ###
variable='chl'
dset1= netcdf_read(dir_data_in_h+variable+'_Omon_'+str(GCM)+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read

start_date_i,end_date_i = dset1.find_time(dset1.times, start_date_cal_h, end_date_cal_h) # This finds which row in the time variable corresponds to the desired start_date_cal and end_date_cal for this specific calculation
dset1.close_ncfile(dset1.fin)

dset2 = MFDataset(dir_data_in_h+variable+'_Omon_'+str(GCM)+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
Var_dset = dset2.variables[variable] # append phyc data to variable

Data_allmonths_orig=np.asarray(Var_dset[start_date_i:end_date_i+1,:,:]) # chl(time, rlat, rlon) - phydiat of all the desired months at the surface
Data_allmonths_orig[Data_allmonths_orig <= 0] = nan
#### Adding RCP data to the Historical data ####
dset1= netcdf_read(dir_data_in_r+variable+'_Omon_'+str(GCM)+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
start_date_i,end_date_i = dset1.find_time(dset1.times, start_date_cal_r, end_date_cal_r)
dset1.close_ncfile(dset1.fin)
dset2 = MFDataset(dir_data_in_r+variable+'_Omon_'+str(GCM)+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
Var_dset = dset2.variables[variable]
Data_allmonths_orig_r=np.asarray(Var_dset[start_date_i:end_date_i+1,:,:]) # chl(time, rlat, rlon) - phydiat of all the desired months at the surface
Data_allmonths_orig_r[Data_allmonths_orig <= 0] = nan
Data_allmonths_orig = np.concatenate((Data_allmonths_orig, Data_allmonths_orig_r), axis=0)
### Regriding data to 180*360 degrees
Chl_allmonths = func_regrid(Data_allmonths_orig, Lat_orig, Lon_orig, Lat_regrid_2D, Lon_regrid_2D) # Chl at surface - all months 
Chl_allmonths[Chl_allmonths > 1E19] = nan # Replacing missing values with NaN - Missing values in the original .nc file are filled with 1E+20 values

########################################
### Phydiat (Large Phyto) percentage ###
Phydiat_prcnt_allmonths = (Phydiat_allmonths / Phyc_allmonths) * 100
### Small phyto ######
Physmall_allmonths = Phyc_allmonths - Phydiat_allmonths

#############################################
### Sea Surface Temperature (SST) reading ###
variable='thetao'

dir_data_in_h_t=(dir_data_in2+ GCM + '/historical/mo/') # Directory to read Historical data
dir_data_in_r_t=(dir_data_in2+ GCM + '/rcp85/mo/') # Directory to read RCP8.5 data

dset1= netcdf_read(dir_data_in_h_t+variable+'_Omon_'+str(GCM)+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read

start_date_i,end_date_i = dset1.find_time(dset1.times, start_date_cal_h, end_date_cal_h) # This finds which row in the time variable corresponds to the desired start_date_cal and end_date_cal for this specific calculation
dset1.close_ncfile(dset1.fin)

dset2 = MFDataset(dir_data_in_h_t+variable+'_Omon_'+str(GCM)+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
Var_dset = dset2.variables[variable] # append phyc data to variable

Data_allmonths_orig=np.asarray(Var_dset[start_date_i:end_date_i+1,0,:,:]) # thetao(time, lev, rlat, rlon) - phydiat of all the desired years, only on the surface
Data_allmonths_orig[Data_allmonths_orig <= 0] = nan
#### Adding RCP data to the Historical data ####
dset1= netcdf_read(dir_data_in_r_t+variable+'_Omon_'+str(GCM)+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
start_date_i,end_date_i = dset1.find_time(dset1.times, start_date_cal_r, end_date_cal_r)
dset1.close_ncfile(dset1.fin)
dset2 = MFDataset(dir_data_in_r_t+variable+'_Omon_'+str(GCM)+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
Var_dset = dset2.variables[variable]
Data_allmonths_orig_r=np.asarray(Var_dset[start_date_i:end_date_i+1,0,:,:]) # thetao(time, lev, rlat, rlon) - phydiat of all the desired years, only on the surface
Data_allmonths_orig_r[Data_allmonths_orig <= 0] = nan
Data_allmonths_orig = np.concatenate((Data_allmonths_orig, Data_allmonths_orig_r), axis=0)
### Regriding data to 180*360 degrees
SST_allmonths = func_regrid(Data_allmonths_orig, Lat_orig, Lon_orig, Lat_regrid_2D, Lon_regrid_2D) # thetao at surface - all months 
SST_allmonths[SST_allmonths > 1E19] = nan # Replacing missing values with NaN - Missing values in the original .nc file are filled with 1E+20 values
SST_allmonths = SST_allmonths - 273.15 # To convert Temperature unit from K to C

############################
###     SpCO2 reading    ###
variable='spco2'
dset1= netcdf_read(dir_data_in_h+variable+'_Omon_'+str(GCM)+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read

start_date_i,end_date_i = dset1.find_time(dset1.times, start_date_cal_h, end_date_cal_h) # This finds which row in the time variable corresponds to the desired start_date_cal and end_date_cal for this specific calculation
dset1.close_ncfile(dset1.fin)

dset2 = MFDataset(dir_data_in_h+variable+'_Omon_'+str(GCM)+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
Var_dset = dset2.variables[variable] # append phyc data to variable

Data_allmonths_orig=np.asarray(Var_dset[start_date_i:end_date_i+1,:,:]) # spco2(time, rlat, rlon) - phydiat of all the desired months at the surface
Data_allmonths_orig[Data_allmonths_orig <= 0] = nan
#### Adding RCP data to the Historical data ####
dir_data_in_r_s=(dir_data_in3+ GCM + '/rcp85/mo/') # Directory to read RCP8.5 data
dset1= netcdf_read(dir_data_in_r_s+variable+'_Omon_'+str(GCM)+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
start_date_i,end_date_i = dset1.find_time(dset1.times, start_date_cal_r, end_date_cal_r)
dset1.close_ncfile(dset1.fin)
dset2 = MFDataset(dir_data_in_r_s+variable+'_Omon_'+str(GCM)+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
Var_dset = dset2.variables[variable]
Data_allmonths_orig_r=np.asarray(Var_dset[start_date_i:end_date_i+1,:,:]) # spco2(time, rlat, rlon) - phydiat of all the desired months at the surface
Data_allmonths_orig_r[Data_allmonths_orig <= 0] = nan
Data_allmonths_orig = np.concatenate((Data_allmonths_orig, Data_allmonths_orig_r), axis=0)
### Regriding data to 180*360 degrees
SpCO2_allmonths = func_regrid(Data_allmonths_orig, Lat_orig, Lon_orig, Lat_regrid_2D, Lon_regrid_2D) # SpCO2 at surface - all months 
SpCO2_allmonths[SpCO2_allmonths > 1E19] = nan # Replacing missing values with NaN - Missing values in the original .nc file are filled with 1E+20 values
SpCO2_allmonths = SpCO2_allmonths / 0.101325 # To Convert spco2 unit from Pa to PPM

################################
###     Heat Flux reading    ###
variable='hfds' # Downward Heat Flux at Sea Water Surface - units: W m-2

dir_data_in_h_t=(dir_data_in4+ GCM + '/historical/mo/') # Directory to read Historical data
dir_data_in_r_t=(dir_data_in4+ GCM + '/rcp85/mo/') # Directory to read RCP8.5 data

dset1= netcdf_read(dir_data_in_h_t+variable+'_Omon_'+str(GCM)+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read

start_date_i,end_date_i = dset1.find_time(dset1.times, start_date_cal_h, end_date_cal_h) # This finds which row in the time variable corresponds to the desired start_date_cal and end_date_cal for this specific calculation
dset1.close_ncfile(dset1.fin)

dset2 = MFDataset(dir_data_in_h_t+variable+'_Omon_'+str(GCM)+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
Var_dset = dset2.variables[variable] # append phyc data to variable

Data_allmonths_orig=np.asarray(Var_dset[start_date_i:end_date_i+1,:,:]) # hfds(time, rlat, rlon) - phydiat of all the desired years, only on the surface
Data_allmonths_orig[Data_allmonths_orig <= 0] = nan
#### Adding RCP data to the Historical data ####
dset1= netcdf_read(dir_data_in_r_t+variable+'_Omon_'+str(GCM)+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
start_date_i,end_date_i = dset1.find_time(dset1.times, start_date_cal_r, end_date_cal_r)
dset1.close_ncfile(dset1.fin)
dset2 = MFDataset(dir_data_in_r_t+variable+'_Omon_'+str(GCM)+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
Var_dset = dset2.variables[variable]
Data_allmonths_orig_r=np.asarray(Var_dset[start_date_i:end_date_i+1,:,:]) # hfds(time, rlat, rlon) - phydiat of all the desired years, only on the surface
Data_allmonths_orig_r[Data_allmonths_orig <= 0] = nan
Data_allmonths_orig = np.concatenate((Data_allmonths_orig, Data_allmonths_orig_r), axis=0)
### Regriding data to 180*360 degrees
HFDS_allmonths = func_regrid(Data_allmonths_orig, Lat_orig, Lon_orig, Lat_regrid_2D, Lon_regrid_2D) # HFDS at surface - all months 
HFDS_allmonths[HFDS_allmonths > 1E19] = nan # Replacing missing values with NaN - Missing values in the original .nc file are filled with 1E+20 values


############################
###     MLD reading     ###
#from Behzadlib2 import func_MLD_Allmonths_bio
#
#variable_thetao='thetao' # Sea Water Temperature
#variable_so='so' # Sea Water Salinity
#
#dir_data_in_h_t=(dir_data_in2+ GCM + '/historical/mo/') # Directory to read Historical data
#dir_data_in_r_t=(dir_data_in2+ GCM + '/rcp85/mo/') # Directory to read RCP8.5 data
#
#dset_thetao= netcdf_read(dir_data_in_h_t+variable_thetao+'_Omon_'+str(GCM)+'*12.nc',str(GCM), variable_thetao)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
#dset_so= netcdf_read(dir_data_in_h_t+variable_so+'_Omon_'+str(GCM)+'*12.nc',str(GCM), variable_so)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
#
#Data_allmonths_orig = func_MLD_Allmonths_bio(dset_thetao, dset_so, start_date_cal_h, end_date_cal_h, lat_n_regrid, lon_n_regrid)
#Data_allmonths_orig[Data_allmonths_orig <= 0] = nan
##### Adding RCP data to the Historical data ####
#dset_thetao= netcdf_read(dir_data_in_r_t+variable_thetao+'_Omon_'+str(GCM)+'*12.nc',str(GCM), variable_thetao)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
#dset_so= netcdf_read(dir_data_in_r_t+variable_so+'_Omon_'+str(GCM)+'*12.nc',str(GCM), variable_so)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
#
#Data_allmonths_orig_r = func_MLD_Allmonths_bio(dset_thetao, dset_so, start_date_cal_r, end_date_cal_r, lat_n_regrid, lon_n_regrid)
#Data_allmonths_orig_r[Data_allmonths_orig_r <= 0] = nan
#Data_allmonths_orig = np.concatenate((Data_allmonths_orig, Data_allmonths_orig_r), axis=0)
#
#MLD_allmonths = copy.deepcopy(Data_allmonths_orig) # MLD - all months 
#MLD_allmonths[MLD_allmonths > 1E19] = nan # Replacing missing values with NaN - Missing values in the original .nc file are filled with 1E+20 values
##np.save('MLD_allmonths.npy',MLD_allmonths)


dir_data_in_m=('/data4/CMIP5/CMIP5_models/descriptors/ocean_physics/') # Directory to read Historical data
variable='MLD'
dset1 = Dataset(dir_data_in_m+'mld_Omon_'+str(GCM)+'_historicalrcp85_r1i1p1_186501-210012.nc') 

Var_dset = dset1.variables[variable] # append phyc data to variable
Data_allmonths_orig=np.asarray(dset1.variables[variable][(start_date_plt-1865)*12:(end_date_plt-1865+1)*12,:,:])

MLD_allmonths  = func_regrid(Data_allmonths_orig, Lat_orig, Lon_orig, Lat_regrid_2D, Lon_regrid_2D) # SpCO2 at surface - all months 
MLD_allmonths[MLD_allmonths > 1E19] = nan # Replacing missing values with NaN - Missing values in the original .nc file are filled with 1E+20 values
MLD_allmonths[MLD_allmonths < -1E19] = nan # Replacing missing values with NaN - Missing values in the original .nc file are filled with 1E+20 values

########################################################
########################################################

SpCO2_monthlymean = SpCO2_allmonths.reshape(( yrs_n ,12,lat_n_regrid,lon_n_regrid)) # Monthly means of SpCO2 over the time period
SpCO2_monthlymean = np.nanmean(SpCO2_monthlymean,axis=0)
SpCO2_monthlymean=np.squeeze(SpCO2_monthlymean)

SST_monthlymean = SST_allmonths.reshape(( yrs_n ,12,lat_n_regrid,lon_n_regrid)) # Monthly means of SST over the time period
SST_monthlymean = np.nanmean(SST_monthlymean,axis=0)
SST_monthlymean=np.squeeze(SST_monthlymean)


SpCO2_allmonths_TEMP_eq1 = empty((SpCO2_allmonths.shape[0], lat_n_regrid,lon_n_regrid)) * nan #SpCO2_allmonths_TEMP_eq1 = np.zeros([SpCO2_allmonths.shape[0], lat_n_regrid,lon_n_regrid])
SpCO2_allmonths_NonTEMP_eq2 = empty((SpCO2_allmonths.shape[0], lat_n_regrid,lon_n_regrid)) * nan #SpCO2_allmonths_NonTEMP_eq2 = np.zeros([SpCO2_allmonths.shape[0], lat_n_regrid,lon_n_regrid])
for ii in range(SpCO2_allmonths.shape[0]):
    SpCO2_allmonths_TEMP_eq1[ii,:,:] = np.nanmean(SpCO2_allmonths,axis=0) * ( np.exp( 0.0423 * ( SST_allmonths[ii,:,:] - np.nanmean(SST_allmonths,axis=0) ) ) );
    SpCO2_allmonths_NonTEMP_eq2[ii,:,:] = SpCO2_allmonths[ii,:,:] * ( np.exp( 0.0423 * ( np.nanmean(SST_allmonths,axis=0) - SST_allmonths[ii,:,:] ) ) );


SpCO2_monthlymean_TEMP_eq1 = empty((12, lat_n_regrid,lon_n_regrid)) * nan
SpCO2_monthlymean_NonTEMP_eq2 = empty((12, lat_n_regrid,lon_n_regrid)) * nan
for ii in range(12):
    SpCO2_monthlymean_TEMP_eq1[ii,:,:] = np.nanmean(SpCO2_allmonths,axis=0) * ( np.exp( 0.0423 * ( SST_monthlymean[ii,:,:] - np.nanmean(SST_allmonths,axis=0) ) ) );
    SpCO2_monthlymean_NonTEMP_eq2[ii,:,:] = SpCO2_monthlymean[ii,:,:] * ( np.exp( 0.0423 * ( np.nanmean(SST_allmonths,axis=0) - SST_monthlymean[ii,:,:] ) ) );    


SpCO2_allmonths_NonTEMP_subtracted = SpCO2_allmonths - SpCO2_allmonths_TEMP_eq1;

###############################################################################
##############          PAPA and NABE Calculations            #################

SST_monthlymean_PAPA = SST_monthlymean[:,135:140,210:220]
SST_monthlymean_PAPA_timeseries=np.nanmean( np.nanmean(SST_monthlymean_PAPA,axis=2) , axis=1); SST_monthlymean_PAPA_timeseries=np.squeeze(SST_monthlymean_PAPA_timeseries)

SpCO2_monthlymean_PAPA = SpCO2_monthlymean[:,135:140,210:220]
SpCO2_monthlymean_PAPA_timeseries=np.nanmean( np.nanmean(SpCO2_monthlymean_PAPA,axis=2) , axis=1); SpCO2_monthlymean_PAPA_timeseries=np.squeeze(SpCO2_monthlymean_PAPA_timeseries)

SpCO2_monthlymean_TEMP_eq1_PAPA = SpCO2_monthlymean_TEMP_eq1[:,135:140,210:220]
SpCO2_monthlymean_TEMP_eq1_PAPA_timeseries=np.nanmean( np.nanmean(SpCO2_monthlymean_TEMP_eq1_PAPA,axis=2) , axis=1); SpCO2_monthlymean_TEMP_eq1_PAPA_timeseries=np.squeeze(SpCO2_monthlymean_TEMP_eq1_PAPA_timeseries)

SpCO2_monthlymean_NonTEMP_eq2_PAPA = SpCO2_monthlymean_NonTEMP_eq2[:,135:140,210:220]
SpCO2_monthlymean_NonTEMP_eq2_PAPA_timeseries=np.nanmean( np.nanmean(SpCO2_monthlymean_NonTEMP_eq2_PAPA,axis=2) , axis=1); SpCO2_monthlymean_NonTEMP_eq2_PAPA_timeseries=np.squeeze(SpCO2_monthlymean_NonTEMP_eq2_PAPA_timeseries)


SST_monthlymean_NABE = SST_monthlymean[:,135:140,325:335]
SST_monthlymean_NABE_timeseries=np.nanmean( np.nanmean(SST_monthlymean_NABE,axis=2) , axis=1); SST_monthlymean_NABE_timeseries=np.squeeze(SST_monthlymean_NABE_timeseries)

SpCO2_monthlymean_NABE = SpCO2_monthlymean[:,135:140,325:335]
SpCO2_monthlymean_NABE_timeseries=np.nanmean( np.nanmean(SpCO2_monthlymean_NABE,axis=2) , axis=1); SpCO2_monthlymean_NABE_timeseries=np.squeeze(SpCO2_monthlymean_NABE_timeseries)

SpCO2_monthlymean_TEMP_eq1_NABE = SpCO2_monthlymean_TEMP_eq1[:,135:140,325:335]
SpCO2_monthlymean_TEMP_eq1_NABE_timeseries=np.nanmean( np.nanmean(SpCO2_monthlymean_TEMP_eq1_NABE,axis=2) , axis=1); SpCO2_monthlymean_TEMP_eq1_NABE_timeseries=np.squeeze(SpCO2_monthlymean_TEMP_eq1_NABE_timeseries)

SpCO2_monthlymean_NonTEMP_eq2_NABE = SpCO2_monthlymean_NonTEMP_eq2[:,135:140,325:335]
SpCO2_monthlymean_NonTEMP_eq2_NABE_timeseries=np.nanmean( np.nanmean(SpCO2_monthlymean_NonTEMP_eq2_NABE,axis=2) , axis=1); SpCO2_monthlymean_NonTEMP_eq2_NABE_timeseries=np.squeeze(SpCO2_monthlymean_NonTEMP_eq2_NABE_timeseries)

###############################################################################
#Corr_SpCO2_Chl = empty((lat_n_regrid,lon_n_regrid)) * nan
#Pvalue_SpCO2_Chl = empty((lat_n_regrid,lon_n_regrid)) * nan
#Bslope_SpCO2_Chl = empty((lat_n_regrid,lon_n_regrid)) * nan
#
#Var_y=SpCO2_allmonths
#Var_x=np.log10(Chl_allmonths)
#for ii in range(0,lat_n_regrid):
#    for jj in range(0,lon_n_regrid):
#        if np.logical_not(np.isnan(Var_x[0,ii,jj])): # If this grid cell has data and is non NaN
#            if Ocean_Index[ii,jj] != 10: # Excluding the Arctic Ocean                
#                xx=Var_x[:,ii,jj] ; yy=Var_y[:,ii,jj] # Creating a time series of the variabel for that grid cell                
#                yy=yy[np.logical_not(np.isnan(xx))] ; xx=xx[np.logical_not(np.isnan(xx))] # Excluding any NaN values in the time series                   
#                xx=xx[np.logical_not(np.isnan(yy))] ; yy=yy[np.logical_not(np.isnan(yy))] # Excluding any NaN values in the time series 
#                slope_ij, intercept_ij, r_value_ij, p_value_ij, std_err_ij = stats.linregress(xx, yy)
#
#                Corr_SpCO2_Chl [ii,jj] = r_value_ij
#                Pvalue_SpCO2_Chl [ii,jj] = p_value_ij
#                Bslope_SpCO2_Chl [ii,jj] = slope_ij          

###############################################################################
##############            Correlation Calculations            #################   
from Behzadlib import func_corr_3d

###############################################################################
##############   Spco2 vs Chlorophyll (Chl) Correlations   ####################
Var_y=SpCO2_allmonths
Var_x=np.log10(Chl_allmonths)
Corr_xy, Pvalue_xy = func_corr_3d(Var_y, Var_x);  Corr_xy[Ocean_Index==10]=nan; Pvalue_xy[Ocean_Index==10]=nan;   #  Excluding the Arctic Ocean

Corr_SpCO2_Chl=Corr_xy
Pvalue_SpCO2_Chl=Pvalue_xy      
##############
Var_y=SpCO2_allmonths_TEMP_eq1
Var_x=np.log10(Chl_allmonths)
Corr_xy, Pvalue_xy = func_corr_3d(Var_y, Var_x);  Corr_xy[Ocean_Index==10]=nan; Pvalue_xy[Ocean_Index==10]=nan;   #  Excluding the Arctic Ocean

Corr_SpCO2_TEMP_eq1_Chl = Corr_xy
Pvalue_SpCO2_TEMP_eq1_Chl = Pvalue_xy
##############
Var_y=SpCO2_allmonths_NonTEMP_eq2
Var_x=np.log10(Chl_allmonths)
Corr_xy, Pvalue_xy = func_corr_3d(Var_y, Var_x);  Corr_xy[Ocean_Index==10]=nan; Pvalue_xy[Ocean_Index==10]=nan;   #  Excluding the Arctic Ocean

Corr_SpCO2_NonTEMP_eq2_Chl = Corr_xy
Pvalue_SpCO2_NonTEMP_eq2_Chl = Pvalue_xy
###############################################################################
##############    Spco2 vs Biomass (Phyc) Correlations   ######################
Var_y=SpCO2_allmonths
Var_x=np.log10(Phyc_allmonths)
Corr_xy, Pvalue_xy = func_corr_3d(Var_y, Var_x);  Corr_xy[Ocean_Index==10]=nan; Pvalue_xy[Ocean_Index==10]=nan;   #  Excluding the Arctic Ocean

Corr_SpCO2_Phyc = Corr_xy
Pvalue_SpCO2_Phyc = empty((lat_n_regrid,lon_n_regrid)) * nan
##############
Var_y=SpCO2_allmonths_TEMP_eq1
Var_x=np.log10(Phyc_allmonths)
Corr_xy, Pvalue_xy = func_corr_3d(Var_y, Var_x);  Corr_xy[Ocean_Index==10]=nan; Pvalue_xy[Ocean_Index==10]=nan;   #  Excluding the Arctic Ocean

Corr_SpCO2_TEMP_eq1_Phyc = Corr_xy
Pvalue_SpCO2_TEMP_eq1_Phyc = Pvalue_xy
##############
Var_y=SpCO2_allmonths_NonTEMP_eq2
Var_x=np.log10(Phyc_allmonths)
Corr_xy, Pvalue_xy = func_corr_3d(Var_y, Var_x);  Corr_xy[Ocean_Index==10]=nan; Pvalue_xy[Ocean_Index==10]=nan;   #  Excluding the Arctic Ocean

Corr_SpCO2_NonTEMP_eq2_Phyc = Corr_xy
Pvalue_SpCO2_NonTEMP_eq2_Phyc = Pvalue_xy

###############################################################################
############## Spco2 vs Biomass (Phydiat_prcnt) Correlations ##################
Var_y=SpCO2_allmonths
Var_x=np.log10(Phydiat_prcnt_allmonths)
Corr_xy, Pvalue_xy = func_corr_3d(Var_y, Var_x);  Corr_xy[Ocean_Index==10]=nan; Pvalue_xy[Ocean_Index==10]=nan;   #  Excluding the Arctic Ocean

Corr_SpCO2_Phydiat_prcnt = Corr_xy
Pvalue_SpCO2_Phydiat_prcnt = Pvalue_xy
##############
Var_y=SpCO2_allmonths_TEMP_eq1
Var_x=np.log10(Phydiat_prcnt_allmonths)
Corr_xy, Pvalue_xy = func_corr_3d(Var_y, Var_x);  Corr_xy[Ocean_Index==10]=nan; Pvalue_xy[Ocean_Index==10]=nan;   #  Excluding the Arctic Ocean

Corr_SpCO2_TEMP_eq1_Phydiat_prcnt = Corr_xy
Pvalue_SpCO2_TEMP_eq1_Phydiat_prcnt = Pvalue_xy
##############
Var_y=SpCO2_allmonths_NonTEMP_eq2
Var_x=np.log10(Phydiat_prcnt_allmonths)
Corr_xy, Pvalue_xy = func_corr_3d(Var_y, Var_x);  Corr_xy[Ocean_Index==10]=nan; Pvalue_xy[Ocean_Index==10]=nan;   #  Excluding the Arctic Ocean

Corr_SpCO2_NonTEMP_eq2_Phydiat_prcnt = Corr_xy
Pvalue_SpCO2_NonTEMP_eq2_Phydiat_prcnt = Pvalue_xy

###############################################################################
##############    Regression Calculations - Deseasonalized    #################

###############         Deseasonalization of the variables         ############ 
help_t= np.tile(SpCO2_monthlymean,(yrs_n,1,1)) # Repeats the matrix for the number of times specified for each axis # Creates a n-years matrix with the monthly mean values repeated N times
SpCO2_allmonths_des = SpCO2_allmonths - help_t # Deseasonalized

help_t = SpCO2_allmonths_TEMP_eq1.reshape(( yrs_n ,12,lat_n_regrid,lon_n_regrid)) # Monthly means of the Variable over the time period
help_t = np.nanmean(help_t,axis=0)
help_t=np.squeeze(help_t)
help_t= np.tile(help_t,(yrs_n,1,1)) # Repeats the matrix for the number of times specified for each axis # Creates a n-years matrix with the monthly mean values repeated N times
SpCO2_allmonths_TEMP_eq1_des = SpCO2_allmonths_TEMP_eq1 - help_t # Deseasonalized

help_t = SpCO2_allmonths_NonTEMP_eq2.reshape(( yrs_n ,12,lat_n_regrid,lon_n_regrid)) # Monthly means of the Variable over the time period
help_t = np.nanmean(help_t,axis=0)
help_t=np.squeeze(help_t)
help_t= np.tile(help_t,(yrs_n,1,1)) # Repeats the matrix for the number of times specified for each axis # Creates a n-years matrix with the monthly mean values repeated N times
SpCO2_allmonths_NonTEMP_eq2_des = SpCO2_allmonths_NonTEMP_eq2 - help_t # Deseasonalized

SpCO2_allmonths_NonTEMP_subtracted_des = SpCO2_allmonths_des - SpCO2_allmonths_TEMP_eq1_des


help_t = Chl_allmonths.reshape(( yrs_n ,12,lat_n_regrid,lon_n_regrid)) # Monthly means of the Variable over the time period
help_t = np.nanmean(help_t,axis=0)
help_t=np.squeeze(help_t)
help_t= np.tile(help_t,(yrs_n,1,1)) # Repeats the matrix for the number of times specified for each axis # Creates a n-years matrix with the monthly mean values repeated N times
Chl_allmonths_des = Chl_allmonths - help_t # Deseasonalized

help_t = Phyc_allmonths.reshape(( yrs_n ,12,lat_n_regrid,lon_n_regrid)) # Monthly means of the Variable over the time period
help_t = np.nanmean(help_t,axis=0)
help_t=np.squeeze(help_t)
help_t= np.tile(help_t,(yrs_n,1,1)) # Repeats the matrix for the number of times specified for each axis # Creates a n-years matrix with the monthly mean values repeated N times
Phyc_allmonths_des = Phyc_allmonths - help_t # Deseasonalized

help_t = Phydiat_prcnt_allmonths.reshape(( yrs_n ,12,lat_n_regrid,lon_n_regrid)) # Monthly means of the Variable over the time period
help_t = np.nanmean(help_t,axis=0)
help_t=np.squeeze(help_t)
help_t= np.tile(help_t,(yrs_n,1,1)) # Repeats the matrix for the number of times specified for each axis # Creates a n-years matrix with the monthly mean values repeated N times
Phydiat_prcnt_allmonths_des = Phydiat_prcnt_allmonths - help_t # Deseasonalized


###############################################################################
##############   Spco2 vs Chlorophyll (Chl) Correlations   ####################
Var_y=SpCO2_allmonths_des
Var_x=np.log10(Chl_allmonths_des)
Corr_xy, Pvalue_xy = func_corr_3d(Var_y, Var_x);  Corr_xy[Ocean_Index==10]=nan; Pvalue_xy[Ocean_Index==10]=nan;   #  Excluding the Arctic Ocean

Corr_SpCO2_Chl_des = Corr_xy
Pvalue_SpCO2_Chl_des = Pvalue_xy
##############
Var_y=SpCO2_allmonths_TEMP_eq1_des
Var_x=np.log10(Chl_allmonths_des)
Corr_xy, Pvalue_xy = func_corr_3d(Var_y, Var_x);  Corr_xy[Ocean_Index==10]=nan; Pvalue_xy[Ocean_Index==10]=nan;   #  Excluding the Arctic Ocean

Corr_SpCO2_TEMP_eq1_Chl_des = Corr_xy
Pvalue_SpCO2_TEMP_eq1_Chl_des = Pvalue_xy
##############
Var_y=SpCO2_allmonths_NonTEMP_eq2_des
Var_x=np.log10(Chl_allmonths_des)
Corr_xy, Pvalue_xy = func_corr_3d(Var_y, Var_x);  Corr_xy[Ocean_Index==10]=nan; Pvalue_xy[Ocean_Index==10]=nan;   #  Excluding the Arctic Ocean

Corr_SpCO2_NonTEMP_eq2_Chl_des = Corr_xy
Pvalue_SpCO2_NonTEMP_eq2_Chl_des = Pvalue_xy

###############################################################################
##############    Spco2 vs Biomass (Phyc) Correlations   ######################
Var_y=SpCO2_allmonths_des
Var_x=np.log10(Phyc_allmonths_des)
Corr_xy, Pvalue_xy = func_corr_3d(Var_y, Var_x);  Corr_xy[Ocean_Index==10]=nan; Pvalue_xy[Ocean_Index==10]=nan;   #  Excluding the Arctic Ocean

Corr_SpCO2_Phyc_des = Corr_xy
Pvalue_SpCO2_Phyc_des = Pvalue_xy
##############
Var_y=SpCO2_allmonths_TEMP_eq1_des
Var_x=np.log10(Phyc_allmonths_des)
Corr_xy, Pvalue_xy = func_corr_3d(Var_y, Var_x);  Corr_xy[Ocean_Index==10]=nan; Pvalue_xy[Ocean_Index==10]=nan;   #  Excluding the Arctic Ocean

Corr_SpCO2_TEMP_eq1_Phyc_des = Corr_xy
Pvalue_SpCO2_TEMP_eq1_Phyc_des = Pvalue_xy
##############
Var_y=SpCO2_allmonths_NonTEMP_eq2_des
Var_x=np.log10(Phyc_allmonths_des)
Corr_xy, Pvalue_xy = func_corr_3d(Var_y, Var_x);  Corr_xy[Ocean_Index==10]=nan; Pvalue_xy[Ocean_Index==10]=nan;   #  Excluding the Arctic Ocean

Corr_SpCO2_NonTEMP_eq2_Phyc_des = Corr_xy
Pvalue_SpCO2_NonTEMP_eq2_Phyc_des = Pvalue_xy

###############################################################################
############## Spco2 vs Biomass (Phydiat_prcnt) Correlations ##################
Var_y=SpCO2_allmonths_des
Var_x=np.log10(Phydiat_prcnt_allmonths_des)
Corr_xy, Pvalue_xy = func_corr_3d(Var_y, Var_x);  Corr_xy[Ocean_Index==10]=nan; Pvalue_xy[Ocean_Index==10]=nan;   #  Excluding the Arctic Ocean

Corr_SpCO2_Phydiat_prcnt_des = Corr_xy
Pvalue_SpCO2_Phydiat_prcnt_des = Pvalue_xy
##############
Var_y=SpCO2_allmonths_TEMP_eq1_des
Var_x=np.log10(Phydiat_prcnt_allmonths_des)
Corr_xy, Pvalue_xy = func_corr_3d(Var_y, Var_x);  Corr_xy[Ocean_Index==10]=nan; Pvalue_xy[Ocean_Index==10]=nan;   #  Excluding the Arctic Ocean

Corr_SpCO2_TEMP_eq1_Phydiat_prcnt_des = Corr_xy
Pvalue_SpCO2_TEMP_eq1_Phydiat_prcnt_des = Pvalue_xy
##############
Var_y=SpCO2_allmonths_NonTEMP_eq2_des
Var_x=np.log10(Phydiat_prcnt_allmonths_des)
Corr_xy, Pvalue_xy = func_corr_3d(Var_y, Var_x);  Corr_xy[Ocean_Index==10]=nan; Pvalue_xy[Ocean_Index==10]=nan;   #  Excluding the Arctic Ocean

Corr_SpCO2_NonTEMP_eq2_Phydiat_prcnt_des = Corr_xy
Pvalue_SpCO2_NonTEMP_eq2_Phydiat_prcnt_des = Pvalue_xy

                
###############################################################################                
###############################################################################                
###############################################################################                
###############################################################################               
###############################################################################
##############            Correlation Calculations            #################          
                
###############################################################################
##############    MLD vs Chlorophyll (Chl) Correlations    ####################
Var_y=MLD_allmonths
Var_x=np.log10(Chl_allmonths)
Corr_xy, Pvalue_xy = func_corr_3d(Var_y, Var_x);  Corr_xy[Ocean_Index==10]=nan; Pvalue_xy[Ocean_Index==10]=nan;   #  Excluding the Arctic Ocean

Corr_MLD_Chl = Corr_xy
Pvalue_MLD_Chl = Pvalue_xy
##############          MLD vs Spco2 Correlations          ####################
Var_y=MLD_allmonths
Var_x=SpCO2_allmonths
Corr_xy, Pvalue_xy = func_corr_3d(Var_y, Var_x);  Corr_xy[Ocean_Index==10]=nan; Pvalue_xy[Ocean_Index==10]=nan;   #  Excluding the Arctic Ocean

Corr_MLD_SpCO2 = Corr_xy
Pvalue_MLD_SpCO2 = Pvalue_xy
##############          MLD vs SST Correlations           ####################
Var_y=MLD_allmonths
Var_x=SST_allmonths
Corr_xy, Pvalue_xy = func_corr_3d(Var_y, Var_x);  Corr_xy[Ocean_Index==10]=nan; Pvalue_xy[Ocean_Index==10]=nan;   #  Excluding the Arctic Ocean

Corr_MLD_SST = Corr_xy
Pvalue_MLD_SST = Pvalue_xy

###############################################################################
############       Biomass (Phyc) vs MLD Correlations          ################
Var_y=np.log10(Phyc_allmonths)
Var_x=MLD_allmonths
Corr_xy, Pvalue_xy = func_corr_3d(Var_y, Var_x);  Corr_xy[Ocean_Index==10]=nan; Pvalue_xy[Ocean_Index==10]=nan;   #  Excluding the Arctic Ocean

Corr_Phyc_MLD = Corr_xy
Pvalue_Phyc_MLD = Pvalue_xy
#############       Biomass (Phyc) vs SST Correlations       ##################
Var_y=np.log10(Phyc_allmonths)
Var_x=SST_allmonths
Corr_xy, Pvalue_xy = func_corr_3d(Var_y, Var_x);  Corr_xy[Ocean_Index==10]=nan; Pvalue_xy[Ocean_Index==10]=nan;   #  Excluding the Arctic Ocean

Corr_Phyc_SST = Corr_xy
Pvalue_Phyc_SST = Pvalue_xy
#########    Biomass (Phyc) vs Chlorophyll (Chl) Correlations    ##############
Var_y=np.log10(Phyc_allmonths)
Var_x=np.log10(Chl_allmonths)
Corr_xy, Pvalue_xy = func_corr_3d(Var_y, Var_x);  Corr_xy[Ocean_Index==10]=nan; Pvalue_xy[Ocean_Index==10]=nan;   #  Excluding the Arctic Ocean

Corr_Phyc_Chl = Corr_xy
Pvalue_Phyc_Chl = Pvalue_xy

###############################################################################
############          Heatflux vs Biomass Correlations         ################
Var_y=np.log10(Phyc_allmonths)
Var_x=HFDS_allmonths
Corr_xy, Pvalue_xy = func_corr_3d(Var_y, Var_x);  Corr_xy[Ocean_Index==10]=nan; Pvalue_xy[Ocean_Index==10]=nan;   #  Excluding the Arctic Ocean

Corr_Phyc_HFDS = Corr_xy
Pvalue_Phyc_HFDS = Pvalue_xy
############            Heatflux vs Spco2 Correlations         ################
Var_y=SpCO2_allmonths
Var_x=HFDS_allmonths
Corr_xy, Pvalue_xy = func_corr_3d(Var_y, Var_x);  Corr_xy[Ocean_Index==10]=nan; Pvalue_xy[Ocean_Index==10]=nan;   #  Excluding the Arctic Ocean

Corr_SpCO2_HFDS = Corr_xy
Pvalue_SpCO2_HFDS = Pvalue_xy
############             Heatflux vs SST Correlations          ################
Var_y=HFDS_allmonths
Var_x=SST_allmonths
Corr_xy, Pvalue_xy = func_corr_3d(Var_y, Var_x);  Corr_xy[Ocean_Index==10]=nan; Pvalue_xy[Ocean_Index==10]=nan;   #  Excluding the Arctic Ocean

Corr_HFDS_SST = Corr_xy
Pvalue_HFDS_SST = Pvalue_xy

###############################################################################
##############    Regression Calculations - Deseasonalized    #################

###############         Deseasonalization of the variables         ############ 

help_t = MLD_allmonths.reshape(( yrs_n ,12,lat_n_regrid,lon_n_regrid)) # Monthly means of the Variable over the time period
help_t = np.nanmean(help_t,axis=0)
help_t=np.squeeze(help_t)
help_t= np.tile(help_t,(yrs_n,1,1)) # Repeats the matrix for the number of times specified for each axis # Creates a n-years matrix with the monthly mean values repeated N times
MLD_allmonths_des = MLD_allmonths - help_t # Deseasonalized

help_t = SST_allmonths.reshape(( yrs_n ,12,lat_n_regrid,lon_n_regrid)) # Monthly means of the Variable over the time period
help_t = np.nanmean(help_t,axis=0)
help_t=np.squeeze(help_t)
help_t= np.tile(help_t,(yrs_n,1,1)) # Repeats the matrix for the number of times specified for each axis # Creates a n-years matrix with the monthly mean values repeated N times
SST_allmonths_des = SST_allmonths - help_t # Deseasonalized

help_t = HFDS_allmonths.reshape(( yrs_n ,12,lat_n_regrid,lon_n_regrid)) # Monthly means of the Variable over the time period
help_t = np.nanmean(help_t,axis=0)
help_t=np.squeeze(help_t)
help_t= np.tile(help_t,(yrs_n,1,1)) # Repeats the matrix for the number of times specified for each axis # Creates a n-years matrix with the monthly mean values repeated N times
HFDS_allmonths_des = HFDS_allmonths - help_t # Deseasonalized


###############################################################################
##############    MLD vs Chlorophyll (Chl) Correlations    ####################
Var_y=MLD_allmonths_des
Var_x=np.log10(Chl_allmonths_des)
Corr_xy, Pvalue_xy = func_corr_3d(Var_y, Var_x);  Corr_xy[Ocean_Index==10]=nan; Pvalue_xy[Ocean_Index==10]=nan;   #  Excluding the Arctic Ocean

Corr_MLD_Chl_des = Corr_xy
Pvalue_MLD_Chl_des = Pvalue_xy
##############          MLD vs Spco2 Correlations          ####################
Var_y=MLD_allmonths_des
Var_x=SpCO2_allmonths_des
Corr_xy, Pvalue_xy = func_corr_3d(Var_y, Var_x);  Corr_xy[Ocean_Index==10]=nan; Pvalue_xy[Ocean_Index==10]=nan;   #  Excluding the Arctic Ocean

Corr_MLD_SpCO2_des = Corr_xy
Pvalue_MLD_SpCO2_des = Pvalue_xy
##############          MLD vs SST Correlations           ####################
Var_y=MLD_allmonths_des
Var_x=SST_allmonths_des
Corr_xy, Pvalue_xy = func_corr_3d(Var_y, Var_x);  Corr_xy[Ocean_Index==10]=nan; Pvalue_xy[Ocean_Index==10]=nan;   #  Excluding the Arctic Ocean

Corr_MLD_SST_des = Corr_xy
Pvalue_MLD_SST_des = Pvalue_xy

###############################################################################
############       Biomass (Phyc) vs MLD Correlations          ################
Var_y=np.log10(Phyc_allmonths_des)
Var_x=MLD_allmonths_des
Corr_xy, Pvalue_xy = func_corr_3d(Var_y, Var_x);  Corr_xy[Ocean_Index==10]=nan; Pvalue_xy[Ocean_Index==10]=nan;   #  Excluding the Arctic Ocean

Corr_Phyc_MLD_des = Corr_xy
Pvalue_Phyc_MLD_des = Pvalue_xy
#############       Biomass (Phyc) vs SST Correlations       ##################
Var_y=np.log10(Phyc_allmonths_des)
Var_x=SST_allmonths_des
Corr_xy, Pvalue_xy = func_corr_3d(Var_y, Var_x);  Corr_xy[Ocean_Index==10]=nan; Pvalue_xy[Ocean_Index==10]=nan;   #  Excluding the Arctic Ocean 

Corr_Phyc_SST_des = Corr_xy
Pvalue_Phyc_SST_des = Pvalue_xy
#########    Biomass (Phyc) vs Chlorophyll (Chl) Correlations    ##############
Var_y=np.log10(Phyc_allmonths_des)
Var_x=np.log10(Chl_allmonths_des)
Corr_xy, Pvalue_xy = func_corr_3d(Var_y, Var_x);  Corr_xy[Ocean_Index==10]=nan; Pvalue_xy[Ocean_Index==10]=nan;   #  Excluding the Arctic Ocean

Corr_Phyc_Chl_des = Corr_xy
Pvalue_Phyc_Chl_des = Pvalue_xy

###############################################################################
############          Heatflux vs Biomass Correlations         ################
Var_y=np.log10(Phyc_allmonths_des)
Var_x=HFDS_allmonths_des
Corr_xy, Pvalue_xy = func_corr_3d(Var_y, Var_x);  Corr_xy[Ocean_Index==10]=nan; Pvalue_xy[Ocean_Index==10]=nan;   #  Excluding the Arctic Ocean

Corr_Phyc_HFDS_des = Corr_xy
Pvalue_Phyc_HFDS_des = Pvalue_xy
############            Heatflux vs Spco2 Correlations         ################
Var_y=SpCO2_allmonths_des
Var_x=HFDS_allmonths_des
Corr_xy, Pvalue_xy = func_corr_3d(Var_y, Var_x);  Corr_xy[Ocean_Index==10]=nan; Pvalue_xy[Ocean_Index==10]=nan;   #  Excluding the Arctic Ocean

Corr_SpCO2_HFDS_des = Corr_xy
Pvalue_SpCO2_HFDS_des = Pvalue_xy
############             Heatflux vs SST Correlations          ################
Var_y=HFDS_allmonths_des
Var_x=SST_allmonths_des
Corr_xy, Pvalue_xy = func_corr_3d(Var_y, Var_x);  Corr_xy[Ocean_Index==10]=nan; Pvalue_xy[Ocean_Index==10]=nan;   #  Excluding the Arctic Ocean

Corr_HFDS_SST_des = Corr_xy
Pvalue_HFDS_SST_des = Pvalue_xy
         

###############################################################################
###########                 Saving Results                #####################
#var_list_1=dir() # You need to edit this list to exclude the libraries and Datastes
#var_list_1=[item for item in var_list_1 if not item.startswith("_") and not item.startswith("__") and not item.startswith("dset")]               
#np.save('var_list_1.npy',var_list_1)
###########          
import os
import shelve

#filename_out = (dir_results + 'AllResults_'+GCM+'pco2_1981_2010.out') # Directory to save processed data
filename_out = (dir_results + 'AllResults_'+GCM+'_pco2_'+str(start_date_plt)+'_'+str(end_date_plt)+'.out') # Directory to save processed data
var_list_1=np.load('var_list_1.npy')

### To save
my_shelf = shelve.open(filename_out,'n') # 'n' for new

for key in var_list_1:
    try:
        my_shelf[key] = globals()[key]
    except TypeError:
        print('ERROR shelving: {0}'.format(key))
my_shelf.close()
##############################################################################
#################        To restore:        ##################################
#import os
#import shelve
#
##filename_out = (dir_results + 'AllResults_'+GCM+'pco2_1981_2010.out') # Directory to save processed data
#filename_out = (dir_results + 'AllResults_'+GCM+'pco2_'+str(start_date_plt)+'_'+str(end_date_plt)+'.out') # Directory to save processed data
#
#my_shelf = shelve.open(filename_out)
#for key in my_shelf:
#    globals()[key]=my_shelf[key]
#my_shelf.close()
##############################################################################



###############################################################################
##############                PAPA - NABE plots               #################
###############################################################################
Time_months = ('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec');
P_Var_x = np.linspace(1,12,12);

P_title= 'Monthly pCO2 and temperature averages - '+str(GCM)+' , '+str(start_date_plt)+'-'+str(end_date_plt)
fig=plt.figure()  

ax = fig.add_subplot(2,1,1)   
plt.title('PAPA [45N-50N, 140W-150W]', fontsize=24)
l = plt.axhline(y=  np.nanmean( SpCO2_monthlymean_PAPA_timeseries, axis=0)    , color='k', linewidth=3.0)
im1=plt.plot(P_Var_x, SpCO2_monthlymean_PAPA_timeseries, c='r',  label='pCO2', marker='o', markersize=15, markerfacecolor='W', linewidth=5.0)
im2=plt.plot(P_Var_x, SpCO2_monthlymean_TEMP_eq1_PAPA_timeseries, c='b',  label='pCO2-T', marker='o', markersize=15, markerfacecolor='W', linewidth=5.0)
im3=plt.plot(P_Var_x, SpCO2_monthlymean_NonTEMP_eq2_PAPA_timeseries, c='g',  label='pCO2-nonT', marker='o', markersize=15, markerfacecolor='W', linewidth=5.0)
plt.xticks(fontsize = 26); plt.yticks(fontsize = 26)
#plt.xlabel('Month', fontsize=26)
plt.ylabel('pCO2', fontsize=26)
plt.ylim(275,450)
plt.xticks(P_Var_x, '', rotation=45) # Setting the X-ticks ## rotation='horizontal' # rotation='vertical'
plt.legend(prop={'size': 12}, loc='best', fancybox=True, framealpha=0.8)

ax = fig.add_subplot(2,1,2)   
plt.title('NABE [45N-50N, 25W-35W]', fontsize=24)
l = plt.axhline(y=  np.nanmean( SpCO2_monthlymean_NABE_timeseries, axis=0)    , color='k', linewidth=3.0)
im1=plt.plot(P_Var_x, SpCO2_monthlymean_NABE_timeseries, c='r',  label='pCO2', marker='o', markersize=15, markerfacecolor='W', linewidth=5.0)
im2=plt.plot(P_Var_x, SpCO2_monthlymean_TEMP_eq1_NABE_timeseries, c='b',  label='pCO2-T', marker='o', markersize=15, markerfacecolor='W', linewidth=5.0)
im3=plt.plot(P_Var_x, SpCO2_monthlymean_NonTEMP_eq2_NABE_timeseries, c='g',  label='pCO2-nonT', marker='o', markersize=15, markerfacecolor='W', linewidth=5.0)
plt.xticks(fontsize = 26); plt.yticks(fontsize = 26)
plt.xlabel('Month', fontsize=26)
plt.ylabel('pCO2', fontsize=26)
plt.ylim(275,450)
plt.xticks(P_Var_x, Time_months, rotation=45, fontsize=24) # Setting the X-ticks ## rotation='horizontal' # rotation='vertical'
#plt.legend(prop={'size': 14}, loc='best', fancybox=True, framealpha=0.8)

#plt.subplots_adjust(left=0.1, bottom=0.11, right=0.9, top=0.9, hspace=0.4, wspace=0.05) # the amount of height/width reserved for space between subplots
plt.subplots_adjust(left=0.33, bottom=0.11, right=0.69, top=0.9, hspace=0.2, wspace=0.05) # the amount of height/width reserved for space between subplots
plt.suptitle(P_title, fontsize=26)    
plt.show()
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full    

fig.savefig(dir_figs+str(GCM)+'_PAPA_NABE_TimeSeries_pCO2_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')



P_title='Monthly pCO2 and temperature averages - '+str(GCM)+' , '+str(start_date_plt)+'-'+str(end_date_plt)
fig=plt.figure()  

ax = fig.add_subplot(2,1,1)   
plt.title('PAPA [45N-50N, 140W-150W]', fontsize=24)
im1=plt.plot(P_Var_x, SST_monthlymean_PAPA_timeseries, c='y', marker='D', markersize=15, markerfacecolor='W', linewidth=5.0)
plt.xticks(fontsize = 26); plt.yticks(fontsize = 26)
#plt.xlabel('Month', fontsize=26)
plt.ylabel('Temp', color='y', fontsize=26)
plt.tick_params(axis='y', labelcolor='y')
plt.ylim(4,20)
plt.xticks(P_Var_x, '', rotation=45) # Setting the X-ticks ## rotation='horizontal' # rotation='vertical'

ax = fig.add_subplot(2,1,2)   
plt.title('NABE [45N-50N, 25W-35W]', fontsize=24)
im1=plt.plot(P_Var_x, SST_monthlymean_NABE_timeseries, c='y', marker='D', markersize=15, markerfacecolor='W', linewidth=5.0)
plt.xticks(fontsize = 26); plt.yticks(fontsize = 26)
plt.xlabel('Month', fontsize=26)
plt.ylabel('Temp', color='y', fontsize=26)
plt.tick_params(axis='y', labelcolor='y')
plt.ylim(4,20)
plt.xticks(P_Var_x, Time_months, rotation=45, fontsize=24) # Setting the X-ticks ## rotation='horizontal' # rotation='vertical'

plt.subplots_adjust(left=0.33, bottom=0.11, right=0.69, top=0.9, hspace=0.2, wspace=0.05) # the amount of height/width reserved for space between subplots
plt.suptitle(P_title, fontsize=26)    
plt.show()
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full      

fig.savefig(dir_figs+str(GCM)+'_PAPA_NABE_TimeSeries_TEMP_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')



###############################################################################
#############              pCo2 Combination Plots          ####################
###############################################################################
var_names=['Corr_SpCO2_Chl','Corr_SpCO2_Phyc','Corr_SpCO2_Phydiat_prcnt','Corr_SpCO2_TEMP_eq1_Chl','Corr_SpCO2_TEMP_eq1_Phyc','Corr_SpCO2_TEMP_eq1_Phydiat_prcnt','Corr_SpCO2_NonTEMP_eq2_Chl','Corr_SpCO2_NonTEMP_eq2_Phyc','Corr_SpCO2_NonTEMP_eq2_Phydiat_prcnt']
var_names_pval=['Pvalue_SpCO2_Chl','Pvalue_SpCO2_Phyc','Pvalue_SpCO2_Phydiat_prcnt','Pvalue_SpCO2_TEMP_eq1_Chl','Pvalue_SpCO2_TEMP_eq1_Phyc','Pvalue_SpCO2_TEMP_eq1_Phydiat_prcnt','Pvalue_SpCO2_NonTEMP_eq2_Chl','Pvalue_SpCO2_NonTEMP_eq2_Phyc','Pvalue_SpCO2_NonTEMP_eq2_Phydiat_prcnt']
plot_titles=['pCO2 vs. chl','pCO2 vs. biomass','pCO2 vs. % large phyto','pCO2-T vs. chl','pCO2-T vs. biomass','pCO2-T vs. % large phyto','pCO2-nonT vs. chl','pCO2-nonT vs. biomass','pCO2-nonT vs. % large phyto']

P_cmap=plt.cm.bwr; P_c_fill='fill'; P_Lon=Lon_regrid_2D; P_Lat=Lat_regrid_2D; P_range=np.linspace(-1.,1.,51);
P_proj='cyl'; P_lon0=210.; P_latN=90.; P_latS=-90.; 

n_r=3 ; n_c=3 ; n_t=9
P_title='Correlations of pCO2 and log Chl, log Biomass, % large Phyto - '+str(GCM)+' , '+str(start_date_plt)+'-'+str(end_date_plt)
fig=plt.figure()    
for ii in range (0, n_t):
    
    Var_plot_ii=copy.deepcopy(eval(var_names[ii]))
    
    ax = fig.add_subplot(n_r,n_c,ii+1)   
    m = Basemap( projection=P_proj, lon_0=P_lon0, llcrnrlon=P_lon0-180, llcrnrlat=P_latS, urcrnrlon=P_lon0+180, urcrnrlat=P_latN)    
    if ii == (n_c*(n_r-1)): # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        m.drawparallels(np.arange(P_latS, P_latN+0.001, 30.),labels=[True,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Latitutes
        m.drawmeridians(np.arange(0,360,90.),labels=[False,False,False,True], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Longitudes    
    elif ii==0 or ii==n_c or ii==n_c*2 or ii==n_c*3 or ii==n_c*4 or ii==n_c*5 or ii==n_c*6 or ii==n_c*7 or ii==n_c*8:
        m.drawparallels(np.arange(P_latS, P_latN+0.001, 30.),labels=[True,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Latitutes
    elif ii >= n_t-n_c and ii != (n_c*(n_r-1)): # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        m.drawmeridians(np.arange(0,360,90.),labels=[False,False,False,True], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Longitudes
    else:
        m.drawparallels(np.arange(P_latS, P_latN+0.001, 30.),labels=[False,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Latitutes
        m.drawmeridians(np.arange(0,360,90.),labels=[False,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Longitudes
    
    if P_c_fill=='fill':
        m.fillcontinents(color='0')
    m.drawcoastlines(linewidth=1.0, linestyle='solid', color='k', antialiased=1, ax=None, zorder=None)
    im=m.contourf(P_Lon, P_Lat, Var_plot_ii,P_range,latlon=True, cmap=P_cmap) #, extend='both')
    plt.title(plot_titles[ii], fontsize=16)
    
plt.subplots_adjust(left=0.07, bottom=0.1, right=0.92, top=0.92, hspace=0.1, wspace=0.02) # the amount of height/width reserved for space between subplots
cbar_ax = fig.add_axes([0.93, 0.11, 0.015, 0.8]) # [right,bottom,width,height] 
fig.colorbar(im, cax=cbar_ax, ticks=[-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1])
cbar_ax.set_yticklabels(['- 1', '- 0.8', '- 0.6', '- 0.4', '- 0.2', '0', '0.2', '0.4', '0.6', '0.8', '1'])
plt.suptitle(P_title, fontsize=20)    
plt.show()
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full   
 
fig.savefig(dir_figs+str(GCM)+'_Correlation_pco2_ALL_cyl.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


P_val_bar=0.05

P_proj='cyl'; P_lon0=210.; P_latN=90.; P_latS=-90.; 

n_r=3 ; n_c=3 ; n_t=9
P_title='Correlations of pCO2 and log Chl, log Biomass, % large Phyto [95% Significance] - '+str(GCM)+' , '+str(start_date_plt)+'-'+str(end_date_plt)
fig=plt.figure()    
for ii in range (0, n_t):
    
    Var_plot_ii=copy.deepcopy(eval(var_names[ii]))
    Var_P_val=copy.deepcopy(eval(var_names_pval[ii]))
    Var_plot_ii [ Var_P_val > P_val_bar ] = nan 
    
    ax = fig.add_subplot(n_r,n_c,ii+1)   
    m = Basemap( projection=P_proj, lon_0=P_lon0, llcrnrlon=P_lon0-180, llcrnrlat=P_latS, urcrnrlon=P_lon0+180, urcrnrlat=P_latN)    
    if ii == (n_c*(n_r-1)): # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        m.drawparallels(np.arange(P_latS, P_latN+0.001, 30.),labels=[True,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Latitutes
        m.drawmeridians(np.arange(0,360,90.),labels=[False,False,False,True], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Longitudes    
    elif ii==0 or ii==n_c or ii==n_c*2 or ii==n_c*3 or ii==n_c*4 or ii==n_c*5 or ii==n_c*6 or ii==n_c*7 or ii==n_c*8:
        m.drawparallels(np.arange(P_latS, P_latN+0.001, 30.),labels=[True,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Latitutes
    elif ii >= n_t-n_c and ii != (n_c*(n_r-1)): # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        m.drawmeridians(np.arange(0,360,90.),labels=[False,False,False,True], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Longitudes
    else:
        m.drawparallels(np.arange(P_latS, P_latN+0.001, 30.),labels=[False,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Latitutes
        m.drawmeridians(np.arange(0,360,90.),labels=[False,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Longitudes
    
    if P_c_fill=='fill':
        m.fillcontinents(color='0')
    m.drawcoastlines(linewidth=1.0, linestyle='solid', color='k', antialiased=1, ax=None, zorder=None)
    im=m.contourf(P_Lon, P_Lat, Var_plot_ii,P_range,latlon=True, cmap=P_cmap) #, extend='both')
    plt.title(plot_titles[ii], fontsize=16)
    
plt.subplots_adjust(left=0.07, bottom=0.1, right=0.92, top=0.92, hspace=0.1, wspace=0.02) # the amount of height/width reserved for space between subplots
cbar_ax = fig.add_axes([0.93, 0.11, 0.015, 0.8]) # [right,bottom,width,height] 
fig.colorbar(im, cax=cbar_ax, ticks=[-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1])
cbar_ax.set_yticklabels(['- 1', '- 0.8', '- 0.6', '- 0.4', '- 0.2', '0', '0.2', '0.4', '0.6', '0.8', '1'])
plt.suptitle(P_title, fontsize=20)    
plt.show()
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full   
 
fig.savefig(dir_figs+str(GCM)+'_Correlation_pco2_ALL_P05_cyl.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


###############################################################################
#############   pCo2 Combination Plots  - Deseasonalized   ####################
###############################################################################
var_names=['Corr_SpCO2_Chl_des','Corr_SpCO2_Phyc_des','Corr_SpCO2_Phydiat_prcnt_des','Corr_SpCO2_TEMP_eq1_Chl_des','Corr_SpCO2_TEMP_eq1_Phyc_des','Corr_SpCO2_TEMP_eq1_Phydiat_prcnt_des','Corr_SpCO2_NonTEMP_eq2_Chl_des','Corr_SpCO2_NonTEMP_eq2_Phyc_des','Corr_SpCO2_NonTEMP_eq2_Phydiat_prcnt_des']
var_names_pval=['Pvalue_SpCO2_Chl_des','Pvalue_SpCO2_Phyc_des','Pvalue_SpCO2_Phydiat_prcnt_des','Pvalue_SpCO2_TEMP_eq1_Chl_des','Pvalue_SpCO2_TEMP_eq1_Phyc_des','Pvalue_SpCO2_TEMP_eq1_Phydiat_prcnt_des','Pvalue_SpCO2_NonTEMP_eq2_Chl_des','Pvalue_SpCO2_NonTEMP_eq2_Phyc_des','Pvalue_SpCO2_NonTEMP_eq2_Phydiat_prcnt_des']
plot_titles=['pCO2 vs. chl','pCO2 vs. biomass','pCO2 vs. % large phyto','pCO2-T vs. chl','pCO2-T vs. biomass','pCO2-T vs. % large phyto','pCO2-nonT vs. chl','pCO2-nonT vs. biomass','pCO2-nonT vs. % large phyto']

P_cmap=plt.cm.bwr; P_c_fill='fill'; P_Lon=Lon_regrid_2D; P_Lat=Lat_regrid_2D; P_range=np.linspace(-0.4,0.4,41);
P_proj='cyl'; P_lon0=210.; P_latN=90.; P_latS=-90.; 

n_r=3 ; n_c=3 ; n_t=9
P_title='Correlations of pCO2 and log Chl, log Biomass, % large Phyto - Deseasonalized - '+str(GCM)+' , '+str(start_date_plt)+'-'+str(end_date_plt)
fig=plt.figure()    
for ii in range (0, n_t):
    
    Var_plot_ii=copy.deepcopy(eval(var_names[ii]))
    
    ax = fig.add_subplot(n_r,n_c,ii+1)   
    m = Basemap( projection=P_proj, lon_0=P_lon0, llcrnrlon=P_lon0-180, llcrnrlat=P_latS, urcrnrlon=P_lon0+180, urcrnrlat=P_latN)    
    if ii == (n_c*(n_r-1)): # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        m.drawparallels(np.arange(P_latS, P_latN+0.001, 30.),labels=[True,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Latitutes
        m.drawmeridians(np.arange(0,360,90.),labels=[False,False,False,True], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Longitudes    
    elif ii==0 or ii==n_c or ii==n_c*2 or ii==n_c*3 or ii==n_c*4 or ii==n_c*5 or ii==n_c*6 or ii==n_c*7 or ii==n_c*8:
        m.drawparallels(np.arange(P_latS, P_latN+0.001, 30.),labels=[True,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Latitutes
    elif ii >= n_t-n_c and ii != (n_c*(n_r-1)): # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        m.drawmeridians(np.arange(0,360,90.),labels=[False,False,False,True], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Longitudes
    else:
        m.drawparallels(np.arange(P_latS, P_latN+0.001, 30.),labels=[False,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Latitutes
        m.drawmeridians(np.arange(0,360,90.),labels=[False,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Longitudes
    
    if P_c_fill=='fill':
        m.fillcontinents(color='0')
    m.drawcoastlines(linewidth=1.0, linestyle='solid', color='k', antialiased=1, ax=None, zorder=None)
    im=m.contourf(P_Lon, P_Lat, Var_plot_ii,P_range,latlon=True, cmap=P_cmap, extend='both')
    plt.title(plot_titles[ii], fontsize=16)
    
plt.subplots_adjust(left=0.07, bottom=0.1, right=0.92, top=0.92, hspace=0.1, wspace=0.02) # the amount of height/width reserved for space between subplots
cbar_ax = fig.add_axes([0.93, 0.11, 0.015, 0.8]) # [right,bottom,width,height] 
fig.colorbar(im, cax=cbar_ax, ticks=[-0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4])
cbar_ax.set_yticklabels(['- 0.4', '- 0.3', '- 0.2', '- 0.1', '0', '0.1', '0.2', '0.3', '0.4'])
plt.suptitle(P_title, fontsize=20)    
plt.show()
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full   
 
fig.savefig(dir_figs+str(GCM)+'_Correlation_pco2_ALL_deseason_cyl.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


P_val_bar=0.05

P_proj='cyl'; P_lon0=210.; P_latN=90.; P_latS=-90.; 

n_r=3 ; n_c=3 ; n_t=9
P_title='Correlations of pCO2 and log Chl, log Biomass, % large Phyto [95% Significance] - Deseasonalized - '+str(GCM)+' , '+str(start_date_plt)+'-'+str(end_date_plt)
fig=plt.figure()    
for ii in range (0, n_t):
    
    Var_plot_ii=copy.deepcopy(eval(var_names[ii]))
    Var_P_val=copy.deepcopy(eval(var_names_pval[ii]))
    Var_plot_ii [ Var_P_val > P_val_bar ] = nan 
    
    ax = fig.add_subplot(n_r,n_c,ii+1)   
    m = Basemap( projection=P_proj, lon_0=P_lon0, llcrnrlon=P_lon0-180, llcrnrlat=P_latS, urcrnrlon=P_lon0+180, urcrnrlat=P_latN)    
    if ii == (n_c*(n_r-1)): # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        m.drawparallels(np.arange(P_latS, P_latN+0.001, 30.),labels=[True,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Latitutes
        m.drawmeridians(np.arange(0,360,90.),labels=[False,False,False,True], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Longitudes    
    elif ii==0 or ii==n_c or ii==n_c*2 or ii==n_c*3 or ii==n_c*4 or ii==n_c*5 or ii==n_c*6 or ii==n_c*7 or ii==n_c*8:
        m.drawparallels(np.arange(P_latS, P_latN+0.001, 30.),labels=[True,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Latitutes
    elif ii >= n_t-n_c and ii != (n_c*(n_r-1)): # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        m.drawmeridians(np.arange(0,360,90.),labels=[False,False,False,True], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Longitudes
    else:
        m.drawparallels(np.arange(P_latS, P_latN+0.001, 30.),labels=[False,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Latitutes
        m.drawmeridians(np.arange(0,360,90.),labels=[False,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Longitudes
    
    if P_c_fill=='fill':
        m.fillcontinents(color='0')
    m.drawcoastlines(linewidth=1.0, linestyle='solid', color='k', antialiased=1, ax=None, zorder=None)
    im=m.contourf(P_Lon, P_Lat, Var_plot_ii,P_range,latlon=True, cmap=P_cmap, extend='both')
    plt.title(plot_titles[ii], fontsize=16)
    
plt.subplots_adjust(left=0.07, bottom=0.1, right=0.92, top=0.92, hspace=0.1, wspace=0.02) # the amount of height/width reserved for space between subplots
cbar_ax = fig.add_axes([0.93, 0.11, 0.015, 0.8]) # [right,bottom,width,height] 
fig.colorbar(im, cax=cbar_ax, ticks=[-0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4])
cbar_ax.set_yticklabels(['- 0.4', '- 0.3', '- 0.2', '- 0.1', '0', '0.1', '0.2', '0.3', '0.4'])
plt.suptitle(P_title, fontsize=20)    
plt.show()
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full   
 
fig.savefig(dir_figs+str(GCM)+'_Correlation_pco2_ALL_deseason_P05_cyl.png', format='png', dpi=300, transparent=True, bbox_inches='tight')




###############################################################################
#############        MLD-SST-SpCO2 Combination Plots       ####################
###############################################################################
#var_names=['Corr_MLD_Chl','Corr_Phyc_MLD','Corr_SpCO2_Chl',    'Corr_MLD_SpCO2','Corr_Phyc_SST','Corr_SpCO2_TEMP_eq1_Chl',       'Corr_MLD_SST','Corr_Phyc_Chl','Corr_SpCO2_NonTEMP_eq2_Chl']
#var_names_pval=['Pvalue_MLD_Chl','Pvalue_Phyc_MLD','Pvalue_SpCO2_Chl',    'Pvalue_MLD_SpCO2','Pvalue_Phyc_SST','Pvalue_SpCO2_TEMP_eq1_Chl',       'Pvalue_MLD_SST','Pvalue_Phyc_Chl','Pvalue_SpCO2_NonTEMP_eq2_Chl']
#plot_titles=['MLD vs. chl','biomass vs. MLD','pCO2 vs. chl',    'MLD vs. pCO2','biomass vs. SST','pCO2-T vs. chl',      'MLD vs. SST','biomass vs. chl','pCO2-nonT vs. chl']
var_names=['Corr_MLD_Chl','Corr_Phyc_MLD','Corr_Phyc_HFDS',    'Corr_MLD_SpCO2','Corr_SpCO2_Phyc','Corr_SpCO2_HFDS',       'Corr_MLD_SST','Corr_Phyc_SST','Corr_HFDS_SST']
var_names_pval=['Pvalue_MLD_Chl','Pvalue_Phyc_MLD','Pvalue_Phyc_HFDS',    'Pvalue_MLD_SpCO2','Pvalue_SpCO2_Phyc','Pvalue_SpCO2_HFDS',       'Pvalue_MLD_SST','Pvalue_Phyc_SST','Pvalue_HFDS_SST']
plot_titles=['MLD vs. chl','biomass vs. MLD','biomass vs. heat flux',    'MLD vs. pCO2','biomass vs. pCO2','heat flux vs. pCO2',      'MLD vs. SST','biomass vs. SST','heat flux vs. SST']

P_cmap=plt.cm.bwr; P_c_fill='fill'; P_Lon=Lon_regrid_2D; P_Lat=Lat_regrid_2D; P_range=np.linspace(-1.,1.,51);
P_proj='cyl'; P_lon0=210.; P_latN=90.; P_latS=-90.; 

n_r=3 ; n_c=3 ; n_t=9
P_title='Monthly Correlations of pCO2, MLD, SST, log Chl, log Biomass - '+str(GCM)+' , '+str(start_date_plt)+'-'+str(end_date_plt)
fig=plt.figure()    
for ii in range (0, n_t):
    
    Var_plot_ii=copy.deepcopy(eval(var_names[ii]))
    
    ax = fig.add_subplot(n_r,n_c,ii+1)   
    m = Basemap( projection=P_proj, lon_0=P_lon0, llcrnrlon=P_lon0-180, llcrnrlat=P_latS, urcrnrlon=P_lon0+180, urcrnrlat=P_latN)    
    if ii == (n_c*(n_r-1)): # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        m.drawparallels(np.arange(P_latS, P_latN+0.001, 30.),labels=[True,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Latitutes
        m.drawmeridians(np.arange(0,360,90.),labels=[False,False,False,True], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Longitudes    
    elif ii==0 or ii==n_c or ii==n_c*2 or ii==n_c*3 or ii==n_c*4 or ii==n_c*5 or ii==n_c*6 or ii==n_c*7 or ii==n_c*8:
        m.drawparallels(np.arange(P_latS, P_latN+0.001, 30.),labels=[True,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Latitutes
    elif ii >= n_t-n_c and ii != (n_c*(n_r-1)): # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        m.drawmeridians(np.arange(0,360,90.),labels=[False,False,False,True], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Longitudes
    else:
        m.drawparallels(np.arange(P_latS, P_latN+0.001, 30.),labels=[False,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Latitutes
        m.drawmeridians(np.arange(0,360,90.),labels=[False,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Longitudes
    
    if P_c_fill=='fill':
        m.fillcontinents(color='0')
    m.drawcoastlines(linewidth=1.0, linestyle='solid', color='k', antialiased=1, ax=None, zorder=None)
    im=m.contourf(P_Lon, P_Lat, Var_plot_ii,P_range,latlon=True, cmap=P_cmap) #, extend='both')
    plt.title(plot_titles[ii], fontsize=16)
    
plt.subplots_adjust(left=0.07, bottom=0.1, right=0.92, top=0.92, hspace=0.1, wspace=0.02) # the amount of height/width reserved for space between subplots
cbar_ax = fig.add_axes([0.93, 0.11, 0.015, 0.8]) # [right,bottom,width,height] 
fig.colorbar(im, cax=cbar_ax, ticks=[-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1])
cbar_ax.set_yticklabels(['- 1', '- 0.8', '- 0.6', '- 0.4', '- 0.2', '0', '0.2', '0.4', '0.6', '0.8', '1'])
plt.suptitle(P_title, fontsize=20)    
plt.show()
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full   
 
fig.savefig(dir_figs+str(GCM)+'_Correlation_MLD_SST_SpCO2_ALL_2.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


P_val_bar=0.05

P_proj='cyl'; P_lon0=210.; P_latN=90.; P_latS=-90.; 

n_r=3 ; n_c=3 ; n_t=9
P_title='Monthly Correlations of pCO2, MLD, SST, log Chl, log Biomass [95% Significance] - '+str(GCM)+' , '+str(start_date_plt)+'-'+str(end_date_plt)
fig=plt.figure()    
for ii in range (0, n_t):
    
    Var_plot_ii=copy.deepcopy(eval(var_names[ii]))
    Var_P_val=copy.deepcopy(eval(var_names_pval[ii]))
    Var_plot_ii [ Var_P_val > P_val_bar ] = nan 
    
    ax = fig.add_subplot(n_r,n_c,ii+1)   
    m = Basemap( projection=P_proj, lon_0=P_lon0, llcrnrlon=P_lon0-180, llcrnrlat=P_latS, urcrnrlon=P_lon0+180, urcrnrlat=P_latN)    
    if ii == (n_c*(n_r-1)): # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        m.drawparallels(np.arange(P_latS, P_latN+0.001, 30.),labels=[True,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Latitutes
        m.drawmeridians(np.arange(0,360,90.),labels=[False,False,False,True], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Longitudes    
    elif ii==0 or ii==n_c or ii==n_c*2 or ii==n_c*3 or ii==n_c*4 or ii==n_c*5 or ii==n_c*6 or ii==n_c*7 or ii==n_c*8:
        m.drawparallels(np.arange(P_latS, P_latN+0.001, 30.),labels=[True,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Latitutes
    elif ii >= n_t-n_c and ii != (n_c*(n_r-1)): # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        m.drawmeridians(np.arange(0,360,90.),labels=[False,False,False,True], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Longitudes
    else:
        m.drawparallels(np.arange(P_latS, P_latN+0.001, 30.),labels=[False,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Latitutes
        m.drawmeridians(np.arange(0,360,90.),labels=[False,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Longitudes
    
    if P_c_fill=='fill':
        m.fillcontinents(color='0')
    m.drawcoastlines(linewidth=1.0, linestyle='solid', color='k', antialiased=1, ax=None, zorder=None)
    im=m.contourf(P_Lon, P_Lat, Var_plot_ii,P_range,latlon=True, cmap=P_cmap) #, extend='both')
    plt.title(plot_titles[ii], fontsize=16)
    
plt.subplots_adjust(left=0.07, bottom=0.1, right=0.92, top=0.92, hspace=0.1, wspace=0.02) # the amount of height/width reserved for space between subplots
cbar_ax = fig.add_axes([0.93, 0.11, 0.015, 0.8]) # [right,bottom,width,height] 
fig.colorbar(im, cax=cbar_ax, ticks=[-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1])
cbar_ax.set_yticklabels(['- 1', '- 0.8', '- 0.6', '- 0.4', '- 0.2', '0', '0.2', '0.4', '0.6', '0.8', '1'])
plt.suptitle(P_title, fontsize=20)    
plt.show()
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full   
 
fig.savefig(dir_figs+str(GCM)+'_Correlation_MLD_SST_SpCO2_ALL_2_P05.png', format='png', dpi=300, transparent=True, bbox_inches='tight')



###############################################################################
########     MLD-SST-SpCO2 Combination Plots  - Deseasonalized     ############
###############################################################################
#var_names=['Corr_MLD_Chl_des','Corr_Phyc_MLD_des','Corr_SpCO2_Chl_des',    'Corr_MLD_SpCO2_des','Corr_Phyc_SST_des','Corr_SpCO2_TEMP_eq1_Chl_des',       'Corr_MLD_SST_des','Corr_Phyc_Chl_des','Corr_SpCO2_NonTEMP_eq2_Chl_des']
#var_names_pval=['Pvalue_MLD_Chl_des','Pvalue_Phyc_MLD_des','Pvalue_SpCO2_Chl_des',    'Pvalue_MLD_SpCO2_des','Pvalue_Phyc_SST_des','Pvalue_SpCO2_TEMP_eq1_Chl_des',       'Pvalue_MLD_SST_des','Pvalue_Phyc_Chl_des','Pvalue_SpCO2_NonTEMP_eq2_Chl_des']
#plot_titles=['MLD vs. chl','biomass vs. MLD','pCO2 vs. chl',    'MLD vs. pCO2','biomass vs. SST','pCO2-T vs. chl',      'MLD vs. SST','biomass vs. chl','pCO2-nonT vs. chl']
var_names=['Corr_MLD_Chl_des','Corr_Phyc_MLD_des','Corr_Phyc_HFDS_des',    'Corr_MLD_SpCO2_des','Corr_SpCO2_Phyc_des','Corr_SpCO2_HFDS_des',       'Corr_MLD_SST_des','Corr_Phyc_SST_des','Corr_HFDS_SST_des']
var_names_pval=['Pvalue_MLD_Chl_des','Pvalue_Phyc_MLD_des','Pvalue_Phyc_HFDS_des',    'Pvalue_MLD_SpCO2_des','Pvalue_SpCO2_Phyc_des','Pvalue_SpCO2_HFDS_des',       'Pvalue_MLD_SST_des','Pvalue_Phyc_SST_des','Pvalue_HFDS_SST_des']
plot_titles=['MLD vs. chl','biomass vs. MLD','biomass vs. heat flux',    'MLD vs. pCO2','biomass vs. pCO2','heat flux vs. pCO2',      'MLD vs. SST','biomass vs. SST','heat flux vs. SST']

P_cmap=plt.cm.bwr; P_c_fill='fill'; P_Lon=Lon_regrid_2D; P_Lat=Lat_regrid_2D; P_range=np.linspace(-0.4,0.4,41);
P_proj='cyl'; P_lon0=210.; P_latN=90.; P_latS=-90.; 

n_r=3 ; n_c=3 ; n_t=9
P_title='Monthly Correlations of pCO2, MLD, SST, log Chl, log Biomass - Deseasonalized - '+str(GCM)+' , '+str(start_date_plt)+'-'+str(end_date_plt)
fig=plt.figure()    
for ii in range (0, n_t):
    
    Var_plot_ii=copy.deepcopy(eval(var_names[ii]))
    
    ax = fig.add_subplot(n_r,n_c,ii+1)   
    m = Basemap( projection=P_proj, lon_0=P_lon0, llcrnrlon=P_lon0-180, llcrnrlat=P_latS, urcrnrlon=P_lon0+180, urcrnrlat=P_latN)    
    if ii == (n_c*(n_r-1)): # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        m.drawparallels(np.arange(P_latS, P_latN+0.001, 30.),labels=[True,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Latitutes
        m.drawmeridians(np.arange(0,360,90.),labels=[False,False,False,True], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Longitudes    
    elif ii==0 or ii==n_c or ii==n_c*2 or ii==n_c*3 or ii==n_c*4 or ii==n_c*5 or ii==n_c*6 or ii==n_c*7 or ii==n_c*8:
        m.drawparallels(np.arange(P_latS, P_latN+0.001, 30.),labels=[True,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Latitutes
    elif ii >= n_t-n_c and ii != (n_c*(n_r-1)): # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        m.drawmeridians(np.arange(0,360,90.),labels=[False,False,False,True], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Longitudes
    else:
        m.drawparallels(np.arange(P_latS, P_latN+0.001, 30.),labels=[False,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Latitutes
        m.drawmeridians(np.arange(0,360,90.),labels=[False,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Longitudes
    
    if P_c_fill=='fill':
        m.fillcontinents(color='0')
    m.drawcoastlines(linewidth=1.0, linestyle='solid', color='k', antialiased=1, ax=None, zorder=None)
    im=m.contourf(P_Lon, P_Lat, Var_plot_ii,P_range,latlon=True, cmap=P_cmap, extend='both')
    plt.title(plot_titles[ii], fontsize=16)
    
plt.subplots_adjust(left=0.07, bottom=0.1, right=0.92, top=0.92, hspace=0.1, wspace=0.02) # the amount of height/width reserved for space between subplots
cbar_ax = fig.add_axes([0.93, 0.11, 0.015, 0.8]) # [right,bottom,width,height] 
fig.colorbar(im, cax=cbar_ax, ticks=[-0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4])
cbar_ax.set_yticklabels(['- 0.4', '- 0.3', '- 0.2', '- 0.1', '0', '0.1', '0.2', '0.3', '0.4'])
plt.suptitle(P_title, fontsize=20)    
plt.show()
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full   
 
fig.savefig(dir_figs+str(GCM)+'_Correlation_MLD_SST_SpCO2_ALL_deseason_2.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


P_val_bar=0.05

P_proj='cyl'; P_lon0=210.; P_latN=90.; P_latS=-90.; 

n_r=3 ; n_c=3 ; n_t=9
P_title='Monthly Correlations of pCO2, MLD, SST, log Chl, log Biomass [95% Significance] - Deseasonalized - '+str(GCM)+' , '+str(start_date_plt)+'-'+str(end_date_plt)
fig=plt.figure()    
for ii in range (0, n_t):
    
    Var_plot_ii=copy.deepcopy(eval(var_names[ii]))
    Var_P_val=copy.deepcopy(eval(var_names_pval[ii]))
    Var_plot_ii [ Var_P_val > P_val_bar ] = nan 
    
    ax = fig.add_subplot(n_r,n_c,ii+1)   
    m = Basemap( projection=P_proj, lon_0=P_lon0, llcrnrlon=P_lon0-180, llcrnrlat=P_latS, urcrnrlon=P_lon0+180, urcrnrlat=P_latN)    
    if ii == (n_c*(n_r-1)): # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        m.drawparallels(np.arange(P_latS, P_latN+0.001, 30.),labels=[True,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Latitutes
        m.drawmeridians(np.arange(0,360,90.),labels=[False,False,False,True], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Longitudes    
    elif ii==0 or ii==n_c or ii==n_c*2 or ii==n_c*3 or ii==n_c*4 or ii==n_c*5 or ii==n_c*6 or ii==n_c*7 or ii==n_c*8:
        m.drawparallels(np.arange(P_latS, P_latN+0.001, 30.),labels=[True,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Latitutes
    elif ii >= n_t-n_c and ii != (n_c*(n_r-1)): # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        m.drawmeridians(np.arange(0,360,90.),labels=[False,False,False,True], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Longitudes
    else:
        m.drawparallels(np.arange(P_latS, P_latN+0.001, 30.),labels=[False,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Latitutes
        m.drawmeridians(np.arange(0,360,90.),labels=[False,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Longitudes
    
    if P_c_fill=='fill':
        m.fillcontinents(color='0')
    m.drawcoastlines(linewidth=1.0, linestyle='solid', color='k', antialiased=1, ax=None, zorder=None)
    im=m.contourf(P_Lon, P_Lat, Var_plot_ii,P_range,latlon=True, cmap=P_cmap, extend='both')
    plt.title(plot_titles[ii], fontsize=16)
    
plt.subplots_adjust(left=0.07, bottom=0.1, right=0.92, top=0.92, hspace=0.1, wspace=0.02) # the amount of height/width reserved for space between subplots
cbar_ax = fig.add_axes([0.93, 0.11, 0.015, 0.8]) # [right,bottom,width,height] 
fig.colorbar(im, cax=cbar_ax, ticks=[-0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4])
cbar_ax.set_yticklabels(['- 0.4', '- 0.3', '- 0.2', '- 0.1', '0', '0.1', '0.2', '0.3', '0.4'])
plt.suptitle(P_title, fontsize=20)    
plt.show()
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full   
 
fig.savefig(dir_figs+str(GCM)+'_Correlation_MLD_SST_SpCO2_ALL_deseason_2_P05.png', format='png', dpi=300, transparent=True, bbox_inches='tight')




















