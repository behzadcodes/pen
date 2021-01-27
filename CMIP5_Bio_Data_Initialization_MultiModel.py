

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
#dir_data_out = (dir_pwd + '/Python_Data_CMIP5/') # Directory to save results
dir_data_out = ('/data5/scratch/Behzad/ocean/ocean_biogeochemistry/processed_data/python_data_CMIP5_bio_1991_2010/')

### Regridding calculations ###
# creating new coordinate grid, same which was used in interpolation in data processing code
lat_n_regrid, lon_n_regrid = 180, 360 # Number of Lat and Lon elements in the regridded data
lon_min_regrid, lon_max_regrid = 0, 360 # Min and Max value of Lon in the regridded data
lat_min_regrid, lat_max_regrid = -90, 90 # Min and Max value of Lat in the regridded data

Lat_regrid_1D, Lon_regrid_1D, Lat_bound_regrid, Lon_bound_regrid = func_latlon_regrid(lat_n_regrid, lon_n_regrid, lat_min_regrid, lat_max_regrid, lon_min_regrid, lon_max_regrid)
Lon_regrid_2D, Lat_regrid_2D = np.meshgrid(Lon_regrid_1D, Lat_regrid_1D)

Ocean_Land_mask = func_oceanlandmask(Lat_regrid_2D, Lon_regrid_2D) # 1= ocean, 0= land

Ocean_Index = func_oceanindex (Lat_regrid_2D, Lon_regrid_2D) # [0=land] [2=Pacific] [3=Indian Ocean] [6=Atlantic] [10=Arctic] [8=Baffin Bay (west of Greenland)] [9=Norwegian Sea (east of Greenland)] [11=Hudson Bay (Canada)] 
###############################################################################

### Variables to be edited by user ###
start_date_plt=1991 # Start Date used for plotting
end_date_plt=2010 # End Date used for plotting

start_date_cal_h=1991 # Start Date for Calculations # Define as string, not number
end_date_cal_h=2005 # End Date for Calculations
start_date_cal_r=2006 # Start Date for Calculations # Define as string, not number
end_date_cal_r=2010 # End Date for Calculations

yrs_n=20;

###############################################################################

### 'GISS-E2-H-CC' and 'GISS-E2-R-CC' have chl data only until 2010
### GISS MODELS have problems with SST data   -   GISS MODELS has problem with NO3 data
### 'MPI-ESM-MR' # has problems with EPC100_allmonths  PPint_allmonths, data are overlaping
### MRI-ESM1 has aconflicting historical and esmHistorical data in it's folder - same for rcp85 and esmrcp85

GCM_Names = ['CESM1-BGC' ,'GFDL-ESM2G','GFDL-ESM2M','IPSL-CM5A-MR','IPSL-CM5A-LR','GISS-E2-H-CC','GISS-E2-R-CC','MRI-ESM1','HadGEM2-ES','MPI-ESM-MR','NorESM1-ME','CanESM2','MPI-ESM-LR']
#GCM_Names = ['GISS-E2-H-CC','GISS-E2-R-CC']

for M_i in range(len(GCM_Names)): # GCM=GCM_Names[0] # M_i=0
    
    GCM=GCM_Names[M_i]
    print (GCM)
    
    dir_data_in_h=(dir_data_in1+ GCM + '/historical/mo/') # Directory to read Historical data
    dir_data_in_r=(dir_data_in1+ GCM + '/rcp85/mo/') # Directory to read RCP8.5 data

    ############################
    ### Phyc (P1+P2) reading ###
    variable='phyc' # Units: mol/m3 # mole_concentration_of_phytoplankton_expressed_as_carbon_in_sea_water
    dset1= netcdf_read(dir_data_in_h+variable+'_Omon_'+str(GCM)+'_historical'+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
    
    start_date_i,end_date_i = dset1.find_time(dset1.times, start_date_cal_h, end_date_cal_h) # This finds which row in the time variable corresponds to the desired start_date_cal and end_date_cal for this specific calculation
    Lat_orig=dset1.y
    Lon_orig=dset1.x    
    dset1.close_ncfile(dset1.fin)
    
    dset2 = MFDataset(dir_data_in_h+variable+'_Omon_'+str(GCM)+'_historical'+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
    Var_dset = dset2.variables[variable] # append phyc data to variable
    
    Data_allmonths_orig=np.asarray(Var_dset[start_date_i:end_date_i+1,:,:]) # phyc(time, rlat, rlon) - phydiat of all the desired months at the surface
    Data_allmonths_orig[Data_allmonths_orig <= 0] = nan
    dset2.close()
    #### Adding RCP data to the Historical data ####
    dset1= netcdf_read(dir_data_in_r+variable+'_Omon_'+str(GCM)+'_rcp85'+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
    start_date_i,end_date_i = dset1.find_time(dset1.times, start_date_cal_r, end_date_cal_r)
    dset1.close_ncfile(dset1.fin)
    dset2 = MFDataset(dir_data_in_r+variable+'_Omon_'+str(GCM)+'_rcp85'+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
    Var_dset = dset2.variables[variable]
    Data_allmonths_orig_r=np.asarray(Var_dset[start_date_i:end_date_i+1,:,:]) # phyc(time, rlat, rlon) - phydiat of all the desired months at the surface
    Data_allmonths_orig_r[Data_allmonths_orig <= 0] = nan
    dset2.close()
    Data_allmonths_orig = np.concatenate((Data_allmonths_orig, Data_allmonths_orig_r), axis=0)
    ### Regriding data to 180*360 degrees
    Phyc_allmonths = func_regrid(Data_allmonths_orig, Lat_orig, Lon_orig, Lat_regrid_2D, Lon_regrid_2D) # Phyc at surface - all months 
    Phyc_allmonths[Phyc_allmonths > 1E19] = nan # Replacing missing values with NaN - Missing values in the original .nc file are filled with 1E+20 values

    if  GCM !='IPSL-CM5A-MR' and GCM !='IPSL-CM5A-LR' and GCM !='IPSL-CM5B-LR' and GCM !='MRI-ESM1' and GCM !='MPI-ESM-MR' and GCM !='NorESM1-ME' and GCM !='CanESM2' and GCM !='MPI-ESM-LR':  # These models don't have phyico data
        #####################################
        ### Phypico (small phyto) reading ###
        variable='phypico' # Units: mol/m3 # mole_concentration_of_diatoms_expressed_as_carbon_in_sea_water
        dset1= netcdf_read(dir_data_in_h+variable+'_Omon_'+str(GCM)+'_historical'+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
        ## dset1.variables  # Shows the variables in the .nc file
        ## This finds which row in the time variable corresponds to the desired start_date_cal and end_date_cal for this specific calculation
        start_date_i,end_date_i = dset1.find_time(dset1.times, start_date_cal_h, end_date_cal_h)
        dset1.close_ncfile(dset1.fin)
        
        dset2 = MFDataset(dir_data_in_h+variable+'_Omon_'+str(GCM)+'_historical'+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
            
        Var_dset = dset2.variables[variable] # append phyico data to variable
        
        Data_allmonths_orig=np.asarray(Var_dset[start_date_i:end_date_i+1,:,:]) # phyico(time, rlat, rlon) - phyico of all the desired months at the surface
        Data_allmonths_orig[Data_allmonths_orig <= 0] = nan
        dset2.close()
        #### Adding RCP data to the Historical data ####
        if GCM=='GFDL-ESM2M':
            dset2 = Dataset(dir_data_in_r+'phypico_Omon_GFDL-ESM2M_historical_r1i1p1_200101-200512.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
            Var_dset = dset2.variables[variable]
            Data_allmonths_orig_r=np.asarray(Var_dset[:,:,:]) # phyico(time, rlat, rlon) - phydiat of all the desired months at the surface
        else:
            dset1= netcdf_read(dir_data_in_r+variable+'_Omon_'+str(GCM)+'_rcp85'+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
            start_date_i,end_date_i = dset1.find_time(dset1.times, start_date_cal_r, end_date_cal_r)
            dset1.close_ncfile(dset1.fin)
            dset2 = MFDataset(dir_data_in_r+variable+'_Omon_'+str(GCM)+'_rcp85'+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
            Var_dset = dset2.variables[variable]
            Data_allmonths_orig_r=np.asarray(Var_dset[start_date_i:end_date_i+1,:,:]) # phyico(time, rlat, rlon) - phydiat of all the desired months at the surface
        Data_allmonths_orig_r[Data_allmonths_orig <= 0] = nan
        dset2.close()
        Data_allmonths_orig = np.concatenate((Data_allmonths_orig, Data_allmonths_orig_r), axis=0)
        ### Regriding data to 180*360 degrees
        Phypico_allmonths = func_regrid(Data_allmonths_orig, Lat_orig, Lon_orig, Lat_regrid_2D, Lon_regrid_2D) # Phypico at surface - all months 
        Phypico_allmonths[Phypico_allmonths > 1E19] = nan # Replacing missing values with NaN - Missing values in the original .nc file are filled with 1E+20 values


    if  GCM !='MRI-ESM1' and GCM !='MPI-ESM-MR' and GCM !='NorESM1-ME' and GCM !='CanESM2' and GCM !='MPI-ESM-LR':  # These models don't have phydiat data
        ############################
        ### Phydiat (P2) reading ###
        variable='phydiat' # Units: mol/m3 # mole_concentration_of_diatoms_expressed_as_carbon_in_sea_water
        dset1= netcdf_read(dir_data_in_h+variable+'_Omon_'+str(GCM)+'_historical'+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
        ## dset1.variables  # Shows the variables in the .nc file
        ## This finds which row in the time variable corresponds to the desired start_date_cal and end_date_cal for this specific calculation
        start_date_i,end_date_i = dset1.find_time(dset1.times, start_date_cal_h, end_date_cal_h)
        dset1.close_ncfile(dset1.fin)
        
        dset2 = MFDataset(dir_data_in_h+variable+'_Omon_'+str(GCM)+'_historical'+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
            
        Var_dset = dset2.variables[variable] # append phydiat data to variable
        
        Data_allmonths_orig=np.asarray(Var_dset[start_date_i:end_date_i+1,:,:]) # phydiat(time, rlat, rlon) - phydiat of all the desired months at the surface
        Data_allmonths_orig[Data_allmonths_orig <= 0] = nan
        dset2.close()
        #### Adding RCP data to the Historical data ####
        if GCM=='GFDL-ESM2M': 
            dir_data_in_r_p=(dir_data_in3+ GCM + '/rcp85/mo/') # Directory to read RCP8.5 data
        else:
            dir_data_in_r_p=(dir_data_in1+ GCM + '/rcp85/mo/') # Directory to read RCP8.5 data  

        dset1= netcdf_read(dir_data_in_r_p+variable+'_Omon_'+str(GCM)+'_rcp85'+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
        start_date_i,end_date_i = dset1.find_time(dset1.times, start_date_cal_r, end_date_cal_r)
        dset1.close_ncfile(dset1.fin)
        dset2 = MFDataset(dir_data_in_r_p+variable+'_Omon_'+str(GCM)+'_rcp85'+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
        Var_dset = dset2.variables[variable]
        Data_allmonths_orig_r=np.asarray(Var_dset[start_date_i:end_date_i+1,:,:]) # phydiat(time, rlat, rlon) - phydiat of all the desired months at the surface
        Data_allmonths_orig_r[Data_allmonths_orig <= 0] = nan
        dset2.close()
        Data_allmonths_orig = np.concatenate((Data_allmonths_orig, Data_allmonths_orig_r), axis=0)
        ### Regriding data to 180*360 degrees
        Phydiat_allmonths = func_regrid(Data_allmonths_orig, Lat_orig, Lon_orig, Lat_regrid_2D, Lon_regrid_2D) # Phydiat at surface - all months 
        Phydiat_allmonths[Phydiat_allmonths > 1E19] = nan # Replacing missing values with NaN - Missing values in the original .nc file are filled with 1E+20 values

        ########################################
        ### Phydiat (Large Phyto) percentage ###
        #Phydiat_prcnt_allmonths = (Phydiat_allmonths / Phyc_allmonths) * 100
        ### Small phyto ######
        #Physmall_allmonths = Phyc_allmonths - Phydiat_allmonths # Units: mol/m3   

    if GCM !='NorESM1-ME':
        ############################
        ### Zooplankton reading  ###
        variable='zooc' # Units: mol/m3 # mole_concentration_of_zooplankton_expressed_as_carbon_in_sea_water
        
        if GCM=='IPSL-CM5A-MR' or GCM=='IPSL-CM5A-LR' or  GCM=='MPI-ESM-LR':
            dir_data_in_h_z=(dir_data_in3+ GCM + '/historical/mo/') # Directory to read RCP8.5 data
        else:
            dir_data_in_h_z=(dir_data_in1+ GCM + '/historical/mo/') # Directory to read RCP8.5 data      
         
        dset1= netcdf_read(dir_data_in_h_z+variable+'_Omon_'+str(GCM)+'_historical'+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
        
        start_date_i,end_date_i = dset1.find_time(dset1.times, start_date_cal_h, end_date_cal_h) # This finds which row in the time variable corresponds to the desired start_date_cal and end_date_cal for this specific calculation
        dset1.close_ncfile(dset1.fin)
        
        dset2 = MFDataset(dir_data_in_h_z+variable+'_Omon_'+str(GCM)+'_historical'+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
        Var_dset = dset2.variables[variable] # append phyc data to variable
        
        Data_allmonths_orig=np.asarray(Var_dset[start_date_i:end_date_i+1,:,:]) # Zooc(time, rlat, rlon) - phydiat of all the desired months at the surface
        Data_allmonths_orig[Data_allmonths_orig <= 0] = nan
        dset2.close()
        #### Adding RCP data to the Historical data ####
        dset1= netcdf_read(dir_data_in_r+variable+'_Omon_'+str(GCM)+'_rcp85'+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
        start_date_i,end_date_i = dset1.find_time(dset1.times, start_date_cal_r, end_date_cal_r)
        dset1.close_ncfile(dset1.fin)
        dset2 = MFDataset(dir_data_in_r+variable+'_Omon_'+str(GCM)+'_rcp85'+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
        Var_dset = dset2.variables[variable]
        Data_allmonths_orig_r=np.asarray(Var_dset[start_date_i:end_date_i+1,:,:]) # Zooc(time, rlat, rlon) - phydiat of all the desired months at the surface
        Data_allmonths_orig_r[Data_allmonths_orig <= 0] = nan
        dset2.close()
        Data_allmonths_orig = np.concatenate((Data_allmonths_orig, Data_allmonths_orig_r), axis=0)
        ### Regriding data to 180*360 degrees
        Zooc_allmonths = func_regrid(Data_allmonths_orig, Lat_orig, Lon_orig, Lat_regrid_2D, Lon_regrid_2D) # Zooc at surface - all months 
        Zooc_allmonths[Zooc_allmonths > 1E19] = nan # Replacing missing values with NaN - Missing values in the original .nc file are filled with 1E+20 values    
    
    #################################
    ### Primary Production reading ###
    variable='intpp' # Units: mol m-2 s-1 # net_primary_mole_productivity_of_carbon_by_phytoplankton - Vertically integrated total primary (organic carbon) production by phytoplankton
    dset1= netcdf_read(dir_data_in_h+variable+'_Omon_'+str(GCM)+'_historical'+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
    
    start_date_i,end_date_i = dset1.find_time(dset1.times, start_date_cal_h, end_date_cal_h) # This finds which row in the time variable corresponds to the desired start_date_cal and end_date_cal for this specific calculation
    dset1.close_ncfile(dset1.fin)
    
    dset2 = MFDataset(dir_data_in_h+variable+'_Omon_'+str(GCM)+'_historical'+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
    Var_dset = dset2.variables[variable] # append phyc data to variable
    
    Data_allmonths_orig=np.asarray(Var_dset[start_date_i:end_date_i+1,:,:]) # intpp(time, rlat, rlon) - phydiat of all the desired months at the surface
    Data_allmonths_orig[Data_allmonths_orig <= 0] = nan
    dset2.close()
    #### Adding RCP data to the Historical data ####
    dset1= netcdf_read(dir_data_in_r+variable+'_Omon_'+str(GCM)+'_rcp85'+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
    start_date_i,end_date_i = dset1.find_time(dset1.times, start_date_cal_r, end_date_cal_r)
    dset1.close_ncfile(dset1.fin)
    dset2 = MFDataset(dir_data_in_r+variable+'_Omon_'+str(GCM)+'_rcp85'+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
    Var_dset = dset2.variables[variable]
    Data_allmonths_orig_r=np.asarray(Var_dset[start_date_i:end_date_i+1,:,:]) # intpp(time, rlat, rlon) - phydiat of all the desired months at the surface
    Data_allmonths_orig_r[Data_allmonths_orig <= 0] = nan
    dset2.close()
    Data_allmonths_orig = np.concatenate((Data_allmonths_orig, Data_allmonths_orig_r), axis=0)
    ### Regriding data to 180*360 degrees
    PPint_allmonths = func_regrid(Data_allmonths_orig, Lat_orig, Lon_orig, Lat_regrid_2D, Lon_regrid_2D) # intpp at surface - all months 
    PPint_allmonths[PPint_allmonths > 1E19] = nan # Replacing missing values with NaN - Missing values in the original .nc file are filled with 1E+20 values 


    if GCM !='GISS-E2-H-CC' and GCM !='GISS-E2-R-CC': # These models don't have epc100 data
        #########################################
        ### Export Production at 100m reading ###
        variable='epc100' # Units: mol m-2 s-1 # Downward Flux of Particle Organic Carbon - at 100 m depth
        dset1= netcdf_read(dir_data_in_h+variable+'_Omon_'+str(GCM)+'_historical'+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
        
        start_date_i,end_date_i = dset1.find_time(dset1.times, start_date_cal_h, end_date_cal_h) # This finds which row in the time variable corresponds to the desired start_date_cal and end_date_cal for this specific calculation
        dset1.close_ncfile(dset1.fin)
        
        dset2 = MFDataset(dir_data_in_h+variable+'_Omon_'+str(GCM)+'_historical'+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
        Var_dset = dset2.variables[variable] # append phyc data to variable
        
        Data_allmonths_orig=np.asarray(Var_dset[start_date_i:end_date_i+1,:,:]) # intpp(time, rlat, rlon) - phydiat of all the desired months at the surface
        Data_allmonths_orig[Data_allmonths_orig <= 0] = nan
        dset2.close()
        #### Adding RCP data to the Historical data ####
        dset1= netcdf_read(dir_data_in_r+variable+'_Omon_'+str(GCM)+'_rcp85'+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
        start_date_i,end_date_i = dset1.find_time(dset1.times, start_date_cal_r, end_date_cal_r)
        dset1.close_ncfile(dset1.fin)
        dset2 = MFDataset(dir_data_in_r+variable+'_Omon_'+str(GCM)+'_rcp85'+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
        Var_dset = dset2.variables[variable]
        Data_allmonths_orig_r=np.asarray(Var_dset[start_date_i:end_date_i+1,:,:]) # intpp(time, rlat, rlon) - phydiat of all the desired months at the surface
        Data_allmonths_orig_r[Data_allmonths_orig <= 0] = nan
        dset2.close()
        Data_allmonths_orig = np.concatenate((Data_allmonths_orig, Data_allmonths_orig_r), axis=0)
        ### Regriding data to 180*360 degrees
        EPC100_allmonths = func_regrid(Data_allmonths_orig, Lat_orig, Lon_orig, Lat_regrid_2D, Lon_regrid_2D) # intpp at surface - all months 
        EPC100_allmonths[EPC100_allmonths > 1E19] = nan # Replacing missing values with NaN - Missing values in the original .nc file are filled with 1E+20 values        
        

    if GCM !='NorESM1-ME' and GCM !='MPI-ESM-LR':      
        if GCM =='GISS-E2-H-CC' or GCM=='GISS-E2-R-CC': # for these two models, the historical data folder for chl has dataup to 2010
            ############################
            ### Chlorophyll reading ###
            variable='chl' # Units: kg m-3 # mass_concentration_of_phytoplankton_expressed_as_chlorophyll_in_sea_water
            dset1= netcdf_read(dir_data_in_h+variable+'_Omon_'+str(GCM)+'_historical'+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
            
            start_date_i,end_date_i = dset1.find_time(dset1.times, start_date_cal_h, end_date_cal_r) # This finds which row in the time variable corresponds to the desired start_date_cal and end_date_cal for this specific calculation
            dset1.close_ncfile(dset1.fin)
            
            dset2 = MFDataset(dir_data_in_h+variable+'_Omon_'+str(GCM)+'_historical'+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
            Var_dset = dset2.variables[variable] # append phyc data to variable
            
            Data_allmonths_orig=np.asarray(Var_dset[start_date_i:end_date_i+1,:,:]) # chl(time, rlat, rlon) - phydiat of all the desired months at the surface
            Data_allmonths_orig[Data_allmonths_orig <= 0] = nan
            dset2.close()
    
            ### Regriding data to 180*360 degrees
            Chl_allmonths = func_regrid(Data_allmonths_orig, Lat_orig, Lon_orig, Lat_regrid_2D, Lon_regrid_2D) # Chl at surface - all months 
            Chl_allmonths[Chl_allmonths > 1E19] = nan # Replacing missing values with NaN - Missing values in the original .nc file are filled with 1E+20 values        
            
        else:
            ############################
            ### Chlorophyll reading ###
            variable='chl' # Units: kg m-3 # mass_concentration_of_phytoplankton_expressed_as_chlorophyll_in_sea_water
            dset1= netcdf_read(dir_data_in_h+variable+'_Omon_'+str(GCM)+'_historical'+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
            
            start_date_i,end_date_i = dset1.find_time(dset1.times, start_date_cal_h, end_date_cal_h) # This finds which row in the time variable corresponds to the desired start_date_cal and end_date_cal for this specific calculation
            dset1.close_ncfile(dset1.fin)
            
            dset2 = MFDataset(dir_data_in_h+variable+'_Omon_'+str(GCM)+'_historical'+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
            Var_dset = dset2.variables[variable] # append phyc data to variable
            
            Data_allmonths_orig=np.asarray(Var_dset[start_date_i:end_date_i+1,:,:]) # chl(time, rlat, rlon) - phydiat of all the desired months at the surface
            Data_allmonths_orig[Data_allmonths_orig <= 0] = nan
            dset2.close()
            #### Adding RCP data to the Historical data ####
            dset1= netcdf_read(dir_data_in_r+variable+'_Omon_'+str(GCM)+'_rcp85'+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
            start_date_i,end_date_i = dset1.find_time(dset1.times, start_date_cal_r, end_date_cal_r)
            dset1.close_ncfile(dset1.fin)
            dset2 = MFDataset(dir_data_in_r+variable+'_Omon_'+str(GCM)+'_rcp85'+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
            Var_dset = dset2.variables[variable]
            Data_allmonths_orig_r=np.asarray(Var_dset[start_date_i:end_date_i+1,:,:]) # chl(time, rlat, rlon) - phydiat of all the desired months at the surface
            Data_allmonths_orig_r[Data_allmonths_orig <= 0] = nan
            dset2.close()
            Data_allmonths_orig = np.concatenate((Data_allmonths_orig, Data_allmonths_orig_r), axis=0)
            ### Regriding data to 180*360 degrees
            Chl_allmonths = func_regrid(Data_allmonths_orig, Lat_orig, Lon_orig, Lat_regrid_2D, Lon_regrid_2D) # Chl at surface - all months 
            Chl_allmonths[Chl_allmonths > 1E19] = nan # Replacing missing values with NaN - Missing values in the original .nc file are filled with 1E+20 values
    
    ############################
    ###  Nitrate reading     ###
    variable='no3' # Units: mol m-3 # mole_concentration_of_nitrate_in_sea_water
    dset1= netcdf_read(dir_data_in_h+variable+'_Omon_'+str(GCM)+'_historical'+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
    
    start_date_i,end_date_i = dset1.find_time(dset1.times, start_date_cal_h, end_date_cal_h) # This finds which row in the time variable corresponds to the desired start_date_cal and end_date_cal for this specific calculation
    dset1.close_ncfile(dset1.fin)
    
    dset2 = MFDataset(dir_data_in_h+variable+'_Omon_'+str(GCM)+'_historical'+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
    Var_dset = dset2.variables[variable] # append phyc data to variable
    
    Data_allmonths_orig=np.asarray(Var_dset[start_date_i:end_date_i+1,:,:]) # NO3(time, rlat, rlon) - phydiat of all the desired months at the surface
    Data_allmonths_orig[Data_allmonths_orig <= 0] = nan
    dset2.close()
    #### Adding RCP data to the Historical data ####
    dset1= netcdf_read(dir_data_in_r+variable+'_Omon_'+str(GCM)+'_rcp85'+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
    start_date_i,end_date_i = dset1.find_time(dset1.times, start_date_cal_r, end_date_cal_r)
    dset1.close_ncfile(dset1.fin)
    dset2 = MFDataset(dir_data_in_r+variable+'_Omon_'+str(GCM)+'_rcp85'+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
    Var_dset = dset2.variables[variable]
    Data_allmonths_orig_r=np.asarray(Var_dset[start_date_i:end_date_i+1,:,:]) # NO3(time, rlat, rlon) - phydiat of all the desired months at the surface
    Data_allmonths_orig_r[Data_allmonths_orig <= 0] = nan
    dset2.close()
    Data_allmonths_orig = np.concatenate((Data_allmonths_orig, Data_allmonths_orig_r), axis=0)
    ### Regriding data to 180*360 degrees
    NO3_allmonths = func_regrid(Data_allmonths_orig, Lat_orig, Lon_orig, Lat_regrid_2D, Lon_regrid_2D) # NO3 at surface - all months 
    NO3_allmonths[NO3_allmonths > 1E19] = nan # Replacing missing values with NaN - Missing values in the original .nc file are filled with 1E+20 values

    if GCM !='MRI-ESM1' and GCM !='CanESM2': # These models don't have dfe and Si data
        ############################
        ###   Iron reading       ### 
        variable='dfe' # Units: mol m-3 # dissolved iron in sea water is meant to include both Fe2+ and Fe3+ ions (but not, e.g., particulate detrital iron)
        dset1= netcdf_read(dir_data_in_h+variable+'_Omon_'+str(GCM)+'_historical'+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
        
        start_date_i,end_date_i = dset1.find_time(dset1.times, start_date_cal_h, end_date_cal_h) # This finds which row in the time variable corresponds to the desired start_date_cal and end_date_cal for this specific calculation
        dset1.close_ncfile(dset1.fin)
        
        dset2 = MFDataset(dir_data_in_h+variable+'_Omon_'+str(GCM)+'_historical'+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
        Var_dset = dset2.variables[variable] # append phyc data to variable
        
        Data_allmonths_orig=np.asarray(Var_dset[start_date_i:end_date_i+1,:,:]) # dfe(time, rlat, rlon) - phydiat of all the desired months at the surface
        Data_allmonths_orig[Data_allmonths_orig <= 0] = nan
        dset2.close()
        #### Adding RCP data to the Historical data ####
        dset1= netcdf_read(dir_data_in_r+variable+'_Omon_'+str(GCM)+'_rcp85'+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
        start_date_i,end_date_i = dset1.find_time(dset1.times, start_date_cal_r, end_date_cal_r)
        dset1.close_ncfile(dset1.fin)
        dset2 = MFDataset(dir_data_in_r+variable+'_Omon_'+str(GCM)+'_rcp85'+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
        Var_dset = dset2.variables[variable]
        Data_allmonths_orig_r=np.asarray(Var_dset[start_date_i:end_date_i+1,:,:]) # dfe(time, rlat, rlon) - phydiat of all the desired months at the surface
        Data_allmonths_orig_r[Data_allmonths_orig <= 0] = nan
        dset2.close()
        Data_allmonths_orig = np.concatenate((Data_allmonths_orig, Data_allmonths_orig_r), axis=0)
        ### Regriding data to 180*360 degrees
        Fe_allmonths = func_regrid(Data_allmonths_orig, Lat_orig, Lon_orig, Lat_regrid_2D, Lon_regrid_2D) # dfe at surface - all months 
        Fe_allmonths[Fe_allmonths > 1E19] = nan # Replacing missing values with NaN - Missing values in the original .nc file are filled with 1E+20 values

        if  GCM!='GFDL-ESM2M':
            ############################
            ###   Silicate reading   ### 
            variable='si' # Units: mol m-3 # Dissolved Silicate Concentration at Surface
            dset1= netcdf_read(dir_data_in_h+variable+'_Omon_'+str(GCM)+'_historical'+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
            
            start_date_i,end_date_i = dset1.find_time(dset1.times, start_date_cal_h, end_date_cal_h) # This finds which row in the time variable corresponds to the desired start_date_cal and end_date_cal for this specific calculation
            dset1.close_ncfile(dset1.fin)
            
            dset2 = MFDataset(dir_data_in_h+variable+'_Omon_'+str(GCM)+'_historical'+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
            Var_dset = dset2.variables[variable] # append phyc data to variable
            
            Data_allmonths_orig=np.asarray(Var_dset[start_date_i:end_date_i+1,:,:]) # Si(time, rlat, rlon) - phydiat of all the desired months at the surface
            Data_allmonths_orig[Data_allmonths_orig <= 0] = nan
            dset2.close()
            #### Adding RCP data to the Historical data ####
            dset1= netcdf_read(dir_data_in_r+variable+'_Omon_'+str(GCM)+'_rcp85'+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
            start_date_i,end_date_i = dset1.find_time(dset1.times, start_date_cal_r, end_date_cal_r)
            dset1.close_ncfile(dset1.fin)
            dset2 = MFDataset(dir_data_in_r+variable+'_Omon_'+str(GCM)+'_rcp85'+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
            Var_dset = dset2.variables[variable]
            Data_allmonths_orig_r=np.asarray(Var_dset[start_date_i:end_date_i+1,:,:]) # Si(time, rlat, rlon) - phydiat of all the desired months at the surface
            Data_allmonths_orig_r[Data_allmonths_orig <= 0] = nan
            dset2.close()
            Data_allmonths_orig = np.concatenate((Data_allmonths_orig, Data_allmonths_orig_r), axis=0)
            ### Regriding data to 180*360 degrees
            Si_allmonths = func_regrid(Data_allmonths_orig, Lat_orig, Lon_orig, Lat_regrid_2D, Lon_regrid_2D) # Si at surface - all months 
            Si_allmonths[Si_allmonths > 1E19] = nan # Replacing missing values with NaN - Missing values in the original .nc file are filled with 1E+20 values    

    if GCM !='IPSL-CM5A-MR' and GCM !='IPSL-CM5A-LR' and GCM !='IPSL-CM5B-LR' : # These models don't have SpCO2 data
        ############################
        ###     SpCO2 reading    ###
        variable='spco2' # Units: Pa # surface_partial_pressure_of_carbon_dioxide_in_sea_water
        dset1= netcdf_read(dir_data_in_h+variable+'_Omon_'+str(GCM)+'_historical'+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
        
        start_date_i,end_date_i = dset1.find_time(dset1.times, start_date_cal_h, end_date_cal_h) # This finds which row in the time variable corresponds to the desired start_date_cal and end_date_cal for this specific calculation
        dset1.close_ncfile(dset1.fin)
        
        dset2 = MFDataset(dir_data_in_h+variable+'_Omon_'+str(GCM)+'_historical'+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
        Var_dset = dset2.variables[variable] # append phyc data to variable
        
        Data_allmonths_orig=np.asarray(Var_dset[start_date_i:end_date_i+1,:,:]) # spco2(time, rlat, rlon) - phydiat of all the desired months at the surface
        Data_allmonths_orig[Data_allmonths_orig <= 0] = nan
        dset2.close()
        #### Adding RCP data to the Historical data ####       
        if GCM=='GFDL-ESM2G' or GCM=='GFDL-ESM2M' or  GCM=='MIROC-ESM' or GCM=='MIROC-ESM-CHEM': 
            dir_data_in_r_s=(dir_data_in3+ GCM + '/rcp85/mo/') # Directory to read RCP8.5 data
        else:
            dir_data_in_r_s=(dir_data_in1+ GCM + '/rcp85/mo/') # Directory to read RCP8.5 data  
    
        dset1= netcdf_read(dir_data_in_r_s+variable+'_Omon_'+str(GCM)+'_rcp85'+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
        start_date_i,end_date_i = dset1.find_time(dset1.times, start_date_cal_r, end_date_cal_r)
        dset1.close_ncfile(dset1.fin)
        dset2 = MFDataset(dir_data_in_r_s+variable+'_Omon_'+str(GCM)+'_rcp85'+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
        Var_dset = dset2.variables[variable]
        Data_allmonths_orig_r=np.asarray(Var_dset[start_date_i:end_date_i+1,:,:]) # spco2(time, rlat, rlon) - phydiat of all the desired months at the surface
        Data_allmonths_orig_r[Data_allmonths_orig <= 0] = nan
        dset2.close()
        Data_allmonths_orig = np.concatenate((Data_allmonths_orig, Data_allmonths_orig_r), axis=0)
        ### Regriding data to 180*360 degrees
        SpCO2_allmonths = func_regrid(Data_allmonths_orig, Lat_orig, Lon_orig, Lat_regrid_2D, Lon_regrid_2D) # SpCO2 at surface - all months 
        SpCO2_allmonths[SpCO2_allmonths > 1E19] = nan # Replacing missing values with NaN - Missing values in the original .nc file are filled with 1E+20 values
        SpCO2_allmonths = SpCO2_allmonths / 0.101325 # To Convert spco2 unit from Pa to PPM
        
        if GCM=='NorESM1-ME':
            SpCO2_allmonths=SpCO2_allmonths/1e6 # This model has a weird conversion factor apparently    
        elif GCM=='MPI-ESM-LR':
            SpCO2_allmonths=SpCO2_allmonths/10
    
    
    ####################################
    ###   Air-sea CO2 flux reading   ### 
    variable='fgco2' # Units: kg m-2 s-1 # Surface Downward CO2 Flux - Gas exchange flux of CO2 (positive into ocean)
    if GCM=='MRI-ESM1': # This model does not have rcp85 for fgco2, it has esmrcp8.5, so we use esmHistorical as well
        dset1= netcdf_read(dir_data_in_h+variable+'_Omon_'+str(GCM)+'_esmHistorical'+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
    else:
        dset1= netcdf_read(dir_data_in_h+variable+'_Omon_'+str(GCM)+'_historical'+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
    
    start_date_i,end_date_i = dset1.find_time(dset1.times, start_date_cal_h, end_date_cal_h) # This finds which row in the time variable corresponds to the desired start_date_cal and end_date_cal for this specific calculation
    dset1.close_ncfile(dset1.fin)

    if GCM=='MRI-ESM1':
        dset2 = MFDataset(dir_data_in_h+variable+'_Omon_'+str(GCM)+'_esmHistorical'+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
    else:
        dset2 = MFDataset(dir_data_in_h+variable+'_Omon_'+str(GCM)+'_historical'+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
    Var_dset = dset2.variables[variable] # append phyc data to variable
    
    Data_allmonths_orig=np.asarray(Var_dset[start_date_i:end_date_i+1,:,:]) # fgco2(time, rlat, rlon) - phydiat of all the desired months at the surface
    dset2.close()
    #### Adding RCP data to the Historical data ####
    if GCM=='MRI-ESM1':    
        dset1= netcdf_read(dir_data_in_r+variable+'_Omon_'+str(GCM)+'_esmrcp85'+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
    else:
        dset1= netcdf_read(dir_data_in_r+variable+'_Omon_'+str(GCM)+'_rcp85'+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
    start_date_i,end_date_i = dset1.find_time(dset1.times, start_date_cal_r, end_date_cal_r)
    dset1.close_ncfile(dset1.fin)
    
    if GCM=='MRI-ESM1':    
        dset2 = MFDataset(dir_data_in_r+variable+'_Omon_'+str(GCM)+'_esmrcp85'+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
    else:
        dset2 = MFDataset(dir_data_in_r+variable+'_Omon_'+str(GCM)+'_rcp85'+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
    Var_dset = dset2.variables[variable]
    Data_allmonths_orig_r=np.asarray(Var_dset[start_date_i:end_date_i+1,:,:]) # fgco2(time, rlat, rlon) - phydiat of all the desired months at the surface
    dset2.close()
    Data_allmonths_orig = np.concatenate((Data_allmonths_orig, Data_allmonths_orig_r), axis=0)
    ### Regriding data to 180*360 degrees
    FgCO2_allmonths = func_regrid(Data_allmonths_orig, Lat_orig, Lon_orig, Lat_regrid_2D, Lon_regrid_2D) # fgco2 at surface - all months 
    FgCO2_allmonths[FgCO2_allmonths > 1E19] = nan # Replacing missing values with NaN - Missing values in the original .nc file are filled with 1E+20 values       
    FgCO2_allmonths = FgCO2_allmonths * 3600 * 24 * 365 * 1000 # To convert from second to year and from kg to g => g/m2/year of co2
    FgCO2_allmonths = FgCO2_allmonths / (44/12) # to convert from gram of co2 to gram of C - 44ra CO2 == 12gr C - Final unit: gram C per m2 per year
    # Final unit: gram C per m2 per year #
    
    
    ###############################################################################
    ###########                 Saving Results                #####################    
    import os
    import shelve 

    filename_out = (dir_data_out + 'AllResults_'+GCM+'_bio_'+str(start_date_plt)+'_'+str(end_date_plt)+'.out') # Directory to save processed data
    
    if GCM=='IPSL-CM5A-MR' or GCM=='IPSL-CM5A-LR' or GCM=='IPSL-CM5B-LR': # These models don't have SpCO2 data
        var_list_1=['Phydiat_allmonths','Phyc_allmonths','Zooc_allmonths','PPint_allmonths','EPC100_allmonths','Chl_allmonths','FgCO2_allmonths','Fe_allmonths','NO3_allmonths','Si_allmonths','Lat_regrid_1D','Lon_regrid_1D','Lat_bound_regrid','Lon_bound_regrid','Lon_regrid_2D','Lat_regrid_2D']    
    elif GCM=='MPI-ESM-MR':  # These models don't have phydiat data
        var_list_1=['Phyc_allmonths','Zooc_allmonths','PPint_allmonths','EPC100_allmonths','Chl_allmonths','FgCO2_allmonths','SpCO2_allmonths','Fe_allmonths','NO3_allmonths','Si_allmonths','Lat_regrid_1D','Lon_regrid_1D','Lat_bound_regrid','Lon_bound_regrid','Lon_regrid_2D','Lat_regrid_2D']        
    elif GCM=='GFDL-ESM2M':  # These models don't have Si data
        var_list_1=['Phydiat_allmonths','Phyc_allmonths','Phypico_allmonths','Zooc_allmonths','PPint_allmonths','EPC100_allmonths','Chl_allmonths','FgCO2_allmonths','SpCO2_allmonths','Fe_allmonths','NO3_allmonths','Lat_regrid_1D','Lon_regrid_1D','Lat_bound_regrid','Lon_bound_regrid','Lon_regrid_2D','Lat_regrid_2D']        
    elif GCM=='NorESM1-ME' or GCM=='MPI-ESM-LR':  # These models don't have phydiat data and Zooc data and Chl data
        var_list_1=['Phyc_allmonths','PPint_allmonths','EPC100_allmonths','FgCO2_allmonths','SpCO2_allmonths','Fe_allmonths','NO3_allmonths','Si_allmonths','Lat_regrid_1D','Lon_regrid_1D','Lat_bound_regrid','Lon_bound_regrid','Lon_regrid_2D','Lat_regrid_2D']        
    elif GCM=='MRI-ESM1' or GCM=='CanESM2':  # These models don't have phydiat data and dfe and Si data
        var_list_1=['Phyc_allmonths','Zooc_allmonths','PPint_allmonths','EPC100_allmonths','Chl_allmonths','FgCO2_allmonths','SpCO2_allmonths','NO3_allmonths','Lat_regrid_1D','Lon_regrid_1D','Lat_bound_regrid','Lon_bound_regrid','Lon_regrid_2D','Lat_regrid_2D']        

    elif GCM =='GISS-E2-H-CC' or GCM=='GISS-E2-R-CC':  # These models don't have EPC100 data
        var_list_1=['Phydiat_allmonths','Phyc_allmonths','Phypico_allmonths','Zooc_allmonths','PPint_allmonths','Chl_allmonths','FgCO2_allmonths','SpCO2_allmonths','Fe_allmonths','NO3_allmonths','Si_allmonths','Lat_regrid_1D','Lon_regrid_1D','Lat_bound_regrid','Lon_bound_regrid','Lon_regrid_2D','Lat_regrid_2D']        
    else:
        var_list_1=['Phydiat_allmonths','Phyc_allmonths','Phypico_allmonths','Zooc_allmonths','PPint_allmonths','EPC100_allmonths','Chl_allmonths','FgCO2_allmonths','SpCO2_allmonths','Fe_allmonths','NO3_allmonths','Si_allmonths','Lat_regrid_1D','Lon_regrid_1D','Lat_bound_regrid','Lon_bound_regrid','Lon_regrid_2D','Lat_regrid_2D']    
    
    ### To save
    my_shelf = shelve.open(filename_out,'n') # 'n' for new
    
    for key in var_list_1:
        try:
            my_shelf[key] = globals()[key]
        except TypeError:
            print('ERROR shelving: {0}'.format(key))
    my_shelf.close()
   


    
#GCM_Names = ['HadGEM2-ES'] # Reading historical thetao problem
#GCM_Names = ['MPI-ESM-MR'] # Reading rcp so problem
GCM_Names = ['HadGEM2-ES','MPI-ESM-MR']

for M_i in range(len(GCM_Names)): # GCM=GCM_Names[0]
    
    GCM=GCM_Names[M_i]
    print (GCM)
    
    dir_data_in_h_t=(dir_data_in2+ GCM + '/historical/mo/') # Directory to read Historical data
    dir_data_in_r_t=(dir_data_in2+ GCM + '/rcp85/mo/') # Directory to read RCP8.5 data

    #############################################
    ### Sea Surface Temperature (SST) reading ###
    variable='thetao'
    
    if GCM=='HadGEM2-ES' or GCM=='MIROC-ESM': # Some models have only concated data which the file name ends with concat.nc, instead of ending in 12.nc
        dset1= netcdf_read(dir_data_in_h_t+variable+'_Omon_'+str(GCM)+'_historical'+'*concat.nc',str(GCM), 'THETAO')# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
    else:
        dset1= netcdf_read(dir_data_in_h_t+variable+'_Omon_'+str(GCM)+'_historical'+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
    
    start_date_i,end_date_i = dset1.find_time(dset1.times, start_date_cal_h, end_date_cal_h) # This finds which row in the time variable corresponds to the desired start_date_cal and end_date_cal for this specific calculation
    Lat_orig=dset1.y
    Lon_orig=dset1.x
    dset1.close_ncfile(dset1.fin)
    
    if GCM=='HadGEM2-ES' or GCM=='MIROC-ESM': # Some models have only concated data which the file name ends with concat.nc, instead of ending in 12.nc
        dset2 = MFDataset(dir_data_in_h_t+variable+'_Omon_'+str(GCM)+'_historical'+'*concat.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
        Var_dset = dset2.variables['THETAO'] # append data to variable
    else:
        dset2 = MFDataset(dir_data_in_h_t+variable+'_Omon_'+str(GCM)+'_historical'+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
        Var_dset = dset2.variables[variable] # append data to variable

    
    Data_allmonths_orig=np.asarray(Var_dset[start_date_i:end_date_i+1,0,:,:]) # thetao(time, lev, rlat, rlon) - phydiat of all the desired years, only on the surface
    Data_allmonths_orig[Data_allmonths_orig <= 0] = nan
    dset2.close()
    #### Adding RCP data to the Historical data ####
    if GCM=='HadGEM2-ES': # This models data names end in 11.nc
        dset1= netcdf_read(dir_data_in_r_t+variable+'_Omon_'+str(GCM)+'_rcp85'+'*11.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
    else:
        dset1= netcdf_read(dir_data_in_r_t+variable+'_Omon_'+str(GCM)+'_rcp85'+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
    start_date_i,end_date_i = dset1.find_time(dset1.times, start_date_cal_r, end_date_cal_r)
    dset1.close_ncfile(dset1.fin)
    
    if GCM=='HadGEM2-ES': # This models data names end in 11.nc
        dset2 = MFDataset(dir_data_in_r_t+variable+'_Omon_'+str(GCM)+'_rcp85'+'*11.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
    else:
        dset2 = MFDataset(dir_data_in_r_t+variable+'_Omon_'+str(GCM)+'_rcp85'+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
    Var_dset = dset2.variables[variable]
    Data_allmonths_orig_r=np.asarray(Var_dset[start_date_i:end_date_i+1,0,:,:]) # thetao(time, lev, rlat, rlon) - phydiat of all the desired years, only on the surface
    Data_allmonths_orig_r[Data_allmonths_orig <= 0] = nan
    dset2.close()
    Data_allmonths_orig = np.concatenate((Data_allmonths_orig, Data_allmonths_orig_r), axis=0)
    ### Regriding data to 180*360 degrees
    SST_allmonths = func_regrid(Data_allmonths_orig, Lat_orig, Lon_orig, Lat_regrid_2D, Lon_regrid_2D) # thetao at surface - all months 
    SST_allmonths[SST_allmonths > 1E19] = nan # Replacing missing values with NaN - Missing values in the original .nc file are filled with 1E+20 values
    SST_allmonths = SST_allmonths - 273.15 # To convert Temperature unit from K to C

    ###########################
    ##     MLD reading     ###
    from Behzadlib2 import func_MLD_Allmonths_bio
    
    variable_thetao='thetao' # Sea Water Temperature
    variable_so='so' # Sea Water Salinity

    if GCM=='HadGEM2-ES' or GCM=='MIROC-ESM': # Some models have only concated data which the file name ends with concat.nc, instead of ending in 12.nc
        dset_thetao= netcdf_read(dir_data_in_h_t+variable_thetao+'_Omon_'+str(GCM)+'_historical'+'*concat.nc',str(GCM), variable_thetao)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
        dset_so= netcdf_read(dir_data_in_h_t+variable_so+'_Omon_'+str(GCM)+'_historical'+'*concat.nc',str(GCM), variable_so)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
    else:
        dset_thetao= netcdf_read(dir_data_in_h_t+variable_thetao+'_Omon_'+str(GCM)+'_historical'+'*12.nc',str(GCM), variable_thetao)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
        dset_so= netcdf_read(dir_data_in_h_t+variable_so+'_Omon_'+str(GCM)+'_historical'+'*12.nc',str(GCM), variable_so)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
    
    Data_allmonths_orig = func_MLD_Allmonths_bio(dset_thetao, dset_so, start_date_cal_h, end_date_cal_h, lat_n_regrid, lon_n_regrid)
    Data_allmonths_orig[Data_allmonths_orig <= 0] = nan
    #### Adding RCP data to the Historical data ####
    if GCM=='HadGEM2-ES': # This models thetao data names end in 11.nc and SO data end in concat.nc
        dset_thetao= netcdf_read(dir_data_in_r_t+variable_thetao+'_Omon_'+str(GCM)+'_rcp85'+'*11.nc',str(GCM), variable_thetao)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
        dset_so= netcdf_read(dir_data_in_r_t+variable_so+'_Omon_'+str(GCM)+'_rcp85'+'*concat.nc',str(GCM), variable_so)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read        
    else:
        dset_thetao= netcdf_read(dir_data_in_r_t+variable_thetao+'_Omon_'+str(GCM)+'_rcp85'+'*12.nc',str(GCM), variable_thetao)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
        dset_so= netcdf_read(dir_data_in_r_t+variable_so+'_Omon_'+str(GCM)+'_rcp85'+'*12.nc',str(GCM), variable_so)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read        

    Data_allmonths_orig_r = func_MLD_Allmonths_bio(dset_thetao, dset_so, start_date_cal_r, end_date_cal_r, lat_n_regrid, lon_n_regrid)
    Data_allmonths_orig_r[Data_allmonths_orig_r <= 0] = nan
    Data_allmonths_orig = np.concatenate((Data_allmonths_orig, Data_allmonths_orig_r), axis=0)
    
    MLD_allmonths = copy.deepcopy(Data_allmonths_orig) # MLD - all months 
    MLD_allmonths[MLD_allmonths > 1E19] = nan # Replacing missing values with NaN - Missing values in the original .nc file are filled with 1E+20 values
    np.save('MLD_allmonths_'+GCM+'_'+str(start_date_plt)+'_'+str(end_date_plt)+'.npy',MLD_allmonths)   
    
#    dir_data_in_m=('/data4/CMIP5/CMIP5_models/descriptors/ocean_physics/') # Directory to read Historical data
#    variable='MLD'
#    dset1 = Dataset(dir_data_in_m+'mld_Omon_'+str(GCM)+'_historicalrcp85_r1i1p1_186501-210012.nc') 
#    
#    Var_dset = dset1.variables[variable] # append phyc data to variable
#    Data_allmonths_orig=np.asarray(dset1.variables[variable][(start_date_plt-1865)*12:(end_date_plt-1865+1)*12,:,:])
#    
#    MLD_allmonths  = func_regrid(Data_allmonths_orig, Lat_orig, Lon_orig, Lat_regrid_2D, Lon_regrid_2D) # SpCO2 at surface - all months 
#    MLD_allmonths[MLD_allmonths > 1E19] = nan # Replacing missing values with NaN - Missing values in the original .nc file are filled with 1E+20 values
#    MLD_allmonths[MLD_allmonths < -1E19] = nan # Replacing missing values with NaN - Missing values in the original .nc file are filled with 1E+20 values    
    
    if GCM=='GFDL-ESM2G' or GCM=='GISS-E2-R' or GCM=='MIROC-ESM' or GCM=='MIROC-ESM-CHEM' or GCM=='MPI-ESM-LR' or GCM=='MPI-ESM-MR' or GCM=='MPI-ESM-P' or GCM=='MRI-CGCM3' or GCM=='NorESM1-ME' or GCM=='CNRM-CM5':   
        ################################
        ###     Heat Flux reading    ###
        variable='hfds' # Units: W m-2 # Downward Heat Flux at Sea Water Surface
        
        if GCM=='GFDL-ESM2G':
            dir_data_in_h_t=(dir_data_in4+ GCM + '/historical/mo/') # Directory to read Historical data
            dir_data_in_r_t=(dir_data_in4+ GCM + '/rcp85/mo/') # Directory to read RCP8.5 data
        elif GCM=='NorESM1-ME':
            dir_data_in_h_t=(dir_data_in4+ 'NorESM1-M/NorESM1-ME/') # Directory to read Historical data
            dir_data_in_r_t=(dir_data_in4+ 'NorESM1-M/NorESM1-ME/') # Directory to read RCP8.5 data            
        else:
            dir_data_in_h_t=(dir_data_in4+ GCM + '/') # Directory to read Historical data
            dir_data_in_r_t=(dir_data_in4+ GCM + '/') # Directory to read RCP8.5 data 
        
        dset1= netcdf_read(dir_data_in_h_t+variable+'_Omon_'+str(GCM)+'_historical'+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
        
        start_date_i,end_date_i = dset1.find_time(dset1.times, start_date_cal_h, end_date_cal_h) # This finds which row in the time variable corresponds to the desired start_date_cal and end_date_cal for this specific calculation
        dset1.close_ncfile(dset1.fin)
        
        dset2 = MFDataset(dir_data_in_h_t+variable+'_Omon_'+str(GCM)+'_historical'+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
        Var_dset = dset2.variables[variable] # append phyc data to variable
        
        Data_allmonths_orig=np.asarray(Var_dset[start_date_i:end_date_i+1,:,:]) # hfds(time, rlat, rlon) - phydiat of all the desired years, only on the surface
        Data_allmonths_orig[Data_allmonths_orig <= 0] = nan
        dset2.close()
        #### Adding RCP data to the Historical data ####
        dset1= netcdf_read(dir_data_in_r_t+variable+'_Omon_'+str(GCM)+'_rcp85'+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
        start_date_i,end_date_i = dset1.find_time(dset1.times, start_date_cal_r, end_date_cal_r)
        dset1.close_ncfile(dset1.fin)
        dset2 = MFDataset(dir_data_in_r_t+variable+'_Omon_'+str(GCM)+'_rcp85'+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
        Var_dset = dset2.variables[variable]
        Data_allmonths_orig_r=np.asarray(Var_dset[start_date_i:end_date_i+1,:,:]) # hfds(time, rlat, rlon) - phydiat of all the desired years, only on the surface
        Data_allmonths_orig_r[Data_allmonths_orig <= 0] = nan
        dset2.close()
        Data_allmonths_orig = np.concatenate((Data_allmonths_orig, Data_allmonths_orig_r), axis=0)
        ### Regriding data to 180*360 degrees
        HFDS_allmonths = func_regrid(Data_allmonths_orig, Lat_orig, Lon_orig, Lat_regrid_2D, Lon_regrid_2D) # HFDS at surface - all months 
        HFDS_allmonths[HFDS_allmonths > 1E19] = nan # Replacing missing values with NaN - Missing values in the original .nc file are filled with 1E+20 values


    if GCM !='GISS-E2-H-CC' and GCM !='GISS-E2-R-CC' and GCM!='HadGEM2-ES': # These models don't have rsntds data
        ######################################################################
        ### Net Downward Shortwave Radiation at Sea Water Surface  reading ###
        variable='rsntds' # units: W m-2
        
        if GCM=='HadGEM2-ES' or GCM=='MIROC-ESM': # Some models have only concated data which the file name ends with concat.nc, instead of ending in 12.nc
            dset1= netcdf_read(dir_data_in_h_t+variable+'_Omon_'+str(GCM)+'_historical'+'*concat.nc',str(GCM), 'THETAO')# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
        elif GCM=='MRI-ESM1':
            dset1= netcdf_read(dir_data_in_h_t+variable+'_Omon_'+str(GCM)+'_esmHistorical'+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
        else:
            dset1= netcdf_read(dir_data_in_h_t+variable+'_Omon_'+str(GCM)+'_historical'+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
        
        start_date_i,end_date_i = dset1.find_time(dset1.times, start_date_cal_h, end_date_cal_h) # This finds which row in the time variable corresponds to the desired start_date_cal and end_date_cal for this specific calculation
        Lat_orig=dset1.y
        Lon_orig=dset1.x
        dset1.close_ncfile(dset1.fin)
        
        if GCM=='HadGEM2-ES' or GCM=='MIROC-ESM': # Some models have only concated data which the file name ends with concat.nc, instead of ending in 12.nc
            dset2 = MFDataset(dir_data_in_h_t+variable+'_Omon_'+str(GCM)+'_historical'+'*concat.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
        elif GCM=='MRI-ESM1':
            dset2 = MFDataset(dir_data_in_h_t+variable+'_Omon_'+str(GCM)+'_esmHistorical'+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read              
        else:
            dset2 = MFDataset(dir_data_in_h_t+variable+'_Omon_'+str(GCM)+'_historical'+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
        Var_dset = dset2.variables[variable] # append data to variable       
        Data_allmonths_orig=np.asarray(Var_dset[start_date_i:end_date_i+1,:,:]) # thetao(time, lev, rlat, rlon) - phydiat of all the desired years, only on the surface
        dset2.close()
        #### Adding RCP data to the Historical data ####
        if GCM=='HadGEM2-ES': # This models data names end in 11.nc
            dset1= netcdf_read(dir_data_in_r_t+variable+'_Omon_'+str(GCM)+'_rcp85'+'*11.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
        elif GCM=='MRI-ESM1':
            dset1= netcdf_read(dir_data_in_r_t+variable+'_Omon_'+str(GCM)+'_esmrcp85'+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
        else:
            dset1= netcdf_read(dir_data_in_r_t+variable+'_Omon_'+str(GCM)+'_rcp85'+'*12.nc',str(GCM), variable)# Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
        start_date_i,end_date_i = dset1.find_time(dset1.times, start_date_cal_r, end_date_cal_r)
        dset1.close_ncfile(dset1.fin)
        
        if GCM=='HadGEM2-ES': # This models data names end in 11.nc
            dset2 = MFDataset(dir_data_in_r_t+variable+'_Omon_'+str(GCM)+'_rcp85'+'*11.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
        elif GCM=='MRI-ESM1':
            dset2 = MFDataset(dir_data_in_r_t+variable+'_Omon_'+str(GCM)+'_esmrcp85'+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
        else:
            dset2 = MFDataset(dir_data_in_r_t+variable+'_Omon_'+str(GCM)+'_rcp85'+'*12.nc') # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
        Var_dset = dset2.variables[variable]
        Data_allmonths_orig_r=np.asarray(Var_dset[start_date_i:end_date_i+1,:,:]) # thetao(time, lev, rlat, rlon) - phydiat of all the desired years, only on the surface
        dset2.close()
        Data_allmonths_orig = np.concatenate((Data_allmonths_orig, Data_allmonths_orig_r), axis=0)
        ### Regriding data to 180*360 degrees
        RSNTDS_allmonths = func_regrid(Data_allmonths_orig, Lat_orig, Lon_orig, Lat_regrid_2D, Lon_regrid_2D) # thetao at surface - all months 
        if GCM=='CESM1-BGC': # for this model the values are negative
            RSNTDS_allmonths = RSNTDS_allmonths *-1 # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
    
        RSNTDS_allmonths[RSNTDS_allmonths > 1E19] = nan # Replacing missing values with NaN - Missing values in the original .nc file are filled with 1E+20 values

    ###############################################################################
    ###########                 Saving Results                #####################    
    import os
    import shelve 

    filename_out = (dir_data_out + 'AllResults_'+GCM+'_physics_'+str(start_date_plt)+'_'+str(end_date_plt)+'.out') # Directory to save processed data
    if GCM=='GFDL-ESM2G' or GCM=='GISS-E2-R' or GCM=='MIROC-ESM' or GCM=='MIROC-ESM-CHEM' or GCM=='MPI-ESM-LR' or GCM=='MPI-ESM-MR' or GCM=='MPI-ESM-P' or GCM=='MRI-CGCM3' or GCM=='NorESM1-ME' or GCM=='CNRM-CM5':  # These models have hfds data  
        var_list_1=['SST_allmonths','MLD_allmonths','HFDS_allmonths','RSNTDS_allmonths','Lat_regrid_1D','Lon_regrid_1D','Lat_bound_regrid','Lon_bound_regrid','Lon_regrid_2D','Lat_regrid_2D']  
    elif GCM =='GISS-E2-H-CC' or GCM =='GISS-E2-R-CC' or GCM=='HadGEM2-ES': # These models don't have rsntds data
        var_list_1=['SST_allmonths','MLD_allmonths','Lat_regrid_1D','Lon_regrid_1D','Lat_bound_regrid','Lon_bound_regrid','Lon_regrid_2D','Lat_regrid_2D']  
    else:
        var_list_1=['SST_allmonths','MLD_allmonths','RSNTDS_allmonths','Lat_regrid_1D','Lon_regrid_1D','Lat_bound_regrid','Lon_bound_regrid','Lon_regrid_2D','Lat_regrid_2D']  
    ### To save
    my_shelf = shelve.open(filename_out,'n') # 'n' for new
    
    for key in var_list_1:
        try:
            my_shelf[key] = globals()[key]
        except TypeError:
            print('ERROR shelving: {0}'.format(key))
    my_shelf.close()









