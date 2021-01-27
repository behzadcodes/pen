### IMPORTANT: The following path needs to be appended so that the neccessary codes can be read from that directory ###
import sys
sys.path.append('/data1/home/basadieh/behzadcodes/behzadlibrary')
#from Behzadlib import * ### This is loaded from '/data1/home/basadieh/behzadcodes/behzadlibrary'
#from CMIP5lib import *  ### This is loaded from '/data1/home/basadieh/behzadcodes/behzadlibrary'
from Behzadlib import func_MLD, func_time_depth_plot, func_stream, func_ENSO, func_NAO, func_latlon_regrid, func_latlon_regrid_eq, func_oceanlandmask, func_oceanindex, func_regrid, func_plot, func_plot_save, func_plot_bounds_save
from BehzadlibPlot import func_plotmap_contourf, func_plot2Dcolor_contourf, func_plotline_1var, func_plotline_2var, func_plot2Dcontour_1var, func_plot2Dcontour_2var, func_plot_lagcor_sig, func_plot_laggedmaps
from CMIP5lib import netcdf_read, runningMeanFast, plot_PSD_welch, lag_cor, lag_cor_data, calc_Atl_Mask
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
dir_figs = (dir_pwd + '/Figures/') # Directory to save figures
dir_results = (dir_pwd + '/Results/') # Directory to save results

lat_n_regrid, lon_n_regrid = 180, 360 # Number of Lat and Lon elements in the regridded data
lon_min_regrid, lon_max_regrid = 0, 360 # Min and Max value of Lon in the regridded data
lat_min_regrid, lat_max_regrid = -90, 90 # Min and Max value of Lat in the regridded data

Lat_regrid_1D, Lon_regrid_1D, Lat_bound_regrid, Lon_bound_regrid = func_latlon_regrid(lat_n_regrid, lon_n_regrid, lat_min_regrid, lat_max_regrid, lon_min_regrid, lon_max_regrid)
Lon_regrid_2D, Lat_regrid_2D = np.meshgrid(Lon_regrid_1D, Lat_regrid_1D)

Ocean_Land_mask = func_oceanlandmask(Lat_regrid_2D, Lon_regrid_2D) # 1= ocean, 0= land
Ocean_Index = func_oceanindex (Lat_regrid_2D, Lon_regrid_2D) # [0=land] [2=Pacific] [3=Indian Ocean] [6=Atlantic] [10=Arctic] [8=Baffin Bay (west of Greenland)] [9=Norwegian Sea (east of Greenland)] [11=Hudson Bay (Canada)] 
##############################################################################
from matplotlib.ticker import MultipleLocator
dir_pwd=os.getcwd() # Gets the current directory (and in which the code is placed)
dir_data_in1 = ('/data6/CMIP6/DCPP') # Directory to raed raw .nc data from

GCM_Institution_ID = ['CCCma', 'MIROC', 'NCC' ]
GCM_Source_ID = ['CanESM5',  'MIROC6', 'NorCPM1' ]

GCM_runs=[ ['r1i1p2f1', 'r2i1p2f1', 'r3i1p2f1', 'r4i1p2f1', 'r5i1p2f1', 'r6i1p2f1', 'r7i1p2f1', 'r8i1p2f1', 'r9i1p2f1', 'r10i1p2f1'], 
         ['r1i1p1f1', 'r2i1p1f1', 'r3i1p1f1', 'r4i1p1f1', 'r5i1p1f1', 'r6i1p1f1', 'r7i1p1f1', 'r8i1p1f1', 'r9i1p1f1', 'r10i1p1f1'],
         ['r1i1p1f1', 'r2i1p1f1', 'r3i1p1f1', 'r4i1p1f1', 'r5i1p1f1', 'r6i1p1f1', 'r7i1p1f1', 'r8i1p1f1', 'r9i1p1f1', 'r10i1p1f1'] ]

Var_names = ['pr', 'tas', 'mrros', 'sos', 'tos']
Var_AOL= ['A', 'A', 'L', 'O', 'O'] # A = Atmo var, O = Ocean var, L = Land var, etc
Var_names_long = ['Precipitation [mm/day]', 'Surface Air Temperature [C]', 'Surface Runoff', 'Sea Surf. Salinity', 'Sea Surf. Temperature']

###############################################################################
################   Reading ERA5 reanalysis data   #############################
dset_era5 = Dataset('/data6/Observations/reanalysis/ERA5/ERA5_TotalPrecipitation_1979_2019_Monthly.nc')
ERA5_PRCP_1979_2018_monthly=np.asarray(dset_era5.variables['tp'][:,:,:]) 
ERA5_PRCP_1979_2018_monthly=ERA5_PRCP_1979_2018_monthly[0:480,:,:] * 1000 ## *1000 is to convert from m/day units to mm/day
ERA5_PRCP_1979_2018_annual=ERA5_PRCP_1979_2018_monthly.reshape(40,12,ERA5_PRCP_1979_2018_monthly.shape[1],ERA5_PRCP_1979_2018_monthly.shape[2])
ERA5_PRCP_1979_2018_annual=np.nanmean(ERA5_PRCP_1979_2018_annual, axis=1)

dset_era5 = Dataset('/data6/Observations/reanalysis/ERA5/ERA5_SurfaceAirTemperature_1979_2019_Monthly.nc')
ERA5_TEMP_1979_2018_monthly=np.asarray(dset_era5.variables['t2m'][:,:,:]) 
ERA5_TEMP_1979_2018_monthly=ERA5_TEMP_1979_2018_monthly[0:480,:,:] - 273.15 ## - 273.15 is to convert from K units to C
ERA5_TEMP_1979_2018_annual=ERA5_TEMP_1979_2018_monthly.reshape(40,12,ERA5_TEMP_1979_2018_monthly.shape[1],ERA5_TEMP_1979_2018_monthly.shape[2])
ERA5_TEMP_1979_2018_annual=np.nanmean(ERA5_TEMP_1979_2018_annual, axis=1)
dset_era5.close

Lat_orig_era5=np.asarray(dset_era5.variables['latitude'][:])
Lon_orig_era5=np.asarray(dset_era5.variables['longitude'][:])

###############################################################################
################   Reading WDFEI reanalysis data   #############################
Days_M=[31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31] # Number of days in each month
Days_M_long=np.tile(Days_M,10)

WFDEI_PRCP_1961_2016_monthly=[]
for yr in range(1961,2015,10): # yr=1961 # ii=0
    
    if yr==2011:
        dset_wfdei = Dataset('/data6/Observations/reanalysis/WDFEI/pr/pr_gpcc_wfdei_{}_{}.nc'.format(yr,yr+5))
    else:
        dset_wfdei = Dataset('/data6/Observations/reanalysis/WDFEI/pr/pr_gpcc_wfdei_{}_{}.nc'.format(yr,yr+9))
    data_daily=np.asarray(dset_wfdei.variables['pr'][:,:,:]) * 24 * 3600; # Converting prcp from kg/m2/s to kg/m2/day == mm/day   (kg/m2 = mm)
    data_daily=data_daily[:-2,:,:] # Omitting the last 2 days - Simplification to avoid leap year calculations
    dset_wfdei.close

    cc=0
    for ii in range(0,np.int(data_daily.shape[0]/365)*12): # Number of months in the daily time series       
        WFDEI_PRCP_1961_2016_monthly.append(np.nanmean( data_daily[cc:cc+Days_M_long[ii],:,:],axis=0 ))
        cc=cc+Days_M_long[ii]
    print('Year {}-{} - appended'.format(yr,yr+9))

WFDEI_PRCP_1961_2016_monthly=np.asarray(WFDEI_PRCP_1961_2016_monthly)
WFDEI_PRCP_1961_2016_monthly [WFDEI_PRCP_1961_2016_monthly>1e19]=nan
###########################################
WFDEI_TEMP_1961_2016_monthly=[]
for yr in range(1961,2015,10): # yr=1961 # ii=0
    
    if yr==2011:
        dset_wfdei = Dataset('/data6/Observations/reanalysis/WDFEI/pr/tas_gpcc_wfdei_{}_{}.nc'.format(yr,yr+5))
    else:
        dset_wfdei = Dataset('/data6/Observations/reanalysis/WDFEI/pr/tas_gpcc_wfdei_{}_{}.nc'.format(yr,yr+9))
    data_daily=np.asarray(dset_wfdei.variables['pr'][:,:,:]) * 24 * 3600; # Converting prcp from kg/m2/s to kg/m2/day == mm/day   (kg/m2 = mm)
    data_daily=data_daily[:-2,:,:] # Omitting the last 2 days - Simplification to avoid leap year calculations
    dset_wfdei.close

    cc=0
    for ii in range(0,np.int(data_daily.shape[0]/365)*12): # Number of months in the daily time series       
        WFDEI_TEMP_1961_2016_monthly.append(np.nanmean( data_daily[cc:cc+Days_M_long[ii],:,:],axis=0 ))
        cc=cc+Days_M_long[ii]
    print('Year {}-{} - appended'.format(yr,yr+9))

WFDEI_TEMP_1961_2016_monthly=np.asarray(WFDEI_TEMP_1961_2016_monthly)
WFDEI_TEMP_1961_2016_monthly [WFDEI_TEMP_1961_2016_monthly>1e19]=nan

###############################################################################          
for M_i in range(len(GCM_Institution_ID)): # M_i=2
   
    GCM_id1 = GCM_Institution_ID[M_i]
    GCM_id2 = GCM_Source_ID[M_i]
    
    for vv in range(0, len(Var_names)): # vv=0
        Var_names_vv=Var_names[vv]
        Var_AOL_vv=Var_AOL[vv]        
        
        r_b_count=0 # a counter to alternate betwwen red and blur line colors            
        fig=plt.figure()

        for yr in range(1960,2019): # yr=1960 
            
            if yr%5 == 0: # Only take the years with 5yr increments, 1960,965,1970,...
                
                Plot_X = np.linspace(yr+1,yr+10,10)
                Plot_var=[]
                for rr in range(0,len(GCM_runs[M_i])): # rr=0 # ensemble number
                    
                    file_dir = dir_data_in1+'/'+GCM_id1+'/'+'dcppA-hindcast'+'/'+GCM_runs[M_i][rr]+'/'+'mon'+'/'+Var_names_vv # Final path where the files will be saved

                    if GCM_id1 == 'CCCma':
                        if os.path.isfile(file_dir+'/'+Var_names_vv+'_'+Var_AOL_vv+'mon_'+GCM_id2+'_dcppA-hindcast_s'+str(yr)+'-'+GCM_runs[M_i][rr]+'_gn_'+str(yr+1)+'01-'+str(yr+10)+'12.nc'):# If this file exists in this directory
                            file_dir_name=file_dir+'/'+Var_names_vv+'_'+Var_AOL_vv+'mon_'+GCM_id2+'_dcppA-hindcast_s'+str(yr)+'-'+GCM_runs[M_i][rr]+'_gn_'+str(yr+1)+'01-'+str(yr+10)+'12.nc'
                            
                    elif GCM_id1 == 'NCC':
                        if os.path.isfile(file_dir+'/'+Var_names_vv+'_'+Var_AOL_vv+'mon_'+GCM_id2+'_dcppA-hindcast_s'+str(yr)+'-'+GCM_runs[M_i][rr]+'_gn_'+str(yr)+'10-'+str(yr+10)+'12.nc'):# If this file does not alredy exist (hasn't been downloaded)
                            file_dir_name=file_dir+'/'+Var_names_vv+'_'+Var_AOL_vv+'mon_'+GCM_id2+'_dcppA-hindcast_s'+str(yr)+'-'+GCM_runs[M_i][rr]+'_gn_'+str(yr)+'10-'+str(yr+10)+'12.nc'
                            
                    elif GCM_id1 == 'MIROC' or GCM_id1 == 'NCAR' or GCM_id1 == 'MPI-M':
                        if os.path.isfile(file_dir+'/'+Var_names_vv+'_'+Var_AOL_vv+'mon_'+GCM_id2+'_dcppA-hindcast_s'+str(yr)+'-'+GCM_runs[M_i][rr]+'_gn_'+str(yr)+'11-'+str(yr+10)+'12.nc'):# If this file does not alredy exist (hasn't been downloaded)
                            file_dir_name=file_dir+'/'+Var_names_vv+'_'+Var_AOL_vv+'mon_'+GCM_id2+'_dcppA-hindcast_s'+str(yr)+'-'+GCM_runs[M_i][rr]+'_gn_'+str(yr)+'11-'+str(yr+10)+'12.nc'
                            
                    else:
                        if os.path.isfile(file_dir+'/'+Var_names_vv+'_'+Var_AOL_vv+'mon_'+GCM_id2+'_dcppA-hindcast_s'+str(yr)+'-'+GCM_runs[M_i][rr]+'_gn_'+str(yr)+'11-'+str(yr+10)+'12.nc'):# If this file does not alredy exist (hasn't been downloaded)
                            file_dir_name=file_dir+'/'+Var_names_vv+'_'+Var_AOL_vv+'mon_'+GCM_id2+'_dcppA-hindcast_s'+str(yr)+'-'+GCM_runs[M_i][rr]+'_gn_'+str(yr)+'11-'+str(yr+10)+'12.nc'

                    dset1 = Dataset(file_dir_name) # Specifies the netCDF dataset from which the data will be read
                    #dset1.variables
                    Lat_orig=np.asarray(dset1.variables['lat'])
                    Lon_orig=np.asarray(dset1.variables['lon'])   
                    
                    Data_orig=np.asarray(dset1.variables[Var_names_vv]) # phyc(time, rlat, rlon) - phydiat of all the desired months at the surface
                    Data_orig[np.abs(Data_orig) > 1e19] = nan # Replacing mising values with nan
                    dset1.close()
                    
                    if GCM_id1 == 'MIROC' or GCM_id1 == 'NCAR' or GCM_id1 == 'MPI-M': # For this model, the first 2 months are from the previous year
                        Data_orig=Data_orig[2:,:,:]
                    if GCM_id1 == 'NCC': # For this model, the first 3 months are from the previous year
                        Data_orig=Data_orig[3:,:,:]
                    
                    #Data_regrid = func_regrid(Data_orig, Lat_orig, Lon_orig, Lat_regrid_2D, Lon_regrid_2D) # Phydiat at surface - all months 
                    Data_regrid = Data_orig # We dont' regrid at this stage
                    
                    if Var_names_vv == 'pr':
                        Data_regrid = Data_regrid * 24 * 3600; # Converting prcp from kg/m2/s to kg/m2/day == mm/day   (kg/m2 = mm)
                    if Var_names_vv == 'tas':
                        Data_regrid = Data_regrid - 273.15; # Converting temp from Kelvin to C                   
                        
                    #############################################################
                    Data_regrid_month_ave=Data_regrid.reshape(10,12,Data_regrid.shape[1],Data_regrid.shape[2])
                    Data_regrid_month_ave=np.nanmean(Data_regrid_month_ave,axis=1)
                    
                    
                    ##########   Gridcell area weighted averaging   ###########
                    Lat_cell=Lat_orig; Lon_cell=Lon_orig; 
                    lat_n_regrid = Lat_cell.shape[0]; lon_n_regrid = Lon_cell.shape[0]
                    earth_R = 6378 # Earth Radius - Unit is kilometer (km)
                    GridCell_Areas = empty((lat_n_regrid, lon_n_regrid )) *nan
                    for ii in range(1,lat_n_regrid-1): # the first and last grid cells are not inlcuded, for the sake of simplification
                        for jj in range(1,lon_n_regrid-1): # the first and last grid cells are not inlcuded, for the sake of simplification
                            GridCell_Areas [ii,jj] = math.fabs( (earth_R**2) * (math.pi/180) * np.abs( (Lon_cell[jj]+Lon_cell[jj+1])/2 - (Lon_cell[jj]+Lon_cell[jj-1])/2 )  * ( math.sin(math.radians((Lat_cell[ii]+Lat_cell[ii+1])/2)) - math.sin(math.radians((Lat_cell[ii]+Lat_cell[ii-1])/2))) )
                    GridCell_Areas=GridCell_Areas / 1e6 # to convert the area to million km2
                    
                    var_area_weighted=[]
                    for ii in range(Data_regrid_month_ave.shape[0]):
                        var_area_weighted.append( np.nansum( np.multiply(Data_regrid_month_ave[ii,:,:], GridCell_Areas) ) / np.nansum(GridCell_Areas) )
                    
                    #Plot_var.append(np.nanmean(np.nanmean(Data_regrid_month_ave,axis=2),axis=1))
                    Plot_var.append(var_area_weighted)
                    print('Model: {}, Var: {}, Year: {}, Run#: {} - appended'.format(GCM_id2,Var_names_vv,yr,GCM_runs[M_i][rr]) )
                          
                Plot_var=np.asarray(Plot_var)
                Plot_var_U=np.nanmean(Plot_var,axis=0)+np.nanstd(Plot_var,axis=0) # Mean + 1 StDev
                Plot_var_D=np.nanmean(Plot_var,axis=0)-np.nanstd(Plot_var,axis=0) # Mean - 1 StDev
                
                r_b_count=r_b_count+1
                if r_b_count%2 == 0:
                    r_b='r'
                else:
                    r_b='b'
                
                im1=plt.fill_between(Plot_X, Plot_var_D, Plot_var_U, facecolor=r_b, color=r_b, alpha=0.2, linewidth=1.0)
                im2=plt.plot( Plot_X, np.nanmean(Plot_var,axis=0), c=r_b, label=GCM_id2, linewidth=3.0)

        ##########   Gridcell area weighted averaging for reanalysis data   ###########
        Lat_cell=Lat_orig_era5; Lon_cell=Lon_orig_era5; 
        lat_n_regrid = Lat_cell.shape[0]; lon_n_regrid = Lon_cell.shape[0]
        earth_R = 6378 # Earth Radius - Unit is kilometer (km)
        GridCell_Areas = empty((lat_n_regrid, lon_n_regrid )) *nan
        for ii in range(1,lat_n_regrid-1): # the first and last grid cells are not inlcuded, for the sake of simplification
            for jj in range(1,lon_n_regrid-1): # the first and last grid cells are not inlcuded, for the sake of simplification
                GridCell_Areas [ii,jj] = math.fabs( (earth_R**2) * (math.pi/180) * np.abs( (Lon_cell[jj]+Lon_cell[jj+1])/2 - (Lon_cell[jj]+Lon_cell[jj-1])/2 )  * ( math.sin(math.radians((Lat_cell[ii]+Lat_cell[ii+1])/2)) - math.sin(math.radians((Lat_cell[ii]+Lat_cell[ii-1])/2))) )
        GridCell_Areas=GridCell_Areas / 1e6 # to convert the area to million km2
    
        if Var_names_vv == 'pr':
            var_re = ERA5_PRCP_1979_2018_annual
        if Var_names_vv == 'tas':
            var_re = ERA5_TEMP_1979_2018_annual
        
        Plot_var_era5=[]
        for ii in range(var_re.shape[0]):
            Plot_var_era5.append( np.nansum( np.multiply(var_re[ii,:,:], GridCell_Areas) ) / np.nansum(GridCell_Areas) )        
        
        im2=plt.plot( np.linspace(1979,2018,40), Plot_var_era5, color='k', label='ERA5 ranalaysis',  linewidth=3.0)
        plt.legend(prop={'size': 24}, loc='best', fancybox=True, framealpha=0.8)
        plt.xticks(fontsize = 26); plt.yticks(fontsize = 26)
        plt.xlabel('Year', fontsize=18)
        plt.ylabel(Var_names_long[vv], fontsize=26)
        plt.title('CMIP6-DCPP hindcasts - Global average (shading= mean +/- st.dev.) - {} - {}  '.format(Var_names_long[vv],GCM_id2,GCM_runs[M_i][rr]), fontsize=18)
        mng = plt.get_current_fig_manager()
        mng.window.showMaximized() # Maximizes the plot window to save figures in full 
        ax = plt.axes()        
        ax.xaxis.grid() # horizontal lines
        ax.xaxis.set_minor_locator(MultipleLocator(5))
        ax.grid(which = 'minor')
        fig.savefig(dir_figs+'CMIP6_DCPP_Hindcasts_GlobalAve_'+GCM_id2+'_'+Var_names_vv+'_AllRuns.png', format='png', dpi=300, transparent=True, bbox_inches='tight')            
            
            























           
################################################################################  
################################################################################
#for M_i in range(len(GCM_Institution_ID)): # M_i=2
#    
#    GCM_id1 = GCM_Institution_ID[M_i]
#    GCM_id2 = GCM_Source_ID[M_i]
#    
#    for rr in range(0,len(GCM_runs[M_i])): # rr=0 # ensemble number
#        
#        
#        for vv in range(0, len(Var_names)): # vv=1
#
#            r_b_count=0 # a counter to alternate betwwen red and blur line colors            
#            fig=plt.figure()
#            
#            Var_names_vv=Var_names[vv]
#            Var_AOL_vv=Var_AOL[vv]
#            
#            file_dir = dir_data_in1+'/'+GCM_id1+'/'+'dcppA-hindcast'+'/'+GCM_runs[M_i][rr]+'/'+'mon'+'/'+Var_names_vv # Final path where the files will be saved
#            
#            for yr in range(1960,2019): # yr=1965 
#                
#                if GCM_id1 == 'CCCma':
#                    if os.path.isfile(file_dir+'/'+Var_names_vv+'_'+Var_AOL_vv+'mon_'+GCM_id2+'_dcppA-hindcast_s'+str(yr)+'-'+GCM_runs[M_i][rr]+'_gn_'+str(yr+1)+'01-'+str(yr+10)+'12.nc'):# If this file exists in this directory
#                        file_dir_name=file_dir+'/'+Var_names_vv+'_'+Var_AOL_vv+'mon_'+GCM_id2+'_dcppA-hindcast_s'+str(yr)+'-'+GCM_runs[M_i][rr]+'_gn_'+str(yr+1)+'01-'+str(yr+10)+'12.nc'
#                        
#                elif GCM_id1 == 'NCC':
#                    if os.path.isfile(file_dir+'/'+Var_names_vv+'_'+Var_AOL_vv+'mon_'+GCM_id2+'_dcppA-hindcast_s'+str(yr)+'-'+GCM_runs[M_i][rr]+'_gn_'+str(yr)+'10-'+str(yr+10)+'12.nc'):# If this file does not alredy exist (hasn't been downloaded)
#                        file_dir_name=file_dir+'/'+Var_names_vv+'_'+Var_AOL_vv+'mon_'+GCM_id2+'_dcppA-hindcast_s'+str(yr)+'-'+GCM_runs[M_i][rr]+'_gn_'+str(yr)+'10-'+str(yr+10)+'12.nc'
#                        
#                elif GCM_id1 == 'MIROC' or GCM_id1 == 'NCAR' or GCM_id1 == 'MPI-M':
#                    if os.path.isfile(file_dir+'/'+Var_names_vv+'_'+Var_AOL_vv+'mon_'+GCM_id2+'_dcppA-hindcast_s'+str(yr)+'-'+GCM_runs[M_i][rr]+'_gn_'+str(yr)+'11-'+str(yr+10)+'12.nc'):# If this file does not alredy exist (hasn't been downloaded)
#                        file_dir_name=file_dir+'/'+Var_names_vv+'_'+Var_AOL_vv+'mon_'+GCM_id2+'_dcppA-hindcast_s'+str(yr)+'-'+GCM_runs[M_i][rr]+'_gn_'+str(yr)+'11-'+str(yr+10)+'12.nc'
#                        
#                else:
#                    if os.path.isfile(file_dir+'/'+Var_names_vv+'_'+Var_AOL_vv+'mon_'+GCM_id2+'_dcppA-hindcast_s'+str(yr)+'-'+GCM_runs[M_i][rr]+'_gn_'+str(yr)+'11-'+str(yr+10)+'12.nc'):# If this file does not alredy exist (hasn't been downloaded)
#                        file_dir_name=file_dir+'/'+Var_names_vv+'_'+Var_AOL_vv+'mon_'+GCM_id2+'_dcppA-hindcast_s'+str(yr)+'-'+GCM_runs[M_i][rr]+'_gn_'+str(yr)+'11-'+str(yr+10)+'12.nc'
#
#                if yr%5 == 0: # Only take the years with 5yr increments, 1960,965,1970,...
#                    
#                    dset1 = Dataset(file_dir_name) # Specifies the netCDF dataset from which the data will be read
#                    #dset1.variables
#                    Lat_orig=np.asarray(dset1.variables['lat'])
#                    Lon_orig=np.asarray(dset1.variables['lon'])   
#                    
#                    Data_orig=np.asarray(dset1.variables[Var_names_vv]) # phyc(time, rlat, rlon) - phydiat of all the desired months at the surface
#                    Data_orig[np.abs(Data_orig) > 1e19] = nan # Replacing mising values with nan
#                    dset1.close()
#                    
#                    if GCM_id1 == 'MIROC' or GCM_id1 == 'NCAR' or GCM_id1 == 'MPI-M': # For this model, the first 2 months are from the previous year
#                        Data_orig=Data_orig[2:,:,:]
#                    if GCM_id1 == 'NCC': # For this model, the first 3 months are from the previous year
#                        Data_orig=Data_orig[3:,:,:]
#                    
#                    Data_regrid = func_regrid(Data_orig, Lat_orig, Lon_orig, Lat_regrid_2D, Lon_regrid_2D) # Phydiat at surface - all months 
#                    
#                    if Var_names_vv == 'pr':
#                        Data_regrid = Data_regrid * 24 * 3600; # Converting prcp from kg/m2/s to kg/m2/day == mm/day   (kg/m2 = mm)
#                    if Var_names_vv == 'tas':
#                        Data_regrid = Data_regrid - 273.15; # Converting temp from Kelvin to C                     
#                        
#                    #############################################################
#                    Data_regrid_month_ave=Data_regrid.reshape(10,12,Data_regrid.shape[1],Data_regrid.shape[2])
#                    Data_regrid_month_ave=np.nanmean(Data_regrid_month_ave,axis=1)
#                    
#                    r_b_count=r_b_count+1
#                    if r_b_count%2 == 0:
#                        r_b='r'
#                    else:
#                        r_b='b'
#                    
#                    im1=plt.plot( np.linspace(yr+1,yr+10,10), np.nanmean(np.nanmean(Data_regrid_month_ave,axis=2),axis=1), c=r_b, linewidth=3.0)
#            
#            
#            plt.xticks(fontsize = 26); plt.yticks(fontsize = 26)
#            plt.xlabel('Year', fontsize=18)
#            plt.ylabel(Var_names_long[vv], fontsize=26)
#            plt.title('CMIP6-DCPP Decadal Hindcasts, Global average - Variable: {} - Model: {} - Run#: {} '.format(Var_names_long[vv],GCM_id2,GCM_runs[M_i][rr]), fontsize=18)
#            mng = plt.get_current_fig_manager()
#            mng.window.showMaximized() # Maximizes the plot window to save figures in full 
#            ax = plt.axes()        
#            ax.xaxis.grid() # horizontal lines
#            ax.xaxis.set_minor_locator(MultipleLocator(5))
#            ax.grid(which = 'minor')
#            fig.savefig(dir_figs+'CMIP6_DCPP_Hindcasts_GlobalAve_'+GCM_id2+'_'+Var_names_vv+'_'+GCM_runs[M_i][rr]+'.png', format='png', dpi=300, transparent=True, bbox_inches='tight')            
#            
#            
#
#    
#        
#Plot_Var = np.nanmean(Data_regrid,axis=0)
##Plot_Var[ Ocean_Land_mask==0 ]=nan # masking over land, so grid cells that fall on land area (value=0) will be deleted
#cmap_limit=np.nanmax(np.abs( np.nanpercentile(Plot_Var, 99)))
#Plot_range=np.linspace(0,cmap_limit,29)
#Plot_unit='[ mm/day ]'; Plot_title= 'Average daily precipitation [mm/day] - '+str(yr)+'-'+str(yr+1)+' - '+str(GCM_id2)
#fig, m = func_plotmap_contourf(Plot_Var, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.jet, 'cyl', 0., 80., -80., '-')
#plt.show()
##fig.savefig(dir_figs+str(GCM)+'_Density_AveMap_1000m.png', format='png', dpi=300, transparent=True, bbox_inches='tight')           
#            
#fig=plt.figure()
#im1=plt.plot( np.linspace(yr+1,yr+10,10), np.nanmean(np.nanmean(Data_regrid_month_ave,axis=2),axis=1),  linewidth=3.0)
#im1=plt.plot( np.linspace(yr+1,yr+10,10)+5, np.nanmean(np.nanmean(Data_regrid_month_ave,axis=2),axis=1),  linewidth=3.0)
#plt.xticks(fontsize = 26); plt.yticks(fontsize = 26)
##plt.xticks(np.linspace(yr+1,yr+10,10))
##plt.ylabel(P_ylable, fontsize=18)
##plt.title(P_title, fontsize=18)
#mng = plt.get_current_fig_manager()
#mng.window.showMaximized() # Maximizes the plot window to save figures in full         
            
            
       
