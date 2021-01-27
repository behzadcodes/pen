### IMPORTANT: The following path needs to be appended so that the neccessary codes can be read from that directory ###
import sys
sys.path.append('/data1/home/basadieh/behzadcodes/behzadlibrary')
#from Behzadlib import * ### This is loaded from '/data1/home/basadieh/behzadcodes/behzadlibrary'
#from CMIP5lib import *  ### This is loaded from '/data1/home/basadieh/behzadcodes/behzadlibrary'
from Behzadlib import func_MLD, func_time_depth_plot, func_stream, func_ENSO, func_NAO, func_latlon_regrid, func_latlon_regrid_eq, func_oceanlandmask, func_oceanindex, func_regrid, func_plot, func_plot_save, func_plot_bounds_save
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
var_list_short=np.load('var_list_short.npy')

dir_pwd = os.getcwd() # Gets the current directory (and in which the code is placed)
dir_figs = (dir_pwd + '/Figures1/') # Directory to save figures
dir_results = (dir_pwd + '/Results/') # Directory to save results

GCM = 'GFDL-ESM2G'
year_start=1
year_end=500

#GCM_Names = ['GFDL-ESM2M', 'CanESM2','CESM1-BGC','CMCC-CESM','CNRM-CM5','GFDL-ESM2G','GFDL-ESM2M','GISS-E2-H-CC','GISS-E2-R-CC','HadGEM2-CC','HadGEM2-ES','IPSL-CM5A-LR','IPSL-CM5A-MR','IPSL-CM5B-LR','MIROC-ESM','MIROC-ESM-CHEM','MPI-ESM-LR','MPI-ESM-MR','MRI-ESM1','IPSL-CM5B-LR','NorESM1-ME']

### Regrdridding calculations ###
# creating new coordinate grid, same which was used in interpolation in data processing code
lat_n_regrid, lon_n_regrid = 180, 360 # Number of Lat and Lon elements in the regridded data
lon_min_regrid, lon_max_regrid = 0, 360 # Min and Max value of Lon in the regridded data
lat_min_regrid, lat_max_regrid = -90, 90 # Min and Max value of Lat in the regridded data

Lat_regrid_1D, Lon_regrid_1D, Lat_bound_regrid, Lon_bound_regrid = func_latlon_regrid_eq(lat_n_regrid, lon_n_regrid, lat_min_regrid, lat_max_regrid, lon_min_regrid, lon_max_regrid)
Lon_regrid_2D, Lat_regrid_2D = np.meshgrid(Lon_regrid_1D, Lat_regrid_1D)

Ocean_Land_mask = func_oceanlandmask(Lat_regrid_2D, Lon_regrid_2D) # 1= ocean, 0= land

Ocean_Index = func_oceanindex (Lat_regrid_2D, Lon_regrid_2D) # [0=land] [2=Pacific] [3=Indian Ocean] [6=Atlantic] [10=Arctic] [8=Baffin Bay (west of Greenland)] [9=Norwegian Sea (east of Greenland)] [11=Hudson Bay (Canada)] 

##################################################
#%% Convection Index Calculations - WS and LAB ###

conv_index_depth_ws=500 # args[7] in func_time_depth_plot code - depth for Convection Index
conv_index_depth_lab=500

dir_area=(dir_pwd + '/areacello_data/')
dset_cellarea = Dataset(dir_area+'areacello_fx_'+str(GCM)+'_historical_r0i0p0.nc')
CellArea=np.asarray(dset_cellarea.variables['areacello'])
dset_cellarea.close()

t_frequency='Omon'
variable_thetao='thetao' # Sea Water Temperature
variable_so='so' # Sea Water Salinity

dir_data_in1 = ('/data5/scratch/CMIP5/CMIP5_models/ocean_physics/') # Directory to raed raw data from
dir_data_in2=(dir_data_in1+ GCM + '/piControl/mo/')

dset_thetao = netcdf_read (dir_data_in2+str(variable_thetao)+'_'+str(t_frequency)+'_'+str(GCM)+'*12.nc',str(GCM), variable_thetao) # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
dset_so = netcdf_read (dir_data_in2+str(variable_so)+'_'+str(t_frequency)+'_'+str(GCM)+'*12.nc',str(GCM), variable_so)
Depths=dset_thetao.lvl[:]

### LAB_index and WS_index are the convection indeces at 500m depth ###
#### Weddel Sea Caclulations ####
month_ws=9 # args[2] in main code # month of the year
hemis_ws=0 # args[3] in main code # 0 means NH / 1 means SH
Conv_area_timeseries_WS, MLD_years_WS, lon_mld_WS, lat_mld_WS, WS_indeces_lonlat = func_MLD(dset_thetao, dset_so, month_ws, hemis_ws, year_start, year_end, CellArea, 90, 180)
MLD_average_WS=np.nanmean(MLD_years_WS,axis=0)

WS_index_time_depth, WS_index=func_time_depth_plot(dset_thetao, year_start, year_end, WS_indeces_lonlat, 90, 180, conv_index_depth_ws) # WS_index is the convection indeces at 500m depth
WS_index_norm=(WS_index-np.nanmean(WS_index))/np.std(WS_index) # Normalized Convection Index

WS_index_norm_rm = runningMeanFast(WS_index_norm, 10)#if data is not smoothed
WS_index_norm_rm=WS_index_norm_rm[:-9]

#### Labrador Sea Caclulations ####
month_lab=3 # args[2] in main code # month of the year 
hemis_lab=1 # args[3] in main code # 0 means NH / 1 means SH
Conv_area_timeseries_LAB, MLD_years_LAB, lon_mld_LAB, lat_mld_LAB, LAB_indeces_lonlat = func_MLD(dset_thetao, dset_so, month_lab, hemis_lab, year_start, year_end, CellArea, 90, 180)
MLD_average_LAB=np.nanmean(MLD_years_LAB,axis=0)

LAB_index_time_depth, LAB_index=func_time_depth_plot(dset_thetao, year_start, year_end, LAB_indeces_lonlat, 90, 180, conv_index_depth_lab) # LAB_index is the convection indeces at 500m depth
LAB_index_norm=(LAB_index-np.nanmean(LAB_index))/np.std(LAB_index) # Normalized Convection Index

LAB_index_norm_rm = runningMeanFast(LAB_index_norm, 10)#if data is not smoothed
LAB_index_norm_rm=LAB_index_norm_rm[:-9]
#########################################
dset_thetao.close_ncfile(dset_thetao.fin)
dset_so.close_ncfile(dset_so.fin)

#####################################
#%% Streamfunction Calculations ####

mask_atl=calc_Atl_Mask()
t_frequency='Omon'
variable='vo'# vo = sea_water_y_velocity - units: m s-1

dir_data_in1 = ('/data5/scratch/CMIP5/CMIP5_models/ocean_physics/') # Directory to raed raw data from
dir_data_in2=(dir_data_in1+ GCM + '/piControl/mo/')

dset_vo = netcdf_read (dir_data_in2+str(variable)+'_'+str(t_frequency)+'_'+str(GCM)+'*12.nc',str(GCM), variable) # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read

#the output structure:
#1st array - this is lon averaged 2d array with lats and depths of transport average over selected time period, can be used to calculate streamfunctions, see example below
#2nd array - time averaged, upper 1000m averaged transport, plotted as a map andd is default output of function of function
#3rd array - same as second butn 2000m-3000m depth average
#4,5,6,7 - indeces of lon/lat of max transport for upper 1000m and min transport for 2000-3000m, see example plots
Transport_lon_final, transport_0_1000, transport_1000_2000, transport_2000_3000, transport_3000_4000, transport_4000_5000, transport_5000_below, latlon_depths, Depth_indx, lon_stream, lat_stream, ii_max,jj_max, ii_min, jj_min, transport_0_1000_mean, transport_2000_3000_mean  = func_stream(dset_vo, year_start, year_end, mask_atl, 180, 360)

transport_1000_2000_mean=np.nanmean(transport_1000_2000,axis=0)
transport_3000_4000_mean=np.nanmean(transport_3000_4000,axis=0)
transport_4000_5000_mean=np.nanmean(transport_4000_5000,axis=0)
transport_5000_below_mean=np.nanmean(transport_5000_below,axis=0)

transport_1000_below_mean = np.mean( np.array([transport_1000_2000_mean, transport_2000_3000_mean, transport_3000_4000_mean, transport_4000_5000_mean, transport_5000_below_mean ]), axis=0 )
############################################
#### Streamfunction calculations Methode 1
#Stream_function=empty((Transport_lon_final.shape[0], Transport_lon_final.shape[1], Transport_lon_final.shape[2]))*nan # streamfunction
#Stream_function[:,0,:]=Transport_lon_final[:,0,:]
#for ii in range(1,Stream_function.shape[1]): # Depths
#    Stream_function[:,ii,:]=Stream_function[:,ii-1,:]+Transport_lon_final[:,ii,:]
##Stream_function=np.nancumsum(Transport_lon_final, axis=1)
#    
#Stream_function_ave=np.nanmean(Stream_function, axis=0) # Stream Function averaged over the years
#AMOC_max=np.nanmax(Stream_function, axis=1)
#SMOC_min=np.nanmin(Stream_function, axis=1)

###############################################
### Streamfunction calculations Methode Anna
Stream_function=empty((Transport_lon_final.shape[0], Transport_lon_final.shape[1], Transport_lon_final.shape[2]))*nan # streamfunction
Stream_function[:,0,:]=Transport_lon_final[:,0,:]
for ii in range(1,Stream_function.shape[1]): # Depths
    Stream_function[:,ii,:]=Stream_function[:,ii-1,:]+Transport_lon_final[:,ii,:]

stfunc_sum=np.nansum(Transport_lon_final,axis=1)
for ii in range(1,Stream_function.shape[1]): # Depths
    Stream_function[:,ii,:]=Stream_function[:,ii,:] - stfunc_sum

Stream_function_ave=np.nanmean(Stream_function, axis=0) # Stream Function averaged over the years
AMOC_max=np.nanmax(Stream_function, axis=1)
SMOC_min=np.nanmin(Stream_function, axis=1)

#######################################################
Transport_lon_final_mean=np.nanmean(Transport_lon_final, axis=0) #  just transport
Transport_lon_final_mean=np.asarray(Transport_lon_final_mean)
Transport_lon_final_mean=np.squeeze(Transport_lon_final_mean)

##################################################
#%%########     ENSO calculations     ############
##################################################
mask=calc_Atl_Mask()
t_frequency='Omon'
variable='thetao'

dir_data_in1 = ('/data5/scratch/CMIP5/CMIP5_models/ocean_physics/') # Directory to raed raw data from
dir_data_in2=(dir_data_in1+ GCM + '/piControl/mo/')

dset_thetao = netcdf_read (dir_data_in2+str(variable)+'_'+str(t_frequency)+'_'+str(GCM)+'*12.nc',str(GCM), variable) # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read

#finding indeces of our time range
start_date_i,end_date_i = dset_thetao.find_time(dset_thetao.times, year_start, year_end)
#order of args in function -  data, start_date_i, end_date_i, 1/0 for plotting yes/no
ENSO_index, ENSO_plot_data, lon_enso, lat_enso = func_ENSO(dset_thetao, year_start, year_end, 1)

dset_thetao.close_ncfile(dset_thetao.fin)
#################################################
#%%########     NAO Calculations    #############
#################################################
mask=calc_Atl_Mask()
t_frequency='Amon'
variable='psl'

dir_data_in1 = ('/data5/scratch/CMIP5/CMIP5_models/atmosphere_physics/') # Directory to raed raw data from
dir_data_in2=(dir_data_in1+ GCM + '/piControl/mo/')

dset_psl = netcdf_read (dir_data_in2+str(variable)+'_'+str(t_frequency)+'_'+str(GCM)+'*12.nc',str(GCM), variable) # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read

#finding indeces of our time range
start_date_i,end_date_i = dset_psl.find_time(dset_psl.times, year_start, year_end)
#order of args in function -  data, start_date_i, end_date_i, 1/0 for plotting yes/no
NAO_index, NAO_spatial_pattern, lon_nao, lat_nao = func_NAO(dset_psl, year_start, year_end, 1) ### NAO is first EOF of PSL

dset_psl.close_ncfile(dset_psl.fin)
####################################################################
#%%####           Westerlies and Trades Calculations        ########
####################################################################
GCM = 'GFDL-ESM2G'
t_frequency='Amon'
variable='tauu'

dir_data_in1 = ('/data5/scratch/CMIP5/CMIP5_models/atmosphere_physics/') # Directory to raed raw data from
dir_data_in2=(dir_data_in1+ GCM + '/piControl/mo/')

dset_tauu = netcdf_read (dir_data_in2+str(variable)+'_'+str(t_frequency)+'_'+str(GCM)+'*12.nc',str(GCM), variable) # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
start_date_i,end_date_i = dset_tauu.find_time(dset_tauu.times, year_start, year_end)

dset_tauu_extracted=dset_tauu.extract_data(dset_tauu.variable,start_date_i,end_date_i)
#here we in order to calculate yearly data from monthly data we first reshape array
dset_tauu_yearly=dset_tauu_extracted.reshape(( (year_end-year_start+1) ,12,len(dset_tauu_extracted[0]),len(dset_tauu_extracted[0][0])))
#second - we calculate mean over each 12 elements, which is second dimension of array
dset_tauu_yearly=np.mean(dset_tauu_yearly,axis=1)
#resulted yealry data we interpolate on regular grid
#lon,lat,Zonal_winds=interpolate_2_reg_grid(dset_tauu.x,dset_tauu.y,dset_tauu_yearly)
Zonal_winds = func_regrid(dset_tauu_yearly, dset_tauu.y, dset_tauu.x, Lat_regrid_2D, Lon_regrid_2D)
#average over all the lons
Zonal_winds=np.mean(Zonal_winds,axis=2)
# You can adjust the selected lats for calculating the wind indeces
North_Westerlies=Zonal_winds[:,135:145]
South_Westerlies=Zonal_winds[:,35:45]
North_Trades=Zonal_winds[:,100:110]
South_Trades=Zonal_winds[:,70:80]

dset_tauu.close_ncfile(dset_tauu.fin)
#################################################################
#%% surface_temperature_composites.py Calculations and Plots ####
t_frequency='Amon'
variable='ts'

dir_data_in1 = ('/data5/scratch/CMIP5/CMIP5_models/atmosphere_physics/') # Directory to raed raw data from
dir_data_in2=(dir_data_in1+ GCM + '/piControl/mo/')

dset_ts= netcdf_read(dir_data_in2+str(variable)+'_'+str(t_frequency)+'_'+str(GCM)+'*12.nc',str(GCM),variable)
start_date_i,end_date_i = dset_ts.find_time(dset_ts.times, year_start, year_end)

ii_LAB = np.where(LAB_index_norm_rm <-0.5)
jj_LAB = np.where(LAB_index_norm_rm >0.5)

ii_WS = np.where(WS_index_norm_rm <-0.5)
jj_WS = np.where(WS_index_norm_rm >0.5)

Airtemp_final_LAB=[]

for i in range(int((end_date_i+1-start_date_i)/12)):
    #print (int((end_date_i+1-start_date_i)/12))
    print('Airtem Composte - LAB - Year: ', i)
    data_vo_extracted=dset_ts.extract_data(dset_ts.variable,start_date_i+12*i,start_date_i+12*i+11)
    data=np.squeeze(data_vo_extracted)
    data=np.mean(data, axis=0)
    #x_i, y_i = np.meshgrid(dset_ts.y,dset_ts.x)
    #lon,lat,data_i=interpolate_2_reg_grid(dset_ts.x,dset_ts.y,data)
    data_i = func_regrid(data, dset_ts.y, dset_ts.x, Lat_regrid_2D, Lon_regrid_2D)
    data_i[data_i>10000000]=np.nan
    Airtemp_final_LAB.append(data_i)

Airtemp_final_LAB=np.asarray(Airtemp_final_LAB)
Airtemp_final_LAB_c=Airtemp_final_LAB[ii_LAB]
Airtemp_final_LAB_n=Airtemp_final_LAB[jj_LAB]
conv_LAB=np.nanmean(Airtemp_final_LAB_c,axis=0)
nonconv_LAB=np.nanmean(Airtemp_final_LAB_n,axis=0)

composite_LAB=np.subtract(conv_LAB,nonconv_LAB)

#######################################################
Airtemp_final_WS=[]

for i in range(int((end_date_i+1-start_date_i)/12)):
    #print (int((end_date_i+1-start_date_i)/12))
    print('Airtem Composte - WS - Year: ', i)
    data_vo_extracted=dset_ts.extract_data(dset_ts.variable,start_date_i+12*i,start_date_i+12*i+11)
    data=np.squeeze(data_vo_extracted)
    data=np.mean(data, axis=0)
    #x_i, y_i = np.meshgrid(dset_ts.y,dset_ts.x)
    #lon,lat,data_i=interpolate_2_reg_grid(dset_ts.x,dset_ts.y,data)
    data_i = func_regrid(data, dset_ts.y, dset_ts.x, Lat_regrid_2D, Lon_regrid_2D)
    data_i[data_i>10000000]=np.nan
    Airtemp_final_WS.append(data_i)

Airtemp_final_WS=np.asarray(Airtemp_final_WS)
Airtemp_final_WS_c=Airtemp_final_WS[ii_WS]
Airtemp_final_WS_n=Airtemp_final_WS[jj_WS]
conv_WS=np.nanmean(Airtemp_final_WS_c,axis=0)
nonconv_WS=np.nanmean(Airtemp_final_WS_n,axis=0)

composite_WS=np.subtract(conv_WS,nonconv_WS)

dset_ts.close_ncfile(dset_ts.fin)
###############################################################################
###########                 SAving Results                #####################
import os
import shelve

dir_pwd = os.getcwd() # Gets the current directory (and in which the code is placed)
filename_out = (dir_pwd + '/Results/AllResults_'+GCM+'_500yr.out') # Directory to save processed data

### To save
my_shelf = shelve.open(filename_out,'n') # 'n' for new

for key in var_list_short:
    try:
        my_shelf[key] = globals()[key]
    except TypeError:
        #
        # __builtins__, my_shelf, and imported modules can not be shelved.
        #
        print('ERROR shelving: {0}'.format(key))
my_shelf.close()

##############################################################################
#################        To restore:        ##################################
import os
import shelve

dir_pwd = os.getcwd() # Gets the current directory (and in which the code is placed)
filename_out = (dir_pwd + '/Results/AllResults_'+GCM+'_500yr.out') # Directory to save processed data

my_shelf = shelve.open(filename_out)
for key in my_shelf:
    globals()[key]=my_shelf[key]
my_shelf.close()
##############################################################################
##############################################################################

#### Composits Plost ####

fig=plt.figure()
m = Basemap( projection='mill',lon_0=210)
#m.fillcontinents(color='0.8')
m.drawmapboundary(fill_color='0.9')
m.drawmapboundary(fill_color='#000099')
m.drawcoastlines(linewidth=1.0, linestyle='solid', color='k', antialiased=1, ax=None, zorder=None)
m.drawparallels(np.arange(-90.,90.001,30.),labels=[True,False,False,False], linewidth=0.01) # labels = [left,right,top,bottom]
m.drawmeridians(np.arange(0.,360.,30.),labels=[False,False,False,True], linewidth=0.01) # labels = [left,right,top,bottom]
cmap_limit=np.nanmax(np.abs( np.nanpercentile(composite_LAB, 99.9)))
levels=np.linspace(-cmap_limit,cmap_limit,200)
im=m.contourf(lon,lat,composite_LAB,levels,latlon=True, cmap=plt.cm.seismic)
cb = m.colorbar(im,"right", size="3%", pad="2%")
#plt.clim(30,40)
#plt.clim(30,40)
#m.scatter(atl_mask[1,:],atl_mask[0,:], latlon=True)
plt.show()
plt.title('Surface air temperature composites - LAB - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_surface_temperature_composites_LAB.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


fig=plt.figure()
m = Basemap( projection='mill',lon_0=210)
#m.fillcontinents(color='0.8')
m.drawmapboundary(fill_color='0.9')
m.drawmapboundary(fill_color='#000099')
m.drawcoastlines(linewidth=1.0, linestyle='solid', color='k', antialiased=1, ax=None, zorder=None)
m.drawparallels(np.arange(-90.,90.001,30.),labels=[True,False,False,False], linewidth=0.01) # labels = [left,right,top,bottom]
m.drawmeridians(np.arange(0.,360.,30.),labels=[False,False,False,True], linewidth=0.01) # labels = [left,right,top,bottom]
cmap_limit=np.nanmax(np.abs( np.nanpercentile(composite_WS, 99.9)))
levels=np.linspace(-cmap_limit,cmap_limit,200)
im=m.contourf(lon,lat,composite_WS,levels,latlon=True, cmap=plt.cm.seismic)
cb = m.colorbar(im,"right", size="3%", pad="2%")
#plt.clim(30,40)
#plt.clim(30,40)
#m.scatter(atl_mask[1,:],atl_mask[0,:], latlon=True)
plt.show()
plt.title('Surface air temperature composites - WS - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_surface_temperature_composites_WS.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

dset_ts.close_ncfile(dset_ts.fin)

##################################################
##################################################
#%% Convection Index Plots - WS and LAB ###
###########################################
### Weddel Sea Plots ####
fig=plt.figure()
years=np.linspace(year_start,year_end,year_end-year_start+1)
#print (years,area)
plt.plot(years, Conv_area_timeseries_WS,'b')
#plt.plot(years[st2],Conv_area_timeseries_WS[st2],'r')
plt.xlabel('Years', fontsize=18)
plt.ylabel('Deep Convection Area (m2) - WS', fontsize=18)
plt.title('Deep Convection Area [m2] (Weddell Sea) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full    
fig.savefig(dir_figs+str(GCM)+'_deep_conv_area_WS.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
if hemis_ws==0:
    m = Basemap( projection='spstere',lon_0=0,boundinglat=-30)
else:
    m = Basemap( projection='npstere',lon_0=0,boundinglat=30)
#m = Basemap(projection='mill',lon_0=180)
m.drawcoastlines(linewidth=1.25)
m.fillcontinents(color='0.8')
#m.drawmapboundary(fill_color='#000099')
m.drawparallels(np.arange(-90,90,20), labels=[1,1,0,1])
m.drawmeridians(np.arange(0,360,30), labels=[1,1,0,1])
#lon[lon > 180] -= 360
im=m.contourf(lon_mld_WS,lat_mld_WS,MLD_average_WS,200,latlon=True, cmap=plt.cm.jet)
plt.colorbar(im)
m.scatter(lon_mld_WS[WS_indeces_lonlat],lat_mld_WS[WS_indeces_lonlat],2,latlon=True)
plt.title('Average MLD [m] (Weddell Sea) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_average_MLD_WS.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


fig=plt.figure()
years=np.linspace(year_start,year_end,year_end-year_start+1)
im = plt.contourf(years, Depths, np.transpose(WS_index_time_depth-273.15), 60, cmap=plt.cm.jet) # contour(X,Y,Z,N) - N shows the number of contour levels
plt.gca().invert_yaxis()
l = plt.axhline(y=conv_index_depth_ws)
plt.xlabel('Years', fontsize=18)
plt.ylabel('Depth [m]', fontsize=18)
plt.title('Temperature timeseries (C) in Convection Area (Weddell Sea) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
cbar = m.colorbar(im,"right", size="3%", pad="2%", extend='max') # extend='both' will extend the colorbar in both sides (upper side and down side)
cbar.set_label('(C)')
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_index_time_depth_WS.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
years=np.linspace(year_start,year_end,year_end-year_start+1)
im = plt.contourf(years, Depths[0:11], np.transpose(WS_index_time_depth[:,0:11]-273.15), 40, cmap=plt.cm.jet)
plt.gca().invert_yaxis()
plt.xlabel('Years', fontsize=18)
plt.ylabel('Depth [m]', fontsize=18)
plt.title('Temperature timeseries (C) in Convection Area (Weddell Sea) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
cbar = m.colorbar(im,"right", size="3%", pad="2%", extend='max') # extend='both' will extend the colorbar in both sides (upper side and down side)
cbar.set_label('(C)')
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_index_time_depth_surface_WS.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

#### Labrador Sea Plots ####
fig=plt.figure()
years=np.linspace(year_start,year_end,year_end-year_start+1)
#print (years,area)
plt.plot(years,Conv_area_timeseries_LAB,'b')
#plt.plot(years[st2],Conv_area_timeseries_LAB[st2],'r')
plt.xlabel('Years', fontsize=18)
plt.ylabel('Deep Convection Area (m2) - LAB', fontsize=18)
plt.title('Deep Convection Area [m2] (Labrador Sea) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full    
fig.savefig(dir_figs+str(GCM)+'_deep_conv_area_LAB.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
if hemis_lab==0:
    m = Basemap( projection='spstere',lon_0=0,boundinglat=-30)
else:
    m = Basemap( projection='npstere',lon_0=0,boundinglat=30)
#m = Basemap(projection='mill',lon_0=180)
m.drawcoastlines(linewidth=1.25)
m.fillcontinents(color='0.8')
#m.drawmapboundary(fill_color='#000099')
m.drawparallels(np.arange(-90,90,20), labels=[1,1,0,1])
m.drawmeridians(np.arange(0,360,30), labels=[1,1,0,1])
#lon[lon > 180] -= 360
im=m.contourf(lon_mld_LAB,lat_mld_LAB,MLD_average_LAB,200,latlon=True, cmap=plt.cm.jet)
plt.colorbar(im)
m.scatter(lon_mld_LAB[LAB_indeces_lonlat],lat_mld_LAB[LAB_indeces_lonlat],2,latlon=True)
plt.title('Average MLD [m] (Labrador Sea) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_average_MLD_LAB.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


fig=plt.figure()
years=np.linspace(year_start,year_end,year_end-year_start+1)
im = plt.contourf(years, Depths[0:42], np.transpose(LAB_index_time_depth[:,0:42]-273.15), 60, cmap=plt.cm.jet) # contour(X,Y,Z,N) - N shows the number of contour levels
plt.gca().invert_yaxis()
l = plt.axhline(y=conv_index_depth_lab)
plt.xlabel('Years', fontsize=18)
plt.ylabel('Depth [m]', fontsize=18)
plt.title('Temperature timeseries (C) in Convection Area (Labrador Sea) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
cbar = m.colorbar(im,"right", size="3%", pad="2%", extend='max') # extend='both' will extend the colorbar in both sides (upper side and down side)
cbar.set_label('(C)')
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_index_time_depth_LAB.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
years=np.linspace(year_start,year_end,year_end-year_start+1)
im = plt.contourf(years, Depths[0:11], np.transpose(LAB_index_time_depth[:,0:11]-273.15), 40, cmap=plt.cm.jet)
plt.gca().invert_yaxis()
plt.xlabel('Years', fontsize=18)
plt.ylabel('Temperature at depth', fontsize=18)
plt.title('Annual Temperature in Convection Area (Labrador Sea) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
cbar = m.colorbar(im,"right", size="3%", pad="2%", extend='max') # extend='both' will extend the colorbar in both sides (upper side and down side)
cbar.set_label('(C)')
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_index_time_depth_surface_LAB.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

###########################################
### Convection Index Time Series Plots ####
fig=plt.figure()
years=np.linspace(year_start,year_end,year_end-year_start+1)
#print (years,area)
plt.plot(years,WS_index_norm,'b')
l = plt.axhline(y=0)
#plt.plot(years[st2],deep_conv_area_WS[st2],'r')
plt.xlabel('Years', fontsize=18)
plt.ylabel('Convection index (°C)', fontsize=18)
plt.title('Weddel Sea convection index (normalized) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full    
fig.savefig(dir_figs+str(GCM)+'_convec_index_WS.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
years=np.linspace(year_start,year_end,year_end-year_start+1)
#print (years,area)
plt.plot(years,LAB_index_norm,'b')
l = plt.axhline(y=0)
#plt.plot(years[st2],deep_conv_area_WS[st2],'r')
plt.xlabel('Years', fontsize=18)
plt.ylabel('Convection index (°C)', fontsize=18)
plt.title('Labrador Sea convection index (normalized) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full    
fig.savefig(dir_figs+str(GCM)+'_convec_index_LAB.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

#######################################
#%% Streamfunction Plots ####

Plot_title=('Ocean depth map (Atlantic) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM))
Plot_save_dir=(dir_figs+str(GCM)+'_depth_map_atl.png')
bounds_max=np.int(np.nanpercentile(latlon_depths, 99.99))
bounds = np.arange(0, bounds_max, bounds_max/40)
#func_plot(latlon_depths, lat_stream, lon_stream, bounds_max, '-', '-', 'mill', 0)
func_plot_bounds_save(latlon_depths, lat_stream, lon_stream, bounds, '(m)', Plot_title, 'mill', 0, Plot_save_dir)

Plot_title=('Mean northward transport upper 1000m [Sv] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM))
Plot_save_dir=(dir_figs+str(GCM)+'_transport_0_1000_mean.png')
#bounds_max=np.nanpercentile(transport_0_1000_mean, 99.99)
bounds_max=float("{0:.02f}".format(np.nanpercentile(transport_0_1000_mean, 98))) # Upper bound of plotted values to be used for colorbar, which is 98th percentile of data, with 2 decimals
bounds = np.arange(-1*bounds_max, bounds_max, bounds_max/20)
#bounds = np.arange(-0.1, 0.1, 0.2/20)
func_plot_bounds_save(transport_0_1000_mean, lat_stream, lon_stream, bounds, '(Sv)', Plot_title, 'mill', 0, Plot_save_dir)  

Plot_title=('Mean northward transport at 2000m-3000m [Sv] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM))
Plot_save_dir=(dir_figs+str(GCM)+'_transport_2000_3000_mean.png')
#bounds_max=np.nanpercentile(transport_0_1000_mean, 99.99)
bounds_max=float("{0:.02f}".format(np.nanpercentile(transport_2000_3000_mean, 98))) # Upper bound of plotted values to be used for colorbar, which is 98th percentile of data, with 2 decimals
bounds = np.arange(-1*bounds_max, bounds_max, bounds_max/20)
#bounds = np.arange(-0.1, 0.1, 0.2/20)
func_plot_bounds_save(transport_2000_3000_mean, lat_stream, lon_stream, bounds, '(Sv)', Plot_title, 'mill', 0, Plot_save_dir)     


levels=np.linspace(-35,35,100)  
#plot #1, just transport
fig=plt.figure()
im=plt.contourf(Lat_regrid_1D, Depths, Transport_lon_final_mean, 100, cmap=plt.cm.RdBu_r, extend = 'both')
plt.gca().invert_yaxis()
plt.xlabel('Latitude', fontsize=18)
plt.ylabel('Depth', fontsize=18)
plt.title('Mean northward transport (Atlantic) [Sv] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
cbar = m.colorbar(im,"right", size="3%", pad="2%") # extend='both' will extend the colorbar in both sides (upper side and down side)
cbar.set_label('(Sv)')
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_Transport_lon_final_mean.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
                
#plot #2, streamfunction
fig=plt.figure()
im=plt.contourf(Lat_regrid_1D, Depths, Stream_function_ave, levels, cmap=plt.cm.RdBu_r, extend = 'both')
plt.gca().invert_yaxis()
plt.xlabel('Latitude', fontsize=18)
plt.ylabel('Depth', fontsize=18)
plt.title('Stream function (Atlantic) [Sv] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
cbar = m.colorbar(im,"right", size="3%", pad="2%") # extend='both' will extend the colorbar in both sides (upper side and down side)
cbar.set_label('(Sv)')
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_Stream_function_ave.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

##plot #3, map of mean transport upper 1000m with indeces of max for each lat
#fig=plt.figure()
#m = Basemap( projection='mill',lon_0=0)
#m.fillcontinents(color='0.8')
#m.drawparallels(np.arange(-90,90,20), labels=[1,1,0,1])
#m.drawmeridians(np.arange(0,360,30), labels=[1,1,0,1])
#im=m.contourf(lon_stream,lat_stream,transport_0_1000_mean,200,latlon=True, cmap=plt.cm.jet)
#m.scatter(lon[ii_max,jj_max],lat[ii_max,jj_max],1,c='k',latlon=True)
#cbar = m.colorbar(im,"right", size="3%", pad="6%", extend = 'both') # extend='both' will extend the colorbar in both sides (upper side and down side)
#cbar.set_label('(Sv)')
#plt.title('Mean transport upper 1000m (Stiplings= max transport locations) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
#mng = plt.get_current_fig_manager()
#mng.window.showMaximized() # Maximizes the plot window to save figures in full
#fig.savefig(dir_figs+str(GCM)+'_transport_0_1000_mean.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

#plot #3, map of mean transport upper 1000m with indeces of max for each lat
fig=plt.figure()
m = Basemap( projection='mill',lon_0=0)
m.fillcontinents(color='0.8')
m.drawparallels(np.arange(-90,90,20), labels=[1,1,0,1])
m.drawmeridians(np.arange(0,360,30), labels=[1,1,0,1])
bounds_max=float("{0:.02f}".format(np.nanpercentile(transport_0_1000_mean, 99.9))) # Upper bound of plotted values to be used for colorbar, which is 98th percentile of data, with 2 decimals
levels=np.linspace(-1*bounds_max,bounds_max,100)  
im=m.contourf(lon_stream,lat_stream,transport_0_1000_mean,levels,latlon=True, cmap=plt.cm.cool)
m.scatter(lon[ii_max,jj_max],lat[ii_max,jj_max],1,c='k',latlon=True)
cbar = m.colorbar(im,"right", size="3%", pad="6%", extend = 'both') # extend='both' will extend the colorbar in both sides (upper side and down side)
cbar.set_label('(Sv)')
plt.title('Mean transport upper 1000m (Stiplings= max transport locations) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_transport_0_1000_mean_stiple.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

#plot #4, map of mean transport 2000m-3000m with indeces of max for each lat
fig=plt.figure()
m = Basemap( projection='mill',lon_0=0)
m.fillcontinents(color='0.8')
m.drawparallels(np.arange(-90,90,20), labels=[1,1,0,1])
m.drawmeridians(np.arange(0,360,30), labels=[1,1,0,1])
bounds_max=float("{0:.02f}".format(np.nanpercentile(transport_2000_3000_mean, 99.9))) # Upper bound of plotted values to be used for colorbar, which is 98th percentile of data, with 2 decimals
levels=np.linspace(-1*bounds_max,bounds_max,100)  
im=m.contourf(lon_stream,lat_stream,transport_2000_3000_mean,levels,latlon=True, cmap=plt.cm.cool)
m.scatter(lon[ii_min,jj_min],lat[ii_min,jj_min],1,c='k',latlon=True)
plt.show()
cbar = m.colorbar(im,"right", size="3%", pad="6%", extend = 'both') # extend='both' will extend the colorbar in both sides (upper side and down side)
cbar.set_label('(Sv)')
plt.title('Mean transport 2000-3000m (Stiplings= min transport locations) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_transport_2000_3000_mean_stiple.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


dset_vo.close_ncfile(dset_vo.fin)
###########################################
#%%########     ENSO Plots     ############
###########################################
m = Basemap( projection='mill',lon_0=210)
m.fillcontinents(color='0.8')
m.drawmapboundary(fill_color='0.9')
m.drawmapboundary(fill_color='#000099')
m.drawparallels(np.arange(-90.,90.001,30.),labels=[True,False,False,False], linewidth=0.01) # labels = [left,right,top,bottom]
m.drawmeridians(np.arange(0.,360.,30.),labels=[False,False,False,True], linewidth=0.01) # labels = [left,right,top,bottom]
plt.title('ENSO map - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
mask_enso=ma.masked_where(ENSO_plot_data == np.nan, ENSO_plot_data)
im=m.contourf(lon_enso,lat_enso,mask_enso,200,latlon=True, cmap=plt.cm.jet)
cb = m.colorbar(im,"right", size="3%", pad="2%")
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_ENSO_map.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
years=np.linspace(year_start, year_end, year_end-year_start+1)
plt.plot(years,ENSO_index, 'k') 
y2=np.zeros(len(ENSO_index))
plt.fill_between(years, ENSO_index, y2, where=ENSO_index >= y2, color = 'r', interpolate=True)
plt.fill_between(years, ENSO_index, y2, where=ENSO_index <= y2, color = 'b', interpolate=True)
plt.axhline(linewidth=1, color='k')
plt.title('ENSO Index - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
plt.xlabel('Years', fontsize=18)
plt.ylabel('ENSO Index', fontsize=18)    
plt.show()
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_ENSO_index.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

##########################################
#%%########     NAO Plots    #############
##########################################
m = Basemap( projection='mill',lon_0=0)
m.fillcontinents(color='0.8')
m.drawmapboundary(fill_color='0.9')
m.drawmapboundary(fill_color='#000099')
m.drawparallels(np.arange(-90.,90.001,30.),labels=[True,False,False,False], linewidth=0.01) # labels = [left,right,top,bottom]
m.drawmeridians(np.arange(0.,360.,30.),labels=[False,False,False,True], linewidth=0.01) # labels = [left,right,top,bottom]
cmap_limit=np.nanmax(np.abs(NAO_spatial_pattern))
levels=np.linspace(-cmap_limit,cmap_limit,200)
plt.title('NAO (1st EOF of SLP) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
mask_nao=ma.masked_where(NAO_spatial_pattern == np.nan, NAO_spatial_pattern)
lon_nao[lon_nao > 180] -= 360
im=m.contourf(lon_nao,lat_nao,mask_nao,levels,latlon=True, cmap=plt.cm.seismic)
cb = m.colorbar(im,"right", size="3%", pad="2%")
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_NAO_map.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


fig=plt.figure()
years=np.linspace(year_start, year_end, year_end-year_start+1)
plt.plot(years,NAO_index, 'k') 
y2=np.zeros(len(NAO_index))
plt.fill_between(years, NAO_index, y2, where=NAO_index >= y2, color = 'r', interpolate=True)
plt.fill_between(years, NAO_index, y2, where=NAO_index <= y2, color = 'b', interpolate=True)
plt.axhline(linewidth=1, color='k')
plt.title('NAO Index - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
plt.xlabel('Years', fontsize=18)
plt.ylabel('NAO Index', fontsize=18)    
plt.show()
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_NAO_index.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

dset_psl.close_ncfile(dset_psl.fin)
############################################
#%% transport_lag_cor_Atlantic.py Plots ####

AMOC_transport_all = copy.deepcopy(transport_0_1000)
SMOC_transport_all= copy.deepcopy(transport_2000_3000)
SMOC_transport_all=SMOC_transport_all*(-1)
#As SMOC is transport averaged over 2000-3000m it is southward flow, thus the sign is negative (for southward flow) so we multiply by -1 to calculate further statistics

LAB_index_rm = copy.deepcopy(LAB_index)
WS_index_rm = copy.deepcopy(WS_index)

LAB_index_rm=runningMeanFast(LAB_index_rm, 10)
#multiplying by -1 as we have convection when temperature is low(opposite to the index which we calculate)
LAB_index_rm=LAB_index_rm*(-1)
WS_index_rm=runningMeanFast(WS_index_rm, 10)
WS_index_rm=WS_index_rm*(-1)

#summing transport over all longtitudes
AMOC_transport=np.nansum(AMOC_transport_all,axis=2)
SMOC_transport=np.nansum(SMOC_transport_all,axis=2)

###LAB plots###
lag_time=40

fig=plt.figure()
data=[]
for i in range(len(AMOC_transport[0][20:161])):
    stream=runningMeanFast(AMOC_transport[:,i+lag_time], 10)
    r=lag_cor_data(LAB_index_rm,stream,lag_time)
    data.append(r)
data=np.asarray(data)
x=np.linspace(-lag_time,lag_time+1, 2*lag_time)
y=np.linspace(-70,70+1, 2*70+1)

x,y=np.meshgrid(x,y)
cmap_limit=np.nanmax(np.abs(data))
levels=np.linspace(-cmap_limit,cmap_limit,200)
im=plt.contourf(x,y,data,levels,latlon=True, cmap=plt.cm.jet)
cb = plt.colorbar(im)
plt.xlabel('Year lag', fontsize=18)
plt.ylabel('Latitude', fontsize=18)
plt.title('Labrador peak convection vs. transport in upper 1000m - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_convec_transport_0_1000_LAB.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


fig=plt.figure()
data=[]
for i in range(len(SMOC_transport[0][20:161])):
    stream=runningMeanFast(SMOC_transport[:,i+lag_time], 10)
    r=lag_cor_data(LAB_index_rm,stream,lag_time)
    data.append(r)
data=np.asarray(data)
x=np.linspace(-lag_time,lag_time+1, 2*lag_time)
y=np.linspace(-70,70+1, 2*70+1)

x,y=np.meshgrid(x,y)
cmap_limit=np.nanmax(np.abs(data))
levels=np.linspace(-cmap_limit,cmap_limit,200)
im=plt.contourf(x,y,data,levels,latlon=True, cmap=plt.cm.jet)
cb = plt.colorbar(im)
plt.xlabel('Year lag', fontsize=18)
plt.ylabel('Latitude', fontsize=18)
plt.title('Labrador peak convection vs. transport in 2000m-3000m - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_convec_transport_2000_3000_LAB.png', format='png', dpi=300, transparent=True, bbox_inches='tight')



fig=plt.figure()
data=[]
for i in range(len(SMOC_transport[0][20:161])):
    stream=runningMeanFast(SMOC_transport[:,i+lag_time], 10)
    r=lag_cor_data(LAB_index_rm,stream,lag_time)
    data.append(r)
data=np.asarray(data)
x=np.linspace(-lag_time,lag_time+1, 2*lag_time)
y=np.linspace(-70,70+1, 2*70+1)

x,y=np.meshgrid(x,y)
cmap_limit=np.nanmax(np.abs(data))
levels=np.linspace(-cmap_limit,cmap_limit,200)
im=plt.contourf(x,y,data,levels,latlon=True, cmap=plt.cm.jet)
cb = plt.colorbar(im)
plt.xlabel('Year lag', fontsize=18)
plt.ylabel('Latitude', fontsize=18)
plt.title('Labrador peak convection vs. transport in 2000m-3000m - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_convec_transport_2000_3000_LAB.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

#### WS plots ####
lag_time=40

fig=plt.figure()
data=[]
for i in range(len(AMOC_transport[0][20:161])):
    stream=runningMeanFast(AMOC_transport[:,i+lag_time], 10)
    r=lag_cor_data(WS_index_rm,stream,lag_time)
    data.append(r)
data=np.asarray(data)
x=np.linspace(-lag_time,lag_time+1, 2*lag_time)
y=np.linspace(-70,70+1, 2*70+1)

x,y=np.meshgrid(x,y)
cmap_limit=np.nanmax(np.abs(data))
levels=np.linspace(-cmap_limit,cmap_limit,200)
im=plt.contourf(x,y,data,levels,latlon=True, cmap=plt.cm.jet)
cb = plt.colorbar(im)
plt.xlabel('Year lag', fontsize=18)
plt.ylabel('Latitude', fontsize=18)
plt.title('WS peak convection vs. transport in upper 1000m - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_convec_transport_0_1000_WS.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


fig=plt.figure()
data=[]
for i in range(len(SMOC_transport[0][20:161])):
    stream=runningMeanFast(SMOC_transport[:,i+lag_time], 10)
    r=lag_cor_data(WS_index_rm,stream,lag_time)
    data.append(r)
data=np.asarray(data)
x=np.linspace(-lag_time,lag_time+1, 2*lag_time)
y=np.linspace(-70,70+1, 2*70+1)

x,y=np.meshgrid(x,y)
cmap_limit=np.nanmax(np.abs(data))
levels=np.linspace(-cmap_limit,cmap_limit,200)
im=plt.contourf(x,y,data,levels,latlon=True, cmap=plt.cm.jet)
cb = plt.colorbar(im)
plt.xlabel('Year lag', fontsize=18)
plt.ylabel('Latitude', fontsize=18)
plt.title('WS peak convection vs. transport in 2000m-3000m - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_convec_transport_2000_3000_WS.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

##############################################################################
##############################################################################
###LAB plots###
lag_time=40

fig=plt.figure()
data=[]
for i in range(len(AMOC_max[0][20:161])):
    stream=runningMeanFast(AMOC_max[:,i+lag_time], 10)
    r=lag_cor_data(LAB_index_rm,stream,lag_time)
    data.append(r)
data=np.asarray(data)
x=np.linspace(-lag_time,lag_time+1, 2*lag_time)
y=np.linspace(-70,70+1, 2*70+1)

x,y=np.meshgrid(x,y)
cmap_limit=np.nanmax(np.abs(data))
levels=np.linspace(-cmap_limit,cmap_limit,200)
im=plt.contourf(x,y,data,levels,latlon=True, cmap=plt.cm.jet)
cb = plt.colorbar(im)
plt.xlabel('Year lag', fontsize=18)
plt.ylabel('Longitude', fontsize=18)
plt.title('Labrador peak convection vs. AMOC (max streamfunction) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_convec_AMOC_max_LAB.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


fig=plt.figure()
data=[]
for i in range(len(SMOC_min[0][20:161])):
    stream=runningMeanFast(SMOC_min[:,i+lag_time], 10)
    r=lag_cor_data(LAB_index_rm,stream,lag_time)
    data.append(r)
data=np.asarray(data)
x=np.linspace(-lag_time,lag_time+1, 2*lag_time)
y=np.linspace(-70,70+1, 2*70+1)

x,y=np.meshgrid(x,y)
cmap_limit=np.nanmax(np.abs(data))
levels=np.linspace(-cmap_limit,cmap_limit,200)
im=plt.contourf(x,y,data,levels,latlon=True, cmap=plt.cm.jet)
cb = plt.colorbar(im)
plt.xlabel('Year lag', fontsize=18)
plt.ylabel('Latitude', fontsize=18)
plt.title('Labrador peak convection vs. SMOC (min streamfunction) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_convec_SMOC_min_LAB.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

#### WS plots ####
lag_time=40

fig=plt.figure()
data=[]
for i in range(len(AMOC_max[0][20:161])):
    stream=runningMeanFast(AMOC_max[:,i+lag_time], 10)
    r=lag_cor_data(WS_index_rm,stream,lag_time)
    data.append(r)
data=np.asarray(data)
x=np.linspace(-lag_time,lag_time+1, 2*lag_time)
y=np.linspace(-70,70+1, 2*70+1)

x,y=np.meshgrid(x,y)
cmap_limit=np.nanmax(np.abs(data))
levels=np.linspace(-cmap_limit,cmap_limit,200)
im=plt.contourf(x,y,data,levels,latlon=True, cmap=plt.cm.jet)
cb = plt.colorbar(im)
plt.xlabel('Year lag', fontsize=18)
plt.ylabel('Latitude', fontsize=18)
plt.title('WS peak convection vs. AMOC (max streamfunction) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_convec_AMOC_max_WS.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


fig=plt.figure()
data=[]
for i in range(len(SMOC_min[0][20:161])):
    stream=runningMeanFast(SMOC_min[:,i+lag_time], 10)
    r=lag_cor_data(WS_index_rm,stream,lag_time)
    data.append(r)
data=np.asarray(data)
x=np.linspace(-lag_time,lag_time+1, 2*lag_time)
y=np.linspace(-70,70+1, 2*70+1)

x,y=np.meshgrid(x,y)
cmap_limit=np.nanmax(np.abs(data))
levels=np.linspace(-cmap_limit,cmap_limit,200)
im=plt.contourf(x,y,data,levels,latlon=True, cmap=plt.cm.jet)
cb = plt.colorbar(im)
plt.xlabel('Year lag', fontsize=18)
plt.ylabel('Latitude', fontsize=18)
plt.title('WS peak convection vs. SMOC (min streamfunction) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_convec_SMOC_min_WS.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

#######################################################
#%%####   Lagged Corrolations and PSD plots    ########
#######################################################

# extracting time series only  fro 30Nlatitude
Transport_lon_30N=Transport_lon_final[:,:,120]
# finding max over depth
Transport_lon_30N=np.nanmax(Transport_lon_30N,axis=1)

### new figure for lagged cors
NAO_m=runningMeanFast(NAO_index, 10)
ENSO_m=runningMeanFast(ENSO_index, 10)
Transport_lon_30N_m=runningMeanFast(Transport_lon_30N, 10)

AMOC_transport_50S=AMOC_transport[:,40] # AMOC at 50S
AMOC_transport_50S_m=runningMeanFast(AMOC_transport_50S, 10)
AMOC_transport_30S=AMOC_transport[:,60] # AMOC at 30S
AMOC_transport_30S_m=runningMeanFast(AMOC_transport_30S, 10)
AMOC_transport_50N=AMOC_transport[:,140] # AMOC at 50N
AMOC_transport_50N_m=runningMeanFast(AMOC_transport_50N, 10)
AMOC_transport_30N=AMOC_transport[:,120] # AMOC at 30N
AMOC_transport_30N_m=runningMeanFast(AMOC_transport_30N, 10)

AMOC_max_50S= AMOC_max[:,40] # AMOC at 50S # Max of streamfunction method
AMOC_max_50S_m=runningMeanFast(AMOC_max_50S, 10)
AMOC_max_30S= AMOC_max[:,60] # AMOC at 30S # Max of streamfunction method
AMOC_max_30S_m=runningMeanFast(AMOC_max_30S, 10)
AMOC_max_50N= AMOC_max[:,140] # AMOC at 50N # Max of streamfunction method
AMOC_max_50N_m=runningMeanFast(AMOC_max_50N, 10)
AMOC_max_30N= AMOC_max[:,120] # AMOC at 30N # Max of streamfunction method
AMOC_max_30N_m=runningMeanFast(AMOC_max_30N, 10)
AMOC_max_eq= AMOC_max[:,90] # AMOC at 30N # Max of streamfunction method
AMOC_max_eq_m=runningMeanFast(AMOC_max_eq, 10)

SMOC_min_50S= SMOC_min[:,40] # SMOC at 50S # Max of streamfunction method
SMOC_min_50S_m=runningMeanFast(SMOC_min_50S, 10)
SMOC_min_30S= SMOC_min[:,60] # SMOC at 30S # Max of streamfunction method
SMOC_min_30S_m=runningMeanFast(SMOC_min_30S, 10)
SMOC_min_50N= SMOC_min[:,140] # SMOC at 50N # Max of streamfunction method
SMOC_min_50N_m=runningMeanFast(SMOC_min_50N, 10)
SMOC_min_30N= SMOC_min[:,120] # SMOC at 30N # Max of streamfunction method
SMOC_min_30N_m=runningMeanFast(SMOC_min_30N, 10)
SMOC_min_eq= SMOC_min[:,90] # SMOC at 30N # Max of streamfunction method
SMOC_min_eq_m=runningMeanFast(SMOC_min_eq, 10)

#### Calulating PSD using welch method ####
fig=plt.figure()
f1,PSD1=plot_PSD_welch(ENSO_index,1,'r','ENSO index')
f2,PSD2=plot_PSD_welch(NAO_index,1,'c','NAO index')
f3,PSD3=plot_PSD_welch(Transport_lon_30N,1,'y','Transport 30N')
plt.legend()
plt.gca().invert_xaxis()
plt.title('Normalized PSD - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
plt.xlabel('Period (years)', fontsize=18)
plt.ylabel('Spectral Density', fontsize=18)    
plt.show()
plt.legend(shadow=True, loc='lower left', prop={'size': 20})
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_PSD_tans30N_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


fig=plt.figure()
lag_cor(Transport_lon_30N,ENSO_m,40,'r','ENSO index')
lag_cor(Transport_lon_30N,NAO_m,40,'b','NAO index')
plt.legend()
plt.show()
plt.title('Lagged Correlation with respect to max streamfunction at 30N - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
plt.xlabel('Years lag', fontsize=18)
plt.ylabel('Correlation coefficient', fontsize=18)    
plt.show()
plt.legend(shadow=True, loc='upper left', prop={'size': 20})
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_tans30N_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


fig=plt.figure()
lag_cor(WS_index_rm,ENSO_m,40,'r','ENSO index')
lag_cor(WS_index_rm,NAO_m,40,'b','NAO index')
plt.legend()
plt.show()
plt.title('Lagged Correlation with respect to WS convection - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
plt.xlabel('Years lag', fontsize=18)
plt.ylabel('Correlation coefficient', fontsize=18)    
plt.show()
plt.legend(shadow=True, loc='upper left', prop={'size': 20})
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_WSconvection_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
lag_cor(LAB_index_rm,ENSO_m,40,'r','ENSO index')
lag_cor(LAB_index_rm,NAO_m,40,'b','NAO index')
plt.legend()
plt.show()
plt.title('Lagged Correlation with respect to LAB convection - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
plt.xlabel('Years lag', fontsize=18)
plt.ylabel('Correlation coefficient', fontsize=18)    
plt.show()
plt.legend(shadow=True, loc='upper left', prop={'size': 20})
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_LABconvection_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

#################################################
### AMOC Lagged Corrolation Plots - Transport ###

fig=plt.figure()
lag_cor(LAB_index_rm[:-9],ENSO_m[:-9],40,'m','ENSO index')
lag_cor(LAB_index_rm[:-9],NAO_m[:-9],40,'c','NAO index')
lag_cor(LAB_index_rm[:-9],AMOC_transport_50S_m[:-9],40,'r','AMOC 50S')
lag_cor(LAB_index_rm[:-9],AMOC_transport_30S_m[:-9],40,'g','AMOC 30S')
lag_cor(LAB_index_rm[:-9],AMOC_transport_30N_m[:-9],40,'y','AMOC 30N')
lag_cor(LAB_index_rm[:-9],AMOC_transport_50N_m[:-9],40,'b','AMOC 50N')
plt.legend()
plt.show()
plt.title('Peak LAB Convection lagged correlation [AMOC=transport upper 1000m] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
plt.xlabel('Years lag', fontsize=18)
plt.ylabel('Correlation coefficient', fontsize=18)    
plt.show()
plt.legend(shadow=True, prop={'size': 20})
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_LABconvection_AMOC_trans_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
lag_cor(WS_index_rm[:-9],ENSO_m[:-9],40,'m','ENSO index')
lag_cor(WS_index_rm[:-9],NAO_m[:-9],40,'c','NAO index')
lag_cor(WS_index_rm[:-9],AMOC_transport_50S_m[:-9],40,'r','AMOC 50S')
lag_cor(WS_index_rm[:-9],AMOC_transport_30S_m[:-9],40,'g','AMOC 30S')
lag_cor(WS_index_rm[:-9],AMOC_transport_30N_m[:-9],40,'y','AMOC 30N')
lag_cor(WS_index_rm[:-9],AMOC_transport_50N_m[:-9],40,'b','AMOC 50N')
plt.legend()
plt.show()
plt.title('Peak WS Convection lagged correlation [AMOC=transport upper 1000m] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
plt.xlabel('Years lag', fontsize=18)
plt.ylabel('Correlation coefficient', fontsize=18)    
plt.show()
plt.legend(shadow=True, prop={'size': 20})
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_WSconvection_AMOC_trans_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

####################################################################
### AMOC/SMOC Lagged Corrolation Plots - Max/Min Stream function ###

fig=plt.figure()
lag_cor(LAB_index_rm[:-9],ENSO_m[:-9],40,'m','ENSO index')
lag_cor(LAB_index_rm[:-9],NAO_m[:-9],40,'c','NAO index')
lag_cor(LAB_index_rm[:-9],AMOC_max_50S_m[:-9],40,'r','AMOC max 50S')
lag_cor(LAB_index_rm[:-9],AMOC_max_30S_m[:-9],40,'g','AMOC max 30S')
lag_cor(LAB_index_rm[:-9],AMOC_max_30N_m[:-9],40,'y','AMOC max 30N')
lag_cor(LAB_index_rm[:-9],AMOC_max_50N_m[:-9],40,'b','AMOC max 50N')
plt.legend()
plt.show()
plt.title('Peak LAB Convection lagged correlation [AMOC=max streamfunction] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
plt.xlabel('Years lag', fontsize=18)
plt.ylabel('Correlation coefficient', fontsize=18)    
plt.show()
plt.legend(shadow=True, prop={'size': 15})
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_LABconvection_AMOC_max_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
lag_cor(WS_index_rm[:-9],ENSO_m[:-9],40,'m','ENSO index')
lag_cor(WS_index_rm[:-9],NAO_m[:-9],40,'c','NAO index')
lag_cor(WS_index_rm[:-9],AMOC_max_50S_m[:-9],40,'r','AMOC max 50S')
lag_cor(WS_index_rm[:-9],AMOC_max_30S_m[:-9],40,'g','AMOC max 30S')
lag_cor(WS_index_rm[:-9],AMOC_max_30N_m[:-9],40,'y','AMOC max 30N')
lag_cor(WS_index_rm[:-9],AMOC_max_50N_m[:-9],40,'b','AMOC max 50N')
plt.legend()
plt.show()
plt.title('Peak WS Convection lagged correlation [AMOC=max streamfunction] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
plt.xlabel('Years lag', fontsize=18)
plt.ylabel('Correlation coefficient', fontsize=18)    
plt.show()
plt.legend(shadow=True, prop={'size': 15})
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_WSconvection_AMOC_max_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
lag_cor(LAB_index_rm[:-9],ENSO_m[:-9],40,'m','ENSO index')
lag_cor(LAB_index_rm[:-9],NAO_m[:-9],40,'c','NAO index')
lag_cor(LAB_index_rm[:-9],SMOC_min_50S_m[:-9],40,'r','SMOC min 50S')
lag_cor(LAB_index_rm[:-9],SMOC_min_30S_m[:-9],40,'g','SMOC min 30S')
lag_cor(LAB_index_rm[:-9],SMOC_min_30N_m[:-9],40,'y','SMOC min 30N')
lag_cor(LAB_index_rm[:-9],SMOC_min_50N_m[:-9],40,'b','SMOC min 50N')
plt.legend()
plt.show()
plt.title('Peak LAB Convection lagged correlation [SMOC=min streamfunction] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
plt.xlabel('Years lag', fontsize=18)
plt.ylabel('Correlation coefficient', fontsize=18)    
plt.show()
plt.legend(shadow=True, prop={'size': 15})
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_LABconvection_SMOC_min_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
lag_cor(WS_index_rm[:-9],ENSO_m[:-9],40,'m','ENSO index')
lag_cor(WS_index_rm[:-9],NAO_m[:-9],40,'c','NAO index')
lag_cor(WS_index_rm[:-9],SMOC_min_50S_m[:-9],40,'r','SMOC min 50S')
lag_cor(WS_index_rm[:-9],SMOC_min_30S_m[:-9],40,'g','SMOC min 30S')
lag_cor(WS_index_rm[:-9],SMOC_min_30N_m[:-9],40,'y','SMOC min 30N')
lag_cor(WS_index_rm[:-9],SMOC_min_50N_m[:-9],40,'b','SMOC min 50N')
plt.legend()
plt.show()
plt.title('Peak WS Convection lagged correlation [SMOC=min streamfunction] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
plt.xlabel('Years lag', fontsize=18)
plt.ylabel('Correlation coefficient', fontsize=18)    
plt.show()
plt.legend(shadow=True, prop={'size': 15})
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_WSconvection_SMOC_min_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


fig=plt.figure()
lag_cor(LAB_index_rm[:-9],ENSO_m[:-9],40,'m','ENSO index')
lag_cor(LAB_index_rm[:-9],NAO_m[:-9],40,'c','NAO index')
lag_cor(LAB_index_rm[:-9],WS_index_rm[:-9],40,'y','WS index')
lag_cor(LAB_index_rm[:-9], np.nanmean(North_Westerlies, axis=1)[:-9],40,'b','North Westerlies')
lag_cor(LAB_index_rm[:-9], np.nanmean(North_Trades, axis=1)[:-9],40,'r','North Trades')
lag_cor(LAB_index_rm[:-9],AMOC_max_30N_m[:-9],40,'slategrey','AMOC max 30N')
lag_cor(LAB_index_rm[:-9],AMOC_max_50N_m[:-9],40,'g','AMOC max 50N')
plt.legend()
plt.show()
plt.title('Peak LAB Convection lagged correlation [AMOC=max streamfunction] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
plt.xlabel('Years lag', fontsize=18)
plt.ylabel('Correlation coefficient', fontsize=18)    
plt.show()
plt.legend(shadow=True, prop={'size': 15})
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_LABconvection_AMOC_max_2.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
lag_cor(WS_index_rm[:-9],ENSO_m[:-9],40,'m','ENSO index')
lag_cor(WS_index_rm[:-9],NAO_m[:-9],40,'c','NAO index')
lag_cor(WS_index_rm[:-9],LAB_index_rm[:-9],40,'y','LAB index')
lag_cor(WS_index_rm[:-9], np.nanmean(South_Westerlies, axis=1)[:-9],40,'b','South Westerlies')
lag_cor(WS_index_rm[:-9], np.nanmean(South_Trades, axis=1)[:-9],40,'r','South Trades')
lag_cor(WS_index_rm[:-9],SMOC_min_30S_m[:-9],40,'slategrey','SMOC min 30S')
lag_cor(WS_index_rm[:-9],SMOC_min_50S_m[:-9],40,'g','SMOC min 50S')
plt.legend()
plt.show()
plt.title('Peak WS Convection lagged correlation [SMOC=min streamfunction] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
plt.xlabel('Years lag', fontsize=18)
plt.ylabel('Correlation coefficient', fontsize=18)    
plt.show()
plt.legend(shadow=True, prop={'size': 15})
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_WSconvection_SMOC_min_2.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


fig=plt.figure()
lag_cor(LAB_index_rm[:-9],WS_index_rm[:-9],40,'c','WS index')
lag_cor(LAB_index_rm[:-9], np.nanmean(North_Westerlies, axis=1)[:-9],40,'m','North Westerlies')
lag_cor(LAB_index_rm[:-9],AMOC_max_30N_m[:-9],40,'r','AMOC max 30N')
lag_cor(LAB_index_rm[:-9],AMOC_max_50N_m[:-9],40,'g','AMOC max 50N')
lag_cor(LAB_index_rm[:-9],AMOC_max_30S_m[:-9],40,'b','AMOC max 30S')
plt.legend()
plt.show()
plt.title('Peak LAB Convection lagged correlation [AMOC=max streamfunction] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
plt.xlabel('Years lag', fontsize=18)
plt.ylabel('Correlation coefficient', fontsize=18)    
plt.show()
plt.legend(shadow=True, prop={'size': 15})
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_LABconvection_AMOC_max_3.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
lag_cor(WS_index_rm[:-9],LAB_index_rm[:-9],40,'c','LAB index')
lag_cor(WS_index_rm[:-9], np.nanmean(South_Westerlies, axis=1)[:-9],40,'m','South Westerlies')
lag_cor(WS_index_rm[:-9],SMOC_min_30S_m[:-9],40,'r','SMOC min 30S')
lag_cor(WS_index_rm[:-9],SMOC_min_50S_m[:-9],40,'g','SMOC min 50S')
plt.legend()
plt.show()
plt.title('Peak WS Convection lagged correlation [SMOC=min streamfunction] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
plt.xlabel('Years lag', fontsize=18)
plt.ylabel('Correlation coefficient', fontsize=18)    
plt.show()
plt.legend(shadow=True, prop={'size': 15})
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_WSconvection_SMOC_min_3.png', format='png', dpi=300, transparent=True, bbox_inches='tight')





