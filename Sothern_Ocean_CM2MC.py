### IMPORTANT: The following path needs to be appended so that the neccessary codes can be read from that directory ###
import sys
sys.path.append('/data1/home/basadieh/behzadcodes/behzadlibrary')
#from Behzadlib import * ### This is loaded from '/data1/home/basadieh/behzadcodes/behzadlibrary'
#from CMIP5lib import *  ### This is loaded from '/data1/home/basadieh/behzadcodes/behzadlibrary'
from Behzadlib import func_MLD, func_time_depth_plot, func_stream, func_ENSO, func_NAO, func_latlon_regrid, func_latlon_regrid_eq, func_oceanlandmask, func_oceanindex, func_regrid, func_plot, func_plot_save, func_plot_bounds_save
from BehzadlibPlot import func_plotmap_contourf
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

GCM2 = 'GFDL-CM2MC'
year_start2=4001
year_end2=5000
###############################################################################
###############################################################################
GCM = 'GFDL-ESM2G'
year_start=1
year_end=500

### Regrdridding calculations ###
# creating new coordinate grid, same which was used in interpolation in data processing code
lat_n_regrid, lon_n_regrid = 180, 360 # Number of Lat and Lon elements in the regridded data
lon_min_regrid, lon_max_regrid = 0, 360 # Min and Max value of Lon in the regridded data
lat_min_regrid, lat_max_regrid = -90, 90 # Min and Max value of Lat in the regridded data

Lat_regrid_1D, Lon_regrid_1D, Lat_bound_regrid, Lon_bound_regrid = func_latlon_regrid_eq(lat_n_regrid, lon_n_regrid, lat_min_regrid, lat_max_regrid, lon_min_regrid, lon_max_regrid)
Lon_regrid_2D, Lat_regrid_2D = np.meshgrid(Lon_regrid_1D, Lat_regrid_1D)

Ocean_Land_mask = func_oceanlandmask(Lat_regrid_2D, Lon_regrid_2D) # 1= ocean, 0= land

Ocean_Index = func_oceanindex (Lat_regrid_2D, Lon_regrid_2D) # [0=land] [2=Pacific] [3=Indian Ocean] [6=Atlantic] [10=Arctic] [8=Baffin Bay (west of Greenland)] [9=Norwegian Sea (east of Greenland)] [11=Hudson Bay (Canada)] 

###############################################################################
dir_data_in1 = (dir_pwd + '/Results_CM2MC/data/') # Directory to raed raw data from
###############################################################################
from numpy import loadtxt

#file_name_dir=(dir_data_in1+ 'composites.dat'); composites = loadtxt(file_name_dir)
file_name_dir=(dir_data_in1+ 'depths.dat'); depths = loadtxt(file_name_dir)
file_name_dir=(dir_data_in1+ 'latitudes.dat'); latitudes = loadtxt(file_name_dir)
file_name_dir=(dir_data_in1+ 'longitudes.dat'); longitudes = loadtxt(file_name_dir)

#file_name_dir=(dir_data_in1+ 'lag_labsst_amoc30n.dat'); lag_labsst_amoc30n = loadtxt(file_name_dir)
#file_name_dir=(dir_data_in1+ 'lag_labsst_amoc30s.dat'); lag_labsst_amoc30s = loadtxt(file_name_dir)
#file_name_dir=(dir_data_in1+ 'lag_labsst_amoc50n.dat'); lag_labsst_amoc50n = loadtxt(file_name_dir)
#file_name_dir=(dir_data_in1+ 'lag_labsst_wedsst.dat'); lag_labsst_wedsst = loadtxt(file_name_dir)
#file_name_dir=(dir_data_in1+ 'lag_wedsst_amoc50s.dat'); lag_wedsst_amoc50s = loadtxt(file_name_dir)
#file_name_dir=(dir_data_in1+ 'lag_wedsst_labsst.dat'); lag_wedsst_labsst = loadtxt(file_name_dir)
#file_name_dir=(dir_data_in1+ 'lag_wedsst_smoc30s.dat'); lag_wedsst_smoc30s = loadtxt(file_name_dir)
#file_name_dir=(dir_data_in1+ 'lag_wedsst_smoc50s.dat'); lag_wedsst_smoc50s = loadtxt(file_name_dir)
#file_name_dir=(dir_data_in1+ 'lag_wedsst_winds_55s_to_45s_zonal_avg.dat'); lag_wedsst_winds_55s_to_45s_zonal_avg = loadtxt(file_name_dir)
#file_name_dir=(dir_data_in1+ 'lag_wedsst_smoc50s.dat'); lag_wedsst_smoc50s = loadtxt(file_name_dir)

file_name_dir=(dir_data_in1+ 'air_sst_composite_labrador_convection.dat'); air_sst_composite_labrador_convection = loadtxt(file_name_dir)
file_name_dir=(dir_data_in1+ 'air_sst_composite_weddell_convection.dat'); air_sst_composite_weddell_convection = loadtxt(file_name_dir)

file_name_dir=(dir_data_in1+ 'sst_composite_labrador_convection.dat'); sst_composite_labrador_convection = loadtxt(file_name_dir)
file_name_dir=(dir_data_in1+ 'sst_composite_weddell_convection.dat'); sst_composite_weddell_convection = loadtxt(file_name_dir)

file_name_dir=(dir_data_in1+ 'lag_labrador_convection_amoc_30n.dat'); lag_labrador_convection_amoc_30n = loadtxt(file_name_dir)
file_name_dir=(dir_data_in1+ 'lag_labrador_convection_amoc_30s.dat'); lag_labrador_convection_amoc_30s = loadtxt(file_name_dir)
file_name_dir=(dir_data_in1+ 'lag_labrador_convection_amoc_50n.dat'); lag_labrador_convection_amoc_50n = loadtxt(file_name_dir)

file_name_dir=(dir_data_in1+ 'lag_labrador_convection_weddell_convection.dat'); lag_labrador_convection_weddell_convection = loadtxt(file_name_dir)
file_name_dir=(dir_data_in1+ 'lag_weddell_convection_labrador_convection.dat'); lag_weddell_convection_labrador_convection = loadtxt(file_name_dir)

file_name_dir=(dir_data_in1+ 'lag_weddell_convection_amoc50s.dat'); lag_weddell_convection_amoc50s = loadtxt(file_name_dir)
file_name_dir=(dir_data_in1+ 'lag_weddell_convection_smoc30s.dat'); lag_weddell_convection_smoc30s = loadtxt(file_name_dir)
file_name_dir=(dir_data_in1+ 'lag_weddell_convection_smoc50s.dat'); lag_weddell_convection_smoc50s = loadtxt(file_name_dir)
file_name_dir=(dir_data_in1+ 'lag_weddell_convection_winds_55s_to_45s.dat'); lag_weddell_convection_winds_55s_to_45s = loadtxt(file_name_dir)

file_name_dir=(dir_data_in1+ 'lag2dcorr_labrador_convection_amoc.dat'); lag2dcorr_labrador_convection_amoc = loadtxt(file_name_dir)
file_name_dir=(dir_data_in1+ 'lag2dcorr_weddell_convection_amoc.dat'); lag2dcorr_weddell_convection_amoc = loadtxt(file_name_dir)
file_name_dir=(dir_data_in1+ 'lag2dslope_labrador_convection_amoc.dat'); lag2dslope_labrador_convection_amoc = loadtxt(file_name_dir)
file_name_dir=(dir_data_in1+ 'lag2dslope_weddell_convection_amoc.dat'); lag2dslope_weddell_convection_amoc = loadtxt(file_name_dir)

file_name_dir=(dir_data_in1+ 'potential_density_labrador_average.dat'); potential_density_labrador_average = loadtxt(file_name_dir)
file_name_dir=(dir_data_in1+ 'potential_density_labrador_convection.dat'); potential_density_labrador_convection = loadtxt(file_name_dir)
file_name_dir=(dir_data_in1+ 'potential_density_labrador_nonconvection.dat'); potential_density_labrador_nonconvection = loadtxt(file_name_dir)

file_name_dir=(dir_data_in1+ 'potential_density_weddell_average.dat'); potential_density_weddell_average = loadtxt(file_name_dir)
file_name_dir=(dir_data_in1+ 'potential_density_weddell_convection.dat'); potential_density_weddell_convection = loadtxt(file_name_dir)
file_name_dir=(dir_data_in1+ 'potential_density_weddell_nonconvection.dat'); potential_density_weddell_nonconvection = loadtxt(file_name_dir)

file_name_dir=(dir_data_in1+ 'salinity_labrador_average.dat'); salinity_labrador_average = loadtxt(file_name_dir)
file_name_dir=(dir_data_in1+ 'salinity_labrador_convection.dat'); salinity_labrador_convection = loadtxt(file_name_dir)
file_name_dir=(dir_data_in1+ 'salinity_labrador_nonconvection.dat'); salinity_labrador_nonconvection = loadtxt(file_name_dir)

file_name_dir=(dir_data_in1+ 'salinity_weddell_average.dat'); salinity_weddell_average = loadtxt(file_name_dir)
file_name_dir=(dir_data_in1+ 'salinity_weddell_convection.dat'); salinity_weddell_convection = loadtxt(file_name_dir)
file_name_dir=(dir_data_in1+ 'salinity_weddell_nonconvection.dat'); salinity_weddell_nonconvection = loadtxt(file_name_dir)

file_name_dir=(dir_data_in1+ 'temperature_labrador_average.dat'); temperature_labrador_average = loadtxt(file_name_dir)
file_name_dir=(dir_data_in1+ 'temperature_labrador_convection.dat'); temperature_labrador_convection = loadtxt(file_name_dir)
file_name_dir=(dir_data_in1+ 'temperature_labrador_nonconvection.dat'); temperature_labrador_nonconvection = loadtxt(file_name_dir)

file_name_dir=(dir_data_in1+ 'temperature_weddell_average.dat'); temperature_weddell_average = loadtxt(file_name_dir)
file_name_dir=(dir_data_in1+ 'temperature_weddell_convection.dat'); temperature_weddell_convection = loadtxt(file_name_dir)
file_name_dir=(dir_data_in1+ 'temperature_weddell_nonconvection.dat'); temperature_weddell_nonconvection = loadtxt(file_name_dir)

file_name_dir=(dir_data_in1+ 'tau_curl_lab_comp.dat'); tau_curl_lab_comp = loadtxt(file_name_dir)
file_name_dir=(dir_data_in1+ 'tau_curl_lab_conv.dat'); tau_curl_lab_conv = loadtxt(file_name_dir)
file_name_dir=(dir_data_in1+ 'tau_curl_lab_ncnv.dat'); tau_curl_lab_ncnv = loadtxt(file_name_dir)

file_name_dir=(dir_data_in1+ 'tau_curl_wed_comp.dat'); tau_curl_wed_comp = loadtxt(file_name_dir)
file_name_dir=(dir_data_in1+ 'tau_curl_wed_conv.dat'); tau_curl_wed_conv = loadtxt(file_name_dir)
file_name_dir=(dir_data_in1+ 'tau_curl_wed_ncnv.dat'); tau_curl_wed_ncnv = loadtxt(file_name_dir)

file_name_dir=(dir_data_in1+ 'tau_x_lab_comp.dat'); tau_x_lab_comp = loadtxt(file_name_dir)
file_name_dir=(dir_data_in1+ 'tau_x_lab_conv.dat'); tau_x_lab_conv = loadtxt(file_name_dir)
file_name_dir=(dir_data_in1+ 'tau_x_lab_ncnv.dat'); tau_x_lab_ncnv = loadtxt(file_name_dir)

file_name_dir=(dir_data_in1+ 'tau_x_wed_comp.dat'); tau_x_wed_comp = loadtxt(file_name_dir)
file_name_dir=(dir_data_in1+ 'tau_x_wed_conv.dat'); tau_x_wed_conv = loadtxt(file_name_dir)
file_name_dir=(dir_data_in1+ 'tau_x_wed_ncnv.dat'); tau_x_wed_ncnv = loadtxt(file_name_dir)

file_name_dir=(dir_data_in1+ 'sst_avg.dat'); sst_avg = loadtxt(file_name_dir)
file_name_dir=(dir_data_in1+ 'tau_x_avg.dat'); tau_x_avg = loadtxt(file_name_dir)



###############################################################################
#file_name_dir=(dir_data_in1+ 'windstress.nc')
#dset_t = Dataset(file_name_dir)
#dset_t.variables


#fig=plt.figure()
#im1=plt.plot(lag_labsst_amoc30n[:,0], lag_labsst_amoc30n[:,1], 'b', label='AMOC max 30N')
#plt.grid(True,which="both",ls="-", color='0.65')
#plt.legend()
#plt.show()
#plt.title('Peak LAB Convection lagged correlation [AMOC=max streamfunction] - '+str(year_start2)+'-'+str(year_end2)+' - '+str(GCM2), fontsize=18)
#plt.xlabel('Years lag', fontsize=18)
#plt.ylabel('Correlation coefficient', fontsize=18)    
#plt.show()
#plt.legend(shadow=True, prop={'size': 15})
#mng = plt.get_current_fig_manager()
#mng.window.showMaximized() # Maximizes the plot window to save figures in full





####################################################################
#%%####           Westerlies and Trades Calculations        ########
####################################################################
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
Tau_X = func_regrid(dset_tauu_yearly, dset_tauu.y, dset_tauu.x, Lat_regrid_2D, Lon_regrid_2D)

Winds_40S60S_0W30W=Tau_X[:,30:51,330:]
Winds_40S60S_0W30W=np.nanmean( Winds_40S60S_0W30W, axis=2)
Winds_40S60S_0W30W=np.nanmean( Winds_40S60S_0W30W, axis=1)
Winds_40S60S_0W30W_m=runningMeanFast(Winds_40S60S_0W30W, 10)

Winds_50N65N_20W65W=Tau_X[:,140:156,295:341]
Winds_50N65N_20W65W=np.nanmean( Winds_50N65N_20W65W, axis=2)
Winds_50N65N_20W65W=np.nanmean( Winds_50N65N_20W65W, axis=1)
Winds_50N65N_20W65W_m=runningMeanFast(Winds_50N65N_20W65W, 10)

#################################################
#%%####           Air temperature        ########
#################################################
t_frequency='Amon'
variable='ts'# Air temperature

dir_data_in1 = ('/data5/scratch/CMIP5/CMIP5_models/atmosphere_physics/') # Directory to raed raw data from
dir_data_in2=(dir_data_in1+ GCM + '/piControl/mo/')

dset_ts = netcdf_read (dir_data_in2+str(variable)+'_'+str(t_frequency)+'_'+str(GCM)+'*12.nc',str(GCM), variable) # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
start_date_i,end_date_i = dset_ts.find_time(dset_ts.times, year_start, year_end)

dset_ts_extracted=dset_ts.extract_data(dset_ts.variable,start_date_i,end_date_i)
#here we in order to calculate yearly data from monthly data we first reshape array
dset_ts_yearly=dset_ts_extracted.reshape(( (year_end-year_start+1) ,12,len(dset_ts_extracted[0]),len(dset_ts_extracted[0][0])))
#second - we calculate mean over each 12 elements, which is second dimension of array
dset_ts_yearly=np.mean(dset_ts_yearly,axis=1)
Air_Temp = func_regrid(dset_ts_yearly, dset_ts.y, dset_ts.x, Lat_regrid_2D, Lon_regrid_2D)
Air_Temp=np.nanmean(Air_Temp, axis=0)

##############################################################################
#################        To restore:        ##################################
import os
import shelve

dir_pwd = os.getcwd() # Gets the current directory (and in which the code is placed)
filename_out = (dir_pwd + '/Results/AllResults_'+GCM+'_500yr.out') # Directory to save processed data

var_list_short=np.load('var_list_short.npy')

my_shelf = shelve.open(filename_out)
#globals()[key]=my_shelf[key]
globals()['WS_index_rm']=my_shelf['WS_index_rm']
globals()['LAB_index_rm']=my_shelf['LAB_index_rm']
globals()['Zonal_winds']=my_shelf['Zonal_winds']
globals()['AMOC_max']=my_shelf['AMOC_max']
globals()['composite_LAB']=my_shelf['composite_LAB']
globals()['composite_WS']=my_shelf['composite_WS']
#my_shelf.close()

North_Westerlies=Zonal_winds[:,135:145]

AMOC_max_50S= AMOC_max[:,40] # AMOC at 50S # Max of streamfunction method
AMOC_max_50S_m=runningMeanFast(AMOC_max_50S, 10)
AMOC_max_30S= AMOC_max[:,60] # AMOC at 30S # Max of streamfunction method
AMOC_max_30S_m=runningMeanFast(AMOC_max_30S, 10)
AMOC_max_50N= AMOC_max[:,140] # AMOC at 50N # Max of streamfunction method
AMOC_max_50N_m=runningMeanFast(AMOC_max_50N, 10)
AMOC_max_30N= AMOC_max[:,120] # AMOC at 30N # Max of streamfunction method
AMOC_max_30N_m=runningMeanFast(AMOC_max_30N, 10)


##### GFDL-ESM2G Plots ######
dir_figs = (dir_pwd + '/Figures_ESM2G/') # Directory to save figures

fig=plt.figure()
lag_cor(LAB_index_rm[:-9],WS_index_rm[:-9],40,'y','WS index')
lag_cor(LAB_index_rm[:-9], Winds_50N65N_20W65W_m[:-9],40,'g','Zonal Winds [50N65N, 20W65W]')
lag_cor(LAB_index_rm[:-9],AMOC_max_50N_m[:-9],40,'b','AMOC max 50N')
plt.legend()
plt.show()
plt.xticks(fontsize = 20); plt.yticks(fontsize = 24)
plt.ylim(-0.7,0.7)
plt.title('Peak LAB Convection lagged correlation - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
plt.xlabel('Years lag', fontsize=18)
plt.ylabel('Correlation coefficient', fontsize=18)    
plt.show()
plt.legend(shadow=True, prop={'size': 15})
l = plt.axhline(y=0, color='k')
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_LABconvection_Final_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
lag_cor(WS_index_rm[:-9],LAB_index_rm[:-9],40,'y','LAB index')
lag_cor(WS_index_rm[:-9], Winds_40S60S_0W30W_m[:-9],40,'g','Zonal Winds [40S60S, 0W30W]')
lag_cor(WS_index_rm[:-9],AMOC_max_50S_m[:-9],40,'b','AMOC max 50S')
lag_cor(WS_index_rm[:-9],AMOC_max_30S_m[:-9],40,'r','AMOC max 30S') 
plt.legend()
plt.show()
plt.xticks(fontsize = 20); plt.yticks(fontsize = 24)
plt.ylim(-0.7,0.7)
plt.title('Peak WS Convection lagged correlation - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
plt.xlabel('Years lag', fontsize=18)
plt.ylabel('Correlation coefficient', fontsize=18)    
plt.show()
plt.legend(shadow=True, prop={'size': 15})
l = plt.axhline(y=0, color='k')
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_WSconvection_Final_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

#### Composits Plost ####

Plot_Var = composite_LAB
#cmap_limit=np.nanmax(np.abs( np.nanpercentile(Plot_Var, 99)))
Plot_range=np.linspace(-2,2,101)
Plot_unit='(°C)'; Plot_title= 'Surface air temperature composites - LAB - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig, m = func_plotmap_contourf(Plot_Var, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.seismic, 'cyl', 210., 80., -80., '-')
fig.savefig(dir_figs+str(GCM)+'_surface_temperature_composites_LAB.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

Plot_Var = composite_WS
#cmap_limit=np.nanmax(np.abs( np.nanpercentile(Plot_Var, 99)))
Plot_range=np.linspace(-1.6,1.6,101)
Plot_unit='(°C)'; Plot_title= 'Surface air temperature composites - WS - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig, m = func_plotmap_contourf(Plot_Var, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.seismic, 'cyl', 210., 80., -80., '-')
fig.savefig(dir_figs+str(GCM)+'_surface_temperature_composites_WS.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

Plot_Var = Air_Temp-273.15
#cmap_limit=np.nanmax(np.abs( np.nanpercentile(Plot_Var, 99)))
Plot_range=np.linspace(-2,34,68)
Plot_unit='(°C)'; Plot_title= 'Surface air temperature [°C] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig, m = func_plotmap_contourf(Plot_Var, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.jet, 'cyl', 210., 80., -80., 'fill')
fig.savefig(dir_figs+str(GCM)+'_air_temperature.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


Var_plot_unit='Unit = °C'
fig=plt.figure()
m = Basemap( projection='mill',lon_0=210, llcrnrlon=30.,llcrnrlat=-80.,urcrnrlon=390.,urcrnrlat=80.)
m.fillcontinents(color='0.95')
m.drawmapboundary(fill_color='0.9')
m.drawmapboundary(fill_color='#000099')
m.drawcoastlines(linewidth=1.0, linestyle='solid', color='k', antialiased=1, ax=None, zorder=None)
m.drawparallels(np.arange(-90.,90.001,30.),labels=[True,False,False,False], linewidth=0.01) # labels = [left,right,top,bottom]
m.drawmeridians(np.arange(0.,360.,30.),labels=[False,False,False,True], linewidth=0.01) # labels = [left,right,top,bottom]
#cmap_limit=np.nanmax(np.abs( np.nanpercentile(sst_avg_regrid, 99.9)))
#levels=np.linspace(-cmap_limit,cmap_limit,200)
levels=np.linspace(-2,34,68)
im=m.contourf(Lon_regrid_2D,Lat_regrid_2D,Air_Temp-273.15 ,levels,latlon=True, cmap=plt.cm.jet)
cbar = m.colorbar(im,"right", size="3%", pad="2%", extend='max') # extend='both' will extend the colorbar in both sides (upper side and down side)
cbar.set_label(Var_plot_unit)
cbar.ax.tick_params(labelsize=18) 
plt.show()
plt.title('Surface air temperature [°C] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_air_temperature.png', format='png', dpi=300, transparent=True, bbox_inches='tight')






##### GFDL-CM2MC Plots ######
dir_figs = (dir_pwd + '/Figures_CM2MC/') # Directory to save figures

fig=plt.figure()
im1=plt.plot(lag_labrador_convection_weddell_convection[:,0], lag_labrador_convection_weddell_convection[:,1], 'y', label='WS index')
im1=plt.plot(lag_labrador_convection_amoc_50n[39:119,0], lag_labrador_convection_amoc_50n[39:119,1], 'b', label='AMOC max 50N')
#im1=plt.plot(lag_labrador_convection_amoc_50n[39:119,0], lag_labrador_convection_amoc_50n[39:119,1], 'g', label='Zonal Winds [40S60S, 0W30W]')
plt.grid(True,which="both",ls="-", color='0.65')
plt.legend()
plt.show()
plt.xticks(fontsize = 20); plt.yticks(fontsize = 24)
plt.ylim(-0.7,0.7)
plt.title('Peak LAB Convection lagged correlation - '+str(year_start2)+'-'+str(year_end2)+' - '+str(GCM2), fontsize=18)
plt.xlabel('Years lag', fontsize=18)
plt.ylabel('Correlation coefficient', fontsize=18)    
plt.show()
plt.legend(shadow=True, prop={'size': 15})
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM2)+'_lagcor_LABconvection_Final_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
im1=plt.plot(lag_weddell_convection_labrador_convection[:,0], lag_weddell_convection_labrador_convection[:,1], 'y', label='LAB index')
im1=plt.plot(lag_weddell_convection_amoc50s[39:119,0], lag_weddell_convection_amoc50s[39:119,1], 'b', label='AMOC max 50S')
im1=plt.plot(lag_weddell_convection_winds_55s_to_45s[39:119,0], lag_weddell_convection_winds_55s_to_45s[39:119,1], 'g', label='Zonal Winds [45S-55S]')
plt.grid(True,which="both",ls="-", color='0.65')
plt.legend()
plt.show()
plt.xticks(fontsize = 20); plt.yticks(fontsize = 24)
plt.ylim(-0.7,0.7)
plt.title('Peak WS Convection lagged correlation - '+str(year_start2)+'-'+str(year_end2)+' - '+str(GCM2), fontsize=18)
plt.xlabel('Years lag', fontsize=18)
plt.ylabel('Correlation coefficient', fontsize=18)    
plt.show()
plt.legend(shadow=True, prop={'size': 15})
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM2)+'_lagcor_WSconvection_Final_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


### Regridding ###

air_sst_composite_labrador_convection_regrid = func_regrid(air_sst_composite_labrador_convection, latitudes, longitudes, Lat_regrid_2D, Lon_regrid_2D)
air_sst_composite_weddell_convection_regrid = func_regrid(air_sst_composite_weddell_convection, latitudes, longitudes, Lat_regrid_2D, Lon_regrid_2D)

sst_avg_mat=sst_avg.reshape((latitudes.shape[0], longitudes.shape[0] )); sst_avg_mat [ sst_avg_mat < -1e10 ] = nan
sst_avg_regrid = func_regrid(sst_avg_mat, latitudes, longitudes, Lat_regrid_2D, Lon_regrid_2D)

tau_x_avg_mat=tau_x_avg.reshape((latitudes.shape[0], longitudes.shape[0] )); tau_x_avg_mat [ tau_x_avg_mat < -1e10 ] = nan
tau_x_avg_regrid = func_regrid(tau_x_avg_mat, latitudes, longitudes, Lat_regrid_2D, Lon_regrid_2D)


tau_curl_lab_comp_mat=tau_curl_lab_comp.reshape((latitudes.shape[0], longitudes.shape[0] )); tau_curl_lab_comp_mat [ tau_curl_lab_comp_mat < -1e10 ] = nan
tau_curl_lab_comp_regrid = func_regrid(tau_curl_lab_comp_mat, latitudes, longitudes, Lat_regrid_2D, Lon_regrid_2D)
tau_curl_lab_conv_mat=tau_curl_lab_conv.reshape((latitudes.shape[0], longitudes.shape[0] )); tau_curl_lab_conv_mat [ tau_curl_lab_conv_mat < -1e10 ] = nan
tau_curl_lab_conv_regrid = func_regrid(tau_curl_lab_conv_mat, latitudes, longitudes, Lat_regrid_2D, Lon_regrid_2D)
tau_curl_lab_ncnv_mat=tau_curl_lab_ncnv.reshape((latitudes.shape[0], longitudes.shape[0] )); tau_curl_lab_ncnv_mat [ tau_curl_lab_ncnv_mat < -1e10 ] = nan
tau_curl_lab_ncnv_regrid = func_regrid(tau_curl_lab_ncnv_mat, latitudes, longitudes, Lat_regrid_2D, Lon_regrid_2D)

tau_curl_wed_comp_mat=tau_curl_wed_comp.reshape((latitudes.shape[0], longitudes.shape[0] )); tau_curl_wed_comp_mat [ tau_curl_wed_comp_mat < -1e10 ] = nan
tau_curl_wed_comp_regrid = func_regrid(tau_curl_wed_comp_mat, latitudes, longitudes, Lat_regrid_2D, Lon_regrid_2D)
tau_curl_wed_conv_mat=tau_curl_wed_conv.reshape((latitudes.shape[0], longitudes.shape[0] )); tau_curl_wed_conv_mat [ tau_curl_wed_conv_mat < -1e10 ] = nan
tau_curl_wed_conv_regrid = func_regrid(tau_curl_wed_conv_mat, latitudes, longitudes, Lat_regrid_2D, Lon_regrid_2D)
tau_curl_wed_ncnv_mat=tau_curl_wed_ncnv.reshape((latitudes.shape[0], longitudes.shape[0] )); tau_curl_wed_ncnv_mat [ tau_curl_wed_ncnv_mat < -1e10 ] = nan
tau_curl_wed_ncnv_regrid = func_regrid(tau_curl_wed_ncnv_mat, latitudes, longitudes, Lat_regrid_2D, Lon_regrid_2D)

tau_curl_lab_comp_regrid = func_regrid(tau_curl_lab_comp, latitudes, longitudes, Lat_regrid_2D, Lon_regrid_2D)
tau_curl_lab_conv_regrid = func_regrid(tau_curl_lab_conv, latitudes, longitudes, Lat_regrid_2D, Lon_regrid_2D)
tau_curl_lab_ncnv_regrid = func_regrid(tau_curl_lab_ncnv, latitudes, longitudes, Lat_regrid_2D, Lon_regrid_2D)
tau_curl_wed_comp_regrid = func_regrid(tau_curl_wed_comp, latitudes, longitudes, Lat_regrid_2D, Lon_regrid_2D)
tau_curl_wed_conv_regrid = func_regrid(tau_curl_wed_conv, latitudes, longitudes, Lat_regrid_2D, Lon_regrid_2D)
tau_curl_wed_ncnv_regrid = func_regrid(tau_curl_wed_ncnv, latitudes, longitudes, Lat_regrid_2D, Lon_regrid_2D)

tau_x_lab_comp_regrid = func_regrid(tau_x_lab_comp, latitudes, longitudes, Lat_regrid_2D, Lon_regrid_2D)
tau_x_lab_conv_regrid = func_regrid(tau_x_lab_conv, latitudes, longitudes, Lat_regrid_2D, Lon_regrid_2D)
tau_x_lab_ncnv_regrid = func_regrid(tau_x_lab_ncnv, latitudes, longitudes, Lat_regrid_2D, Lon_regrid_2D)
tau_x_wed_comp_regrid = func_regrid(tau_x_wed_comp, latitudes, longitudes, Lat_regrid_2D, Lon_regrid_2D)
tau_x_wed_conv_regrid = func_regrid(tau_x_wed_conv, latitudes, longitudes, Lat_regrid_2D, Lon_regrid_2D)
tau_x_wed_ncnv_regrid = func_regrid(tau_x_wed_ncnv, latitudes, longitudes, Lat_regrid_2D, Lon_regrid_2D)


#### Composits Plost ####

Plot_Var = air_sst_composite_labrador_convection_regrid
#cmap_limit=np.nanmax(np.abs( np.nanpercentile(Plot_Var, 99)))
Plot_range=np.linspace(-2,2,101)
Plot_unit='(°C)'; Plot_title= 'Surface air temperature composites - LAB - '+str(year_start2)+'-'+str(year_end2)+' - '+str(GCM2)
fig, m = func_plotmap_contourf(Plot_Var, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.seismic, 'cyl', 210., 80., -80., 'fill')
fig.savefig(dir_figs+str(GCM2)+'_surface_temperature_composites_LAB.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

Plot_Var = air_sst_composite_weddell_convection_regrid
#cmap_limit=np.nanmax(np.abs( np.nanpercentile(Plot_Var, 99)))
Plot_range=np.linspace(-1.6,1.6,101)
Plot_unit='(°C)'; Plot_title= 'Surface air temperature composites - WS - '+str(year_start2)+'-'+str(year_end2)+' - '+str(GCM2)
fig, m = func_plotmap_contourf(Plot_Var, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.seismic, 'cyl', 210., 80., -80., 'fill')
fig.savefig(dir_figs+str(GCM2)+'_surface_temperature_composites_WS.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

Plot_Var = sst_avg_regrid
#cmap_limit=np.nanmax(np.abs( np.nanpercentile(Plot_Var, 99)))
Plot_range=np.linspace(-2,34,68)
Plot_unit='(°C)'; Plot_title= 'Surface air temperature [°C] '+str(year_start2)+'-'+str(year_end2)+' - '+str(GCM2)
fig, m = func_plotmap_contourf(Plot_Var, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.jet, 'cyl', 210., 80., -80., 'fill')
fig.savefig(dir_figs+str(GCM2)+'_air_temperature.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

Plot_Var = np.subtract(sst_avg_regrid, Air_Temp-273.15)
#cmap_limit=np.nanmax(np.abs( np.nanpercentile(Plot_Var, 99)))
Plot_range=np.linspace(-8,8,101)
Plot_unit='(°C)'; Plot_title= 'Surface air temperature [°C] - CM2MC minus ESM2G'
fig, m = func_plotmap_contourf(Plot_Var, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.seismic, 'cyl', 210., 80., -80., 'fill')
fig.savefig(dir_figs+'air_temperature_CM2MCminusESM2G.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


#####################
### Wind Plots ###
Plot_Var = tau_x_avg_regrid
Plot_Var_l=copy.deepcopy(Plot_Var)  ; Plot_Var_l[ Ocean_Land_mask==0 ]=nan
#cmap_limit=np.nanmax(np.abs( np.nanpercentile(Plot_Var, 99)))
Plot_range=np.linspace(-0.412,0.412,27)
Plot_unit='(N/m2)'; Plot_title= 'Wind Stress X (N/m2) - '+str(year_start2)+'-'+str(year_end2)+' - '+str(GCM2)
fig, m = func_plotmap_contourf(Plot_Var, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.seismic, 'cyl', 210., 80., -80., 'fill')
im3=m.contour(Lon_regrid_2D[0:70,:], Lat_regrid_2D[0:70,:],Plot_Var_l[0:70,:], levels = [0], latlon=True, colors='darkgreen')
plt.clabel(im3, fontsize=8, inline=1)
fig.savefig(dir_figs+'Wind_Tau_X_'+str(GCM2)+'.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


###############################################################################
########################       Density Plots       ############################
fig=plt.figure()
im1=plt.plot(potential_density_labrador_average[0:21], depths[0:21])
plt.gca().invert_yaxis()
plt.xlim(26.8,28)
plt.xticks(fontsize = 20); plt.yticks(fontsize = 20)
plt.xlabel('Potential Density [kg/m3]', fontsize=18)
plt.ylabel('Depth [m]', fontsize=18)
plt.title('Potential Density Ave - Labrador Sea - '+str(year_start2)+'-'+str(year_end2)+' - '+str(GCM2), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM2)+'_density_depth_LAB.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
im1=plt.plot(potential_density_labrador_convection[0:21], depths[0:21], 'r')
im1=plt.plot(potential_density_labrador_nonconvection[0:21], depths[0:21], 'b')
plt.gca().invert_yaxis()
plt.xlim(26.8,28)
plt.xticks(fontsize = 20); plt.yticks(fontsize = 20)
plt.xlabel('Potential Density [kg/m3]', fontsize=18)
plt.ylabel('Depth [m]', fontsize=18)
plt.title('Potential Density Ave - Labrador Sea (red=Lab. Sea conv, blue=nonconv) - '+str(year_start2)+'-'+str(year_end2)+' - '+str(GCM2), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM2)+'_density_depth_LAB_conv.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
im1=plt.plot(potential_density_weddell_average[0:21], depths[0:21])
plt.gca().invert_yaxis()
plt.xlim(27.1,27.8)
plt.xticks(fontsize = 20); plt.yticks(fontsize = 20)
plt.xlabel('Potential Density [kg/m3]', fontsize=18)
plt.ylabel('Depth [m]', fontsize=18)
plt.title('Potential Density Ave - Weddell Sea - '+str(year_start2)+'-'+str(year_end2)+' - '+str(GCM2), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM2)+'_density_depth_WS.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
im1=plt.plot(potential_density_weddell_convection[0:21], depths[0:21], 'r')
im1=plt.plot(potential_density_weddell_nonconvection[0:21], depths[0:21], 'b')
plt.gca().invert_yaxis()
plt.xlim(27.1,27.8)
plt.xticks(fontsize = 20); plt.yticks(fontsize = 20)
plt.xlabel('Potential Density [kg/m3]', fontsize=18)
plt.ylabel('Depth [m]', fontsize=18)
plt.title('Potential Density Ave - Weddell Sea (red=Wedd. Sea conv, blue=nonconv) - '+str(year_start2)+'-'+str(year_end2)+' - '+str(GCM2), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM2)+'_density_depth_WS_conv.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


###############################################################################
######################       Temperature Plots       ##########################
fig=plt.figure()
im1=plt.plot(temperature_labrador_average[0:21], depths[0:21])
plt.gca().invert_yaxis()
#plt.xlim(4,6.5)
plt.xticks(fontsize = 20); plt.yticks(fontsize = 20)
plt.xlabel('Potential Temperature [C]', fontsize=18)
plt.ylabel('Depth [m]', fontsize=18)
plt.title('Potential Temperature Ave - Labrador Sea - '+str(year_start2)+'-'+str(year_end2)+' - '+str(GCM2), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM2)+'_temp_depth_LAB.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
im1=plt.plot(temperature_labrador_convection[0:21], depths[0:21], 'r')
im1=plt.plot(temperature_labrador_nonconvection[0:21], depths[0:21], 'b')
plt.gca().invert_yaxis()
#plt.xlim(4,6.5)
plt.xticks(fontsize = 20); plt.yticks(fontsize = 20)
plt.xlabel('Potential Temperature [C]', fontsize=18)
plt.ylabel('Depth [m]', fontsize=18)
plt.title('Potential Temperature Ave - Labrador Sea (red=Lab. Sea conv, blue=nonconv) - '+str(year_start2)+'-'+str(year_end2)+' - '+str(GCM2), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM2)+'_temp_depth_LAB_conv.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
im1=plt.plot(temperature_weddell_average[0:21], depths[0:21])
plt.gca().invert_yaxis()
#plt.xlim(4,6.5)
plt.xticks(fontsize = 20); plt.yticks(fontsize = 20)
plt.xlabel('Potential Temperature [C]', fontsize=18)
plt.ylabel('Depth [m]', fontsize=18)
plt.title('Potential Temperature Ave - Weddell Sea - '+str(year_start2)+'-'+str(year_end2)+' - '+str(GCM2), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM2)+'_temp_depth_WS.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
im1=plt.plot(temperature_weddell_convection[0:21], depths[0:21], 'r')
im1=plt.plot(temperature_weddell_nonconvection[0:21], depths[0:21], 'b')
plt.gca().invert_yaxis()
#plt.xlim(4,6.5)
plt.xticks(fontsize = 20); plt.yticks(fontsize = 20)
plt.xlabel('Potential Temperature [C]', fontsize=18)
plt.ylabel('Depth [m]', fontsize=18)
plt.title('Potential Temperature Ave - Weddell Sea (red=Wedd. Sea conv, blue=nonconv) - '+str(year_start2)+'-'+str(year_end2)+' - '+str(GCM2), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM2)+'_temp_depth_WS_conv.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

###############################################################################
#######################       Salinity Plots       ###########################
fig=plt.figure()
im1=plt.plot(salinity_labrador_average[0:21], depths[0:21])
plt.gca().invert_yaxis()
plt.xlim(34.2,35.3)
plt.xticks(fontsize = 20); plt.yticks(fontsize = 20)
plt.xlabel('Salinity [kg/m3]', fontsize=18)
plt.ylabel('Depth [m]', fontsize=18)
plt.title('Salinity Ave - Labrador Sea - '+str(year_start2)+'-'+str(year_end2)+' - '+str(GCM2), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM2)+'_sal_depth_LAB.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
im1=plt.plot(salinity_labrador_convection[0:21], depths[0:21], 'r')
im1=plt.plot(salinity_labrador_nonconvection[0:21], depths[0:21], 'b')
plt.gca().invert_yaxis()
plt.xlim(34.2,35.3)
plt.xticks(fontsize = 20); plt.yticks(fontsize = 20)
plt.xlabel('Salinity [kg/m3]', fontsize=18)
plt.ylabel('Depth [m]', fontsize=18)
plt.title('Salinity Ave - Labrador Sea (red=Lab. Sea conv, blue=nonconv) - '+str(year_start2)+'-'+str(year_end2)+' - '+str(GCM2), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM2)+'_sal_depth_LAB_conv.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
im1=plt.plot(salinity_weddell_average[0:21], depths[0:21])
plt.gca().invert_yaxis()
plt.xlim(33.95,34.65)
plt.xticks(fontsize = 20); plt.yticks(fontsize = 20)
plt.xlabel('Salinity [kg/m3]', fontsize=18)
plt.ylabel('Depth [m]', fontsize=18)
plt.title('Salinity Ave - Weddell Sea - '+str(year_start2)+'-'+str(year_end2)+' - '+str(GCM2), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM2)+'_sal_depth_WS.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
im1=plt.plot(salinity_weddell_convection[0:21], depths[0:21], 'r')
im1=plt.plot(salinity_weddell_nonconvection[0:21], depths[0:21], 'b')
plt.gca().invert_yaxis()
plt.xlim(33.8,34.65)
plt.xticks(fontsize = 20); plt.yticks(fontsize = 20)
plt.xlabel('Salinity [kg/m3]', fontsize=18)
plt.ylabel('Depth [m]', fontsize=18)
plt.title('Salinity Ave - Weddell Sea (red=Wedd. Sea conv, blue=nonconv) - '+str(year_start2)+'-'+str(year_end2)+' - '+str(GCM2), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM2)+'_sal_depth_WS_conv.png', format='png', dpi=300, transparent=True, bbox_inches='tight')














