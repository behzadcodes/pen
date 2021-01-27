### IMPORTANT: The following path needs to be appended so that the neccessary codes can be read from that directory ###
import sys
sys.path.append('/data1/home/basadieh/behzadcodes/behzadlibrary')
#from Behzadlib import * ### This is loaded from '/data1/home/basadieh/behzadcodes/behzadlibrary'
#from CMIP5lib import *  ### This is loaded from '/data1/home/basadieh/behzadcodes/behzadlibrary'
from Behzadlib import func_MLD, func_time_depth_plot, func_stream, func_ENSO, func_NAO, func_latlon_regrid, func_latlon_regrid_eq, func_oceanlandmask, func_oceanindex, func_regrid, func_plot, func_plot_save, func_plot_bounds_save, func_EOF
from BehzadlibPlot import func_plot2Dcolor_contourf, func_plotline_1var, func_plotline_2var, func_plot2Dcontour_1var, func_plot2Dcontour_2var, func_plotmap_contourf
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
from scipy import stats
########################################

dir_pwd = os.getcwd() # Gets the current directory (and in which the code is placed)
dir_figs = (dir_pwd + '/Figures1/') # Directory to save figures
dir_results = (dir_pwd + '/Results/') # Directory to save results

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

##############################################################################
#################        To restore:        ##################################
import os
import shelve

dir_pwd = os.getcwd() # Gets the current directory (and in which the code is placed)
filename_out = (dir_pwd + '/Results/AllResults_'+GCM+'_500yr.out') # Directory to save processed data

my_shelf = shelve.open(filename_out)
#globals()[key]=my_shelf[key]
globals()['Depths']=my_shelf['Depths']
my_shelf.close()

##############################################################################
##############################################################################
Density_allyears=np.load('Density_allyears_GFDL-ESM2G_500yr.npy')

Density_allyears_WS = copy.deepcopy(Density_allyears[:,:,20:35,320:]) # [40W-0W,70S-55S]
Density_allyears_WS_ave = np.nanmean( np.nanmean (Density_allyears_WS, axis=3), axis=2)

Delta_rho_WS = ( (Density_allyears_WS_ave[:,30]*53.04 + Density_allyears_WS_ave[:,29]*24.17 ) / (53.04+24.17) ) - Density_allyears_WS_ave[:,0] # Rho_500 minus Rho_0 # Interpolating between Depth[30]=524.17 and Depth[29]=446.94 to calculate density at depth 500m
Delta_rho_WS_CoeffVar = np.std(Delta_rho_WS) / np.nanmean(Delta_rho_WS)

R_value_depth_WS = empty((Density_allyears_WS_ave.shape[1])) * nan
Slope_depth_WS = empty((Density_allyears_WS_ave.shape[1])) * nan
for ii in range(0,R_value_depth_WS.shape[0]):
    slope_ij, intercept_ij, r_value_ij, p_value_ij, std_err_ij = stats.linregress(Delta_rho_WS, Density_allyears_WS_ave[:,ii])
    R_value_depth_WS[ii] = r_value_ij
    Slope_depth_WS[ii] = slope_ij


Density_allyears_WS_conS=np.load('Density_allyears_70W_0W_75S_50S_WS_conS.npy') ## [70W-0W,75S-50S] ## [15:41,290:360]
Density_allyears_WS_conT=np.load('Density_allyears_70W_0W_75S_50S_WS_conT.npy') ## [70W-0W,75S-50S] ## [15:41,290:360]
Density_allyears_WS_conS=Density_allyears_WS_conS [:,:,5:20,30:] # selecting [40W-0W,70S-55S] box from the original [70W-0W,75S-50S] box
Density_allyears_WS_conT=Density_allyears_WS_conT [:,:,5:20,30:] # selecting [40W-0W,70S-55S] box from the original [70W-0W,75S-50S] box

Density_allyears_WS_conS_ave = np.nanmean( np.nanmean (Density_allyears_WS_conS, axis=3), axis=2)
Density_allyears_WS_conT_ave = np.nanmean( np.nanmean (Density_allyears_WS_conT, axis=3), axis=2)

Delta_rho_WS_conS = ( (Density_allyears_WS_conS_ave[:,30]*53.04 + Density_allyears_WS_conS_ave[:,29]*24.17 ) / (53.04+24.17) ) - Density_allyears_WS_conS_ave[:,0] # Rho_500 minus Rho_0 # Interpolating between Depth[30]=524.17 and Depth[29]=446.94 to calculate density at depth 500m
Delta_rho_WS_conT = ( (Density_allyears_WS_conT_ave[:,30]*53.04 + Density_allyears_WS_conT_ave[:,29]*24.17 ) / (53.04+24.17) ) - Density_allyears_WS_conT_ave[:,0] # Rho_500 minus Rho_0 # Interpolating between Depth[30]=524.17 and Depth[29]=446.94 to calculate density at depth 500m

R_value_depth_WS_conS = empty((Density_allyears_WS_conS_ave.shape[1])) * nan
Slope_depth_WS_conS = empty((Density_allyears_WS_conS_ave.shape[1])) * nan
R_value_depth_WS_conT = empty((Density_allyears_WS_conT_ave.shape[1])) * nan
Slope_depth_WS_conT = empty((Density_allyears_WS_conT_ave.shape[1])) * nan
for ii in range(0,R_value_depth_WS_conS.shape[0]):
    slope_ij, intercept_ij, r_value_ij, p_value_ij, std_err_ij = stats.linregress(Delta_rho_WS, Density_allyears_WS_conS_ave[:,ii])
    R_value_depth_WS_conS[ii] = r_value_ij
    Slope_depth_WS_conS[ii] = slope_ij
    slope_ij, intercept_ij, r_value_ij, p_value_ij, std_err_ij = stats.linregress(Delta_rho_WS, Density_allyears_WS_conT_ave[:,ii])
    R_value_depth_WS_conT[ii] = r_value_ij
    Slope_depth_WS_conT[ii] = slope_ij    
    

P_Var_x=Slope_depth_WS[0:30]  ;P_Var_y=Depths[0:30]
P_xlable='Regression Slope'  ;P_ylable='Depth [m]'
P_title='Regress of Rho_500 minus Rho_0 vs Rho - Weddell Sea [40W-0W,70S-55S] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plotline_1var(P_Var_x, P_Var_y, P_xlable, P_ylable, P_title, 'darkcyan', 'invert_yaxis')
#fig.savefig(dir_figs+str(GCM)+'_density_regress_depth_WS_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


P_Var_x1=Slope_depth_WS[0:31] ;  P_Var_y=Depths[0:31]
P_Var_x2=Slope_depth_WS_conT[0:31] ; P_Var_x3=Slope_depth_WS_conS[0:31] 
P_xlable='Regression Slope'  ;P_ylable='Depth [m]'
P_title='Regress of (Rho_500 - Rho_0) vs. (Rho at each drpth) - Weddell Sea [40W-0W,70S-55S] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)

fig=plt.figure()
im1=plt.plot(P_Var_x1, P_Var_y, c='k', label='Density', linewidth=3.0)
im2=plt.plot(P_Var_x2, P_Var_y, c='r', label='Density - (varying S, constant T)', linewidth=3.0)
im3=plt.plot(P_Var_x3, P_Var_y, c='g', label='Density - (varying T, constant S)', linewidth=3.0)
plt.gca().invert_yaxis()
plt.legend(prop={'size': 24}, loc='lower left', fancybox=True, framealpha=0.8)
plt.xticks(fontsize = 26); plt.yticks(fontsize = 26)
plt.xlabel(P_xlable, fontsize=18) ; plt.ylabel(P_ylable, fontsize=18)
plt.title(P_title, fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full

fig.savefig(dir_figs+str(GCM)+'_density_regress_constant_S_T_depth_WS_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


##############################################################################

Delta_rho = ( (Density_allyears[:,30,:,:]*53.04 + Density_allyears[:,29,:,:]*24.17 ) / (53.04+24.17) ) - Density_allyears[:,0,:,:] # Rho_500 minus Rho_0 # Interpolating between Depth[30]=524.17 and Depth[29]=446.94 to calculate density at depth 500m
Delta_rho_CoeffVar = np.std(Delta_rho ,axis=0) / np.nanmean(Delta_rho ,axis=0)


Plot_Var = copy.deepcopy(Delta_rho_CoeffVar)
#Plot_Var[ Ocean_Land_mask==0 ]=nan # masking over land, so grid cells that fall on land area (value=0) will be deleted
Plot_range=51 #np.linspace(1018,1028,41)
Plot_unit='[ - ]'; Plot_title= 'Coefficient of variation (st.dev / mean) in Stratification [i.e. Rho_500 minus Rho_0] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)

fig, m = func_plotmap_contourf(Plot_Var, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.jet, 'cyl', 210., 80., -80., '-')
plt.show()
#fig.savefig(dir_figs+str(GCM)+'_Coeff_of_var_R500minusR0_Map.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


#################################################################################
## Weddell Sea [70W-0W,75S-50S] 

t_frequency='Omon'
variable_thetao='thetao' # Sea Water Temperature
variable_so='so' # Sea Water Salinity

dir_data_in1 = ('/data5/scratch/CMIP5/CMIP5_models/ocean_physics/') # Directory to raed raw data from
dir_data_in2=(dir_data_in1+ GCM + '/piControl/mo/')

dset_thetao = netcdf_read (dir_data_in2+str(variable_thetao)+'_'+str(t_frequency)+'_'+str(GCM)+'*12.nc',str(GCM), variable_thetao) # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
dset_so = netcdf_read (dir_data_in2+str(variable_so)+'_'+str(t_frequency)+'_'+str(GCM)+'*12.nc',str(GCM), variable_so)
Depths=dset_thetao.lvl[:]

import seawater as sw
start_date_i,end_date_i = dset_thetao.find_time(dset_thetao.times, year_start, year_end) # start_date_i=0 , end_date_i=5999

Lat_orig=dset_thetao.y
Lon_orig=dset_thetao.x

#### Density with constant temperature - varying salinity

#data_thetao_extracted_ws_con=[]
#for t in range(int((end_date_i+1-start_date_i)/12)):
#    print('Temp - Year: ', year_start+t)
#    data_thetao_extracted=dset_thetao.extract_data(dset_thetao.variable,start_date_i+12*t+1-1,start_date_i+12*t+12-1)
#    data_thetao_extracted=np.nanmean(data_thetao_extracted, axis=0)
#    data_thetao_extracted=np.squeeze(data_thetao_extracted)
#    data_thetao_extracted = data_thetao_extracted[:,0:34,210:280] # [70W-0W,75S-50S]
#
#    data_thetao_extracted_ws_con.append(data_thetao_extracted)   
#
#data_thetao_extracted_ws_con=np.asarray(data_thetao_extracted_ws_con)
#data_thetao_extracted_ws_con=np.nanmean(data_thetao_extracted_ws_con, axis=0)

Temp_allyears=np.load('Temp_allyears_GFDL-ESM2G_500yr.npy')
data_thetao_extracted_ws_con=np.nanmean(Temp_allyears[:,:,0:34,210:280], axis=0)
#np.save('data_thetao_extracted_ws_con.npy',data_thetao_extracted_ws_con)

Density_allyears_WS_conT=[]
for t in range(int((end_date_i+1-start_date_i)/12)):
    print('Density calc - Year: ', year_start+t)
    data_so_extracted=dset_so.extract_data(dset_so.variable,start_date_i+12*t+1-1,start_date_i+12*t+12-1)
    data_so_extracted=np.nanmean(data_so_extracted, axis=0)
    data_so_extracted=np.squeeze(data_so_extracted)
     
    data_so_extracted = data_so_extracted [:,0:34,210:280] # [70W-0W,75S-50S]
    
    data_dens=sw.dens0(data_so_extracted, data_thetao_extracted_ws_con-273.15)

    data_i = func_regrid(data_dens, Lat_orig[0:34,210:280], Lon_orig[0:34,210:280], Lat_regrid_2D[15:41,290:360], Lon_regrid_2D[15:41,290:360])
    data_i[data_i>100000]=np.nan
    Density_allyears_WS_conT.append(data_i)   
#    data_dens=np.asarray(data_dens)
#    Density_allyears_WS_conT.append(data_dens)

Density_allyears_WS_conT=np.asarray(Density_allyears_WS_conT)

np.save('Density_allyears_GFDL-ESM2G_500yr_WS_conT.npy',Density_allyears_WS_conT) ## [70W-0W,75S-50S] ## [15:41,290:360]

Density_allyears_WS_conT=np.load('Density_allyears_70W_0W_75S_50S_WS_conT.npy') ## [70W-0W,75S-50S] ## [15:41,290:360]


#### Density with constant salinity - varying stemperature

#data_so_extracted_ws_con=[]
#for t in range(int((end_date_i+1-start_date_i)/12)):
#    print('Temp - Year: ', year_start+t)
#    data_so_extracted=dset_so.extract_data(dset_so.variable,start_date_i+12*t+1-1,start_date_i+12*t+12-1)
#    data_so_extracted=np.nanmean(data_so_extracted, axis=0)
#    data_so_extracted=np.squeeze(data_so_extracted)
#    data_so_extracted = data_so_extracted[:,0:34,210:280] # [70W-0W,75S-50S]
#
#    data_so_extracted_ws_con.append(data_so_extracted)   
#
#data_so_extracted_ws_con=np.asarray(data_so_extracted_ws_con)
#data_so_extracted_ws_con=np.nanmean(data_so_extracted_ws_con, axis=0)

Salinity_allyears=np.load('Salinity_allyears_GFDL-ESM2G_500yr.npy')
data_so_extracted_ws_con=np.nanmean(Salinity_allyears[:,:,0:34,210:280], axis=0)
#np.save('data_so_extracted_ws_con.npy',data_so_extracted_ws_con)

Density_allyears_WS_conS=[]
for t in range(int((end_date_i+1-start_date_i)/12)):
    print('Density calc - Year: ', year_start+t)
    data_thetao_extracted=dset_thetao.extract_data(dset_thetao.variable,start_date_i+12*t+1-1,start_date_i+12*t+12-1)
    data_thetao_extracted=np.nanmean(data_thetao_extracted, axis=0)
    data_thetao_extracted=np.squeeze(data_thetao_extracted)
     
    data_thetao_extracted = data_thetao_extracted [:,0:34,210:280] # [70W-0W,75S-50S]
    
    data_dens=sw.dens0(data_so_extracted_ws_con, data_thetao_extracted-273.15)

    data_i = func_regrid(data_dens, Lat_orig[0:34,210:280], Lon_orig[0:34,210:280], Lat_regrid_2D[15:41,290:360], Lon_regrid_2D[15:41,290:360])
    data_i[data_i>100000]=np.nan
    Density_allyears_WS_conS.append(data_i)   
#    data_dens=np.asarray(data_dens)
#    Density_allyears_WS_conS.append(data_dens)

Density_allyears_WS_conS=np.asarray(Density_allyears_WS_conS)

np.save('Density_allyears_GFDL-ESM2G_500yr_WS_conS.npy',Density_allyears_WS_conS) ## [70W-0W,75S-50S] ## [15:41,290:360]

Density_allyears_WS_conS=np.load('Density_allyears_70W_0W_75S_50S_WS_conS.npy') ## [70W-0W,75S-50S] ## [15:41,290:360]

########################################################################
























