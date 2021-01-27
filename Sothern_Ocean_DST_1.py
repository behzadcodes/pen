### IMPORTANT: The following path needs to be appended so that the neccessary codes can be read from that directory ###
import sys
sys.path.append('/data1/home/basadieh/behzadcodes/behzadlibrary')
#from Behzadlib import * ### This is loaded from '/data1/home/basadieh/behzadcodes/behzadlibrary'
#from CMIP5lib import *  ### This is loaded from '/data1/home/basadieh/behzadcodes/behzadlibrary'
from Behzadlib import func_MLD, func_time_depth_plot, func_stream, func_ENSO, func_NAO, func_latlon_regrid, func_latlon_regrid_eq, func_oceanlandmask, func_oceanindex, func_regrid, func_plot, func_plot_save, func_plot_bounds_save
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
########################################

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

##############################################################################
#################        To restore:        ##################################
#import os
#import shelve
#
#dir_pwd = os.getcwd() # Gets the current directory (and in which the code is placed)
#filename_out = (dir_pwd + '/Results/AllResults_'+GCM+'_500yr.out') # Directory to save processed data
#
#my_shelf = shelve.open(filename_out)
#for key in my_shelf:
#    globals()[key]=my_shelf[key]
#my_shelf.close()
##############################################################################
##############################################################################
import os
import shelve

dir_pwd = os.getcwd() # Gets the current directory (and in which the code is placed)
filename_out = (dir_pwd + '/Results/AllResults_'+GCM+'_500yr.out') # Directory to save processed data

my_shelf = shelve.open(filename_out)
#globals()[key]=my_shelf[key]
globals()['WS_index_norm_rm']=my_shelf['WS_index_norm_rm']
globals()['LAB_index_norm_rm']=my_shelf['LAB_index_norm_rm']
globals()['Depths']=my_shelf['Depths']
my_shelf.close()

R_ii_WS = np.where(WS_index_norm_rm <-0.5)
R_jj_WS = np.where(WS_index_norm_rm >0.5)

R_ii_LAB = np.where(LAB_index_norm_rm <-0.5)
R_jj_LAB = np.where(LAB_index_norm_rm >0.5)

##############################################################################
##############################################################################

conv_index_depth_ws=500 # args[7] in func_time_depth_plot code - depth for Convection Index
conv_index_depth_lab=500

t_frequency='Omon'
variable_thetao='thetao' # Sea Water Temperature
variable_so='so' # Sea Water Salinity

dir_data_in1 = ('/data5/scratch/CMIP5/CMIP5_models/ocean_physics/') # Directory to raed raw data from
dir_data_in2=(dir_data_in1+ GCM + '/piControl/mo/')

dset_thetao = netcdf_read (dir_data_in2+str(variable_thetao)+'_'+str(t_frequency)+'_'+str(GCM)+'*12.nc',str(GCM), variable_thetao) # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
dset_so = netcdf_read (dir_data_in2+str(variable_so)+'_'+str(t_frequency)+'_'+str(GCM)+'*12.nc',str(GCM), variable_so)
Depths=dset_thetao.lvl[:]

R_ii_WS = np.where(WS_index_norm_rm <-0.5)
R_jj_WS = np.where(WS_index_norm_rm >0.5)

import seawater as sw
start_date_i,end_date_i = dset_thetao.find_time(dset_thetao.times, year_start, year_end)

Density_allyears=[]

for t in range(int((end_date_i+1-start_date_i)/12)):
    print('Density calc - Year: ', year_start+t)
    data_thetao_extracted=dset_thetao.extract_data(dset_thetao.variable,start_date_i+12*t+1-1,start_date_i+12*t+12-1)
    data_thetao_extracted=np.nanmean(data_thetao_extracted, axis=0)
    data_so_extracted=dset_so.extract_data(dset_so.variable,start_date_i+12*t+1-1,start_date_i+12*t+12-1)
    data_so_extracted=np.nanmean(data_so_extracted, axis=0)
    data_thetao_extracted=np.squeeze(data_thetao_extracted)
    data_so_extracted=np.squeeze(data_so_extracted)
    data_dens=sw.dens0(data_so_extracted, data_thetao_extracted-273.15)

    data_i = func_regrid(data_dens, dset_thetao.y, dset_thetao.x, Lat_regrid_2D, Lon_regrid_2D)
    data_i[data_i>100000]=np.nan
    
    Density_allyears.append(data_i)

Density_allyears=np.asarray(Density_allyears)
Density_allyears_ws_c=Density_allyears[R_ii_WS]
Density_allyears_ws_n=Density_allyears[R_jj_WS]
R_conv_WS_All=np.nanmean(Density_allyears_ws_c,axis=0)
R_nonconv_WS_All=np.nanmean(Density_allyears_ws_n,axis=0)

R_conv_WS_60W0E=R_conv_WS_All[:,:,300:]
R_nonconv_WS_60W0E=R_nonconv_WS_All[:,:,300:]

R_conv_WS_All=np.nanmean(R_conv_WS_All,axis=2)
R_nonconv_WS_All=np.nanmean(R_nonconv_WS_All,axis=2)
R_conv_WS_60W0E=np.nanmean(R_conv_WS_60W0E,axis=2)
R_nonconv_WS_60W0E=np.nanmean(R_nonconv_WS_60W0E,axis=2)

np.save('Density_allyears_GFDL-ESM2G_500yr.npy',Density_allyears)
#Density_allyears=np.load('Density_allyears_GFDL-ESM2G_500yr.npy')


P_Var_x1=P_Var_x2=Lat_regrid_1D[20:36]
P_Var_y1=P_Var_y2=Depths[0:31]
P_Var_z1=R_conv_WS_All[0:31,20:36]-1000
P_Var_z2=R_nonconv_WS_All[0:31,20:36]-1000
P_xlable='Latitude'  ;P_ylable='Depth [m]'
P_title='Potential Density Contours (red=conv, blue=nonconv) - All Southern Ocean - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plot2Dcontour_2var(P_Var_x1, P_Var_y1, P_Var_z1, P_Var_x2, P_Var_y2, P_Var_z2, P_xlable, P_ylable, P_title, 'r', 'b', 10, 'invert_yaxis')
#fig.savefig(dir_figs+str(GCM)+'_density_composites_WS_All.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

P_Var_x1=P_Var_x2=Lat_regrid_1D[20:36]
P_Var_y1=P_Var_y2=Depths[0:31]
P_Var_z1=R_conv_WS_60W0E[0:31,20:36]-1000
P_Var_z2=R_nonconv_WS_60W0E[0:31,20:36]-1000
P_xlable='Latitude'  ;P_ylable='Depth [m]'
P_title='Potential Density Contours (red=conv, blue=nonconv) - 60W-0E - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plot2Dcontour_2var(P_Var_x1, P_Var_y1, P_Var_z1, P_Var_x2, P_Var_y2, P_Var_z2, P_xlable, P_ylable, P_title, 'r', 'b', 10, 'invert_yaxis')
#fig.savefig(dir_figs+str(GCM)+'_density_composites_WS_60W0E.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

P_Var_x1=P_Var_x2=Lat_regrid_1D[20:30]
P_Var_y1=P_Var_y2=Depths[0:31]
P_Var_z1=R_conv_WS_60W0E[0:31,20:30]-1000
P_Var_z2=R_nonconv_WS_60W0E[0:31,20:30]-1000
P_xlable='Latitude'  ;P_ylable='Depth [m]'
P_title='Potential Density Contours (red=conv, blue=nonconv) - 60W-0E - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plot2Dcontour_2var(P_Var_x1, P_Var_y1, P_Var_z1, P_Var_x2, P_Var_y2, P_Var_z2, P_xlable, P_ylable, P_title, 'r', 'b', 10, 'invert_yaxis')
#fig.savefig(dir_figs+str(GCM)+'_density_composites_WS_60W0E.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


#############################
### Average Density Plots ###

#Density_allyears=np.load('Density_allyears_GFDL-ESM2G_500yr.npy')

P_Var_x1=Lat_regrid_1D
P_Var_y1=Depths[0:35]
P_Var_z1=np.nanmean(Density_allyears,axis=0)
for ii in range(P_Var_z1.shape[0]):
    P_Var_z1[ii,:,:][Ocean_Index!=6]=nan
P_Var_z1=np.nanmean(P_Var_z1,axis=2)
P_Var_z1=P_Var_z1[0:35,:]-1000
P_xlable='Latitude'  ;P_ylable='Depth [m]'
P_title='Potential Density Contours - Atlantic - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plot2Dcontour_1var(P_Var_x1, P_Var_y1, P_Var_z1, P_xlable, P_ylable, P_title, 'b', 40, 'invert_yaxis')
P_Var_z2=copy.deepcopy(P_Var_z1)
P_Var_z2[ P_Var_z2 > 27.301]=nan ; P_Var_z2[ P_Var_z2 < 27.199 ]=nan
#fig = func_plot2Dcontour_1var(P_Var_x1, P_Var_y1, P_Var_z2, P_xlable, P_ylable, P_title, 'r', 40, 'invert_yaxis')
im2=m.contour(P_Var_x1, P_Var_y1, P_Var_z2, latlon=True, colors='darkgreen')
fig.savefig(dir_figs+str(GCM)+'_Density_Atlantic_c_2.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


#P_Var_z2=np.nanmean(Density_allyears,axis=0)
#for ii in range(P_Var_z2.shape[0]):
#    P_Var_z2[ii,:,:][Ocean_Index!=6]=nan
#    
#P_Var_z2[ P_Var_z2 > 1027.301]=nan ; P_Var_z2[ P_Var_z2 < 1027.199 ]=nan
#    
#P_Var_z2=np.nanmean(P_Var_z2,axis=2)
#P_Var_z2=P_Var_z1[0:35,:]-1000
#
#P_Var_z2[ P_Var_z2 > 27.301]=nan ; P_Var_z2[ P_Var_z2 < 27.199 ]=nan
#
#fig = func_plot2Dcontour_1var(P_Var_x1, P_Var_y1, P_Var_z2, P_xlable, P_ylable, P_title, 'b', 40, 'invert_yaxis')


P_Var_x1=Lat_regrid_1D
P_Var_y1=Depths[0:35]
P_Var_z1=np.nanmean(Density_allyears,axis=0)
for ii in range(P_Var_z1.shape[0]):
    P_Var_z1[ii,:,:][Ocean_Index!=6]=nan
P_Var_z1=np.nanmean(P_Var_z1,axis=2)
P_Var_z1=P_Var_z1[0:35,:]-1000
P_Var_z2=copy.deepcopy(P_Var_z1)
P_Var_z2[ P_Var_z2 > 27.301]=nan ; P_Var_z2[ P_Var_z2 < 27.199 ]=nan
P_xlable='Latitude'  ;P_ylable='Depth [m]'
P_title='Potential Density Contours - Atlantic - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plot2Dcontour_2var(P_Var_x1, P_Var_y1, P_Var_z1, P_Var_x1, P_Var_y1, P_Var_z2, P_xlable, P_ylable, P_title, 'r', 'b', 10, 'invert_yaxis')


fig, m = func_plot2Dcontour_1var(P_Var_x1, P_Var_y1, P_Var_z1, P_xlable, P_ylable, P_title, 'b', 40, 'invert_yaxis')
P_Var_z2=copy.deepcopy(P_Var_z1)
P_Var_z2[ P_Var_z2 > 27.301]=nan ; P_Var_z2[ P_Var_z2 < 27.199 ]=nan
#fig = func_plot2Dcontour_1var(P_Var_x1, P_Var_y1, P_Var_z2, P_xlable, P_ylable, P_title, 'r', 40, 'invert_yaxis')
im2=m.contour(P_Var_x1, P_Var_y1, P_Var_z2, latlon=True, colors='darkgreen')
fig.savefig(dir_figs+str(GCM)+'_Density_Atlantic_c_2.png', format='png', dpi=300, transparent=True, bbox_inches='tight')




P_Var_x1=Lat_regrid_1D
P_Var_y1=Depths[0:35]
P_Var_z1=np.nanmean(Density_allyears,axis=0)
for ii in range(P_Var_z1.shape[0]):
    P_Var_z1[ii,:,:][Ocean_Index!=2]=nan
P_Var_z1=np.nanmean(P_Var_z1,axis=2)
P_Var_z1=P_Var_z1[0:35,:]-1000
P_xlable='Latitude'  ;P_ylable='Depth [m]'
P_title='Potential Density Contours - Pacific - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plot2Dcontour_1var(P_Var_x1, P_Var_y1, P_Var_z1, P_xlable, P_ylable, P_title, 'b', 40, 'invert_yaxis')
fig.savefig(dir_figs+str(GCM)+'_Density_Pacific_c_2.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


Plot_Var = np.nanmean(Density_allyears[400:500,0,:,:],axis=0)
Plot_Var[ Ocean_Land_mask==0 ]=nan # masking over land, so grid cells that fall on land area (value=0) will be deleted
Plot_range=np.linspace(1018,1028,41)
Plot_unit='[ kg/m3 ]'; Plot_title= 'Average potential density at surface [kg/m3] - '+str(401)+'-'+str(500)+' - '+str(GCM)

fig, m = func_plotmap_contourf(Plot_Var, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.jet, 'cyl', 210., 80., -80., '-')
plt.show()
#fig.savefig(dir_figs+str(GCM)+'_Density_AveMap_0m.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

Plot_Var = np.nanmean(Density_allyears[400:500,10,:,:],axis=0)
Plot_Var[ Ocean_Land_mask==0 ]=nan # masking over land, so grid cells that fall on land area (value=0) will be deleted
Plot_range=np.linspace(1018,1028,41)
Plot_unit='[ kg/m3 ]'; Plot_title= 'Average potential density at 105m depth [kg/m3] - '+str(401)+'-'+str(500)+' - '+str(GCM)

fig, m = func_plotmap_contourf(Plot_Var, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.jet, 'cyl', 210., 80., -80., '-')
plt.show()
#fig.savefig(dir_figs+str(GCM)+'_Density_AveMap_100m.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

Plot_Var = np.nanmean(Density_allyears[400:500,31,:,:],axis=0)
Plot_Var[ Ocean_Land_mask==0 ]=nan # masking over land, so grid cells that fall on land area (value=0) will be deleted
Plot_range=np.linspace(1018,1028,41)
Plot_unit='[ kg/m3 ]'; Plot_title= 'Average potential density at 600m depth [kg/m3] - '+str(401)+'-'+str(500)+' - '+str(GCM)

fig, m = func_plotmap_contourf(Plot_Var, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.jet, 'cyl', 210., 80., -80., '-')
plt.show()
#fig.savefig(dir_figs+str(GCM)+'_Density_AveMap_1000m.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

##############################
### Average Salinity Plots ###

#Salinity_allyears=np.load('Salinity_allyears_GFDL-ESM2G_500yr.npy')

P_Var_x1=Lat_regrid_1D
P_Var_y1=Depths[0:35]
P_Var_z1=np.nanmean(Salinity_allyears,axis=0)
for ii in range(P_Var_z1.shape[0]):
    P_Var_z1[ii,:,:][Ocean_Index!=6]=nan
P_Var_z1=np.nanmean(P_Var_z1,axis=2)
P_Var_z1=P_Var_z1[0:35,:]
P_xlable='Latitude'  ;P_ylable='Depth [m]'
P_title='Salinity Contours - Atlantic - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plot2Dcontour_1var(P_Var_x1, P_Var_y1, P_Var_z1, P_xlable, P_ylable, P_title, 'r', 30, 'invert_yaxis')
fig.savefig(dir_figs+str(GCM)+'_Salinity_Atlantic_c_3.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

P_Var_x1=Lat_regrid_1D
P_Var_y1=Depths[0:35]
P_Var_z1=np.nanmean(Salinity_allyears,axis=0)
for ii in range(P_Var_z1.shape[0]):
    P_Var_z1[ii,:,:][Ocean_Index!=2]=nan
P_Var_z1=np.nanmean(P_Var_z1,axis=2)
P_Var_z1=P_Var_z1[0:35,:]
P_xlable='Latitude'  ;P_ylable='Depth [m]'
P_title='Salinity Contours - Pacific - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plot2Dcontour_1var(P_Var_x1, P_Var_y1, P_Var_z1, P_xlable, P_ylable, P_title, 'r', 30, 'invert_yaxis')
fig.savefig(dir_figs+str(GCM)+'_Salinity_Pacific_c_3.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


Plot_Var = np.nanmean(Salinity_allyears[400:500,0,:,:],axis=0)
Plot_Var[ Ocean_Land_mask==0 ]=nan # masking over land, so grid cells that fall on land area (value=0) will be deleted
Plot_range=np.linspace(23,38,46)
Plot_unit='[ psu ]'; Plot_title= 'Average salinity at surface [psu] - '+str(401)+'-'+str(500)+' - '+str(GCM)

fig, m = func_plotmap_contourf(Plot_Var, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.jet, 'cyl', 210., 80., -80., '-')
plt.show()
#fig.savefig(dir_figs+str(GCM)+'_Salinity_AveMap_0m.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

Plot_Var = np.nanmean(Salinity_allyears[400:500,10,:,:],axis=0)
Plot_Var[ Ocean_Land_mask==0 ]=nan # masking over land, so grid cells that fall on land area (value=0) will be deleted
Plot_range=np.linspace(23,38,46)
Plot_unit='[ psu ]'; Plot_title= 'Average salinity at 105m depth [psu] - '+str(401)+'-'+str(500)+' - '+str(GCM)

fig, m = func_plotmap_contourf(Plot_Var, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.jet, 'cyl', 210., 80., -80., '-')
plt.show()
#fig.savefig(dir_figs+str(GCM)+'_Salinity_AveMap_100m.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

Plot_Var = np.nanmean(Salinity_allyears[400:500,31,:,:],axis=0)
Plot_Var[ Ocean_Land_mask==0 ]=nan # masking over land, so grid cells that fall on land area (value=0) will be deleted
Plot_range=np.linspace(23,38,46)
Plot_unit='[ psu ]'; Plot_title= 'Average salinity at 600m depth [psu] - '+str(401)+'-'+str(500)+' - '+str(GCM)

fig, m = func_plotmap_contourf(Plot_Var, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.jet, 'cyl', 210., 80., -80., '-')
plt.show()
#fig.savefig(dir_figs+str(GCM)+'_Salinity_AveMap_600m.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

###########################################################
###############    Labrador Sea Plots    ##################
###########################################################

R_ii_LAB = np.where(LAB_index_norm_rm <-0.5)
R_jj_LAB = np.where(LAB_index_norm_rm >0.5)

Density_allyears_lab_c=Density_allyears[R_ii_LAB]
Density_allyears_lab_n=Density_allyears[R_jj_LAB]

Density_allyears_LAB=Density_allyears[:,:,140:156,295:341]
Density_allyears_LAB=np.nanmean(Density_allyears_LAB,axis=2)
Density_allyears_LAB=np.nanmean(Density_allyears_LAB,axis=0)

P_Var_x1=Lon_regrid_1D[295:341]-360
P_Var_y1=Depths[0:31]
P_Var_z1=Density_allyears_LAB[0:31,:]-1000
P_xlable='Longitude'  ;P_ylable='Depth [m]'
P_title='Potential Density Contours - 20W-65W - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plot2Dcontour_1var(P_Var_x1, P_Var_y1, P_Var_z1, P_xlable, P_ylable, P_title, 'b', 40, 'invert_yaxis')
fig.savefig(dir_figs+'Density_'+str(GCM)+'_LAB_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

Density_allyears_LAB_ave=np.nanmean(Density_allyears_LAB,axis=1)

P_Var_x=Density_allyears_LAB_ave[0:40]-1000  ;P_Var_y=Depths[0:40]
P_xlable='Potential Density [kg/m3]'  ;P_ylable='Depth [m]'
P_title='Potential Density Average- 20W-65W , 50N-65N - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plotline_1var(P_Var_x, P_Var_y, P_xlable, P_ylable, P_title, 'darkcyan', 'invert_yaxis')
#fig.savefig(dir_figs+'Density_'+str(GCM)+'_LAB_3.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


Density_allyears_LAB_c=Density_allyears_lab_c[:,:,140:156,295:341]
Density_allyears_LAB_c=np.nanmean(Density_allyears_LAB_c,axis=2)
Density_allyears_LAB_c=np.nanmean(Density_allyears_LAB_c,axis=0)
Density_allyears_LAB_n=Density_allyears_lab_n[:,:,140:156,295:341]
Density_allyears_LAB_n=np.nanmean(Density_allyears_LAB_n,axis=2)
Density_allyears_LAB_n=np.nanmean(Density_allyears_LAB_n,axis=0)

P_Var_x1=P_Var_x2=Lon_regrid_1D[295:341]-360
P_Var_y1=P_Var_y2=Depths[0:40]
P_Var_z1=Density_allyears_LAB_c[0:40,:]-1000
P_Var_z2=Density_allyears_LAB_n[0:40,:]-1000
P_xlable='Longitude'  ;P_ylable='Depth [m]'
P_title='Potential Density Contours - Labrador Sea (50N-65N ave.) (red=conv, blue=nonconv) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plot2Dcontour_2var(P_Var_x1, P_Var_y1, P_Var_z1, P_Var_x2, P_Var_y2, P_Var_z2, P_xlable, P_ylable, P_title, 'r', 'b', 40, 'invert_yaxis')
#fig.savefig(dir_figs+'Density_'+str(GCM)+'_LAB_C1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

P_Var_x1=P_Var_x2=Lon_regrid_1D[295:341]-360
P_Var_y1=P_Var_y2=Depths[0:31]
P_Var_z1=Density_allyears_LAB_c[0:31,:]-1000
P_Var_z2=Density_allyears_LAB_n[0:31,:]-1000
P_xlable='Longitude'  ;P_ylable='Depth [m]'
P_title='Potential Density Contours - Labrador Sea (50N-65N ave.) (red=conv, blue=nonconv) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plot2Dcontour_2var(P_Var_x1, P_Var_y1, P_Var_z1, P_Var_x2, P_Var_y2, P_Var_z2, P_xlable, P_ylable, P_title, 'r', 'b', 40, 'invert_yaxis')
#fig.savefig(dir_figs+'Density_'+str(GCM)+'_LAB_C2.png', format='png', dpi=300, transparent=True, bbox_inches='tight')



P_Var_x1=P_Var_x2=Lon_regrid_1D[295:341]-360
P_Var_y1=P_Var_y2=Depths[0:31]
P_Var_z1=Density_allyears_LAB_n[0:31,:]-1000
P_Var_z2=Density_allyears_LAB_c[0:31,:]-1000
P_range=np.linspace(25,28,11) 
P_xlable='Longitude'  ;P_ylable='Depth [m]';  P_unit = ' - '
P_title='Potential Density profile - Labrador Sea (50N-65N ave.) (Color=nonvonv, Contours=conv) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plot2Dcolor_contourf(P_Var_x1, P_Var_y1, P_Var_z1, P_range, P_xlable, P_ylable, P_title, P_unit, plt.cm.jet, 'invert_yaxis')
im1=plt.contour(P_Var_x2, P_Var_y2, P_Var_z2, P_range, colors='k')
plt.clabel(im1, fontsize=8, inline=1)
#fig.savefig(dir_figs+'Density_'+str(GCM)+'_LAB_C22.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

Density_allyears_LAB_c_ave=np.nanmean(Density_allyears_LAB_c,axis=1)
Density_allyears_LAB_n_ave=np.nanmean(Density_allyears_LAB_n,axis=1)

P_Var_x1=Density_allyears_LAB_c_ave[0:31]-1000
P_Var_x2=Density_allyears_LAB_n_ave[0:31]-1000
P_Var_y1=P_Var_y2=Depths[0:31]
P_xlable='Potential Density [kg/m3]'  ;P_ylable='Depth [m]'
P_title='Potential Density Average (red=conv, blue=nonconv) - 20W-65W , 50N-65N - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plotline_2var(P_Var_x1, P_Var_y1, P_Var_x2, P_Var_y2, P_xlable, P_ylable, P_title, 'r', 'b', 'Lab. Sea conv', 'Lab. Sea nonconv', 'best', 'invert_yaxis')
#fig.savefig(dir_figs+'Density_'+str(GCM)+'_LAB_C4.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

#### Temperature ####

Temp_allyears=[]

for t in range(int((end_date_i+1-start_date_i)/12)):
    print('Temp calc - Year: ', year_start+t)
    data_thetao_extracted=dset_thetao.extract_data(dset_thetao.variable,start_date_i+12*t+1-1,start_date_i+12*t+12-1)
    data_thetao_extracted=np.nanmean(data_thetao_extracted, axis=0)
    data_thetao_extracted=np.squeeze(data_thetao_extracted)

    data_i = func_regrid(data_thetao_extracted, dset_thetao.y, dset_thetao.x, Lat_regrid_2D, Lon_regrid_2D)
    data_i[data_i>100000]=np.nan
    
    Temp_allyears.append(data_i)

Temp_allyears=np.asarray(Temp_allyears)
Temp_allyears_c=Temp_allyears[R_ii_WS]
Temp_allyears_n=Temp_allyears[R_jj_WS]
T_conv_WS_All=np.nanmean(Temp_allyears_c,axis=0)
T_nonconv_WS_All=np.nanmean(Temp_allyears_n,axis=0)

T_conv_WS_60W0E=T_conv_WS_All[:,:,300:]
T_nonconv_WS_60W0E=T_nonconv_WS_All[:,:,300:]

T_conv_WS_All=np.nanmean(T_conv_WS_All,axis=2)
T_nonconv_WS_All=np.nanmean(T_nonconv_WS_All,axis=2)
T_conv_WS_60W0E=np.nanmean(T_conv_WS_60W0E,axis=2)
T_nonconv_WS_60W0E=np.nanmean(T_nonconv_WS_60W0E,axis=2)

np.save('Temp_allyears_GFDL-ESM2G_500yr.npy',Temp_allyears)
#Temp_allyears=np.load('Temp_allyears_GFDL-ESM2G_500yr.npy')

P_Var_x1=P_Var_x2=Lat_regrid_1D[20:36]
P_Var_y1=P_Var_y2=Depths[0:31]
P_Var_z1=T_conv_WS_All[0:31,20:36]-273.15
P_Var_z2=T_nonconv_WS_All[0:31,20:36]-273.15
P_xlable='Latitude'  ;P_ylable='Depth [m]'
P_title='Potential Temperature Contours (red=conv, blue=nonconv) - All Southern Ocean - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plot2Dcontour_2var(P_Var_x1, P_Var_y1, P_Var_z1, P_Var_x2, P_Var_y2, P_Var_z2, P_xlable, P_ylable, P_title, 'r', 'b', 10, 'invert_yaxis')
#fig.savefig(dir_figs+str(GCM)+'_temp_composites_WS_All.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

P_Var_x1=P_Var_x2=Lat_regrid_1D[20:36]
P_Var_y1=P_Var_y2=Depths[0:31]
P_Var_z1=T_conv_WS_60W0E[0:31,20:36]-273.15
P_Var_z2=T_nonconv_WS_60W0E[0:31,20:36]-273.15
P_xlable='Latitude'  ;P_ylable='Depth [m]'
P_title='Potential Temperature Contours (red=conv, blue=nonconv) - 60W-0E - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plot2Dcontour_2var(P_Var_x1, P_Var_y1, P_Var_z1, P_Var_x2, P_Var_y2, P_Var_z2, P_xlable, P_ylable, P_title, 'r', 'b', 10, 'invert_yaxis')
#fig.savefig(dir_figs+str(GCM)+'_temp_composites_WS_60W0E.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

P_Var_x1=P_Var_x2=Lat_regrid_1D[20:30]
P_Var_y1=P_Var_y2=Depths[0:31]
P_Var_z1=T_conv_WS_60W0E[0:31,20:30]-273.15
P_Var_z2=T_nonconv_WS_60W0E[0:31,20:30]-273.15
P_xlable='Latitude'  ;P_ylable='Depth [m]'
P_title='Potential Temperature Contours (red=conv, blue=nonconv) - 60W-0E - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
#fig = func_plot2Dcontour_2var(P_Var_x1, P_Var_y1, P_Var_z1, P_Var_x2, P_Var_y2, P_Var_z2, P_xlable, P_ylable, P_title, 'r', 'b', 10, 'invert_yaxis')


R_ii_LAB = np.where(LAB_index_norm_rm <-0.5)
R_jj_LAB = np.where(LAB_index_norm_rm >0.5)

Temp_allyears_lab_c=Temp_allyears[R_ii_LAB]
Temp_allyears_lab_n=Temp_allyears[R_jj_LAB]

Temp_allyears_LAB_c=Temp_allyears_lab_c[:,:,140:156,295:341]
Temp_allyears_LAB_c=np.nanmean(Temp_allyears_LAB_c,axis=2)
Temp_allyears_LAB_c=np.nanmean(Temp_allyears_LAB_c,axis=0)
Temp_allyears_LAB_n=Temp_allyears_lab_n[:,:,140:156,295:341]
Temp_allyears_LAB_n=np.nanmean(Temp_allyears_LAB_n,axis=2)
Temp_allyears_LAB_n=np.nanmean(Temp_allyears_LAB_n,axis=0)


P_Var_x1=P_Var_x2=Lon_regrid_1D[295:341]-360
P_Var_y1=P_Var_y2=Depths[0:40]
P_Var_z1=Temp_allyears_LAB_c[0:40,:]-273.15
P_Var_z2=Temp_allyears_LAB_n[0:40,:]-273.15
P_xlable='Longitude'  ;P_ylable='Depth [m]'
P_title='Potential Temperature Contours - Labrador Sea (50N-65N ave.) (red=conv, blue=nonconv) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plot2Dcontour_2var(P_Var_x1, P_Var_y1, P_Var_z1, P_Var_x2, P_Var_y2, P_Var_z2, P_xlable, P_ylable, P_title, 'r', 'b', 40, 'invert_yaxis')
fig.savefig(dir_figs+'Temp_'+str(GCM)+'_LAB_C1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

P_Var_x1=P_Var_x2=Lon_regrid_1D[295:341]-360
P_Var_y1=P_Var_y2=Depths[0:31]
P_Var_z1=Temp_allyears_LAB_c[0:31,:]-273.15
P_Var_z2=Temp_allyears_LAB_n[0:31,:]-273.15
P_xlable='Longitude'  ;P_ylable='Depth [m]'
P_title='Potential Temperature Contours - Labrador Sea (50N-65N ave.) (red=conv, blue=nonconv) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plot2Dcontour_2var(P_Var_x1, P_Var_y1, P_Var_z1, P_Var_x2, P_Var_y2, P_Var_z2, P_xlable, P_ylable, P_title, 'r', 'b', 40, 'invert_yaxis')
fig.savefig(dir_figs+'Temp_'+str(GCM)+'_LAB_C2.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


#### Salinity ####

Salinity_allyears=[]

for t in range(int((end_date_i+1-start_date_i)/12)):
    print('Salinity calc - Year: ', year_start+t)
    data_so_extracted=dset_so.extract_data(dset_so.variable,start_date_i+12*t+1-1,start_date_i+12*t+12-1)
    data_so_extracted=np.nanmean(data_so_extracted, axis=0)
    data_so_extracted=np.squeeze(data_so_extracted)

    data_i = func_regrid(data_so_extracted, dset_so.y, dset_so.x, Lat_regrid_2D, Lon_regrid_2D)
    data_i[data_i>100000]=np.nan
    
    Salinity_allyears.append(data_i)

Salinity_allyears=np.asarray(Salinity_allyears)
Salinity_allyears_c=Salinity_allyears[R_ii_WS]
Salinity_allyears_n=Salinity_allyears[R_jj_WS]
S_conv_WS_All=np.nanmean(Salinity_allyears_c,axis=0)
S_nonconv_WS_All=np.nanmean(Salinity_allyears_n,axis=0)

S_conv_WS_60W0E=S_conv_WS_All[:,:,300:]
S_nonconv_WS_60W0E=S_nonconv_WS_All[:,:,300:]

S_conv_WS_All=np.nanmean(S_conv_WS_All,axis=2)
S_nonconv_WS_All=np.nanmean(S_nonconv_WS_All,axis=2)
S_conv_WS_60W0E=np.nanmean(S_conv_WS_60W0E,axis=2)
S_nonconv_WS_60W0E=np.nanmean(S_nonconv_WS_60W0E,axis=2)

np.save('Salinity_allyears_GFDL-ESM2G_500yr.npy',Salinity_allyears)
#Salinity_allyears=np.load('Salinity_allyears_GFDL-ESM2G_500yr.npy')

P_Var_x1=P_Var_x2=Lat_regrid_1D[20:36]
P_Var_y1=P_Var_y2=Depths[0:31]
P_Var_z1=S_conv_WS_All[0:31,20:36]
P_Var_z2=S_nonconv_WS_All[0:31,20:36]
P_xlable='Latitude'  ;P_ylable='Depth [m]'
P_title='Salinity Contours (red=conv, blue=nonconv) - All Southern Ocean - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plot2Dcontour_2var(P_Var_x1, P_Var_y1, P_Var_z1, P_Var_x2, P_Var_y2, P_Var_z2, P_xlable, P_ylable, P_title, 'r', 'b', 10, 'invert_yaxis')
#fig.savefig(dir_figs+str(GCM)+'_salt_composites_WS_All.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

P_Var_x1=P_Var_x2=Lat_regrid_1D[20:36]
P_Var_y1=P_Var_y2=Depths[0:31]
P_Var_z1=S_conv_WS_60W0E[0:31,20:36]
P_Var_z2=S_nonconv_WS_60W0E[0:31,20:36]
P_xlable='Latitude'  ;P_ylable='Depth [m]'
P_title='Salinity Contours (red=conv, blue=nonconv) - 60W-0E - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plot2Dcontour_2var(P_Var_x1, P_Var_y1, P_Var_z1, P_Var_x2, P_Var_y2, P_Var_z2, P_xlable, P_ylable, P_title, 'r', 'b', 10, 'invert_yaxis')
#fig.savefig(dir_figs+str(GCM)+'_salt_composites_WS_60W0E.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

P_Var_x1=P_Var_x2=Lat_regrid_1D[20:30]
P_Var_y1=P_Var_y2=Depths[0:31]
P_Var_z1=S_conv_WS_60W0E[0:31,20:30]
P_Var_z2=S_nonconv_WS_60W0E[0:31,20:30]
P_xlable='Latitude'  ;P_ylable='Depth [m]'
P_title='Salinity Contours (red=conv, blue=nonconv) - 60W-0E - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plot2Dcontour_2var(P_Var_x1, P_Var_y1, P_Var_z1, P_Var_x2, P_Var_y2, P_Var_z2, P_xlable, P_ylable, P_title, 'r', 'b', 10, 'invert_yaxis')


R_ii_LAB = np.where(LAB_index_norm_rm <-0.5)
R_jj_LAB = np.where(LAB_index_norm_rm >0.5)

Salinity_allyears_lab_c=Salinity_allyears[R_ii_LAB]
Salinity_allyears_lab_n=Salinity_allyears[R_jj_LAB]

Salinity_allyears_LAB_c=Salinity_allyears_lab_c[:,:,140:156,295:341]
Salinity_allyears_LAB_c=np.nanmean(Salinity_allyears_LAB_c,axis=2)
Salinity_allyears_LAB_c=np.nanmean(Salinity_allyears_LAB_c,axis=0)
Salinity_allyears_LAB_n=Salinity_allyears_lab_n[:,:,140:156,295:341]
Salinity_allyears_LAB_n=np.nanmean(Salinity_allyears_LAB_n,axis=2)
Salinity_allyears_LAB_n=np.nanmean(Salinity_allyears_LAB_n,axis=0)

P_Var_x1=P_Var_x2=Lon_regrid_1D[295:341]-360
P_Var_y1=P_Var_y2=Depths[0:40]
P_Var_z1=Salinity_allyears_LAB_c[0:40,:]
P_Var_z2=Salinity_allyears_LAB_n[0:40,:]
P_xlable='Longitude'  ;P_ylable='Depth [m]'
P_title='Salinity Contours - Labrador Sea (50N-65N ave.) (red=conv, blue=nonconv) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plot2Dcontour_2var(P_Var_x1, P_Var_y1, P_Var_z1, P_Var_x2, P_Var_y2, P_Var_z2, P_xlable, P_ylable, P_title, 'r', 'b', 40, 'invert_yaxis')
#fig.savefig(dir_figs+'Salinity_'+str(GCM)+'_LAB_C1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

P_Var_x1=P_Var_x2=Lon_regrid_1D[295:341]-360
P_Var_y1=P_Var_y2=Depths[0:31]
P_Var_z1=Salinity_allyears_LAB_c[0:31,:]
P_Var_z2=Salinity_allyears_LAB_n[0:31,:]
P_xlable='Longitude'  ;P_ylable='Depth [m]'
P_title='Salinity Contours - Labrador Sea (50N-65N ave.) (red=conv, blue=nonconv) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plot2Dcontour_2var(P_Var_x1, P_Var_y1, P_Var_z1, P_Var_x2, P_Var_y2, P_Var_z2, P_xlable, P_ylable, P_title, 'r', 'b', 40, 'invert_yaxis')
#fig.savefig(dir_figs+'Salinity_'+str(GCM)+'_LAB_C2.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

###############################################################################
###############################################################################
import os
import shelve
dir_pwd = os.getcwd() # Gets the current directory (and in which the code is placed)
filename_out = (dir_pwd + '/Results/AllResults_'+GCM+'_500yr.out') # Directory to save processed data

my_shelf = shelve.open(filename_out)
#globals()[key]=my_shelf[key]
globals()['WS_index_norm_rm']=my_shelf['WS_index_norm_rm']
globals()['LAB_index_norm_rm']=my_shelf['LAB_index_norm_rm']
my_shelf.close()

R_ii_WS = np.where(WS_index_norm_rm <-0.5)
R_jj_WS = np.where(WS_index_norm_rm >0.5)
R_ii_LAB = np.where(LAB_index_norm_rm <-0.5)
R_jj_LAB = np.where(LAB_index_norm_rm >0.5)

###############################################################################
########################       Density Plots       ############################
###############################################################################
Density_allyears=np.load('Density_allyears_GFDL-ESM2G_500yr.npy')
##############   Labrador Sea   ##############
Density_allyears_lab_c=Density_allyears[R_ii_LAB]
Density_allyears_lab_n=Density_allyears[R_jj_LAB]

Density_allyears_LAB=Density_allyears[:,:,140:156,295:341]
Density_allyears_LAB=np.nanmean(Density_allyears_LAB,axis=2)
Density_allyears_LAB=np.nanmean(Density_allyears_LAB,axis=0)

Density_allyears_LAB_c=Density_allyears_lab_c[:,:,140:156,295:341]
Density_allyears_LAB_c=np.nanmean(Density_allyears_LAB_c,axis=2)
Density_allyears_LAB_c=np.nanmean(Density_allyears_LAB_c,axis=0)
Density_allyears_LAB_n=Density_allyears_lab_n[:,:,140:156,295:341]
Density_allyears_LAB_n=np.nanmean(Density_allyears_LAB_n,axis=2)
Density_allyears_LAB_n=np.nanmean(Density_allyears_LAB_n,axis=0)

Density_allyears_LAB_ave=np.nanmean(Density_allyears_LAB,axis=1)
Density_allyears_LAB_c_ave=np.nanmean(Density_allyears_LAB_c,axis=1)
Density_allyears_LAB_n_ave=np.nanmean(Density_allyears_LAB_n,axis=1)

dRho_dZ_LAB=np.zeros(( Depths.shape[0]-1))
dRho_dZ_LAB_c=np.zeros(( Depths.shape[0]-1))
dRho_dZ_LAB_n=np.zeros(( Depths.shape[0]-1))
for ii in range(0,dRho_dZ_LAB.shape[0]):
    dRho_dZ_LAB[ii]=(Density_allyears_LAB_ave[ii]-Density_allyears_LAB_ave[ii+1])/(Depths[ii]-Depths[ii+1])
    dRho_dZ_LAB_c[ii]=(Density_allyears_LAB_c_ave[ii]-Density_allyears_LAB_c_ave[ii+1])/(Depths[ii]-Depths[ii+1])
    dRho_dZ_LAB_n[ii]=(Density_allyears_LAB_n_ave[ii]-Density_allyears_LAB_n_ave[ii+1])/(Depths[ii]-Depths[ii+1])
    

P_Var_x=Density_allyears_LAB_ave[0:40]-1000  ;P_Var_y=Depths[0:40]
P_xlable='Potential Density [kg/m3]'  ;P_ylable='Depth [m]'
P_title='Potential Density Ave - Labrador Sea [20W-65W,50N-65N] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plotline_1var(P_Var_x, P_Var_y, P_xlable, P_ylable, P_title, 'darkcyan', 'invert_yaxis')
#fig.savefig(dir_figs+str(GCM)+'_density_depth_LAB_20W65W_50N65N.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

P_Var_x1=Density_allyears_LAB_c_ave[0:40]-1000
P_Var_x2=Density_allyears_LAB_n_ave[0:40]-1000
P_Var_y1=P_Var_y2=Depths[0:40]
P_xlable='Potential Density [kg/m3]'  ;P_ylable='Depth [m]'
P_title='Potential Density Ave - Lab. Sea [20W-65W,50N-65N] (red=Lab. Sea conv, blue=nonconv) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plotline_2var(P_Var_x1, P_Var_y1, P_Var_x2, P_Var_y2, P_xlable, P_ylable, P_title, 'r', 'b', 'Lab. Sea conv', 'Lab. Sea nonconv', 'best', 'invert_yaxis')
#fig.savefig(dir_figs+str(GCM)+'_density_depth_LAB_20W65W_50N65N_conv.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

P_Var_x=dRho_dZ_LAB[0:40]  ;P_Var_y=Depths[0:40]
P_xlable='dRho/dZ [kg/m3/m]'  ;P_ylable='Depth [m]'
P_title='dRho/dZ - Labrador Sea [20W-65W,50N-65N] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plotline_1var(P_Var_x, P_Var_y, P_xlable, P_ylable, P_title, 'darkcyan', 'invert_yaxis')
#fig.savefig(dir_figs+str(GCM)+'_dRho_dZ_LAB_20W65W_50N65N.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

P_Var_x1=dRho_dZ_LAB_c[0:40]
P_Var_x2=dRho_dZ_LAB_n[0:40]
P_Var_y1=P_Var_y2=Depths[0:40]
P_xlable='dRho/dZ [kg/m3/m]'  ;P_ylable='Depth [m]'
P_title='dRho/dZ - Lab. Sea [20W-65W,50N-65N] (red=Lab. Sea conv, blue=nonconv) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plotline_2var(P_Var_x1, P_Var_y1, P_Var_x2, P_Var_y2, P_xlable, P_ylable, P_title, 'r', 'b', 'Lab. Sea conv', 'Lab. Sea nonconv', 'best', 'invert_yaxis')
#fig.savefig(dir_figs+str(GCM)+'_dRho_dZ_LAB_20W65W_50N65N_conv.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

#######################################

dRho_dZ_LAB_Box=np.zeros(( Density_allyears_LAB.shape[0]-1, Density_allyears_LAB.shape[1]))
dRho_dZ_LAB_c_Box=np.zeros(( Density_allyears_LAB.shape[0]-1, Density_allyears_LAB.shape[1]))
dRho_dZ_LAB_n_Box=np.zeros(( Density_allyears_LAB.shape[0]-1, Density_allyears_LAB.shape[1]))
for ii in range(0,Density_allyears_LAB.shape[0]-1):
    for jj in range(0,Density_allyears_LAB.shape[1]):
        dRho_dZ_LAB_Box[ii,jj]=(Density_allyears_LAB[ii,jj] - Density_allyears_LAB[ii+1,jj])/(Depths[ii]-Depths[ii+1])
        dRho_dZ_LAB_c_Box[ii,jj]=(Density_allyears_LAB_c[ii,jj] - Density_allyears_LAB_c[ii+1,jj])/(Depths[ii]-Depths[ii+1])
        dRho_dZ_LAB_n_Box[ii,jj]=(Density_allyears_LAB_n[ii,jj] - Density_allyears_LAB_n[ii+1,jj])/(Depths[ii]-Depths[ii+1])


P_Var_x1=Lon_regrid_1D[295:341]-360  ; P_Var_y1=Depths[0:31]
P_Var_z1=dRho_dZ_LAB_Box[0:31,:]
P_range=np.linspace(0,0.03) 
P_xlable='Longitude'  ;P_ylable='Depth [m]';  P_unit = '(C)'
P_title='dRho/dZ - All Years - Labrador Sea [20W-65W,50N-65N] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plot2Dcolor_contourf(P_Var_x1, P_Var_y1, P_Var_z1, P_range, P_xlable, P_ylable, P_title, P_unit, plt.cm.jet, 'invert_yaxis')
fig.savefig(dir_figs+str(GCM)+'_dRho_dZ_LAB_20W65W_50N65N_Box.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

P_Var_x1=Lon_regrid_1D[295:341]-360  ; P_Var_y1=Depths[0:31]
P_Var_z1=dRho_dZ_LAB_c_Box[0:31,:]
P_range=np.linspace(0,0.03) 
P_xlable='Longitude'  ;P_ylable='Depth [m]';  P_unit = '(C)'
P_title='dRho/dZ - Convective Years - Labrador Sea [20W-65W,50N-65N] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plot2Dcolor_contourf(P_Var_x1, P_Var_y1, P_Var_z1, P_range, P_xlable, P_ylable, P_title, P_unit, plt.cm.jet, 'invert_yaxis')
fig.savefig(dir_figs+str(GCM)+'_dRho_dZ_LAB_20W65W_50N65N_conv_Box.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

P_Var_x1=Lon_regrid_1D[295:341]-360  ; P_Var_y1=Depths[0:31]
P_Var_z1=dRho_dZ_LAB_n_Box[0:31,:]
P_range=np.linspace(0,0.03) 
P_xlable='Longitude'  ;P_ylable='Depth [m]';  P_unit = '(C)'
P_title='dRho/dZ - Non-Convective Years - Labrador Sea [20W-65W,50N-65N] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plot2Dcolor_contourf(P_Var_x1, P_Var_y1, P_Var_z1, P_range, P_xlable, P_ylable, P_title, P_unit, plt.cm.jet, 'invert_yaxis')
fig.savefig(dir_figs+str(GCM)+'_dRho_dZ_LAB_20W65W_50N65N_nonconv_Box.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

P_Var_x1=Lon_regrid_1D[295:341]-360  ; P_Var_y1=Depths[0:31]
P_Var_z1=dRho_dZ_LAB_c_Box[0:31,:] - dRho_dZ_LAB_n_Box[0:31,:]
P_range=np.linspace(-0.0024,0.0024,41)
P_xlable='Longitude'  ;P_ylable='Depth [m]';  P_unit = '(C)'
P_title='dRho/dZ - Convective minus Non-Convective - Labrador Sea [20W-65W,50N-65N] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plot2Dcolor_contourf(P_Var_x1, P_Var_y1, P_Var_z1, P_range, P_xlable, P_ylable, P_title, P_unit, plt.cm.seismic, 'invert_yaxis')
fig.savefig(dir_figs+str(GCM)+'_dRho_dZ_LAB_20W65W_50N65N_composite_Box.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


Density_allyears_LAB_B1=Density_allyears[:,:,140:156,295:320]
Density_allyears_LAB_B1=np.nanmean(Density_allyears_LAB_B1,axis=2); Density_allyears_LAB_B1=np.nanmean(Density_allyears_LAB_B1,axis=0)

Density_allyears_LAB_c_B1=Density_allyears_lab_c[:,:,140:156,295:320]
Density_allyears_LAB_c_B1=np.nanmean(Density_allyears_LAB_c_B1,axis=2); Density_allyears_LAB_c_B1=np.nanmean(Density_allyears_LAB_c_B1,axis=0)

Density_allyears_LAB_n_B1=Density_allyears_lab_n[:,:,140:156,295:320]
Density_allyears_LAB_n_B1=np.nanmean(Density_allyears_LAB_n_B1,axis=2); Density_allyears_LAB_n_B1=np.nanmean(Density_allyears_LAB_n_B1,axis=0)

Density_allyears_LAB_c_ave_B1=np.nanmean(Density_allyears_LAB_c_B1,axis=1); Density_allyears_LAB_n_ave_B1=np.nanmean(Density_allyears_LAB_n_B1,axis=1)

P_Var_x1=Density_allyears_LAB_c_ave_B1[0:40]-1000
P_Var_x2=Density_allyears_LAB_n_ave_B1[0:40]-1000
P_Var_y1=P_Var_y2=Depths[0:40]
P_xlable='Potential Density [kg/m3]'  ;P_ylable='Depth [m]'
P_title='Potential Density Ave - Labrador Box1 [65W-41W, 50N-65N] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plotline_2var(P_Var_x1, P_Var_y1, P_Var_x2, P_Var_y2, P_xlable, P_ylable, P_title, 'r', 'b', 'Lab. Sea conv', 'Lab. Sea nonconv', 'lower left', 'invert_yaxis')
plt.xlim(26.2,28)
#fig.savefig(dir_figs+str(GCM)+'_density_depth_LAB_41W65W_50N65N_conv.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


Density_allyears_LAB_B2=Density_allyears[:,:,140:156,320:341]
Density_allyears_LAB_B2=np.nanmean(Density_allyears_LAB_B2,axis=2); Density_allyears_LAB_B2=np.nanmean(Density_allyears_LAB_B2,axis=0)

Density_allyears_LAB_c_B2=Density_allyears_lab_c[:,:,140:156,320:341]
Density_allyears_LAB_c_B2=np.nanmean(Density_allyears_LAB_c_B2,axis=2); Density_allyears_LAB_c_B2=np.nanmean(Density_allyears_LAB_c_B2,axis=0)

Density_allyears_LAB_n_B2=Density_allyears_lab_n[:,:,140:156,320:341]
Density_allyears_LAB_n_B2=np.nanmean(Density_allyears_LAB_n_B2,axis=2); Density_allyears_LAB_n_B2=np.nanmean(Density_allyears_LAB_n_B2,axis=0)

Density_allyears_LAB_c_ave_B2=np.nanmean(Density_allyears_LAB_c_B2,axis=1); Density_allyears_LAB_n_ave_B2=np.nanmean(Density_allyears_LAB_n_B2,axis=1)

P_Var_x1=Density_allyears_LAB_c_ave_B2[0:40]-1000
P_Var_x2=Density_allyears_LAB_n_ave_B2[0:40]-1000
P_Var_y1=P_Var_y2=Depths[0:40]
P_xlable='Potential Density [kg/m3]'  ;P_ylable='Depth [m]'
P_title='Potential Density Ave - Labrador Box2 [41W-20W, 50N-65N] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plotline_2var(P_Var_x1, P_Var_y1, P_Var_x2, P_Var_y2, P_xlable, P_ylable, P_title, 'r', 'b', 'Lab. Sea conv', 'Lab. Sea nonconv', 'lower left', 'invert_yaxis')
plt.xlim(26.2,28)
#fig.savefig(dir_figs+str(GCM)+'_density_depth_LAB_20W41W_50N65N_conv.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

##############   Weddell Sea   ##############
Density_allyears_ws_c=Density_allyears[R_ii_WS]
Density_allyears_ws_n=Density_allyears[R_jj_WS]

Density_allyears_WS=Density_allyears[:,:,20:35,300:]
Density_allyears_WS=np.nanmean(Density_allyears_WS,axis=2)
Density_allyears_WS=np.nanmean(Density_allyears_WS,axis=0)

Density_allyears_WS_c=Density_allyears_ws_c[:,:,20:35,300:]
Density_allyears_WS_c=np.nanmean(Density_allyears_WS_c,axis=2)
Density_allyears_WS_c=np.nanmean(Density_allyears_WS_c,axis=0)
Density_allyears_WS_n=Density_allyears_ws_n[:,:,20:35,300:]
Density_allyears_WS_n=np.nanmean(Density_allyears_WS_n,axis=2)
Density_allyears_WS_n=np.nanmean(Density_allyears_WS_n,axis=0)

Density_allyears_WS_ave=np.nanmean(Density_allyears_WS,axis=1)
Density_allyears_WS_c_ave=np.nanmean(Density_allyears_WS_c,axis=1)
Density_allyears_WS_n_ave=np.nanmean(Density_allyears_WS_n,axis=1)

dRho_dZ_WS=np.zeros(( Depths.shape[0]-1))
dRho_dZ_WS_c=np.zeros(( Depths.shape[0]-1))
dRho_dZ_WS_n=np.zeros(( Depths.shape[0]-1))
for ii in range(0,dRho_dZ_LAB.shape[0]):
    dRho_dZ_WS[ii]=(Density_allyears_WS_ave[ii]-Density_allyears_WS_ave[ii+1])/(Depths[ii]-Depths[ii+1])
    dRho_dZ_WS_c[ii]=(Density_allyears_WS_c_ave[ii]-Density_allyears_WS_c_ave[ii+1])/(Depths[ii]-Depths[ii+1])
    dRho_dZ_WS_n[ii]=(Density_allyears_WS_n_ave[ii]-Density_allyears_WS_n_ave[ii+1])/(Depths[ii]-Depths[ii+1])

P_Var_x=Density_allyears_WS_ave[0:40]-1000  ;P_Var_y=Depths[0:40]
P_xlable='Potential Density [kg/m3]'  ;P_ylable='Depth [m]'
P_title='Potential Density Ave- Weddell Sea [60W-0W,70S-55S] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plotline_1var(P_Var_x, P_Var_y, P_xlable, P_ylable, P_title, 'darkcyan', 'invert_yaxis')
#fig.savefig(dir_figs+str(GCM)+'_density_depth_WS_60W0W_70S55S.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

P_Var_x1=Density_allyears_WS_c_ave[0:40]-1000
P_Var_x2=Density_allyears_WS_n_ave[0:40]-1000
P_Var_y1=P_Var_y2=Depths[0:40]
P_xlable='Potential Density [kg/m3]'  ;P_ylable='Depth [m]'
P_title='Potential Density Ave - Wedd. Sea [60W-0W,70S-55S] (red=Wedd. Sea conv, blue=nonconv) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plotline_2var(P_Var_x1, P_Var_y1, P_Var_x2, P_Var_y2, P_xlable, P_ylable, P_title, 'r', 'b', 'Wedd. Sea conv', 'Wedd. Sea nonconv', 'lower left', 'invert_yaxis')
#fig.savefig(dir_figs+str(GCM)+'_density_depth_WS_60W0W_70S55S_conv.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

P_Var_x=dRho_dZ_WS[0:40]  ;P_Var_y=Depths[0:40]
P_xlable='dRho/dZ [kg/m3/m]'  ;P_ylable='Depth [m]'
P_title='dRho/dZ - Weddell Sea [60W-0W,70S-55S] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plotline_1var(P_Var_x, P_Var_y, P_xlable, P_ylable, P_title, 'darkcyan', 'invert_yaxis')
#fig.savefig(dir_figs+str(GCM)+'_dRho_dZ_WS_60W0W_70S55S.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

P_Var_x1=dRho_dZ_WS_c[0:40]
P_Var_x2=dRho_dZ_WS_n[0:40]
P_Var_y1=P_Var_y2=Depths[0:40]
P_xlable='dRho/dZ [kg/m3/m]'  ;P_ylable='Depth [m]'
P_title='dRho/dZ - Weddell Sea [60W-0W,70S-55S] (red=Wedd. Sea conv, blue=nonconv) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plotline_2var(P_Var_x1, P_Var_y1, P_Var_x2, P_Var_y2, P_xlable, P_ylable, P_title, 'r', 'b', 'Wedd. Sea conv', 'Wedd. Sea nonconv', 'lower right', 'invert_yaxis')
#fig.savefig(dir_figs+str(GCM)+'_density_depth_WS_60W0W_70S55S_conv.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

###############################################################################
######################       Temperature Plots       ##########################
###############################################################################
Temp_allyears=np.load('Temp_allyears_GFDL-ESM2G_500yr.npy')
##############   Labrador Sea   ##############
Temp_allyears_lab_c=Temp_allyears[R_ii_LAB]
Temp_allyears_lab_n=Temp_allyears[R_jj_LAB]

Temp_allyears_LAB=Temp_allyears[:,:,140:156,295:341]
Temp_allyears_LAB=np.nanmean(Temp_allyears_LAB,axis=2)
Temp_allyears_LAB=np.nanmean(Temp_allyears_LAB,axis=0)

Temp_allyears_LAB_c=Temp_allyears_lab_c[:,:,140:156,295:341]
Temp_allyears_LAB_c=np.nanmean(Temp_allyears_LAB_c,axis=2)
Temp_allyears_LAB_c=np.nanmean(Temp_allyears_LAB_c,axis=0)
Temp_allyears_LAB_n=Temp_allyears_lab_n[:,:,140:156,295:341]
Temp_allyears_LAB_n=np.nanmean(Temp_allyears_LAB_n,axis=2)
Temp_allyears_LAB_n=np.nanmean(Temp_allyears_LAB_n,axis=0)

Temp_allyears_LAB_ave=np.nanmean(Temp_allyears_LAB,axis=1)
Temp_allyears_LAB_c_ave=np.nanmean(Temp_allyears_LAB_c,axis=1)
Temp_allyears_LAB_n_ave=np.nanmean(Temp_allyears_LAB_n,axis=1)

dTheta_dZ_LAB=np.zeros(( Depths.shape[0]-1))
dTheta_dZ_LAB_c=np.zeros(( Depths.shape[0]-1))
dTheta_dZ_LAB_n=np.zeros(( Depths.shape[0]-1))
for ii in range(0,dTheta_dZ_LAB.shape[0]):
    dTheta_dZ_LAB[ii]=(Temp_allyears_LAB_ave[ii]-Temp_allyears_LAB_ave[ii+1])/(Depths[ii]-Depths[ii+1])
    dTheta_dZ_LAB_c[ii]=(Temp_allyears_LAB_c_ave[ii]-Temp_allyears_LAB_c_ave[ii+1])/(Depths[ii]-Depths[ii+1])
    dTheta_dZ_LAB_n[ii]=(Temp_allyears_LAB_n_ave[ii]-Temp_allyears_LAB_n_ave[ii+1])/(Depths[ii]-Depths[ii+1])

P_Var_x=Temp_allyears_LAB_ave[0:40]-273.15  ;P_Var_y=Depths[0:40]
P_xlable='Potential Temperature [C]'  ;P_ylable='Depth [m]'
P_title='Potential Temperature Ave - Labrador Sea [20W-65W,50N-65N] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plotline_1var(P_Var_x, P_Var_y, P_xlable, P_ylable, P_title, 'darkcyan', 'invert_yaxis')
#fig.savefig(dir_figs+str(GCM)+'_temp_depth_LAB_20W65W_50N65N.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

P_Var_x1=Temp_allyears_LAB_c_ave[0:40]-273.15
P_Var_x2=Temp_allyears_LAB_n_ave[0:40]-273.15
P_Var_y1=P_Var_y2=Depths[0:40]
P_xlable='Potential Temperature [C]'  ;P_ylable='Depth [m]'
P_title='Potential Temp. Ave - Lab. Sea [20W-65W,50N-65N] (red=Lab. Sea conv, blue=nonconv) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plotline_2var(P_Var_x1, P_Var_y1, P_Var_x2, P_Var_y2, P_xlable, P_ylable, P_title, 'r', 'b', 'Wedd. Sea conv', 'Wedd. Sea nonconv', 'best', 'invert_yaxis')
#fig.savefig(dir_figs+str(GCM)+'_temp_depth_LAB_20W65W_50N65N_conv.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

P_Var_x=dTheta_dZ_LAB[0:40]  ;P_Var_y=Depths[0:40]
P_xlable='dTheta/dZ [C/m]'  ;P_ylable='Depth [m]'
P_title='dTheta/dZ - Labrador Sea [20W-65W,50N-65N] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plotline_1var(P_Var_x, P_Var_y, P_xlable, P_ylable, P_title, 'darkcyan', 'invert_yaxis')
#fig.savefig(dir_figs+str(GCM)+'_dTheta_dZ_LAB_20W65W_50N65N.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

P_Var_x1=dTheta_dZ_LAB_c[0:40]
P_Var_x2=dTheta_dZ_LAB_n[0:40]
P_Var_y1=P_Var_y2=Depths[0:40]
P_xlable='dTheta/dZ [C/m]'  ;P_ylable='Depth [m]'
P_title='dTheta/dZ - Labrador Sea [20W-65W,50N-65N] (red=Lab. Sea conv, blue=nonconv) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plotline_2var(P_Var_x1, P_Var_y1, P_Var_x2, P_Var_y2, P_xlable, P_ylable, P_title, 'r', 'b', 'Wedd. Sea conv', 'Wedd. Sea nonconv', 'best', 'invert_yaxis')
#fig.savefig(dir_figs+str(GCM)+'_dTheta_dZ_LAB_20W65W_50N65N_conv.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

#######################################

Temp_allyears_LAB_B1=Temp_allyears[:,:,140:156,295:320]
Temp_allyears_LAB_B1=np.nanmean(Temp_allyears_LAB_B1,axis=2); Temp_allyears_LAB_B1=np.nanmean(Temp_allyears_LAB_B1,axis=0)

Temp_allyears_LAB_c_B1=Temp_allyears_lab_c[:,:,140:156,295:320]
Temp_allyears_LAB_c_B1=np.nanmean(Temp_allyears_LAB_c_B1,axis=2); Temp_allyears_LAB_c_B1=np.nanmean(Temp_allyears_LAB_c_B1,axis=0)

Temp_allyears_LAB_n_B1=Temp_allyears_lab_n[:,:,140:156,295:320]
Temp_allyears_LAB_n_B1=np.nanmean(Temp_allyears_LAB_n_B1,axis=2); Temp_allyears_LAB_n_B1=np.nanmean(Temp_allyears_LAB_n_B1,axis=0)

Temp_allyears_LAB_c_ave_B1=np.nanmean(Temp_allyears_LAB_c_B1,axis=1); Temp_allyears_LAB_n_ave_B1=np.nanmean(Temp_allyears_LAB_n_B1,axis=1)


P_Var_x1=Temp_allyears_LAB_c_ave_B1[0:40]-273.15
P_Var_x2=Temp_allyears_LAB_n_ave_B1[0:40]-273.15
P_Var_y1=P_Var_y2=Depths[0:40]
P_xlable='Potential Temperature [C]'  ;P_ylable='Depth [m]'
P_title='Potential Temperature Ave - Labrador Box1 [20W-41W, 50N-65N] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plotline_2var(P_Var_x1, P_Var_y1, P_Var_x2, P_Var_y2, P_xlable, P_ylable, P_title, 'r', 'b', 'Lab. Sea conv', 'Lab. Sea nonconv', 'lower right', 'invert_yaxis')
plt.xlim(1,9.8)
fig.savefig(dir_figs+str(GCM)+'_temp_depth_LAB_20W41W_50N65N_conv.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


Temp_allyears_LAB_B2=Temp_allyears[:,:,140:156,320:341]
Temp_allyears_LAB_B2=np.nanmean(Temp_allyears_LAB_B2,axis=2); Temp_allyears_LAB_B2=np.nanmean(Temp_allyears_LAB_B2,axis=0)

Temp_allyears_LAB_c_B2=Temp_allyears_lab_c[:,:,140:156,320:341]
Temp_allyears_LAB_c_B2=np.nanmean(Temp_allyears_LAB_c_B2,axis=2); Temp_allyears_LAB_c_B2=np.nanmean(Temp_allyears_LAB_c_B2,axis=0)

Temp_allyears_LAB_n_B2=Temp_allyears_lab_n[:,:,140:156,320:341]
Temp_allyears_LAB_n_B2=np.nanmean(Temp_allyears_LAB_n_B2,axis=2); Temp_allyears_LAB_n_B2=np.nanmean(Temp_allyears_LAB_n_B2,axis=0)

Temp_allyears_LAB_c_ave_B2=np.nanmean(Temp_allyears_LAB_c_B2,axis=1); Temp_allyears_LAB_n_ave_B2=np.nanmean(Temp_allyears_LAB_n_B2,axis=1)


P_Var_x1=Temp_allyears_LAB_c_ave_B2[0:40]-273.15
P_Var_x2=Temp_allyears_LAB_n_ave_B2[0:40]-273.15
P_Var_y1=P_Var_y2=Depths[0:40]
P_xlable='Potential Temperature [C]'  ;P_ylable='Depth [m]'
P_title='Potential Temperature Ave - Labrador Box2 [41W-65W, 50N-65N] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plotline_2var(P_Var_x1, P_Var_y1, P_Var_x2, P_Var_y2, P_xlable, P_ylable, P_title, 'r', 'b', 'Lab. Sea conv', 'Lab. Sea nonconv', 'lower right', 'invert_yaxis')
plt.xlim(1,9.8)
fig.savefig(dir_figs+str(GCM)+'_temp_depth_LAB_41W65W_50N65N_conv.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


##############   Weddell Sea   ##############
Temp_allyears_ws_c=Temp_allyears[R_ii_WS]
Temp_allyears_ws_n=Temp_allyears[R_jj_WS]

Temp_allyears_WS=Temp_allyears[:,:,20:35,300:]
Temp_allyears_WS=np.nanmean(Temp_allyears_WS,axis=2)
Temp_allyears_WS=np.nanmean(Temp_allyears_WS,axis=0)

Temp_allyears_WS_c=Temp_allyears_ws_c[:,:,20:35,300:]
Temp_allyears_WS_c=np.nanmean(Temp_allyears_WS_c,axis=2)
Temp_allyears_WS_c=np.nanmean(Temp_allyears_WS_c,axis=0)
Temp_allyears_WS_n=Temp_allyears_ws_n[:,:,20:35,300:]
Temp_allyears_WS_n=np.nanmean(Temp_allyears_WS_n,axis=2)
Temp_allyears_WS_n=np.nanmean(Temp_allyears_WS_n,axis=0)

Temp_allyears_WS_ave=np.nanmean(Temp_allyears_WS,axis=1)
Temp_allyears_WS_c_ave=np.nanmean(Temp_allyears_WS_c,axis=1)
Temp_allyears_WS_n_ave=np.nanmean(Temp_allyears_WS_n,axis=1)

dTheta_dZ_WS=np.zeros(( Depths.shape[0]-1))
dTheta_dZ_WS_c=np.zeros(( Depths.shape[0]-1))
dTheta_dZ_WS_n=np.zeros(( Depths.shape[0]-1))
for ii in range(0,dRho_dZ_LAB.shape[0]):
    dTheta_dZ_WS[ii]=(Temp_allyears_WS_ave[ii]-Temp_allyears_WS_ave[ii+1])/(Depths[ii]-Depths[ii+1])
    dTheta_dZ_WS_c[ii]=(Temp_allyears_WS_c_ave[ii]-Temp_allyears_WS_c_ave[ii+1])/(Depths[ii]-Depths[ii+1])
    dTheta_dZ_WS_n[ii]=(Temp_allyears_WS_n_ave[ii]-Temp_allyears_WS_n_ave[ii+1])/(Depths[ii]-Depths[ii+1])

fig=plt.figure()
im1=plt.plot(Temp_allyears_WS_ave[0:40]-273.15, Depths[0:40])
plt.gca().invert_yaxis()
plt.xticks(fontsize = 20); plt.yticks(fontsize = 20)
plt.xlabel('Potential Temperature [C]', fontsize=18)
plt.ylabel('Depth [m]', fontsize=18)
plt.title('Potential Temperature Ave - Weddell Sea [60W-0W,70S-55S] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_temp_depth_WS_60W0W_70S55S.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
im1=plt.plot(Temp_allyears_WS_c_ave[0:40]-273.15, Depths[0:40], 'r')
im1=plt.plot(Temp_allyears_WS_n_ave[0:40]-273.15, Depths[0:40], 'b')
plt.gca().invert_yaxis()
plt.xticks(fontsize = 20); plt.yticks(fontsize = 20)
plt.xlabel('Potential Temperature [C]', fontsize=18)
plt.ylabel('Depth', fontsize=18)
plt.title('Potential Temp. Ave - Wedd. Sea [60W-0W,70S-55S] (red=Wedd. Sea conv, blue=nonconv) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_temp_depth_WS_60W0W_70S55S_conv.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


fig=plt.figure()
im1=plt.plot(dTheta_dZ_WS[0:40], Depths[0:40], 'darkgreen')
plt.gca().invert_yaxis()
plt.xticks(fontsize = 20); plt.yticks(fontsize = 20)
plt.xlabel('dTheta/dZ [C/m]', fontsize=18)
plt.ylabel('Depth [m]', fontsize=18)
plt.title('dTheta/dZ - Weddell Sea [60W-0W,70S-55S] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_dTheta_dZ_WS_60W0W_70S55S.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
im1=plt.plot(dTheta_dZ_WS_c[0:40], Depths[0:40], 'r')
im1=plt.plot(dTheta_dZ_WS_n[0:40], Depths[0:40], 'b')
plt.gca().invert_yaxis()
plt.xticks(fontsize = 20); plt.yticks(fontsize = 20)
plt.xlabel('dTheta [C/m]', fontsize=18)
plt.ylabel('Depth [m]', fontsize=18)
plt.title('dTheta/dZ - Weddell Sea [60W-0W,70S-55S] (red=Wedd. Sea conv, blue=nonconv) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_dTheta_dZ_WS_60W0W_70S55S_conv.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


###############################################################################
#######################       Salinity Plots       ###########################
###############################################################################
##############   Labrador Sea   ##############
Salinity_allyears=np.load('Salinity_allyears_GFDL-ESM2G_500yr.npy')
##############   Labrador Sea   ##############
Salinity_allyears_lab_c=Salinity_allyears[R_ii_WS]
Salinity_allyears_lab_n=Salinity_allyears[R_jj_WS]

Salinity_allyears_LAB=Salinity_allyears[:,:,140:156,295:341]
Salinity_allyears_LAB=np.nanmean(Salinity_allyears_LAB,axis=2)
Salinity_allyears_LAB=np.nanmean(Salinity_allyears_LAB,axis=0)

Salinity_allyears_LAB_c=Salinity_allyears_lab_c[:,:,140:156,295:341]
Salinity_allyears_LAB_c=np.nanmean(Salinity_allyears_LAB_c,axis=2)
Salinity_allyears_LAB_c=np.nanmean(Salinity_allyears_LAB_c,axis=0)
Salinity_allyears_LAB_n=Salinity_allyears_lab_n[:,:,140:156,295:341]
Salinity_allyears_LAB_n=np.nanmean(Salinity_allyears_LAB_n,axis=2)
Salinity_allyears_LAB_n=np.nanmean(Salinity_allyears_LAB_n,axis=0)

Salinity_allyears_LAB_ave=np.nanmean(Salinity_allyears_LAB,axis=1)
Salinity_allyears_LAB_c_ave=np.nanmean(Salinity_allyears_LAB_c,axis=1)
Salinity_allyears_LAB_n_ave=np.nanmean(Salinity_allyears_LAB_n,axis=1)

dSal_dZ_LAB=np.zeros(( Depths.shape[0]-1))
dSal_dZ_LAB_c=np.zeros(( Depths.shape[0]-1))
dSal_dZ_LAB_n=np.zeros(( Depths.shape[0]-1))
for ii in range(0,dSal_dZ_LAB.shape[0]):
    dSal_dZ_LAB[ii]=(Salinity_allyears_LAB_ave[ii]-Salinity_allyears_LAB_ave[ii+1])/(Depths[ii]-Depths[ii+1])
    dSal_dZ_LAB_c[ii]=(Salinity_allyears_LAB_c_ave[ii]-Salinity_allyears_LAB_c_ave[ii+1])/(Depths[ii]-Depths[ii+1])
    dSal_dZ_LAB_n[ii]=(Salinity_allyears_LAB_n_ave[ii]-Salinity_allyears_LAB_n_ave[ii+1])/(Depths[ii]-Depths[ii+1])


P_Var_x=Salinity_allyears_LAB_ave[0:40]  ;P_Var_y=Depths[0:40]
P_xlable='Salinity [PSU]'  ;P_ylable='Depth [m]'
P_title='Salinity Ave - Labrador Sea [20W-65W,50N-65N] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plotline_1var(P_Var_x, P_Var_y, P_xlable, P_ylable, P_title, 'darkcyan', 'invert_yaxis')
fig.savefig(dir_figs+str(GCM)+'_salinity_depth_LAB_20W65W_50N65N.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


P_Var_x1=Salinity_allyears_LAB_c_ave[0:40]
P_Var_x2=Salinity_allyears_LAB_n_ave[0:40]
P_Var_y1=P_Var_y2=Depths[0:40]
P_xlable='Salinity [PSU]'  ;P_ylable='Depth [m]'
P_title='Salinity Ave - Labrador Sea [20W-65W,50N-65N] (red=Lab. Sea conv, blue=nonconv) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plotline_2var(P_Var_x1, P_Var_y1, P_Var_x2, P_Var_y2, P_xlable, P_ylable, P_title, 'r', 'b', 'Wedd. Sea conv', 'Wedd. Sea nonconv', 'best', 'invert_yaxis')
fig.savefig(dir_figs+str(GCM)+'_salinity_depth_LAB_20W65W_50N65N_conv.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


P_Var_x=Salinity_allyears_LAB_ave[0:40]  ;P_Var_y=Depths[0:40]
P_xlable='dSal/dZ [PSU/m]'  ;P_ylable='Depth [m]'
P_title='dSal/dZ - Labrador Sea [20W-65W,50N-65N] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plotline_1var(P_Var_x, P_Var_y, P_xlable, P_ylable, P_title, 'darkcyan', 'invert_yaxis')
fig.savefig(dir_figs+str(GCM)+'_dSal_dZ_LAB_20W65W_50N65N.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


P_Var_x1=dSal_dZ_LAB_c[0:40]
P_Var_x2=dSal_dZ_LAB_n[0:40]
P_Var_y1=P_Var_y2=Depths[0:40]
P_xlable='Salinity [PSU]'  ;P_ylable='Depth [m]'
P_title='dSal/dZ - Labrador Sea [20W-65W,50N-65N] (red=Lab. Sea conv, blue=nonconv) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plotline_2var(P_Var_x1, P_Var_y1, P_Var_x2, P_Var_y2, P_xlable, P_ylable, P_title, 'r', 'b', 'Wedd. Sea conv', 'Wedd. Sea nonconv', 'best', 'invert_yaxis')
fig.savefig(dir_figs+str(GCM)+'_dSal_dZ_LAB_20W65W_50N65N_conv.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


#######################################

Salinity_allyears_LAB_B1=Salinity_allyears[:,:,140:156,295:320]
Salinity_allyears_LAB_B1=np.nanmean(Salinity_allyears_LAB_B1,axis=2); Salinity_allyears_LAB_B1=np.nanmean(Salinity_allyears_LAB_B1,axis=0)

Salinity_allyears_LAB_c_B1=Salinity_allyears_lab_c[:,:,140:156,295:320]
Salinity_allyears_LAB_c_B1=np.nanmean(Salinity_allyears_LAB_c_B1,axis=2); Salinity_allyears_LAB_c_B1=np.nanmean(Salinity_allyears_LAB_c_B1,axis=0)

Salinity_allyears_LAB_n_B1=Salinity_allyears_lab_n[:,:,140:156,295:320]
Salinity_allyears_LAB_n_B1=np.nanmean(Salinity_allyears_LAB_n_B1,axis=2); Salinity_allyears_LAB_n_B1=np.nanmean(Salinity_allyears_LAB_n_B1,axis=0)

Salinity_allyears_LAB_c_ave_B1=np.nanmean(Salinity_allyears_LAB_c_B1,axis=1); Salinity_allyears_LAB_n_ave_B1=np.nanmean(Salinity_allyears_LAB_n_B1,axis=1)


P_Var_x1=Salinity_allyears_LAB_c_ave_B1[0:40]
P_Var_x2=Salinity_allyears_LAB_n_ave_B1[0:40]
P_Var_y1=P_Var_y2=Depths[0:40]
P_xlable='Salinity [PSU]'  ;P_ylable='Depth [m]'
P_title='Salinity Ave - Labrador Box1 [65W-41W, 50N-65N] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plotline_2var(P_Var_x1, P_Var_y1, P_Var_x2, P_Var_y2, P_xlable, P_ylable, P_title, 'r', 'b', 'Lab. Sea conv', 'Lab. Sea nonconv', 'lower left', 'invert_yaxis')
plt.xlim(32.8,35.6)
fig.savefig(dir_figs+str(GCM)+'_salinity_depth_LAB_41W65W_50N65N_conv.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


Salinity_allyears_LAB_B2=Salinity_allyears[:,:,140:156,320:341]
Salinity_allyears_LAB_B2=np.nanmean(Salinity_allyears_LAB_B2,axis=2); Salinity_allyears_LAB_B2=np.nanmean(Salinity_allyears_LAB_B2,axis=0)

Salinity_allyears_LAB_c_B2=Salinity_allyears_lab_c[:,:,140:156,320:341]
Salinity_allyears_LAB_c_B2=np.nanmean(Salinity_allyears_LAB_c_B2,axis=2); Salinity_allyears_LAB_c_B2=np.nanmean(Salinity_allyears_LAB_c_B2,axis=0)

Salinity_allyears_LAB_n_B2=Salinity_allyears_lab_n[:,:,140:156,320:341]
Salinity_allyears_LAB_n_B2=np.nanmean(Salinity_allyears_LAB_n_B2,axis=2); Salinity_allyears_LAB_n_B2=np.nanmean(Salinity_allyears_LAB_n_B2,axis=0)

Salinity_allyears_LAB_c_ave_B2=np.nanmean(Salinity_allyears_LAB_c_B2,axis=1); Salinity_allyears_LAB_n_ave_B2=np.nanmean(Salinity_allyears_LAB_n_B2,axis=1)


P_Var_x1=Salinity_allyears_LAB_c_ave_B2[0:40]
P_Var_x2=Salinity_allyears_LAB_n_ave_B2[0:40]
P_Var_y1=P_Var_y2=Depths[0:40]
P_xlable='Salinity [PSU]'  ;P_ylable='Depth [m]'
P_title='Salinity Ave - Labrador Box2 [41W-20W, 50N-65N] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plotline_2var(P_Var_x1, P_Var_y1, P_Var_x2, P_Var_y2, P_xlable, P_ylable, P_title, 'r', 'b', 'Lab. Sea conv', 'Lab. Sea nonconv', 'lower left', 'invert_yaxis')
plt.xlim(32.8,35.6)
fig.savefig(dir_figs+str(GCM)+'_salinity_depth_LAB_20W41W_50N65N_conv.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

##############   Weddell Sea   ##############
Salinity_allyears_ws_c=Salinity_allyears[R_ii_WS]
Salinity_allyears_ws_n=Salinity_allyears[R_jj_WS]

Salinity_allyears_WS=Salinity_allyears[:,:,20:35,300:]
Salinity_allyears_WS=np.nanmean(Salinity_allyears_WS,axis=2)
Salinity_allyears_WS=np.nanmean(Salinity_allyears_WS,axis=0)

Salinity_allyears_WS_c=Salinity_allyears_ws_c[:,:,20:35,300:]
Salinity_allyears_WS_c=np.nanmean(Salinity_allyears_WS_c,axis=2)
Salinity_allyears_WS_c=np.nanmean(Salinity_allyears_WS_c,axis=0)
Salinity_allyears_WS_n=Salinity_allyears_ws_n[:,:,20:35,300:]
Salinity_allyears_WS_n=np.nanmean(Salinity_allyears_WS_n,axis=2)
Salinity_allyears_WS_n=np.nanmean(Salinity_allyears_WS_n,axis=0)

Salinity_allyears_WS_ave=np.nanmean(Salinity_allyears_WS,axis=1)
Salinity_allyears_WS_c_ave=np.nanmean(Salinity_allyears_WS_c,axis=1)
Salinity_allyears_WS_n_ave=np.nanmean(Salinity_allyears_WS_n,axis=1)

dSal_dZ_WS=np.zeros(( Depths.shape[0]-1))
dSal_dZ_WS_c=np.zeros(( Depths.shape[0]-1))
dSal_dZ_WS_n=np.zeros(( Depths.shape[0]-1))
for ii in range(0,dRho_dZ_LAB.shape[0]):
    dSal_dZ_WS[ii]=(Salinity_allyears_WS_ave[ii]-Salinity_allyears_WS_ave[ii+1])/(Depths[ii]-Depths[ii+1])
    dSal_dZ_WS_c[ii]=(Salinity_allyears_WS_c_ave[ii]-Salinity_allyears_WS_c_ave[ii+1])/(Depths[ii]-Depths[ii+1])
    dSal_dZ_WS_n[ii]=(Salinity_allyears_WS_n_ave[ii]-Salinity_allyears_WS_n_ave[ii+1])/(Depths[ii]-Depths[ii+1])

fig=plt.figure()
im1=plt.plot(Salinity_allyears_WS_ave[0:40], Depths[0:40])
plt.gca().invert_yaxis()
plt.xticks(fontsize = 20); plt.yticks(fontsize = 20)
plt.xlabel('Salinity [kg/m3]', fontsize=18)
plt.ylabel('Depth [m]', fontsize=18)
plt.title('Salinity Ave - Weddell Sea [60W-0W,70S-55S] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_salinity_depth_WS_60W0W_70S55S.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
im1=plt.plot(Salinity_allyears_WS_c_ave[0:40], Depths[0:40], 'r')
im1=plt.plot(Salinity_allyears_WS_n_ave[0:40], Depths[0:40], 'b')
plt.gca().invert_yaxis()
plt.xticks(fontsize = 20); plt.yticks(fontsize = 20)
plt.xlabel('Salinity [kg/m3]', fontsize=18)
plt.ylabel('Depth', fontsize=18)
plt.title('Salinity Ave - Weddell Sea [60W-0W,70S-55S] (red=Wedd. Sea conv, blue=nonconv) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_salinity_depth_WS_60W0W_70S55S_conv.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


fig=plt.figure()
im1=plt.plot(dSal_dZ_WS[0:40], Depths[0:40], 'darkgreen')
plt.gca().invert_yaxis()
plt.xticks(fontsize = 20); plt.yticks(fontsize = 20)
plt.xlabel('dSal/dZ [kg/m3/m]', fontsize=18)
plt.ylabel('Depth [m]', fontsize=18)
plt.title('dSal/dZ - Weddell Sea [60W-0W,70S-55S] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_dSal_dZ_WS_60W0W_70S55S.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
im1=plt.plot(dSal_dZ_WS_c[0:40], Depths[0:40], 'r')
im1=plt.plot(dSal_dZ_WS_n[0:40], Depths[0:40], 'b')
plt.gca().invert_yaxis()
plt.xticks(fontsize = 20); plt.yticks(fontsize = 20)
plt.xlabel('dSal [kg/m3/m]', fontsize=18)
plt.ylabel('Depth [m]', fontsize=18)
plt.title('dSal/dZ - Weddell Sea [60W-0W,70S-55S] (red=Wedd. Sea conv, blue=nonconv) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_dSal_dZ_WS_60W0W_70S55S_conv.png', format='png', dpi=300, transparent=True, bbox_inches='tight')










###############################################################################
########################       ESM2G + CM2MC       ############################
dir_data_in1 = (dir_pwd + '/Results_CM2MC/data/') # Directory to raed raw data from
###############################################################################
from numpy import loadtxt

file_name_dir=(dir_data_in1+ 'depths.dat'); depths = loadtxt(file_name_dir)
file_name_dir=(dir_data_in1+ 'potential_density_labrador_average.dat'); potential_density_labrador_average = loadtxt(file_name_dir)
file_name_dir=(dir_data_in1+ 'potential_density_weddell_average.dat'); potential_density_weddell_average = loadtxt(file_name_dir)
file_name_dir=(dir_data_in1+ 'salinity_labrador_average.dat'); salinity_labrador_average = loadtxt(file_name_dir)
file_name_dir=(dir_data_in1+ 'salinity_weddell_average.dat'); salinity_weddell_average = loadtxt(file_name_dir)
file_name_dir=(dir_data_in1+ 'temperature_labrador_average.dat'); temperature_labrador_average = loadtxt(file_name_dir)
file_name_dir=(dir_data_in1+ 'temperature_weddell_average.dat'); temperature_weddell_average = loadtxt(file_name_dir)

GCM2 = 'GFDL-CM2MC'
dir_figs = (dir_pwd + '/Figures_CM2MC/') # Directory to save figures

P_Var_x1=Density_allyears_LAB_ave[0:40]-1000
P_Var_x2=potential_density_labrador_average[0:21]
P_Var_y1=Depths[0:40]; P_Var_y2=depths[0:21]
P_xlable='Potential Density [kg/m3]'  ;P_ylable='Depth [m]'
P_title='Potential Density Ave - Labrador Sea'
fig = func_plotline_2var(P_Var_x1, P_Var_y1, P_Var_x2, P_Var_y2, P_xlable, P_ylable, P_title, 'r', 'b', 'ESM2G', 'CM2MC', 'best', 'invert_yaxis')
#fig.savefig(dir_figs+str(GCM2)+str(GCM)+'_density_depth_LAB.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

P_Var_x1=Density_allyears_WS_ave[0:40]-1000
P_Var_x2=potential_density_weddell_average[0:21]
P_Var_y1=Depths[0:40]; P_Var_y2=depths[0:21]
P_xlable='Potential Density [kg/m3]'  ;P_ylable='Depth [m]'
P_title='Potential Density Ave - Weddell Sea'
fig = func_plotline_2var(P_Var_x1, P_Var_y1, P_Var_x2, P_Var_y2, P_xlable, P_ylable, P_title, 'r', 'b', 'ESM2G', 'CM2MC', 'best', 'invert_yaxis')
#fig.savefig(dir_figs+str(GCM2)+str(GCM)+'_density_depth_WS.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

P_Var_x1=Temp_allyears_LAB_ave[0:40]-273.15
P_Var_x2=temperature_labrador_average[0:21]
P_Var_y1=Depths[0:40]; P_Var_y2=depths[0:21]
P_xlable='Potential Temperature [C]'  ;P_ylable='Depth [m]'
P_title='Potential Temperature - Labrador Sea'
fig = func_plotline_2var(P_Var_x1, P_Var_y1, P_Var_x2, P_Var_y2, P_xlable, P_ylable, P_title, 'g', 'm', 'ESM2G', 'CM2MC', 'best', 'invert_yaxis')
#fig.savefig(dir_figs+str(GCM2)+str(GCM)+'_temp_depth_LAB.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

P_Var_x1=Temp_allyears_WS_ave[0:40]-273.15
P_Var_x2=temperature_weddell_average[0:21]
P_Var_y1=Depths[0:40]; P_Var_y2=depths[0:21]
P_xlable='Potential Temperature [C]'  ;P_ylable='Depth [m]'
P_title='Potential Temperature - Weddell Sea'
fig = func_plotline_2var(P_Var_x1, P_Var_y1, P_Var_x2, P_Var_y2, P_xlable, P_ylable, P_title, 'g', 'm', 'ESM2G', 'CM2MC', 'best', 'invert_yaxis')
#fig.savefig(dir_figs+str(GCM2)+str(GCM)+'_temp_depth_WS.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

P_Var_x1=Salinity_allyears_LAB_ave[0:40]
P_Var_x2=salinity_labrador_average[0:21]
P_Var_y1=Depths[0:40]; P_Var_y2=depths[0:21]
P_xlable='Salinity [kg/m3]'  ;P_ylable='Depth [m]'
P_title='Salinity - Labrador Sea'
fig = func_plotline_2var(P_Var_x1, P_Var_y1, P_Var_x2, P_Var_y2, P_xlable, P_ylable, P_title, 'g', 'm', 'ESM2G', 'CM2MC', 'best', 'invert_yaxis')
#fig.savefig(dir_figs+str(GCM2)+str(GCM)+'_sal_depth_LAB.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

P_Var_x1=Salinity_allyears_WS_ave[0:40]
P_Var_x2=salinity_weddell_average[0:21]
P_Var_y1=Depths[0:40]; P_Var_y2=depths[0:21]
P_xlable='Salinity [kg/m3]'  ;P_ylable='Depth [m]'
P_title='Salinity - Weddell Sea'
fig = func_plotline_2var(P_Var_x1, P_Var_y1, P_Var_x2, P_Var_y2, P_xlable, P_ylable, P_title, 'g', 'm', 'ESM2G', 'CM2MC', 'best', 'invert_yaxis')
#fig.savefig(dir_figs+str(GCM2)+str(GCM)+'_sal_depth_WS.png', format='png', dpi=300, transparent=True, bbox_inches='tight')



###############################################################################
from numpy import loadtxt

file_name_dir=(dir_data_in1+ 'latitudes.dat'); latitudes_CM2MC = loadtxt(file_name_dir)
file_name_dir=(dir_data_in1+ 'longitudes.dat'); longitudes_CM2MC = loadtxt(file_name_dir)

file_dire=dir_data_in1+'amoc.nc'
dset_t = Dataset(file_dire)  #  dset_t.variables
AMOC_CM2MC=np.asarray(dset_t.variables['AMOC'][:])
#AMOC_CM2MC=np.transpose(AMOC_CM2MC)

AMOC_CM2MC=np.float64(AMOC_CM2MC)

file_dire=dir_data_in1+'lab_sst_index.nc'
dset_t = Dataset(file_dire)  #  dset_t.variables
LAB_index_CM2MC=np.asarray(dset_t.variables['SSTLAB'][:])

LAB_index_CM2MC_rm = copy.deepcopy(LAB_index_CM2MC)
LAB_index_CM2MC_rm=runningMeanFast(LAB_index_CM2MC_rm, 10)

###LAB plots###
lag_time=40

fig=plt.figure()
data=[]
for i in range(len(AMOC_CM2MC[0][3:73])):
    stream=runningMeanFast(AMOC_CM2MC[:,i+lag_time], 10)
    r=lag_cor_data(LAB_index_CM2MC_rm,stream,lag_time)
    data.append(r)
data=np.asarray(data)
data=np.asarray(data)
xxx=np.linspace(-lag_time,lag_time+1, 2*lag_time)
yyy=latitudes_CM2MC[3:73]

xxx,yyy=np.meshgrid(xxx,yyy)
cmap_limit=np.nanmax(np.abs(data))
levels=np.linspace(-cmap_limit,cmap_limit,200)
im=plt.contourf(xxx,yyy,data,levels,latlon=True, cmap=plt.cm.jet)
cb = plt.colorbar(im)
plt.xlabel('Year lag', fontsize=18)
plt.ylabel('Longitude', fontsize=18)
plt.title('Labrador peak convection vs. AMOC (max streamfunction) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM2), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
#fig.savefig(dir_figs+str(GCM2)+'_lagcor_convec_AMOC_max_LAB.png', format='png', dpi=300, transparent=True, bbox_inches='tight')












#fig=plt.figure()
#im1=plt.contour(Lat_regrid_1D[10:60], Depths[0:31], R_conv_WS_60W0E[0:31,10:60]-1000, 30, colors='r')
#plt.clabel(im1, fontsize=8, inline=1)
#plt.gca().invert_yaxis()
#im2=plt.contour(Lat_regrid_1D[10:60], Depths[0:31], R_nonconv_WS_60W0E[0:31,10:60]-1000, 30, colors='b')
#plt.clabel(im2, fontsize=8, inline=1)
#plt.xlabel('Latitude', fontsize=18)
#plt.ylabel('Depth', fontsize=18)
#plt.title('Potential Density Contours (red=conv, blue=nonconv) - 60W-0E - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
#mng = plt.get_current_fig_manager()
#mng.window.showMaximized() # Maximizes the plot window to save figures in full
#
#fig=plt.figure()
#im1=plt.contour(Lat_regrid_1D[10:60], Depths[0:31], R_conv_WS_All[0:31,10:60]-1000, 30, colors='r')
#plt.clabel(im1, fontsize=8, inline=1)
#plt.gca().invert_yaxis()
#im2=plt.contour(Lat_regrid_1D[10:60], Depths[0:31], R_nonconv_WS_All[0:31,10:60]-1000, 30, colors='b')
#plt.clabel(im2, fontsize=8, inline=1)
#plt.xlabel('Latitude', fontsize=18)
#plt.ylabel('Depth', fontsize=18)
#plt.title('Potential Density Contours (red=conv, blue=nonconv) - All Southern Ocean - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
#mng = plt.get_current_fig_manager()
#mng.window.showMaximized() # Maximizes the plot window to save figures in full
#
#fig=plt.figure()
#im1=plt.contour(Lat_regrid_1D, Depths[0:31], R_conv_WS_All[0:31,:]-1000, 40, colors='b')
#plt.clabel(im1, fontsize=8, inline=1)
#plt.gca().invert_yaxis()
#plt.xlabel('Latitude', fontsize=18)
#plt.ylabel('Depth', fontsize=18)
#plt.title('Potential Density Contours - All Ocean - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
#mng = plt.get_current_fig_manager()
#mng.window.showMaximized() # Maximizes the plot window to save figures in full
#
#fig=plt.figure()
#im1=plt.contour(Lat_regrid_1D, Depths[0:31], R_conv_WS_60W0E[0:31,:]-1000, 40, colors='b')
#plt.clabel(im1, fontsize=8, inline=1)
#plt.gca().invert_yaxis()
#plt.xlabel('Latitude', fontsize=18)
#plt.ylabel('Depth', fontsize=18)
#plt.title('Potential Density Contours - 60W-0E - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
#mng = plt.get_current_fig_manager()
#mng.window.showMaximized() # Maximizes the plot window to save figures in full






