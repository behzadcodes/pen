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
import os
import shelve
########################################
var_list_long=np.load('var_list_long.npy')
var_list_short=np.load('var_list_short.npy')

dir_pwd = os.getcwd() # Gets the current directory (and in which the code is placed)
dir_figs = (dir_pwd + '/Figures_MultiGCM/') # Directory to save figures
dir_results = (dir_pwd + '/Results/') # Directory to save results

GCM_Names = ['GFDL-ESM2G', 'GFDL-ESM2M','HadGEM2-ES','IPSL-CM5A-MR','IPSL-CM5A-LR','MPI-ESM-LR','MRI-ESM1','CESM1-BGC','CanESM2','MIROC-ESM-CHEM'] ##  ,'CMCC-CESM','CNRM-CM5','GFDL-ESM2G','GFDL-ESM2M','GISS-E2-H-CC','GISS-E2-R-CC','HadGEM2-CC','HadGEM2-ES','IPSL-CM5B-LR','MIROC-ESM','MIROC-ESM-CHEM','MPI-ESM-MR','MRI-ESM1','IPSL-CM5B-LR','NorESM1-ME']
#GCM_Names = ['GFDL-ESM2G', 'GFDL-ESM2M','IPSL-CM5A-MR']

### Regrdridding calculations ###
# creating new coordinate grid, same which was used in interpolation in data processing code
lat_n_regrid, lon_n_regrid = 180, 360 # Number of Lat and Lon elements in the regridded data
lon_min_regrid, lon_max_regrid = 0, 360 # Min and Max value of Lon in the regridded data
lat_min_regrid, lat_max_regrid = -90, 90 # Min and Max value of Lat in the regridded data

Lat_regrid_1D, Lon_regrid_1D, Lat_bound_regrid, Lon_bound_regrid = func_latlon_regrid_eq(lat_n_regrid, lon_n_regrid, lat_min_regrid, lat_max_regrid, lon_min_regrid, lon_max_regrid)
Lon_regrid_2D, Lat_regrid_2D = np.meshgrid(Lon_regrid_1D, Lat_regrid_1D)

##################################################
#%% Convection Index Plots - WS and LAB ###
###########################################
### Weddel Sea Plots ####

n_r=4 # Number of rows for subplot
n_c=3 # Number of columns for subplot
n_range=list(range(len(GCM_Names)))

fig=plt.figure()
levels=np.linspace(2000,5000,100)
for ii in n_range:
    
    GCM = GCM_Names[ii]    
    filename_in = (dir_results + 'AllResults_'+GCM+'_500yr.out') # Directory to save processed data
    my_shelf = shelve.open(filename_in)
    globals()['year_start']=my_shelf['year_start']
    globals()['year_end']=my_shelf['year_end']    
    globals()['MLD_average_WS']=my_shelf['MLD_average_WS']
    globals()['lon_mld_WS']=my_shelf['lon_mld_WS']
    globals()['lat_mld_WS']=my_shelf['lat_mld_WS']
    globals()['WS_indeces_lonlat']=my_shelf['WS_indeces_lonlat']
    
    ax = fig.add_subplot(n_r,n_c,ii+1)
    m = Basemap( projection='spstere',lon_0=0,boundinglat=-30)
    m.drawcoastlines(linewidth=1.25)
    m.fillcontinents(color='0.8')
    m.drawparallels(np.arange(-90,90,20), labels=[1,1,0,1])
    m.drawmeridians(np.arange(0,360,30), labels=[1,1,0,1])
    im=m.contourf(lon_mld_WS,lat_mld_WS,MLD_average_WS, levels,latlon=True, cmap=plt.cm.jet)
    m.scatter(lon_mld_WS[WS_indeces_lonlat],lat_mld_WS[WS_indeces_lonlat], 0.05, c="k", latlon=True)
    plt.title(GCM_Names[ii]+' (yr '+str(year_start)+'-'+str(year_end)+')', fontsize=12)

plt.subplots_adjust(left=0.25, bottom=0.05, right=0.75, top=0.9, hspace=0.3, wspace=0.01) # the amount of height/width reserved for space between subplots
cbar_ax = fig.add_axes([0.78, 0.1, 0.015, 0.8]) # cax = [left position, bottom postion, width, height] 
fig.colorbar(im, cax=cbar_ax)
plt.suptitle('Average MLD [m] (Weddell Sea)', fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+'AllModels_average_MLD_WS.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
for ii in n_range:
    
    GCM = GCM_Names[ii]    
    filename_in = (dir_results + 'AllResults_'+GCM+'_500yr.out') # Directory to save processed data
    my_shelf = shelve.open(filename_in)
    globals()['year_start']=my_shelf['year_start']
    globals()['year_end']=my_shelf['year_end']    
    globals()['Conv_area_timeseries_WS']=my_shelf['Conv_area_timeseries_WS']
    
    ax = fig.add_subplot(n_r,n_c,ii+1)
    years=np.linspace(year_start,year_end,year_end-year_start+1)
    #print (years,area)
    plt.plot(years, Conv_area_timeseries_WS / 1e6,'b')
    if ii+1 >= n_range[-n_c+1]: # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        plt.xlabel('Years', fontsize=18)
    if ii+1 == n_c+1: # Adds longitude ranges only to the last subplots that appear at the bottom of plot     
        plt.ylabel('Deep Convection Area (km2)', fontsize=18)
    plt.title(GCM_Names[ii]+' (yr '+str(year_start)+'-'+str(year_end)+')', fontsize=12)

plt.subplots_adjust(hspace=0.4, wspace=0.2) # the amount of height/width reserved for space between subplots
plt.suptitle('Deep Convection Area [km2] (MLD>2000m) - Weddell Sea', fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+'AllModels_deep_conv_area_WS.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


#fig=plt.figure()
#levels=np.linspace(-1.5, 2,100)
#for ii in n_range:
#    
#    GCM = GCM_Names[ii]    
#    filename_in = (dir_results + 'AllResults_'+GCM+'_500yr.out') # Directory to save processed data
#    my_shelf = shelve.open(filename_in)
#    globals()['year_start']=my_shelf['year_start']
#    globals()['year_end']=my_shelf['year_end']    
#    globals()['WS_index_time_depth']=my_shelf['WS_index_time_depth']
#    globals()['conv_index_depth_ws']=my_shelf['conv_index_depth_ws']
#    globals()['Depths']=my_shelf['Depths']
#    
#    ax = fig.add_subplot(n_r,n_c,ii+1)
#    years=np.linspace(year_start,year_end,year_end-year_start+1)
#    #print (years,area)
#    im = plt.contourf(years, Depths, np.transpose(WS_index_time_depth-273.15), levels, cmap=plt.cm.jet) # contour(X,Y,Z,N) - N shows the number of contour levels
#    plt.gca().invert_yaxis()
#    l = plt.axhline(y=conv_index_depth_ws)
#    if ii+1 >= n_range[-n_c+1]: # Adds longitude ranges only to the last subplots that appear at the bottom of plot
#        plt.xlabel('Years', fontsize=18)
#    if ii+1 == n_c+1: # Adds longitude ranges only to the last subplots that appear at the bottom of plot     
#        plt.ylabel('Depth [m]', fontsize=18)
#    plt.title(GCM_Names[ii]+' (yr '+str(year_start)+'-'+str(year_end)+')', fontsize=12)
#
#plt.subplots_adjust(hspace=0.3, wspace=0.2) # the amount of height/width reserved for space between subplots
#cbar_ax = fig.add_axes([0.93, 0.1, 0.015, 0.8])
#fig.colorbar(im, cax=cbar_ax)
#plt.suptitle('Temperature timeseries (C) in Convection Area (Weddell Sea)', fontsize=18)
#mng = plt.get_current_fig_manager()
#mng.window.showMaximized() # Maximizes the plot window to save figures in full
#fig.savefig(dir_figs+'AllModels_index_time_depth_WS_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
for ii in n_range:
    
    GCM = GCM_Names[ii]    
    filename_in = (dir_results + 'AllResults_'+GCM+'_500yr.out') # Directory to save processed data
    my_shelf = shelve.open(filename_in)
    globals()['year_start']=my_shelf['year_start']
    globals()['year_end']=my_shelf['year_end']    
    globals()['WS_index_time_depth']=my_shelf['WS_index_time_depth']
    globals()['conv_index_depth_ws']=my_shelf['conv_index_depth_ws']
    globals()['Depths']=my_shelf['Depths']
    
    ax = fig.add_subplot(n_r,n_c,ii+1)
    years=np.linspace(year_start,year_end,year_end-year_start+1)
    #print (years,area)
    im = plt.contourf(years, Depths, np.transpose(WS_index_time_depth-273.15), 60, cmap=plt.cm.jet) # contour(X,Y,Z,N) - N shows the number of contour levels
    plt.gca().invert_yaxis()
    l = plt.axhline(y=conv_index_depth_ws)
    cbar = m.colorbar(im,"right", size="3%", pad="2%", extend='max') # extend='both' will extend the colorbar in both sides (upper side and down side)
    if ii+1 >= n_range[-n_c+1]: # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        plt.xlabel('Years', fontsize=18)
    if ii+1 == n_c+1: # Adds longitude ranges only to the last subplots that appear at the bottom of plot     
        plt.ylabel('Depth [m]', fontsize=18)
    plt.title(GCM_Names[ii]+' (yr '+str(year_start)+'-'+str(year_end)+')', fontsize=12)

plt.subplots_adjust(hspace=0.4, wspace=0.35) # the amount of height/width reserved for space between subplots
plt.suptitle('Temperature timeseries (C) in Convection Area (Weddell Sea)', fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+'AllModels_index_time_depth_WS.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
for ii in n_range:
    
    GCM = GCM_Names[ii]    
    filename_in = (dir_results + 'AllResults_'+GCM+'_500yr.out') # Directory to save processed data
    my_shelf = shelve.open(filename_in)
    globals()['year_start']=my_shelf['year_start']
    globals()['year_end']=my_shelf['year_end']    
    globals()['WS_index_time_depth']=my_shelf['WS_index_time_depth']
    globals()['conv_index_depth_ws']=my_shelf['conv_index_depth_ws']
    globals()['Depths']=my_shelf['Depths']
    
    ax = fig.add_subplot(n_r,n_c,ii+1)
    years=np.linspace(year_start,year_end,year_end-year_start+1)
    #print (years,area)
    im = plt.contourf(years, Depths[0:11], np.transpose(WS_index_time_depth[:,0:11]-273.15), 40, cmap=plt.cm.jet)
    plt.gca().invert_yaxis()
    cbar = m.colorbar(im,"right", size="3%", pad="2%", extend='max') # extend='both' will extend the colorbar in both sides (upper side and down side)
    if ii+1 >= n_range[-n_c+1]: # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        plt.xlabel('Years', fontsize=18)
    if ii+1 == n_c+1: # Adds longitude ranges only to the last subplots that appear at the bottom of plot     
        plt.ylabel('Depth [m]', fontsize=18)
    plt.title(GCM_Names[ii]+' (yr '+str(year_start)+'-'+str(year_end)+')', fontsize=12)

plt.subplots_adjust(hspace=0.4, wspace=0.35) # the amount of height/width reserved for space between subplots
plt.suptitle('Temperature timeseries (C) in Convection Area (Weddell Sea)', fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+'AllModels_index_time_depth_surface_WS.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

############################
#### Labrador Sea Plots ####
fig=plt.figure()
levels=np.linspace(1000,3000,100)
for ii in n_range:
    
    GCM = GCM_Names[ii]    
    filename_in = (dir_results + 'AllResults_'+GCM+'_500yr.out') # Directory to save processed data
    my_shelf = shelve.open(filename_in)
    globals()['year_start']=my_shelf['year_start']
    globals()['year_end']=my_shelf['year_end']    
    globals()['MLD_average_LAB']=my_shelf['MLD_average_LAB']
    globals()['lon_mld_LAB']=my_shelf['lon_mld_LAB']
    globals()['lat_mld_LAB']=my_shelf['lat_mld_LAB']
    globals()['LAB_indeces_lonlat']=my_shelf['LAB_indeces_lonlat']
    
    ax = fig.add_subplot(n_r,n_c,ii+1)
    m = Basemap( projection='npstere',lon_0=0,boundinglat=30)
    m.drawcoastlines(linewidth=1.25)
    m.fillcontinents(color='0.8')
    m.drawparallels(np.arange(-90,90,20), labels=[1,1,0,1])
    m.drawmeridians(np.arange(0,360,30), labels=[1,1,0,1])
    im=m.contourf(lon_mld_LAB,lat_mld_LAB,MLD_average_LAB, levels,latlon=True, cmap=plt.cm.jet)
    m.scatter(lon_mld_LAB[LAB_indeces_lonlat],lat_mld_LAB[LAB_indeces_lonlat], 0.05, c="k", latlon=True)
    plt.title(GCM_Names[ii]+' (yr '+str(year_start)+'-'+str(year_end)+')', fontsize=12)

plt.subplots_adjust(left=0.25, bottom=0.05, right=0.75, top=0.9, hspace=0.3, wspace=0.01) # the amount of height/width reserved for space between subplots
cbar_ax = fig.add_axes([0.78, 0.1, 0.015, 0.8]) # cax = [left position, bottom postion, width, height] 
fig.colorbar(im, cax=cbar_ax)
plt.suptitle('Average MLD [m] (Labrador Sea)', fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+'AllModels_average_MLD_LAB.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
for ii in n_range:
    
    GCM = GCM_Names[ii]    
    filename_in = (dir_results + 'AllResults_'+GCM+'_500yr.out') # Directory to save processed data
    my_shelf = shelve.open(filename_in)
    globals()['year_start']=my_shelf['year_start']
    globals()['year_end']=my_shelf['year_end']    
    globals()['Conv_area_timeseries_LAB']=my_shelf['Conv_area_timeseries_LAB']
    
    ax = fig.add_subplot(n_r,n_c,ii+1)
    years=np.linspace(year_start,year_end,year_end-year_start+1)
    #print (years,area)
    plt.plot(years, Conv_area_timeseries_LAB / 1e6,'b')
    if ii+1 >= n_range[-n_c+1]: # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        plt.xlabel('Years', fontsize=18)
    if ii+1 == n_c+1: # Adds longitude ranges only to the last subplots that appear at the bottom of plot     
        plt.ylabel('Deep Convection Area (km2)', fontsize=18)
    plt.title(GCM_Names[ii]+' (yr '+str(year_start)+'-'+str(year_end)+')', fontsize=12)

plt.subplots_adjust(hspace=0.4, wspace=0.2) # the amount of height/width reserved for space between subplots
plt.suptitle('Deep Convection Area [km2] (MLD>1000m) - Labrador Sea', fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+'AllModels_deep_conv_area_LAB.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

#fig=plt.figure()
#levels=np.linspace(-1.5, 2,100)
#norm = matplotlib.colors.BoundaryNorm(boundaries=bounds, ncolors=256)
#for ii in n_range:
#    
#    GCM = GCM_Names[ii]    
#    filename_in = (dir_results + 'AllResults_'+GCM+'_500yr.out') # Directory to save processed data
#    my_shelf = shelve.open(filename_in)
#    globals()['year_start']=my_shelf['year_start']
#    globals()['year_end']=my_shelf['year_end']    
#    globals()['LAB_index_time_depth']=my_shelf['LAB_index_time_depth']
#    globals()['conv_index_depth_lab']=my_shelf['conv_index_depth_lab']
#    globals()['Depths']=my_shelf['Depths']
#    
#    ax = fig.add_subplot(n_r,n_c,ii+1)
#    years=np.linspace(year_start,year_end,year_end-year_start+1)
#    #print (years,area)
#    im = plt.contourf(years, Depths, np.transpose(LAB_index_time_depth-273.15), levels, cmap=plt.cm.jet) # contour(X,Y,Z,N) - N shows the number of contour levels
#    plt.gca().invert_yaxis()
#    l = plt.axhline(y=conv_index_depth_ws)
#    if ii+1 >= n_range[-n_c+1]: # Adds longitude ranges only to the last subplots that appear at the bottom of plot
#        plt.xlabel('Years', fontsize=18)
#    if ii+1 == n_c+1: # Adds longitude ranges only to the last subplots that appear at the bottom of plot     
#        plt.ylabel('Depth [m]', fontsize=18)
#    plt.title(GCM_Names[ii]+' (yr '+str(year_start)+'-'+str(year_end)+')', fontsize=12)
#
#plt.subplots_adjust(hspace=0.3, wspace=0.2) # the amount of height/width reserved for space between subplots
#cbar_ax = fig.add_axes([0.93, 0.1, 0.015, 0.8])
#fig.colorbar(im, cax=cbar_ax)
#plt.suptitle('Temperature timeseries (C) in Convection Area (Labrador Sea)', fontsize=18)
#mng = plt.get_current_fig_manager()
#mng.window.showMaximized() # Maximizes the plot window to save figures in full
#fig.savefig(dir_figs+'AllModels_index_time_depth_LAB_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
for ii in n_range:
    
    GCM = GCM_Names[ii]    
    filename_in = (dir_results + 'AllResults_'+GCM+'_500yr.out') # Directory to save processed data
    my_shelf = shelve.open(filename_in)
    globals()['year_start']=my_shelf['year_start']
    globals()['year_end']=my_shelf['year_end']    
    globals()['LAB_index_time_depth']=my_shelf['LAB_index_time_depth']
    globals()['conv_index_depth_lab']=my_shelf['conv_index_depth_lab']
    globals()['Depths']=my_shelf['Depths']
    
    ax = fig.add_subplot(n_r,n_c,ii+1)
    years=np.linspace(year_start,year_end,year_end-year_start+1)
    #print (years,area)
    im = plt.contourf(years, Depths, np.transpose(LAB_index_time_depth-273.15), 60, cmap=plt.cm.jet) # contour(X,Y,Z,N) - N shows the number of contour levels
    plt.gca().invert_yaxis()
    l = plt.axhline(y=conv_index_depth_lab)
    cbar = m.colorbar(im,"right", size="3%", pad="2%", extend='max') # extend='both' will extend the colorbar in both sides (upper side and down side)
    if ii+1 >= n_range[-n_c+1]: # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        plt.xlabel('Years', fontsize=18)
    if ii+1 == n_c+1: # Adds longitude ranges only to the last subplots that appear at the bottom of plot     
        plt.ylabel('Depth [m]', fontsize=18)
    plt.title(GCM_Names[ii]+' (yr '+str(year_start)+'-'+str(year_end)+')', fontsize=12)

plt.subplots_adjust(hspace=0.4, wspace=0.35) # the amount of height/width reserved for space between subplots
plt.suptitle('Temperature timeseries (C) in Convection Area (Labrador Sea)', fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+'AllModels_index_time_depth_LAB.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
for ii in n_range:
    
    GCM = GCM_Names[ii]    
    filename_in = (dir_results + 'AllResults_'+GCM+'_500yr.out') # Directory to save processed data
    my_shelf = shelve.open(filename_in)
    globals()['year_start']=my_shelf['year_start']
    globals()['year_end']=my_shelf['year_end']    
    globals()['LAB_index_time_depth']=my_shelf['LAB_index_time_depth']
    globals()['conv_index_depth_lab']=my_shelf['conv_index_depth_lab']
    globals()['Depths']=my_shelf['Depths']
    
    ax = fig.add_subplot(n_r,n_c,ii+1)
    years=np.linspace(year_start,year_end,year_end-year_start+1)
    #print (years,area)
    im = plt.contourf(years, Depths[0:11], np.transpose(LAB_index_time_depth[:,0:11]-273.15), 40, cmap=plt.cm.jet)
    plt.gca().invert_yaxis()
    cbar = m.colorbar(im,"right", size="3%", pad="2%", extend='max') # extend='both' will extend the colorbar in both sides (upper side and down side)
    if ii+1 >= n_range[-n_c+1]: # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        plt.xlabel('Years', fontsize=18)
    if ii+1 == n_c+1: # Adds longitude ranges only to the last subplots that appear at the bottom of plot     
        plt.ylabel('Depth [m]', fontsize=18)
    plt.title(GCM_Names[ii]+' (yr '+str(year_start)+'-'+str(year_end)+')', fontsize=12)

plt.subplots_adjust(hspace=0.4, wspace=0.35) # the amount of height/width reserved for space between subplots
plt.suptitle('Temperature timeseries (C) in Convection Area (Labrador Sea)', fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+'AllModels_index_time_depth_surface_LAB.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

###########################################
### Convection Index Time Series Plots ####
fig=plt.figure()
for ii in n_range:
    
    GCM = GCM_Names[ii]    
    filename_in = (dir_results + 'AllResults_'+GCM+'_500yr.out') # Directory to save processed data
    my_shelf = shelve.open(filename_in)
    globals()['year_start']=my_shelf['year_start']
    globals()['year_end']=my_shelf['year_end']    
    globals()['WS_index_norm']=my_shelf['WS_index_norm']
    
    ax = fig.add_subplot(n_r,n_c,ii+1)
    years=np.linspace(year_start,year_end,year_end-year_start+1)
    plt.plot(years,WS_index_norm,'b')
    l = plt.axhline(y=0)
    if ii+1 >= n_range[-n_c+1]: # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        plt.xlabel('Years', fontsize=18)
    if ii+1 == n_c+1: # Adds longitude ranges only to the last subplots that appear at the bottom of plot     
        plt.ylabel('Convection index (°C)', fontsize=18)
    plt.title(GCM_Names[ii]+' (yr '+str(year_start)+'-'+str(year_end)+')', fontsize=12)

plt.subplots_adjust(hspace=0.4, wspace=0.2) # the amount of height/width reserved for space between subplots
plt.suptitle('Weddel Sea convection index (normalized)', fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+'AllModels_convec_index_WS.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
for ii in n_range:
    
    GCM = GCM_Names[ii]    
    filename_in = (dir_results + 'AllResults_'+GCM+'_500yr.out') # Directory to save processed data
    my_shelf = shelve.open(filename_in)
    globals()['year_start']=my_shelf['year_start']
    globals()['year_end']=my_shelf['year_end']    
    globals()['LAB_index_norm']=my_shelf['LAB_index_norm']
    
    ax = fig.add_subplot(n_r,n_c,ii+1)
    years=np.linspace(year_start,year_end,year_end-year_start+1)
    plt.plot(years,LAB_index_norm,'b')
    l = plt.axhline(y=0)
    if ii+1 >= n_range[-n_c+1]: # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        plt.xlabel('Years', fontsize=18)
    if ii+1 == n_c+1: # Adds longitude ranges only to the last subplots that appear at the bottom of plot     
        plt.ylabel('Convection index (°C)', fontsize=18)
    plt.title(GCM_Names[ii]+' (yr '+str(year_start)+'-'+str(year_end)+')', fontsize=12)

plt.subplots_adjust(hspace=0.4, wspace=0.2) # the amount of height/width reserved for space between subplots
plt.suptitle('Labrador Sea convection index (normalized)', fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+'AllModels_convec_index_LAB.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

#######################################
#%% Streamfunction Plots ####
#plot #1, just transport
fig=plt.figure()
for ii in n_range:
    
    GCM = GCM_Names[ii]    
    filename_in = (dir_results + 'AllResults_'+GCM+'_500yr.out') # Directory to save processed data
    my_shelf = shelve.open(filename_in)
    globals()['year_start']=my_shelf['year_start']
    globals()['year_end']=my_shelf['year_end']    
    globals()['Transport_lon_final_mean']=my_shelf['Transport_lon_final_mean']
    globals()['Depths']=my_shelf['Depths']
    
    ax = fig.add_subplot(n_r,n_c,ii+1)
    years=np.linspace(year_start,year_end,year_end-year_start+1)
    varr=Transport_lon_final_mean
    varr=np.nanmax([ np.absolute(float("{0:.0f}".format(np.nanpercentile(varr, 99)))) ,  np.absolute(float("{0:.0f}".format(np.nanpercentile(varr, 1)))) ])
    levels=np.linspace(-varr,varr,101) 
    #print (years,area)
    im=plt.contourf(Lat_regrid_1D, Depths, Transport_lon_final_mean, levels, cmap=plt.cm.RdBu_r, extend = 'both')
    plt.gca().invert_yaxis()
    cbar = m.colorbar(im,"right", size="3%", pad="2%", extend='max') # extend='both' will extend the colorbar in both sides (upper side and down side)
    if ii+1 >= n_range[-n_c+1]: # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        plt.xlabel('Latitude', fontsize=18)
    if ii+1 == n_c+1: # Adds longitude ranges only to the last subplots that appear at the bottom of plot     
        plt.ylabel('Depth [m]', fontsize=18)
    plt.title(GCM_Names[ii]+' (yr '+str(year_start)+'-'+str(year_end)+')', fontsize=12)

plt.subplots_adjust(hspace=0.4, wspace=0.35) # the amount of height/width reserved for space between subplots
plt.suptitle('Mean northward transport (Atlantic) [Sv]', fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+'AllModels_Transport_lon_final_mean.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
levels=np.linspace(-12,12,101)  
for ii in n_range:
    
    GCM = GCM_Names[ii]    
    filename_in = (dir_results + 'AllResults_'+GCM+'_500yr.out') # Directory to save processed data
    my_shelf = shelve.open(filename_in)
    globals()['year_start']=my_shelf['year_start']
    globals()['year_end']=my_shelf['year_end']    
    globals()['Transport_lon_final_mean']=my_shelf['Transport_lon_final_mean']
    globals()['Depths']=my_shelf['Depths']
    
    ax = fig.add_subplot(n_r,n_c,ii+1)
    im=plt.contourf(Lat_regrid_1D, Depths, Transport_lon_final_mean, levels, cmap=plt.cm.RdBu_r, extend = 'both')
    plt.gca().invert_yaxis()
    if ii+1 >= n_range[-n_c+1]: # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        plt.xlabel('Latitude', fontsize=18)
    if ii+1 == n_c+1: # Adds longitude ranges only to the last subplots that appear at the bottom of plot     
        plt.ylabel('Depth [m]', fontsize=18)
    plt.title(GCM_Names[ii]+' (yr '+str(year_start)+'-'+str(year_end)+')', fontsize=12)

plt.subplots_adjust(hspace=0.4, wspace=0.2) # the amount of height/width reserved for space between subplots
cbar_ax = fig.add_axes([0.93, 0.1, 0.015, 0.8])
fig.colorbar(im, cax=cbar_ax)
plt.suptitle('Mean northward transport (Atlantic) [Sv]', fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+'AllModels_Transport_lon_final_mean_2.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

#plot #2, streamfunction
fig=plt.figure()
for ii in n_range:
    
    GCM = GCM_Names[ii]    
    filename_in = (dir_results + 'AllResults_'+GCM+'_500yr.out') # Directory to save processed data
    my_shelf = shelve.open(filename_in)
    globals()['year_start']=my_shelf['year_start']
    globals()['year_end']=my_shelf['year_end']    
    globals()['Stream_function_ave']=my_shelf['Stream_function_ave']
    globals()['Depths']=my_shelf['Depths']
    
    ax = fig.add_subplot(n_r,n_c,ii+1)
    years=np.linspace(year_start,year_end,year_end-year_start+1)
    varr=Stream_function_ave
    varr=np.nanmax([ np.absolute(float("{0:.0f}".format(np.nanpercentile(varr, 95)))) ,  np.absolute(float("{0:.0f}".format(np.nanpercentile(varr, 5)))) ])
    levels=np.linspace(-varr,varr,101) 
    #print (years,area)
    im=plt.contourf(Lat_regrid_1D, Depths, Stream_function_ave, levels, cmap=plt.cm.RdBu_r, extend = 'both')
    plt.gca().invert_yaxis()
    cbar = m.colorbar(im,"right", size="3%", pad="2%", extend='max') # extend='both' will extend the colorbar in both sides (upper side and down side)
    if ii+1 >= n_range[-n_c+1]: # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        plt.xlabel('Latitude', fontsize=18)
    if ii+1 == n_c+1: # Adds longitude ranges only to the last subplots that appear at the bottom of plot     
        plt.ylabel('Depth [m]', fontsize=18)
    plt.title(GCM_Names[ii]+' (yr '+str(year_start)+'-'+str(year_end)+')', fontsize=12)

plt.subplots_adjust(hspace=0.4, wspace=0.35) # the amount of height/width reserved for space between subplots
plt.suptitle('Stream function (Atlantic) [Sv]', fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+'AllModels_Stream_function_ave.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
levels=np.linspace(-35,35,101)  
for ii in n_range:
    
    GCM = GCM_Names[ii]    
    filename_in = (dir_results + 'AllResults_'+GCM+'_500yr.out') # Directory to save processed data
    my_shelf = shelve.open(filename_in)
    globals()['year_start']=my_shelf['year_start']
    globals()['year_end']=my_shelf['year_end']    
    globals()['Stream_function_ave']=my_shelf['Stream_function_ave']
    globals()['Depths']=my_shelf['Depths']
    
    ax = fig.add_subplot(n_r,n_c,ii+1)
    im=plt.contourf(Lat_regrid_1D, Depths, Stream_function_ave, levels, cmap=plt.cm.RdBu_r, extend = 'both')
    plt.gca().invert_yaxis()
    if ii+1 >= n_range[-n_c+1]: # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        plt.xlabel('Latitude', fontsize=18)
    if ii+1 == n_c+1: # Adds longitude ranges only to the last subplots that appear at the bottom of plot     
        plt.ylabel('Depth [m]', fontsize=18)
    plt.title(GCM_Names[ii]+' (yr '+str(year_start)+'-'+str(year_end)+')', fontsize=12)

plt.subplots_adjust(hspace=0.4, wspace=0.2) # the amount of height/width reserved for space between subplots
cbar_ax = fig.add_axes([0.93, 0.1, 0.015, 0.8])
fig.colorbar(im, cax=cbar_ax)
plt.suptitle('Stream function (Atlantic) [Sv]', fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+'AllModels_Stream_function_ave_2.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

#plot #3, map of mean transport upper 1000m with indeces of max for each lat
fig=plt.figure()
for ii in n_range:
    
    GCM = GCM_Names[ii]    
    filename_in = (dir_results + 'AllResults_'+GCM+'_500yr.out') # Directory to save processed data
    my_shelf = shelve.open(filename_in)
    globals()['year_start']=my_shelf['year_start']
    globals()['year_end']=my_shelf['year_end']    
    globals()['transport_0_1000_mean']=my_shelf['transport_0_1000_mean']
    globals()['lon_stream']=my_shelf['lon_stream']
    globals()['lat_stream']=my_shelf['lat_stream']
    globals()['ii_max']=my_shelf['ii_max']
    globals()['jj_max']=my_shelf['jj_max']

    ax = fig.add_subplot(n_r,n_c,ii+1)
    m = Basemap(llcrnrlon=-120.,llcrnrlat=-90.,urcrnrlon=60.,urcrnrlat=90., projection='mill',lon_0=0)
    m.fillcontinents(color='0.8')
    m.drawparallels(np.arange(-90,90,30), labels=[1,1,0,1])
    m.drawmeridians(np.arange(0,360,60), labels=[1,1,0,1])
    bounds_max=float("{0:.01f}".format(np.nanpercentile(transport_0_1000_mean, 99.9))) # Upper bound of plotted values to be used for colorbar, which is 98th percentile of data, with 2 decimals
    levels=np.linspace(-1*bounds_max,bounds_max,101)  
    im=m.contourf(lon_stream,lat_stream,transport_0_1000_mean,levels,latlon=True, cmap=plt.cm.cool)
    m.scatter(lon_stream[ii_max,jj_max],lat_stream[ii_max,jj_max],0.05,c='k',latlon=True)
    cbar = m.colorbar(im,"right", size="8%", pad="40%", extend = 'both') # extend='both' will extend the colorbar in both sides (upper side and down side)
    plt.title(GCM_Names[ii]+' (yr '+str(year_start)+'-'+str(year_end)+')', fontsize=12)

plt.subplots_adjust(left=0.25, bottom=0.05, right=0.75, top=0.9, hspace=0.3, wspace=0.3) # the amount of height/width reserved for space between subplots
plt.suptitle('Mean transport upper 1000m [Sv] (Stippling= max transport locations)', fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+'AllModels_transport_0_1000_mean_stiple.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

#plot #4, map of mean transport 2000m-3000m with indeces of max for each lat
fig=plt.figure()
for ii in n_range:
    
    GCM = GCM_Names[ii]    
    filename_in = (dir_results + 'AllResults_'+GCM+'_500yr.out') # Directory to save processed data
    my_shelf = shelve.open(filename_in)
    globals()['year_start']=my_shelf['year_start']
    globals()['year_end']=my_shelf['year_end']    
    globals()['transport_2000_3000_mean']=my_shelf['transport_2000_3000_mean']
    globals()['lon_stream']=my_shelf['lon_stream']
    globals()['lat_stream']=my_shelf['lat_stream']
    globals()['ii_min']=my_shelf['ii_min']
    globals()['jj_min']=my_shelf['jj_min']

    ax = fig.add_subplot(n_r,n_c,ii+1)
    m = Basemap(llcrnrlon=-120.,llcrnrlat=-90.,urcrnrlon=60.,urcrnrlat=90., projection='mill',lon_0=0)
    m.fillcontinents(color='0.8')
    m.drawparallels(np.arange(-90,90,30), labels=[1,1,0,1])
    m.drawmeridians(np.arange(0,360,60), labels=[1,1,0,1])
    bounds_max=float("{0:.01f}".format(np.nanpercentile(transport_2000_3000_mean, 99.9))) # Upper bound of plotted values to be used for colorbar, which is 98th percentile of data, with 2 decimals
    levels=np.linspace(-1*bounds_max,bounds_max,101)  
    im=m.contourf(lon_stream,lat_stream,transport_2000_3000_mean,levels,latlon=True, cmap=plt.cm.cool)
    m.scatter(lon_stream[ii_min,jj_min],lat_stream[ii_min,jj_min],0.05,c='k',latlon=True)
    cbar = m.colorbar(im,"right", size="8%", pad="40%", extend = 'both') # extend='both' will extend the colorbar in both sides (upper side and down side)
    plt.title(GCM_Names[ii]+' (yr '+str(year_start)+'-'+str(year_end)+')', fontsize=12)

plt.subplots_adjust(left=0.25, bottom=0.05, right=0.75, top=0.9, hspace=0.3, wspace=0.3) # the amount of height/width reserved for space between subplots
plt.suptitle('Mean transport 2000-3000m [Sv] (Stippling= min transport locations)', fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+'AllModels_transport_2000_3000_mean_stiple.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

#########################
#### Composits Plost ####

fig=plt.figure()
for ii in n_range:
    
    GCM = GCM_Names[ii]    
    filename_in = (dir_results + 'AllResults_'+GCM+'_500yr.out') # Directory to save processed data
    my_shelf = shelve.open(filename_in)
    globals()['year_start']=my_shelf['year_start']
    globals()['year_end']=my_shelf['year_end']    
    globals()['composite_LAB']=my_shelf['composite_LAB']

    ax = fig.add_subplot(n_r,n_c,ii+1)
    m = Basemap( projection='mill',lon_0=210)
    m.drawmapboundary(fill_color='0.9')
    m.drawmapboundary(fill_color='#000099')
    m.drawparallels(np.arange(-90,90,30), labels=[1,1,0,1])
    m.drawmeridians(np.arange(0,360,60), labels=[1,1,0,1])                      
    m.drawcoastlines(linewidth=1.0, linestyle='solid', color='k', antialiased=1, ax=None, zorder=None)
    varr=composite_LAB
    varr=np.nanmax([ np.absolute(float("{0:.01f}".format(np.nanpercentile(varr, 99)))) ,  np.absolute(float("{0:.0f}".format(np.nanpercentile(varr, 1)))) ])    
    levels=np.linspace(-varr,varr,101) 
    im=m.contourf(Lon_regrid_2D,Lat_regrid_2D,composite_LAB,levels,latlon=True, cmap=plt.cm.seismic)
    cbar = m.colorbar(im,"right", size="5%", pad="20%", extend = 'both') # extend='both' will extend the colorbar in both sides (upper side and down side)
    plt.title(GCM_Names[ii]+' (yr '+str(year_start)+'-'+str(year_end)+')', fontsize=12)

plt.subplots_adjust(left=0.15, bottom=0.05, right=0.85, top=0.9, hspace=0.3, wspace=0.1) # the amount of height/width reserved for space between subplots
plt.suptitle('Surface air temperature composites (degree C) - Labrador Sea', fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+'AllModels_surface_temperature_composites_LAB.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
levels=np.linspace(-1,1,101)
for ii in n_range:
    
    GCM = GCM_Names[ii]    
    filename_in = (dir_results + 'AllResults_'+GCM+'_500yr.out') # Directory to save processed data
    my_shelf = shelve.open(filename_in)
    globals()['year_start']=my_shelf['year_start']
    globals()['year_end']=my_shelf['year_end']    
    globals()['composite_LAB']=my_shelf['composite_LAB']

    ax = fig.add_subplot(n_r,n_c,ii+1)
    m = Basemap( projection='mill',lon_0=210)
    m.drawmapboundary(fill_color='0.9')
    m.drawmapboundary(fill_color='#000099')
    m.drawparallels(np.arange(-90,90,30), labels=[1,1,0,1])
    m.drawmeridians(np.arange(0,360,60), labels=[1,1,0,1])                      
    m.drawcoastlines(linewidth=1.0, linestyle='solid', color='k', antialiased=1, ax=None, zorder=None)
    im=m.contourf(Lon_regrid_2D,Lat_regrid_2D,composite_LAB,levels,latlon=True, cmap=plt.cm.seismic)
    plt.title(GCM_Names[ii]+' (yr '+str(year_start)+'-'+str(year_end)+')', fontsize=12)

plt.subplots_adjust(left=0.2, bottom=0.1, right=0.8, top=0.9, hspace=0.3, wspace=0.05) # the amount of height/width reserved for space between subplots
cbar_ax = fig.add_axes([0.8, 0.1, 0.015, 0.8])
fig.colorbar(im, cax=cbar_ax)
plt.suptitle('Surface air temperature composites (degree C) - Labrador Sea', fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+'AllModels_surface_temperature_composites_LAB_2.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


fig=plt.figure()
for ii in n_range:
    
    GCM = GCM_Names[ii]    
    filename_in = (dir_results + 'AllResults_'+GCM+'_500yr.out') # Directory to save processed data
    my_shelf = shelve.open(filename_in)
    globals()['year_start']=my_shelf['year_start']
    globals()['year_end']=my_shelf['year_end']    
    globals()['composite_WS']=my_shelf['composite_WS']

    ax = fig.add_subplot(n_r,n_c,ii+1)
    m = Basemap( projection='mill',lon_0=210)
    m.drawmapboundary(fill_color='0.9')
    m.drawmapboundary(fill_color='#000099')
    m.drawparallels(np.arange(-90,90,30), labels=[1,1,0,1])
    m.drawmeridians(np.arange(0,360,60), labels=[1,1,0,1])                      
    m.drawcoastlines(linewidth=1.0, linestyle='solid', color='k', antialiased=1, ax=None, zorder=None)
    varr=composite_WS
    varr=np.nanmax([ np.absolute(float("{0:.01f}".format(np.nanpercentile(varr, 99)))) ,  np.absolute(float("{0:.0f}".format(np.nanpercentile(varr, 1)))) ])    
    levels=np.linspace(-varr,varr,101) 
    im=m.contourf(Lon_regrid_2D,Lat_regrid_2D,composite_WS,levels,latlon=True, cmap=plt.cm.seismic)
    cbar = m.colorbar(im,"right", size="5%", pad="20%", extend = 'both') # extend='both' will extend the colorbar in both sides (upper side and down side)
    plt.title(GCM_Names[ii]+' (yr '+str(year_start)+'-'+str(year_end)+')', fontsize=12)

plt.subplots_adjust(left=0.15, bottom=0.05, right=0.85, top=0.9, hspace=0.3, wspace=0.1) # the amount of height/width reserved for space between subplots
plt.suptitle('Surface air temperature composites (degree C) - Weddell Sea', fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+'AllModels_surface_temperature_composites_WS.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
levels=np.linspace(-1,1,101)
for ii in n_range:
    
    GCM = GCM_Names[ii]    
    filename_in = (dir_results + 'AllResults_'+GCM+'_500yr.out') # Directory to save processed data
    my_shelf = shelve.open(filename_in)
    globals()['year_start']=my_shelf['year_start']
    globals()['year_end']=my_shelf['year_end']    
    globals()['composite_WS']=my_shelf['composite_WS']

    ax = fig.add_subplot(n_r,n_c,ii+1)
    m = Basemap( projection='mill',lon_0=210)
    m.drawmapboundary(fill_color='0.9')
    m.drawmapboundary(fill_color='#000099')
    m.drawparallels(np.arange(-90,90,30), labels=[1,1,0,1])
    m.drawmeridians(np.arange(0,360,60), labels=[1,1,0,1])                      
    m.drawcoastlines(linewidth=1.0, linestyle='solid', color='k', antialiased=1, ax=None, zorder=None)
    im=m.contourf(Lon_regrid_2D,Lat_regrid_2D,composite_WS,levels,latlon=True, cmap=plt.cm.seismic)
    plt.title(GCM_Names[ii]+' (yr '+str(year_start)+'-'+str(year_end)+')', fontsize=12)

plt.subplots_adjust(left=0.2, bottom=0.1, right=0.8, top=0.9, hspace=0.3, wspace=0.05) # the amount of height/width reserved for space between subplots
cbar_ax = fig.add_axes([0.8, 0.1, 0.015, 0.8])
fig.colorbar(im, cax=cbar_ax)
plt.suptitle('Surface air temperature composites (degree C) - Weddell Sea', fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+'AllModels_surface_temperature_composites_WS_2.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

############################################
#%% transport_lag_cor_Atlantic.py Plots ####
###LAB plots###
lag_time=40

fig=plt.figure()
for ii in n_range:
    
    GCM = GCM_Names[ii]    
    filename_in = (dir_results + 'AllResults_'+GCM+'_500yr.out') # Directory to save processed data
    my_shelf = shelve.open(filename_in)
    globals()['year_start']=my_shelf['year_start']
    globals()['year_end']=my_shelf['year_end']    
    globals()['AMOC_max']=my_shelf['AMOC_max']
    globals()['LAB_index']=my_shelf['LAB_index']
    
    LAB_index_rm = copy.deepcopy(LAB_index)
    LAB_index_rm=runningMeanFast(LAB_index_rm, 10)
    LAB_index_rm=LAB_index_rm*(-1) #multiplying by -1 as we have convection when temperature is low(opposite to the index which we calculate)

    ax = fig.add_subplot(n_r,n_c,ii+1)
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
    levels=np.linspace(-cmap_limit,cmap_limit,201)
    im=plt.contourf(x,y,data,levels,latlon=True, cmap=plt.cm.jet)
    cb = plt.colorbar(im)
    if ii+1 >= n_range[-n_c+1]: # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        plt.xlabel('Year lag', fontsize=18)
    if ii+1 == n_c+1: # Adds longitude ranges only to the last subplots that appear at the bottom of plot     
        plt.ylabel('Latitude', fontsize=18)
    plt.title(GCM_Names[ii]+' (yr '+str(year_start)+'-'+str(year_end)+')', fontsize=14)

plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, hspace=0.4, wspace=0.15) # the amount of height/width reserved for space between subplots
plt.suptitle('Labrador peak convection vs. AMOC (max streamfunction)', fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+'AllModels_lagcor_convec_AMOC_max_LAB.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
levels=np.linspace(-0.8,0.8,201)
for ii in n_range:
    
    GCM = GCM_Names[ii]    
    filename_in = (dir_results + 'AllResults_'+GCM+'_500yr.out') # Directory to save processed data
    my_shelf = shelve.open(filename_in)
    globals()['year_start']=my_shelf['year_start']
    globals()['year_end']=my_shelf['year_end']    
    globals()['AMOC_max']=my_shelf['AMOC_max']
    globals()['LAB_index']=my_shelf['LAB_index']
    
    LAB_index_rm = copy.deepcopy(LAB_index)
    LAB_index_rm=runningMeanFast(LAB_index_rm, 10)
    LAB_index_rm=LAB_index_rm*(-1) #multiplying by -1 as we have convection when temperature is low(opposite to the index which we calculate)

    ax = fig.add_subplot(n_r,n_c,ii+1)
    data=[]
    for i in range(len(AMOC_max[0][20:161])):
        stream=runningMeanFast(AMOC_max[:,i+lag_time], 10)
        r=lag_cor_data(LAB_index_rm,stream,lag_time)
        data.append(r)
    data=np.asarray(data)
    x=np.linspace(-lag_time,lag_time+1, 2*lag_time)
    y=np.linspace(-70,70+1, 2*70+1)
    x,y=np.meshgrid(x,y)
    im=plt.contourf(x,y,data,levels,latlon=True, cmap=plt.cm.jet)
    if ii+1 >= n_range[-n_c+1]: # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        plt.xlabel('Year lag', fontsize=18)
    if ii+1 == n_c+1: # Adds longitude ranges only to the last subplots that appear at the bottom of plot     
        plt.ylabel('Latitude', fontsize=18)
    plt.title(GCM_Names[ii]+' (yr '+str(year_start)+'-'+str(year_end)+')', fontsize=14)

plt.subplots_adjust(left=0.15, bottom=0.1, right=0.85, top=0.9, hspace=0.4, wspace=0.15) # the amount of height/width reserved for space between subplots
cbar_ax = fig.add_axes([0.87, 0.1, 0.015, 0.8])
fig.colorbar(im, cax=cbar_ax)
plt.suptitle('Labrador peak convection vs. AMOC (max streamfunction)', fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+'AllModels_lagcor_convec_AMOC_max_LAB_2.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
for ii in n_range:
    
    GCM = GCM_Names[ii]    
    filename_in = (dir_results + 'AllResults_'+GCM+'_500yr.out') # Directory to save processed data
    my_shelf = shelve.open(filename_in)
    globals()['year_start']=my_shelf['year_start']
    globals()['year_end']=my_shelf['year_end']    
    globals()['SMOC_min']=my_shelf['SMOC_min']
    globals()['LAB_index']=my_shelf['LAB_index']
    
    LAB_index_rm = copy.deepcopy(LAB_index)
    LAB_index_rm=runningMeanFast(LAB_index_rm, 10)
    LAB_index_rm=LAB_index_rm*(-1) #multiplying by -1 as we have convection when temperature is low(opposite to the index which we calculate)

    ax = fig.add_subplot(n_r,n_c,ii+1)
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
    levels=np.linspace(-cmap_limit,cmap_limit,201)
    im=plt.contourf(x,y,data,levels,latlon=True, cmap=plt.cm.jet)
    cb = plt.colorbar(im)
    if ii+1 >= n_range[-n_c+1]: # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        plt.xlabel('Year lag', fontsize=18)
    if ii+1 == n_c+1: # Adds longitude ranges only to the last subplots that appear at the bottom of plot     
        plt.ylabel('Latitude', fontsize=18)
    plt.title(GCM_Names[ii]+' (yr '+str(year_start)+'-'+str(year_end)+')', fontsize=14)

plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, hspace=0.4, wspace=0.15) # the amount of height/width reserved for space between subplots
plt.suptitle('Labrador peak convection vs. SMOC (min streamfunction)', fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+'AllModels_lagcor_convec_SMOC_min_LAB.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
levels=np.linspace(-0.8,0.8,201)
for ii in n_range:
    
    GCM = GCM_Names[ii]    
    filename_in = (dir_results + 'AllResults_'+GCM+'_500yr.out') # Directory to save processed data
    my_shelf = shelve.open(filename_in)
    globals()['year_start']=my_shelf['year_start']
    globals()['year_end']=my_shelf['year_end']    
    globals()['SMOC_min']=my_shelf['SMOC_min']
    globals()['LAB_index']=my_shelf['LAB_index']
    
    LAB_index_rm = copy.deepcopy(LAB_index)
    LAB_index_rm=runningMeanFast(LAB_index_rm, 10)
    LAB_index_rm=LAB_index_rm*(-1) #multiplying by -1 as we have convection when temperature is low(opposite to the index which we calculate)

    ax = fig.add_subplot(n_r,n_c,ii+1)
    data=[]
    for i in range(len(SMOC_min[0][20:161])):
        stream=runningMeanFast(SMOC_min[:,i+lag_time], 10)
        r=lag_cor_data(LAB_index_rm,stream,lag_time)
        data.append(r)
    data=np.asarray(data)
    x=np.linspace(-lag_time,lag_time+1, 2*lag_time)
    y=np.linspace(-70,70+1, 2*70+1)
    x,y=np.meshgrid(x,y)
    im=plt.contourf(x,y,data,levels,latlon=True, cmap=plt.cm.jet)
    if ii+1 >= n_range[-n_c+1]: # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        plt.xlabel('Year lag', fontsize=18)
    if ii+1 == n_c+1: # Adds longitude ranges only to the last subplots that appear at the bottom of plot     
        plt.ylabel('Latitude', fontsize=18)
    plt.title(GCM_Names[ii]+' (yr '+str(year_start)+'-'+str(year_end)+')', fontsize=14)

plt.subplots_adjust(left=0.15, bottom=0.1, right=0.85, top=0.9, hspace=0.4, wspace=0.15) # the amount of height/width reserved for space between subplots
cbar_ax = fig.add_axes([0.87, 0.1, 0.015, 0.8])
fig.colorbar(im, cax=cbar_ax)
plt.suptitle('Labrador peak convection vs. SMOC (min streamfunction)', fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+'AllModels_lagcor_convec_SMOC_min_LAB_2.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

#### WS plots ####
fig=plt.figure()
for ii in n_range:
    
    GCM = GCM_Names[ii]    
    filename_in = (dir_results + 'AllResults_'+GCM+'_500yr.out') # Directory to save processed data
    my_shelf = shelve.open(filename_in)
    globals()['year_start']=my_shelf['year_start']
    globals()['year_end']=my_shelf['year_end']    
    globals()['AMOC_max']=my_shelf['AMOC_max']
    globals()['WS_index']=my_shelf['WS_index']
    
    WS_index_rm = copy.deepcopy(WS_index)
    WS_index_rm=runningMeanFast(WS_index_rm, 10)
    WS_index_rm=WS_index_rm*(-1) #multiplying by -1 as we have convection when temperature is low(opposite to the index which we calculate)

    ax = fig.add_subplot(n_r,n_c,ii+1)
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
    levels=np.linspace(-cmap_limit,cmap_limit,201)
    im=plt.contourf(x,y,data,levels,latlon=True, cmap=plt.cm.jet)
    cb = plt.colorbar(im)
    if ii+1 >= n_range[-n_c+1]: # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        plt.xlabel('Year lag', fontsize=18)
    if ii+1 == n_c+1: # Adds longitude ranges only to the last subplots that appear at the bottom of plot     
        plt.ylabel('Latitude', fontsize=18)
    plt.title(GCM_Names[ii]+' (yr '+str(year_start)+'-'+str(year_end)+')', fontsize=14)

plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, hspace=0.4, wspace=0.15) # the amount of height/width reserved for space between subplots
plt.suptitle('Weddell peak convection vs. AMOC (max streamfunction)', fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+'AllModels_lagcor_convec_AMOC_max_WS.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
levels=np.linspace(-0.8,0.8,201)
for ii in n_range:
    
    GCM = GCM_Names[ii]    
    filename_in = (dir_results + 'AllResults_'+GCM+'_500yr.out') # Directory to save processed data
    my_shelf = shelve.open(filename_in)
    globals()['year_start']=my_shelf['year_start']
    globals()['year_end']=my_shelf['year_end']    
    globals()['AMOC_max']=my_shelf['AMOC_max']
    globals()['WS_index']=my_shelf['WS_index']
    
    WS_index_rm = copy.deepcopy(WS_index)
    WS_index_rm=runningMeanFast(WS_index_rm, 10)
    WS_index_rm=WS_index_rm*(-1) #multiplying by -1 as we have convection when temperature is low(opposite to the index which we calculate)

    ax = fig.add_subplot(n_r,n_c,ii+1)
    data=[]
    for i in range(len(AMOC_max[0][20:161])):
        stream=runningMeanFast(AMOC_max[:,i+lag_time], 10)
        r=lag_cor_data(WS_index_rm,stream,lag_time)
        data.append(r)
    data=np.asarray(data)
    x=np.linspace(-lag_time,lag_time+1, 2*lag_time)
    y=np.linspace(-70,70+1, 2*70+1)
    x,y=np.meshgrid(x,y)
    im=plt.contourf(x,y,data,levels,latlon=True, cmap=plt.cm.jet)
    if ii+1 >= n_range[-n_c+1]: # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        plt.xlabel('Year lag', fontsize=18)
    if ii+1 == n_c+1: # Adds longitude ranges only to the last subplots that appear at the bottom of plot     
        plt.ylabel('Latitude', fontsize=18)
    plt.title(GCM_Names[ii]+' (yr '+str(year_start)+'-'+str(year_end)+')', fontsize=14)

plt.subplots_adjust(left=0.15, bottom=0.1, right=0.85, top=0.9, hspace=0.4, wspace=0.15) # the amount of height/width reserved for space between subplots
cbar_ax = fig.add_axes([0.87, 0.1, 0.015, 0.8])
fig.colorbar(im, cax=cbar_ax)
plt.suptitle('Weddell peak convection vs. AMOC (max streamfunction)', fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+'AllModels_lagcor_convec_AMOC_max_WS_2.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
for ii in n_range:
    
    GCM = GCM_Names[ii]    
    filename_in = (dir_results + 'AllResults_'+GCM+'_500yr.out') # Directory to save processed data
    my_shelf = shelve.open(filename_in)
    globals()['year_start']=my_shelf['year_start']
    globals()['year_end']=my_shelf['year_end']    
    globals()['SMOC_min']=my_shelf['SMOC_min']
    globals()['WS_index']=my_shelf['WS_index']
    
    WS_index_rm = copy.deepcopy(WS_index)
    WS_index_rm=runningMeanFast(WS_index_rm, 10)
    WS_index_rm=WS_index_rm*(-1) #multiplying by -1 as we have convection when temperature is low(opposite to the index which we calculate)

    ax = fig.add_subplot(n_r,n_c,ii+1)
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
    levels=np.linspace(-cmap_limit,cmap_limit,201)
    im=plt.contourf(x,y,data,levels,latlon=True, cmap=plt.cm.jet)
    cb = plt.colorbar(im)
    if ii+1 >= n_range[-n_c+1]: # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        plt.xlabel('Year lag', fontsize=18)
    if ii+1 == n_c+1: # Adds longitude ranges only to the last subplots that appear at the bottom of plot     
        plt.ylabel('Latitude', fontsize=18)
    plt.title(GCM_Names[ii]+' (yr '+str(year_start)+'-'+str(year_end)+')', fontsize=14)

plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, hspace=0.4, wspace=0.15) # the amount of height/width reserved for space between subplots
plt.suptitle('Weddell peak convection vs. SMOC (min streamfunction)', fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+'AllModels_lagcor_convec_SMOC_min_WS.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
levels=np.linspace(-0.8,0.8,201)
for ii in n_range:
    
    GCM = GCM_Names[ii]    
    filename_in = (dir_results + 'AllResults_'+GCM+'_500yr.out') # Directory to save processed data
    my_shelf = shelve.open(filename_in)
    globals()['year_start']=my_shelf['year_start']
    globals()['year_end']=my_shelf['year_end']    
    globals()['SMOC_min']=my_shelf['SMOC_min']
    globals()['WS_index']=my_shelf['WS_index']
    
    WS_index_rm = copy.deepcopy(WS_index)
    WS_index_rm=runningMeanFast(WS_index_rm, 10)
    WS_index_rm=WS_index_rm*(-1) #multiplying by -1 as we have convection when temperature is low(opposite to the index which we calculate)

    ax = fig.add_subplot(n_r,n_c,ii+1)
    data=[]
    for i in range(len(SMOC_min[0][20:161])):
        stream=runningMeanFast(SMOC_min[:,i+lag_time], 10)
        r=lag_cor_data(WS_index_rm,stream,lag_time)
        data.append(r)
    data=np.asarray(data)
    x=np.linspace(-lag_time,lag_time+1, 2*lag_time)
    y=np.linspace(-70,70+1, 2*70+1)
    x,y=np.meshgrid(x,y)
    im=plt.contourf(x,y,data,levels,latlon=True, cmap=plt.cm.jet)
    if ii+1 >= n_range[-n_c+1]: # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        plt.xlabel('Year lag', fontsize=18)
    if ii+1 == n_c+1: # Adds longitude ranges only to the last subplots that appear at the bottom of plot     
        plt.ylabel('Latitude', fontsize=18)
    plt.title(GCM_Names[ii]+' (yr '+str(year_start)+'-'+str(year_end)+')', fontsize=14)

plt.subplots_adjust(left=0.15, bottom=0.1, right=0.85, top=0.9, hspace=0.4, wspace=0.15) # the amount of height/width reserved for space between subplots
cbar_ax = fig.add_axes([0.87, 0.1, 0.015, 0.8])
fig.colorbar(im, cax=cbar_ax)
plt.suptitle('Weddell peak convection vs. SMOC (min streamfunction)', fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+'AllModels_lagcor_convec_SMOC_min_WS_2.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

####################################################################
### AMOC/SMOC Lagged Corrolation Plots - Max/Min Stream function ###

fig=plt.figure()
for ii in n_range:
    
    GCM = GCM_Names[ii]    
    filename_in = (dir_results + 'AllResults_'+GCM+'_500yr.out') # Directory to save processed data
    my_shelf = shelve.open(filename_in)
    globals()['year_start']=my_shelf['year_start']
    globals()['year_end']=my_shelf['year_end']    
    globals()['AMOC_max']=my_shelf['AMOC_max']
    globals()['LAB_index']=my_shelf['LAB_index']
    globals()['WS_index']=my_shelf['WS_index'] 
    globals()['Zonal_winds']=my_shelf['Zonal_winds']
    globals()['NAO_index']=my_shelf['NAO_index']    
    globals()['ENSO_index']=my_shelf['ENSO_index']    

    WS_index_rm = copy.deepcopy(WS_index)
    WS_index_rm=runningMeanFast(WS_index_rm, 10)
    WS_index_rm=WS_index_rm*(-1) #multiplying by -1 as we have convection when temperature is low(opposite to the index which we calculate)    
    LAB_index_rm = copy.deepcopy(LAB_index)
    LAB_index_rm=runningMeanFast(LAB_index_rm, 10)
    LAB_index_rm=LAB_index_rm*(-1) #multiplying by -1 as we have convection when temperature is low(opposite to the index which we calculate)
    North_Westerlies=Zonal_winds[:,135:145]
    North_Trades=Zonal_winds[:,100:110]
    NAO_m=runningMeanFast(NAO_index, 10)
    ENSO_m=runningMeanFast(ENSO_index, 10)
    AMOC_max_50S= AMOC_max[:,40] # AMOC at 50S # Max of streamfunction method
    AMOC_max_50S_m=runningMeanFast(AMOC_max_50S, 10)
    AMOC_max_30S= AMOC_max[:,60] # AMOC at 30S # Max of streamfunction method
    AMOC_max_30S_m=runningMeanFast(AMOC_max_30S, 10)
    AMOC_max_50N= AMOC_max[:,140] # AMOC at 50N # Max of streamfunction method
    AMOC_max_50N_m=runningMeanFast(AMOC_max_50N, 10)
    AMOC_max_30N= AMOC_max[:,120] # AMOC at 30N # Max of streamfunction method
    AMOC_max_30N_m=runningMeanFast(AMOC_max_30N, 10)

    ax = fig.add_subplot(n_r,n_c,ii+1)
    lag_cor(LAB_index_rm[:-9],WS_index_rm[:-9],40,'y','WS index')
    lag_cor(LAB_index_rm[:-9],ENSO_m[:-9],40,'m','ENSO index')
    lag_cor(LAB_index_rm[:-9],NAO_m[:-9],40,'c','NAO index')
    lag_cor(LAB_index_rm[:-9], np.nanmean(North_Westerlies, axis=1)[:-9],40,'slategrey','North Westerlies')
    lag_cor(LAB_index_rm[:-9], np.nanmean(North_Trades, axis=1)[:-9],40,'darkslategrey','North Trades')
    lag_cor(LAB_index_rm[:-9],AMOC_max_50N_m[:-9],40,'b','AMOC max 50N')
    lag_cor(LAB_index_rm[:-9],AMOC_max_30N_m[:-9],40,'g','AMOC max 30N')
    lag_cor(LAB_index_rm[:-9],AMOC_max_50S_m[:-9],40,'r','AMOC max 50S')
    lag_cor(LAB_index_rm[:-9],AMOC_max_30S_m[:-9],40,'darkorange','AMOC max 30S')    
    #plt.legend(shadow=True, prop={'size': 15})
    if ii+1 >= n_range[-n_c+1]: # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        plt.xlabel('Years lag', fontsize=18)
    else:
        plt.xlabel('')
    if ii+1 == n_c+1: # Adds longitude ranges only to the last subplots that appear at the bottom of plot     
        plt.ylabel('Correlation coefficient', fontsize=18)  
    else:
        plt.ylabel('')  
    plt.title(GCM_Names[ii]+' (yr '+str(year_start)+'-'+str(year_end)+')', fontsize=14)

plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, hspace=0.4, wspace=0.15) # the amount of height/width reserved for space between subplots
#cbar_ax = fig.add_axes([0.87, 0.1, 0.015, 0.8])
fig.legend(shadow=True, loc='lower right', bbox_to_anchor=(0.91, 0.05), ncol=2 , prop={'size': 13})
#fig.legend(shadow=True, bbox_to_anchor=(0., 1.02, 1., .102), loc='lower right', ncol=2, mode="expand", borderaxespad=0.)
plt.suptitle('Lagged correlations - Labrador peak convection ', fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+'AllModels_lagcor_LABconvection_AMOC_max_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
for ii in n_range:
    
    GCM = GCM_Names[ii]    
    filename_in = (dir_results + 'AllResults_'+GCM+'_500yr.out') # Directory to save processed data
    my_shelf = shelve.open(filename_in)
    globals()['year_start']=my_shelf['year_start']
    globals()['year_end']=my_shelf['year_end']    
    globals()['SMOC_min']=my_shelf['SMOC_min']
    globals()['LAB_index']=my_shelf['LAB_index']
    globals()['WS_index']=my_shelf['WS_index'] 
    globals()['Zonal_winds']=my_shelf['Zonal_winds']
    globals()['NAO_index']=my_shelf['NAO_index']    
    globals()['ENSO_index']=my_shelf['ENSO_index']    

    WS_index_rm = copy.deepcopy(WS_index)
    WS_index_rm=runningMeanFast(WS_index_rm, 10)
    WS_index_rm=WS_index_rm*(-1) #multiplying by -1 as we have convection when temperature is low(opposite to the index which we calculate)    
    LAB_index_rm = copy.deepcopy(LAB_index)
    LAB_index_rm=runningMeanFast(LAB_index_rm, 10)
    LAB_index_rm=LAB_index_rm*(-1) #multiplying by -1 as we have convection when temperature is low(opposite to the index which we calculate)
    South_Westerlies=Zonal_winds[:,35:45]
    South_Trades=Zonal_winds[:,70:80]
    NAO_m=runningMeanFast(NAO_index, 10)
    ENSO_m=runningMeanFast(ENSO_index, 10)
    SMOC_min_50S= SMOC_min[:,40] # SMOC at 50S # Max of streamfunction method
    SMOC_min_50S_m=runningMeanFast(SMOC_min_50S, 10)
    SMOC_min_30S= SMOC_min[:,60] # SMOC at 30S # Max of streamfunction method
    SMOC_min_30S_m=runningMeanFast(SMOC_min_30S, 10)
    SMOC_min_50N= SMOC_min[:,140] # SMOC at 50N # Max of streamfunction method
    SMOC_min_50N_m=runningMeanFast(SMOC_min_50N, 10)
    SMOC_min_30N= SMOC_min[:,120] # SMOC at 30N # Max of streamfunction method
    SMOC_min_30N_m=runningMeanFast(SMOC_min_30N, 10)

    ax = fig.add_subplot(n_r,n_c,ii+1)
    lag_cor(LAB_index_rm[:-9],WS_index_rm[:-9],40,'y','WS index')
    lag_cor(LAB_index_rm[:-9],ENSO_m[:-9],40,'m','ENSO index')
    lag_cor(LAB_index_rm[:-9],NAO_m[:-9],40,'c','NAO index')
    lag_cor(LAB_index_rm[:-9], np.nanmean(South_Westerlies, axis=1)[:-9],40,'slategrey','South Westerlies')
    lag_cor(LAB_index_rm[:-9], np.nanmean(South_Trades, axis=1)[:-9],40,'darkslategrey','South Trades')
    lag_cor(LAB_index_rm[:-9],SMOC_min_50N_m[:-9],40,'b','SMOC min 50N')
    lag_cor(LAB_index_rm[:-9],SMOC_min_30N_m[:-9],40,'g','SMOC min 30N')
    lag_cor(LAB_index_rm[:-9],SMOC_min_50S_m[:-9],40,'r','SMOC min 50S')
    lag_cor(LAB_index_rm[:-9],SMOC_min_30S_m[:-9],40,'darkorange','SMOC min 30S')    
    #plt.legend(shadow=True, prop={'size': 15})
    if ii+1 >= n_range[-n_c+1]: # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        plt.xlabel('Years lag', fontsize=18)
    else:
        plt.xlabel('')
    if ii+1 == n_c+1: # Adds longitude ranges only to the last subplots that appear at the bottom of plot     
        plt.ylabel('Correlation coefficient', fontsize=18)  
    else:
        plt.ylabel('')  
    plt.title(GCM_Names[ii]+' (yr '+str(year_start)+'-'+str(year_end)+')', fontsize=14)

plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, hspace=0.4, wspace=0.15) # the amount of height/width reserved for space between subplots
#cbar_ax = fig.add_axes([0.87, 0.1, 0.015, 0.8])
fig.legend(shadow=True, loc='lower right', bbox_to_anchor=(0.91, 0.05), ncol=2 , prop={'size': 13})
#fig.legend(shadow=True, bbox_to_anchor=(0., 1.02, 1., .102), loc='lower right', ncol=2, mode="expand", borderaxespad=0.)
plt.suptitle('Lagged correlations - Labrador peak convection ', fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+'AllModels_lagcor_LABconvection_SMOC_min_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


fig=plt.figure()
for ii in n_range:
    
    GCM = GCM_Names[ii]    
    filename_in = (dir_results + 'AllResults_'+GCM+'_500yr.out') # Directory to save processed data
    my_shelf = shelve.open(filename_in)
    globals()['year_start']=my_shelf['year_start']
    globals()['year_end']=my_shelf['year_end']    
    globals()['AMOC_max']=my_shelf['AMOC_max']
    globals()['LAB_index']=my_shelf['LAB_index']
    globals()['WS_index']=my_shelf['WS_index'] 
    globals()['Zonal_winds']=my_shelf['Zonal_winds']
    globals()['NAO_index']=my_shelf['NAO_index']    
    globals()['ENSO_index']=my_shelf['ENSO_index']    

    WS_index_rm = copy.deepcopy(WS_index)
    WS_index_rm=runningMeanFast(WS_index_rm, 10)
    WS_index_rm=WS_index_rm*(-1) #multiplying by -1 as we have convection when temperature is low(opposite to the index which we calculate)    
    LAB_index_rm = copy.deepcopy(LAB_index)
    LAB_index_rm=runningMeanFast(LAB_index_rm, 10)
    LAB_index_rm=LAB_index_rm*(-1) #multiplying by -1 as we have convection when temperature is low(opposite to the index which we calculate)
    North_Westerlies=Zonal_winds[:,135:145]
    North_Trades=Zonal_winds[:,100:110]
    NAO_m=runningMeanFast(NAO_index, 10)
    ENSO_m=runningMeanFast(ENSO_index, 10)
    AMOC_max_50S= AMOC_max[:,40] # AMOC at 50S # Max of streamfunction method
    AMOC_max_50S_m=runningMeanFast(AMOC_max_50S, 10)
    AMOC_max_30S= AMOC_max[:,60] # AMOC at 30S # Max of streamfunction method
    AMOC_max_30S_m=runningMeanFast(AMOC_max_30S, 10)
    AMOC_max_50N= AMOC_max[:,140] # AMOC at 50N # Max of streamfunction method
    AMOC_max_50N_m=runningMeanFast(AMOC_max_50N, 10)
    AMOC_max_30N= AMOC_max[:,120] # AMOC at 30N # Max of streamfunction method
    AMOC_max_30N_m=runningMeanFast(AMOC_max_30N, 10)

    ax = fig.add_subplot(n_r,n_c,ii+1)
    lag_cor(WS_index_rm[:-9],LAB_index_rm[:-9],40,'y','LAB index')
    lag_cor(WS_index_rm[:-9],ENSO_m[:-9],40,'m','ENSO index')
    lag_cor(WS_index_rm[:-9],NAO_m[:-9],40,'c','NAO index')
    lag_cor(WS_index_rm[:-9], np.nanmean(North_Westerlies, axis=1)[:-9],40,'slategrey','North Westerlies')
    lag_cor(WS_index_rm[:-9], np.nanmean(North_Trades, axis=1)[:-9],40,'darkslategrey','North Trades')
    lag_cor(WS_index_rm[:-9],AMOC_max_50N_m[:-9],40,'b','AMOC max 50N')
    lag_cor(WS_index_rm[:-9],AMOC_max_30N_m[:-9],40,'g','AMOC max 30N')
    lag_cor(WS_index_rm[:-9],AMOC_max_50S_m[:-9],40,'r','AMOC max 50S')
    lag_cor(WS_index_rm[:-9],AMOC_max_30S_m[:-9],40,'darkorange','AMOC max 30S')    
    #plt.legend(shadow=True, prop={'size': 15})
    if ii+1 >= n_range[-n_c+1]: # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        plt.xlabel('Years lag', fontsize=18)
    else:
        plt.xlabel('')
    if ii+1 == n_c+1: # Adds longitude ranges only to the last subplots that appear at the bottom of plot     
        plt.ylabel('Correlation coefficient', fontsize=18)  
    else:
        plt.ylabel('')  
    plt.title(GCM_Names[ii]+' (yr '+str(year_start)+'-'+str(year_end)+')', fontsize=14)

plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, hspace=0.4, wspace=0.15) # the amount of height/width reserved for space between subplots
#cbar_ax = fig.add_axes([0.87, 0.1, 0.015, 0.8])
fig.legend(shadow=True, loc='lower right', bbox_to_anchor=(0.91, 0.05), ncol=2 , prop={'size': 13})
#fig.legend(shadow=True, bbox_to_anchor=(0., 1.02, 1., .102), loc='lower right', ncol=2, mode="expand", borderaxespad=0.)
plt.suptitle('Lagged correlations - Weddell peak convection ', fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+'AllModels_lagcor_WSconvection_AMOC_max_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
for ii in n_range:
    
    GCM = GCM_Names[ii]    
    filename_in = (dir_results + 'AllResults_'+GCM+'_500yr.out') # Directory to save processed data
    my_shelf = shelve.open(filename_in)
    globals()['year_start']=my_shelf['year_start']
    globals()['year_end']=my_shelf['year_end']    
    globals()['SMOC_min']=my_shelf['SMOC_min']
    globals()['LAB_index']=my_shelf['LAB_index']
    globals()['WS_index']=my_shelf['WS_index'] 
    globals()['Zonal_winds']=my_shelf['Zonal_winds']
    globals()['NAO_index']=my_shelf['NAO_index']    
    globals()['ENSO_index']=my_shelf['ENSO_index']    

    WS_index_rm = copy.deepcopy(WS_index)
    WS_index_rm=runningMeanFast(WS_index_rm, 10)
    WS_index_rm=WS_index_rm*(-1) #multiplying by -1 as we have convection when temperature is low(opposite to the index which we calculate)    
    LAB_index_rm = copy.deepcopy(LAB_index)
    LAB_index_rm=runningMeanFast(LAB_index_rm, 10)
    LAB_index_rm=LAB_index_rm*(-1) #multiplying by -1 as we have convection when temperature is low(opposite to the index which we calculate)
    South_Westerlies=Zonal_winds[:,35:45]
    South_Trades=Zonal_winds[:,70:80]
    NAO_m=runningMeanFast(NAO_index, 10)
    ENSO_m=runningMeanFast(ENSO_index, 10)
    SMOC_min_50S= SMOC_min[:,40] # SMOC at 50S # Max of streamfunction method
    SMOC_min_50S_m=runningMeanFast(SMOC_min_50S, 10)
    SMOC_min_30S= SMOC_min[:,60] # SMOC at 30S # Max of streamfunction method
    SMOC_min_30S_m=runningMeanFast(SMOC_min_30S, 10)
    SMOC_min_50N= SMOC_min[:,140] # SMOC at 50N # Max of streamfunction method
    SMOC_min_50N_m=runningMeanFast(SMOC_min_50N, 10)
    SMOC_min_30N= SMOC_min[:,120] # SMOC at 30N # Max of streamfunction method
    SMOC_min_30N_m=runningMeanFast(SMOC_min_30N, 10)

    ax = fig.add_subplot(n_r,n_c,ii+1)
    lag_cor(WS_index_rm[:-9],LAB_index_rm[:-9],40,'y','LAB index')
    lag_cor(WS_index_rm[:-9],ENSO_m[:-9],40,'m','ENSO index')
    lag_cor(WS_index_rm[:-9],NAO_m[:-9],40,'c','NAO index')
    lag_cor(WS_index_rm[:-9], np.nanmean(South_Westerlies, axis=1)[:-9],40,'slategrey','South Westerlies')
    lag_cor(WS_index_rm[:-9], np.nanmean(South_Trades, axis=1)[:-9],40,'darkslategrey','South Trades')
    lag_cor(WS_index_rm[:-9],SMOC_min_50N_m[:-9],40,'b','SMOC min 50N')
    lag_cor(WS_index_rm[:-9],SMOC_min_30N_m[:-9],40,'g','SMOC min 30N')
    lag_cor(WS_index_rm[:-9],SMOC_min_50S_m[:-9],40,'r','SMOC min 50S')
    lag_cor(WS_index_rm[:-9],SMOC_min_30S_m[:-9],40,'darkorange','SMOC min 30S')    
    #plt.legend(shadow=True, prop={'size': 15})
    if ii+1 >= n_range[-n_c+1]: # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        plt.xlabel('Years lag', fontsize=18)
    else:
        plt.xlabel('')
    if ii+1 == n_c+1: # Adds longitude ranges only to the last subplots that appear at the bottom of plot     
        plt.ylabel('Correlation coefficient', fontsize=18)  
    else:
        plt.ylabel('')  
    plt.title(GCM_Names[ii]+' (yr '+str(year_start)+'-'+str(year_end)+')', fontsize=14)

plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, hspace=0.4, wspace=0.15) # the amount of height/width reserved for space between subplots
#cbar_ax = fig.add_axes([0.87, 0.1, 0.015, 0.8])
fig.legend(shadow=True, loc='lower right', bbox_to_anchor=(0.91, 0.05), ncol=2 , prop={'size': 13})
#fig.legend(shadow=True, bbox_to_anchor=(0., 1.02, 1., .102), loc='lower right', ncol=2, mode="expand", borderaxespad=0.)
plt.suptitle('Lagged correlations - Weddell peak convection ', fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+'AllModels_lagcor_WSconvection_SMOC_min_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


































########################################
##%% Streamfunction Plots ####
#
#Plot_title=('Ocean depth map (Atlantic) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM))
#Plot_save_dir=(dir_figs+str(GCM)+'_depth_map_atl.png')
#bounds_max=np.int(np.nanpercentile(latlon_depths, 99.99))
#bounds = np.arange(0, bounds_max, bounds_max/40)
##func_plot(latlon_depths, lat_stream, lon_stream, bounds_max, '-', '-', 'mill', 0)
#func_plot_bounds_save(latlon_depths, lat_stream, lon_stream, bounds, '(m)', Plot_title, 'mill', 0, Plot_save_dir)
#
#Plot_title=('Mean northward transport upper 1000m [Sv] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM))
#Plot_save_dir=(dir_figs+str(GCM)+'_transport_0_1000_mean.png')
##bounds_max=np.nanpercentile(transport_0_1000_mean, 99.99)
#bounds_max=float("{0:.02f}".format(np.nanpercentile(transport_0_1000_mean, 98))) # Upper bound of plotted values to be used for colorbar, which is 98th percentile of data, with 2 decimals
#bounds = np.arange(-1*bounds_max, bounds_max, bounds_max/20)
##bounds = np.arange(-0.1, 0.1, 0.2/20)
#func_plot_bounds_save(transport_0_1000_mean, lat_stream, lon_stream, bounds, '(Sv)', Plot_title, 'mill', 0, Plot_save_dir)  
#
#Plot_title=('Mean northward transport at 2000m-3000m [Sv] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM))
#Plot_save_dir=(dir_figs+str(GCM)+'_transport_2000_3000_mean.png')
##bounds_max=np.nanpercentile(transport_0_1000_mean, 99.99)
#bounds_max=float("{0:.02f}".format(np.nanpercentile(transport_2000_3000_mean, 98))) # Upper bound of plotted values to be used for colorbar, which is 98th percentile of data, with 2 decimals
#bounds = np.arange(-1*bounds_max, bounds_max, bounds_max/20)
##bounds = np.arange(-0.1, 0.1, 0.2/20)
#func_plot_bounds_save(transport_2000_3000_mean, lat_stream, lon_stream, bounds, '(Sv)', Plot_title, 'mill', 0, Plot_save_dir)     
#
#############################################
##%% transport_lag_cor_Atlantic.py Plots ####
#
#AMOC_transport_all = copy.deepcopy(transport_0_1000)
#SMOC_transport_all= copy.deepcopy(transport_2000_3000)
#SMOC_transport_all=SMOC_transport_all*(-1)
##As SMOC is transport averaged over 2000-3000m it is southward flow, thus the sign is negative (for southward flow) so we multiply by -1 to calculate further statistics
#
#LAB_index_rm = copy.deepcopy(LAB_index)
#WS_index_rm = copy.deepcopy(WS_index)
#
#LAB_index_rm=runningMeanFast(LAB_index_rm, 10)
##multiplying by -1 as we have convection when temperature is low(opposite to the index which we calculate)
#LAB_index_rm=LAB_index_rm*(-1)
#WS_index_rm=runningMeanFast(WS_index_rm, 10)
#WS_index_rm=WS_index_rm*(-1)
#
##summing transport over all longtitudes
#AMOC_transport=np.nansum(AMOC_transport_all,axis=2)
#SMOC_transport=np.nansum(SMOC_transport_all,axis=2)
#
####LAB plots###
#lag_time=40
#
#fig=plt.figure()
#data=[]
#for i in range(len(AMOC_transport[0][20:161])):
#    stream=runningMeanFast(AMOC_transport[:,i+lag_time], 10)
#    r=lag_cor_data(LAB_index_rm,stream,lag_time)
#    data.append(r)
#data=np.asarray(data)
#x=np.linspace(-lag_time,lag_time+1, 2*lag_time)
#y=np.linspace(-70,70+1, 2*70+1)
#x,y=np.meshgrid(x,y)
#cmap_limit=np.nanmax(np.abs(data))
#levels=np.linspace(-cmap_limit,cmap_limit,200)
#im=plt.contourf(x,y,data,levels,latlon=True, cmap=plt.cm.jet)
#cb = plt.colorbar(im)
#plt.xlabel('Year lag', fontsize=18)
#plt.ylabel('Latitude', fontsize=18)
#plt.title('Labrador peak convection vs. transport in upper 1000m - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
#mng = plt.get_current_fig_manager()
#mng.window.showMaximized() # Maximizes the plot window to save figures in full
#fig.savefig(dir_figs+str(GCM)+'_lagcor_convec_transport_0_1000_LAB.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
#
#
#fig=plt.figure()
#data=[]
#for i in range(len(SMOC_transport[0][20:161])):
#    stream=runningMeanFast(SMOC_transport[:,i+lag_time], 10)
#    r=lag_cor_data(LAB_index_rm,stream,lag_time)
#    data.append(r)
#data=np.asarray(data)
#x=np.linspace(-lag_time,lag_time+1, 2*lag_time)
#y=np.linspace(-70,70+1, 2*70+1)
#x,y=np.meshgrid(x,y)
#cmap_limit=np.nanmax(np.abs(data))
#levels=np.linspace(-cmap_limit,cmap_limit,200)
#im=plt.contourf(x,y,data,levels,latlon=True, cmap=plt.cm.jet)
#cb = plt.colorbar(im)
#plt.xlabel('Year lag', fontsize=18)
#plt.ylabel('Latitude', fontsize=18)
#plt.title('Labrador peak convection vs. transport in 2000m-3000m - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
#mng = plt.get_current_fig_manager()
#mng.window.showMaximized() # Maximizes the plot window to save figures in full
#fig.savefig(dir_figs+str(GCM)+'_lagcor_convec_transport_2000_3000_LAB.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
#
##### WS plots ####
#lag_time=40
#
#fig=plt.figure()
#data=[]
#for i in range(len(AMOC_transport[0][20:161])):
#    stream=runningMeanFast(AMOC_transport[:,i+lag_time], 10)
#    r=lag_cor_data(WS_index_rm,stream,lag_time)
#    data.append(r)
#data=np.asarray(data)
#x=np.linspace(-lag_time,lag_time+1, 2*lag_time)
#y=np.linspace(-70,70+1, 2*70+1)
#x,y=np.meshgrid(x,y)
#cmap_limit=np.nanmax(np.abs(data))
#levels=np.linspace(-cmap_limit,cmap_limit,200)
#im=plt.contourf(x,y,data,levels,latlon=True, cmap=plt.cm.jet)
#cb = plt.colorbar(im)
#plt.xlabel('Year lag', fontsize=18)
#plt.ylabel('Latitude', fontsize=18)
#plt.title('WS peak convection vs. transport in upper 1000m - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
#mng = plt.get_current_fig_manager()
#mng.window.showMaximized() # Maximizes the plot window to save figures in full
#fig.savefig(dir_figs+str(GCM)+'_lagcor_convec_transport_0_1000_WS.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
#
#
#fig=plt.figure()
#data=[]
#for i in range(len(SMOC_transport[0][20:161])):
#    stream=runningMeanFast(SMOC_transport[:,i+lag_time], 10)
#    r=lag_cor_data(WS_index_rm,stream,lag_time)
#    data.append(r)
#data=np.asarray(data)
#x=np.linspace(-lag_time,lag_time+1, 2*lag_time)
#y=np.linspace(-70,70+1, 2*70+1)
#x,y=np.meshgrid(x,y)
#cmap_limit=np.nanmax(np.abs(data))
#levels=np.linspace(-cmap_limit,cmap_limit,200)
#im=plt.contourf(x,y,data,levels,latlon=True, cmap=plt.cm.jet)
#cb = plt.colorbar(im)
#plt.xlabel('Year lag', fontsize=18)
#plt.ylabel('Latitude', fontsize=18)
#plt.title('WS peak convection vs. transport in 2000m-3000m - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
#mng = plt.get_current_fig_manager()
#mng.window.showMaximized() # Maximizes the plot window to save figures in full
#fig.savefig(dir_figs+str(GCM)+'_lagcor_convec_transport_2000_3000_WS.png', format='png', dpi=300, transparent=True, bbox_inches='tight'
#











