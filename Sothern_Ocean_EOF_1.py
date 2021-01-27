### IMPORTANT: The following path needs to be appended so that the neccessary codes can be read from that directory ###
import sys
sys.path.append('/data1/home/basadieh/behzadcodes/behzadlibrary')
#from Behzadlib import * ### This is loaded from '/data1/home/basadieh/behzadcodes/behzadlibrary'
#from CMIP5lib import *  ### This is loaded from '/data1/home/basadieh/behzadcodes/behzadlibrary'
from Behzadlib import func_MLD, func_time_depth_plot, func_stream, func_ENSO, func_NAO, func_latlon_regrid, func_latlon_regrid_eq, func_oceanlandmask, func_oceanindex, func_regrid, func_plot, func_plot_save, func_plot_bounds_save, func_EOF
from BehzadlibPlot import func_plot2Dcolor_contourf, func_plotline_1var, func_plotline_2var, func_plot2Dcontour_1var, func_plot2Dcontour_2var, func_plotmap_contourf, plot_PSD_welch_conf
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
from scipy import stats, signal
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
Density_allyears=np.load('Density_allyears_GFDL-ESM2G_500yr.npy')

###############################################################################
##########    Labrador Sea Density EOFs     ###################################
###############################################################################

c_a=[140,156,295,320] # Cell area range [lat_S, lat_N, lon_E, lon_W]

Density_allyears_LAB=copy.deepcopy(Density_allyears[:,:,c_a[0]:c_a[1],c_a[2]:c_a[3]]) ## [65W-41W, 50N-65N] - Western Labrador box only
Density_allyears_LAB500m  = ( (Density_allyears_LAB[:,30,:,:]*53.04 + Density_allyears_LAB[:,29,:,:]*24.17 ) / (53.04+24.17) ) # Rho_500 - Interpolating between Depth[30]=524.17m and Depth[29]=446.94m to calculate density at depth 500m

Calc_Var = Density_allyears_LAB500m
Calc_Lat = Lat_regrid_2D [c_a[0]:c_a[1],c_a[2]:c_a[3]]
Calc_Lon = Lon_regrid_2D [c_a[0]:c_a[1],c_a[2]:c_a[3]]

EOF_spatial_pattern_LAB500m, EOF_time_series_LAB500m, EOF_variance_prcnt_LAB500m = func_EOF (Calc_Var, Calc_Lat)

#########################

Plot_Var = EOF_spatial_pattern_LAB500m
Plot_prcnt=EOF_variance_prcnt_LAB500m

P_cmap=plt.cm.seismic; P_proj='cyl'; P_lon0=330.; P_lon_range=90.; P_latN=90.; P_latS=0.; P_range=51; P_Lon=Lon_regrid_2D; P_Lat=Lat_regrid_2D;

n_r=2 ; n_c=2 ; n_t=4
fig=plt.figure()
for M_i in range(n_t):
    ax = fig.add_subplot(n_r,n_c,M_i+1) 

    Var_plot_ii = empty((Lat_regrid_2D.shape[0],Lat_regrid_2D.shape[1]))*nan    
    Var_plot_ii[c_a[0]:c_a[1],c_a[2]:c_a[3]]=Plot_Var[M_i,:,:]    
    
    m = Basemap( projection=P_proj, lon_0=P_lon0, llcrnrlon=P_lon0-P_lon_range, llcrnrlat=P_latS, urcrnrlon=P_lon0+P_lon_range, urcrnrlat=P_latN)    
    if M_i == (n_c*(n_r-1)): # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        m.drawparallels(np.arange(P_latS, P_latN+0.001, 30.),labels=[True,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Latitutes
        m.drawmeridians(np.arange(0,360,90.),labels=[False,False,False,True], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Longitudes    
    elif M_i==0 or M_i==n_c or M_i==n_c*2 or M_i==n_c*3 or M_i==n_c*4 or M_i==n_c*5 or M_i==n_c*6 or M_i==n_c*7 or M_i==n_c*8:
        m.drawparallels(np.arange(P_latS, P_latN+0.001, 30.),labels=[True,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Latitutes
    elif M_i >= n_t-n_c and M_i != (n_c*(n_r-1)): # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        m.drawmeridians(np.arange(0,360,90.),labels=[False,False,False,True], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Longitudes
    else:
        m.drawparallels(np.arange(P_latS, P_latN+0.001, 30.),labels=[False,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Latitutes
        m.drawmeridians(np.arange(0,360,90.),labels=[False,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Longitudes

    m.fillcontinents(color='0')
    m.drawcoastlines(linewidth=1.0, linestyle='solid', antialiased=1, ax=None, zorder=None)
    m.fillcontinents(color='0.95')
    im=m.contourf(P_Lon, P_Lat, Var_plot_ii, P_range, latlon=True, cmap=P_cmap, extend='both')
    plt.title('EOF #'+str(M_i+1)+'  ,  '+str(round(Plot_prcnt[M_i], 2))+' % of the variance', fontsize=14)
        
plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, hspace=0.05, wspace=0.1) # the amount of height/width reserved for space between subplots
cbar_ax = fig.add_axes([0.93, 0.1, 0.015, 0.8]) # [right,bottom,width,height] 
fig.colorbar(im, cax=cbar_ax)
plt.suptitle('EOF spatial pattern of potential density at 500m - Labrador Sea [65W-41W, 50N-65N] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=20)    
plt.show()
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full   

fig.savefig(dir_figs+str(GCM)+'_EOF_map_Rho500m_LAB1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


Plot_Var = EOF_time_series_LAB500m
Plot_prcnt=EOF_variance_prcnt_LAB500m
fig, ax = plt.subplots(nrows=2, ncols=2)
for ii in range(4):
    EOF_time_series_plot_norm=(Plot_Var[ii,:]-np.nanmean(Plot_Var[ii,:]))/np.std(Plot_Var[ii,:])    
    EOF_time_series_plot_norm_rm=runningMeanFast(EOF_time_series_plot_norm, 10)
    plt.subplot(2, 2, ii+1)
    n_l=EOF_time_series_plot_norm.shape[0]
    years=np.linspace(0, n_l, n_l)
    plt.plot(years,EOF_time_series_plot_norm_rm, 'k') 
    y2=np.zeros(len(EOF_time_series_plot_norm_rm))
    plt.fill_between(years, EOF_time_series_plot_norm_rm, y2, where=EOF_time_series_plot_norm_rm >= y2, color = 'r', interpolate=True)
    plt.fill_between(years, EOF_time_series_plot_norm_rm, y2, where=EOF_time_series_plot_norm_rm <= y2, color = 'b', interpolate=True)
    plt.axhline(linewidth=1, color='k')
    plt.title('EOF # '+str(ii+1)+'  ,  '+str(round(Plot_prcnt[ii], 2))+' % of the variance', fontsize=18)
    plt.xticks(fontsize = 18); plt.yticks(fontsize = 18)
    plt.suptitle('EOF normalized/smoothed index of potential density at 500m - Labrador Sea [65W-41W, 50N-65N] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)  
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full

fig.savefig(dir_figs+str(GCM)+'_EOF_timeseries_Rho500m_LAB1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

###############################################################################
##########    Weddell Sea Density EOFs     ####################################
###############################################################################
c_a=[20,35,320,360] # Cell area range [lat_S, lat_N, lon_E, lon_W] # [40W-0W,70S-55S]

Density_allyears_WS=copy.deepcopy(Density_allyears[:,:,c_a[0]:c_a[1],c_a[2]:c_a[3]]) ## [40W-0W,70S-55S] -Weddell Sea
Density_allyears_WS500m  = ( (Density_allyears_WS[:,30,:,:]*53.04 + Density_allyears_WS[:,29,:,:]*24.17 ) / (53.04+24.17) ) # Rho_500 - Interpolating between Depth[30]=524.17m and Depth[29]=446.94m to calculate density at depth 500m

Calc_Var = Density_allyears_WS500m
Calc_Lat = Lat_regrid_2D [c_a[0]:c_a[1],c_a[2]:c_a[3]]
Calc_Lon = Lon_regrid_2D [c_a[0]:c_a[1],c_a[2]:c_a[3]]

EOF_spatial_pattern_WS500m, EOF_time_series_WS500m, EOF_variance_prcnt_WS500m = func_EOF (Calc_Var, Calc_Lat)

#########################

Plot_Var = EOF_spatial_pattern_WS500m
Plot_prcnt=EOF_variance_prcnt_WS500m

P_cmap=plt.cm.seismic; P_proj='cyl'; P_lon0=330.; P_lon_range=90.; P_latN=0.; P_latS=-90.; P_range=51; P_Lon=Lon_regrid_2D; P_Lat=Lat_regrid_2D;

n_r=2 ; n_c=2 ; n_t=4
fig=plt.figure()
for M_i in range(n_t):
    ax = fig.add_subplot(n_r,n_c,M_i+1) 

    Var_plot_ii = empty((Lat_regrid_2D.shape[0],Lat_regrid_2D.shape[1]))*nan    
    Var_plot_ii[c_a[0]:c_a[1],c_a[2]:c_a[3]]=Plot_Var[M_i,:,:]    
    
    m = Basemap( projection=P_proj, lon_0=P_lon0, llcrnrlon=P_lon0-P_lon_range, llcrnrlat=P_latS, urcrnrlon=P_lon0+P_lon_range, urcrnrlat=P_latN)    
    if M_i == (n_c*(n_r-1)): # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        m.drawparallels(np.arange(P_latS, P_latN+0.001, 30.),labels=[True,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Latitutes
        m.drawmeridians(np.arange(0,360,90.),labels=[False,False,False,True], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Longitudes    
    elif M_i==0 or M_i==n_c or M_i==n_c*2 or M_i==n_c*3 or M_i==n_c*4 or M_i==n_c*5 or M_i==n_c*6 or M_i==n_c*7 or M_i==n_c*8:
        m.drawparallels(np.arange(P_latS, P_latN+0.001, 30.),labels=[True,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Latitutes
    elif M_i >= n_t-n_c and M_i != (n_c*(n_r-1)): # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        m.drawmeridians(np.arange(0,360,90.),labels=[False,False,False,True], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Longitudes
    else:
        m.drawparallels(np.arange(P_latS, P_latN+0.001, 30.),labels=[False,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Latitutes
        m.drawmeridians(np.arange(0,360,90.),labels=[False,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Longitudes

    m.fillcontinents(color='0')
    m.drawcoastlines(linewidth=1.0, linestyle='solid', antialiased=1, ax=None, zorder=None)
    m.fillcontinents(color='0.95')
    im=m.contourf(P_Lon, P_Lat, Var_plot_ii, P_range, latlon=True, cmap=P_cmap, extend='both')
    plt.title('EOF #'+str(M_i+1)+'  ,  '+str(round(Plot_prcnt[M_i], 2))+' % of the variance', fontsize=14)
        
plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, hspace=0.05, wspace=0.1) # the amount of height/width reserved for space between subplots
cbar_ax = fig.add_axes([0.93, 0.1, 0.015, 0.8]) # [right,bottom,width,height] 
fig.colorbar(im, cax=cbar_ax)
plt.suptitle('EOF spatial pattern of potential density at 500m - Weddell Sea [40W-0W,70S-55S] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=20)    
plt.show()
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full   

fig.savefig(dir_figs+str(GCM)+'_EOF_map_Rho500m_WS1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


Plot_Var = EOF_time_series_WS500m
Plot_prcnt=EOF_variance_prcnt_WS500m
fig, ax = plt.subplots(nrows=2, ncols=2)
for ii in range(4):
    EOF_time_series_plot_norm=(Plot_Var[ii,:]-np.nanmean(Plot_Var[ii,:]))/np.std(Plot_Var[ii,:])    
    EOF_time_series_plot_norm_rm=runningMeanFast(EOF_time_series_plot_norm, 10)
    plt.subplot(2, 2, ii+1)
    n_l=EOF_time_series_plot_norm.shape[0]
    years=np.linspace(0, n_l, n_l)
    plt.plot(years,EOF_time_series_plot_norm_rm, 'k') 
    y2=np.zeros(len(EOF_time_series_plot_norm_rm))
    plt.fill_between(years, EOF_time_series_plot_norm_rm, y2, where=EOF_time_series_plot_norm_rm >= y2, color = 'r', interpolate=True)
    plt.fill_between(years, EOF_time_series_plot_norm_rm, y2, where=EOF_time_series_plot_norm_rm <= y2, color = 'b', interpolate=True)
    plt.axhline(linewidth=1, color='k')
    plt.title('EOF # '+str(ii+1)+'  ,  '+str(round(Plot_prcnt[ii], 2))+' % of the variance', fontsize=18)
    plt.xticks(fontsize = 18); plt.yticks(fontsize = 18)
    plt.suptitle('EOF normalized/smoothed index of potential density at 500m - Weddell Sea [40W-0W,70S-55S] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)  
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full

fig.savefig(dir_figs+str(GCM)+'_EOF_timeseries_Rho500m_WS1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

###############################################################################
##########    Atlantic Ocean Density EOFs     #################################
###############################################################################

Density_allyears_Atl=copy.deepcopy(Density_allyears)
Density_allyears_Atl_500m  = ( (Density_allyears_Atl[:,30,:,:]*53.04 + Density_allyears_Atl[:,29,:,:]*24.17 ) / (53.04+24.17) ) # Rho_500 - Interpolating between Depth[30]=524.17m and Depth[29]=446.94m to calculate density at depth 500m
for ii in range(Density_allyears_Atl_500m.shape[0]): # Only keeping Atlantic cells
    Density_allyears_Atl_500m[ii,:,:][ np.where( np.logical_and(   np.logical_and( Ocean_Index != 6 , Ocean_Index != 7, Ocean_Index != 8 ) , np.logical_and( Ocean_Index != 8 , Ocean_Index != 9 )   ) ) ] = nan

Calc_Var = Density_allyears_Atl_500m
Calc_Lat = Lat_regrid_2D
Calc_Lon = Lon_regrid_2D

EOF_spatial_pattern_Atl500m, EOF_time_series_Atl500m, EOF_variance_prcnt_Atl500m = func_EOF (Calc_Var, Calc_Lat)

######################

Plot_Var = EOF_spatial_pattern_Atl500m
Plot_prcnt=EOF_variance_prcnt_Atl500m

P_cmap=plt.cm.seismic; P_proj='cyl'; P_lon0=0.; P_latN=90.; P_latS=-90.; P_range=51; P_Lon=Lon_regrid_2D; P_Lat=Lat_regrid_2D;

n_r=2 ; n_c=2 ; n_t=4
fig=plt.figure()
for M_i in range(n_t):
    ax = fig.add_subplot(n_r,n_c,M_i+1)     
    Var_plot_ii=Plot_Var[M_i,:,:]    
    
    m = Basemap( projection=P_proj, lon_0=P_lon0, llcrnrlon=P_lon0-180, llcrnrlat=P_latS, urcrnrlon=P_lon0+180, urcrnrlat=P_latN)    
    if M_i == (n_c*(n_r-1)): # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        m.drawparallels(np.arange(P_latS, P_latN+0.001, 30.),labels=[True,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Latitutes
        m.drawmeridians(np.arange(0,360,90.),labels=[False,False,False,True], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Longitudes    
    elif M_i==0 or M_i==n_c or M_i==n_c*2 or M_i==n_c*3 or M_i==n_c*4 or M_i==n_c*5 or M_i==n_c*6 or M_i==n_c*7 or M_i==n_c*8:
        m.drawparallels(np.arange(P_latS, P_latN+0.001, 30.),labels=[True,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Latitutes
    elif M_i >= n_t-n_c and M_i != (n_c*(n_r-1)): # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        m.drawmeridians(np.arange(0,360,90.),labels=[False,False,False,True], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Longitudes
    else:
        m.drawparallels(np.arange(P_latS, P_latN+0.001, 30.),labels=[False,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Latitutes
        m.drawmeridians(np.arange(0,360,90.),labels=[False,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Longitudes

    m.fillcontinents(color='0')
    m.drawcoastlines(linewidth=1.0, linestyle='solid', antialiased=1, ax=None, zorder=None)
    m.fillcontinents(color='0.95')
    im=m.contourf(P_Lon, P_Lat, Var_plot_ii, P_range, latlon=True, cmap=P_cmap, extend='both')
    plt.title('EOF #'+str(M_i+1)+'  ,  '+str(round(Plot_prcnt[M_i], 2))+' % of the variance', fontsize=14)
        
plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, hspace=0.05, wspace=0.1) # the amount of height/width reserved for space between subplots
cbar_ax = fig.add_axes([0.93, 0.1, 0.015, 0.8]) # [right,bottom,width,height] 
fig.colorbar(im, cax=cbar_ax)
plt.suptitle('EOF spatial pattern of potential density at 500m - Atlantic Ocean '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=20)    
plt.show()
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full   

fig.savefig(dir_figs+str(GCM)+'_EOF_map_Rho500m_atlantic.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


Plot_Var = EOF_time_series_Atl500m
Plot_prcnt=EOF_variance_prcnt_Atl500m
fig, ax = plt.subplots(nrows=2, ncols=2)
for ii in range(4):
    EOF_time_series_plot_norm=(Plot_Var[ii,:]-np.nanmean(Plot_Var[ii,:]))/np.std(Plot_Var[ii,:])    
    EOF_time_series_plot_norm_rm=runningMeanFast(EOF_time_series_plot_norm, 10)
    plt.subplot(2, 2, ii+1)
    n_l=EOF_time_series_plot_norm.shape[0]
    years=np.linspace(0, n_l, n_l)
    plt.plot(years,EOF_time_series_plot_norm_rm, 'k') 
    y2=np.zeros(len(EOF_time_series_plot_norm_rm))
    plt.fill_between(years, EOF_time_series_plot_norm_rm, y2, where=EOF_time_series_plot_norm_rm >= y2, color = 'r', interpolate=True)
    plt.fill_between(years, EOF_time_series_plot_norm_rm, y2, where=EOF_time_series_plot_norm_rm <= y2, color = 'b', interpolate=True)
    plt.axhline(linewidth=1, color='k')
    plt.title('EOF # '+str(ii+1)+'  ,  '+str(round(Plot_prcnt[ii], 2))+' % of the variance', fontsize=18)
    plt.xticks(fontsize = 18); plt.yticks(fontsize = 18)
    plt.suptitle('EOF normalized/smoothed index of potential density at 500m - Atlantic Ocean - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)  
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full

fig.savefig(dir_figs+str(GCM)+'_EOF_timeseries_Rho500m_atlantic.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


###############################################################################
###############################################################################
from scipy import stats, signal
import os
import shelve

dir_pwd = os.getcwd() # Gets the current directory (and in which the code is placed)
filename_out = (dir_pwd + '/Results/AllResults_'+GCM+'_500yr.out') # Directory to save processed data

my_shelf = shelve.open(filename_out)
#globals()[key]=my_shelf[key]
globals()['WS_index']=my_shelf['WS_index']
globals()['LAB_index']=my_shelf['LAB_index']
globals()['WS_index_norm']=my_shelf['WS_index_norm']
globals()['WS_index_norm_rm']=my_shelf['WS_index_norm_rm']
globals()['LAB_index_norm']=my_shelf['LAB_index_norm']
globals()['LAB_index_norm_rm']=my_shelf['LAB_index_norm_rm']
globals()['AMOC_max_50S']=my_shelf['AMOC_max_50S']
globals()['AMOC_max_30S']=my_shelf['AMOC_max_30S']
globals()['AMOC_max_30N']=my_shelf['AMOC_max_30N']
globals()['AMOC_max_50N']=my_shelf['AMOC_max_50N']
globals()['NAO_index']=my_shelf['NAO_index']
globals()['ENSO_index']=my_shelf['ENSO_index']
my_shelf.close()

CutOff_T = 5 # Cut-off period 
n_order = 4 # Order of filtering

fs = 1  # Sampling frequency, equal to 1 year in our case
fc = 1/CutOff_T  # Cut-off frequency of the filter
ww = fc / (fs / 2) # Normalize the frequency
bb, aa = signal.butter(n_order, ww, 'low')

WS_index_BWfilt = signal.filtfilt(bb, aa, WS_index)
LAB_index_BWfilt = signal.filtfilt(bb, aa, LAB_index)

AMOC_max_50S_BWfilt = signal.filtfilt(bb, aa, AMOC_max_50S)
AMOC_max_30S_BWfilt = signal.filtfilt(bb, aa, AMOC_max_30S)
AMOC_max_50N_BWfilt = signal.filtfilt(bb, aa, AMOC_max_50N)
AMOC_max_30N_BWfilt = signal.filtfilt(bb, aa, AMOC_max_30N)
ENSO_index_BWfilt = signal.filtfilt(bb, aa, ENSO_index)
NAO_index_BWfilt = signal.filtfilt(bb, aa, NAO_index)


##################################################
#######   Power Spectral Density plots    ########
##################################################

P_title='Power Spectral Density - Labrador Sea [65W-41W, 50N-65N] EOF of density at 500m - BW low-pass filt., '+str(CutOff_T)+'yr cut-off - '+str(GCM)
fig=plt.figure(); ax = fig.add_axes([0.12,0.1,0.78,0.8]) # ax [left, bottom, width, height]
ff,PSD=plot_PSD_welch_conf(signal.filtfilt(bb, aa, EOF_time_series_LAB500m[0,:]), 1, '-', '-', '-', P_title, 'b', 'LAB 500m Density 1st EOF', 'lower right')
ff,PSD=plot_PSD_welch_conf(signal.filtfilt(bb, aa, EOF_time_series_LAB500m[1,:]), 1, '-', '-', '-', P_title, 'g', 'LAB 500m Density 2nd EOF', 'lower right')
ff,PSD=plot_PSD_welch_conf(signal.filtfilt(bb, aa, EOF_time_series_LAB500m[2,:]), 1, '-', '-', '-', P_title, 'darkorange', 'LAB 500m Density 3rd EOF', 'lower right')
ff,PSD=plot_PSD_welch_conf(signal.filtfilt(bb, aa, EOF_time_series_LAB500m[3,:]), 1, '-', '-', '-', P_title, 'r', 'LAB 500m Density 4th EOF', 'lower right')
#plt.ylim(1e-6,9.8) #;plt.xlim(1,500)
fig.savefig(dir_figs+str(GCM)+'_PSD_LAB_Rho500m_EOF_BW'+str(CutOff_T)+'.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

P_title='Power Spectral Density - Weddell Sea [40W-0W,70S-55S] EOF of density at 500m - BW low-pass filt., '+str(CutOff_T)+'yr cut-off - '+str(GCM)
fig=plt.figure(); ax = fig.add_axes([0.12,0.1,0.78,0.8]) # ax [left, bottom, width, height]
ff,PSD=plot_PSD_welch_conf(signal.filtfilt(bb, aa, EOF_time_series_WS500m[0,:]), 1, '-', '-', '-', P_title, 'b', 'WS 500m Density 1st EOF', 'lower right')
ff,PSD=plot_PSD_welch_conf(signal.filtfilt(bb, aa, EOF_time_series_WS500m[1,:]), 1, '-', '-', '-', P_title, 'g', 'WS 500m Density 2nd EOF', 'lower right')
ff,PSD=plot_PSD_welch_conf(signal.filtfilt(bb, aa, EOF_time_series_WS500m[2,:]), 1, '-', '-', '-', P_title, 'darkorange', 'WS 500m Density 3rd EOF', 'lower right')
ff,PSD=plot_PSD_welch_conf(signal.filtfilt(bb, aa, EOF_time_series_WS500m[3,:]), 1, '-', '-', '-', P_title, 'r', 'WS 500m Density 4th EOF', 'lower right')
#plt.ylim(1e-6,9.8) #;plt.xlim(1,500)
fig.savefig(dir_figs+str(GCM)+'_PSD_WS_Rho500m_EOF_BW'+str(CutOff_T)+'.png', format='png', dpi=300, transparent=True, bbox_inches='tight')



P_Var_x1=P_Var_x2=years
P_Var_y1=copy.deepcopy(EOF_time_series_LAB500m[2,:])*-1
P_Var_y2=copy.deepcopy(EOF_time_series_LAB500m[2,:])*-1
P_xlable='Years' ;  P_ylable='Convection index (°C)'
P_title='Labrador Sea conv. index - Butterworth low-pass filt. '+str(CutOff_T)+'yr cut-off, '+str(n_order)+'th order (normalized, NOT-smoothed) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plotline_2var(P_Var_x1, P_Var_y1, P_Var_x2, P_Var_y2, P_xlable, P_ylable, P_title, 'b', 'r', '1st EOF', '2nd EOF', 'best', '-')
#fig.savefig(dir_figs+str(GCM)+'_convec_index_LAB_BW'+str(CutOff_T)+'.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig = func_plotline_1var(P_Var_x1, P_Var_y1, P_xlable, P_ylable, P_title, 'darkcyan', 'invert_yaxis')

###############################################################################
###########           Lagged Correlations          ############################
###############################################################################

fig=plt.figure()
lag_cor(LAB_index_BWfilt, signal.filtfilt(bb, aa, EOF_time_series_LAB500m[0,:]), 40,'b','LAB 500m Density 1st EOF')
lag_cor(LAB_index_BWfilt, signal.filtfilt(bb, aa, EOF_time_series_LAB500m[1,:]), 40,'g','LAB 500m Density 2nd EOF')
lag_cor(LAB_index_BWfilt, signal.filtfilt(bb, aa, EOF_time_series_LAB500m[2,:]), 40,'darkorange','LAB 500m Density 3rd EOF')
lag_cor(LAB_index_BWfilt, signal.filtfilt(bb, aa, EOF_time_series_LAB500m[3,:]), 40,'r','LAB 500m Density 4th EOF')
plt.legend()
plt.show()
plt.title('Peak LAB Conv. lagged corr. with EOF of density at 500m - BW low-pass filt., '+str(CutOff_T)+'yr cut-off) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
plt.xlabel('Years lag', fontsize=18); plt.xticks(fontsize=18)
plt.ylabel('Correlation coefficient', fontsize=18); plt.yticks(fontsize=18)    
plt.show()
plt.legend(shadow=True, prop={'size': 15})
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_LABconvection_Rho500m_EOF_BW'+str(CutOff_T)+'_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
lag_cor(WS_index_BWfilt, signal.filtfilt(bb, aa, EOF_time_series_WS500m[0,:]), 40,'b','WS 500m Density 1st EOF')
lag_cor(WS_index_BWfilt, signal.filtfilt(bb, aa, EOF_time_series_WS500m[1,:]), 40,'g','WS 500m Density 2nd EOF')
lag_cor(WS_index_BWfilt, signal.filtfilt(bb, aa, EOF_time_series_WS500m[2,:]), 40,'darkorange','WS 500m Density 3rd EOF')
lag_cor(WS_index_BWfilt, signal.filtfilt(bb, aa, EOF_time_series_WS500m[3,:]), 40,'r','WS 500m Density 4th EOF')
plt.legend()
plt.show()
plt.title('Peak WS Conv. lagged corr. with EOF of density at 500m  (BW low-pass filt., '+str(CutOff_T)+'yr cut-off) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
plt.xlabel('Years lag', fontsize=18); plt.xticks(fontsize=18)
plt.ylabel('Correlation coefficient', fontsize=18); plt.yticks(fontsize=18)    
plt.show()
plt.legend(shadow=True, prop={'size': 15})
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_WSconvection_Rho500m_EOF_BW'+str(CutOff_T)+'_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

######################################

fig=plt.figure()
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_LAB500m[0,:]), WS_index_BWfilt, 40,'y','WS index')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_LAB500m[0,:]), LAB_index_BWfilt, 40,'k','LAB index')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_LAB500m[0,:]), AMOC_max_50N_BWfilt, 40,'b','AMOC max 50N')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_LAB500m[0,:]), AMOC_max_30N_BWfilt, 40,'g','AMOC max 30N')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_LAB500m[0,:]), AMOC_max_30S_BWfilt, 40,'darkorange','AMOC max 30S')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_LAB500m[0,:]), AMOC_max_50S_BWfilt, 40,'r','AMOC max 50S')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_LAB500m[0,:]), ENSO_index_BWfilt, 40,'m','ENSO index')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_LAB500m[0,:]), NAO_index_BWfilt, 40,'c','NAO index')
plt.legend()
plt.show()
plt.title('1st EOF of Labrador Sea density at 500m lagged correlation (BW low-pass filt., '+str(CutOff_T)+'yr cut-off) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
plt.xlabel('Years lag', fontsize=18); plt.xticks(fontsize=18)
plt.ylabel('Correlation coefficient', fontsize=18); plt.yticks(fontsize=18)    
plt.show()
plt.legend(shadow=True, prop={'size': 15})
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_LAB_Rho500m_1st_EOF_BW'+str(CutOff_T)+'_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_LAB500m[1,:]), WS_index_BWfilt, 40,'y','WS index')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_LAB500m[1,:]), LAB_index_BWfilt, 40,'k','LAB index')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_LAB500m[1,:]), AMOC_max_50N_BWfilt, 40,'b','AMOC max 50N')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_LAB500m[1,:]), AMOC_max_30N_BWfilt, 40,'g','AMOC max 30N')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_LAB500m[1,:]), AMOC_max_30S_BWfilt, 40,'darkorange','AMOC max 30S')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_LAB500m[1,:]), AMOC_max_50S_BWfilt, 40,'r','AMOC max 50S')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_LAB500m[1,:]), ENSO_index_BWfilt, 40,'m','ENSO index')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_LAB500m[1,:]), NAO_index_BWfilt, 40,'c','NAO index')
plt.legend()
plt.show()
plt.title('2nd EOF of Labrador Sea density at 500m lagged correlation (BW low-pass filt., '+str(CutOff_T)+'yr cut-off) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
plt.xlabel('Years lag', fontsize=18); plt.xticks(fontsize=18)
plt.ylabel('Correlation coefficient', fontsize=18); plt.yticks(fontsize=18)    
plt.show()
plt.legend(shadow=True, prop={'size': 15})
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_LAB_Rho500m_2nd_EOF_BW'+str(CutOff_T)+'_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


fig=plt.figure()
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_LAB500m[2,:]), WS_index_BWfilt, 40,'y','WS index')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_LAB500m[2,:]), LAB_index_BWfilt, 40,'k','LAB index')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_LAB500m[2,:]), AMOC_max_50N_BWfilt, 40,'b','AMOC max 50N')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_LAB500m[2,:]), AMOC_max_30N_BWfilt, 40,'g','AMOC max 30N')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_LAB500m[2,:]), AMOC_max_30S_BWfilt, 40,'darkorange','AMOC max 30S')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_LAB500m[2,:]), AMOC_max_50S_BWfilt, 40,'r','AMOC max 50S')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_LAB500m[2,:]), ENSO_index_BWfilt, 40,'m','ENSO index')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_LAB500m[2,:]), NAO_index_BWfilt, 40,'c','NAO index')
plt.legend()
plt.show()
plt.title('3rd EOF of Labrador Sea density at 500m lagged correlation (BW low-pass filt., '+str(CutOff_T)+'yr cut-off) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
plt.xlabel('Years lag', fontsize=18); plt.xticks(fontsize=18)
plt.ylabel('Correlation coefficient', fontsize=18); plt.yticks(fontsize=18)    
plt.show()
plt.legend(shadow=True, prop={'size': 15})
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_LAB_Rho500m_3rd_EOF_BW'+str(CutOff_T)+'_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

######################################

fig=plt.figure()
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_WS500m[0,:]), WS_index_BWfilt, 40,'y','WS index')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_WS500m[0,:]), LAB_index_BWfilt, 40,'k','LAB index')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_WS500m[0,:]), AMOC_max_50N_BWfilt, 40,'b','AMOC max 50N')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_WS500m[0,:]), AMOC_max_30N_BWfilt, 40,'g','AMOC max 30N')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_WS500m[0,:]), AMOC_max_30S_BWfilt, 40,'darkorange','AMOC max 30S')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_WS500m[0,:]), AMOC_max_50S_BWfilt, 40,'r','AMOC max 50S')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_WS500m[0,:]), ENSO_index_BWfilt, 40,'m','ENSO index')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_WS500m[0,:]), NAO_index_BWfilt, 40,'c','NAO index')
plt.legend()
plt.show()
plt.title('1st EOF of Weddell Sea density at 500m lagged correlation (BW low-pass filt., '+str(CutOff_T)+'yr cut-off) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
plt.xlabel('Years lag', fontsize=18); plt.xticks(fontsize=18)
plt.ylabel('Correlation coefficient', fontsize=18); plt.yticks(fontsize=18)    
plt.show()
plt.legend(shadow=True, prop={'size': 15})
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_WS_Rho500m_1st_EOF_BW'+str(CutOff_T)+'_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_WS500m[1,:]), WS_index_BWfilt, 40,'y','WS index')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_WS500m[1,:]), LAB_index_BWfilt, 40,'k','LAB index')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_WS500m[1,:]), AMOC_max_50N_BWfilt, 40,'b','AMOC max 50N')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_WS500m[1,:]), AMOC_max_30N_BWfilt, 40,'g','AMOC max 30N')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_WS500m[1,:]), AMOC_max_30S_BWfilt, 40,'darkorange','AMOC max 30S')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_WS500m[1,:]), AMOC_max_50S_BWfilt, 40,'r','AMOC max 50S')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_WS500m[1,:]), ENSO_index_BWfilt, 40,'m','ENSO index')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_WS500m[1,:]), NAO_index_BWfilt, 40,'c','NAO index')
plt.legend()
plt.show()
plt.title('2nd EOF of Weddell Sea density at 500m lagged correlation (BW low-pass filt., '+str(CutOff_T)+'yr cut-off) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
plt.xlabel('Years lag', fontsize=18); plt.xticks(fontsize=18)
plt.ylabel('Correlation coefficient', fontsize=18); plt.yticks(fontsize=18)    
plt.show()
plt.legend(shadow=True, prop={'size': 15})
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_WS_Rho500m_2nd_EOF_BW'+str(CutOff_T)+'_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_WS500m[2,:]), WS_index_BWfilt, 40,'y','WS index')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_WS500m[2,:]), LAB_index_BWfilt, 40,'k','LAB index')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_WS500m[2,:]), AMOC_max_50N_BWfilt, 40,'b','AMOC max 50N')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_WS500m[2,:]), AMOC_max_30N_BWfilt, 40,'g','AMOC max 30N')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_WS500m[2,:]), AMOC_max_30S_BWfilt, 40,'darkorange','AMOC max 30S')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_WS500m[2,:]), AMOC_max_50S_BWfilt, 40,'r','AMOC max 50S')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_WS500m[2,:]), ENSO_index_BWfilt, 40,'m','ENSO index')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_WS500m[2,:]), NAO_index_BWfilt, 40,'c','NAO index')
plt.legend()
plt.show()
plt.title('3rd EOF of Weddell Sea density at 500m lagged correlation (BW low-pass filt., '+str(CutOff_T)+'yr cut-off) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
plt.xlabel('Years lag', fontsize=18); plt.xticks(fontsize=18)
plt.ylabel('Correlation coefficient', fontsize=18); plt.yticks(fontsize=18)    
plt.show()
plt.legend(shadow=True, prop={'size': 15})
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_WS_Rho500m_3rd_EOF_BW'+str(CutOff_T)+'_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')



###############################################################################
###############################################################################
Temp_allyears=np.load('Temp_allyears_GFDL-ESM2G_500yr.npy')

###############################################################################
##########    Labrador Sea Density EOFs     ###################################
###############################################################################

c_a=[140,156,295,320] # Cell area range [lat_S, lat_N, lon_E, lon_W]

Temp_allyears_LAB=copy.deepcopy(Temp_allyears[:,:,c_a[0]:c_a[1],c_a[2]:c_a[3]]) ## [65W-41W, 50N-65N] - Western Labrador box only
Temp_allyears_LAB500m  = ( (Temp_allyears_LAB[:,30,:,:]*53.04 + Temp_allyears_LAB[:,29,:,:]*24.17 ) / (53.04+24.17) ) # Temp_500 - Interpolating between Depth[30]=524.17m and Depth[29]=446.94m to calculate density at depth 500m

Calc_Var = Temp_allyears_LAB500m
Calc_Lat = Lat_regrid_2D [c_a[0]:c_a[1],c_a[2]:c_a[3]]
Calc_Lon = Lon_regrid_2D [c_a[0]:c_a[1],c_a[2]:c_a[3]]

EOF_spatial_pattern_T_LAB500m, EOF_time_series_T_LAB500m, EOF_variance_prcnt_T_LAB500m = func_EOF (Calc_Var, Calc_Lat)

EOF_spatial_pattern_T_LAB0m, EOF_time_series_T_LAB0m, EOF_variance_prcnt_T_LAB0m = func_EOF (Temp_allyears_LAB[:,0,:,:], Calc_Lat)

#########################

Plot_Var = EOF_spatial_pattern_T_LAB0m
Plot_prcnt=EOF_variance_prcnt_T_LAB0m

P_cmap=plt.cm.seismic; P_proj='cyl'; P_lon0=330.; P_lon_range=90.; P_latN=90.; P_latS=0.; P_range=51; P_Lon=Lon_regrid_2D; P_Lat=Lat_regrid_2D;

n_r=2 ; n_c=2 ; n_t=4
fig=plt.figure()
for M_i in range(n_t):
    ax = fig.add_subplot(n_r,n_c,M_i+1) 

    Var_plot_ii = empty((Lat_regrid_2D.shape[0],Lat_regrid_2D.shape[1]))*nan    
    Var_plot_ii[c_a[0]:c_a[1],c_a[2]:c_a[3]]=Plot_Var[M_i,:,:]    
    
    m = Basemap( projection=P_proj, lon_0=P_lon0, llcrnrlon=P_lon0-P_lon_range, llcrnrlat=P_latS, urcrnrlon=P_lon0+P_lon_range, urcrnrlat=P_latN)    
    if M_i == (n_c*(n_r-1)): # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        m.drawparallels(np.arange(P_latS, P_latN+0.001, 30.),labels=[True,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Latitutes
        m.drawmeridians(np.arange(0,360,90.),labels=[False,False,False,True], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Longitudes    
    elif M_i==0 or M_i==n_c or M_i==n_c*2 or M_i==n_c*3 or M_i==n_c*4 or M_i==n_c*5 or M_i==n_c*6 or M_i==n_c*7 or M_i==n_c*8:
        m.drawparallels(np.arange(P_latS, P_latN+0.001, 30.),labels=[True,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Latitutes
    elif M_i >= n_t-n_c and M_i != (n_c*(n_r-1)): # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        m.drawmeridians(np.arange(0,360,90.),labels=[False,False,False,True], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Longitudes
    else:
        m.drawparallels(np.arange(P_latS, P_latN+0.001, 30.),labels=[False,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Latitutes
        m.drawmeridians(np.arange(0,360,90.),labels=[False,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Longitudes

    m.fillcontinents(color='0')
    m.drawcoastlines(linewidth=1.0, linestyle='solid', antialiased=1, ax=None, zorder=None)
    m.fillcontinents(color='0.95')
    im=m.contourf(P_Lon, P_Lat, Var_plot_ii, P_range, latlon=True, cmap=P_cmap, extend='both')
    plt.title('EOF #'+str(M_i+1)+'  ,  '+str(round(Plot_prcnt[M_i], 2))+' % of the variance', fontsize=14)
        
plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, hspace=0.05, wspace=0.1) # the amount of height/width reserved for space between subplots
cbar_ax = fig.add_axes([0.93, 0.1, 0.015, 0.8]) # [right,bottom,width,height] 
fig.colorbar(im, cax=cbar_ax)
plt.suptitle('EOF spatial pattern of potential temperature at 500m - Labrador Sea [65W-41W, 50N-65N] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=20)    
plt.show()
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full   

fig.savefig(dir_figs+str(GCM)+'_EOF_map_Temp500m_LAB1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


Plot_Var = EOF_time_series_T_LAB500m
Plot_prcnt=EOF_variance_prcnt_T_LAB500m
fig, ax = plt.subplots(nrows=2, ncols=2)
for ii in range(4):
    EOF_time_series_plot_norm=(Plot_Var[ii,:]-np.nanmean(Plot_Var[ii,:]))/np.std(Plot_Var[ii,:])    
    EOF_time_series_plot_norm_rm=runningMeanFast(EOF_time_series_plot_norm, 10)
    plt.subplot(2, 2, ii+1)
    n_l=EOF_time_series_plot_norm.shape[0]
    years=np.linspace(0, n_l, n_l)
    plt.plot(years,EOF_time_series_plot_norm_rm, 'k') 
    y2=np.zeros(len(EOF_time_series_plot_norm_rm))
    plt.fill_between(years, EOF_time_series_plot_norm_rm, y2, where=EOF_time_series_plot_norm_rm >= y2, color = 'r', interpolate=True)
    plt.fill_between(years, EOF_time_series_plot_norm_rm, y2, where=EOF_time_series_plot_norm_rm <= y2, color = 'b', interpolate=True)
    plt.axhline(linewidth=1, color='k')
    plt.title('EOF # '+str(ii+1)+'  ,  '+str(round(Plot_prcnt[ii], 2))+' % of the variance', fontsize=18)
    plt.xticks(fontsize = 18); plt.yticks(fontsize = 18)
    plt.suptitle('EOF normalized/smoothed index of potential temperature at 500m - Labrador Sea [65W-41W, 50N-65N] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)  
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full

fig.savefig(dir_figs+str(GCM)+'_EOF_timeseries_Temp500m_LAB1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


###################################################################

CutOff_T = 5 # Cut-off period 
n_order = 4 # Order of filtering

fs = 1  # Sampling frequency, equal to 1 year in our case
fc = 1/CutOff_T  # Cut-off frequency of the filter
ww = fc / (fs / 2) # Normalize the frequency
bb, aa = signal.butter(n_order, ww, 'low')


P_title='Power Spectral Density - Labrador Sea [65W-41W, 50N-65N] EOF of temperature at 500m - BW low-pass filt., '+str(CutOff_T)+'yr cut-off - '+str(GCM)
fig=plt.figure(); ax = fig.add_axes([0.12,0.1,0.78,0.8]) # ax [left, bottom, width, height]
ff,PSD=plot_PSD_welch_conf(signal.filtfilt(bb, aa, EOF_time_series_T_LAB0m[0,:]), 1, '-', '-', '-', P_title, 'b', 'LAB 500m Density 1st EOF', 'lower right')
ff,PSD=plot_PSD_welch_conf(signal.filtfilt(bb, aa, EOF_time_series_T_LAB0m[1,:]), 1, '-', '-', '-', P_title, 'g', 'LAB 500m Density 2nd EOF', 'lower right')
ff,PSD=plot_PSD_welch_conf(signal.filtfilt(bb, aa, EOF_time_series_T_LAB0m[2,:]), 1, '-', '-', '-', P_title, 'darkorange', 'LAB 500m Density 3rd EOF', 'lower right')
ff,PSD=plot_PSD_welch_conf(signal.filtfilt(bb, aa, EOF_time_series_T_LAB0m[3,:]), 1, '-', '-', '-', P_title, 'r', 'LAB 500m Density 4th EOF', 'lower right')
#plt.ylim(1e-6,9.8) #;plt.xlim(1,500)
fig.savefig(dir_figs+str(GCM)+'_PSD_LAB_Temp500m_EOF_BW'+str(CutOff_T)+'.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

###############################################################################
###########           Lagged Correlations          ############################
###############################################################################

fig=plt.figure()
lag_cor(LAB_index_BWfilt, signal.filtfilt(bb, aa, EOF_time_series_T_LAB500m[0,:]), 40,'b','LAB 500m Temp 1st EOF')
lag_cor(LAB_index_BWfilt, signal.filtfilt(bb, aa, EOF_time_series_T_LAB500m[1,:]), 40,'g','LAB 500m Temp 2nd EOF')
lag_cor(LAB_index_BWfilt, signal.filtfilt(bb, aa, EOF_time_series_T_LAB500m[2,:]), 40,'darkorange','LAB 500m Temp 3rd EOF')
lag_cor(LAB_index_BWfilt, signal.filtfilt(bb, aa, EOF_time_series_T_LAB500m[3,:]), 40,'r','LAB 500m Temp 4th EOF')
plt.legend()
plt.show()
plt.title('Peak LAB Conv. lagged corr. with EOF of temperature at 500m - BW low-pass filt., '+str(CutOff_T)+'yr cut-off) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
plt.xlabel('Years lag', fontsize=18); plt.xticks(fontsize=18)
plt.ylabel('Correlation coefficient', fontsize=18); plt.yticks(fontsize=18)    
plt.show()
plt.legend(shadow=True, prop={'size': 15})
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_LABconvection_Temp500m_EOF_BW'+str(CutOff_T)+'_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')



fig=plt.figure()
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_T_LAB0m[0,:]), WS_index_BWfilt, 40,'y','WS index')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_T_LAB0m[0,:]), LAB_index_BWfilt, 40,'k','LAB index')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_T_LAB0m[0,:]), AMOC_max_50N_BWfilt, 40,'b','AMOC max 50N')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_T_LAB0m[0,:]), AMOC_max_30N_BWfilt, 40,'g','AMOC max 30N')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_T_LAB0m[0,:]), AMOC_max_30S_BWfilt, 40,'darkorange','AMOC max 30S')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_T_LAB0m[0,:]), AMOC_max_50S_BWfilt, 40,'r','AMOC max 50S')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_T_LAB0m[0,:]), ENSO_index_BWfilt, 40,'m','ENSO index')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_T_LAB0m[0,:]), NAO_index_BWfilt, 40,'c','NAO index')
plt.legend()
plt.show()
plt.title('1st EOF of Labrador Sea temperature at 500m lagged correlation (BW low-pass filt., '+str(CutOff_T)+'yr cut-off) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
plt.xlabel('Years lag', fontsize=18); plt.xticks(fontsize=18)
plt.ylabel('Correlation coefficient', fontsize=18); plt.yticks(fontsize=18)    
plt.show()
plt.legend(shadow=True, prop={'size': 15})
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_LAB_Temp500m_1st_EOF_BW'+str(CutOff_T)+'_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_T_LAB500m[1,:]), WS_index_BWfilt, 40,'y','WS index')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_T_LAB500m[1,:]), LAB_index_BWfilt, 40,'k','LAB index')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_T_LAB500m[1,:]), AMOC_max_50N_BWfilt, 40,'b','AMOC max 50N')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_T_LAB500m[1,:]), AMOC_max_30N_BWfilt, 40,'g','AMOC max 30N')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_T_LAB500m[1,:]), AMOC_max_30S_BWfilt, 40,'darkorange','AMOC max 30S')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_T_LAB500m[1,:]), AMOC_max_50S_BWfilt, 40,'r','AMOC max 50S')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_T_LAB500m[1,:]), ENSO_index_BWfilt, 40,'m','ENSO index')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_T_LAB500m[1,:]), NAO_index_BWfilt, 40,'c','NAO index')
plt.legend()
plt.show()
plt.title('2nd EOF of Labrador Sea temperature at 500m lagged correlation (BW low-pass filt., '+str(CutOff_T)+'yr cut-off) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
plt.xlabel('Years lag', fontsize=18); plt.xticks(fontsize=18)
plt.ylabel('Correlation coefficient', fontsize=18); plt.yticks(fontsize=18)    
plt.show()
plt.legend(shadow=True, prop={'size': 15})
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_LAB_Temp500m_2nd_EOF_BW'+str(CutOff_T)+'_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


fig=plt.figure()
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_T_LAB500m[2,:]), WS_index_BWfilt, 40,'y','WS index')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_T_LAB500m[2,:]), LAB_index_BWfilt, 40,'k','LAB index')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_T_LAB500m[2,:]), AMOC_max_50N_BWfilt, 40,'b','AMOC max 50N')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_T_LAB500m[2,:]), AMOC_max_30N_BWfilt, 40,'g','AMOC max 30N')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_T_LAB500m[2,:]), AMOC_max_30S_BWfilt, 40,'darkorange','AMOC max 30S')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_T_LAB500m[2,:]), AMOC_max_50S_BWfilt, 40,'r','AMOC max 50S')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_T_LAB500m[2,:]), ENSO_index_BWfilt, 40,'m','ENSO index')
lag_cor(signal.filtfilt(bb, aa, EOF_time_series_T_LAB500m[2,:]), NAO_index_BWfilt, 40,'c','NAO index')
plt.legend()
plt.show()
plt.title('3rd EOF of Labrador Sea temperature at 500m lagged correlation (BW low-pass filt., '+str(CutOff_T)+'yr cut-off) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
plt.xlabel('Years lag', fontsize=18); plt.xticks(fontsize=18)
plt.ylabel('Correlation coefficient', fontsize=18); plt.yticks(fontsize=18)    
plt.show()
plt.legend(shadow=True, prop={'size': 15})
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_LAB_Temp500m_3rd_EOF_BW'+str(CutOff_T)+'_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')




















