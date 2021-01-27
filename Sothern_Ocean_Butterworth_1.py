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
dir_figs = (dir_pwd + '/Figures1/') # Directory to save figures
dir_results = (dir_pwd + '/Results/') # Directory to save results

GCM = 'GFDL-ESM2G'
year_start=1
year_end=500

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
for key in my_shelf:
    globals()[key]=my_shelf[key]
my_shelf.close()
##############################################################################
#fs = 1000  # Sampling frequency
## Generate the time vector properly
#t = np.arange(1000) / fs
#fc = 35  # Cut-off frequency of the filter
#w = fc / (fs / 2) # Normalize the frequency
from scipy import stats, signal
from BehzadlibPlot import plot_PSD_welch_conf

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

AMOC_max_BWfilt = copy.deepcopy(AMOC_max)
for ii in range(AMOC_max.shape[1]):
    AMOC_max_BWfilt[:,ii] = signal.filtfilt(bb, aa, AMOC_max[:,ii])
    

###############################################################################
##################################################
#######   Power Spectral Density plots    ########
##################################################

P_title='Power Spectral Density - Weddell Sea Conv. Index Butterworth low-pass filtered ('+str(CutOff_T)+'yr cut-off, '+str(n_order)+'th order) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig=plt.figure(); ax = fig.add_axes([0.12,0.1,0.78,0.8]) # ax [left, bottom, width, height]
ff,PSD=plot_PSD_welch_conf(WS_index, 1, '-', '-', '-', P_title, 'k', 'WS Conv. Index Original', 'lower right')
plt.ylim(1e-6,9.8) #;plt.xlim(1,500)
ff,PSD=plot_PSD_welch_conf(WS_index_BWfilt, 1, '-', '-', '-', P_title, 'r', 'WS Conv. Index Butterworth filtered', 'lower right')
fig.savefig(dir_figs+str(GCM)+'_PSD_WS_BWfiltered'+str(CutOff_T)+'.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


P_title='Power Spectral Density - Labrador Sea Conv. Index Butterworth low-pass filtered ('+str(CutOff_T)+'yr cut-off, '+str(n_order)+'th order) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig=plt.figure(); ax = fig.add_axes([0.12,0.1,0.78,0.8]) # ax [left, bottom, width, height]
ff,PSD=plot_PSD_welch_conf(LAB_index, 1, '-', '-', '-', P_title, 'k', 'LAB Conv. Index Original', 'lower right')
plt.ylim(1e-6,9.8) #;plt.xlim(1,500)
ff,PSD=plot_PSD_welch_conf(LAB_index_BWfilt, 1, '-', '-', '-', P_title, 'r', 'LAB Conv. Index Butterworth filtered', 'lower right')
fig.savefig(dir_figs+str(GCM)+'_PSD_LAB_BWfiltered'+str(CutOff_T)+'.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

###############################################################################
WS_index_BWfilt_norm=(WS_index_BWfilt-np.nanmean(WS_index_BWfilt))/np.std(WS_index_BWfilt) # Normalized Convection Index
LAB_index_BWfilt_norm=(LAB_index_BWfilt-np.nanmean(LAB_index_BWfilt))/np.std(LAB_index_BWfilt) # Normalized Convection Index


P_Var_x1=P_Var_x2=years
P_Var_y1=copy.deepcopy(WS_index_norm)
P_Var_y2=copy.deepcopy(WS_index_BWfilt_norm)
P_xlable='Years' ;  P_ylable='Convection index (°C)'
P_title='Weddell Sea conv. index - Butterworth low-pass filt. '+str(CutOff_T)+'yr cut-off, '+str(n_order)+'th order (normalized, NOT-smoothed) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plotline_2var(P_Var_x1, P_Var_y1, P_Var_x2, P_Var_y2, P_xlable, P_ylable, P_title, 'b', 'r', 'Conv. Index Original', 'Conv. Index Butterworth filtered', 'best', '-')
for l in fig.gca().lines:
    l.set_alpha(0.7)
l = plt.axhline(y=-0.5, color='k', linestyle ='--')
l = plt.axhline(y=0, color='k')
l = plt.axhline(y=0.5, color='k', linestyle ='--')
plt.text(10, 1.1, 'Non-Conv.', horizontalalignment='center', verticalalignment='center', fontsize=20, color='b')
plt.text(0, -1.1, 'Conv.', horizontalalignment='center', verticalalignment='center', fontsize=20, color='b')
fig.savefig(dir_figs+str(GCM)+'_convec_index_WS_BW'+str(CutOff_T)+'.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

P_Var_x1=P_Var_x2=years
P_Var_y1=copy.deepcopy(LAB_index_norm)
P_Var_y2=copy.deepcopy(LAB_index_BWfilt_norm)
P_xlable='Years' ;  P_ylable='Convection index (°C)'
P_title='Labrador Sea conv. index - Butterworth low-pass filt. '+str(CutOff_T)+'yr cut-off, '+str(n_order)+'th order (normalized, NOT-smoothed) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plotline_2var(P_Var_x1, P_Var_y1, P_Var_x2, P_Var_y2, P_xlable, P_ylable, P_title, 'b', 'r', 'Conv. Index Original', 'Conv. Index Butterworth filtered', 'best', '-')
for l in fig.gca().lines:
    l.set_alpha(0.7)
l = plt.axhline(y=-0.5, color='k', linestyle ='--')
l = plt.axhline(y=0, color='k')
l = plt.axhline(y=0.5, color='k', linestyle ='--')
plt.text(10, 1.1, 'Non-Conv.', horizontalalignment='center', verticalalignment='center', fontsize=20, color='b')
plt.text(0, -1.1, 'Conv.', horizontalalignment='center', verticalalignment='center', fontsize=20, color='b')
fig.savefig(dir_figs+str(GCM)+'_convec_index_LAB_BW'+str(CutOff_T)+'.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

###############################################################################
lag_time=40

fig=plt.figure()
data=[]
for i in range(len(AMOC_max_BWfilt[0][20:161])):
    stream=signal.filtfilt(bb, aa, AMOC_max_BWfilt[:,i+lag_time])
    r=lag_cor_data(WS_index_BWfilt,stream,lag_time)
    data.append(r)
data=np.asarray(data)
x=np.linspace(-lag_time,lag_time+1, 2*lag_time)
y=np.linspace(-70,70+1, 2*70+1)

x,y=np.meshgrid(x,y)
cmap_limit=np.nanmax(np.abs(data))
levels=np.linspace(-cmap_limit,cmap_limit,200)
im=plt.contourf(x,y,data,levels,latlon=True, cmap=plt.cm.jet)
cb = plt.colorbar(im); cb.ax.tick_params(labelsize=18) 
plt.xlabel('Year lag', fontsize=18); plt.xticks(fontsize=18)
plt.ylabel('Latitude', fontsize=18); plt.yticks(fontsize=18)
plt.title('WS peak convection vs. AMOC - Butterworth filtered ('+str(CutOff_T)+'yr cut-off, '+str(n_order)+'th order) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_convec_AMOC_max_WS_BW'+str(CutOff_T)+'.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


fig=plt.figure()
data=[]
for i in range(len(AMOC_max_BWfilt[0][20:161])):
    stream=signal.filtfilt(bb, aa, AMOC_max_BWfilt[:,i+lag_time])
    r=lag_cor_data(LAB_index_BWfilt,stream,lag_time)
    data.append(r)
data=np.asarray(data)
x=np.linspace(-lag_time,lag_time+1, 2*lag_time)
y=np.linspace(-70,70+1, 2*70+1)

x,y=np.meshgrid(x,y)
cmap_limit=np.nanmax(np.abs(data))
levels=np.linspace(-cmap_limit,cmap_limit,200)
im=plt.contourf(x,y,data,levels,latlon=True, cmap=plt.cm.jet)
cb = plt.colorbar(im); cb.ax.tick_params(labelsize=18)
plt.xlabel('Year lag', fontsize=18); plt.xticks(fontsize=18)
plt.ylabel('Latitude', fontsize=18); plt.yticks(fontsize=18)
plt.title('LAB peak convection vs. AMOC - Butterworth filtered ('+str(CutOff_T)+'yr cut-off, '+str(n_order)+'th order) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_convec_AMOC_max_LAB_BW'+str(CutOff_T)+'.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

###########################
#fig=plt.figure()
#data=[]
#for i in range(len(AMOC_max[0][20:161])):
#    stream=AMOC_max[:,i+lag_time]
#    r=lag_cor_data(WS_index,stream,lag_time)
#    data.append(r)
#data=np.asarray(data)
#x=np.linspace(-lag_time,lag_time+1, 2*lag_time)
#y=np.linspace(-70,70+1, 2*70+1)
#
#x,y=np.meshgrid(x,y)
#cmap_limit=np.nanmax(np.abs(data))
#levels=np.linspace(-cmap_limit,cmap_limit,200)
#im=plt.contourf(x,y,data,levels,latlon=True, cmap=plt.cm.jet)
#cb = plt.colorbar(im); cb.ax.tick_params(labelsize=18) 
#plt.xlabel('Year lag', fontsize=18); plt.xticks(fontsize=18)
#plt.ylabel('Latitude', fontsize=18); plt.yticks(fontsize=18)
#plt.title('WS peak convection vs. AMOC - Original Time Series, No Smoothing At All - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
#mng = plt.get_current_fig_manager()
#mng.window.showMaximized() # Maximizes the plot window to save figures in full
#fig.savefig(dir_figs+str(GCM)+'_lagcor_convec_AMOC_max_WS_NOsmoothingAtAll.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
#
#
#fig=plt.figure()
#data=[]
#for i in range(len(AMOC_max[0][20:161])):
#    stream=AMOC_max[:,i+lag_time]
#    r=lag_cor_data(LAB_index,stream,lag_time)
#    data.append(r)
#data=np.asarray(data)
#x=np.linspace(-lag_time,lag_time+1, 2*lag_time)
#y=np.linspace(-70,70+1, 2*70+1)
#
#x,y=np.meshgrid(x,y)
#cmap_limit=np.nanmax(np.abs(data))
#levels=np.linspace(-cmap_limit,cmap_limit,200)
#im=plt.contourf(x,y,data,levels,latlon=True, cmap=plt.cm.jet)
#cb = plt.colorbar(im); cb.ax.tick_params(labelsize=18) 
#plt.xlabel('Year lag', fontsize=18); plt.xticks(fontsize=18)
#plt.ylabel('Latitude', fontsize=18); plt.yticks(fontsize=18)
#plt.title('LAB peak convection vs. AMOC - Original Time Series, No Smoothing At All - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
#mng = plt.get_current_fig_manager()
#mng.window.showMaximized() # Maximizes the plot window to save figures in full
#fig.savefig(dir_figs+str(GCM)+'_lagcor_convec_AMOC_max_LAB_NOsmoothingAtAll.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
###############################################################################

fig=plt.figure()
lag_cor(WS_index_BWfilt, ENSO_index_BWfilt, 40,'m','ENSO index')
lag_cor(WS_index_BWfilt, NAO_index_BWfilt, 40,'c','NAO index')
lag_cor(WS_index_BWfilt, AMOC_max_50N_BWfilt, 40,'b','AMOC max 50N')
lag_cor(WS_index_BWfilt, AMOC_max_30N_BWfilt, 40,'g','AMOC max 30N')
lag_cor(WS_index_BWfilt, AMOC_max_30S_BWfilt, 40,'darkorange','AMOC max 30S')
lag_cor(WS_index_BWfilt, AMOC_max_50S_BWfilt, 40,'r','AMOC max 50S')
plt.legend()
plt.show()
plt.title('Peak WS Conv. lagged corr. Butterworth low-pass filt. ('+str(CutOff_T)+'yr cut-off, '+str(n_order)+'th order) [AMOC=max streamfunc.] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
plt.xlabel('Years lag', fontsize=18); plt.xticks(fontsize=18)
plt.ylabel('Correlation coefficient', fontsize=18) ; plt.yticks(fontsize=18)   
plt.show()
plt.legend(shadow=True, prop={'size': 15})
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_WSconvection_AMOC_max_BW'+str(CutOff_T)+'_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


fig=plt.figure()
lag_cor(LAB_index_BWfilt, ENSO_index_BWfilt, 40,'m','ENSO index')
lag_cor(LAB_index_BWfilt, NAO_index_BWfilt, 40,'c','NAO index')
lag_cor(LAB_index_BWfilt, AMOC_max_50N_BWfilt, 40,'b','AMOC max 50N')
lag_cor(LAB_index_BWfilt, AMOC_max_30N_BWfilt, 40,'g','AMOC max 30N')
lag_cor(LAB_index_BWfilt, AMOC_max_30S_BWfilt, 40,'darkorange','AMOC max 30S')
lag_cor(LAB_index_BWfilt, AMOC_max_50S_BWfilt, 40,'r','AMOC max 50S')
plt.legend()
plt.show()
plt.title('Peak LAB Conv. lagged corr. Butterworth low-pass filt. ('+str(CutOff_T)+'yr cut-off, '+str(n_order)+'th order) [AMOC=max streamfunc.] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
plt.xlabel('Years lag', fontsize=18); plt.xticks(fontsize=18)
plt.ylabel('Correlation coefficient', fontsize=18); plt.yticks(fontsize=18)    
plt.show()
plt.legend(shadow=True, prop={'size': 15})
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_LABconvection_AMOC_max_BW'+str(CutOff_T)+'_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

########################################################
########    LAB Gyre Strength Time Series    ###########
Stream_function_Barotropic_NAtl_VyDx_allyears=np.load('Stream_function_Barotropic_NAtl_VyDx_allyears_GFDL-ESM2G_500yr.npy')
Stream_function_Barotropic_NAtl_VyDx_allyears_LAB_B1_max = np.nanmax (Stream_function_Barotropic_NAtl_VyDx_allyears [:,140:156,295:320] , axis =2)
LAB_index_gyre = np.nanmax (Stream_function_Barotropic_NAtl_VyDx_allyears_LAB_B1_max , axis =1)
LAB_index_gyre_norm=(LAB_index_gyre-np.nanmean(LAB_index_gyre))/np.std(LAB_index_gyre) # Normalized Convection Index

LAB_index_gyre_BWfilt = signal.filtfilt(bb, aa, LAB_index_gyre)
LAB_index_gyre_BWfilt_norm=(LAB_index_gyre_BWfilt-np.nanmean(LAB_index_gyre_BWfilt))/np.std(LAB_index_gyre_BWfilt) # Normalized Convection Index


P_title='Power Spectral Density - Labrador Sea GYRE Index Butterworth low-pass filtered ('+str(CutOff_T)+'yr cut-off, '+str(n_order)+'th order) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig=plt.figure(); ax = fig.add_axes([0.12,0.1,0.78,0.8]) # ax [left, bottom, width, height]
ff,PSD=plot_PSD_welch_conf(LAB_index_gyre, 1, '-', '-', '-', P_title, 'k', 'LAB Gyre Index Original', 'lower right')
plt.ylim(1e-2,980) #;plt.xlim(1,500)
ff,PSD=plot_PSD_welch_conf(LAB_index_gyre_BWfilt, 1, '-', '-', '-', P_title, 'r', 'LAB Gyre Index Butterworth filtered', 'lower right')
fig.savefig(dir_figs+str(GCM)+'_PSD_LAB_gyre_BWfiltered'+str(CutOff_T)+'.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


P_Var_x1=P_Var_x2=years
P_Var_y1=copy.deepcopy(LAB_index_gyre_norm)
P_Var_y2=copy.deepcopy(LAB_index_gyre_BWfilt_norm)
P_xlable='Years' ;  P_ylable='Convection/Gyre index - Barometric Streamfunction Maxima (Sv)'
P_title='Labrador Sea GYRE index - Butterworth low-pass filt. '+str(CutOff_T)+'yr cut-off, '+str(n_order)+'th order (normalized, NOT-smoothed) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plotline_2var(P_Var_x1, P_Var_y1, P_Var_x2, P_Var_y2, P_xlable, P_ylable, P_title, 'b', 'r', 'Gyre Index Original', 'Gyre Index Butterworth filtered', 'best', '-')
for l in fig.gca().lines:
    l.set_alpha(0.7)
l = plt.axhline(y=-0.5, color='k', linestyle ='--')
l = plt.axhline(y=0, color='k')
l = plt.axhline(y=0.5, color='k', linestyle ='--')
plt.text(10, 1.1, 'Non-Conv.', horizontalalignment='center', verticalalignment='center', fontsize=20, color='b')
plt.text(0, -1.1, 'Conv.', horizontalalignment='center', verticalalignment='center', fontsize=20, color='b')
fig.savefig(dir_figs+str(GCM)+'_convec_index_LAB_gyre_BW'+str(CutOff_T)+'.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


fig=plt.figure()
data=[]
for i in range(len(AMOC_max[0][20:161])):
    stream=signal.filtfilt(bb, aa, AMOC_max[:,i+lag_time])
    r=lag_cor_data(LAB_index_gyre_BWfilt,stream,lag_time)
    data.append(r)
data=np.asarray(data)
x=np.linspace(-lag_time,lag_time+1, 2*lag_time)
y=np.linspace(-70,70+1, 2*70+1)

x,y=np.meshgrid(x,y)
cmap_limit=np.nanmax(np.abs(data))
levels=np.linspace(-cmap_limit,cmap_limit,200)
im=plt.contourf(x,y,data,levels,latlon=True, cmap=plt.cm.jet)
cb = plt.colorbar(im); cb.ax.tick_params(labelsize=18)
plt.xlabel('Year lag', fontsize=18); plt.xticks(fontsize=18)
plt.ylabel('Latitude', fontsize=18); plt.yticks(fontsize=18)
plt.title('LAB peak GYRE strength vs. AMOC - Butterworth filtered ('+str(CutOff_T)+'yr cut-off, '+str(n_order)+'th order) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_convec_AMOC_max_LAB_gyre_BW'+str(CutOff_T)+'.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


fig=plt.figure()
lag_cor(LAB_index_gyre_BWfilt, ENSO_index_BWfilt, 40,'m','ENSO index')
lag_cor(LAB_index_gyre_BWfilt, NAO_index_BWfilt, 40,'c','NAO index')
lag_cor(LAB_index_gyre_BWfilt, AMOC_max_50N_BWfilt, 40,'b','AMOC max 50N')
lag_cor(LAB_index_gyre_BWfilt, AMOC_max_30N_BWfilt, 40,'g','AMOC max 30N')
lag_cor(LAB_index_gyre_BWfilt, AMOC_max_30S_BWfilt, 40,'darkorange','AMOC max 30S')
lag_cor(LAB_index_gyre_BWfilt, AMOC_max_50S_BWfilt, 40,'r','AMOC max 50S')
plt.legend()
plt.show()
plt.title('Peak LAB GYRE strength lagged corr. Butterworth low-pass filt. ('+str(CutOff_T)+'yr cut-off, '+str(n_order)+'th order) [AMOC=max streamfunc.] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
plt.xlabel('Years lag', fontsize=18); plt.xticks(fontsize=18)
plt.ylabel('Correlation coefficient', fontsize=18); plt.yticks(fontsize=18)    
plt.show()
plt.legend(shadow=True, prop={'size': 15})
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_LAB_gyre_AMOC_max_BW'+str(CutOff_T)+'_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')






######################################################
##### Plots with NO SMOOTHING / FILTERING at all #####

fig=plt.figure()
data=[]
for i in range(len(AMOC_max[0][20:161])):
    stream=AMOC_max[:,i+lag_time]
    r=lag_cor_data(WS_index,stream,lag_time)
    data.append(r)
data=np.asarray(data)
x=np.linspace(-lag_time,lag_time+1, 2*lag_time)
y=np.linspace(-70,70+1, 2*70+1)

x,y=np.meshgrid(x,y)
cmap_limit=np.nanmax(np.abs(data))
levels=np.linspace(-cmap_limit,cmap_limit,200)
im=plt.contourf(x,y,data,levels,latlon=True, cmap=plt.cm.jet)
cb = plt.colorbar(im); cb.ax.tick_params(labelsize=18) 
plt.xlabel('Year lag', fontsize=18); plt.xticks(fontsize=18)
plt.ylabel('Latitude', fontsize=18); plt.yticks(fontsize=18)
plt.title('WS peak convection vs. AMOC - Original Time Series, No Smoothing At All - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_convec_AMOC_max_WS_NOsmoothingAtAll.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


fig=plt.figure()
data=[]
for i in range(len(AMOC_max[0][20:161])):
    stream=AMOC_max[:,i+lag_time]
    r=lag_cor_data(LAB_index,stream,lag_time)
    data.append(r)
data=np.asarray(data)
x=np.linspace(-lag_time,lag_time+1, 2*lag_time)
y=np.linspace(-70,70+1, 2*70+1)

x,y=np.meshgrid(x,y)
cmap_limit=np.nanmax(np.abs(data))
levels=np.linspace(-cmap_limit,cmap_limit,200)
im=plt.contourf(x,y,data,levels,latlon=True, cmap=plt.cm.jet)
cb = plt.colorbar(im); cb.ax.tick_params(labelsize=18) 
plt.xlabel('Year lag', fontsize=18); plt.xticks(fontsize=18)
plt.ylabel('Latitude', fontsize=18); plt.yticks(fontsize=18)
plt.title('LAB peak convection vs. AMOC - Original Time Series, No Smoothing At All - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_convec_AMOC_max_LAB_NOsmoothingAtAll.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

fig=plt.figure()
lag_cor(WS_index, ENSO_index, 40,'m','ENSO index')
lag_cor(WS_index, NAO_index, 40,'c','NAO index')
lag_cor(WS_index, AMOC_max_50N, 40,'b','AMOC max 50N')
lag_cor(WS_index, AMOC_max_30N, 40,'g','AMOC max 30N')
lag_cor(WS_index, AMOC_max_30S, 40,'darkorange','AMOC max 30S')
lag_cor(WS_index, AMOC_max_50S, 40,'r','AMOC max 50S')
plt.legend()
plt.show()
plt.title('Peak WS Conv. lagged corr. - NO smoothing/filtering at all [AMOC=max streamfunc.] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
plt.xlabel('Years lag', fontsize=18); plt.xticks(fontsize=18)
plt.ylabel('Correlation coefficient', fontsize=18) ; plt.yticks(fontsize=18)   
plt.show()
plt.legend(shadow=True, prop={'size': 15})
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_WSconvection_AMOC_max_NOsmoothingAtAll.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


fig=plt.figure()
lag_cor(LAB_index, ENSO_index, 40,'m','ENSO index')
lag_cor(LAB_index, NAO_index, 40,'c','NAO index')
lag_cor(LAB_index, AMOC_max_50N, 40,'b','AMOC max 50N')
lag_cor(LAB_index, AMOC_max_30N, 40,'g','AMOC max 30N')
lag_cor(LAB_index, AMOC_max_30S, 40,'darkorange','AMOC max 30S')
lag_cor(LAB_index, AMOC_max_50S, 40,'r','AMOC max 50S')
plt.legend()
plt.show()
plt.title('Peak LAB Conv. lagged corr. - NO smoothing/filtering at all [AMOC=max streamfunc.] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
plt.xlabel('Years lag', fontsize=18); plt.xticks(fontsize=18)
plt.ylabel('Correlation coefficient', fontsize=18) ; plt.yticks(fontsize=18)   
plt.show()
plt.legend(shadow=True, prop={'size': 15})
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_LABconvection_AMOC_max_NOsmoothingAtAll.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


fig=plt.figure()
data=[]
for i in range(len(AMOC_max[0][20:161])):
    stream=signal.filtfilt(bb, aa, AMOC_max[:,i+lag_time])
    r=lag_cor_data(LAB_index_gyre,stream,lag_time)
    data.append(r)
data=np.asarray(data)
x=np.linspace(-lag_time,lag_time+1, 2*lag_time)
y=np.linspace(-70,70+1, 2*70+1)

x,y=np.meshgrid(x,y)
cmap_limit=np.nanmax(np.abs(data))
levels=np.linspace(-cmap_limit,cmap_limit,200)
im=plt.contourf(x,y,data,levels,latlon=True, cmap=plt.cm.jet)
cb = plt.colorbar(im); cb.ax.tick_params(labelsize=18)
plt.xlabel('Year lag', fontsize=18); plt.xticks(fontsize=18)
plt.ylabel('Latitude', fontsize=18); plt.yticks(fontsize=18)
plt.title('LAB peak GYRE strength vs. AMOC - Original Time Series, No Smoothing At All - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_convec_AMOC_max_LAB_gyre_NOsmoothingAtAll.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


fig=plt.figure()
lag_cor(LAB_index_gyre, ENSO_index, 40,'m','ENSO index')
lag_cor(LAB_index_gyre, NAO_index, 40,'c','NAO index')
lag_cor(LAB_index_gyre, AMOC_max_50N, 40,'b','AMOC max 50N')
lag_cor(LAB_index_gyre, AMOC_max_30N, 40,'g','AMOC max 30N')
lag_cor(LAB_index_gyre, AMOC_max_30S, 40,'darkorange','AMOC max 30S')
lag_cor(LAB_index_gyre, AMOC_max_50S, 40,'r','AMOC max 50S')
plt.legend()
plt.show()
plt.title('Peak LAB GYRE strength lagged corr. - NO smoothing/filtering at all [AMOC=max streamfunc.] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
plt.xlabel('Years lag', fontsize=18); plt.xticks(fontsize=18)
plt.ylabel('Correlation coefficient', fontsize=18); plt.yticks(fontsize=18)    
plt.show()
plt.legend(shadow=True, prop={'size': 15})
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_lagcor_LAB_gyre_AMOC_max_NOsmoothingAtAll.png', format='png', dpi=300, transparent=True, bbox_inches='tight')






























