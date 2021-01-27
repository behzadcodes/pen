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

##################################################
#%%####   Power Spectral Density plots    ########
##################################################
from scipy import stats

#P_Var=WS_index ; P_legend='WS Convection Index' ; P_color='r'; P_c_probability=0.95; P_rho='yes'; P_smooth=5
#ff,PSD=plot_PSD_welch_conf(P_Var, 1, P_c_probability, P_rho, P_smooth, P_title, P_color, P_legend, 'lower right')

slope, intercept, r_value, p_value, std_err = stats.linregress(WS_index,LAB_index)
P_title='Power Spectral Density - Weddell Sea Convection Index - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig=plt.figure(); ax = fig.add_axes([0.12,0.1,0.78,0.8]) # ax [left, bottom, width, height]
ff,PSD=plot_PSD_welch_conf(WS_index, 1, '-', '-', '-', P_title, 'k', 'WS Convection Index', 'lower right')
ff,PSD=plot_PSD_welch_conf(LAB_index, 1, '-', '-', '-', P_title, 'r', 'LAB Convection Index', 'lower right')
plt.ylim(1e-6,9.8) #;plt.xlim(1,500)
plt.text(0.5, 0.97, 'Corrolation (R)= '+str(float("{0:.4f}".format(r_value))), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=20, color='blue')
#fig.savefig(dir_figs+str(GCM)+'_PSD_WS_LAB.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

slope, intercept, r_value, p_value, std_err = stats.linregress(WS_index_norm,LAB_index_norm)
P_title='Power Spectral Density - Weddell Sea Convection Index - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig=plt.figure(); ax = fig.add_axes([0.12,0.1,0.78,0.8]) # ax [left, bottom, width, height]
ff,PSD=plot_PSD_welch_conf(WS_index_norm, 1, '-', '-', '-', P_title, 'k', 'WS Convection Index Normalized', 'lower right')
ff,PSD=plot_PSD_welch_conf(LAB_index_norm, 1, '-', '-', '-', P_title, 'r', 'LAB Convection Index Normalized', 'lower right')
plt.ylim(1e-6,100) #;plt.xlim(1,500)
plt.text(0.5, 0.97, 'Corrolation (R)= '+str(float("{0:.4f}".format(r_value))), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=20, color='blue')
#fig.savefig(dir_figs+str(GCM)+'_PSD_WS_LAB_norm.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

slope, intercept, r_value, p_value, std_err = stats.linregress(WS_index,np.nanmean(South_Westerlies, axis=1))
P_title='Power Spectral Density - Weddell Sea Convection Index - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig=plt.figure(); ax = fig.add_axes([0.12,0.1,0.78,0.8]) # ax [left, bottom, width, height]
ff,PSD=plot_PSD_welch_conf(WS_index, 1, '-', '-', '-', P_title, 'k', 'WS Convection Index', 'lower right')
ff,PSD=plot_PSD_welch_conf(np.nanmean(South_Westerlies, axis=1), 1, '-', '-', '-', P_title, 'r', 'South Westerlies (S45-S55)', 'lower right')
plt.ylim(1e-6,9.8) #;plt.xlim(1,500)
plt.text(0.5, 0.97, 'Corrolation (R)= '+str(float("{0:.4f}".format(r_value))), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=20, color='blue')
#fig.savefig(dir_figs+str(GCM)+'_PSD_WS_SouthWesterlies.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

slope, intercept, r_value, p_value, std_err = stats.linregress(LAB_index,np.nanmean(South_Westerlies, axis=1))
P_title='Power Spectral Density - Labrador Sea Convection Index - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig=plt.figure(); ax = fig.add_axes([0.12,0.1,0.78,0.8]) # ax [left, bottom, width, height]
ff,PSD=plot_PSD_welch_conf(LAB_index, 1, '-', '-', '-', P_title, 'k', 'LAB Convection Index', 'lower right')
ff,PSD=plot_PSD_welch_conf(np.nanmean(South_Westerlies, axis=1), 1, '-', '-', '-', P_title, 'r', 'South Westerlies (S45-S55)', 'lower right')
plt.ylim(1e-6,9.8) #;plt.xlim(1,500)
plt.text(0.5, 0.97, 'Corrolation (R)= '+str(float("{0:.4f}".format(r_value))), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=20, color='blue')
#fig.savefig(dir_figs+str(GCM)+'_PSD_LAB_SouthWesterlies.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

slope, intercept, r_value, p_value, std_err = stats.linregress(WS_index,SMOC_min_50S)
P_title='Power Spectral Density - Weddell Sea Convection Index - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig=plt.figure(); ax = fig.add_axes([0.12,0.1,0.78,0.8]) # ax [left, bottom, width, height]
ff,PSD=plot_PSD_welch_conf(WS_index, 1, '-', '-', '-', P_title, 'k', 'WS Convection Index', 'lower right')
ff,PSD=plot_PSD_welch_conf(SMOC_min_50S, 1, '-', '-', '-', P_title, 'r', 'SMOC min 50S', 'lower right')
plt.ylim(1e-5,1e3) #;plt.xlim(1,500)
plt.text(0.5, 0.97, 'Corrolation (R)= '+str(float("{0:.4f}".format(r_value))), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=20, color='blue')
#fig.savefig(dir_figs+str(GCM)+'_PSD_WS_SMOC50S.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

slope, intercept, r_value, p_value, std_err = stats.linregress(WS_index,AMOC_max_30S)
P_title='Power Spectral Density - Weddell Sea Convection Index - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig=plt.figure(); ax = fig.add_axes([0.12,0.1,0.78,0.8]) # ax [left, bottom, width, height]
ff,PSD=plot_PSD_welch_conf(WS_index, 1, '-', '-', '-', P_title, 'k', 'WS Convection Index', 'lower right')
ff,PSD=plot_PSD_welch_conf(AMOC_max_30S, 1, '-', '-', '-', P_title, 'r', 'AMOC max 30S', 'lower right')
plt.ylim(1e-5,1e3) #;plt.xlim(1,500)
plt.text(0.5, 0.97, 'Corrolation (R)= '+str(float("{0:.4f}".format(r_value))), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=20, color='blue')
#fig.savefig(dir_figs+str(GCM)+'_PSD_WS_AMOC30S.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

slope, intercept, r_value, p_value, std_err = stats.linregress(WS_index,AMOC_max_50S)
P_title='Power Spectral Density - Weddell Sea Convection Index - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig=plt.figure(); ax = fig.add_axes([0.12,0.1,0.78,0.8]) # ax [left, bottom, width, height]
ff,PSD=plot_PSD_welch_conf(WS_index, 1, '-', '-', '-', P_title, 'k', 'WS Convection Index', 'lower right')
ff,PSD=plot_PSD_welch_conf(AMOC_max_50S, 1, '-', '-', '-', P_title, 'r', 'AMOC max 50S', 'lower right')
plt.ylim(1e-5,1e3) #;plt.xlim(1,500)
plt.text(0.5, 0.97, 'Corrolation (R)= '+str(float("{0:.4f}".format(r_value))), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=20, color='blue')
#fig.savefig(dir_figs+str(GCM)+'_PSD_WS_AMOC50S.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

slope, intercept, r_value, p_value, std_err = stats.linregress(WS_index,AMOC_max_30N)
P_title='Power Spectral Density - Weddell Sea Convection Index - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig=plt.figure(); ax = fig.add_axes([0.12,0.1,0.78,0.8]) # ax [left, bottom, width, height]
ff,PSD=plot_PSD_welch_conf(WS_index, 1, '-', '-', '-', P_title, 'k', 'WS Convection Index', 'lower right')
ff,PSD=plot_PSD_welch_conf(AMOC_max_30N, 1, '-', '-', '-', P_title, 'r', 'AMOC max 30N', 'lower right')
plt.ylim(1e-5,1e3) #;plt.xlim(1,500)
plt.text(0.5, 0.97, 'Corrolation (R)= '+str(float("{0:.4f}".format(r_value))), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=20, color='blue')
#fig.savefig(dir_figs+str(GCM)+'_PSD_WS_AMOC30N.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

slope, intercept, r_value, p_value, std_err = stats.linregress(WS_index,AMOC_max_50N)
P_title='Power Spectral Density - Weddell Sea Convection Index - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig=plt.figure(); ax = fig.add_axes([0.12,0.1,0.78,0.8]) # ax [left, bottom, width, height]
ff,PSD=plot_PSD_welch_conf(WS_index, 1, '-', '-', '-', P_title, 'k', 'WS Convection Index', 'lower right')
ff,PSD=plot_PSD_welch_conf(AMOC_max_50N, 1, '-', '-', '-', P_title, 'r', 'AMOC max 50N', 'lower right')
plt.ylim(1e-5,1e3) #;plt.xlim(1,500)
plt.text(0.5, 0.97, 'Corrolation (R)= '+str(float("{0:.4f}".format(r_value))), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=20, color='blue')
#fig.savefig(dir_figs+str(GCM)+'_PSD_WS_AMOC50N.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

slope, intercept, r_value, p_value, std_err = stats.linregress(LAB_index,SMOC_min_50S)
P_title='Power Spectral Density - Labrador Sea Convection Index - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig=plt.figure(); ax = fig.add_axes([0.12,0.1,0.78,0.8]) # ax [left, bottom, width, height]
ff,PSD=plot_PSD_welch_conf(LAB_index, 1, '-', '-', '-', P_title, 'k', 'LAB Convection Index', 'lower right')
ff,PSD=plot_PSD_welch_conf(SMOC_min_50S, 1, '-', '-', '-', P_title, 'r', 'SMOC min 50S', 'lower right')
plt.ylim(1e-5,1e3) #;plt.xlim(1,500)
plt.text(0.5, 0.97, 'Corrolation (R)= '+str(float("{0:.4f}".format(r_value))), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=20, color='blue')
#fig.savefig(dir_figs+str(GCM)+'_PSD_LAB_SMOC50S.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

slope, intercept, r_value, p_value, std_err = stats.linregress(LAB_index,AMOC_max_30S)
P_title='Power Spectral Density - Labrador Sea Convection Index - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig=plt.figure(); ax = fig.add_axes([0.12,0.1,0.78,0.8]) # ax [left, bottom, width, height]
ff,PSD=plot_PSD_welch_conf(LAB_index, 1, '-', '-', '-', P_title, 'k', 'LAB Convection Index', 'lower right')
ff,PSD=plot_PSD_welch_conf(AMOC_max_30S, 1, '-', '-', '-', P_title, 'r', 'AMOC max 30S', 'lower right')
plt.ylim(1e-5,1e3) #;plt.xlim(1,500)
plt.text(0.5, 0.97, 'Corrolation (R)= '+str(float("{0:.4f}".format(r_value))), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=20, color='blue')
#fig.savefig(dir_figs+str(GCM)+'_PSD_LAB_AMOC30S.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

slope, intercept, r_value, p_value, std_err = stats.linregress(LAB_index,AMOC_max_50S)
P_title='Power Spectral Density - Labrador Sea Convection Index - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig=plt.figure(); ax = fig.add_axes([0.12,0.1,0.78,0.8]) # ax [left, bottom, width, height]
ff,PSD=plot_PSD_welch_conf(LAB_index, 1, '-', '-', '-', P_title, 'k', 'LAB Convection Index', 'lower right')
ff,PSD=plot_PSD_welch_conf(AMOC_max_50S, 1, '-', '-', '-', P_title, 'r', 'AMOC max 50S', 'lower right')
plt.ylim(1e-5,1e3) #;plt.xlim(1,500)
plt.text(0.5, 0.97, 'Corrolation (R)= '+str(float("{0:.4f}".format(r_value))), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=20, color='blue')
#fig.savefig(dir_figs+str(GCM)+'_PSD_LAB_AMOC50S.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

slope, intercept, r_value, p_value, std_err = stats.linregress(LAB_index,AMOC_max_30N)
P_title='Power Spectral Density - Labrador Sea Convection Index - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig=plt.figure(); ax = fig.add_axes([0.12,0.1,0.78,0.8]) # ax [left, bottom, width, height]
ff,PSD=plot_PSD_welch_conf(LAB_index, 1, '-', '-', '-', P_title, 'k', 'LAB Convection Index', 'lower right')
ff,PSD=plot_PSD_welch_conf(AMOC_max_30N, 1, '-', '-', '-', P_title, 'r', 'AMOC max 30N', 'lower right')
plt.ylim(1e-5,1e3) #;plt.xlim(1,500)
plt.text(0.5, 0.97, 'Corrolation (R)= '+str(float("{0:.4f}".format(r_value))), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=20, color='blue')
#fig.savefig(dir_figs+str(GCM)+'_PSD_LAB_AMOC30N.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

slope, intercept, r_value, p_value, std_err = stats.linregress(LAB_index,AMOC_max_50N)
P_title='Power Spectral Density - Labrador Sea Convection Index - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig=plt.figure(); ax = fig.add_axes([0.12,0.1,0.78,0.8]) # ax [left, bottom, width, height]
ff,PSD=plot_PSD_welch_conf(LAB_index, 1, '-', '-', '-', P_title, 'k', 'LAB Convection Index', 'lower right')
ff,PSD=plot_PSD_welch_conf(AMOC_max_50N, 1, '-', '-', '-', P_title, 'r', 'AMOC max 50N', 'lower right')
plt.ylim(1e-5,1e3) #;plt.xlim(1,500)
plt.text(0.5, 0.97, 'Corrolation (R)= '+str(float("{0:.4f}".format(r_value))), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=20, color='blue')
#fig.savefig(dir_figs+str(GCM)+'_PSD_LAB_AMOC50N.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

###############################################################################
###############################################################################

P_Var=WS_index ; P_legend='WS Convection Index' ; P_color='r'; P_c_probability=0.95; P_rho='yes'; P_smooth=5
P_title='Power Spectral Density - WS Convection Index - yrs '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig=plt.figure()
ff,PSD=plot_PSD_welch_conf(P_Var, 1, P_c_probability, P_rho, P_smooth, P_title, P_color, P_legend, 'lower right')
fig.savefig(dir_figs+str(GCM)+'_PSD_conf_WS_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


P_Var=LAB_index ; P_legend='LAB Convection Index' ; P_color='r'; P_c_probability=0.95; P_rho='yes'; P_smooth=5
P_title='Power Spectral Density - LAB Convection Index - yrs '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig=plt.figure()
ff,PSD=plot_PSD_welch_conf(P_Var, 1, P_c_probability, P_rho, P_smooth, P_title, P_color, P_legend, 'lower right')
fig.savefig(dir_figs+str(GCM)+'_PSD_conf_LAB_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


P_Var=AMOC_max_30N ; P_legend='AMOC max 30N' ; P_color='r'; P_c_probability=0.95; P_rho='yes'; P_smooth=5
P_title='Power Spectral Density - AMOC at 30N - yrs '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig=plt.figure()
ff,PSD=plot_PSD_welch_conf(P_Var, 1, P_c_probability, P_rho, P_smooth, P_title, P_color, P_legend, 'lower right')
fig.savefig(dir_figs+str(GCM)+'_PSD_conf_AMOC30N_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


P_Var=AMOC_max_30S ; P_legend='AMOC max 30S' ; P_color='r'; P_c_probability=0.95; P_rho='yes'; P_smooth=5
P_title='Power Spectral Density - AMOC at 30S - yrs '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig=plt.figure()
ff,PSD=plot_PSD_welch_conf(P_Var, 1, P_c_probability, P_rho, P_smooth, P_title, P_color, P_legend, 'lower right')
fig.savefig(dir_figs+str(GCM)+'_PSD_conf_AMOC30S_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


P_Var=AMOC_max_50S ; P_legend='AMOC max 50S' ; P_color='r'; P_c_probability=0.95; P_rho='yes'; P_smooth=5
P_title='Power Spectral Density - AMOC at 50S - yrs '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig=plt.figure()
ff,PSD=plot_PSD_welch_conf(P_Var, 1, P_c_probability, P_rho, P_smooth, P_title, P_color, P_legend, 'lower right')
fig.savefig(dir_figs+str(GCM)+'_PSD_conf_AMOC50S_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


P_Var=SMOC_min_50S ; P_legend='SMOC min 50S' ; P_color='r'; P_c_probability=0.95; P_rho='yes'; P_smooth=5
P_title='Power Spectral Density - AMOC at 50S - yrs '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig=plt.figure()
ff,PSD=plot_PSD_welch_conf(P_Var, 1, P_c_probability, P_rho, P_smooth, P_title, P_color, P_legend, 'lower right')
fig.savefig(dir_figs+str(GCM)+'_PSD_conf_SMOC50S_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')






def plot_PSD_welch_conf(P_Var, P_fq, P_c_probability, P_rho, P_smooth, P_title, P_color, P_legend, P_legend_loc):

### P_c_probability= 0.95 or 0.05 or '-' = confidence interval levels, if the input is not a number then teh confidence intervals won't be plptted
### P_rho = 'yes' or 0.7 or '-' = Red noise significance line, if 'yes' then the Rho and line will be calculated, if P_rho is a number then it will be given as the Rho value, else NO Red noise significance line will be plotted
### P_smooth = 9 (should be even number) or '-' = Number of years for smoothing the PSD line and confidence intervals, if no number is given then NO smoothing will be applied
### P_Var= Plotting variable at X-axis, 1D(X) || P_title=Plot title || P_color=Plot line color || P_legend=Plot variable name to be shown in Plot Legend 
### P_fq=PSD Welch method's sampling frequency of the x time series in units of Hz (Defaults to 1.0) || P_legend_loc=location of legend, 'best' or 'lower left' or 'right' or 'center' or ...
#%% Example :

#P_Var=WS_index ; P_legend='WS Convection Index' ; P_color='r'
#P_title='Power Spectral Density - WS Convection Index with 95% confidence intervals - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
#fig=plt.figure()
#ff,PSD=plot_PSD_welch_conf(P_Var, 1, 0.95, 'yes', 9, P_title, P_color, P_legend, 'lower right')
#plt.ylim(1e-5,1e3) #;plt.xlim(1,500)
#fig.savefig(dir_figs+'name.png', format='png', dpi=300, transparent=True, bbox_inches='tight') 
### 

    from scipy.stats import chi2    
    from scipy import signal

    ff,PSD = signal.welch(P_Var,P_fq) # ff= Array of sample frequencies  ,  PSD = Pxx, Power spectral density or power spectrum of x (which is P_Var)
    X_var=np.linspace(1,len(P_Var)/2+1,len(P_Var)/2+1)
    X_var=ff**(-1); X_var=X_var
    
    if P_rho=='yes':
        Rho = np.nansum( (P_Var[:-1] - np.nanmean(P_Var,axis=0)) * (P_Var[1:] - np.nanmean(P_Var,axis=0)) ,axis=0) / np.nansum( (P_Var[:-1] - np.nanmean(P_Var,axis=0)) * (P_Var[:-1] - np.nanmean(P_Var,axis=0)) ,axis=0)
    elif type(P_rho)==float:
        Rho=P_rho # Rho is the memory parameter
    
    if type(P_smooth)==float or type(P_smooth)==int:
        P_smooth=np.int(P_smooth)
        v = 2*P_smooth # P is the number of estimates in welch function and also the degree of freedom.
        
        sm= np.int( (P_smooth-1)/2)
        P_smooth = np.int( (sm*2)+1 ) # P_smoothhs to be even number for centered smoothing; in case odd number was given, it's changed to an even number by subtracting 1
        PSD_m=copy.deepcopy(PSD) ## Smoothing the variable
        for ii in range(sm,PSD.shape[0]-sm+1):
            PSD_m[ii]=np.nanmean(PSD[ii-sm:ii+sm])
        
        PSD_m=PSD_m[sm:-sm]
        X_var_m=X_var[sm:-sm]
        P_legend=P_legend+' ('+str(np.int(P_smooth))+'yr smoothed)'
        
    else:
        v=2

        PSD_m=copy.deepcopy(PSD)
        X_var_m=copy.deepcopy(X_var)
        
    if type(P_c_probability)==float:
        if P_c_probability < 0.5:
            P_c_probability=1-P_c_probability # In case the P_c_probability is input 0.05 instead of 0.95 for example
        alfa = 1 - P_c_probability

    if P_rho=='yes' or type(P_rho)==float or type(P_rho)==int:
        if type(P_c_probability)!=float: # In case the P_c_probability is not given since confidence interval calculation is not necessary, but red noise significance line is needed
            alfa=0.05
            
        F_x_v = (1-Rho**2) / (1 + Rho**2 - 2*Rho*np.cos(2*np.pi*ff ) )  #  F_x_v is the power spectraum 
        F_x_v_star=np.float( np.real( np.nanmean(PSD,axis=0) / np.nanmean(F_x_v,axis=0) ) ) * F_x_v 
        Pr_alpha = (1/v) * F_x_v_star * np.float( chi2.ppf([1 - alfa], v) )
    
    plt.grid(True,which="both",ls="-", color='0.65')
    plt.loglog(X_var_m,PSD_m, color=P_color, label=P_legend)
    plt.legend(loc='best')
    plt.xlabel('Period (years)', fontsize=18)
    plt.ylabel('Spectral Density', fontsize=18) 
    plt.xticks(fontsize = 20); plt.yticks(fontsize = 20)
    #plt.gca().invert_xaxis()
    if type(P_c_probability)==float:
        Chi = chi2.ppf([1 - alfa / 2, alfa / 2], v)
        
        PSDc_lower = PSD_m * ( v / Chi[0] )
        PSDc_upper = PSD_m * ( v / Chi[1] ) 
        plt.loglog(X_var_m,PSDc_lower, color='g', ls='--', label=str(np.int( (1 - alfa) *100))+'% confidence intervals')
        plt.loglog(X_var_m,PSDc_upper, color='g', ls='--')
    if P_rho=='yes' or type(P_rho)==float or type(P_rho)==int:
        plt.loglog(X_var,Pr_alpha , color='b', ls='--', label=str(np.int( (1 - alfa) *100))+'% Red Noise Significance Level')
    plt.legend(shadow=True, loc=P_legend_loc, prop={'size': 20})
    plt.title(P_title, fontsize=18) 
    plt.show()
    mng = plt.get_current_fig_manager()
    mng.window.showMaximized() # Maximizes the plot window to save figures in full
    return ff,PSD





























