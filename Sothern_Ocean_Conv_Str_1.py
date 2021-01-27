### IMPORTANT: The following path needs to be appended so that the neccessary codes can be read from that directory ###
import sys
sys.path.append('/data1/home/basadieh/behzadcodes/behzadlibrary')
#from Behzadlib import * ### This is loaded from '/data1/home/basadieh/behzadcodes/behzadlibrary'
#from CMIP5lib import *  ### This is loaded from '/data1/home/basadieh/behzadcodes/behzadlibrary'
from Behzadlib import func_MLD, func_time_depth_plot, func_stream, func_ENSO, func_NAO, func_latlon_regrid, func_latlon_regrid_eq, func_oceanlandmask, func_oceanindex, func_regrid, func_plot, func_plot_save, func_plot_bounds_save
from BehzadlibPlot import func_plotmap_contourf, func_plot2Dcolor_contourf, func_plotline_1var, func_plotline_2var, func_plot2Dcontour_1var, func_plot2Dcontour_2var, func_plot_lagcor_sig, func_plot_laggedmaps
from Behzadlib2 import func_barotropicstream_Vx_DyIntegrated, func_barotropicstream_Vy_DxIntegrated, func_MLD_AllYears_annual_ave, func_MLD_AllYears_months_ave, func_MLD_AllYears_1month
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
Plot_unit='(m)'; Plot_title= 'Average MLD [m] (Weddell Sea) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig, m = func_plotmap_contourf(MLD_average_WS, lon_mld_WS,lat_mld_WS, 200, Plot_title, Plot_unit, plt.cm.jet, 'spstere', 210., 80., -80., 'fill')
m.scatter(lon_mld_WS[WS_indeces_lonlat],lat_mld_WS[WS_indeces_lonlat],2,latlon=True)
fig.savefig(dir_figs+str(GCM)+'_average_MLD_WS.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

Plot_unit='(m)'; Plot_title= 'Average MLD [m] (Labrador Sea) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig, m = func_plotmap_contourf(MLD_average_LAB, lon_mld_LAB,lat_mld_LAB, 200, Plot_title, Plot_unit, plt.cm.jet, 'npstere', 210., 80., -80., 'fill')
m.scatter(lon_mld_LAB[LAB_indeces_lonlat],lat_mld_LAB[LAB_indeces_lonlat],2,latlon=True)
fig.savefig(dir_figs+str(GCM)+'_average_MLD_LAB.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
###############################################################################
###############################################################################
###############################################################################
########               Wind Stress/Curl Calculations                ###########
###############################################################################
###############################################################################
###############################################################################

###############################################################################
#################        To restore:        ###################################
#import os
#import shelve
#
#dir_pwd = os.getcwd() # Gets the current directory (and in which the code is placed)
#filename_out = (dir_pwd + '/Results/AllResults_'+GCM+'_500yr.out') # Directory to save processed data
#
#my_shelf = shelve.open(filename_out)
##globals()[key]=my_shelf[key]
#globals()['WS_index_norm_rm']=my_shelf['WS_index_norm_rm']
#globals()['LAB_index_norm_rm']=my_shelf['LAB_index_norm_rm']
#my_shelf.close()

R_ii_WS = np.where(WS_index_norm_rm <-0.5)
R_jj_WS = np.where(WS_index_norm_rm >0.5)

R_ii_LAB = np.where(LAB_index_norm_rm <-0.5)
R_jj_LAB = np.where(LAB_index_norm_rm >0.5)

####################################################################
#######           Westerlies and Trades Calculations        ########
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

variable='tauv'

dset_tauv = netcdf_read (dir_data_in2+str(variable)+'_'+str(t_frequency)+'_'+str(GCM)+'*12.nc',str(GCM), variable) # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
start_date_i,end_date_i = dset_tauv.find_time(dset_tauv.times, year_start, year_end)

dset_tauv_extracted=dset_tauv.extract_data(dset_tauv.variable,start_date_i,end_date_i)
#here we in order to calculate yearly data from monthly data we first reshape array
dset_tauv_yearly=dset_tauv_extracted.reshape(( (year_end-year_start+1) ,12,len(dset_tauv_extracted[0]),len(dset_tauv_extracted[0][0])))
#second - we calculate mean over each 12 elements, which is second dimension of array
dset_tauv_yearly=np.mean(dset_tauv_yearly,axis=1)
Tau_Y = func_regrid(dset_tauv_yearly, dset_tauv.y, dset_tauv.x, Lat_regrid_2D, Lon_regrid_2D)

#dset_tauu.close_ncfile(dset_tauu.fin)
#dset_tauv.close_ncfile(dset_tauv.fin)

Wind_Curl = np.zeros(( Tau_X.shape[0], Tau_X.shape[1], Tau_X.shape[2]))

#for tt in range (0,Tau_X.shape[0]):
#    
#    for ii in range (1,Tau_X.shape[1]-1):
#        for jj in range (1,Tau_X.shape[2]-1): # Wind_Curl = ( D_Tau_Y / D_X ) - ( D_Tau_X / D_Y ) # D_X = (Lon_1 - Lon_2) * COS(Lat)
#            
#            Wind_Curl[tt,ii,jj] = (  ( Tau_Y[tt, ii-1,jj] - Tau_Y[tt, ii+1,jj] ) /  np.absolute(  ( ( Lon_regrid_2D[ii,jj-1] -  Lon_regrid_2D[ii,jj+1] ) * 111321  *  math.cos(math.radians( ( Lat_regrid_2D[ii,jj])))   )  )     )   -   (  ( Tau_X[tt, ii,jj-1] - Tau_X[tt, ii,jj+1] ) / np.absolute( ( ( Lat_regrid_2D[ii-1,jj] -  Lat_regrid_2D[ii+1,jj] ) * 111321 ) )  )

for tt in range (0,Tau_X.shape[0]):
    
    for ii in range (1,Tau_X.shape[1]-1):
        for jj in range (1,Tau_X.shape[2]-1): # Wind_Curl = ( D_Tau_Y / D_X ) - ( D_Tau_X / D_Y ) # D_X = (Lon_1 - Lon_2) * COS(Lat)
            
            Wind_Curl[tt,ii,jj] = (  ( Tau_Y[tt, ii,jj+1] - Tau_Y[tt, ii,jj-1] ) /  np.absolute(  ( ( Lon_regrid_2D[ii,jj+1] -  Lon_regrid_2D[ii,jj-1] ) * 111321  *  math.cos(math.radians( ( Lat_regrid_2D[ii,jj])))   )  )     )   -   (  ( Tau_X[tt, ii+1,jj] - Tau_X[tt, ii-1,jj] ) / np.absolute( ( ( Lat_regrid_2D[ii+1,jj] -  Lat_regrid_2D[ii-1,jj] ) * 111321 ) )  )

        Wind_Curl[tt,ii,0] = (  ( Tau_Y[tt, ii,1] - Tau_Y[tt, ii,-1] ) /  np.absolute(  ( ( Lon_regrid_2D[ii,1] -  Lon_regrid_2D[ii,-1] ) * 111321  *  math.cos(math.radians( ( Lat_regrid_2D[ii,0])))   )  )     )   -   (  ( Tau_X[tt, ii+1,0] - Tau_X[tt, ii-1,0] ) / np.absolute( ( ( Lat_regrid_2D[ii+1,0] -  Lat_regrid_2D[ii-1,0] ) * 111321 ) )  )
        Wind_Curl[tt,ii,-1] = (  ( Tau_Y[tt, ii,0] - Tau_Y[tt, ii,-2] ) /  np.absolute(  ( ( Lon_regrid_2D[ii,0] -  Lon_regrid_2D[ii,-2] ) * 111321  *  math.cos(math.radians( ( Lat_regrid_2D[ii,-1])))   )  )     )   -   (  ( Tau_X[tt, ii+1,-1] - Tau_X[tt, ii-1,-1] ) / np.absolute( ( ( Lat_regrid_2D[ii+1,-1] -  Lat_regrid_2D[ii-1,-1] ) * 111321 ) )  )


Wind_Curl_f = np.zeros(( Tau_X.shape[0], Tau_X.shape[1], Tau_X.shape[2])) # Wind_Crul / f # f = coriolis parameter = 2Wsin(LAT) , W = 7.292E-5 rad/s

for tt in range (0,Tau_X.shape[0]):
    
    for ii in range (1,Tau_X.shape[1]-1):
        if np.absolute( Lat_regrid_2D[ii,0] ) >= 5: # Only calulate for Lats > 5N and Lats < 5S, to avoid infinit numbers in equator where f is zero
            for jj in range (1,Tau_X.shape[2]-1): # Wind_Curl = ( D_Tau_Y / D_X ) - ( D_Tau_X / D_Y ) # D_X = (Lon_1 - Lon_2) * COS(Lat)
            
                Wind_Curl_f[tt,ii,jj] = (  ( ( Tau_Y[tt, ii,jj+1] - Tau_Y[tt, ii,jj-1] ) / ( 2*7.292E-5 *  math.sin(math.radians( ( Lat_regrid_2D[ii,jj]))) ) ) /  np.absolute(  ( ( Lon_regrid_2D[ii,jj+1] -  Lon_regrid_2D[ii,jj-1] ) * 111321  *  math.cos(math.radians( ( Lat_regrid_2D[ii,jj])))   )  )     )   -   (  ( ( Tau_X[tt, ii+1,jj] - Tau_X[tt, ii-1,jj] ) / ( 2*7.292E-5 *  math.sin(math.radians( ( Lat_regrid_2D[ii,jj]))) ) ) / np.absolute( ( ( Lat_regrid_2D[ii+1,jj] -  Lat_regrid_2D[ii-1,jj] ) * 111321 ) )  )

            Wind_Curl_f[tt,ii,0] = (  ( ( Tau_Y[tt, ii,1] - Tau_Y[tt, ii,-1] ) / ( 2*7.292E-5 *  math.sin(math.radians( ( Lat_regrid_2D[ii,0]))) ) ) /  np.absolute(  ( ( Lon_regrid_2D[ii,1] -  Lon_regrid_2D[ii,-1] ) * 111321  *  math.cos(math.radians( ( Lat_regrid_2D[ii,0])))   )  )     )   -   (  ( ( Tau_X[tt, ii+1,jj] - Tau_X[tt, ii-1,0] ) / ( 2*7.292E-5 *  math.sin(math.radians( ( Lat_regrid_2D[ii,0]))) ) ) / np.absolute( ( ( Lat_regrid_2D[ii+1,0] -  Lat_regrid_2D[ii-1,0] ) * 111321 ) )  )
            Wind_Curl_f[tt,ii,-1] = (  ( ( Tau_Y[tt, ii,0] - Tau_Y[tt, ii,-2] ) / ( 2*7.292E-5 *  math.sin(math.radians( ( Lat_regrid_2D[ii,-1]))) ) ) /  np.absolute(  ( ( Lon_regrid_2D[ii,0] -  Lon_regrid_2D[ii,-2] ) * 111321  *  math.cos(math.radians( ( Lat_regrid_2D[ii,-1])))   )  )     )   -   (  ( ( Tau_X[tt, ii+1,-1] - Tau_X[tt, ii-1,-1] ) / ( 2*7.292E-5 *  math.sin(math.radians( ( Lat_regrid_2D[ii,-1]))) ) ) / np.absolute( ( ( Lat_regrid_2D[ii+1,-1] -  Lat_regrid_2D[ii-1,-1] ) * 111321 ) )  )


Wind_Curl_c=Wind_Curl[R_ii_WS]
Wind_Curl_n=Wind_Curl[R_jj_WS]
Wind_Curl_conv_WS=np.nanmean(Wind_Curl_c,axis=0)
Wind_Curl_nonconv_WS=np.nanmean(Wind_Curl_n,axis=0)
Wind_Curl_composite_WS=np.subtract(Wind_Curl_conv_WS, Wind_Curl_nonconv_WS)

Wind_Curl_f_c=Wind_Curl_f[R_ii_WS]
Wind_Curl_f_n=Wind_Curl_f[R_jj_WS]
Wind_Curl_f_conv_WS=np.nanmean(Wind_Curl_f_c,axis=0)
Wind_Curl_f_nonconv_WS=np.nanmean(Wind_Curl_f_n,axis=0)
Wind_Curl_f_composite_WS=np.subtract(Wind_Curl_f_conv_WS, Wind_Curl_f_nonconv_WS)

Wind_Curl_f_c=Wind_Curl_f[R_ii_LAB]
Wind_Curl_f_n=Wind_Curl_f[R_jj_LAB]
Wind_Curl_f_conv_LAB=np.nanmean(Wind_Curl_f_c,axis=0)
Wind_Curl_f_nonconv_LAB=np.nanmean(Wind_Curl_f_n,axis=0)
Wind_Curl_f_composite_LAB=np.subtract(Wind_Curl_f_conv_LAB, Wind_Curl_f_nonconv_LAB)

Tau_X_c=Tau_X[R_ii_WS]
Tau_X_n=Tau_X[R_jj_WS]
Tau_X_conv_WS=np.nanmean(Tau_X_c,axis=0)
Tau_X_nonconv_WS=np.nanmean(Tau_X_n,axis=0)
Tau_X_composite_WS=np.subtract(Tau_X_conv_WS, Tau_X_nonconv_WS)

Tau_X_c=Tau_X[R_ii_LAB]
Tau_X_n=Tau_X[R_jj_LAB]
Tau_X_conv_LAB=np.nanmean(Tau_X_c,axis=0)
Tau_X_nonconv_LAB=np.nanmean(Tau_X_n,axis=0)
Tau_X_composite_LAB=np.subtract(Tau_X_conv_LAB, Tau_X_nonconv_LAB)

###############################################################################
filename_out = (dir_pwd + '/Results/'+GCM+'_500yr_Winds.out') # Directory to save processed data
#################        To Save:        #####################################
import os
import shelve
my_shelf = shelve.open(filename_out,'n') # 'n' for new
var_list_1=['Tau_X', 'Tau_Y','Wind_Curl','Wind_Curl_f', 'var_list_1']
for key in var_list_1:
    try:
        my_shelf[key] = globals()[key]
    except TypeError:
        print('ERROR shelving: {0}'.format(key))
my_shelf.close()
#################        To restore:        ##################################
filename_out = (dir_pwd + '/Results/'+GCM+'_500yr_Winds.out') # Directory to save processed data
import os
import shelve

my_shelf = shelve.open(filename_out)
for key in my_shelf:
    globals()[key]=my_shelf[key]
my_shelf.close()
###############################################################################
###############################################################################
Plot_Var = np.nanmean(Wind_Curl,axis=0) * 1E7
Plot_Var[ Ocean_Land_mask==0 ]=nan # masking over land, so grid cells that fall on land area (value=0) will be deleted
Plot_Var2 = np.nanmean(Tau_X,axis=0)
Plot_Var2 [ Ocean_Land_mask==0 ]=nan # masking over land, so grid cells that fall on land area (value=0) will be deleted

cmap_limit=np.nanmax(np.abs( np.nanpercentile(Plot_Var, 99)))
Plot_range=np.linspace(-cmap_limit,cmap_limit,27)
Plot_unit='(1E-7 N/m3)'; Plot_title= 'Wind Curl (1E-7 N/m3) - (contour lines = Tau_x) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)

fig, m = func_plotmap_contourf(Plot_Var, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.seismic, 'cyl', 210., 80., -80., '-')
im2=m.contour(Lon_regrid_2D, Lat_regrid_2D,Plot_Var2, 20, latlon=True, colors='k')
plt.clabel(im2, fontsize=8, inline=1)
plt.show()
#fig.savefig(dir_figs+'Wind_Curl_'+str(GCM)+'.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

###############################################################################
Plot_Var_f = np.nanmean(Wind_Curl_f,axis=0) * 1E3
Plot_Var2 = np.nanmean(Tau_X,axis=0)
Plot_Var2 [ Ocean_Land_mask==0 ]=nan # masking over land, so grid cells that fall on land area (value=0) will be deleted

cmap_limit=np.nanmax(np.abs( np.nanpercentile(Plot_Var_f, 99)))
Plot_range=np.linspace(-cmap_limit,cmap_limit,27)
Plot_unit='(1E-3 N.S/m3.rad)'; Plot_title= 'Curl of (Wind/f) (Ekman upwelling) (1E-3 N.S/m3.rad) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)+'\n(contour lines: dark = Tau_x, green = Curl(wind/f) )'

fig, m = func_plotmap_contourf(Plot_Var_f, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.seismic, 'cyl', 210., 80., -80., '-')
im2=m.contour(Lon_regrid_2D, Lat_regrid_2D,Plot_Var2, 20, latlon=True, colors='k')
plt.clabel(im2, fontsize=8, inline=1)
im3=m.contour(Lon_regrid_2D[25:50,:], Lat_regrid_2D[25:50,:],Plot_Var_f[25:50,:], levels = [0], latlon=True, colors='darkgreen')
plt.clabel(im3, fontsize=8, inline=1)
plt.show()
#fig.savefig(dir_figs+'Wind_Curl_f_'+str(GCM)+'.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

###############################################################################
Plot_Var3 = np.nanmean(Tau_X,axis=0)
Plot_Var3_l=copy.deepcopy(Plot_Var3)  ; Plot_Var3_l[ Ocean_Land_mask==0 ]=nan

cmap_limit=np.nanmax(np.abs( np.nanpercentile(Plot_Var3, 99.99)))
Plot_range=np.linspace(-cmap_limit,cmap_limit,27)
Plot_unit='(N/m2)'; Plot_title= 'Wind Stress X (N/m2) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)

fig, m = func_plotmap_contourf(Plot_Var3, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.seismic, 'cyl', 210., 80., -80., '-')
im2=m.contour(Lon_regrid_2D[0:70,:], Lat_regrid_2D[0:70,:],Plot_Var3_l[0:70,:], levels = [0], latlon=True, colors='darkgreen')
plt.clabel(im2, fontsize=8, inline=1)
plt.show()
#fig.savefig(dir_figs+'Wind_Tau_X_'+str(GCM)+'.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

###############################################################################
Plot_Var4 = np.nanmean(Tau_Y,axis=0)

cmap_limit=np.nanmax(np.abs( np.nanpercentile(Plot_Var4, 99.99)))
Plot_range=np.linspace(-cmap_limit,cmap_limit,27)
Plot_unit='(N/m2)'; Plot_title= 'Wind Stress Y (N/m2) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)

fig, m = func_plotmap_contourf(Plot_Var4, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.seismic, 'cyl', 210., 80., -80., '-')
#fig.savefig(dir_figs+'Wind_Tau_Y_'+str(GCM)+'.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

###############################################################################
Lat_regrid_1D_4, Lon_regrid_1D_4, Lat_bound_regrid_4, Lon_bound_regrid_4 = func_latlon_regrid_eq(45, 90, lat_min_regrid, lat_max_regrid, lon_min_regrid, lon_max_regrid)
Lon_regrid_2D_4, Lat_regrid_2D_4 = np.meshgrid(Lon_regrid_1D_4, Lat_regrid_1D_4)
Tau_X_4 = func_regrid(np.nanmean(Tau_X,axis=0), Lat_regrid_2D, Lon_regrid_2D, Lat_regrid_2D_4, Lon_regrid_2D_4)
Tau_Y_4 = func_regrid(np.nanmean(Tau_Y,axis=0), Lat_regrid_2D, Lon_regrid_2D, Lat_regrid_2D_4, Lon_regrid_2D_4)

Plot_Var_f = np.nanmean(Wind_Curl_f,axis=0) * 1E3

cmap_limit=np.nanmax(np.abs( np.nanpercentile(Plot_Var_f, 99)))
Plot_range=np.linspace(-cmap_limit,cmap_limit,27)
Plot_unit='(1E-3 N.S/m3.rad)'; Plot_title= 'Curl of (Wind/f) (Ekman upwelling) (1E-3 N.S/m3.rad) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)+'\n(Arrows: wind direction) (contour line: Curl(wind/f)=0)'

fig, m = func_plotmap_contourf(Plot_Var_f, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.seismic, 'cyl', 210., 80., -80., '-')
im2=m.quiver(Lon_regrid_2D_4, Lat_regrid_2D_4, Tau_X_4, Tau_Y_4, latlon=True, pivot='middle')
plt.show()
im3=m.contour(Lon_regrid_2D[25:50,:], Lat_regrid_2D[25:50,:],Plot_Var_f[25:50,:], levels = [0], latlon=True, colors='darkgreen')
plt.clabel(im3, fontsize=8, inline=1)
#fig.savefig(dir_figs+'Wind_Curl_f_WQuiver_'+str(GCM)+'.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

###############################################################################
########      Convective/non-Convective decades composites     ################
Plot_Var_f = copy.deepcopy(Wind_Curl_f_composite_WS) * 1E3

cmap_limit=np.nanmax(np.abs( np.nanpercentile(Plot_Var_f, 99.7)))
Plot_range=np.linspace(-cmap_limit,cmap_limit,27)
Plot_unit='(1E-3 N.S/m3.rad)'; Plot_title= 'Curl of (Wind/f) composites (1E-3 N.S/m3.rad) - Weddell Sea Convection - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)+'\n(contour lines: Curl(wind/f)=0 , yellow=nonconv, green=conv)'

fig, m = func_plotmap_contourf(Plot_Var_f, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.seismic, 'cyl', 210., 80., -80., '-')
im2=m.contour(Lon_regrid_2D[25:50,:], Lat_regrid_2D[25:50,:],Wind_Curl_f_nonconv_WS[25:50,:], levels = [0], latlon=True, colors='gold')
plt.clabel(im2, fontsize=8, inline=1)
im3=m.contour(Lon_regrid_2D[25:50,:], Lat_regrid_2D[25:50,:],Wind_Curl_f_conv_WS[25:50,:], levels = [0], latlon=True, colors='darkgreen')
plt.clabel(im3, fontsize=8, inline=1)
plt.show()
#fig.savefig(dir_figs+'Wind_Curl_f_composites_WS_'+str(GCM)+'.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

###############################################################################
Plot_Var_f = copy.deepcopy(Wind_Curl_f_composite_LAB) * 1E3

cmap_limit=np.nanmax(np.abs( np.nanpercentile(Plot_Var_f, 99.7)))
Plot_range=np.linspace(-cmap_limit,cmap_limit,27)
Plot_unit='(1E-3 N.S/m3.rad)'; Plot_title= 'Curl of (Wind/f) composites (1E-3 N.S/m3.rad) - Labrador Sea Convection - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)+'\n(contour lines: Curl(wind/f)=0 , yellow=nonconv, green=conv)'

fig, m = func_plotmap_contourf(Plot_Var_f, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.seismic, 'cyl', 210., 80., -80., '-')
im2=m.contour(Lon_regrid_2D[25:50,:], Lat_regrid_2D[25:50,:],Wind_Curl_f_nonconv_LAB[25:50,:], levels = [0], latlon=True, colors='gold')
plt.clabel(im2, fontsize=8, inline=1)
im3=m.contour(Lon_regrid_2D[25:50,:], Lat_regrid_2D[25:50,:],Wind_Curl_f_conv_LAB[25:50,:], levels = [0], latlon=True, colors='darkgreen')
plt.clabel(im3, fontsize=8, inline=1)
plt.show()
#fig.savefig(dir_figs+'Wind_Curl_f_composites_LAB_'+str(GCM)+'.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

###############################################################################
Plot_Var_f = copy.deepcopy(Tau_X_composite_WS)
#Plot_Var_f[ Ocean_Land_mask==0 ]=nan # masking over land, so grid cells that fall on land area (value=0) will be deleted
Plot_Var1_l=copy.deepcopy(Tau_X_nonconv_WS)  ; Plot_Var1_l[ Ocean_Land_mask==0 ]=nan
Plot_Var2_l=copy.deepcopy(Tau_X_conv_WS)  ; Plot_Var2_l[ Ocean_Land_mask==0 ]=nan

cmap_limit=np.nanmax(np.abs( np.nanpercentile(Plot_Var_f, 99.99)))
Plot_range=np.linspace(-cmap_limit,cmap_limit,27)
Plot_unit='(1E-3 N.S/m3.rad)'; Plot_title= 'Wind Stress X composites (N/m2) - Weddell Sea Convection - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)+'\n(contour lines: Tau_x=0 , yellow=nonconv, green=conv)'

fig, m = func_plotmap_contourf(Plot_Var_f, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.seismic, 'cyl', 210., 80., -80., '-')
im2=m.contour(Lon_regrid_2D[0:70,:], Lat_regrid_2D[0:70,:],Plot_Var1_l[0:70,:], levels = [0], latlon=True, colors='gold')
plt.clabel(im2, fontsize=8, inline=1)
im3=m.contour(Lon_regrid_2D[0:70,:], Lat_regrid_2D[0:70,:],Plot_Var2_l[0:70,:], levels = [0], latlon=True, colors='darkgreen')
plt.clabel(im3, fontsize=8, inline=1)
plt.show()
#fig.savefig(dir_figs+'Wind_Tau_X_composites_WS_'+str(GCM)+'.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

###############################################################################
Plot_Var_f = copy.deepcopy(Tau_X_composite_LAB)
#Plot_Var_f[ Ocean_Land_mask==0 ]=nan # masking over land, so grid cells that fall on land area (value=0) will be deleted
Plot_Var1_l=copy.deepcopy(Tau_X_nonconv_LAB)  ; Plot_Var1_l[ Ocean_Land_mask==0 ]=nan
Plot_Var2_l=copy.deepcopy(Tau_X_conv_LAB)  ; Plot_Var2_l[ Ocean_Land_mask==0 ]=nan

cmap_limit=np.nanmax(np.abs( np.nanpercentile(Plot_Var_f, 99.99)))
Plot_range=np.linspace(-cmap_limit,cmap_limit,27)
Plot_unit='(1E-3 N.S/m3.rad)'; Plot_title= 'Wind Stress X composites (N/m2) - Labrador Sea Convection - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)+'\n(contour lines: Tau_x=0 , yellow=nonconv, green=conv)'

fig, m = func_plotmap_contourf(Plot_Var_f, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.seismic, 'cyl', 210., 80., -80., '-')
im2=m.contour(Lon_regrid_2D[0:70,:], Lat_regrid_2D[0:70,:],Plot_Var1_l[0:70,:], levels = [0], latlon=True, colors='gold')
plt.clabel(im2, fontsize=8, inline=1)
im3=m.contour(Lon_regrid_2D[0:70,:], Lat_regrid_2D[0:70,:],Plot_Var2_l[0:70,:], levels = [0], latlon=True, colors='darkgreen')
plt.clabel(im3, fontsize=8, inline=1)
plt.show()
#fig.savefig(dir_figs+'Wind_Tau_X_composites_LAB_'+str(GCM)+'.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

###############################################################################

###############################################################################
########################       ESM2G + CM2MC       ############################
dir_data_in1 = (dir_pwd + '/Results_CM2MC/data/') # Directory to raed raw data from
###############################################################################
from numpy import loadtxt

file_name_dir=(dir_data_in1+ 'latitudes.dat'); latitudes = loadtxt(file_name_dir)
file_name_dir=(dir_data_in1+ 'longitudes.dat'); longitudes = loadtxt(file_name_dir)
file_name_dir=(dir_data_in1+ 'tau_x_avg.dat'); tau_x_avg = loadtxt(file_name_dir)
tau_x_avg_mat=tau_x_avg.reshape((latitudes.shape[0], longitudes.shape[0] )); tau_x_avg_mat [ tau_x_avg_mat < -1e10 ] = nan
tau_x_avg_regrid = func_regrid(tau_x_avg_mat, latitudes, longitudes, Lat_regrid_2D, Lon_regrid_2D)

GCM2 = 'GFDL-CM2MC'
dir_figs = (dir_pwd + '/Figures_CM2MC/') # Directory to save figures

Plot_Var3 = np.subtract(tau_x_avg_regrid, np.nanmean(Tau_X,axis=0))
Plot_Var3 [ Ocean_Land_mask==0 ]=nan # masking over land, so grid cells that fall on land area (value=0) will be deleted
Plot_range=np.linspace(-0.2,0.2,81)
Plot_unit='(N/m2)'; Plot_title= 'Wind Stress X (N/m2) - CM2MC minus ESM2G'
fig, m = func_plotmap_contourf(Plot_Var3, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.seismic, 'cyl', 210., 80., -80., '-')
fig.savefig(dir_figs+'Wind_Tau_X__CM2MCminusESM2G.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

###############################################################################
###############################################################################
from numpy import savetxt

dir_pwd = os.getcwd() # Gets the current directory (and in which the code is placed)
dir_results = (dir_pwd + '/ESM2G_Wind/') # Directory to save results
#file_name_dir=(dir_results+ 'tau_x_wed_comp.dat'); tau_x_wed_comp = loadtxt(file_name_dir)

Plot_Var=Lat_regrid_1D
file_name_dir=(dir_results+ 'latitude.dat'); savetxt(file_name_dir, Plot_Var)

Plot_Var=Lon_regrid_1D
file_name_dir=(dir_results+ 'longitude.dat'); savetxt(file_name_dir, Plot_Var)

Plot_Var=np.nanmean(Wind_Curl,axis=0) * 1E7
file_name_dir=(dir_results+ 'tau_x_ave.dat'); savetxt(file_name_dir, Plot_Var)

Plot_Var=np.nanmean(Wind_Curl_f,axis=0) * 1E3
file_name_dir=(dir_results+ 'tau_x_dividedby_f_curl_ave.dat'); savetxt(file_name_dir, Plot_Var)

Plot_Var=Wind_Curl_f_composite_LAB * 1E3
file_name_dir=(dir_results+ 'tau_x_lab_comp.dat'); savetxt(file_name_dir, Plot_Var)

Plot_Var=Wind_Curl_f_composite_WS * 1E3
file_name_dir=(dir_results+ 'tau_x_wed_comp.dat'); savetxt(file_name_dir, Plot_Var)

Plot_Var=Wind_Curl_f_composite_LAB * 1E3
file_name_dir=(dir_results+ 'tau_x_dividedby_f_curl_lab_comp.dat'); savetxt(file_name_dir, Plot_Var)

Plot_Var=Wind_Curl_f_composite_WS * 1E3
file_name_dir=(dir_results+ 'tau_x_dividedby_f_curl_wed_comp.dat'); savetxt(file_name_dir, Plot_Var)



###############################################################################
###############################################################################
###############################################################################
########                 Barotropic Streamfunction                  ###########
###############################################################################
###############################################################################
###############################################################################

t_frequency='Omon'
variable='uo'# vo = sea_water_x_velocity - units: m s-1

dir_data_in1 = ('/data5/scratch/CMIP5/CMIP5_models/ocean_physics/') # Directory to raed raw data from
dir_data_in2=(dir_data_in1+ GCM + '/piControl/mo/')

dset_uo = netcdf_read (dir_data_in2+str(variable)+'_'+str(t_frequency)+'_'+str(GCM)+'*12.nc',str(GCM), variable) # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read

Lat_refrence_H = -30 # Higher bound of lat
Lat_refrence_L = -90 # Lower bound of lat
Stream_function_Barotropic_SO_VxDy_allyears = func_barotropicstream_Vx_DyIntegrated(dset_uo, year_start, year_end, Lat_refrence_H, Lat_refrence_L, 180, 360)
Stream_function_Barotropic_SO_VxDy_ave = np.nanmean (Stream_function_Barotropic_SO_VxDy_allyears, axis=0)

#np.save('Stream_function_Barotropic_SO_VxDy_allyears_GFDL-ESM2G_500yr.npy',Stream_function_Barotropic_SO_VxDy_allyears)
#Stream_function_Barotropic_SO_VxDy_allyears=np.load('Stream_function_Barotropic_SO_VxDy_allyears_GFDL-ESM2G_500yr.npy')

Stream_function_Barotropic_NAtl_VxDy_allyears = func_barotropicstream_Vx_DyIntegrated(dset_uo, year_start, year_end, Lat_refrence_H, Lat_refrence_L, 180, 360)
Stream_function_Barotropic_NAtl_VxDy_ave = np.nanmean (Stream_function_Barotropic_NAtl_VxDy_allyears, axis=0)

##############################################################
variable='vo'# vo = sea_water_y_velocity - units: m s-1

dir_data_in1 = ('/data5/scratch/CMIP5/CMIP5_models/ocean_physics/') # Directory to raed raw data from
dir_data_in2=(dir_data_in1+ GCM + '/piControl/mo/')

dset_vo = netcdf_read (dir_data_in2+str(variable)+'_'+str(t_frequency)+'_'+str(GCM)+'*12.nc',str(GCM), variable) # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read

Lat_refrence_H = 90 # Higher bound of lat
Lat_refrence_L = 0 # Lower bound of lat
Stream_function_Barotropic_NAtl_VyDx_allyears = func_barotropicstream_Vy_DxIntegrated(dset_vo, year_start, year_end, Lat_refrence_H, Lat_refrence_L, 180, 360)
Stream_function_Barotropic_NAtl_VyDx_ave = np.nanmean (Stream_function_Barotropic_NAtl_VyDx_allyears, axis=0)

#np.save('Stream_function_Barotropic_NAtl_VyDx_allyears_GFDL-ESM2G_500yr.npy',Stream_function_Barotropic_NAtl_VyDx_allyears)
#Stream_function_Barotropic_NAtl_VyDx_allyears=np.load('Stream_function_Barotropic_NAtl_VyDx_allyears_GFDL-ESM2G_500yr.npy')
################################################################

Plot_Var = copy.deepcopy(Stream_function_Barotropic_SO_VxDy_ave)
Plot_Var[ Ocean_Land_mask==0 ]=nan # masking over land, so grid cells that fall on land area (value=0) will be deleted
cmap_limit=np.nanmax(np.abs( np.nanpercentile(Plot_Var, 99.999999)))
Plot_range=np.linspace(-cmap_limit,cmap_limit,101) ### Or:  Plot_range=100
Plot_unit='(Sv)'; Plot_title= 'Barotropic Streamfunction (Sv) [Integral of (Vx*dZ*dY) - 30S-90S]\nyears '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig, m = func_plotmap_contourf(Plot_Var[0:61,:], Lon_regrid_2D[0:61,:], Lat_regrid_2D[0:61,:], Plot_range, Plot_title, Plot_unit, plt.cm.seismic, 'spstere', 210., 80., -80., 'fill')
#fig.savefig(dir_figs+str(GCM)+'_Stream_function_SO_barotropic_ave_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

Plot_Var = copy.deepcopy(Stream_function_Barotropic_SO_VxDy_ave)
Plot_Var[ Ocean_Land_mask==0 ]=nan # masking over land, so grid cells that fall on land area (value=0) will be deleted
cmap_limit=np.nanmax(np.abs( np.nanpercentile(Plot_Var, 99.999999)))
Plot_range=np.linspace(-cmap_limit,cmap_limit,101) ### Or:  Plot_range=100
Plot_unit='(Sv)'; Plot_title= 'Barotropic Streamfunction (Sv) [Integral of (Vx*dZ*dY) - 30S-90S] - years '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig, m = func_plotmap_contourf(Plot_Var[0:61,:], Lon_regrid_2D[0:61,:], Lat_regrid_2D[0:61,:], Plot_range, Plot_title, Plot_unit, plt.cm.seismic, 'cyl', 210., 80., -80., 'fill')
#fig.savefig(dir_figs+str(GCM)+'_Stream_function_SO_barotropic_ave_2.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

Plot_Var = copy.deepcopy(Stream_function_Barotropic_NAtl_VxDy_ave)
Plot_Var[ Ocean_Land_mask==0 ]=nan # masking over land, so grid cells that fall on land area (value=0) will be deleted
cmap_limit=np.nanmax(np.abs( np.nanpercentile(Plot_Var, 99.999999)))
Plot_range=np.linspace(-cmap_limit,cmap_limit,101) ### Or:  Plot_range=100
Plot_unit='(Sv)'; Plot_title= 'Barotropic Streamfunction (Sv) [Integral of (Vx*dZ*dY) - 5N-90N]\nyears '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig, m = func_plotmap_contourf(Plot_Var[95:180,:], Lon_regrid_2D[95:180,:], Lat_regrid_2D[95:180,:], Plot_range, Plot_title, Plot_unit, plt.cm.seismic, 'npstere', 210., 80., -80., 'fill')
#fig.savefig(dir_figs+str(GCM)+'_Stream_function_NAtl_barotropic_ave_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

Plot_Var = copy.deepcopy(Stream_function_Barotropic_NAtl_VxDy_ave)
Plot_Var[ Ocean_Land_mask==0 ]=nan # masking over land, so grid cells that fall on land area (value=0) will be deleted
cmap_limit=np.nanmax(np.abs( np.nanpercentile(Plot_Var, 99.999999)))
Plot_range=np.linspace(-cmap_limit,cmap_limit,101) ### Or:  Plot_range=100
Plot_unit='(Sv)'; Plot_title= 'Barotropic Streamfunction (Sv) [Integral of (Vx*dZ*dY) - 5N-90N] - years '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig, m = func_plotmap_contourf(Plot_Var[95:180,:], Lon_regrid_2D[95:180,:], Lat_regrid_2D[95:180,:], Plot_range, Plot_title, Plot_unit, plt.cm.seismic, 'cyl', 210., 80., -80., 'fill')
#fig.savefig(dir_figs+str(GCM)+'_Stream_function_NAtl_barotropic_ave_2.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

Plot_Var = copy.deepcopy(Stream_function_Barotropic_NAtl_VyDx_ave)
Plot_Var[ Ocean_Land_mask==0 ]=nan # masking over land, so grid cells that fall on land area (value=0) will be deleted
cmap_limit=np.nanmax(np.abs( np.nanpercentile(Plot_Var, 99.999999)))
Plot_range=np.linspace(-cmap_limit,cmap_limit,101) ### Or:  Plot_range=100
Plot_unit='(Sv)'; Plot_title= 'Barotropic Streamfunction (Sv) [Integral of (Vy*dZ*dX) - E to W integration]\nyears '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig, m = func_plotmap_contourf(Plot_Var[95:180,:], Lon_regrid_2D[95:180,:], Lat_regrid_2D[95:180,:], Plot_range, Plot_title, Plot_unit, plt.cm.seismic, 'npstere', 210., 80., -80., 'fill')
#fig.savefig(dir_figs+str(GCM)+'_Stream_function_NAtl_barotropic_ave_3.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

Plot_Var = copy.deepcopy(Stream_function_Barotropic_NAtl_VyDx_ave)
Plot_Var[ Ocean_Land_mask==0 ]=nan # masking over land, so grid cells that fall on land area (value=0) will be deleted
cmap_limit=np.nanmax(np.abs( np.nanpercentile(Plot_Var, 99.999999)))
Plot_range=np.linspace(-cmap_limit,cmap_limit,101) ### Or:  Plot_range=100
Plot_unit='(Sv)'; Plot_title= 'Barotropic Streamfunction (Sv) [Integral of (Vy*dZ*dX) - E to W integration] - years '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig, m = func_plotmap_contourf(Plot_Var[95:180,:], Lon_regrid_2D[95:180,:], Lat_regrid_2D[95:180,:], Plot_range, Plot_title, Plot_unit, plt.cm.seismic, 'cyl', 210., 80., -80., 'fill')
#fig.savefig(dir_figs+str(GCM)+'_Stream_function_NAtl_barotropic_ave_4.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


########                   Gyre Strength Time Series                ###########

Stream_function_Barotropic_SO_VxDy_allyears_WS_max = np.nanmax (Stream_function_Barotropic_SO_VxDy_allyears [:,0:61,291:] , axis =2)
Stream_function_Barotropic_SO_VxDy_allyears_WS_max = np.nanmax (Stream_function_Barotropic_SO_VxDy_allyears_WS_max , axis =1)
Stream_function_Barotropic_SO_VxDy_allyears_WS_max_rm = copy.deepcopy(Stream_function_Barotropic_SO_VxDy_allyears_WS_max)
Stream_function_Barotropic_SO_VxDy_allyears_WS_max_rm = runningMeanFast(Stream_function_Barotropic_SO_VxDy_allyears_WS_max_rm, 10)#if data is not smoothed

years=np.linspace(year_start,year_end,year_end-year_start+1)
P_Var_x=years[:-9]
P_Var_y=Stream_function_Barotropic_SO_VxDy_allyears_WS_max_rm[:-9]
P_xlable='Years'
P_ylable='Barometric Streamfunction Maxima (Sv)'
P_title='SO Barometric Streamfunction Maxima (Weddell Sea [70W-0W, 90S-30S] gyre strength, smoothed) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plotline_1var(P_Var_x, P_Var_y, P_xlable, P_ylable, P_title, 'darkcyan', '-')
fig.savefig(dir_figs+str(GCM)+'_Stream_function_SO_barotropic_max.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


Stream_function_Barotropic_NAtl_VyDx_allyears_LAB_B1_max = np.nanmax (Stream_function_Barotropic_NAtl_VyDx_allyears [:,140:156,295:320] , axis =2)
Stream_function_Barotropic_NAtl_VyDx_allyears_LAB_B1_max = np.nanmax (Stream_function_Barotropic_NAtl_VyDx_allyears_LAB_B1_max , axis =1)
Stream_function_Barotropic_NAtl_VyDx_allyears_LAB_B1_max_rm = copy.deepcopy(Stream_function_Barotropic_NAtl_VyDx_allyears_LAB_B1_max)
Stream_function_Barotropic_NAtl_VyDx_allyears_LAB_B1_max_rm = runningMeanFast(Stream_function_Barotropic_NAtl_VyDx_allyears_LAB_B1_max_rm, 10)#if data is not smoothed

years=np.linspace(year_start,year_end,year_end-year_start+1)
P_Var_x=years[:-9]
P_Var_y=Stream_function_Barotropic_NAtl_VyDx_allyears_LAB_B1_max_rm[:-9]
P_xlable='Years'
P_ylable='Barometric Streamfunction Maxima (Sv)'
P_title='N.Atl. Barometric Streamfunction Maxima (Labrador Sea [65W-41W, 50N-65N] gyre strength, smoothed) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plotline_1var(P_Var_x, P_Var_y, P_xlable, P_ylable, P_title, 'darkcyan', '-')
fig.savefig(dir_figs+str(GCM)+'_Stream_function_LAB_B1_barotropic_max.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


from BehzadlibPlot import plot_PSD_welch_conf
P_Var=Stream_function_Barotropic_SO_VxDy_allyears_WS_max ; P_legend='SO Bar. Streamfunction Maxima' ; P_color='r'; P_c_probability='-'; P_rho='-'; P_smooth='-'
P_title='Power Spectral Density - SO Barometric Streamfunction Maxima - yrs '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig=plt.figure()
ff,PSD=plot_PSD_welch_conf(P_Var, 1, P_c_probability, P_rho, P_smooth, P_title, P_color, P_legend, 'lower right')
fig.savefig(dir_figs+str(GCM)+'_PSD_BarStr_WSgyre_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

P_Var=Stream_function_Barotropic_NAtl_VyDx_allyears_LAB_B1_max ; P_legend='N.Atl. Bar. Streamfunction Maxima' ; P_color='r'; P_c_probability='-'; P_rho='-'; P_smooth='-'
P_title='Power Spectral Density - N.Atl. Barometric Streamfunction Maxima - yrs '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig=plt.figure()
ff,PSD=plot_PSD_welch_conf(P_Var, 1, P_c_probability, P_rho, P_smooth, P_title, P_color, P_legend, 'lower right')
fig.savefig(dir_figs+str(GCM)+'_PSD_BarStr_LABgyre_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

#np.save('Stream_function_Barotropic_NAtl_VyDx_allyears_LAB_B1_max_rm_GFDL-ESM2G_500yr.npy',Stream_function_Barotropic_NAtl_VyDx_allyears_LAB_B1_max_rm)
#np.save('Stream_function_Barotropic_SO_VxDy_allyears_WS_max_rm_GFDL-ESM2G_500yr.npy',Stream_function_Barotropic_SO_VxDy_allyears_WS_max_rm)

import os
import shelve

dir_pwd = os.getcwd() # Gets the current directory (and in which the code is placed)
filename_out = (dir_pwd + '/Results/'+GCM+'_500yr_Winds.out') # Directory to save processed data

my_shelf = shelve.open(filename_out)
#globals()[key]=my_shelf[key]
globals()['Tau_X']=my_shelf['Tau_X']
globals()['Wind_Curl_f']=my_shelf['Wind_Curl_f']
my_shelf.close()

Wind_Curl_f_40S60S_0W30W=Wind_Curl_f[:,30:51,330:]
Wind_Curl_f_40S60S_0W30W=np.nanmean( Wind_Curl_f_40S60S_0W30W, axis=2)
Wind_Curl_f_40S60S_0W30W=np.nanmean( Wind_Curl_f_40S60S_0W30W, axis=1)
Wind_Curl_f_40S60S_0W30W_m=runningMeanFast(Wind_Curl_f_40S60S_0W30W, 10)

Wind_Curl_f_50N65N_65W20W=Wind_Curl_f[:,140:156,295:341]
Wind_Curl_f_50N65N_65W20W=np.nanmean( Wind_Curl_f_50N65N_65W20W, axis=2)
Wind_Curl_f_50N65N_65W20W=np.nanmean( Wind_Curl_f_50N65N_65W20W, axis=1)
Wind_Curl_f_50N65N_65W20W_m=runningMeanFast(Wind_Curl_f_50N65N_65W20W, 10)

Wind_Curl_f_50N65N_65W41W=Wind_Curl_f[:,140:156,295:320]
Wind_Curl_f_50N65N_65W41W=np.nanmean( Wind_Curl_f_50N65N_65W41W, axis=2)
Wind_Curl_f_50N65N_65W41W=np.nanmean( Wind_Curl_f_50N65N_65W41W, axis=1)
Wind_Curl_f_50N65N_65W41W_m=runningMeanFast(Wind_Curl_f_50N65N_65W41W, 10)

Wind_Curl_f_65N80N_5W20E_1=Wind_Curl_f[:,155:171,355:360]
Wind_Curl_f_65N80N_5W20E_2=Wind_Curl_f[:,155:171,0:21]
Wind_Curl_f_65N80N_5W20E=np.concatenate((Wind_Curl_f_65N80N_5W20E_1, Wind_Curl_f_65N80N_5W20E_2), axis=2)
Wind_Curl_f_65N80N_5W20E=np.nanmean(Wind_Curl_f_65N80N_5W20E, axis=2)
Wind_Curl_f_65N80N_5W20E=np.nanmean(Wind_Curl_f_65N80N_5W20E, axis=1)
Wind_Curl_f_65N80N_5W20E_m=runningMeanFast(Wind_Curl_f_65N80N_5W20E, 10)

North_Westerlies=Zonal_winds[:,135:145]

AMOC_max_50S= AMOC_max[:,40] # AMOC at 50S # Max of streamfunction method
AMOC_max_50S_m=runningMeanFast(AMOC_max_50S, 10)
AMOC_max_30S= AMOC_max[:,60] # AMOC at 30S # Max of streamfunction method
AMOC_max_30S_m=runningMeanFast(AMOC_max_30S, 10)
AMOC_max_50N= AMOC_max[:,140] # AMOC at 50N # Max of streamfunction method
AMOC_max_50N_m=runningMeanFast(AMOC_max_50N, 10)
AMOC_max_30N= AMOC_max[:,120] # AMOC at 30N # Max of streamfunction method
AMOC_max_30N_m=runningMeanFast(AMOC_max_30N, 10)


P_title='Peak Labrador Sea Convection lagged correlation - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig=plt.figure()
func_plot_lagcor_sig(LAB_index_rm[:-9], Stream_function_Barotropic_NAtl_VyDx_allyears_LAB_B1_max_rm[:-9], 40, P_title, 'y', 'LAB Gyre Strength ([65W-41W, 50N-65N] Bar.Str.Func.Max)', 'best', '-')
func_plot_lagcor_sig(LAB_index_rm[:-9], Stream_function_Barotropic_SO_VxDy_allyears_WS_max_rm[:-9], 40, P_title, 'm', 'WS Gyre Strength ([70W-0W, 90S-30S] Bar.Str.Func.Max)', 'best', '-')
func_plot_lagcor_sig(LAB_index_rm[:-9], AMOC_max_50N_m[:-9], 40, P_title, 'b', 'AMOC max 50N', 'best', '-')
func_plot_lagcor_sig(LAB_index_rm[:-9], AMOC_max_30N_m[:-9], 40, P_title, 'g', 'AMOC max 30N', 'best', '-')
func_plot_lagcor_sig(LAB_index_rm[:-9], AMOC_max_30S_m[:-9], 40, P_title, 'darkorange', 'AMOC max 30S', 'best', '-')
func_plot_lagcor_sig(LAB_index_rm[:-9], AMOC_max_50S_m[:-9], 40, P_title, 'r', 'AMOC max 50S', 'best', '-')
plt.ylim(-1,1)
fig.savefig(dir_figs+str(GCM)+'_lagcor_LABconvection_BarStr_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


P_title='Peak LAB Gyre Strength ([65W-41W, 50N-65N] Bar.Str.Func.Max) lagged correlation - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig=plt.figure()
func_plot_lagcor_sig(Stream_function_Barotropic_NAtl_VyDx_allyears_LAB_B1_max_rm[:-9], LAB_index_rm[:-9], 40, P_title, 'y', 'LAB index', 'best', '-')
func_plot_lagcor_sig(Stream_function_Barotropic_NAtl_VyDx_allyears_LAB_B1_max_rm[:-9], WS_index_rm[:-9], 40, P_title, 'm', 'WS index', 'best', '-')
func_plot_lagcor_sig(Stream_function_Barotropic_NAtl_VyDx_allyears_LAB_B1_max_rm[:-9], AMOC_max_50N_m[:-9], 40, P_title, 'b', 'AMOC max 50N', 'best', '-')
func_plot_lagcor_sig(Stream_function_Barotropic_NAtl_VyDx_allyears_LAB_B1_max_rm[:-9], AMOC_max_30N_m[:-9], 40, P_title, 'g', 'AMOC max 30N', 'best', '-')
func_plot_lagcor_sig(Stream_function_Barotropic_NAtl_VyDx_allyears_LAB_B1_max_rm[:-9], AMOC_max_30S_m[:-9], 40, P_title, 'darkorange', 'AMOC max 30S', 'best', '-')
func_plot_lagcor_sig(Stream_function_Barotropic_NAtl_VyDx_allyears_LAB_B1_max_rm[:-9], AMOC_max_50S_m[:-9], 40, P_title, 'r', 'AMOC max 50S', 'best', '-')
plt.ylim(-1,1)
fig.savefig(dir_figs+str(GCM)+'_lagcor_LAB_BarStr_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


P_title='Peak LAB Gyre Strength ([65W-41W, 50N-65N] Bar.Str.Func.Max) lagged correlation - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig=plt.figure()
func_plot_lagcor_sig(Stream_function_Barotropic_NAtl_VyDx_allyears_LAB_B1_max_rm[:-9], LAB_index_rm[:-9], 40, P_title, 'y', 'LAB index', 'best', '-')
func_plot_lagcor_sig(Stream_function_Barotropic_NAtl_VyDx_allyears_LAB_B1_max_rm[:-9], WS_index_rm[:-9], 40, P_title, 'm', 'WS index', 'best', '-')
func_plot_lagcor_sig(Stream_function_Barotropic_NAtl_VyDx_allyears_LAB_B1_max_rm[:-9], Wind_Curl_f_40S60S_0W30W_m[:-9], 40, P_title, 'r', 'Curl of (Wind/f) WS [40S60S,0W30W]', 'best', '-')
func_plot_lagcor_sig(Stream_function_Barotropic_NAtl_VyDx_allyears_LAB_B1_max_rm[:-9], Wind_Curl_f_50N65N_65W20W_m[:-9], 40, P_title, 'b', 'Curl of (Wind/f) LAB [50N65N_65W20W]', 'best', '-')
func_plot_lagcor_sig(Stream_function_Barotropic_NAtl_VyDx_allyears_LAB_B1_max_rm[:-9], Wind_Curl_f_50N65N_65W41W_m[:-9], 40, P_title, 'k', 'Curl of (Wind/f) LAB West [50N65N_65W41W]', 'best', '-')
func_plot_lagcor_sig(Stream_function_Barotropic_NAtl_VyDx_allyears_LAB_B1_max_rm[:-9], Wind_Curl_f_65N80N_5W20E_m[:-9], 40, P_title, 'g', 'Curl of (Wind/f) NOR [65N80N_5W20E]', 'best', '-')
func_plot_lagcor_sig(Stream_function_Barotropic_NAtl_VyDx_allyears_LAB_B1_max_rm[:-9], NAO_m[:-9], 40, P_title, 'c', 'NAO index', 'best', '-')
plt.ylim(-1,1)
fig.savefig(dir_figs+str(GCM)+'_lagcor_LAB_BarStr_2.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


P_title='Peak WS Gyre Strength ([70W-0W, 90S-30S] Bar.Str.Func.Max) lagged correlation - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig=plt.figure()
func_plot_lagcor_sig(Stream_function_Barotropic_SO_VxDy_allyears_WS_max_rm[:-9], LAB_index_rm[:-9], 40, P_title, 'y', 'LAB index', 'best', '-')
func_plot_lagcor_sig(Stream_function_Barotropic_SO_VxDy_allyears_WS_max_rm[:-9], WS_index_rm[:-9], 40, P_title, 'm', 'WS index', 'best', '-')
func_plot_lagcor_sig(Stream_function_Barotropic_SO_VxDy_allyears_WS_max_rm[:-9], AMOC_max_50N_m[:-9], 40, P_title, 'b', 'AMOC max 50N', 'best', '-')
func_plot_lagcor_sig(Stream_function_Barotropic_SO_VxDy_allyears_WS_max_rm[:-9], AMOC_max_30N_m[:-9], 40, P_title, 'g', 'AMOC max 30N', 'best', '-')
func_plot_lagcor_sig(Stream_function_Barotropic_SO_VxDy_allyears_WS_max_rm[:-9], AMOC_max_30S_m[:-9], 40, P_title, 'darkorange', 'AMOC max 30S', 'best', '-')
func_plot_lagcor_sig(Stream_function_Barotropic_SO_VxDy_allyears_WS_max_rm[:-9], AMOC_max_50S_m[:-9], 40, P_title, 'r', 'AMOC max 50S', 'best', '-')
plt.ylim(-1,1)
fig.savefig(dir_figs+str(GCM)+'_lagcor_WS_BarStr_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


P_title='Peak WS Gyre Strength ([70W-0W, 90S-30S] Bar.Str.Func.Max) lagged correlation - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig=plt.figure()
func_plot_lagcor_sig(Stream_function_Barotropic_SO_VxDy_allyears_WS_max_rm[:-9], LAB_index_rm[:-9], 40, P_title, 'y', 'LAB index', 'best', '-')
func_plot_lagcor_sig(Stream_function_Barotropic_SO_VxDy_allyears_WS_max_rm[:-9], WS_index_rm[:-9], 40, P_title, 'm', 'WS index', 'best', '-')
func_plot_lagcor_sig(Stream_function_Barotropic_SO_VxDy_allyears_WS_max_rm[:-9], Wind_Curl_f_40S60S_0W30W_m[:-9], 40, P_title, 'r', 'Curl of (Wind/f) WS [40S60S,0W30W]', 'best', '-')
func_plot_lagcor_sig(Stream_function_Barotropic_SO_VxDy_allyears_WS_max_rm[:-9], Wind_Curl_f_50N65N_65W20W_m[:-9], 40, P_title, 'b', 'Curl of (Wind/f) LAB [50N65N_65W20W]', 'best', '-')
func_plot_lagcor_sig(Stream_function_Barotropic_SO_VxDy_allyears_WS_max_rm[:-9], Wind_Curl_f_50N65N_65W41W_m[:-9], 40, P_title, 'k', 'Curl of (Wind/f) LAB West [50N65N_65W41W]', 'best', '-')
func_plot_lagcor_sig(Stream_function_Barotropic_SO_VxDy_allyears_WS_max_rm[:-9], Wind_Curl_f_65N80N_5W20E_m[:-9], 40, P_title, 'g', 'Curl of (Wind/f) NOR [65N80N_5W20E]', 'best', '-')
func_plot_lagcor_sig(Stream_function_Barotropic_SO_VxDy_allyears_WS_max_rm[:-9], NAO_m[:-9], 40, P_title, 'c', 'NAO index', 'best', '-')
plt.ylim(-1,1)
fig.savefig(dir_figs+str(GCM)+'_lagcor_WS_BarStr_2.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


########            Barotropic Streamfunction Composites            ###########

Stream_function_Barotropic_NAtl_VyDx_lab_c=Stream_function_Barotropic_NAtl_VyDx_allyears[R_ii_LAB]
Stream_function_Barotropic_NAtl_VyDx_lab_n=Stream_function_Barotropic_NAtl_VyDx_allyears[R_jj_LAB]

composite_Stream_function_Barotropic_NAtl_LAB=np.subtract( np.nanmean(Stream_function_Barotropic_NAtl_VyDx_lab_c,axis=0) , np.nanmean(Stream_function_Barotropic_NAtl_VyDx_lab_n,axis=0))

Plot_Var = composite_Stream_function_Barotropic_NAtl_LAB
Plot_Var[ Ocean_Land_mask==0 ]=nan # masking over land, so grid cells that fall on land area (value=0) will be deleted
cmap_limit=6  #np.nanmax(np.abs( np.nanpercentile(Plot_Var, 99.999999999999)))
Plot_range=np.linspace(-cmap_limit,cmap_limit,101) ### Or:  Plot_range=100
Plot_unit='(Sv)'; Plot_title= 'Labrador Sea Convection composites for Barotropic Streamfunction (Sv) (conv - nonconve)\n[Integral of (Vy*dZ*dX) - E to W integration] - years '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig, m = func_plotmap_contourf(Plot_Var, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.seismic, 'cyl', 210., 80., -80., 'fill')
fig.savefig(dir_figs+str(GCM)+'_Stream_function_NAtl_barotropic_composites_LAB.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


Plot_Var1_l=copy.deepcopy( np.nanmean(Stream_function_Barotropic_NAtl_VyDx_lab_c,axis=0) )  ; Plot_Var1_l[ Ocean_Land_mask==0 ]=nan
Plot_Var2_l=copy.deepcopy( np.nanmean(Stream_function_Barotropic_NAtl_VyDx_lab_n,axis=0) )  ; Plot_Var2_l[ Ocean_Land_mask==0 ]=nan

Plot_Var = Stream_function_Barotropic_NAtl_VyDx_ave
Plot_Var[ Ocean_Land_mask==0 ]=nan # masking over land, so grid cells that fall on land area (value=0) will be deleted
cmap_limit=np.nanmax(np.abs( np.nanpercentile(Plot_Var, 99.999999999999)))
Plot_range=np.linspace(-cmap_limit,cmap_limit,101) ### Or:  Plot_range=100
Plot_unit='(Sv)'; Plot_title= 'Barotropic Streamfunction Average (Sv) and Labrador Sea convection contours - '+str(GCM)+'\n(contour lines: green=LAB conve gyre, yellow=LAB nonconv gyre, black=average)'

fig=plt.figure()
m = Basemap( projection='cyl',lon_0=210., llcrnrlon=260.,llcrnrlat=0.,urcrnrlon=390.,urcrnrlat=80.) 
m.drawparallels(np.arange(0, 80+0.001, 40.),labels=[True,False,False,False], linewidth=0.01, color='k', fontsize=20) # labels = [left,right,top,bottom] # Latitutes
m.drawmeridians(np.arange(210-180,210+180,60.),labels=[False,False,False,True], linewidth=0.01, color='k', fontsize=20) # labels = [left,right,top,bottom] # Longitudes        
m.fillcontinents(color='0.8')
m.drawcoastlines(linewidth=1.0, linestyle='solid', color='k', antialiased=1, ax=None, zorder=None)
im=m.contourf(Lon_regrid_2D, Lat_regrid_2D, Plot_Var,Plot_range,latlon=True, cmap=plt.cm.seismic, extend='both')
cbar = m.colorbar(im,"right", size="3%", pad="2%")
cbar.ax.tick_params(labelsize=20) 
cbar.set_label(Plot_unit)
plt.show()
plt.title(Plot_title, fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
im2=m.contour(Lon_regrid_2D[125:155,:], Lat_regrid_2D[125:155,:],Plot_Var1_l[125:155,:], levels = [0], latlon=True, colors='darkgreen')
plt.clabel(im2, fontsize=2, inline=1)
im3=m.contour(Lon_regrid_2D[125:155,:], Lat_regrid_2D[125:155,:],Plot_Var2_l[125:155,:], levels = [0], latlon=True, colors='gold')
plt.clabel(im3, fontsize=2, inline=1)
im4=m.contour(Lon_regrid_2D[125:155,:], Lat_regrid_2D[125:155,:],Plot_Var[125:155,:], levels = [0], latlon=True, colors='k', linestyles='--', linewidths=1)
plt.clabel(im3, fontsize=2, inline=1)
fig.savefig(dir_figs+str(GCM)+'_Stream_function_NAtl_barotropic_composites_LAB_2.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


########                 Barotropic Streamfunction                  ###########



###############################################################################
###############################################################################
########                      MLD calculations                      ###########
###############################################################################
###############################################################################

#Density_allyears=np.load('Density_allyears_GFDL-ESM2G_500yr.npy')
#MLD_allyears_2= np.zeros(( Density_allyears.shape[0], Density_allyears.shape[2], Density_allyears.shape[3]))
#
#for tt in range (0,Density_allyears.shape[0]): # Time loop
#    for ii in range (0, Density_allyears.shape[2]): # Latitude loop
#        for jj in range (0, Density_allyears.shape[3]): # Longitude loop
#            
#            if Ocean_Land_mask[ii,jj] == 1:
#            
#                density_10m = ( Density_allyears[tt,0,ii,jj] + Density_allyears[tt,1,ii,jj] ) / 2
#                
#                for dd in range (0, Density_allyears.shape[1]): # Depth loop
#                    if Density_allyears[tt,dd,ii,jj] - density_10m > 0.03:
#                        if dd==0:
#                            MLD_allyears_2[tt,ii,jj] =  Depths [0]
#                            break
#                        else:
#                            MLD_allyears_2[tt,ii,jj] =  Depths [dd-1]
#                            break
###############################################################################

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

MLD_years_monthlyave = func_MLD_AllYears_months_ave(dset_thetao, dset_so, year_start, year_end, 180, 360)

np.save('MLD_years_monthlyave_GFDL-ESM2G_500yr.npy',MLD_years_monthlyave)
#MLD_years_monthlyave=np.load('MLD_years_monthlyave_GFDL-ESM2G_500yr.npy')

MLD_years_annualave = func_MLD_AllYears_annual_ave(dset_thetao, dset_so, year_start, year_end, 180, 360)

month_no=9 # args[2] in main code # month of the year, 9=September, 3=March
MLD_years_september = func_MLD_AllYears_1month(dset_thetao, dset_so, month_no, year_start, year_end, 180, 360)
month_no=3 # args[2] in main code # month of the year, 9=September, 3=March
MLD_years_march = func_MLD_AllYears_1month(dset_thetao, dset_so, month_no, year_start, year_end, 180, 360)


Plot_Var = np.nanmean( MLD_years_september, axis=0)
cmap_limit=np.nanmax(np.abs( np.nanpercentile(Plot_Var, 99.9999)))
Plot_range=np.linspace(0,cmap_limit,27) ### Or:  Plot_range=100
Plot_unit='(m)'; Plot_title= 'Average MLD (meters) - September only - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig, m = func_plotmap_contourf(Plot_Var, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.jet, 'cyl', 210., 80., -80., 'fill')
fig.savefig(dir_figs+'MLD_ave_SeptemberOnly_5_'+str(GCM)+'.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

Plot_Var = np.nanmean( MLD_years_march, axis=0)
cmap_limit=np.nanmax(np.abs( np.nanpercentile(Plot_Var, 99.9999)))
Plot_range=np.linspace(0,cmap_limit,27) ### Or:  Plot_range=100
Plot_unit='(m)'; Plot_title= 'Average MLD (meters) - March only - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig, m = func_plotmap_contourf(Plot_Var, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.jet, 'cyl', 210., 80., -80., 'fill')
fig.savefig(dir_figs+'MLD_ave_MarchOnly_5_'+str(GCM)+'.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

Plot_Var = np.nanmean( MLD_years_annualave, axis=0)
cmap_limit=np.nanmax(np.abs( np.nanpercentile(Plot_Var, 99.9999)))
Plot_range=np.linspace(0,cmap_limit,27) ### Or:  Plot_range=100
Plot_unit='(m)'; Plot_title= 'Average MLD (meters) - Annual average - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig, m = func_plotmap_contourf(Plot_Var, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.jet, 'cyl', 210., 80., -80., 'fill')
fig.savefig(dir_figs+'MLD_ave_Annual_5_'+str(GCM)+'.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

Plot_Var = np.nanmean( MLD_years_monthlyave, axis=0)
cmap_limit=np.nanmax(np.abs( np.nanpercentile(Plot_Var, 99.9999)))
Plot_range=np.linspace(0,cmap_limit,27) ### Or:  Plot_range=100
Plot_unit='(m)'; Plot_title= 'Average MLD (meters) - Monthly average - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig, m = func_plotmap_contourf(Plot_Var, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.jet, 'cyl', 210., 80., -80., 'fill')
fig.savefig(dir_figs+'MLD_ave_Monthly_500yrs_'+str(GCM)+'.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

Plot_Var = np.nanstd( MLD_years_monthlyave, axis=0)
cmap_limit=np.nanmax(np.abs( np.nanpercentile(Plot_Var, 99.9999)))
Plot_range=np.linspace(0,cmap_limit,27) ### Or:  Plot_range=100
Plot_unit='(m)'; Plot_title= 'St.Dev. of monthly averaged MLD (meters) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)+'\n( contour line: St.Dev.(MLD)=200m )'
fig, m = func_plotmap_contourf(Plot_Var, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.jet, 'cyl', 210., 80., -80., 'fill')
im2=m.contour(Lon_regrid_2D, Lat_regrid_2D,Plot_Var, levels = [200], latlon=True, colors='firebrick')
plt.clabel(im2, fontsize=8, inline=5)
fig.savefig(dir_figs+'MLD_StDev_Monthly_500yrs_'+str(GCM)+'.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

###############################################################################
###############################################################################
########                      MLD calculations                      ###########
###############################################################################
###############################################################################
import os
import shelve

dir_pwd = os.getcwd() # Gets the current directory (and in which the code is placed)
filename_out = (dir_pwd + '/Results/AllResults_'+GCM+'_500yr.out') # Directory to save processed data

my_shelf = shelve.open(filename_out)
#globals()[key]=my_shelf[key]
globals()['WS_indeces_lonlat']=my_shelf['WS_indeces_lonlat']
globals()['LAB_indeces_lonlat']=my_shelf['LAB_indeces_lonlat']
globals()['WS_index_norm']=my_shelf['WS_index_norm']
globals()['WS_index_norm_rm']=my_shelf['WS_index_norm_rm']
globals()['LAB_index_norm']=my_shelf['LAB_index_norm']
globals()['LAB_index_norm_rm']=my_shelf['LAB_index_norm_rm']
my_shelf.close()

##################################################
#%% Convection Index Calculations - WS and LAB ###

conv_index_depth_ws=0 # args[7] in func_time_depth_plot code - depth for Convection Index
conv_index_depth_lab=0

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

WS_index_time_depth_surface, WS_index_surface=func_time_depth_plot(dset_thetao, year_start, year_end, WS_indeces_lonlat, 90, 180, conv_index_depth_ws) # WS_index is the convection indeces at surface
WS_index_surface_norm=(WS_index_surface-np.nanmean(WS_index_surface))/np.std(WS_index_surface) # Normalized Convection Index

WS_index_surface_rm = copy.deepcopy(WS_index_surface)
WS_index_surface_rm = runningMeanFast(WS_index_surface_rm, 10)#if data is not smoothed
#WS_index_surface_rm=WS_index_surface_rm*(-1) #multiplying by -1 as we have convection when temperature is low(opposite to the index which we calculate)

WS_index_surface_norm_rm = runningMeanFast(WS_index_surface_norm, 10)#if data is not smoothed

#### Labrador Sea Caclulations ####

LAB_index_time_depth_surface, LAB_index_surface=func_time_depth_plot(dset_thetao, year_start, year_end, LAB_indeces_lonlat, 90, 180, conv_index_depth_lab) # LAB_index is the convection indeces at surface
LAB_index_surface_norm=(LAB_index_surface-np.nanmean(LAB_index_surface))/np.std(LAB_index_surface) # Normalized Convection Index

LAB_index_surface_rm = copy.deepcopy(LAB_index_surface)
LAB_index_surface_rm = runningMeanFast(LAB_index_surface_rm, 10)#if data is not smoothed
#LAB_index_surface_rm=LAB_index_surface_rm*(-1) #multiplying by -1 as we have convection when temperature is low(opposite to the index which we calculate)

LAB_index_surface_norm_rm = runningMeanFast(LAB_index_surface_norm, 10)#if data is not smoothed

dset_thetao.close_ncfile(dset_thetao.fin)
dset_so.close_ncfile(dset_so.fin)

###############################################################################
###############################################################################
###############################################################################
MLD_years_LAB_timeseries=copy.deepcopy(MLD_years_LAB)
MLD_years_LAB_timeseries=np.nanmean( np.nanmean(MLD_years_LAB_timeseries,axis=2) ,axis=1)
MLD_years_LAB_timeseries_rm = copy.deepcopy(MLD_years_LAB_timeseries)
MLD_years_LAB_timeseries_rm = runningMeanFast(MLD_years_LAB_timeseries_rm, 10)#if data is not smoothed

P_title='Peak Labrador Sea Convection lagged correlation - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig=plt.figure()
func_plot_lagcor_sig(LAB_index_rm[:-9], MLD_years_LAB_timeseries_rm[:-9], 40, P_title, 'b', 'MLD', 'best', '-')
plt.ylim(-1,1)
#fig.savefig(dir_figs+str(GCM)+'_lagcor_LABconvection_MLD_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


MLD_years_WS_timeseries=copy.deepcopy(MLD_years_WS)
MLD_years_WS_timeseries=np.nanmean( np.nanmean(MLD_years_WS_timeseries,axis=2) ,axis=1)
MLD_years_WS_timeseries_rm = copy.deepcopy(MLD_years_WS_timeseries)
MLD_years_WS_timeseries_rm = runningMeanFast(MLD_years_WS_timeseries_rm, 10)#if data is not smoothed

P_title='Peak Weddell Sea Convection lagged correlation - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig=plt.figure()
func_plot_lagcor_sig(WS_index_rm[:-9], MLD_years_WS_timeseries_rm[:-9], 40, P_title, 'b', 'MLD', 'best', '-')
plt.ylim(-1,1)

#fig.savefig(dir_figs+str(GCM)+'_lagcor_WSconvection_MLD_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')



###############################################################################
filename_out = (dir_pwd + '/Results/'+GCM+'_SurfaceConv_500yr.out') # Directory to save processed data
#################        To Save:        #####################################
import os
import shelve
my_shelf = shelve.open(filename_out,'n') # 'n' for new
var_list_1=['WS_index_time_depth_surface','WS_index_surface','WS_index_surface_norm','WS_index_surface_rm','WS_index_surface_norm_rm']
for key in var_list_1:
    try:
        my_shelf[key] = globals()[key]
    except TypeError:
        print('ERROR shelving: {0}'.format(key))
my_shelf.close()
#################        To restore:        ##################################
filename_out = (dir_pwd + '/Results/'+GCM+'_SurfaceConv_500yr.out') # Directory to save processed data

my_shelf = shelve.open(filename_out)
for key in my_shelf:
    globals()[key]=my_shelf[key]
my_shelf.close()
###############################################################################
#########################################

years=np.linspace(year_start,year_end,year_end-year_start+1)

P_Var_x1=P_Var_x2=years
P_Var_y1=copy.deepcopy(WS_index_norm)
P_Var_y2=copy.deepcopy(WS_index_surface_norm_rm)
P_xlable='Years' ;  P_ylable='Convection index (°C)'
P_title='Weddel Sea convection index (normalized, smoothed) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plotline_2var(P_Var_x1, P_Var_y1, P_Var_x2, P_Var_y2, P_xlable, P_ylable, P_title, 'b', 'r', 'Conv. Index at 500m depth', 'Conv. Index at surface', 'best', '-')
for l in fig.gca().lines:
    l.set_alpha(0.7)
l = plt.axhline(y=-0.5, color='k', linestyle ='--')
l = plt.axhline(y=0, color='k')
l = plt.axhline(y=0.5, color='k', linestyle ='--')
plt.text(10, 1.1, 'Non-Conv.', horizontalalignment='center', verticalalignment='center', fontsize=20, color='b')
plt.text(0, -1.1, 'Conv.', horizontalalignment='center', verticalalignment='center', fontsize=20, color='b')
fig.savefig(dir_figs+str(GCM)+'_convec_index_WS_2.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


P_Var_x1=P_Var_x2=years
P_Var_y1=copy.deepcopy(WS_index_norm)
P_Var_y2=copy.deepcopy(WS_index_surface_norm_rm*-1)
P_xlable='Years' ;  P_ylable='Convection index (°C)'
P_title='Weddel Sea convection index (normalized, smoothed) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plotline_2var(P_Var_x1, P_Var_y1, P_Var_x2, P_Var_y2, P_xlable, P_ylable, P_title, 'b', 'g', 'Conv. Index at 500m depth', 'Conv. Index at surface - Inverted', 'best', '-')
for l in fig.gca().lines:
    l.set_alpha(0.7)
l = plt.axhline(y=-0.5, color='k', linestyle ='--')
l = plt.axhline(y=0, color='k')
l = plt.axhline(y=0.5, color='k', linestyle ='--')
plt.text(10, 1.1, 'Non-Conv.', horizontalalignment='center', verticalalignment='center', fontsize=20, color='b')
plt.text(0, -1.1, 'Conv.', horizontalalignment='center', verticalalignment='center', fontsize=20, color='b')
fig.savefig(dir_figs+str(GCM)+'_convec_index_WS_3.png', format='png', dpi=300, transparent=True, bbox_inches='tight')



P_Var_x1=P_Var_x2=years
P_Var_y1=copy.deepcopy(LAB_index_norm)
P_Var_y2=copy.deepcopy(LAB_index_surface_norm_rm)
P_xlable='Years' ;  P_ylable='Convection index (°C)'
P_title='Labrador Sea convection index (normalized, smoothed) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plotline_2var(P_Var_x1, P_Var_y1, P_Var_x2, P_Var_y2, P_xlable, P_ylable, P_title, 'b', 'r', 'Conv. Index at 500m depth', 'Conv. Index at surface', 'best', '-')
for l in fig.gca().lines:
    l.set_alpha(0.7)
l = plt.axhline(y=-0.5, color='k', linestyle ='--')
l = plt.axhline(y=0, color='k')
l = plt.axhline(y=0.5, color='k', linestyle ='--')
plt.text(10, 1, 'Non-Conv.', horizontalalignment='center', verticalalignment='center', fontsize=20, color='b')
plt.text(0, -1, 'Conv.', horizontalalignment='center', verticalalignment='center', fontsize=20, color='b')
fig.savefig(dir_figs+str(GCM)+'_convec_index_LAB_2.png', format='png', dpi=300, transparent=True, bbox_inches='tight')



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


North_Westerlies=Zonal_winds[:,135:145]

AMOC_max_50S= AMOC_max[:,40] # AMOC at 50S # Max of streamfunction method
AMOC_max_50S_m=runningMeanFast(AMOC_max_50S, 10)
AMOC_max_30S= AMOC_max[:,60] # AMOC at 30S # Max of streamfunction method
AMOC_max_30S_m=runningMeanFast(AMOC_max_30S, 10)
AMOC_max_50N= AMOC_max[:,140] # AMOC at 50N # Max of streamfunction method
AMOC_max_50N_m=runningMeanFast(AMOC_max_50N, 10)
AMOC_max_30N= AMOC_max[:,120] # AMOC at 30N # Max of streamfunction method
AMOC_max_30N_m=runningMeanFast(AMOC_max_30N, 10)


P_title='Peak WS Convection lagged correlation (Convection index=Temp at Suface) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig=plt.figure()
func_plot_lagcor_sig(WS_index_surface_rm[:-9], LAB_index_rm[:-9], 40, P_title, 'y', 'LAB index (Temp at 500m depth)', 'best', 'yes')
func_plot_lagcor_sig(WS_index_surface_rm[:-9], Winds_40S60S_0W30W_m[:-9], 40, P_title, 'g', 'Zonal Winds [40S60S, 0W30W]', 'best', 'yes')
func_plot_lagcor_sig(WS_index_surface_rm[:-9], AMOC_max_50S_m[:-9], 40, P_title, 'b', 'AMOC max 50S', 'best', 'yes')
func_plot_lagcor_sig(WS_index_surface_rm[:-9], AMOC_max_30S_m[:-9], 40, P_title, 'r', 'AMOC max 30S', 'best', 'yes')
plt.ylim(-0.7,0.7)
#fig.savefig(dir_figs+str(GCM)+'_lagcor_WSconvection_Final_Surface_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

P_title='Peak WS Convection lagged correlation (Convection index=Temp at Suface) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig=plt.figure()
func_plot_lagcor_sig(WS_index_surface_rm[:-9], LAB_index_rm[:-9], 40, P_title, 'y', 'LAB index (Temp at 500m depth)', 'best', '-')
func_plot_lagcor_sig(WS_index_surface_rm[:-9], Winds_40S60S_0W30W_m[:-9], 40, P_title, 'g', 'Zonal Winds [40S60S, 0W30W]', 'best', '-')
func_plot_lagcor_sig(WS_index_surface_rm[:-9], AMOC_max_50S_m[:-9], 40, P_title, 'b', 'AMOC max 50S', 'best', '-')
func_plot_lagcor_sig(WS_index_surface_rm[:-9], AMOC_max_30S_m[:-9], 40, P_title, 'r', 'AMOC max 30S', 'best', '-')
plt.ylim(-0.7,0.7)
#fig.savefig(dir_figs+str(GCM)+'_lagcor_WSconvection_Final_Surface_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


P_title='Peak WS Convection index at 500m depth lagged correlation - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig=plt.figure()
func_plot_lagcor_sig(WS_index_norm_rm,WS_index_surface_norm_rm[:-9], 40, P_title, 'b', 'WS Convection index at Surface', 'best', '-')
plt.ylim(-0.7,0.7)
#fig.savefig(dir_figs+str(GCM)+'_lagcor_WSconvection_Surf_500m.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


P_title='Peak LAB Convection index at 500m depth lagged correlation - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig=plt.figure()
func_plot_lagcor_sig(LAB_index_norm_rm,LAB_index_surface_norm_rm[:-9], 40, P_title, 'b', 'LAB Convection index at Surface', 'best', '-')
plt.ylim(-0.4,1)
#fig.savefig(dir_figs+str(GCM)+'_lagcor_LABconvection_Surf_500m.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


############################################################
########       Sea level pressure composites      ##########
t_frequency='Amon'
variable='psl'

dir_data_in1 = ('/data5/scratch/CMIP5/CMIP5_models/atmosphere_physics/') # Directory to raed raw data from
dir_data_in2=(dir_data_in1+ GCM + '/piControl/mo/')

dset_psl= netcdf_read(dir_data_in2+str(variable)+'_'+str(t_frequency)+'_'+str(GCM)+'*12.nc',str(GCM),variable)
start_date_i,end_date_i = dset_psl.find_time(dset_psl.times, year_start, year_end)

PSL_final=[]

for i in range(int((end_date_i+1-start_date_i)/12)):
    #print (int((end_date_i+1-start_date_i)/12))
    print('PSL Composte - Year: ', i)
    data_vo_extracted=dset_psl.extract_data(dset_psl.variable,start_date_i+12*i,start_date_i+12*i+11)
    data=np.squeeze(data_vo_extracted)
    data=np.mean(data, axis=0)
    #x_i, y_i = np.meshgrid(dset_ts.y,dset_ts.x)
    #lon,lat,data_i=interpolate_2_reg_grid(dset_ts.x,dset_ts.y,data)
    data_i = func_regrid(data, dset_psl.y, dset_psl.x, Lat_regrid_2D, Lon_regrid_2D)
    data_i[data_i>10000000]=np.nan
    PSL_final.append(data_i)

PSL_final=np.asarray(PSL_final)
dset_psl.close_ncfile(dset_psl.fin)

R_ii_LAB = np.where(LAB_index_norm_rm <-0.5)
R_jj_LAB = np.where(LAB_index_norm_rm >0.5)

PSL_final_LAB_c=PSL_final[R_ii_LAB]
PSL_final_LAB_n=PSL_final[R_jj_LAB]
conv_LAB=np.nanmean(PSL_final_LAB_c,axis=0)
nonconv_LAB=np.nanmean(PSL_final_LAB_n,axis=0)

composite_LAB=np.subtract(conv_LAB,nonconv_LAB)

Plot_Var = composite_LAB
cmap_limit=140 #np.nanmax(np.abs( np.nanpercentile(Plot_Var, 99.9)))
Plot_range=np.linspace(-cmap_limit,cmap_limit,141) ### Or:  Plot_range=100
Plot_unit='(Pa)'; Plot_title= 'Sea Level Pressure composites - Labrador Sea convection - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig, m = func_plotmap_contourf(Plot_Var, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.seismic, 'cyl', 210., 80., -80., '-')
fig.savefig(dir_figs+str(GCM)+'_Sea_level_pressure_composites_LAB.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


R_ii_WS = np.where(WS_index_norm_rm <-0.5)
R_jj_WS = np.where(WS_index_norm_rm >0.5)

PSL_final_WS_c=PSL_final[R_ii_WS]
PSL_final_WS_n=PSL_final[R_jj_WS]
conv_WS=np.nanmean(PSL_final_WS_c,axis=0)
nonconv_WS=np.nanmean(PSL_final_WS_n,axis=0)

composite_WS=np.subtract(conv_WS,nonconv_WS)

Plot_Var = composite_WS
cmap_limit=140 #np.nanmax(np.abs( np.nanpercentile(Plot_Var, 99.9)))
Plot_range=np.linspace(-cmap_limit,cmap_limit,141) ### Or:  Plot_range=100
Plot_unit='(Pa)'; Plot_title= 'Sea Level Pressure composites - Weddell Sea convection - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig, m = func_plotmap_contourf(Plot_Var, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.seismic, 'cyl', 210., 80., -80., '-')
fig.savefig(dir_figs+str(GCM)+'_Sea_level_pressure_composites_WS.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


#########################################
########       Heat Flux       ##########
t_frequency='Omon'
variable='hfds'

dir_data_in1 = ('/data5/scratch/CMIP5/CMIP5_models/ocean_physics/') # Directory to raed raw data from
dir_data_in2=(dir_data_in1+ GCM + '/piControl/mo/')

dset_hfds= netcdf_read(dir_data_in2+str(variable)+'_'+str(t_frequency)+'_'+str(GCM)+'*12.nc',str(GCM),variable)
start_date_i,end_date_i = dset_hfds.find_time(dset_hfds.times, year_start, year_end)

HFDS_final=[]

for i in range(int((end_date_i+1-start_date_i)/12)):
    #print (int((end_date_i+1-start_date_i)/12))
    print('HFDS Composte - Year: ', i)
    data_vo_extracted=dset_hfds.extract_data(dset_hfds.variable,start_date_i+12*i,start_date_i+12*i+11)
    data=np.squeeze(data_vo_extracted)
    data=np.mean(data, axis=0)
    #x_i, y_i = np.meshgrid(dset_ts.y,dset_ts.x)
    #lon,lat,data_i=interpolate_2_reg_grid(dset_ts.x,dset_ts.y,data)
    data_i = func_regrid(data, dset_hfds.y, dset_hfds.x, Lat_regrid_2D, Lon_regrid_2D)
    data_i[data_i>10000000]=np.nan
    HFDS_final.append(data_i)

HFDS_final=np.asarray(HFDS_final)
dset_hfds.close_ncfile(dset_hfds.fin)


Plot_Var = copy.deepcopy(np.nanmean(HFDS_final, axis=0))
#Plot_Var[ Ocean_Land_mask==0 ]=nan # masking over land, so grid cells that fall on land area (value=0) will be deleted
Plot_range=51 #np.linspace(1018,1028,41)
Plot_unit='(W/m2)'; Plot_title= 'Heat Flux into Ocean average map (W/m2) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig, m = func_plotmap_contourf(Plot_Var, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.jet, 'cyl', 210., 80., -80., '-')
plt.show()
#fig.savefig(dir_figs+str(GCM)+'_heat_flux_cave_map.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


R_ii_LAB = np.where(LAB_index_norm_rm <-0.5)
R_jj_LAB = np.where(LAB_index_norm_rm >0.5)

HFDS_final_LAB_c=HFDS_final[R_ii_LAB]
HFDS_final_LAB_n=HFDS_final[R_jj_LAB]
conv_LAB=np.nanmean(HFDS_final_LAB_c,axis=0)
nonconv_LAB=np.nanmean(HFDS_final_LAB_n,axis=0)

composite_LAB=np.subtract(conv_LAB,nonconv_LAB)


Plot_Var = composite_LAB
cmap_limit=30 #np.nanmax(np.abs( np.nanpercentile(Plot_Var, 99.9)))
Plot_range=np.linspace(-cmap_limit,cmap_limit,61) ### Or:  Plot_range=100
Plot_unit='(W/m2)'; Plot_title= 'Heat Flux into Ocean composites - Labrador Sea convection - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig, m = func_plotmap_contourf(Plot_Var, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.seismic, 'cyl', 210., 80., -80., '-')
fig.savefig(dir_figs+str(GCM)+'_heat_flux_composites_LAB.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


R_ii_WS = np.where(WS_index_norm_rm <-0.5)
R_jj_WS = np.where(WS_index_norm_rm >0.5)

HFDS_final_WS_c=HFDS_final[R_ii_WS]
HFDS_final_WS_n=HFDS_final[R_jj_WS]
conv_WS=np.nanmean(HFDS_final_WS_c,axis=0)
nonconv_WS=np.nanmean(HFDS_final_WS_n,axis=0)

composite_WS=np.subtract(conv_WS,nonconv_WS)

Plot_Var = composite_WS
cmap_limit=30 #np.nanmax(np.abs( np.nanpercentile(Plot_Var, 99.9)))
Plot_range=np.linspace(-cmap_limit,cmap_limit,61) ### Or:  Plot_range=100
Plot_unit='(W/m2)'; Plot_title= 'Heat Flux into Ocean composites - Weddell Sea convection - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig, m = func_plotmap_contourf(Plot_Var, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.seismic, 'cyl', 210., 80., -80., '-')
fig.savefig(dir_figs+str(GCM)+'_heat_flux_composites_WS.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

###############
HFDS_final_WS_mean = np.nanmean (HFDS_final [:,0:61,291:] , axis =2)
HFDS_final_WS_mean = np.nanmean (HFDS_final_WS_mean , axis =1)
HFDS_final_WS_mean_rm = copy.deepcopy(HFDS_final_WS_mean)
HFDS_final_WS_mean_rm = runningMeanFast(HFDS_final_WS_mean_rm, 10)#if data is not smoothed

years=np.linspace(year_start,year_end,year_end-year_start+1)
P_Var_x=years[:-9]
P_Var_y=HFDS_final_WS_mean_rm[:-9]
P_xlable='Years'
P_ylable='Heat Flux into Ocean (W/m2)'
P_title='Heat Flux into Ocean average (Weddell Sea [70W-0W, 90S-30S], smoothed) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plotline_1var(P_Var_x, P_Var_y, P_xlable, P_ylable, P_title, 'darkcyan', '-')
fig.savefig(dir_figs+str(GCM)+'_heat_flux_ave_timeseries_WS.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

P_title='Peak Weddel Sea Convection lagged correlation - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig=plt.figure()
func_plot_lagcor_sig(WS_index_rm[:-9], HFDS_final_WS_mean_rm[:-9], 40, P_title, 'r', 'Heat Flux into Ocean over Weddell Sea [70W-0W, 90S-30S]', 'best', '-')
plt.ylim(-1,1)
fig.savefig(dir_figs+str(GCM)+'_lagcor_WSconvection_HFDS_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

HFDS_final_LAB_mean = np.nanmean (HFDS_final [:,140:156,295:320] , axis =2)
HFDS_final_LAB_mean = np.nanmean (HFDS_final_LAB_mean , axis =1)
HFDS_final_LAB_mean_rm = copy.deepcopy(HFDS_final_LAB_mean)
HFDS_final_LAB_mean_rm = runningMeanFast(HFDS_final_LAB_mean_rm, 10)#if data is not smoothed

years=np.linspace(year_start,year_end,year_end-year_start+1)
P_Var_x=years[:-9]
P_Var_y=HFDS_final_LAB_mean_rm[:-9]
P_xlable='Years'
P_ylable='Heat Flux into Ocean (W/m2)'
P_title='Heat Flux into Ocean average (Labrador Sea [65W-41W, 50N-65N], smoothed) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plotline_1var(P_Var_x, P_Var_y, P_xlable, P_ylable, P_title, 'darkcyan', '-')
fig.savefig(dir_figs+str(GCM)+'_heat_flux_ave_timeseries_LAB.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

P_title='Peak Labrador Sea Convection lagged correlation - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig=plt.figure()
func_plot_lagcor_sig(LAB_index_rm[:-9], HFDS_final_LAB_mean_rm[:-9], 40, P_title, 'r', 'Heat Flux into Ocean over Labrador Sea [65W-41W, 50N-65N]', 'best', '-')
plt.ylim(-1,1)
fig.savefig(dir_figs+str(GCM)+'_lagcor_LABconvection_HFDS_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
##############

P_Var=copy.deepcopy(HFDS_final)
P_peaktime=387  ; P_lag=40; P_lag_period=5; P_anomalies='yes'; P_range=np.linspace(-60,60,61) 
P_title='Heat Flux into Ocean anomalies during the WS convection at year '+str(P_peaktime)+' (W/m2) - '+str(GCM)
fig = func_plot_laggedmaps(P_Var, Lon_regrid_2D, Lat_regrid_2D, P_peaktime, P_lag, P_lag_period, P_anomalies, P_range, P_title, plt.cm.seismic, 'mill', 210., 80., -80., 'fill')
#fig.savefig(dir_figs+str(GCM)+'_laggedmaps_WSconvection_HFDS_yr387.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

P_Var=copy.deepcopy(HFDS_final)
P_peaktime=387  ; P_lag=40; P_lag_period=5; P_anomalies='yes'; P_range=np.linspace(-60,60,61)
P_title='Heat Flux into Ocean anomalies during the WS convection at year '+str(P_peaktime)+' - '+str(GCM)
fig = func_plot_laggedmaps(P_Var, Lon_regrid_2D, Lat_regrid_2D, P_peaktime, P_lag, P_lag_period, P_anomalies, P_range, P_title, plt.cm.seismic, 'spstere', 210., 80., -80., 'fill')
#fig.savefig(dir_figs+str(GCM)+'_laggedmaps_WSconvection_HFDS_yr387_S.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


P_Var=copy.deepcopy(HFDS_final)
P_peaktime=367  ; P_lag=40; P_lag_period=5; P_anomalies='yes'; P_range=np.linspace(-60,60,61) 
P_title='Heat Flux into Ocean anomalies during the LAB convection at year '+str(P_peaktime)+' (W/m2) - '+str(GCM)
fig = func_plot_laggedmaps(P_Var, Lon_regrid_2D, Lat_regrid_2D, P_peaktime, P_lag, P_lag_period, P_anomalies, P_range, P_title, plt.cm.seismic, 'mill', 210., 80., -80., 'fill')
#fig.savefig(dir_figs+str(GCM)+'_laggedmaps_LABconvection_HFDS_yr367.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

P_Var=copy.deepcopy(HFDS_final)
P_peaktime=367  ; P_lag=40; P_lag_period=5; P_anomalies='yes'; P_range=np.linspace(-60,60,61)
P_title='Heat Flux into Ocean anomalies during the LAB convection at year '+str(P_peaktime)+' - '+str(GCM)
fig = func_plot_laggedmaps(P_Var, Lon_regrid_2D, Lat_regrid_2D, P_peaktime, P_lag, P_lag_period, P_anomalies, P_range, P_title, plt.cm.seismic, 'npstere', 210., 80., -80., 'fill')
#fig.savefig(dir_figs+str(GCM)+'_laggedmaps_LABconvection_HFDS_yr367_N.png', format='png', dpi=300, transparent=True, bbox_inches='tight')



###############################################################################
###############################################################################
########                  Lagged maps - annomalies                  ###########
###############################################################################
###############################################################################
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
my_shelf.close()

filename_out = (dir_pwd + '/Results/'+GCM+'_SurfaceConv_500yr.out') # Directory to save processed data

my_shelf = shelve.open(filename_out)
for key in my_shelf:
    globals()[key]=my_shelf[key]
my_shelf.close()
########                  Temperature                  ###########
Temp_allyears=np.load('Temp_allyears_GFDL-ESM2G_500yr.npy')

P_Var=copy.deepcopy(Temp_allyears[:,0,:,:]) ; P_Var=P_Var - 273.15
P_peaktime=387  ; P_lag=40; P_lag_period=5; P_anomalies='yes'; P_range=np.linspace(-2.2,2.2,101) 
P_title='Sea Surface Temperature anomalies during the WS convection at year '+str(P_peaktime)+' (degree C) - '+str(GCM)
fig = func_plot_laggedmaps(P_Var, Lon_regrid_2D, Lat_regrid_2D, P_peaktime, P_lag, P_lag_period, P_anomalies, P_range, P_title, plt.cm.seismic, 'mill', 210., 80., -80., 'fill')
#fig.savefig(dir_figs+str(GCM)+'_laggedmaps_WSconvection_SST_yr387.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


P_Var=copy.deepcopy(Temp_allyears[:,0,:,:]) ; P_Var=P_Var - 273.15
P_peaktime=367  ; P_lag=40; P_lag_period=5; P_anomalies='yes'; P_range=np.linspace(-2.2,2.2,101) 
P_title='Sea Surface Temperature anomalies during the LAB convection at year '+str(P_peaktime)+' (degree C) - '+str(GCM)
fig = func_plot_laggedmaps(P_Var, Lon_regrid_2D, Lat_regrid_2D, P_peaktime, P_lag, P_lag_period, P_anomalies, P_range, P_title, plt.cm.seismic, 'mill', 210., 80., -80., 'fill')
#fig.savefig(dir_figs+str(GCM)+'_laggedmaps_LABconvection_SST_yr367.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


#P_Var=copy.deepcopy(Temp_allyears[:,0:35,:,:])
P_Var=np.nanmean(Temp_allyears[:,0:35,:,:], axis=1); P_Var=P_Var - 273.15
P_peaktime=387  ; P_lag=40; P_lag_period=5; P_anomalies='yes'; P_range=np.linspace(-2.2,2.2,101) 
P_title='First 1000m Temperature anomalies during the WS convection at year '+str(P_peaktime)+' (degree C) - '+str(GCM)
fig = func_plot_laggedmaps(P_Var, Lon_regrid_2D, Lat_regrid_2D, P_peaktime, P_lag, P_lag_period, P_anomalies, P_range, P_title, plt.cm.seismic, 'mill', 210., 80., -80., 'fill')
#fig.savefig(dir_figs+str(GCM)+'_laggedmaps_WSconvection_Temp_0_1000_yr387.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


#P_Var=copy.deepcopy(Temp_allyears[:,0:35,:,:])
P_Var=np.nanmean(Temp_allyears[:,0:35,:,:], axis=1); P_Var=P_Var - 273.15
P_peaktime=367  ; P_lag=40; P_lag_period=5; P_anomalies='yes'; P_range=np.linspace(-2.2,2.2,101) 
P_title='First 1000m Temperature anomalies during the LAB convection at year '+str(P_peaktime)+' (degree C) - '+str(GCM)
fig = func_plot_laggedmaps(P_Var, Lon_regrid_2D, Lat_regrid_2D, P_peaktime, P_lag, P_lag_period, P_anomalies, P_range, P_title, plt.cm.seismic, 'mill', 210., 80., -80., 'fill')
#fig.savefig(dir_figs+str(GCM)+'_laggedmaps_LABconvection_Temp_0_1000_yr367.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


########                  Salinity                  ###########
Salinity_allyears=np.load('Salinity_allyears_GFDL-ESM2G_500yr.npy')

P_Var=copy.deepcopy(Salinity_allyears[:,0,:,:])
P_peaktime=387  ; P_lag=40; P_lag_period=5; P_anomalies='yes'; P_range=np.linspace(-1.4,1.4,101) 
P_title='Sea Surface Salinity anomalies during the WS convection at year '+str(P_peaktime)+' (PSU) - '+str(GCM)
fig = func_plot_laggedmaps(P_Var, Lon_regrid_2D, Lat_regrid_2D, P_peaktime, P_lag, P_lag_period, P_anomalies, P_range, P_title, plt.cm.seismic, 'mill', 210., 80., -80., 'fill')
fig.savefig(dir_figs+str(GCM)+'_laggedmaps_WSconvection_Salinity_yr387.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


P_Var=copy.deepcopy(Salinity_allyears[:,0,:,:])
P_peaktime=367  ; P_lag=40; P_lag_period=5; P_anomalies='yes'; P_range=np.linspace(-1.4,1.4,101) 
P_title='Sea Surface Salinity anomalies during the LAB convection at year '+str(P_peaktime)+' (PSU) - '+str(GCM)
fig = func_plot_laggedmaps(P_Var, Lon_regrid_2D, Lat_regrid_2D, P_peaktime, P_lag, P_lag_period, P_anomalies, P_range, P_title, plt.cm.seismic, 'mill', 210., 80., -80., 'fill')
fig.savefig(dir_figs+str(GCM)+'_laggedmaps_LABconvection_Salinity_yr367.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


#P_Var=copy.deepcopy(Salinity_allyears[:,0:35,:,:])
P_Var=np.nanmean(Salinity_allyears[:,0:35,:,:], axis=1)
P_peaktime=387  ; P_lag=40; P_lag_period=5; P_anomalies='yes'; P_range=np.linspace(-0.6,0.6,101) 
P_title='First 1000m Salinity anomalies during the WS convection at year '+str(P_peaktime)+' (PSU) - '+str(GCM)
fig = func_plot_laggedmaps(P_Var, Lon_regrid_2D, Lat_regrid_2D, P_peaktime, P_lag, P_lag_period, P_anomalies, P_range, P_title, plt.cm.seismic, 'mill', 210., 80., -80., 'fill')
#fig.savefig(dir_figs+str(GCM)+'_laggedmaps_WSconvection_Salinity_0_1000_yr387.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


#P_Var=copy.deepcopy(Salinity_allyears[:,0:35,:,:])
P_Var=np.nanmean(Salinity_allyears[:,0:35,:,:], axis=1)
P_peaktime=367  ; P_lag=40; P_lag_period=5; P_anomalies='yes'; P_range=np.linspace(-0.6,0.6,101) 
P_title='First 1000m Salinity anomalies during the LAB convection at year '+str(P_peaktime)+' (PSU) - '+str(GCM)
fig = func_plot_laggedmaps(P_Var, Lon_regrid_2D, Lat_regrid_2D, P_peaktime, P_lag, P_lag_period, P_anomalies, P_range, P_title, plt.cm.seismic, 'mill', 210., 80., -80., 'fill')
#fig.savefig(dir_figs+str(GCM)+'_laggedmaps_LABconvection_Salinity_0_1000_yr367.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


########                  Sea Ice Concentration                  ###########
t_frequency='OImon'
variable_sic='sic' # Sea Ice Concentration

dir_data_in1 = ('/data5/scratch/CMIP5/CMIP5_models/ocean_physics/') # Directory to raed raw data from
dir_data_in2=(dir_data_in1+ GCM + '/piControl/mo/')

dset_sic = netcdf_read (dir_data_in2+str(variable_sic)+'_'+str(t_frequency)+'_'+str(GCM)+'*12.nc',str(GCM), variable_sic) # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
start_date_i,end_date_i = dset_sic.find_time(dset_sic.times, year_start, year_end)

#SIC_Annual_allyears = np.zeros(( int((end_date_i+1-start_date_i)/12), Lat_regrid_2D.shape[0], Lon_regrid_2D.shape[1]))
SIC_Annual_allyears=np.full([int((end_date_i+1-start_date_i)/12),Lat_regrid_2D.shape[0],Lon_regrid_2D.shape[1]], np.nan)
for t in range(int((end_date_i+1-start_date_i)/12)):
    if t%10==0: # Just showing every 10 year
        print('SIC calc - Year: ', year_start+t)
    data_sic_extracted=dset_sic.extract_data(dset_sic.variable, start_date_i+12*t+1-1,start_date_i+12*t+12-1)
    data_sic_extracted=np.nanmean(data_sic_extracted, axis=0)  
    data_sic_extracted=np.squeeze(data_sic_extracted)
    data_i = func_regrid(data_sic_extracted, dset_sic.y, dset_sic.x, Lat_regrid_2D, Lon_regrid_2D)
    data_i[data_i>1e19]=np.nan
    
    SIC_Annual_allyears[t,:,:]=data_i


month_no=9
SIC_September_allyears=np.full([int((end_date_i+1-start_date_i)/12),Lat_regrid_2D.shape[0],Lon_regrid_2D.shape[1]], np.nan)
for t in range(int((end_date_i+1-start_date_i)/12)):
    if t%10==0: # Just showing every 10 year
        print('SIC calc - Year: ', year_start+t)
    data_sic_extracted=dset_sic.extract_data(dset_sic.variable, start_date_i+12*t+month_no-1,start_date_i+12*t+month_no-1)
    data_sic_extracted=np.squeeze(data_sic_extracted)
    data_i = func_regrid(data_sic_extracted, dset_sic.y, dset_sic.x, Lat_regrid_2D, Lon_regrid_2D)
    data_i[data_i>1e19]=np.nan
    
    SIC_September_allyears[t,:,:]=data_i


month_no=3
SIC_March_allyears=np.full([int((end_date_i+1-start_date_i)/12),Lat_regrid_2D.shape[0],Lon_regrid_2D.shape[1]], np.nan)
for t in range(int((end_date_i+1-start_date_i)/12)):
    if t%10==0: # Just showing every 10 year
        print('SIC calc - Year: ', year_start+t)
    data_sic_extracted=dset_sic.extract_data(dset_sic.variable, start_date_i+12*t+month_no-1,start_date_i+12*t+month_no-1)
    data_sic_extracted=np.squeeze(data_sic_extracted)
    data_i = func_regrid(data_sic_extracted, dset_sic.y, dset_sic.x, Lat_regrid_2D, Lon_regrid_2D)
    data_i[data_i>1e19]=np.nan
    
    SIC_March_allyears[t,:,:]=data_i


P_Var=copy.deepcopy(SIC_September_allyears)
P_peaktime=387  ; P_lag=40; P_lag_period=5; P_anomalies='no'; P_range=np.linspace(0,100,51) 
P_title='September Sea Ice Area Fraction (SIC) during the WS convection at year '+str(P_peaktime)+' - '+str(GCM)
fig = func_plot_laggedmaps(P_Var, Lon_regrid_2D, Lat_regrid_2D, P_peaktime, P_lag, P_lag_period, P_anomalies, P_range, P_title, plt.cm.jet, 'mill', 210., 80., -80., 'fill')
fig.savefig(dir_figs+str(GCM)+'_laggedmaps_WSconvection_SIC_September_yr387.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


P_Var=copy.deepcopy(SIC_September_allyears)
P_peaktime=387  ; P_lag=40; P_lag_period=5; P_anomalies='no'; P_range=np.linspace(0,100,51) 
P_title='September Sea Ice Area Fraction (SIC) during the WS convection at year '+str(P_peaktime)+' - '+str(GCM)
fig = func_plot_laggedmaps(P_Var, Lon_regrid_2D, Lat_regrid_2D, P_peaktime, P_lag, P_lag_period, P_anomalies, P_range, P_title, plt.cm.jet, 'spstere', 210., 80., -80., 'fill')
fig.savefig(dir_figs+str(GCM)+'_laggedmaps_WSconvection_SIC_September_yr387_S.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


P_Var=copy.deepcopy(SIC_September_allyears)
P_peaktime=387  ; P_lag=34; P_lag_period=2; P_anomalies='no'; P_range=np.linspace(0,100,51) 
P_title='September Sea Ice Area Fraction (SIC) during the WS convection at year '+str(P_peaktime)+' - '+str(GCM)
fig = func_plot_laggedmaps(P_Var, Lon_regrid_2D, Lat_regrid_2D, P_peaktime, P_lag, P_lag_period, P_anomalies, P_range, P_title, plt.cm.jet, 'spstere', 210., 80., -80., 'fill')
fig.savefig(dir_figs+str(GCM)+'_laggedmaps_WSconvection_SIC_September_yr387_S_2.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


P_Var=copy.deepcopy(SIC_March_allyears)
P_peaktime=367  ; P_lag=40; P_lag_period=5; P_anomalies='no'; P_range=np.linspace(0,100,51) 
P_title='March Sea Ice Area Fraction (SIC) during the LAB convection at year '+str(P_peaktime)+' - '+str(GCM)
fig = func_plot_laggedmaps(P_Var, Lon_regrid_2D, Lat_regrid_2D, P_peaktime, P_lag, P_lag_period, P_anomalies, P_range, P_title, plt.cm.jet, 'mill', 210., 80., -80., 'fill')
fig.savefig(dir_figs+str(GCM)+'_laggedmaps_LABconvection_SIC_March_yr367.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


P_Var=copy.deepcopy(SIC_March_allyears)
P_peaktime=367  ; P_lag=40; P_lag_period=5; P_anomalies='no'; P_range=np.linspace(0,100,51) 
P_title='March Sea Ice Area Fraction (SIC) during the LAB convection at year '+str(P_peaktime)+' - '+str(GCM)
fig = func_plot_laggedmaps(P_Var, Lon_regrid_2D, Lat_regrid_2D, P_peaktime, P_lag, P_lag_period, P_anomalies, P_range, P_title, plt.cm.jet, 'npstere', 210., 80., -80., 'fill')
fig.savefig(dir_figs+str(GCM)+'_laggedmaps_LABconvection_SIC_March_yr367_N.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


P_Var=copy.deepcopy(SIC_Annual_allyears)
P_peaktime=387  ; P_lag=40; P_lag_period=5; P_anomalies='no'; P_range=np.linspace(0,100,51) 
P_title='Annual average Sea Ice Area Fraction (SIC) during the WS convection at year '+str(P_peaktime)+' - '+str(GCM)
fig = func_plot_laggedmaps(P_Var, Lon_regrid_2D, Lat_regrid_2D, P_peaktime, P_lag, P_lag_period, P_anomalies, P_range, P_title, plt.cm.jet, 'mill', 210., 80., -80., 'fill')
fig.savefig(dir_figs+str(GCM)+'_laggedmaps_WSconvection_SIC_Annual_yr387.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


P_Var=copy.deepcopy(SIC_Annual_allyears)
P_peaktime=367  ; P_lag=40; P_lag_period=5; P_anomalies='no'; P_range=np.linspace(0,100,51) 
P_title='Annual average Sea Ice Area Fraction (SIC) during the LAB convection at year '+str(P_peaktime)+' - '+str(GCM)
fig = func_plot_laggedmaps(P_Var, Lon_regrid_2D, Lat_regrid_2D, P_peaktime, P_lag, P_lag_period, P_anomalies, P_range, P_title, plt.cm.jet, 'mill', 210., 80., -80., 'fill')
fig.savefig(dir_figs+str(GCM)+'_laggedmaps_LABconvection_SIC_Annual_yr367.png', format='png', dpi=300, transparent=True, bbox_inches='tight')




P_Var=copy.deepcopy(SIC_September_allyears)
P_peaktime=387  ; P_lag=30; P_lag_period=2; P_anomalies='no'; P_range=np.linspace(0,100,51) 
P_title='September Sea Ice Area Fraction (SIC) during the WS convection at year '+str(P_peaktime)+' - '+str(GCM)
fig = func_plot_laggedmaps(P_Var, Lon_regrid_2D, Lat_regrid_2D, P_peaktime, P_lag, P_lag_period, P_anomalies, P_range, P_title, plt.cm.jet, 'spstere', 210., 80., -80., 'fill')
#fig.savefig(dir_figs+str(GCM)+'_laggedmaps_WSconvection_SIC_September_yr387_2.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


###############
SIC_Annual_allyears_WS_mean = np.nanmean (SIC_Annual_allyears [:,0:61,291:] , axis =2)
SIC_Annual_allyears_WS_mean = np.nanmean (SIC_Annual_allyears_WS_mean , axis =1)
SIC_Annual_allyears_WS_mean_rm = copy.deepcopy(SIC_Annual_allyears_WS_mean)
SIC_Annual_allyears_WS_mean_rm = runningMeanFast(SIC_Annual_allyears_WS_mean_rm, 10)#if data is not smoothed

SIC_September_allyears_WS_mean = np.nanmean (SIC_September_allyears [:,0:61,291:] , axis =2)
SIC_September_allyears_WS_mean = np.nanmean (SIC_September_allyears_WS_mean , axis =1)
SIC_September_allyears_WS_mean_rm = copy.deepcopy(SIC_September_allyears_WS_mean)
SIC_September_allyears_WS_mean_rm = runningMeanFast(SIC_September_allyears_WS_mean_rm, 10)#if data is not smoothed

years=np.linspace(year_start,year_end,year_end-year_start+1)
P_Var_x1=P_Var_x2=years[:-9]
P_Var_y1=SIC_Annual_allyears_WS_mean_rm[:-9]
P_Var_y2=SIC_September_allyears_WS_mean_rm[:-9]
P_xlable='Years'
P_ylable='Sea Ice Area Fraction (%)'
P_title='Sea Ice Area Fraction (%) (Weddell Sea [70W-0W, 90S-30S], smoothed) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig = func_plotline_2var(P_Var_x1, P_Var_y1, P_Var_x2, P_Var_y2, P_xlable, P_ylable, P_title, 'darkcyan', 'm', 'Annual average', 'September', 'lower left', '-')
#fig.savefig(dir_figs+str(GCM)+'_SIC_ave_timeseries_WS.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

P_title='Peak Weddel Sea Convection lagged correlation - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
fig=plt.figure()
func_plot_lagcor_sig(WS_index_rm[:-9], SIC_Annual_allyears_WS_mean_rm[:-9], 40, P_title, 'r', 'Annual average Sea Ice Area Fraction over Weddell Sea [70W-0W, 90S-30S]', 'best', '-')
func_plot_lagcor_sig(WS_index_rm[:-9], SIC_September_allyears_WS_mean_rm[:-9], 40, P_title, 'b', 'September Sea Ice Area Fraction over Weddell Sea [70W-0W, 90S-30S]', 'best', '-')
plt.ylim(-1,1)
fig.savefig(dir_figs+str(GCM)+'_lagcor_WSconvection_SIC_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


#############################








