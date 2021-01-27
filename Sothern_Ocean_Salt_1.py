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

filename_out = (dir_pwd + '/Results/AllResults_'+GCM+'_500yr.out') # Directory to save processed data

my_shelf = shelve.open(filename_out)
for key in my_shelf:
    globals()[key]=my_shelf[key]
my_shelf.close()

###############################################################################

def func_salt_transport(*args):
    ##### order of args :
    ##### 0) data (netcdf_read object from CMIP5lib.py),
    ##### 1) data (netcdf_read object from CMIP5lib.py),
    ##### 2) start year,
    ##### 3) end year,
    ##### 4) Mask,
    ##### 5) lat
    ##### 6) lon    
    depths=args[0].extract_depth() # Units: 'm'
    ### these are upper and lower depths of each cell in an ocean grid
    depths_b=args[0].extract_depth_bounds()
    ### calculate the depth of each cell in an ocean grid
    depths_r=depths_b[:,1]-depths_b[:,0]
    ### find timeindeces
    start,end = args[0].find_time(args[0].times, args[2], args[3])
    transport_lon_final=[] # SUM(V_y * dX * dZ) over longitudes, for all years
    transport_final=[]
    transport_0_1000=[]
    transport_1000_2000=[]
    transport_2000_3000=[]
    transport_3000_4000=[]
    transport_4000_5000=[]
    transport_5000_below=[]
    
    Depth_indx=np.zeros((5,4)) # Rows: Depths of 1000, 2000, 3000, 4000 and 5000 meters
                               # Columns: 0=row number, 1=depths, 2=depths lower range, 3=depths upper range    
    for dd in range(depths.shape[0]):
        if (depths_b[dd,0] < 1000) and (depths_b[dd,1] > 1000):
            Depth_indx[0,0]=dd; Depth_indx[0,1]=depths[dd]
            Depth_indx[0,2]=depths_b[dd,0]; Depth_indx[0,3]=depths_b[dd,1]
        if (depths_b[dd,0] < 2000) and (depths_b[dd,1] > 2000):
            Depth_indx[1,0]=dd; Depth_indx[1,1]=depths[dd]
            Depth_indx[1,2]=depths_b[dd,0]; Depth_indx[1,3]=depths_b[dd,1]
        if (depths_b[dd,0] < 3000) and (depths_b[dd,1] > 3000):
            Depth_indx[2,0]=dd; Depth_indx[2,1]=depths[dd]
            Depth_indx[2,2]=depths_b[dd,0]; Depth_indx[2,3]=depths_b[dd,1]
        if (depths_b[dd,0] < 4000) and (depths_b[dd,1] > 4000):
            Depth_indx[3,0]=dd; Depth_indx[3,1]=depths[dd]
            Depth_indx[3,2]=depths_b[dd,0]; Depth_indx[3,3]=depths_b[dd,1]
        if (depths_b[dd,0] < 5000) and (depths_b[dd,1] > 5000):
            Depth_indx[4,0]=dd; Depth_indx[4,1]=depths[dd]
            Depth_indx[4,2]=depths_b[dd,0]; Depth_indx[4,3]=depths_b[dd,1]

    Lat_regrid_1D, Lon_regrid_1D, Lat_bound_regrid, Lon_bound_regrid = func_latlon_regrid_eq(args[5], args[6], -90, 90, 0, 360)
    lon, lat = np.meshgrid(Lon_regrid_1D, Lat_regrid_1D)

    ### we averaging velocities over the year for monthly data
    for t in range(int((end+1-start)/12)):
        print ('Stream calc - Year: ', t)
        data_extracted_vo=args[0].extract_data(args[0].variable,start+12*t,start+12*t+11)
        data_vo=np.squeeze(data_extracted_vo)
        data_vo=np.mean(data_vo, axis=0) # unit: m/s
        data_vo[data_vo > 1e18] = np.nan

        data_extracted_so=args[1].extract_data(args[1].variable,start+12*t,start+12*t+11)
        data_so=np.squeeze(data_extracted_so)
        data_so=np.mean(data_so, axis=0) # unit: gr/kg = gr/l = 1000 gr/m3 = kg/m3    
        data_so[data_so > 1e18] = np.nan
        
        data = np.multiply(data_vo, data_so) # [Unit: kg/m3 * m/s = kg/m2.s ]
        
#        lon,lat,data_i=interpolate_2_reg_grid(args[0].x,args[0].y,data)
        data_i = func_regrid(data, args[0].y, args[0].x, lat, lon)
        
        ##### converting 1e+20 to nan ######
        data_i[data_i>1000]=np.nan
        if t==0: 
            data_depth=np.full([len(lon),len(lon[0])], np.nan)
            data_depth_ranges=np.full([len(data_i),len(data_i[0]),len(data_i[0][0])], np.nan)
        [ii,jj] = args[4]
        for k in range(len(ii)):
            #### I calculate the depth by looking how many nans is in the depth column
            if sum(~np.isnan(data_i[:,ii[k],jj[k]]))>0:
                if t==0:
                    data_depth[ii[k],jj[k]]=depths[sum(~np.isnan(data_i[:,ii[k],jj[k]]))-1]
                    for l in range(sum(~np.isnan(data_i[:,ii[k],jj[k]]))):
                        data_depth_ranges[l,ii[k],jj[k]]=depths_r[l]
                        
        #### calculating volume transport
        #### first multiplying by 111km*cos(lat)
        mul_by_lat=data_i*(np.cos(np.deg2rad(lat))*111321) # S * V_y * dX [Unit: kg/m2.s * m = kg/m.s]
        #### second multiplying by depth
        transport=mul_by_lat*data_depth_ranges/1000000 # S * V_y * dX * dZ [Unit: kg/m.s *m * 1e-6  = kg/s * 10e6] 
        #### calculating integral over dz
        transport_lon=np.nansum(transport,axis=2) # SUM(S * V_y * dX * dZ) over longitudes [Unit: kg/m.s *m * 1e-6  = kg/s * 10e6]
        #### calculating cum integral over dz
        #stream=np.nancumsum(transport_lon,axis=0)
        transport_lon_final.append(transport_lon) # SUM(S * V_y * dX * dZ) over longitudes, for all years
        transport_0_1000.append(np.nanmean(transport[0:int(Depth_indx[0,0])+1,:,:],axis=0))
        transport_1000_2000.append(np.nanmean(transport[int(Depth_indx[0,0])+1:int(Depth_indx[1,0])+1,:,:],axis=0))        
        transport_2000_3000.append(np.nanmean(transport[int(Depth_indx[1,0])+1:int(Depth_indx[2,0])+1,:,:],axis=0))
        transport_3000_4000.append(np.nanmean(transport[int(Depth_indx[2,0])+1:int(Depth_indx[3,0])+1,:,:],axis=0))
        transport_4000_5000.append(np.nanmean(transport[int(Depth_indx[3,0])+1:int(Depth_indx[4,0])+1,:,:],axis=0))
        if depths.shape[0] - (int(Depth_indx[4,0])+1) ==1:
            transport_5000_below.append(transport[int(Depth_indx[4,0])+1,:,:])
        elif depths.shape[0] - (int(Depth_indx[4,0])+1) >=2:
            transport_5000_below.append(np.nanmean(transport[int(Depth_indx[4,0])+1:,:,:],axis=0))

        #transport_final.append(transport)
    transport_lon_final=np.asarray(transport_lon_final)

    transport_0_1000=np.asarray(transport_0_1000)
    transport_1000_2000=np.asarray(transport_1000_2000)      
    transport_2000_3000=np.asarray(transport_2000_3000)
    transport_3000_4000=np.asarray(transport_3000_4000)
    transport_4000_5000=np.asarray(transport_4000_5000)
    transport_5000_below=np.asarray(transport_5000_below)

    transport_0_1000_mean=np.nanmean(transport_0_1000,axis=0)
    transport_2000_3000_mean=np.nanmean(transport_2000_3000,axis=0)

    maxvals=np.nanmax(transport_0_1000_mean,axis=1)
    minvals=np.nanmin(transport_2000_3000_mean,axis=1)

    ind_max = np.array([np.argwhere(transport_0_1000_mean == [x]) for x in maxvals])
    ind_max=np.concatenate(ind_max).astype(None)
    ind_max=np.concatenate(ind_max).astype(None)
    ii_max=ind_max[0::2].astype(int)
    jj_max=ind_max[1::2].astype(int)
    ind_min = np.array([np.argwhere(transport_2000_3000_mean == [x]) for x in minvals])
    ind_min=np.concatenate(ind_min).astype(None)
    ind_min=np.concatenate(ind_min).astype(None)
    ii_min=ind_min[0::2].astype(int)
    jj_min=ind_min[1::2].astype(int)

    return transport_lon_final, transport_0_1000, transport_1000_2000, transport_2000_3000, transport_3000_4000, transport_4000_5000, transport_5000_below, data_depth, Depth_indx, lon, lat, ii_max,jj_max, ii_min, jj_min, transport_0_1000_mean, transport_2000_3000_mean


###############################################################################
mask_atl=calc_Atl_Mask()
t_frequency='Omon'
variable_vo='vo'# vo = sea_water_y_velocity - units: m s-1
variable_so='so' # Sea Water Salinity - units: psu (gr/kg)

dir_data_in1 = ('/data5/scratch/CMIP5/CMIP5_models/ocean_physics/') # Directory to raed raw data from
dir_data_in2=(dir_data_in1+ GCM + '/piControl/mo/')

dset_vo = netcdf_read (dir_data_in2+str(variable_vo)+'_'+str(t_frequency)+'_'+str(GCM)+'*12.nc',str(GCM), variable_vo) # Specifies the netCDF dataset (consisting of multiple files) from which the data will be read
dset_so = netcdf_read (dir_data_in2+str(variable_so)+'_'+str(t_frequency)+'_'+str(GCM)+'*12.nc',str(GCM), variable_so)
Depths=dset_so.lvl[:]


Salt_Transport_lon_final, Salt_transport_0_1000, Salt_transport_1000_2000, Salt_transport_2000_3000, Salt_transport_3000_4000, Salt_transport_4000_5000, Salt_transport_5000_below, Salt_latlon_depths, Salt_Depth_indx, Salt_lon_stream, Salt_lat_stream, Salt_ii_max, Salt_jj_max, Salt_ii_min, Salt_jj_min, Salt_transport_0_1000_mean, Salt_transport_2000_3000_mean  = func_salt_transport(dset_vo, dset_so, year_start, year_end, mask_atl, 180, 360)


Salt_transport_1000_2000_mean=np.nanmean(Salt_transport_1000_2000,axis=0)
Salt_transport_3000_4000_mean=np.nanmean(Salt_transport_3000_4000,axis=0)
Salt_transport_4000_5000_mean=np.nanmean(Salt_transport_4000_5000,axis=0)
Salt_transport_5000_below_mean=np.nanmean(Salt_transport_5000_below,axis=0)

###############################################
### Salt Streamfunction calculations Methode Anna
Salt_Stream_function=empty((Salt_Transport_lon_final.shape[0], Salt_Transport_lon_final.shape[1], Salt_Transport_lon_final.shape[2]))*nan # streamfunction
Salt_Stream_function[:,0,:]=Salt_Transport_lon_final[:,0,:]
for ii in range(1,Salt_Stream_function.shape[1]): # Depths
    Salt_Stream_function[:,ii,:]=Salt_Stream_function[:,ii-1,:]+Salt_Transport_lon_final[:,ii,:]

stfunc_sum=np.nansum(Salt_Transport_lon_final,axis=1)
for ii in range(1,Salt_Stream_function.shape[1]): # Depths
    Salt_Stream_function[:,ii,:]=Salt_Stream_function[:,ii,:] - stfunc_sum

Salt_Stream_function_ave=np.nanmean(Salt_Stream_function, axis=0) # Stream Function averaged over the years
Salt_AMOC_max=np.nanmax(Salt_Stream_function, axis=1)
Salt_SMOC_min=np.nanmin(Salt_Stream_function, axis=1)

#######################################################
Salt_Transport_lon_final_mean=np.nanmean(Salt_Transport_lon_final, axis=0) #  just transport
Salt_Transport_lon_final_mean=np.asarray(Salt_Transport_lon_final_mean)
Salt_Transport_lon_final_mean=np.squeeze(Salt_Transport_lon_final_mean)


###############################################################################
###########                 Saving Results                #####################
var_list_salt=['Salt_Transport_lon_final', 'Salt_transport_0_1000', 'Salt_transport_1000_2000', 'Salt_transport_2000_3000', 'Salt_transport_3000_4000', 'Salt_transport_4000_5000',
                'Salt_transport_5000_below', 'Salt_latlon_depths', 'Salt_Depth_indx', 'Salt_lon_stream', 'Salt_lat_stream', 'Salt_ii_max', 'Salt_jj_max', 'Salt_ii_min', 'Salt_jj_min',
                'Salt_transport_0_1000_mean', 'Salt_transport_2000_3000_mean','Salt_transport_1000_2000_mean','Salt_transport_3000_4000_mean','Salt_transport_4000_5000_mean','Salt_transport_5000_below_mean',
                'Salt_Stream_function','Salt_Stream_function_ave','Salt_AMOC_max','Salt_SMOC_min','Salt_Transport_lon_final_mean','Depths','GCM']
import os
import shelve

dir_pwd = os.getcwd() # Gets the current directory (and in which the code is placed)
filename_out = (dir_pwd + '/Results/AllResults_'+GCM+'_500yr_Salt.out') # Directory to save processed data

### To save
my_shelf = shelve.open(filename_out,'n') # 'n' for new

for key in var_list_salt:
    try:
        my_shelf[key] = globals()[key]
    except TypeError:
        print('ERROR shelving: {0}'.format(key))
my_shelf.close()

##############################################################################
#################        To restore:        ##################################
import os
import shelve

dir_pwd = os.getcwd() # Gets the current directory (and in which the code is placed)
filename_out = (dir_pwd + '/Results/AllResults_'+GCM+'_500yr_Salt.out') # Directory to save processed data

my_shelf = shelve.open(filename_out)
for key in my_shelf:
    globals()[key]=my_shelf[key]
my_shelf.close()
##############################################################################


#######################################
#%% Salt Streamfunction Plots ####
dir_figs = (dir_pwd + '/Figures1/') # Directory to save figures
m = Basemap( projection='mill',lon_0=210)


Plot_title=('Mean northward salt transport upper 1000m [1e6 kg/s] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM))
Plot_save_dir=(dir_figs+str(GCM)+'_salt_transport_0_1000_mean.png')
#bounds_max=np.nanpercentile(transport_0_1000_mean, 99.99)
bounds_max=float("{0:.02f}".format(np.nanpercentile(Salt_transport_0_1000_mean, 98))) # Upper bound of plotted values to be used for colorbar, which is 98th percentile of data, with 2 decimals
bounds = np.arange(-1*bounds_max, bounds_max, bounds_max/20)
#bounds = np.arange(-0.1, 0.1, 0.2/20)
func_plot_bounds_save(Salt_transport_0_1000_mean, Salt_lat_stream, Salt_lon_stream, bounds, '(Sv)', Plot_title, 'mill', 0, Plot_save_dir)  

Plot_title=('Mean northward transport at 2000m-3000m [1e6 kg/s] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM))
Plot_save_dir=(dir_figs+str(GCM)+'_transport_2000_3000_mean.png')
#bounds_max=np.nanpercentile(transport_0_1000_mean, 99.99)
bounds_max=float("{0:.02f}".format(np.nanpercentile(Salt_transport_2000_3000_mean, 98))) # Upper bound of plotted values to be used for colorbar, which is 98th percentile of data, with 2 decimals
bounds = np.arange(-1*bounds_max, bounds_max, bounds_max/20)
#bounds = np.arange(-0.1, 0.1, 0.2/20)
func_plot_bounds_save(Salt_transport_2000_3000_mean, Salt_lat_stream, Salt_lon_stream, bounds, '(Sv)', Plot_title, 'mill', 0, Plot_save_dir)     


levels=np.linspace(-220,220,45)  
#plot #1, just transport
fig=plt.figure()
im=plt.contourf(Lat_regrid_1D, Depths, Salt_Transport_lon_final_mean, levels, cmap=plt.cm.RdBu_r, extend = 'both')
plt.gca().invert_yaxis()
plt.xlabel('Latitude', fontsize=18)
plt.ylabel('Depth', fontsize=18)
plt.title('Mean northward salt transport (Atlantic) [1e6 kg/s] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
cbar = m.colorbar(im,"right", size="3%", pad="2%") # extend='both' will extend the colorbar in both sides (upper side and down side)
cbar.set_label('(1e6 kg/s)')
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_Salt_Transport_lon_final_mean.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
            
levels=np.linspace(-1000,1000,51)     
#plot #2, streamfunction
fig=plt.figure()
im=plt.contourf(Lat_regrid_1D, Depths, Salt_Stream_function_ave, levels, cmap=plt.cm.RdBu_r, extend = 'both')
plt.gca().invert_yaxis()
plt.xlabel('Latitude', fontsize=18)
plt.ylabel('Depth', fontsize=18)
plt.title('Salt Stream function (Atlantic) [1e6 kg/s] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
cbar = m.colorbar(im,"right", size="3%", pad="2%") # extend='both' will extend the colorbar in both sides (upper side and down side)
cbar.set_label('(1e6 kg/s)')
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_Salt_Stream_function_ave.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

#plot #3, map of mean salt transport upper 1000m with indeces of max for each lat
fig=plt.figure()
m = Basemap( projection='mill',lon_0=0)
m.fillcontinents(color='0.8')
m.drawparallels(np.arange(-90,90,20), labels=[1,1,0,1])
m.drawmeridians(np.arange(0,360,30), labels=[1,1,0,1])
bounds_max=float("{0:.02f}".format(np.nanpercentile(Salt_transport_0_1000_mean, 99.9))) # Upper bound of plotted values to be used for colorbar, which is 98th percentile of data, with 2 decimals
levels=np.linspace(-1*bounds_max,bounds_max,100)  
im=m.contourf(Salt_lon_stream,Salt_lat_stream,Salt_transport_0_1000_mean,levels,latlon=True, cmap=plt.cm.seismic)
m.scatter(lon[Salt_ii_max,Salt_jj_max],lat[Salt_ii_max,Salt_jj_max],1,c='k',latlon=True)
cbar = m.colorbar(im,"right", size="3%", pad="6%", extend = 'both') # extend='both' will extend the colorbar in both sides (upper side and down side)
cbar.set_label('(1e6 kg/s)')
plt.title('Mean salt transport upper 1000m [1e6 kg/s] (Stiplings= max transport locations) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_Salt_transport_0_1000_mean_stiple.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

#plot #4, map of mean transport 2000m-3000m with indeces of max for each lat
fig=plt.figure()
m = Basemap( projection='mill',lon_0=0)
m.fillcontinents(color='0.8')
m.drawparallels(np.arange(-90,90,20), labels=[1,1,0,1])
m.drawmeridians(np.arange(0,360,30), labels=[1,1,0,1])
bounds_max=float("{0:.02f}".format(np.nanpercentile(Salt_transport_2000_3000_mean, 99.9))) # Upper bound of plotted values to be used for colorbar, which is 98th percentile of data, with 2 decimals
levels=np.linspace(-1*bounds_max,bounds_max,100)  
im=m.contourf(Salt_lon_stream, Salt_lat_stream, Salt_transport_2000_3000_mean,levels,latlon=True, cmap=plt.cm.seismic)
m.scatter(lon[Salt_ii_min,Salt_jj_min],lat[Salt_ii_min,Salt_jj_min],1,c='k',latlon=True)
plt.show()
cbar = m.colorbar(im,"right", size="3%", pad="6%", extend = 'both') # extend='both' will extend the colorbar in both sides (upper side and down side)
cbar.set_label('(1e6 kg/s)')
plt.title('Mean salt transport 2000-3000m [1e6 kg/s] (Stiplings= min transport locations) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_Salt_transport_2000_3000_mean_stiple.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

############################################
#%% transport_lag_cor_Atlantic.py Plots ####

Salt_AMOC_transport_all = copy.deepcopy(Salt_transport_0_1000)
Salt_SMOC_transport_all= copy.deepcopy(Salt_transport_2000_3000)
Salt_SMOC_transport_all=Salt_SMOC_transport_all*(-1)
#As SMOC is transport averaged over 2000-3000m it is southward flow, thus the sign is negative (for southward flow) so we multiply by -1 to calculate further statistics

#summing transport over all longtitudes
Salt_AMOC_transport=np.nansum(Salt_AMOC_transport_all,axis=2)
Salt_SMOC_transport=np.nansum(Salt_SMOC_transport_all,axis=2)

#### WS plots ####
lag_time=40

fig=plt.figure()
data=[]
for i in range(len(Salt_AMOC_transport[0][20:161])):
    stream=runningMeanFast(Salt_AMOC_transport[:,i+lag_time], 10)
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
plt.title('WS peak convection vs. Salt transport in upper 1000m - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_Salt_lagcor_convec_transport_0_1000_WS.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


fig=plt.figure()
data=[]
for i in range(len(Salt_AMOC_max[0][20:161])):
    stream=runningMeanFast(Salt_AMOC_max[:,i+lag_time], 10)
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
plt.title('WS peak convection vs. Salt * AMOC (max [salt * streamfunction]) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_Salt_lagcor_convec_AMOC_max_WS.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


###############################################################################
###############################################################################

Salinity_allyears=np.load('Salinity_allyears_GFDL-ESM2G_500yr.npy')
for ii in range(Salinity_allyears.shape[0]): # Only keeping Atlantic cells
    for jj in range(Salinity_allyears.shape[1]):
        Salinity_allyears[ii,jj,:,:][ np.where( np.logical_and(   np.logical_and( Ocean_Index != 6 , Ocean_Index != 7, Ocean_Index != 8 ) , np.logical_and( Ocean_Index != 8 , Ocean_Index != 9 )   ) ) ] = np.nan

Salinity_allyears_150_600_atl_ave = np.nanmean( np.nanmean(Salinity_allyears[:,15:32,:,:], axis=1) , axis=2)
Salinity_allyears_150_1000_atl_ave = np.nanmean( np.nanmean(Salinity_allyears[:,15:35,:,:], axis=1) , axis=2)
Salinity_allyears_1000_2000_atl_ave = np.nanmean( np.nanmean(Salinity_allyears[:,34:40,:,:], axis=1) , axis=2) 
Salinity_allyears_1000_4200_atl_ave = np.nanmean( np.nanmean(Salinity_allyears[:,34:47,:,:], axis=1) , axis=2) 


levels=np.linspace(-0.1,0.1,21)     
#plot #2, streamfunction
fig=plt.figure()
im=plt.contourf(np.arange(1,501,1) , Lat_regrid_1D, np.transpose(Salinity_allyears_150_600_atl_ave - np.matlib.repmat(np.nanmean(Salinity_allyears_150_600_atl_ave, axis=0), 500,1)), levels, cmap=plt.cm.RdBu_r, extend = 'both')
plt.gca().invert_yaxis()
plt.xlabel('Years', fontsize=18)
plt.ylabel('Latitude', fontsize=18)
plt.title('Salt anomalies (Atlantic) 150-600m [PSU] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
cbar = m.colorbar(im,"right", size="3%", pad="2%") # extend='both' will extend the colorbar in both sides (upper side and down side)
cbar.set_label('(PSU)')
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_Salinity_anomaly_150_600.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

levels=np.linspace(-0.1,0.1,21)     
#plot #2, streamfunction
fig=plt.figure()
im=plt.contourf(np.arange(1,501,1) , Lat_regrid_1D, np.transpose(Salinity_allyears_150_1000_atl_ave - np.matlib.repmat(np.nanmean(Salinity_allyears_150_1000_atl_ave, axis=0), 500,1)), levels, cmap=plt.cm.RdBu_r, extend = 'both')
plt.gca().invert_yaxis()
plt.xlabel('Years', fontsize=18)
plt.ylabel('Latitude', fontsize=18)
plt.title('Salt anomalies (Atlantic) 150-1000m [PSU] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
cbar = m.colorbar(im,"right", size="3%", pad="2%") # extend='both' will extend the colorbar in both sides (upper side and down side)
cbar.set_label('(PSU)')
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_Salinity_anomaly_150_1000.png', format='png', dpi=300, transparent=True, bbox_inches='tight')



Salinity_allyears_atl_ave = np.nanmean( np.nanmean(Salinity_allyears[0:10,:,:,:], axis=0) , axis=2)


levels=np.linspace(33.6,36,21)     
#plot #2, streamfunction
fig=plt.figure()
im=plt.contourf(Lat_regrid_1D, Depths, Salinity_allyears_atl_ave, levels, cmap=plt.cm.RdBu_r, extend = 'both')
plt.gca().invert_yaxis()
plt.xlabel('Latitude', fontsize=18)
plt.ylabel('Depths (m)', fontsize=18)
plt.title('Salinity (Atlantic) [PSU] - 10 years', fontsize=18)
cbar = m.colorbar(im,"right", size="3%", pad="2%") # extend='both' will extend the colorbar in both sides (upper side and down side)
cbar.set_label('(PSU)')
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_Salinity_atl_10years.png', format='png', dpi=300, transparent=True, bbox_inches='tight')



Salinity_allyears_300_1000_atl_ave = np.nanmean( np.nanmean(Salinity_allyears[:,26:35,:,:], axis=1) , axis=2)

levels=np.linspace(-0.1,0.1,21)     
#plot #2, streamfunction
fig=plt.figure()
im=plt.contourf(np.arange(1,501,1) , Lat_regrid_1D, np.transpose(Salinity_allyears_300_1000_atl_ave - np.matlib.repmat(np.nanmean(Salinity_allyears_300_1000_atl_ave, axis=0), 500,1)), levels, cmap=plt.cm.RdBu_r, extend = 'both')
plt.gca().invert_yaxis()
plt.xlabel('Years', fontsize=18)
plt.ylabel('Latitude', fontsize=18)
plt.title('Salt anomalies (Atlantic) 300-1000m [PSU] - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
cbar = m.colorbar(im,"right", size="3%", pad="2%") # extend='both' will extend the colorbar in both sides (upper side and down side)
cbar.set_label('(PSU)')
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_figs+str(GCM)+'_Salinity_anomaly_300_1000.png', format='png', dpi=300, transparent=True, bbox_inches='tight')






###############################################################################
###############################################################################
mask_atl=calc_Atl_Mask()
t_frequency='Omon'
variable_so='so' # Sea Water Salinity - units: psu (gr/kg)

dir_data_in1 = ('/data5/scratch/CMIP5/CMIP5_models/ocean_physics/') # Directory to raed raw data from
dir_data_in2=(dir_data_in1+ GCM + '/piControl/mo/')

dset_so = netcdf_read (dir_data_in2+str(variable_so)+'_'+str(t_frequency)+'_'+str(GCM)+'*12.nc',str(GCM), variable_so)
Depths=dset_so.lvl[:]

Salt_Transport_lon_final, Salt_transport_0_1000, Salt_transport_1000_2000, Salt_transport_2000_3000, Salt_transport_3000_4000, Salt_transport_4000_5000, Salt_transport_5000_below, Salt_latlon_depths, Salt_Depth_indx, Salt_lon_stream, Salt_lat_stream, Salt_ii_max, Salt_jj_max, Salt_ii_min, Salt_jj_min, Salt_transport_0_1000_mean, Salt_transport_2000_3000_mean 
 = func_salt_transport(dset_vo, dset_so, year_start, year_end, mask_atl, 180, 360)

start,end = dset_so.find_time(dset_so.times, year_start, year_end)

for t in range(int((end+1-start)/12)):
    print ('Salinity - Year: ', t+1)

    data_extracted_so=dset_so.extract_data(dset_so.variable,start+12*t,start+12*t+11)
    data_so=np.squeeze(data_extracted_so)
    data_so=np.mean(data_so, axis=0) # unit: gr/kg = gr/l = 1000 gr/m3 = kg/m3    
    data_so[data_so > 1e18] = np.nan
    
    data_i = func_regrid(data_so, dset_so.y, dset_so.x, Lat_regrid_2D, Lon_regrid_2D)
    
    data_i[ np.where( np.logical_and(   np.logical_and( Ocean_Index != 6 , Ocean_Index != 7, Ocean_Index != 8 ) , np.logical_and( Ocean_Index != 8 , Ocean_Index != 9 )   ) ) ] = np.nan
    
    
    #### calculating volume transport
    #### first multiplying by 111km*cos(lat)
    mul_by_lat=data_i*(np.cos(np.deg2rad(lat))*111321) # S * V_y * dX [Unit: kg/m2.s * m = kg/m.s]
    #### second multiplying by depth
    transport=mul_by_lat*data_depth_ranges/1000000 # S * V_y * dX * dZ [Unit: kg/m.s *m * 1e-6  = kg/s * 10e6] 
    #### calculating integral over dz
    transport_lon=np.nansum(transport,axis=2) # SUM(S * V_y * dX * dZ) over longitudes [Unit: kg/m.s *m * 1e-6  = kg/s * 10e6]




















depths=dset_so.extract_depth() # Units: 'm'
### these are upper and lower depths of each cell in an ocean grid
depths_b=dset_so.extract_depth_bounds()
### calculate the depth of each cell in an ocean grid
depths_r=depths_b[:,1]-depths_b[:,0]
### find timeindeces
start,end = dset_so.find_time(dset_so.times, year_start, year_end)
transport_lon_final=[] # SUM(V_y * dX * dZ) over longitudes, for all years
transport_final=[]
transport_0_1000=[]
transport_1000_2000=[]
transport_2000_3000=[]
transport_3000_4000=[]
transport_4000_5000=[]
transport_5000_below=[]

Depth_indx=np.zeros((5,4)) # Rows: Depths of 1000, 2000, 3000, 4000 and 5000 meters
                           # Columns: 0=row number, 1=depths, 2=depths lower range, 3=depths upper range    
for dd in range(depths.shape[0]):
    if (depths_b[dd,0] < 1000) and (depths_b[dd,1] > 1000):
        Depth_indx[0,0]=dd; Depth_indx[0,1]=depths[dd]
        Depth_indx[0,2]=depths_b[dd,0]; Depth_indx[0,3]=depths_b[dd,1]
    if (depths_b[dd,0] < 2000) and (depths_b[dd,1] > 2000):
        Depth_indx[1,0]=dd; Depth_indx[1,1]=depths[dd]
        Depth_indx[1,2]=depths_b[dd,0]; Depth_indx[1,3]=depths_b[dd,1]
    if (depths_b[dd,0] < 3000) and (depths_b[dd,1] > 3000):
        Depth_indx[2,0]=dd; Depth_indx[2,1]=depths[dd]
        Depth_indx[2,2]=depths_b[dd,0]; Depth_indx[2,3]=depths_b[dd,1]
    if (depths_b[dd,0] < 4000) and (depths_b[dd,1] > 4000):
        Depth_indx[3,0]=dd; Depth_indx[3,1]=depths[dd]
        Depth_indx[3,2]=depths_b[dd,0]; Depth_indx[3,3]=depths_b[dd,1]
    if (depths_b[dd,0] < 5000) and (depths_b[dd,1] > 5000):
        Depth_indx[4,0]=dd; Depth_indx[4,1]=depths[dd]
        Depth_indx[4,2]=depths_b[dd,0]; Depth_indx[4,3]=depths_b[dd,1]

Lat_regrid_1D, Lon_regrid_1D, Lat_bound_regrid, Lon_bound_regrid = func_latlon_regrid_eq(180, 360, -90, 90, 0, 360)
lon, lat = np.meshgrid(Lon_regrid_1D, Lat_regrid_1D)

### we averaging velocities over the year for monthly data
for t in range(int((end+1-start)/12)):
    print ('Stream calc - Year: ', t)

    data_extracted_so=dset_so.extract_data(dset_so.variable,start+12*t,start+12*t+11)
    data_so=np.squeeze(data_extracted_so)
    data_so=np.mean(data_so, axis=0) # unit: gr/kg = gr/l = 1000 gr/m3 = kg/m3    
    data_so[data_so > 1e18] = np.nan
    
    data = np.multiply(data_vo, data_so) # [Unit: kg/m3 * m/s = kg/m2.s ]
    
#        lon,lat,data_i=interpolate_2_reg_grid(args[0].x,args[0].y,data)
    data_i = func_regrid(data, dset_so.y, dset_so.x, lat, lon)
    
    ##### converting 1e+20 to nan ######
    data_i[data_i>1000]=np.nan
    if t==0: 
        data_depth=np.full([len(lon),len(lon[0])], np.nan)
        data_depth_ranges=np.full([len(data_i),len(data_i[0]),len(data_i[0][0])], np.nan)
    [ii,jj] = mask_atl
    for k in range(len(ii)):
        #### I calculate the depth by looking how many nans is in the depth column
        if sum(~np.isnan(data_i[:,ii[k],jj[k]]))>0:
            if t==0:
                data_depth[ii[k],jj[k]]=depths[sum(~np.isnan(data_i[:,ii[k],jj[k]]))-1]
                for l in range(sum(~np.isnan(data_i[:,ii[k],jj[k]]))):
                    data_depth_ranges[l,ii[k],jj[k]]=depths_r[l]
                    
    #### calculating volume transport
    #### first multiplying by 111km*cos(lat)
    mul_by_lat=data_i*(np.cos(np.deg2rad(lat))*111321) # S * V_y * dX [Unit: kg/m2.s * m = kg/m.s]
    #### second multiplying by depth
    transport=mul_by_lat*data_depth_ranges/1000000 # S * V_y * dX * dZ [Unit: kg/m.s *m * 1e-6  = kg/s * 10e6] 
    #### calculating integral over dz
    transport_lon=np.nansum(transport,axis=2) # SUM(S * V_y * dX * dZ) over longitudes [Unit: kg/m.s *m * 1e-6  = kg/s * 10e6]
    #### calculating cum integral over dz
    #stream=np.nancumsum(transport_lon,axis=0)
    transport_lon_final.append(transport_lon) # SUM(S * V_y * dX * dZ) over longitudes, for all years
    transport_0_1000.append(np.nanmean(transport[0:int(Depth_indx[0,0])+1,:,:],axis=0))
    transport_1000_2000.append(np.nanmean(transport[int(Depth_indx[0,0])+1:int(Depth_indx[1,0])+1,:,:],axis=0))        
    transport_2000_3000.append(np.nanmean(transport[int(Depth_indx[1,0])+1:int(Depth_indx[2,0])+1,:,:],axis=0))
    transport_3000_4000.append(np.nanmean(transport[int(Depth_indx[2,0])+1:int(Depth_indx[3,0])+1,:,:],axis=0))
    transport_4000_5000.append(np.nanmean(transport[int(Depth_indx[3,0])+1:int(Depth_indx[4,0])+1,:,:],axis=0))
    if depths.shape[0] - (int(Depth_indx[4,0])+1) ==1:
        transport_5000_below.append(transport[int(Depth_indx[4,0])+1,:,:])
    elif depths.shape[0] - (int(Depth_indx[4,0])+1) >=2:
        transport_5000_below.append(np.nanmean(transport[int(Depth_indx[4,0])+1:,:,:],axis=0))

    #transport_final.append(transport)
transport_lon_final=np.asarray(transport_lon_final)

transport_0_1000=np.asarray(transport_0_1000)
transport_1000_2000=np.asarray(transport_1000_2000)      
transport_2000_3000=np.asarray(transport_2000_3000)
transport_3000_4000=np.asarray(transport_3000_4000)
transport_4000_5000=np.asarray(transport_4000_5000)
transport_5000_below=np.asarray(transport_5000_below)

transport_0_1000_mean=np.nanmean(transport_0_1000,axis=0)
transport_2000_3000_mean=np.nanmean(transport_2000_3000,axis=0)

maxvals=np.nanmax(transport_0_1000_mean,axis=1)
minvals=np.nanmin(transport_2000_3000_mean,axis=1)

ind_max = np.array([np.argwhere(transport_0_1000_mean == [x]) for x in maxvals])
ind_max=np.concatenate(ind_max).astype(None)
ind_max=np.concatenate(ind_max).astype(None)
ii_max=ind_max[0::2].astype(int)
jj_max=ind_max[1::2].astype(int)
ind_min = np.array([np.argwhere(transport_2000_3000_mean == [x]) for x in minvals])
ind_min=np.concatenate(ind_min).astype(None)
ind_min=np.concatenate(ind_min).astype(None)
ii_min=ind_min[0::2].astype(int)
jj_min=ind_min[1::2].astype(int)



###############################################################################
###############################################################################
###############################################################################
###############################################################################



#################################################
#%% Convection Index Calculations - WS  Salt  ###

import os
import shelve

dir_pwd = os.getcwd() # Gets the current directory (and in which the code is placed)
filename_out = (dir_pwd + '/Results/AllResults_'+GCM+'_500yr.out') # Directory to save processed data

my_shelf = shelve.open(filename_out)
for key in my_shelf:
    globals()[key]=my_shelf[key]
my_shelf.close()

######################################
GCM = 'GFDL-ESM2G'
year_start=1
year_end=500

conv_index_depth_ws=500 # args[7] in func_time_depth_plot code - depth for Convection Index

t_frequency='Omon'
variable_so='so' # Sea Water Salinity

dir_data_in1 = ('/data5/scratch/CMIP5/CMIP5_models/ocean_physics/') # Directory to raed raw data from
dir_data_in2=(dir_data_in1+ GCM + '/piControl/mo/')

dset_so = netcdf_read (dir_data_in2+str(variable_so)+'_'+str(t_frequency)+'_'+str(GCM)+'*12.nc',str(GCM), variable_so)
Depths=dset_so.lvl[:]


WS_salt_time_depth, WS_salt_index=func_time_depth_plot(dset_so, year_start, year_end, WS_indeces_lonlat, 90, 180, conv_index_depth_ws) # WS_index is the convection indeces at 500m depth



#############################################################################
###########                 Saving Results                ###################
import os
import shelve

dir_pwd = os.getcwd() # Gets the current directory (and in which the code is placed)
filename_out = (dir_pwd + '/AllResults_'+GCM+'-WS_salt_index_500yr.out') # Directory to save processed data

### To save
my_shelf = shelve.open(filename_out,'n') # 'n' for new
 
var_list_short=['WS_salt_time_depth', 'WS_salt_index', 'Depths']

for key in var_list_short:
    try:
        my_shelf[key] = globals()[key]
    except TypeError:
        print('ERROR shelving: {0}'.format(key))
my_shelf.close()


import gsw as g

WS_T = WS_index_time_depth
WS_S = WS_salt_time_depth

WS_p = g.p_from_z(-Depths,-50) #pressure at all depths assuming Ross sea is on avg ~50S

WS_SA = g.SA_from_SP(WS_S,WS_p,195,-50)
WS_CT = g.CT_from_pt(WS_SA,WS_T)
WS_alpha= g.alpha(WS_SA,WS_CT,WS_p)
WS_beta= g.beta(WS_SA,WS_CT,WS_p)
WS_rho = g.rho(WS_SA,WS_CT,WS_p)
WS_sigma= g.sigma0(WS_SA,WS_CT)


plt.suptitle('Density in Convection Area (Weddell Sea) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
plt.subplot(211)
plt.title('depth = 2.5m')
plt.plot(WS_rho[:,0]-1000)
plt.subplot(212)
plt.title('depth = 500m')
plt.plot(ROSS_rho[:,13]-1000)
plt.ylabel('(rho - 1000) [kg/m^3]',fontsize=18)
plt.xlabel('Years', fontsize=18)


#########################
s_term_500 = WS_beta[:,13]*(WS_SA[:,13]-np.mean(WS_SA[:,13]))
t_term_500 = WS_alpha[:,13]*(WS_CT[:,13]-np.mean(WS_CT[:,13]))

s_term_0 = WS_beta[:,0]*(WS_SA[:,0]-np.mean(WS_SA[:,0]))
t_term_0 = WS_alpha[:,0]*(WS_CT[:,0]-np.mean(WS_CT[:,0]))

plt.plot(WS_alpha[:,13]*(WS_T[:,13]-np.mean(WS_T[:,13])),'r')
plt.yticks(fontsize = 15)
plt.ylabel('alpha * (CT - CT_0), beta * (SA - SA_0)',fontsize = 18)
plt.plot(WS_beta[:,13]*(WS_S[:,13]-np.mean(WS_S[:,13])))
plt.xlabel('Years', fontsize=18)
plt.xticks(fontsize = 15); plt.yticks(fontsize = 15)
plt.title('Density separation in Convection Area (Weddell Sea), 500m - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
plt.legend(['alpha * (CT - CT_0)','beta * (SA - SA_0)'],loc ='lower left')



plt.plot((1 + s_term_500-t_term_500 + WS_p[13]/10.1325/1000/1500**2)*np.mean(WS_rho[:,13]))
plt.plot(WS_rho[:,13])
plt.ylabel('rho [kg/m^3]',fontsize=18)
plt.suptitle('Density in Convection Area (Weddell Sea), 500m - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
plt.xlabel('Years',fontsize = 18)
plt.legend(['rho','rho computed with alpha, beta'])

#######################
WS_rho_constS= g.rho(np.mean(WS_SA,0),WS_CT,WS_p)
WS_rho_constT= g.rho(WS_SA,np.mean(WS_CT,0),WS_p)


plt.plot(WS_rho_constS[:,13]-1000,'r')
plt.plot(WS_rho_constT[:,13]-1000,'b')
plt.plot(WS_rho[:,13]-1000,'g')
plt.xlabel('Years', fontsize=18)
plt.ylabel('(rho - 1000) [kg/m^3]', fontsize=18)
plt.xticks(fontsize = 15); plt.yticks(fontsize = 15)
plt.legend(['Rho(mean SA, variable CT)','Rho(variable SA, mean CT)','Rho(SA,CT)'])
plt.title('rho calculated with constant SA, CT in Convection Area (Weddell Sea), 500m - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)

plt.plot(WS_rho_constS[:,0]-1000,'r')
plt.plot(WS_rho_constT[:,0]-1000,'b')
plt.plot(WS_rho[:,0]-1000,'g')
plt.xlabel('Years', fontsize=18)
plt.ylabel('(Rho - 1000) [kg/m^3]', fontsize=18)
plt.xticks(fontsize = 15); plt.yticks(fontsize = 15)
plt.legend(['Rho(mean SA, variable CT)','Rho(variable SA, mean CT)','Rho(SA,CT)'])
plt.title('rho calculated with constant SA, CT in Convection Area (Weddell Sea), ' + str(Depths[0]) + 'm - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)









