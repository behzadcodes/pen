#################################
####   CMIP5 - CO2 & Phyto   ####
########################################
####     Behzad Asadieh, Ph.D.      ####
####  University of Pennsylvania    ####
####    basadieh@sas.upenn.edu      ####
########################################

### IMPORTANT: The following path needs to be appended so that the neccessary codes can be read from that directory ###
import sys
sys.path.append('/data1/home/basadieh/behzadcodes/behzadlibrary')
from Behzadlib import func_latlon_regrid, func_latlon_regrid_eq, func_oceanlandmask, func_oceanindex, func_regrid, func_plot, func_plot_save, func_plot_bounds_save
from BehzadlibPlot import func_plot2Dcolor_contourf, func_plotline_1var, func_plotline_2var, func_plot2Dcontour_1var, func_plot2Dcontour_2var
from CMIP5lib import netcdf_read
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
from mpl_toolkits import mplot3d
from Behzadlib import func_corr_3d
from Behzadlib import func_detrend_3d
########################################

dir_pwd = os.getcwd() # Gets the current directory (and in which the code is placed)
dir_data_in = ('/data5/scratch/Behzad/ocean/ocean_biogeochemistry/processed_data/python_data_CMIP5_bio_1991_2010/')
dir_figs = (dir_pwd + '/Figures_MultiGCM/') # Directory to save processed data
dir_results = (dir_pwd + '/Results/') # Directory to save results

### Regrdridding calculations ###
# creating new coordinate grid, same which was used in interpolation in data processing code
lat_n_regrid, lon_n_regrid = 180, 360 # Number of Lat and Lon elements in the regridded data
lon_min_regrid, lon_max_regrid = 0, 360 # Min and Max value of Lon in the regridded data
lat_min_regrid, lat_max_regrid = -90, 90 # Min and Max value of Lat in the regridded data

Lat_regrid_1D, Lon_regrid_1D, Lat_bound_regrid, Lon_bound_regrid = func_latlon_regrid(lat_n_regrid, lon_n_regrid, lat_min_regrid, lat_max_regrid, lon_min_regrid, lon_max_regrid)
Lon_regrid_2D, Lat_regrid_2D = np.meshgrid(Lon_regrid_1D, Lat_regrid_1D)

Ocean_Land_mask = func_oceanlandmask(Lat_regrid_2D, Lon_regrid_2D) # 1= ocean, 0= land

Ocean_Index = func_oceanindex (Lat_regrid_2D, Lon_regrid_2D) # [0=land] [2=Pacific] [3=Indian Ocean] [6=Atlantic] [10=Arctic] [8=Baffin Bay (west of Greenland)] [9=Norwegian Sea (east of Greenland)] [11=Hudson Bay (Canada)] 
###############################################################################

### Variables to be edited by user ###
start_date_plt=1991 # Start Date used for plotting
end_date_plt=2010 # End Date used for plotting

yrs_n=20;

GCM_Names = ['CESM1-BGC' ,'GFDL-ESM2G','GFDL-ESM2M','IPSL-CM5A-MR','IPSL-CM5A-LR','MRI-ESM1','GISS-E2-H-CC','GISS-E2-R-CC','HadGEM2-ES','CanESM2','NorESM1-ME','MPI-ESM-MR','MPI-ESM-LR']
Colors    =  [[1, 1, 0],   [1, 0.8, 0],  [1, 0.6, 0], [1, 0.4, 0],   [1, 0, 0], [1, 0.8, 1],   [1, 0.5, 1],   [1, 0, 1],   [0, 0, 0.5], [0, 0, 1], [0, 0.5, 1],[0, 0.5, 0],  [0, 1, 0]]

############################################
##   Multi-GCM SingleVariable Plotting  ####
############################################
Time_months = ('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec');
P_Var_x = np.linspace(1,12,12);

d_trend='yes'


fig=plt.figure() 
P_title= 'alpha value [ Ln(PP) / Ln(P) = alpha+1 ]'
if d_trend=='yes':
    plt.title('Global - '+str(start_date_plt)+'-'+str(end_date_plt)+ ' - detrended', fontsize=24)
else:
    plt.title('Global - '+str(start_date_plt)+'-'+str(end_date_plt), fontsize=24)

for M_i in range(len(GCM_Names)): # M_i=1     # GCM=GCM_Names[0]
    
    GCM=GCM_Names[M_i]
    print (GCM)     
    
    filename_in = (dir_data_in + 'AllResults_'+GCM+'_bio_1991_2010.out') # Directory to save processed data
    my_shelf = shelve.open(filename_in)
    if 'Phyc_allmonths' in my_shelf and 'PPint_allmonths' in my_shelf: # If that variable exists in the saved data    
        globals()['Phyc_allmonths']=my_shelf['Phyc_allmonths']
        globals()['PPint_allmonths']=my_shelf['PPint_allmonths']
        if d_trend=='yes': # If it is requested to detrend the data before calculations
            Phyc_allmonths= func_detrend_3d(Phyc_allmonths)
            PPint_allmonths= func_detrend_3d(PPint_allmonths)
                
        Phyc_monthly = Phyc_allmonths.reshape(( np.int(Phyc_allmonths.shape[0]/12) ,12,lat_n_regrid,lon_n_regrid)) # Monthly means of SST over the time period
        Phyc_monthly = np.nanmean(Phyc_monthly,axis=0)
        PPint_monthly = PPint_allmonths.reshape(( np.int(PPint_allmonths.shape[0]/12) ,12,lat_n_regrid,lon_n_regrid)) # Monthly means of SST over the time period
        PPint_monthly = np.nanmean(PPint_monthly,axis=0)  
        
        alpha_monthly=np.empty(12)*nan
        for ii in range (0,12):
            #alpha_1= np.log(PPint_monthly[ii,:,:]) / np.log(Phyc_monthly[ii,:,:])
            #alpha_monthly[ii] = np.nanmean(np.nanmean(alpha_1,axis=0),axis=0)
            alpha_monthly[ii] = np.log( np.nanmean(np.nanmean(PPint_monthly[ii,:,:],axis=0),axis=0) ) / np.log( np.nanmean(np.nanmean(Phyc_monthly[ii,:,:],axis=0),axis=0) )
        alpha_monthly = alpha_monthly - 1
        
        
        plt.plot(P_Var_x, alpha_monthly, c=Colors[M_i], label=GCM, marker='o', markersize=7, markerfacecolor='W', linewidth=2.0)
        plt.xticks(fontsize = 26); plt.yticks(fontsize = 26)
        #plt.xlabel('Month', fontsize=26)
        plt.ylabel('alpha', fontsize=26)
        #plt.ylim(275,450)
    
plt.xticks(P_Var_x, '', rotation=45) # Setting the X-ticks ## rotation='horizontal' # rotation='vertical'
plt.legend(prop={'size': 14}, loc='center right', bbox_to_anchor=(1.4, 0.5), fancybox=True, framealpha=0.8)
plt.xticks(P_Var_x, Time_months, rotation=45, fontsize=24) # Setting the X-ticks ## rotation='horizontal' # rotation='vertical'
plt.subplots_adjust(left=0.33, bottom=0.11, right=0.69, top=0.9, hspace=0.2, wspace=0.05) # the amount of height/width reserved for space between subplots
plt.suptitle(P_title, fontsize=26)    
plt.show()
plt.grid()
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full    

if d_trend=='yes':
    fig.savefig(dir_figs+'AllModels_Alpha_value_detrended_Global_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
else:
    fig.savefig(dir_figs+'AllModels_Alpha_value_Global_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')




d_trend='yes'

if d_trend=='yes':
    P_title= 'alpha value [ Ln(PP) / Ln(P) = alpha+1 ] - N.P. and N.A. - detrended - '+str(start_date_plt)+'-'+str(end_date_plt)
else:
    P_title= 'alpha value [ Ln(PP) / Ln(P) = alpha+1 ] - N.P. and N.A. - '+str(start_date_plt)+'-'+str(end_date_plt)

fig=plt.figure() 
for iii in range(0,2):
    ax = fig.add_subplot(1,2,iii+1)

    for M_i in range(len(GCM_Names)): # M_i=1     # GCM=GCM_Names[0]
        
        GCM=GCM_Names[M_i]
        print (GCM)     
        
        filename_in = (dir_data_in + 'AllResults_'+GCM+'_bio_1991_2010.out') # Directory to save processed data
        my_shelf = shelve.open(filename_in)
        if 'Phyc_allmonths' in my_shelf and 'PPint_allmonths' in my_shelf: # If that variable exists in the saved data    
            globals()['Phyc_allmonths']=my_shelf['Phyc_allmonths']
            globals()['PPint_allmonths']=my_shelf['PPint_allmonths']
                
            if iii==0: # 'PAPA'
                Phyc_allmonths_BOX = Phyc_allmonths[:,135:140,210:220]
                PPint_allmonths_BOX = PPint_allmonths[:,135:140,210:220]
            elif iii==1: # 'NABE'
                Phyc_allmonths_BOX = Phyc_allmonths[:,135:140,325:335]
                PPint_allmonths_BOX = PPint_allmonths[:,135:140,325:335]

            if d_trend=='yes': # If it is requested to detrend the data before calculations
                Phyc_allmonths_BOX= func_detrend_3d(Phyc_allmonths_BOX)
                PPint_allmonths_BOX= func_detrend_3d(PPint_allmonths_BOX)
                    
            Phyc_monthly = Phyc_allmonths_BOX.reshape(( np.int(Phyc_allmonths_BOX.shape[0]/12) ,12,Phyc_allmonths_BOX.shape[1],Phyc_allmonths_BOX.shape[2])) # Monthly means of SST over the time period
            Phyc_monthly = np.nanmean(Phyc_monthly,axis=0)
            PPint_monthly = PPint_allmonths_BOX.reshape(( np.int(PPint_allmonths_BOX.shape[0]/12) ,12,PPint_allmonths_BOX.shape[1],PPint_allmonths_BOX.shape[2])) # Monthly means of SST over the time period
            PPint_monthly = np.nanmean(PPint_monthly,axis=0)  

            alpha_monthly=np.empty(12)*nan
            for ii in range (0,12):
                alpha_monthly[ii] = np.log( np.nanmean(np.nanmean(PPint_monthly[ii,:,:],axis=0),axis=0) ) / np.log( np.nanmean(np.nanmean(Phyc_monthly[ii,:,:],axis=0),axis=0) )
            alpha_monthly = alpha_monthly - 1
            
            plt.plot(P_Var_x, alpha_monthly, c=Colors[M_i], label=GCM, marker='o', markersize=7, markerfacecolor='W', linewidth=2.0)
            plt.xticks(fontsize = 26); plt.yticks(fontsize = 26)
            plt.xticks(P_Var_x, '', rotation=45) # Setting the X-ticks ## rotation='horizontal' # rotation='vertical'
            plt.xticks(P_Var_x, Time_months, rotation=45, fontsize=24) # Setting the X-ticks ## rotation='horizontal' # rotation='vertical'
            plt.xlabel('Month', fontsize=26)
            plt.ylabel('alpha', fontsize=26)
            if iii==0: # 'PAPA'
                plt.title('North Pacific [45N-50N, 140W-150W]', fontsize=24)
            elif iii==1: # 'NABE'
                plt.title('North Atlantic [45N-50N, 25W-35W]', fontsize=24)
            plt.ylim(0.7,2.5)
            plt.grid()

plt.legend(prop={'size': 14}, loc='center right', bbox_to_anchor=(1.45, 0.5), fancybox=True, framealpha=0.8)
plt.subplots_adjust(left=0.1, bottom=0.15, right=0.8, top=0.85, hspace=0.2, wspace=0.25) # the amount of height/width reserved for space between subplots
plt.suptitle(P_title, fontsize=26)   
plt.show()
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full    

if d_trend=='yes':
    fig.savefig(dir_figs+'AllModels_Alpha_value_detrended_NPNA_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
else:
    fig.savefig(dir_figs+'AllModels_Alpha_value_NPNA_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')












from scipy import stats
d_trend='no'

fig=plt.figure() 
P_title= 'alpha value [ Ln(PP) = (alpha+1) * Ln(P) + beta ]'
if d_trend=='yes':
    plt.title('Global - '+str(start_date_plt)+'-'+str(end_date_plt)+ ' - detrended', fontsize=24)
else:
    plt.title('Global - '+str(start_date_plt)+'-'+str(end_date_plt), fontsize=24)

for M_i in range(len(GCM_Names)): # M_i=1     # GCM=GCM_Names[0]
    
    GCM=GCM_Names[M_i]
    print (GCM)     
    
    filename_in = (dir_data_in + 'AllResults_'+GCM+'_bio_1991_2010.out') # Directory to save processed data
    my_shelf = shelve.open(filename_in)
    if 'Phyc_allmonths' in my_shelf and 'PPint_allmonths' in my_shelf: # If that variable exists in the saved data    
        globals()['Phyc_allmonths']=my_shelf['Phyc_allmonths']
        globals()['PPint_allmonths']=my_shelf['PPint_allmonths']
        if d_trend=='yes': # If it is requested to detrend the data before calculations
            Phyc_allmonths= func_detrend_3d(Phyc_allmonths)
            PPint_allmonths= func_detrend_3d(PPint_allmonths)
                
        Phyc_monthly = Phyc_allmonths.reshape(( np.int(Phyc_allmonths.shape[0]/12) ,12,lat_n_regrid,lon_n_regrid)) # Monthly means of SST over the time period
        Phyc_monthly = np.nanmean(Phyc_monthly,axis=0)
        PPint_monthly = PPint_allmonths.reshape(( np.int(PPint_allmonths.shape[0]/12) ,12,lat_n_regrid,lon_n_regrid)) # Monthly means of SST over the time period
        PPint_monthly = np.nanmean(PPint_monthly,axis=0)  


        alpha_monthly=np.empty(12)*nan
        beta_monthly=np.empty(12)*nan
        pvalue_monthly=np.empty(12)*nan
        
        for tt in range (0,12):
            
            Var_y=np.log(PPint_monthly[tt,:,:])
            Var_x=np.log(Phyc_monthly[tt,:,:])
            Var_y[ np.abs(Var_y)>1e19 ]=nan;  Var_x[ np.abs(Var_x)>1e19 ]=nan
            xx = Var_x.reshape( np.int(lat_n_regrid*lon_n_regrid)) # Creating a time series of the variabel for that month 
            yy = Var_y.reshape( np.int(lat_n_regrid*lon_n_regrid))             
            yy=yy[np.logical_not(np.isnan(xx))] ; xx=xx[np.logical_not(np.isnan(xx))] # Excluding any NaN values in the time series                   
            xx=xx[np.logical_not(np.isnan(yy))] ; yy=yy[np.logical_not(np.isnan(yy))] # Excluding any NaN values in the time series 
            
            slope_ij, intercept_ij, r_value_ij, p_value_ij, std_err_ij = stats.linregress(xx, yy)

            alpha_monthly[tt] = slope_ij
            beta_monthly[tt] = intercept_ij
            pvalue_monthly[tt] = p_value_ij

        alpha_monthly = alpha_monthly - 1
        
        plt.plot(P_Var_x, alpha_monthly, c=Colors[M_i], label=GCM, marker='o', markersize=7, markerfacecolor='W', linewidth=2.0)
        plt.xticks(fontsize = 26); plt.yticks(fontsize = 26)
        #plt.xlabel('Month', fontsize=26)
        plt.ylabel('alpha', fontsize=26)
        #plt.ylim(275,450)
    
plt.xticks(P_Var_x, '', rotation=45) # Setting the X-ticks ## rotation='horizontal' # rotation='vertical'
plt.legend(prop={'size': 14}, loc='center right', bbox_to_anchor=(1.4, 0.5), fancybox=True, framealpha=0.8)
plt.xticks(P_Var_x, Time_months, rotation=45, fontsize=24) # Setting the X-ticks ## rotation='horizontal' # rotation='vertical'
plt.subplots_adjust(left=0.33, bottom=0.11, right=0.69, top=0.9, hspace=0.2, wspace=0.05) # the amount of height/width reserved for space between subplots
plt.suptitle(P_title, fontsize=26)    
plt.show()
plt.grid()
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full    

if d_trend=='yes':
    fig.savefig(dir_figs+'AllModels_Alpha_value_regress_detrended_Global_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
else:
    fig.savefig(dir_figs+'AllModels_Alpha_value_regress_Global_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')



from scipy import stats
d_trend='no'

fig=plt.figure() 
P_title= 'alpha value [ Ln(PP) = (alpha+1) * Ln(P) + beta ] - N.P. and N.A. - '+str(start_date_plt)+'-'+str(end_date_plt)

for iii in range(0,2):
    ax = fig.add_subplot(1,2,iii+1)

    for M_i in range(len(GCM_Names)): # M_i=1     # GCM=GCM_Names[0]
        
        GCM=GCM_Names[M_i]
        print (GCM)     
        
        filename_in = (dir_data_in + 'AllResults_'+GCM+'_bio_1991_2010.out') # Directory to save processed data
        my_shelf = shelve.open(filename_in)
        if 'Phyc_allmonths' in my_shelf and 'PPint_allmonths' in my_shelf: # If that variable exists in the saved data    
            globals()['Phyc_allmonths']=my_shelf['Phyc_allmonths']
            globals()['PPint_allmonths']=my_shelf['PPint_allmonths']
            if d_trend=='yes': # If it is requested to detrend the data before calculations
                Phyc_allmonths= func_detrend_3d(Phyc_allmonths)
                PPint_allmonths= func_detrend_3d(PPint_allmonths)
                
            if iii==0: # 'PAPA'
                Phyc_allmonths_BOX = Phyc_allmonths[:,135:140,210:220]
                PPint_allmonths_BOX = PPint_allmonths[:,135:140,210:220]
            elif iii==1: # 'NABE'
                Phyc_allmonths_BOX = Phyc_allmonths[:,135:140,325:335]
                PPint_allmonths_BOX = PPint_allmonths[:,135:140,325:335]
                    
            Phyc_monthly = Phyc_allmonths_BOX.reshape(( np.int(Phyc_allmonths_BOX.shape[0]/12) ,12,Phyc_allmonths_BOX.shape[1],Phyc_allmonths_BOX.shape[2])) # Monthly means of SST over the time period
            Phyc_monthly = np.nanmean(Phyc_monthly,axis=0)
            PPint_monthly = PPint_allmonths_BOX.reshape(( np.int(PPint_allmonths_BOX.shape[0]/12) ,12,PPint_allmonths_BOX.shape[1],PPint_allmonths_BOX.shape[2])) # Monthly means of SST over the time period
            PPint_monthly = np.nanmean(PPint_monthly,axis=0)  
    
            alpha_monthly=np.empty(12)*nan
            beta_monthly=np.empty(12)*nan
            pvalue_monthly=np.empty(12)*nan
            
            for tt in range (0,12):
                
                Var_y=np.log(PPint_monthly[tt,:,:])
                Var_x=np.log(Phyc_monthly[tt,:,:])
                Var_y[ np.abs(Var_y)>1e19 ]=nan;  Var_x[ np.abs(Var_x)>1e19 ]=nan
                xx = Var_x.reshape( np.int(Var_x.shape[0]*Var_x.shape[1])) # Creating a time series of the variabel for that month 
                yy = Var_y.reshape( np.int(Var_y.shape[0]*Var_y.shape[1]))             
                yy=yy[np.logical_not(np.isnan(xx))] ; xx=xx[np.logical_not(np.isnan(xx))] # Excluding any NaN values in the time series                   
                xx=xx[np.logical_not(np.isnan(yy))] ; yy=yy[np.logical_not(np.isnan(yy))] # Excluding any NaN values in the time series 
                
                slope_ij, intercept_ij, r_value_ij, p_value_ij, std_err_ij = stats.linregress(xx, yy)
    
                alpha_monthly[tt] = slope_ij
                beta_monthly[tt] = intercept_ij
                pvalue_monthly[tt] = p_value_ij
    
            alpha_monthly = alpha_monthly - 1
            
            plt.plot(P_Var_x, alpha_monthly, c=Colors[M_i], label=GCM, marker='o', markersize=7, markerfacecolor='W', linewidth=2.0)
            plt.xticks(fontsize = 26); plt.yticks(fontsize = 26)
            plt.xticks(P_Var_x, '', rotation=45) # Setting the X-ticks ## rotation='horizontal' # rotation='vertical'
            plt.xticks(P_Var_x, Time_months, rotation=45, fontsize=24) # Setting the X-ticks ## rotation='horizontal' # rotation='vertical'
            plt.xlabel('Month', fontsize=26)
            plt.ylabel('alpha', fontsize=26)
            if iii==0: # 'PAPA'
                plt.title('North Pacific [45N-50N, 140W-150W]', fontsize=24)
            elif iii==1: # 'NABE'
                plt.title('North Atlantic [45N-50N, 25W-35W]', fontsize=24)
            plt.ylim(-2,3)
            plt.grid()

plt.legend(prop={'size': 14}, loc='center right', bbox_to_anchor=(1.45, 0.5), fancybox=True, framealpha=0.8)
plt.subplots_adjust(left=0.1, bottom=0.15, right=0.8, top=0.85, hspace=0.2, wspace=0.25) # the amount of height/width reserved for space between subplots
plt.suptitle(P_title, fontsize=26)   
plt.show()
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full    

if d_trend=='yes':
    fig.savefig(dir_figs+'AllModels_Alpha_value_regress_detrended_NPNA_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
else:
    fig.savefig(dir_figs+'AllModels_Alpha_value_regress_NPNA_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')



from scipy import stats
d_trend='no'

fig=plt.figure() 
P_title= 'beta value [ Ln(PP) = (alpha+1) * Ln(P) + beta ]'
if d_trend=='yes':
    plt.title('Global - '+str(start_date_plt)+'-'+str(end_date_plt)+ ' - detrended', fontsize=24)
else:
    plt.title('Global - '+str(start_date_plt)+'-'+str(end_date_plt), fontsize=24)

for M_i in range(len(GCM_Names)): # M_i=1     # GCM=GCM_Names[0]
    
    GCM=GCM_Names[M_i]
    print (GCM)     
    
    filename_in = (dir_data_in + 'AllResults_'+GCM+'_bio_1991_2010.out') # Directory to save processed data
    my_shelf = shelve.open(filename_in)
    if 'Phyc_allmonths' in my_shelf and 'PPint_allmonths' in my_shelf: # If that variable exists in the saved data    
        globals()['Phyc_allmonths']=my_shelf['Phyc_allmonths']
        globals()['PPint_allmonths']=my_shelf['PPint_allmonths']
        if d_trend=='yes': # If it is requested to detrend the data before calculations
            Phyc_allmonths= func_detrend_3d(Phyc_allmonths)
            PPint_allmonths= func_detrend_3d(PPint_allmonths)
                
        Phyc_monthly = Phyc_allmonths.reshape(( np.int(Phyc_allmonths.shape[0]/12) ,12,lat_n_regrid,lon_n_regrid)) # Monthly means of SST over the time period
        Phyc_monthly = np.nanmean(Phyc_monthly,axis=0)
        PPint_monthly = PPint_allmonths.reshape(( np.int(PPint_allmonths.shape[0]/12) ,12,lat_n_regrid,lon_n_regrid)) # Monthly means of SST over the time period
        PPint_monthly = np.nanmean(PPint_monthly,axis=0)  


        alpha_monthly=np.empty(12)*nan
        beta_monthly=np.empty(12)*nan
        pvalue_monthly=np.empty(12)*nan
        
        for tt in range (0,12):
            
            Var_y=np.log(PPint_monthly[tt,:,:])
            Var_x=np.log(Phyc_monthly[tt,:,:])
            Var_y[ np.abs(Var_y)>1e19 ]=nan;  Var_x[ np.abs(Var_x)>1e19 ]=nan
            xx = Var_x.reshape( np.int(lat_n_regrid*lon_n_regrid)) # Creating a time series of the variabel for that month 
            yy = Var_y.reshape( np.int(lat_n_regrid*lon_n_regrid))             
            yy=yy[np.logical_not(np.isnan(xx))] ; xx=xx[np.logical_not(np.isnan(xx))] # Excluding any NaN values in the time series                   
            xx=xx[np.logical_not(np.isnan(yy))] ; yy=yy[np.logical_not(np.isnan(yy))] # Excluding any NaN values in the time series 
            
            slope_ij, intercept_ij, r_value_ij, p_value_ij, std_err_ij = stats.linregress(xx, yy)

            alpha_monthly[tt] = slope_ij
            beta_monthly[tt] = intercept_ij
            pvalue_monthly[tt] = p_value_ij

        alpha_monthly = alpha_monthly - 1
        
        plt.plot(P_Var_x, beta_monthly, c=Colors[M_i], label=GCM, marker='o', markersize=7, markerfacecolor='W', linewidth=2.0)
        plt.xticks(fontsize = 26); plt.yticks(fontsize = 26)
        #plt.xlabel('Month', fontsize=26)
        plt.ylabel('beta', fontsize=26)
        #plt.ylim(275,450)
    
plt.xticks(P_Var_x, '', rotation=45) # Setting the X-ticks ## rotation='horizontal' # rotation='vertical'
plt.legend(prop={'size': 14}, loc='center right', bbox_to_anchor=(1.4, 0.5), fancybox=True, framealpha=0.8)
plt.xticks(P_Var_x, Time_months, rotation=45, fontsize=24) # Setting the X-ticks ## rotation='horizontal' # rotation='vertical'
plt.subplots_adjust(left=0.33, bottom=0.11, right=0.69, top=0.9, hspace=0.2, wspace=0.05) # the amount of height/width reserved for space between subplots
plt.suptitle(P_title, fontsize=26)    
plt.show()
plt.grid()
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full    

if d_trend=='yes':
    fig.savefig(dir_figs+'AllModels_Beta_value_regress_detrended_Global_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
else:
    fig.savefig(dir_figs+'AllModels_Beta_value_regress_Global_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


#######################################################################################################################################################
#######################################################################################################################################################

from scipy import stats
d_trend='no'

fig=plt.figure() 
P_title= 'alpha value [ Ln(PP) = (alpha+1) * Ln(P) + beta ] - Subtropical Gyres of N.P. and N.A. - '+str(start_date_plt)+'-'+str(end_date_plt)

for iii in range(0,2):
    ax = fig.add_subplot(1,2,iii+1)

    for M_i in range(len(GCM_Names)): # M_i=1     # GCM=GCM_Names[0]
        
        GCM=GCM_Names[M_i]
        print (GCM)     
        
        filename_in = (dir_data_in + 'AllResults_'+GCM+'_bio_1991_2010.out') # Directory to save processed data
        my_shelf = shelve.open(filename_in)
        if 'Phyc_allmonths' in my_shelf and 'PPint_allmonths' in my_shelf: # If that variable exists in the saved data    
            globals()['Phyc_allmonths']=my_shelf['Phyc_allmonths']
            globals()['PPint_allmonths']=my_shelf['PPint_allmonths']
                
            if iii==0: # 'PAPA'
                Phyc_allmonths_BOX = Phyc_allmonths[:,110:131,160:230]
                PPint_allmonths_BOX = PPint_allmonths[:,110:131,160:230]
            elif iii==1: # 'NABE'
                Phyc_allmonths_BOX = Phyc_allmonths[:,110:131,290:350]
                PPint_allmonths_BOX = PPint_allmonths[:,110:131,290:350]

            if d_trend=='yes': # If it is requested to detrend the data before calculations
                Phyc_allmonths_BOX= func_detrend_3d(Phyc_allmonths_BOX)
                PPint_allmonths_BOX= func_detrend_3d(PPint_allmonths_BOX)                               
                    
            Phyc_monthly = Phyc_allmonths_BOX.reshape(( np.int(Phyc_allmonths_BOX.shape[0]/12) ,12,Phyc_allmonths_BOX.shape[1],Phyc_allmonths_BOX.shape[2])) # Monthly means of SST over the time period
            Phyc_monthly = np.nanmean(Phyc_monthly,axis=0)
            PPint_monthly = PPint_allmonths_BOX.reshape(( np.int(PPint_allmonths_BOX.shape[0]/12) ,12,PPint_allmonths_BOX.shape[1],PPint_allmonths_BOX.shape[2])) # Monthly means of SST over the time period
            PPint_monthly = np.nanmean(PPint_monthly,axis=0)  
    
            alpha_monthly=np.empty(12)*nan
            beta_monthly=np.empty(12)*nan
            pvalue_monthly=np.empty(12)*nan
            
            for tt in range (0,12):
                
                Var_y=np.log(PPint_monthly[tt,:,:])
                Var_x=np.log(Phyc_monthly[tt,:,:])
                Var_y[ np.abs(Var_y)>1e19 ]=nan;  Var_x[ np.abs(Var_x)>1e19 ]=nan
                xx = Var_x.reshape( np.int(Var_x.shape[0]*Var_x.shape[1])) # Creating a time series of the variabel for that month 
                yy = Var_y.reshape( np.int(Var_y.shape[0]*Var_y.shape[1]))             
                yy=yy[np.logical_not(np.isnan(xx))] ; xx=xx[np.logical_not(np.isnan(xx))] # Excluding any NaN values in the time series                   
                xx=xx[np.logical_not(np.isnan(yy))] ; yy=yy[np.logical_not(np.isnan(yy))] # Excluding any NaN values in the time series 
                
                slope_ij, intercept_ij, r_value_ij, p_value_ij, std_err_ij = stats.linregress(xx, yy)
    
                alpha_monthly[tt] = slope_ij
                beta_monthly[tt] = intercept_ij
                pvalue_monthly[tt] = p_value_ij
    
            alpha_monthly = alpha_monthly - 1
            
            plt.plot(P_Var_x, alpha_monthly, c=Colors[M_i], label=GCM, marker='o', markersize=7, markerfacecolor='W', linewidth=2.0)
            plt.xticks(fontsize = 26); plt.yticks(fontsize = 26)
            plt.xticks(P_Var_x, '', rotation=45) # Setting the X-ticks ## rotation='horizontal' # rotation='vertical'
            plt.xticks(P_Var_x, Time_months, rotation=45, fontsize=24) # Setting the X-ticks ## rotation='horizontal' # rotation='vertical'
            plt.xlabel('Month', fontsize=26)
            plt.ylabel('alpha', fontsize=26)
            if iii==0: # 'PAPA'
                plt.title('SubTrop.Gyre N.Pac.[20N-40N,160E-130W]', fontsize=24)
            elif iii==1: # 'NABE'
                plt.title('SubTrop.Gyre N.Atl.[20N-40N,10W-70W]', fontsize=24)
            plt.ylim(-2,3)
            plt.grid()

plt.legend(prop={'size': 14}, loc='center right', bbox_to_anchor=(1.45, 0.5), fancybox=True, framealpha=0.8)
plt.subplots_adjust(left=0.1, bottom=0.15, right=0.8, top=0.85, hspace=0.2, wspace=0.25) # the amount of height/width reserved for space between subplots
plt.suptitle(P_title, fontsize=26)   
plt.show()
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full    

if d_trend=='yes':
    fig.savefig(dir_figs+'AllModels_Alpha_value_regress_detrended_SubTropGyreNPNA_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
else:
    fig.savefig(dir_figs+'AllModels_Alpha_value_regress_SubTropGyreNPNA_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')



d_trend='no'

fig=plt.figure() 
P_title= 'beta value [ Ln(PP) = (alpha+1) * Ln(P) + beta ] - Subtropical Gyres of N.P. and N.A. - '+str(start_date_plt)+'-'+str(end_date_plt)

for iii in range(0,2):
    ax = fig.add_subplot(1,2,iii+1)

    for M_i in range(len(GCM_Names)): # M_i=1     # GCM=GCM_Names[0]
        
        GCM=GCM_Names[M_i]
        print (GCM)     
        
        filename_in = (dir_data_in + 'AllResults_'+GCM+'_bio_1991_2010.out') # Directory to save processed data
        my_shelf = shelve.open(filename_in)
        if 'Phyc_allmonths' in my_shelf and 'PPint_allmonths' in my_shelf: # If that variable exists in the saved data    
            globals()['Phyc_allmonths']=my_shelf['Phyc_allmonths']
            globals()['PPint_allmonths']=my_shelf['PPint_allmonths']
                
            if iii==0: # 'PAPA'
                Phyc_allmonths_BOX = Phyc_allmonths[:,110:131,160:230]
                PPint_allmonths_BOX = PPint_allmonths[:,110:131,160:230]
            elif iii==1: # 'NABE'
                Phyc_allmonths_BOX = Phyc_allmonths[:,110:131,290:350]
                PPint_allmonths_BOX = PPint_allmonths[:,110:131,290:350]

            if d_trend=='yes': # If it is requested to detrend the data before calculations
                Phyc_allmonths_BOX= func_detrend_3d(Phyc_allmonths_BOX)
                PPint_allmonths_BOX= func_detrend_3d(PPint_allmonths_BOX)                               
                    
            Phyc_monthly = Phyc_allmonths_BOX.reshape(( np.int(Phyc_allmonths_BOX.shape[0]/12) ,12,Phyc_allmonths_BOX.shape[1],Phyc_allmonths_BOX.shape[2])) # Monthly means of SST over the time period
            Phyc_monthly = np.nanmean(Phyc_monthly,axis=0)
            PPint_monthly = PPint_allmonths_BOX.reshape(( np.int(PPint_allmonths_BOX.shape[0]/12) ,12,PPint_allmonths_BOX.shape[1],PPint_allmonths_BOX.shape[2])) # Monthly means of SST over the time period
            PPint_monthly = np.nanmean(PPint_monthly,axis=0)  
    
            alpha_monthly=np.empty(12)*nan
            beta_monthly=np.empty(12)*nan
            pvalue_monthly=np.empty(12)*nan
            
            for tt in range (0,12):
                
                Var_y=np.log(PPint_monthly[tt,:,:])
                Var_x=np.log(Phyc_monthly[tt,:,:])
                Var_y[ np.abs(Var_y)>1e19 ]=nan;  Var_x[ np.abs(Var_x)>1e19 ]=nan
                xx = Var_x.reshape( np.int(Var_x.shape[0]*Var_x.shape[1])) # Creating a time series of the variabel for that month 
                yy = Var_y.reshape( np.int(Var_y.shape[0]*Var_y.shape[1]))             
                yy=yy[np.logical_not(np.isnan(xx))] ; xx=xx[np.logical_not(np.isnan(xx))] # Excluding any NaN values in the time series                   
                xx=xx[np.logical_not(np.isnan(yy))] ; yy=yy[np.logical_not(np.isnan(yy))] # Excluding any NaN values in the time series 
                
                slope_ij, intercept_ij, r_value_ij, p_value_ij, std_err_ij = stats.linregress(xx, yy)
    
                alpha_monthly[tt] = slope_ij
                beta_monthly[tt] = intercept_ij
                pvalue_monthly[tt] = p_value_ij
    
            alpha_monthly = alpha_monthly - 1
            
            plt.plot(P_Var_x, beta_monthly, c=Colors[M_i], label=GCM, marker='o', markersize=7, markerfacecolor='W', linewidth=2.0)
            plt.xticks(fontsize = 26); plt.yticks(fontsize = 26)
            plt.xticks(P_Var_x, '', rotation=45) # Setting the X-ticks ## rotation='horizontal' # rotation='vertical'
            plt.xticks(P_Var_x, Time_months, rotation=45, fontsize=24) # Setting the X-ticks ## rotation='horizontal' # rotation='vertical'
            plt.xlabel('Month', fontsize=26)
            plt.ylabel('beta', fontsize=26)
            if iii==0: # 'PAPA'
                plt.title('SubTrop.Gyre N.Pac.[20N-40N,160E-130W]', fontsize=24)
            elif iii==1: # 'NABE'
                plt.title('SubTrop.Gyre N.Atl.[20N-40N,10W-70W]', fontsize=24)
            #plt.ylim(-2,3)
            plt.grid()

plt.legend(prop={'size': 14}, loc='center right', bbox_to_anchor=(1.45, 0.5), fancybox=True, framealpha=0.8)
plt.subplots_adjust(left=0.1, bottom=0.15, right=0.8, top=0.85, hspace=0.2, wspace=0.25) # the amount of height/width reserved for space between subplots
plt.suptitle(P_title, fontsize=26)   
plt.show()
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full    

if d_trend=='yes':
    fig.savefig(dir_figs+'AllModels_Beta_value_regress_detrended_SubTropGyreNPNA_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
else:
    fig.savefig(dir_figs+'AllModels_Beta_value_regress_SubTropGyreNPNA_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')



###############################################################################
###############################################################################
#######    Alpha values using grazing data - GFDL-ESM2G model only    #########
###############################################################################
###############################################################################
from scipy import stats
d_trend='no'

dir_data_in_f=('/data5/scratch/Behzad/ocean/ocean_biogeochemistry/processed_data/python_data_CMIP5_bio_npfiles_1991_2010/')

GCM='GFDL-ESM2G'

Phyc_allmonths = np.load( dir_data_in_f + 'Phyc_allmonths_'+GCM+'_'+str(start_date_plt)+'_'+str(end_date_plt)+'.npy' )
Graz_allmonths = np.load( dir_data_in_f + 'Graz_allmonths_'+GCM+'_'+str(start_date_plt)+'_'+str(end_date_plt)+'.npy' )
PP_allmonths = np.load( dir_data_in_f + 'PP_allmonths_'+GCM+'_'+str(start_date_plt)+'_'+str(end_date_plt)+'.npy' )
NO3_allmonths = np.load( dir_data_in_f + 'NO3_allmonths_'+GCM+'_'+str(start_date_plt)+'_'+str(end_date_plt)+'.npy' )
Phypico_allmonths = np.load( dir_data_in_f + 'Phypico_allmonths_'+GCM+'_'+str(start_date_plt)+'_'+str(end_date_plt)+'.npy' )
Zooc_allmonths = np.load( dir_data_in_f + 'Zooc_allmonths_'+GCM+'_'+str(start_date_plt)+'_'+str(end_date_plt)+'.npy' )

NO3_allmonths_ave=np.nanmean(NO3_allmonths,axis=0)
Nit_mask= empty((NO3_allmonths_ave.shape[0], NO3_allmonths_ave.shape[1])) * nan
Nit_mask [NO3_allmonths_ave*1000 >= 1.18 ] = 1
Nit_mask [NO3_allmonths_ave*1000 < 1.18 ] = 0
#Nit_mask [ Ocean_Land_mask == 0 ] = nan # masking over land, so grid cells that fall on land area (value=0) will be deleted

# iii=3     #     'navy', 'blue', 'royalblue', 'lime', 'limegreen', 'green', 'yellow', 'gold', 'orange', 'red', 'magenta','darkviolet'
Colors    =  [ 'lime', 'green', 'yellow', 'gold', 'navy', 'royalblue', 'red', 'magenta','darkviolet', 'indigo', 'cyan', 'steelblue'   ]
plot_regions=[ 'Subtropical N Pacific - [20N-35N,160E-130W]', 'Subtropical N Atlantic - [20N-35N,70W-20W]', 'Subtropical S Pacific - [20S-35S,170W-90W]', 'Subtropical S Atlantic - [20S-35S,30W-10W]', 'Mid-Lat. N. Pacific PAPA [45N-55N, 150W-140W]', 'Mid-Lat. N. Atlantic NABE [45N-55N, 35W-25W]', 'Pacific Equatorial tongue - [30S-30N, 140W-65W ]', 'Peru Coast, Pacific - [55S-20S, 80W-69W ]', 'Argentina Coast, Atlantic - [55S-30S, 69W-50W ]', 'SouthWest Africa Coast, Atlantic - [35S-22S, 8E-20E ]', 'Subtpolar S Pacific - [60S-40S,170W-90W]', 'Subpolar S Atlantic - [60S-40S,20W-0W]'    ]                                                        
plot_lats=[ [110,126],  [110,126], [55,71],    [55,71],    [135,140], [135,140], [60,121],  [35,70],   [35,60],   [55,68], [30,51],    [30,51] ]
plot_lons=[ [160,230],  [290,340], [190,270],  [330,350],  [210,220], [325,335], [140,295], [280,291], [291,311], [8,21],  [190,270],  [330,350] ]

fig=plt.figure()
for iii in range(0,12):
    # iii==0  Subtropical Gyre  N Pacific - [20N-40N,160E-130W]
    # iii==1  Subtropical Gyre N Atlantic - [20N-40N,10W-70W]
    # iii==2  North Pacific Mid-Latitudes PAPA [45N-55N, 140W-150W]
    # iii==3  North Atlantic Mid-Latitudes NABE [45N-55N, 25W-35W]    
    if iii!=6:
        Phyc_allmonths_BOX = Phyc_allmonths[:,plot_lats[iii][0]:plot_lats[iii][1], plot_lons[iii][0]:plot_lons[iii][1]]
        Graz_allmonths_BOX = Graz_allmonths[:,plot_lats[iii][0]:plot_lats[iii][1], plot_lons[iii][0]:plot_lons[iii][1]]
        PP_allmonths_BOX = PP_allmonths[:,plot_lats[iii][0]:plot_lats[iii][1], plot_lons[iii][0]:plot_lons[iii][1]]
        Phypico_allmonths_BOX = Phypico_allmonths[:,plot_lats[iii][0]:plot_lats[iii][1], plot_lons[iii][0]:plot_lons[iii][1]]
        Zooc_allmonths_BOX = Zooc_allmonths[:,plot_lats[iii][0]:plot_lats[iii][1], plot_lons[iii][0]:plot_lons[iii][1]]
        
    elif iii==6: # Pacific Equatorial tongue - [30S-30N, 140E-65W ]      
        Phyc_ave_EqTongue = copy.deepcopy( np.nanmean(Phyc_allmonths,axis=0) )
        Phyc_ave_EqTongue [Nit_mask == 0] = nan # Excluding Nitrate < 1.18 regions
        Phyc_ave_EqTongue [Ocean_Index != 2] = nan # Excluding all but Pacific Ocean
        Phyc_allmonths_BOX_ave=Phyc_ave_EqTongue[60:121,220:295] # Pacific Equatorial tongue [30S-30N, 140W-65W ]
        
        PP_ave_EqTongue = copy.deepcopy( np.nanmean(PP_allmonths,axis=0) )
        PP_ave_EqTongue [Nit_mask == 0] = nan ; PP_ave_EqTongue [Ocean_Index != 2] = nan # Excluding all but Pacific Ocean
        PP_allmonths_BOX_ave=PP_ave_EqTongue[60:121,220:295] # Pacific Equatorial tongue [30S-30N, 140W-65W ]
        
        Graz_ave_EqTongue = copy.deepcopy( np.nanmean(Graz_allmonths,axis=0) )
        Graz_ave_EqTongue [Nit_mask == 0] = nan ; Graz_ave_EqTongue [Ocean_Index != 2] = nan # Excluding all but Pacific Ocean
        Graz_allmonths_BOX_ave=Graz_ave_EqTongue[60:121,220:295] # Pacific Equatorial tongue [30S-30N, 140W-65W ]    

        Phypico_ave_EqTongue = copy.deepcopy( np.nanmean(Phypico_allmonths,axis=0) )
        Phypico_ave_EqTongue [Nit_mask == 0] = nan ; Phypico_ave_EqTongue [Ocean_Index != 2] = nan # Excluding all but Pacific Ocean
        Phypico_allmonths_BOX_ave=Phypico_ave_EqTongue[60:121,220:295] # Pacific Equatorial tongue [30S-30N, 140W-65W ]

        Zooc_ave_EqTongue = copy.deepcopy( np.nanmean(Zooc_allmonths,axis=0) )
        Zooc_ave_EqTongue [Nit_mask == 0] = nan ; Zooc_ave_EqTongue [Ocean_Index != 2] = nan # Excluding all but Pacific Ocean
        Zooc_allmonths_BOX_ave=Phypico_ave_EqTongue[60:121,220:295] # Pacific Equatorial tongue [30S-30N, 140W-65W ]         
    
    if d_trend=='yes' and iii !=6 : # If it is requested to detrend the data before calculations
        Phyc_allmonths_BOX= func_detrend_3d(Phyc_allmonths_BOX)
        Graz_allmonths_BOX= func_detrend_3d(Graz_allmonths_BOX)
        PP_allmonths_BOX= func_detrend_3d(PP_allmonths_BOX)
        Phypico_allmonths_BOX= func_detrend_3d(Phypico_allmonths_BOX)
        Zooc_allmonths_BOX= func_detrend_3d(Zooc_allmonths_BOX)

    PhytoLarge_allmonths_BOX = Phyc_allmonths_BOX - Phypico_allmonths_BOX
    PhytoLarge_prcnt_allmonths_BOX = ( (Phyc_allmonths_BOX - Phypico_allmonths_BOX) / Phyc_allmonths_BOX ) *100
        
    if iii!=6:
        Phyc_allmonths_BOX_ave=np.nanmean(Phyc_allmonths_BOX,axis=0)
        Graz_allmonths_BOX_ave=np.nanmean(Graz_allmonths_BOX,axis=0)
        PP_allmonths_BOX_ave=np.nanmean(PP_allmonths_BOX,axis=0)
        Phypico_allmonths_BOX_ave=np.nanmean(Phypico_allmonths_BOX,axis=0)
        Zooc_allmonths_BOX_ave=np.nanmean(Zooc_allmonths_BOX,axis=0)
        PhytoLarge_allmonths_BOX_ave=np.nanmean(PhytoLarge_allmonths_BOX,axis=0)
        PhytoLarge_prcnt_allmonths_BOX_ave=np.nanmean(PhytoLarge_prcnt_allmonths_BOX,axis=0)
    
    #######################################
    Var_y=np.log(Graz_allmonths_BOX_ave)
    Var_x=np.log(Phyc_allmonths_BOX_ave)
    Var_y[ np.abs(Var_y)>1e19 ]=nan;  Var_x[ np.abs(Var_x)>1e19 ]=nan
    xx = Var_x.reshape( np.int(Var_x.shape[0]*Var_x.shape[1])) # Creating a time series of the variabel for that month 
    yy = Var_y.reshape( np.int(Var_y.shape[0]*Var_y.shape[1]))             
    yy=yy[np.logical_not(np.isnan(xx))] ; xx=xx[np.logical_not(np.isnan(xx))] # Excluding any NaN values in the time series                   
    xx=xx[np.logical_not(np.isnan(yy))] ; yy=yy[np.logical_not(np.isnan(yy))] # Excluding any NaN values in the time series 
    
    alpha_G_P_ave, beta_G_P_ave, r_value_ij, pvalue_G_P_ave, std_err_ij = stats.linregress(xx, yy)
    alpha_G_P_ave = alpha_G_P_ave - 1 ; print(alpha_G_P_ave)
    
    #######################################
    Var_y=np.log(PP_allmonths_BOX_ave)
    Var_x=np.log(Phyc_allmonths_BOX_ave)
    Var_y[ np.abs(Var_y)>1e19 ]=nan;  Var_x[ np.abs(Var_x)>1e19 ]=nan
    xx = Var_x.reshape( np.int(Var_x.shape[0]*Var_x.shape[1])) # Creating a time series of the variabel for that month 
    yy = Var_y.reshape( np.int(Var_y.shape[0]*Var_y.shape[1]))             
    yy=yy[np.logical_not(np.isnan(xx))] ; xx=xx[np.logical_not(np.isnan(xx))] # Excluding any NaN values in the time series                   
    xx=xx[np.logical_not(np.isnan(yy))] ; yy=yy[np.logical_not(np.isnan(yy))] # Excluding any NaN values in the time series 
    
    alpha_PP_P_ave, beta_PP_P_ave, r_value_ij, pvalue_PP_P_ave, std_err_ij = stats.linregress(xx, yy)
    alpha_PP_P_ave = alpha_PP_P_ave - 1 ; print(alpha_PP_P_ave)
    #######################################
    
    print( np.nanmean(np.nanmean(Graz_allmonths_BOX_ave,axis=0),axis=0) )
    print( np.nanmean(np.nanmean(PP_allmonths_BOX_ave,axis=0),axis=0) )
    
    
    n_r=3 ; n_c=3 ;

    ax = fig.add_subplot(n_r,n_c,1) 
    plt.scatter( (np.nanmean(np.nanmean(Graz_allmonths_BOX_ave,axis=0),axis=0))*1E9 , alpha_G_P_ave, s=100, marker='o', c=Colors[iii], label=plot_regions[iii]) # Only the first pannel gets the lables so the legend location be alligned with that
    plt.scatter( (np.nanmean(np.nanmean(Graz_allmonths_BOX_ave,axis=0),axis=0))*1E9 , alpha_G_P_ave, s=100, marker='o', facecolors='none', edgecolors='k')
    plt.scatter( (np.nanmean(np.nanmean(Graz_allmonths_BOX_ave,axis=0),axis=0))*1E9 , alpha_PP_P_ave, s=100, marker='d', c=Colors[iii]) # Only the first pannel gets the lables so the legend location be alligned with that
    plt.scatter( (np.nanmean(np.nanmean(Graz_allmonths_BOX_ave,axis=0),axis=0))*1E9 , alpha_PP_P_ave, s=100, marker='d', facecolors='none', edgecolors='k')
    plt.title('alpha vs. total grazing', fontsize=14)
    plt.xlabel('total grazing (nano mol . m-2 . s-1)',fontsize=14); plt.ylabel('alpha',fontsize=14)  
    plt.xlim(0,35) ; plt.ylim(0,1.2)

    ax = fig.add_subplot(n_r,n_c,2) 
    plt.scatter( (np.nanmean(np.nanmean(PP_allmonths_BOX_ave,axis=0),axis=0))*1E9, alpha_G_P_ave, s=100, marker='o', c=Colors[iii], label=plot_regions[iii]) # Only the first pannel gets the lables so the legend location be alligned with that
    plt.scatter((np.nanmean(np.nanmean(PP_allmonths_BOX_ave,axis=0),axis=0))*1E9, alpha_G_P_ave, s=100, marker='o', facecolors='none', edgecolors='k')    
    plt.scatter( (np.nanmean(np.nanmean(PP_allmonths_BOX_ave,axis=0),axis=0))*1E9, alpha_PP_P_ave, s=100, marker='d', c=Colors[iii]) # Only the first pannel gets the lables so the legend location be alligned with that
    plt.scatter((np.nanmean(np.nanmean(PP_allmonths_BOX_ave,axis=0),axis=0))*1E9, alpha_PP_P_ave, s=100, marker='d', facecolors='none', edgecolors='k')
    plt.title('alpha vs. primary production', fontsize=14)   
    plt.xlabel('Promary Production (nano mol . m-2 . s-1)',fontsize=14); plt.ylabel('alpha',fontsize=14)  
    plt.xlim(0,35) ; plt.ylim(0,1.2)

    ax = fig.add_subplot(n_r,n_c,3) 
    plt.scatter( (np.nanmean(np.nanmean(Phyc_allmonths_BOX_ave,axis=0),axis=0))*1E3, alpha_G_P_ave, s=100, marker='o', c=Colors[iii], label=plot_regions[iii]) # Only the first pannel gets the lables so the legend location be alligned with that
    plt.scatter((np.nanmean(np.nanmean(Phyc_allmonths_BOX_ave,axis=0),axis=0))*1E3, alpha_G_P_ave, s=100, marker='o', facecolors='none', edgecolors='k')
    plt.scatter( (np.nanmean(np.nanmean(Phyc_allmonths_BOX_ave,axis=0),axis=0))*1E3, alpha_PP_P_ave, s=100, marker='d', c=Colors[iii]) # Only the first pannel gets the lables so the legend location be alligned with that
    plt.scatter((np.nanmean(np.nanmean(Phyc_allmonths_BOX_ave,axis=0),axis=0))*1E3, alpha_PP_P_ave, s=100, marker='d', facecolors='none', edgecolors='k')
    plt.title('alpha vs. phyto', fontsize=14)   
    plt.xlabel('Phytoplankton biomass (mili mol . m-3)',fontsize=14); plt.ylabel('alpha',fontsize=14)  
    plt.xlim(0.8,4) ; plt.ylim(0,1.2)

    ax = fig.add_subplot(n_r,n_c,4) 
    plt.scatter( (np.nanmean(np.nanmean(PhytoLarge_prcnt_allmonths_BOX_ave,axis=0),axis=0)), alpha_G_P_ave, s=100, marker='o', c=Colors[iii], label=plot_regions[iii]) # Only the first pannel gets the lables so the legend location be alligned with that
    plt.scatter((np.nanmean(np.nanmean(PhytoLarge_prcnt_allmonths_BOX_ave,axis=0),axis=0)), alpha_G_P_ave, s=100, marker='o', facecolors='none', edgecolors='k')
    plt.scatter( (np.nanmean(np.nanmean(PhytoLarge_prcnt_allmonths_BOX_ave,axis=0),axis=0)), alpha_PP_P_ave, s=100, marker='d', c=Colors[iii]) # Only the first pannel gets the lables so the legend location be alligned with that
    plt.scatter((np.nanmean(np.nanmean(PhytoLarge_prcnt_allmonths_BOX_ave,axis=0),axis=0)), alpha_PP_P_ave, s=100, marker='d', facecolors='none', edgecolors='k')
    plt.title('alpha vs. large phyto %', fontsize=14)   
    plt.xlabel('Large Phytoplankton percentage (%)',fontsize=14); plt.ylabel('alpha',fontsize=14)  
    plt.xlim(5,35) ; plt.ylim(0,1.2)
    
    ax = fig.add_subplot(n_r,n_c,5) 
    plt.scatter( (np.nanmean(np.nanmean(PhytoLarge_allmonths_BOX_ave,axis=0),axis=0))*1E3, alpha_G_P_ave, s=100, marker='o', c=Colors[iii], label=plot_regions[iii]) # Only the first pannel gets the lables so the legend location be alligned with that
    plt.scatter((np.nanmean(np.nanmean(PhytoLarge_allmonths_BOX_ave,axis=0),axis=0))*1E3, alpha_G_P_ave, s=100, marker='o', facecolors='none', edgecolors='k')
    plt.scatter( (np.nanmean(np.nanmean(PhytoLarge_allmonths_BOX_ave,axis=0),axis=0))*1E3, alpha_PP_P_ave, s=100, marker='d', c=Colors[iii]) # Only the first pannel gets the lables so the legend location be alligned with that
    plt.scatter((np.nanmean(np.nanmean(PhytoLarge_allmonths_BOX_ave,axis=0),axis=0))*1E3, alpha_PP_P_ave, s=100, marker='d', facecolors='none', edgecolors='k')
    plt.title('alpha vs. large phyto', fontsize=14)   
    plt.xlabel('Large Phytoplankton (mili mol . m-3)',fontsize=14); plt.ylabel('alpha',fontsize=14)
    plt.xlim(0,1.8) ; plt.ylim(0,1.2)

    ax = fig.add_subplot(n_r,n_c,6) 
    plt.scatter( (np.nanmean(np.nanmean(Phypico_allmonths_BOX_ave,axis=0),axis=0))*1E3, alpha_G_P_ave, s=100, marker='o', c=Colors[iii], label=plot_regions[iii]) # Only the first pannel gets the lables so the legend location be alligned with that
    plt.scatter((np.nanmean(np.nanmean(Phypico_allmonths_BOX_ave,axis=0),axis=0))*1E3, alpha_G_P_ave, s=100, marker='o', facecolors='none', edgecolors='k')
    plt.scatter( (np.nanmean(np.nanmean(Phypico_allmonths_BOX_ave,axis=0),axis=0))*1E3, alpha_PP_P_ave, s=100, marker='d', c=Colors[iii]) # Only the first pannel gets the lables so the legend location be alligned with that
    plt.scatter((np.nanmean(np.nanmean(Phypico_allmonths_BOX_ave,axis=0),axis=0))*1E3, alpha_PP_P_ave, s=100, marker='d', facecolors='none', edgecolors='k')
    plt.title('alpha vs. small phyto', fontsize=14)   
    plt.xlabel('Small Phytoplankton (mili mol . m-3)',fontsize=14); plt.ylabel('alpha',fontsize=14)
    plt.xlim(0.5,2.5) ; plt.ylim(0,1.2)

    ax = fig.add_subplot(n_r,n_c,7) 
    plt.scatter( (np.nanmean(np.nanmean(Zooc_allmonths_BOX_ave,axis=0),axis=0))*1E3, alpha_G_P_ave, s=100, marker='o', c=Colors[iii], label=plot_regions[iii]) # Only the first pannel gets the lables so the legend location be alligned with that
    plt.scatter((np.nanmean(np.nanmean(Zooc_allmonths_BOX_ave,axis=0),axis=0))*1E3, alpha_G_P_ave, s=100, marker='o', facecolors='none', edgecolors='k')
    plt.scatter( (np.nanmean(np.nanmean(Zooc_allmonths_BOX_ave,axis=0),axis=0))*1E3, alpha_PP_P_ave, s=100, marker='d', c=Colors[iii]) # Only the first pannel gets the lables so the legend location be alligned with that
    plt.scatter((np.nanmean(np.nanmean(Zooc_allmonths_BOX_ave,axis=0),axis=0))*1E3, alpha_PP_P_ave, s=100, marker='d', facecolors='none', edgecolors='k')
    plt.title('alpha vs. zooplankton', fontsize=14)   
    plt.xlabel('Zooplankton biomass (mili mol . m-3)',fontsize=14); plt.ylabel('alpha',fontsize=14)  
    plt.xlim(0.8,8) ; plt.ylim(0,1.2)
        
    #plt.plot([-1000,1000],[-1000,1000], 'k--', linewidth=0.75)
    #plt.xlim(plot_ranges[V_i][0],plot_ranges[V_i][1]) ; plt.ylim(plot_ranges[V_i][0],plot_ranges[V_i][1])
    plt.xticks(fontsize=10); plt.yticks(fontsize = 10) # Setting the X-ticks ## rotation='horizontal' # rotation='vertical'   

plt.subplots_adjust(left=0.1, bottom=0.07, right=0.85, top=0.9, hspace=0.45, wspace=0.25) # the amount of height/width reserved for space between subplots
plt.legend(shadow=True, loc='right', bbox_to_anchor=(4, 0.35), markerscale=1 , ncol=2, prop={'size': 14})
if d_trend=='yes':
    plt.suptitle('Circles: Ln(Grazing) = (alpha+1)*Ln(P) + beta  ,  Diamonds: Ln(PP) = (alpha+1)*Ln(P) + beta  - '+str(start_date_plt)+'-'+str(end_date_plt)+' - '+str(GCM)+' - detrended' , fontsize=18)
else:
    plt.suptitle('Circles: Ln(Grazing) = (alpha+1)*Ln(P) + beta  ,  Diamonds: Ln(PP) = (alpha+1)*Ln(P) + beta - '+str(start_date_plt)+'-'+str(end_date_plt)+' - '+str(GCM), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full

if d_trend=='yes':
    fig.savefig(dir_figs+str(GCM)+'_alpha_regress_scatter_detrended_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
else:
    fig.savefig(dir_figs+str(GCM)+'_alpha_regress_scatter_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')





###############################################################################
from BehzadlibPlot import func_plotmap_contourf

Plot_Var = np.nanmean(Graz_allmonths,axis=0)
cmap_limit=np.nanmax(np.abs( np.nanpercentile(Plot_Var, 99.99)))
#Plot_range=np.linspace(0,cmap_limit,27)
Plot_range=np.linspace(0,8E-8,41)
Plot_unit='(mol m-3 s-1)'; Plot_title= 'Surface Total Grazing of Phytoplankton (mol m-3 s-1) - '+str(start_date_plt)+'-'+str(end_date_plt)+' - '+str(GCM)
fig, m = func_plotmap_contourf(Plot_Var, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.jet, 'cyl', 210., 80., -80., '-')
fig.savefig(dir_figs+str(GCM)+'_ave_map_Grazing.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

Plot_Var = np.nanmean(PP_allmonths,axis=0)
cmap_limit=np.nanmax(np.abs( np.nanpercentile(Plot_Var, 99.99)))
#Plot_range=np.linspace(0,cmap_limit,27)
Plot_range=np.linspace(0,8E-8,41)
Plot_unit='(mol m-2 s-1)'; Plot_title= 'Surface Primary Production (mol m-3 s-1) - '+str(start_date_plt)+'-'+str(end_date_plt)+' - '+str(GCM)
fig, m = func_plotmap_contourf(Plot_Var, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.jet, 'cyl', 210., 80., -80., '-')
fig.savefig(dir_figs+str(GCM)+'_ave_map_PP.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

Plot_Var = np.nanmean(PPint_allmonths,axis=0)
cmap_limit=np.nanmax(np.abs( np.nanpercentile(Plot_Var, 99.99)))
#Plot_range=np.linspace(0,cmap_limit,41)
Plot_range=np.linspace(0,2E-6,41)
Plot_unit='(mol m-2 s-1)'; Plot_title= 'Integrated Primary Production (mol m-2 s-1) - '+str(start_date_plt)+'-'+str(end_date_plt)+' - '+str(GCM)
fig, m = func_plotmap_contourf(Plot_Var, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.jet, 'cyl', 210., 80., -80., '-')
fig.savefig(dir_figs+str(GCM)+'_ave_map_PPint.png', format='png', dpi=300, transparent=True, bbox_inches='tight')




Plot_Var = np.nanmean(Phyc_allmonths,axis=0)
cmap_limit=np.nanmax(np.abs( np.nanpercentile(Plot_Var, 99.99)))
#Plot_range=np.linspace(0,cmap_limit,41)
Plot_range=np.linspace(0,0.0055,56)
Plot_unit='(mol carbon / m3)'; Plot_title= 'Total Phytoplankton at Surface (mol m-3) - '+str(start_date_plt)+'-'+str(end_date_plt)+' - '+str(GCM)
fig, m = func_plotmap_contourf(Plot_Var, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.jet, 'cyl', 210., 80., -80., '-')
fig.savefig(dir_figs+str(GCM)+'_ave_map_Phyc.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

Plot_Var = np.nanmean(Phypico_allmonths,axis=0)
cmap_limit=np.nanmax(np.abs( np.nanpercentile(Plot_Var, 99.99)))
#Plot_range=np.linspace(0,cmap_limit,41)
Plot_range=np.linspace(0,0.0055,56)
Plot_unit='(mol carbon / m3)'; Plot_title= 'Small Phytoplankton (phypico) at Surface (mol m-3) - '+str(start_date_plt)+'-'+str(end_date_plt)+' - '+str(GCM)
fig, m = func_plotmap_contourf(Plot_Var, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.jet, 'cyl', 210., 80., -80., '-')
fig.savefig(dir_figs+str(GCM)+'_ave_map_PhytoSmall.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

Plot_Var = np.nanmean(Phyc_allmonths - Phypico_allmonths,axis=0)
cmap_limit=np.nanmax(np.abs( np.nanpercentile(Plot_Var, 99.99)))
#Plot_range=np.linspace(0,cmap_limit,41)
Plot_range=np.linspace(0,0.0055,56)
Plot_unit='(mol carbon / m3)'; Plot_title= 'Large Phytoplankton (phyc minus phypico) at Surface (mol m-3) - '+str(start_date_plt)+'-'+str(end_date_plt)+' - '+str(GCM)
fig, m = func_plotmap_contourf(Plot_Var, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.jet, 'cyl', 210., 80., -80., '-')
fig.savefig(dir_figs+str(GCM)+'_ave_map_PhytoLarge.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

Plot_Var = np.nanmean((Phyc_allmonths - Phypico_allmonths) / Phyc_allmonths ,axis=0) * 100
cmap_limit=np.nanmax(np.abs( np.nanpercentile(Plot_Var, 99.99)))
#Plot_range=np.linspace(0,cmap_limit,41)
Plot_range=np.linspace(0,60,61)
Plot_unit=' % '; Plot_title= 'Large Phytoplankton percentage [(phyc-phypico)/phyc] at Surface (mol m-3) - '+str(start_date_plt)+'-'+str(end_date_plt)+' - '+str(GCM)
fig, m = func_plotmap_contourf(Plot_Var, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.jet, 'cyl', 210., 80., -80., '-')
fig.savefig(dir_figs+str(GCM)+'_ave_map_PhytoLargeprcnt.png', format='png', dpi=300, transparent=True, bbox_inches='tight')











