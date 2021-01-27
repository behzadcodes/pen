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

### Regridding calculations ###
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

#############################################
###   Multi-GCM SingleVariable Plotting  ####
#############################################
#Time_months = ('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec');
#P_Var_x = np.linspace(1,12,12);
#
#P_title= 'Monthly Biomass averages - '+str(start_date_plt)+'-'+str(end_date_plt)
#fig=plt.figure()  
#plt.title('PAPA [45N-50N, 140W-150W]', fontsize=24)
#
#for M_i in range(len(GCM_Names)): # M_i=0     # GCM=GCM_Names[0]
#    
#    GCM=GCM_Names[M_i]
#    print (GCM)     
#    
#    filename_in = (dir_data_in + 'AllResults_'+GCM+'_bio_1991_2010.out') # Directory to save processed data
#    my_shelf = shelve.open(filename_in)
#    if 'Phyc_allmonths' in my_shelf: # If that variable exists in the saved data    
#        globals()['Phyc_allmonths']=my_shelf['Phyc_allmonths']
#        
#        Var_1_allmonths=copy.deepcopy(Phyc_allmonths); 
#        Var_1_monthlymean = Var_1_allmonths.reshape(( np.int(Var_1_allmonths.shape[0]/12) ,12,lat_n_regrid,lon_n_regrid)) # Monthly means of SST over the time period
#        Var_1_monthlymean = np.nanmean(Var_1_monthlymean,axis=0)
#        Var_1_monthlymean=np.squeeze(Var_1_monthlymean)
#        
#        Var_1_monthlymean_PAPA = Var_1_monthlymean[:,135:140,210:220]
#        Var_1_monthlymean_PAPA_timeseries=np.nanmean( np.nanmean(Var_1_monthlymean_PAPA,axis=2) , axis=1); Var_1_monthlymean_PAPA_timeseries=np.squeeze(Var_1_monthlymean_PAPA_timeseries)
#        
#        plt.plot(P_Var_x, Var_1_monthlymean_PAPA_timeseries, c=Colors[M_i], label=GCM, marker='o', markersize=7, markerfacecolor='W', linewidth=2.0)
#        plt.xticks(fontsize = 26); plt.yticks(fontsize = 26)
#        #plt.xlabel('Month', fontsize=26)
#        plt.ylabel('Phyc', fontsize=26)
#        #plt.ylim(275,450)
#    
#plt.xticks(P_Var_x, '', rotation=45) # Setting the X-ticks ## rotation='horizontal' # rotation='vertical'
#plt.legend(prop={'size': 12}, loc='best', fancybox=True, framealpha=0.8)
#plt.xticks(P_Var_x, Time_months, rotation=45, fontsize=24) # Setting the X-ticks ## rotation='horizontal' # rotation='vertical'
#plt.subplots_adjust(left=0.33, bottom=0.11, right=0.69, top=0.9, hspace=0.2, wspace=0.05) # the amount of height/width reserved for space between subplots
#plt.suptitle(P_title, fontsize=26)    
#plt.show()
#mng = plt.get_current_fig_manager()
#mng.window.showMaximized() # Maximizes the plot window to save figures in full    
#
#fig.savefig(dir_figs+str(GCM)+'_PAPA_TimeSeries_phyc_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')    
    
#############################################
###   Multi-GCM MultiVariable Plotting  #####
#############################################
#Time_months = ('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec');
#P_Var_x = np.linspace(1,12,12);
#
#var_names=['Phydiat_allmonths','Phyc_allmonths','Zooc_allmonths','PPint_allmonths','EPC100_allmonths','Chl_allmonths','FgCO2_allmonths','SpCO2_allmonths','Fe_allmonths','NO3_allmonths','Si_allmonths']
#plot_titles=['Diatom (mol/m3)','Biomass (mol/m3)', 'Zooplankton (mol/m3)','Primary Production (mol.m-2.s-1)','Export Production at 100m (mol.m-2.s-1)','Chlorophyll (kg.m-3)','Air-sea CO2 flux (donward +) (gram C /m2/year)','SpCO2 (ppm)', 'Iron - Fe (mol/m3)', 'Nitrate - NO3 (mol/m3)','Silica - Si (mol/m3)']  
#
#for V_i in range(len(var_names)): # V_i=0
#    
#    fig=plt.figure()  
#    P_title= 'PAPA [45N-50N, 140W-150W] - '+str(start_date_plt)+'-'+str(end_date_plt)   
#    plt.title(plot_titles[V_i], fontsize=24)    
#    for M_i in range(len(GCM_Names)): # M_i=0     # GCM=GCM_Names[0]
#        
#        GCM=GCM_Names[M_i]
#        print (GCM)     
#        
#        filename_in = (dir_data_in + 'AllResults_'+GCM+'_bio_1991_2010.out') # Directory to save processed data
#        my_shelf = shelve.open(filename_in)
#        if var_names[V_i] in my_shelf: # If that variable exists in the saved data    
#            globals()['Variable_1']=my_shelf[var_names[V_i]]
#
#            Var_1_allmonths=copy.deepcopy(Variable_1); 
#            Var_1_monthlymean = Var_1_allmonths.reshape(( np.int(Var_1_allmonths.shape[0]/12) ,12,lat_n_regrid,lon_n_regrid)) # Monthly means of SST over the time period
#            Var_1_monthlymean = np.nanmean(Var_1_monthlymean,axis=0)
#            Var_1_monthlymean=np.squeeze(Var_1_monthlymean)
#            
#            Var_1_monthlymean_PAPA = Var_1_monthlymean[:,135:140,210:220]
#            Var_1_monthlymean_PAPA_timeseries=np.nanmean( np.nanmean(Var_1_monthlymean_PAPA,axis=2) , axis=1); Var_1_monthlymean_PAPA_timeseries=np.squeeze(Var_1_monthlymean_PAPA_timeseries)
#            
#            plt.plot(P_Var_x, Var_1_monthlymean_PAPA_timeseries, c=Colors[M_i], label=GCM, marker='o', markersize=7, markerfacecolor='W', linewidth=2.0)
#            plt.xticks(fontsize = 26); plt.yticks(fontsize = 26)
#            #plt.xlabel('Month', fontsize=26)
#            #plt.ylabel('Phyc', fontsize=26)
#            #plt.ylim(275,450)
#        
#    plt.xticks(P_Var_x, '', rotation=45) # Setting the X-ticks ## rotation='horizontal' # rotation='vertical'
#    plt.legend(prop={'size': 12}, loc='best', fancybox=True, framealpha=0.8)
#    plt.xticks(P_Var_x, Time_months, rotation=45, fontsize=24) # Setting the X-ticks ## rotation='horizontal' # rotation='vertical'
#    plt.subplots_adjust(left=0.33, bottom=0.11, right=0.69, top=0.9, hspace=0.2, wspace=0.05) # the amount of height/width reserved for space between subplots
#    plt.suptitle(P_title, fontsize=26)    
#    plt.show()
#    mng = plt.get_current_fig_manager()
#    mng.window.showMaximized() # Maximizes the plot window to save figures in full    
#    
#    fig.savefig(dir_figs+str(GCM)+'_PAPA_'+var_names[V_i]+'_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')     
#    #plt.close()
 
#######################################################
###   Multi-GCM MultiVariable Plotting - Subplots  ####
####################################################### 
GCM_Names = ['CESM1-BGC' ,'GFDL-ESM2G','GFDL-ESM2M','IPSL-CM5A-MR','IPSL-CM5A-LR','MRI-ESM1','GISS-E2-H-CC','GISS-E2-R-CC','HadGEM2-ES','CanESM2','NorESM1-ME','MPI-ESM-MR','MPI-ESM-LR']
Colors    =  [[1, 1, 0],   [1, 0.8, 0],  [1, 0.6, 0], [1, 0.4, 0],   [1, 0, 0], [1, 0.8, 1],   [1, 0.5, 1],   [1, 0, 1],   [0, 0, 0.5], [0, 0, 1], [0, 0.5, 1],[0, 0.5, 0],  [0, 1, 0]]

Time_months = ('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec');
P_Var_x = np.linspace(1,12,12);

plot_ranges=[[-0.001,0.02], [-0.1e-6,5e-6],  [-0.00005, 0.0009],  [-0.0005,0.013], [-0.05e-6,1.4e-6], [-0.1e-7,2.7e-7],   [-20,40],        [270,460],    [-0.1e-6,3.8e-6], [-0.001,0.018], [-0.005,0.08], [0,250]]
var_names=['Phyc_allmonths','Chl_allmonths','Phydiat_allmonths','Zooc_allmonths','PPint_allmonths','EPC100_allmonths','FgCO2_allmonths','SpCO2_allmonths','Fe_allmonths','NO3_allmonths','Si_allmonths','RSNTDS_allmonths']
plot_titles=['Biomass (mol/m3)','Chlorophyll (kg.m-3)','Diatom (mol/m3)', 'Zooplankton (mol/m3)','Primary Production (mol.m-2.s-1)','Export Production at 100m (mol.m-2.s-1)','Air-sea CO2 flux (donward +) (gram C /m2/year)','SpCO2 (ppm)', 'Iron - Fe (mol/m3)', 'Nitrate - NO3 (mol/m3)','Silica - Si (mol/m3)', 'Light - Surface Downward SW Radiation (W/m2)']  

d_trend='yes'

n_r=3 # Number of rows for subplot
n_c=4 # Number of columns for subplot
n_range=list(range(len(var_names)))

fig=plt.figure()
for V_i in range(len(var_names)): # V_i=0
    
    ax = fig.add_subplot(n_r,n_c,V_i+1)   
    for M_i in range(len(GCM_Names)): # M_i=0     # GCM=GCM_Names[0]
        
        GCM=GCM_Names[M_i]
        print (GCM)     
        
        if var_names[V_i]=='RSNTDS_allmonths':
            filename_in = (dir_data_in + 'AllResults_'+GCM+'_physics_1991_2010.out') # Directory to save processed data
        else:
            filename_in = (dir_data_in + 'AllResults_'+GCM+'_bio_1991_2010.out') # Directory to save processed data
        my_shelf = shelve.open(filename_in)
        if var_names[V_i] in my_shelf: # If that variable exists in the saved data    
            globals()['Variable_1']=my_shelf[var_names[V_i]]

            Var_1_allmonths=copy.deepcopy(Variable_1); 
            Var_1_allmonths_PAPA=Var_1_allmonths[:,135:140,210:220]

            if d_trend=='yes': # If it is requested to detrend the data before calculations
                Var_1_allmonths_PAPA= func_detrend_3d(Var_1_allmonths_PAPA)
  
            Var_1_monthlymean_PAPA = Var_1_allmonths_PAPA.reshape(( np.int(Var_1_allmonths_PAPA.shape[0]/12) ,12,Var_1_allmonths_PAPA.shape[1],Var_1_allmonths_PAPA.shape[2])) # Monthly means of SST over the time period
            Var_1_monthlymean_PAPA = np.nanmean(Var_1_monthlymean_PAPA,axis=0)
            Var_1_monthlymean_PAPA=np.squeeze(Var_1_monthlymean_PAPA)           
            Var_1_monthlymean_PAPA_timeseries=np.nanmean( np.nanmean(Var_1_monthlymean_PAPA,axis=2) , axis=1); Var_1_monthlymean_PAPA_timeseries=np.squeeze(Var_1_monthlymean_PAPA_timeseries)
       
            plt.plot(P_Var_x, Var_1_monthlymean_PAPA_timeseries, c=Colors[M_i], label=GCM, marker='o', markersize=7, markerfacecolor='W', linewidth=2.0)
            #plt.xticks(fontsize = 26); plt.yticks(fontsize = 26) 
            plt.ylim(plot_ranges[V_i])

    plt.xticks(P_Var_x, '', rotation=45) # Setting the X-ticks ## rotation='horizontal' # rotation='vertical'
    #plt.legend(prop={'size': 12}, loc='best', fancybox=True, framealpha=0.8)
    plt.xticks(P_Var_x, Time_months, rotation=45, fontsize=10) # Setting the X-ticks ## rotation='horizontal' # rotation='vertical'    
    plt.title(plot_titles[V_i], fontsize=12)

plt.subplots_adjust(left=0.08, bottom=0.07, right=0.87, top=0.9, hspace=0.27, wspace=0.28) # the amount of height/width reserved for space between subplots
#fig.legend(shadow=True, loc='lower right', bbox_to_anchor=(0.96, 0.07), ncol=2 , prop={'size': 13}) # bbox_to_anchor adjusts the location of the legend
fig.legend(shadow=True, loc='lower right', bbox_to_anchor=(0.99, 0.35), ncol=1 , prop={'size': 13}) # bbox_to_anchor adjusts the location of the legend
if d_trend=='yes':
    plt.suptitle('CMIP5 seasonal cylce - North Pacific [45N-50N, 140W-150W] - '+str(start_date_plt)+'-'+str(end_date_plt)+ ' - detrended' , fontsize=18)
else:
    plt.suptitle('CMIP5 seasonal cylce - North Pacific [45N-50N, 140W-150W] - '+str(start_date_plt)+'-'+str(end_date_plt) , fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full

if d_trend=='yes':
    fig.savefig(dir_figs+'AllModels_SesonalCycle_detrended_PAPA_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
else:
    fig.savefig(dir_figs+'AllModels_SesonalCycle_PAPA_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
    

plot_ranges=[[-0.001,0.02], [-0.1e-6,5e-6],  [-0.0005, 0.009],  [-0.0005,0.013], [-0.05e-6,1.4e-6], [-0.1e-7,2.7e-7],   [-20,40],        [270,460],    [-0.1e-6,3.8e-6], [-0.001,0.018], [-0.005,0.08], [0,250]]

fig=plt.figure()
for V_i in range(len(var_names)): # V_i=0
    
    ax = fig.add_subplot(n_r,n_c,V_i+1)   
    for M_i in range(len(GCM_Names)): # M_i=0     # GCM=GCM_Names[0]
        
        GCM=GCM_Names[M_i]
        print (GCM)     
        
        if var_names[V_i]=='RSNTDS_allmonths':
            filename_in = (dir_data_in + 'AllResults_'+GCM+'_physics_1991_2010.out') # Directory to save processed data
        else:
            filename_in = (dir_data_in + 'AllResults_'+GCM+'_bio_1991_2010.out') # Directory to save processed data
        my_shelf = shelve.open(filename_in)
        if var_names[V_i] in my_shelf: # If that variable exists in the saved data    
            globals()['Variable_1']=my_shelf[var_names[V_i]]

            Var_1_allmonths=copy.deepcopy(Variable_1); 
            Var_1_allmonths_NABE=Var_1_allmonths[:,135:140,325:335]

            if d_trend=='yes': # If it is requested to detrend the data before calculations
                Var_1_allmonths_NABE= func_detrend_3d(Var_1_allmonths_NABE)
  
            Var_1_monthlymean_NABE = Var_1_allmonths_NABE.reshape(( np.int(Var_1_allmonths_NABE.shape[0]/12) ,12,Var_1_allmonths_NABE.shape[1],Var_1_allmonths_NABE.shape[2])) # Monthly means of SST over the time period
            Var_1_monthlymean_NABE = np.nanmean(Var_1_monthlymean_NABE,axis=0)
            Var_1_monthlymean_NABE=np.squeeze(Var_1_monthlymean_NABE)           
            Var_1_monthlymean_NABE_timeseries=np.nanmean( np.nanmean(Var_1_monthlymean_NABE,axis=2) , axis=1); Var_1_monthlymean_NABE_timeseries=np.squeeze(Var_1_monthlymean_NABE_timeseries)

            plt.plot(P_Var_x, Var_1_monthlymean_NABE_timeseries, c=Colors[M_i], label=GCM, marker='o', markersize=7, markerfacecolor='W', linewidth=2.0)
            #plt.xticks(fontsize = 26); plt.yticks(fontsize = 26) 
            plt.ylim(plot_ranges[V_i])

    plt.xticks(P_Var_x, '', rotation=45) # Setting the X-ticks ## rotation='horizontal' # rotation='vertical'
    #plt.legend(prop={'size': 12}, loc='best', fancybox=True, framealpha=0.8)
    plt.xticks(P_Var_x, Time_months, rotation=45, fontsize=10) # Setting the X-ticks ## rotation='horizontal' # rotation='vertical'    
    plt.title(plot_titles[V_i], fontsize=12)

plt.subplots_adjust(left=0.08, bottom=0.07, right=0.87, top=0.9, hspace=0.27, wspace=0.28) # the amount of height/width reserved for space between subplots
fig.legend(shadow=True, loc='lower right', bbox_to_anchor=(0.99, 0.35), ncol=1 , prop={'size': 13}) # bbox_to_anchor adjusts the location of the legend
if d_trend=='yes':
    plt.suptitle('CMIP5 seasonal cylce - North Atlantic [45N-50N, 25W-35W] - '+str(start_date_plt)+'-'+str(end_date_plt)+ ' - detrended' , fontsize=18)
else:
    plt.suptitle('CMIP5 seasonal cylce - North Atlantic [45N-50N, 25W-35W] - '+str(start_date_plt)+'-'+str(end_date_plt) , fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full

if d_trend=='yes':
    fig.savefig(dir_figs+'AllModels_SesonalCycle_detrended_NABE_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
else:
    fig.savefig(dir_figs+'AllModels_SesonalCycle_NABE_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


##########################################################
###   SPCO2 Separation Multi-GCM Plotting - Subplots  ####
##########################################################

GCM_Names = ['CESM1-BGC' ,'GFDL-ESM2G','GFDL-ESM2M','MRI-ESM1','GISS-E2-H-CC','GISS-E2-R-CC','HadGEM2-ES','CanESM2','NorESM1-ME','MPI-ESM-MR','MPI-ESM-LR'] # Models that have SpCO2 data
Colors    =  [[1, 1, 0],   [1, 0.8, 0],  [1, 0.6, 0], [1, 0.8, 1],   [1, 0.5, 1],   [1, 0, 1],   [0, 0, 0.5], [0, 0, 1], [0, 0.5, 1],[0, 0.5, 0],  [0, 1, 0]]

Time_months = ('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec');
P_Var_x = np.linspace(1,12,12);

d_trend='yes'

n_r=3 # Number of rows for subplot
n_c=4 # Number of columns for subplot
n_range=list(range(len(GCM_Names)))

####################
### PAPA Plots #####
####################

fig=plt.figure()
for M_i in range(len(GCM_Names)): # M_i=0     # GCM=GCM_Names[0]
    
    ax = fig.add_subplot(n_r,n_c,M_i+1)  

    GCM=GCM_Names[M_i]
    print (GCM)     
    
    filename_in_s = (dir_data_in + 'AllResults_'+GCM+'_bio_1991_2010.out') # Directory to save processed data
    filename_in_t = (dir_data_in + 'AllResults_'+GCM+'_physics_1991_2010.out') # Directory to save processed data
    my_shelf_s = shelve.open(filename_in_s)
    if 'SpCO2_allmonths' in my_shelf_s: # If that variable exists in the saved data    
        globals()['SpCO2_allmonths']=my_shelf_s['SpCO2_allmonths']
            
        my_shelf_t = shelve.open(filename_in_t)
        globals()['SST_allmonths']=my_shelf_t['SST_allmonths']      

        SpCO2_allmonths_PAPA = SpCO2_allmonths[:,135:140,210:220]
        SST_allmonths_PAPA = SST_allmonths[:,135:140,210:220]

        if d_trend=='yes': # If it is requested to detrend the data before calculations
            SpCO2_allmonths_PAPA= func_detrend_3d(SpCO2_allmonths_PAPA)
            SST_allmonths_PAPA= func_detrend_3d(SST_allmonths_PAPA)

#        #### In Takahashi et al formula for SpCO2 separation, average SpCO2 and SST are considerd as average of corresponding year #####
#        SpCO2_anuualave_PAPA=SpCO2_allmonths_PAPA.reshape(( yrs_n ,12,SpCO2_allmonths_PAPA.shape[1],SpCO2_allmonths_PAPA.shape[2])) # Annual means of SpCO2 over the time period
#        SpCO2_anuualave_PAPA = np.nanmean(SpCO2_anuualave_PAPA,axis=1) ; SpCO2_anuualave_PAPA=np.squeeze(SpCO2_anuualave_PAPA)
#        SST_anuualave_PAPA=SST_allmonths_PAPA.reshape(( yrs_n ,12,SST_allmonths_PAPA.shape[1],SST_allmonths_PAPA.shape[2])) # Annual means of SpCO2 over the time period
#        SST_anuualave_PAPA = np.nanmean(SST_anuualave_PAPA,axis=1) ; SST_anuualave_PAPA=np.squeeze(SST_anuualave_PAPA)
#
#        SpCO2_allmonths_TEMP_eq1_PAPA = empty((SpCO2_allmonths_PAPA.shape[0], SpCO2_allmonths_PAPA.shape[1],SpCO2_allmonths_PAPA.shape[2])) * nan #SpCO2_allmonths_TEMP_eq1 = np.zeros([SpCO2_allmonths.shape[0], lat_n_regrid,lon_n_regrid])
#        SpCO2_allmonths_NonTEMP_eq2_PAPA = empty((SpCO2_allmonths_PAPA.shape[0], SpCO2_allmonths_PAPA.shape[1],SpCO2_allmonths_PAPA.shape[2])) * nan #SpCO2_allmonths_NonTEMP_eq2 = np.zeros([SpCO2_allmonths.shape[0], lat_n_regrid,lon_n_regrid])
#        for ii in range(SpCO2_allmonths_PAPA.shape[0]):
#            SpCO2_allmonths_TEMP_eq1_PAPA[ii,:,:] = SpCO2_anuualave_PAPA[ np.int(np.floor(ii/12)) ,:,:] * ( np.exp( 0.0423 * ( SST_allmonths_PAPA[ii,:,:] - SST_anuualave_PAPA[ np.int(np.floor(ii/12)) ,:,:] ) ) );
#            SpCO2_allmonths_NonTEMP_eq2_PAPA[ii,:,:] = SpCO2_allmonths_PAPA[ii,:,:] * ( np.exp( 0.0423 * ( SST_anuualave_PAPA[ np.int(np.floor(ii/12)) ,:,:] - SST_allmonths_PAPA[ii,:,:] ) ) );   

        #### In Takahashi et al formula for SpCO2 separation, average SpCO2 and SST are considerd as average of the whole study period #####
        SpCO2_allmonths_TEMP_eq1_PAPA = empty((SpCO2_allmonths_PAPA.shape[0], SpCO2_allmonths_PAPA.shape[1],SpCO2_allmonths_PAPA.shape[2])) * nan #SpCO2_allmonths_TEMP_eq1 = np.zeros([SpCO2_allmonths.shape[0], lat_n_regrid,lon_n_regrid])
        SpCO2_allmonths_NonTEMP_eq2_PAPA = empty((SpCO2_allmonths_PAPA.shape[0], SpCO2_allmonths_PAPA.shape[1],SpCO2_allmonths_PAPA.shape[2])) * nan #SpCO2_allmonths_NonTEMP_eq2 = np.zeros([SpCO2_allmonths.shape[0], lat_n_regrid,lon_n_regrid])
        for ii in range(SpCO2_allmonths_PAPA.shape[0]):
            SpCO2_allmonths_TEMP_eq1_PAPA[ii,:,:] = np.nanmean(SpCO2_allmonths_PAPA,axis=0) * ( np.exp( 0.0423 * ( SST_allmonths_PAPA[ii,:,:] - np.nanmean(SST_allmonths_PAPA,axis=0) ) ) );
            SpCO2_allmonths_NonTEMP_eq2_PAPA[ii,:,:] = SpCO2_allmonths_PAPA[ii,:,:] * ( np.exp( 0.0423 * ( np.nanmean(SST_allmonths_PAPA,axis=0) - SST_allmonths_PAPA[ii,:,:] ) ) );

        SpCO2_monthlymean_PAPA = SpCO2_allmonths_PAPA.reshape(( yrs_n ,12,SpCO2_allmonths_PAPA.shape[1],SpCO2_allmonths_PAPA.shape[2])) # Monthly means of SpCO2 over the time period
        SpCO2_monthlymean_PAPA = np.nanmean(SpCO2_monthlymean_PAPA,axis=0)
        SpCO2_monthlymean_PAPA=np.squeeze(SpCO2_monthlymean_PAPA)
        SpCO2_monthlymean_PAPA_timeseries=np.nanmean( np.nanmean(SpCO2_monthlymean_PAPA,axis=2) , axis=1); SpCO2_monthlymean_PAPA_timeseries=np.squeeze(SpCO2_monthlymean_PAPA_timeseries)

        SpCO2_monthlymean_TEMP_eq1_PAPA = SpCO2_allmonths_TEMP_eq1_PAPA.reshape(( yrs_n ,12,SpCO2_allmonths_TEMP_eq1_PAPA.shape[1],SpCO2_allmonths_TEMP_eq1_PAPA.shape[2])) # Monthly means of SpCO2 over the time period
        SpCO2_monthlymean_TEMP_eq1_PAPA = np.nanmean(SpCO2_monthlymean_TEMP_eq1_PAPA,axis=0)
        SpCO2_monthlymean_TEMP_eq1_PAPA=np.squeeze(SpCO2_monthlymean_TEMP_eq1_PAPA)
        SpCO2_monthlymean_TEMP_eq1_PAPA_timeseries=np.nanmean( np.nanmean(SpCO2_monthlymean_TEMP_eq1_PAPA,axis=2) , axis=1); SpCO2_monthlymean_TEMP_eq1_PAPA_timeseries=np.squeeze(SpCO2_monthlymean_TEMP_eq1_PAPA_timeseries)

        SpCO2_monthlymean_NonTEMP_eq2_PAPA = SpCO2_allmonths_NonTEMP_eq2_PAPA.reshape(( yrs_n ,12,SpCO2_allmonths_NonTEMP_eq2_PAPA.shape[1],SpCO2_allmonths_NonTEMP_eq2_PAPA.shape[2])) # Monthly means of SpCO2 over the time period
        SpCO2_monthlymean_NonTEMP_eq2_PAPA = np.nanmean(SpCO2_monthlymean_NonTEMP_eq2_PAPA,axis=0)
        SpCO2_monthlymean_NonTEMP_eq2_PAPA=np.squeeze(SpCO2_monthlymean_NonTEMP_eq2_PAPA)
        SpCO2_monthlymean_NonTEMP_eq2_PAPA_timeseries=np.nanmean( np.nanmean(SpCO2_monthlymean_NonTEMP_eq2_PAPA,axis=2) , axis=1); SpCO2_monthlymean_NonTEMP_eq2_PAPA_timeseries=np.squeeze(SpCO2_monthlymean_NonTEMP_eq2_PAPA_timeseries)

        ### Ploting ###
        l = plt.axhline(y=  np.nanmean( SpCO2_monthlymean_PAPA_timeseries, axis=0)    , color='k', linewidth=3.0)
        plt.plot(P_Var_x, SpCO2_monthlymean_PAPA_timeseries, c='r',  label='pCO2', marker='o', markersize=10, markerfacecolor='W', linewidth=4.0)
        plt.plot(P_Var_x, SpCO2_monthlymean_TEMP_eq1_PAPA_timeseries, c='b',  label='pCO2-T', marker='o', markersize=10, markerfacecolor='W', linewidth=4.0)
        plt.plot(P_Var_x, SpCO2_monthlymean_NonTEMP_eq2_PAPA_timeseries, c='g',  label='pCO2-nonT', marker='o', markersize=10, markerfacecolor='W', linewidth=4.0)
        plt.yticks(fontsize = 14)
        #plt.ylabel('pCO2', fontsize=26)
        plt.ylim(260,460)
        plt.xticks(P_Var_x, Time_months, rotation=45, fontsize=10) # Setting the X-ticks ## rotation='horizontal' # rotation='vertical'    
        plt.title(GCM, fontsize=12)

plt.subplots_adjust(left=0.08, bottom=0.06, right=0.96, top=0.9, hspace=0.27, wspace=0.25) # the amount of height/width reserved for space between subplots
fig.legend(shadow=True, loc='lower right', bbox_to_anchor=(0.96, 0.07), ncol=1 , prop={'size': 13}) # bbox_to_anchor adjusts the location of the legend
if d_trend=='yes':
    plt.suptitle('CMIP5 SpCO2 sepration (ppm) - North Pacific [45N-50N,140W-150W] - '+str(start_date_plt)+'-'+str(end_date_plt)+ ' - detrended' , fontsize=18)
else:
    plt.suptitle('CMIP5 SpCO2 sepration (ppm) - North Pacific [45N-50N,140W-150W] - '+str(start_date_plt)+'-'+str(end_date_plt) , fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full

if d_trend=='yes':
    fig.savefig(dir_figs+'AllModels_SpCO2Separation_detrended_PAPA_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
else:
    fig.savefig(dir_figs+'AllModels_SpCO2Separation_PAPA_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
    

fig=plt.figure()
for M_i in range(len(GCM_Names)): # M_i=0     # GCM=GCM_Names[0]
    
    ax = fig.add_subplot(n_r,n_c,M_i+1)  

    GCM=GCM_Names[M_i]
    print (GCM)     
    
    filename_in_s = (dir_data_in + 'AllResults_'+GCM+'_bio_1991_2010.out') # Directory to save processed data
    filename_in_t = (dir_data_in + 'AllResults_'+GCM+'_physics_1991_2010.out') # Directory to save processed data
    my_shelf_s = shelve.open(filename_in_s)
    if 'SpCO2_allmonths' in my_shelf_s: # If that variable exists in the saved data    
            
        my_shelf_t = shelve.open(filename_in_t)
        globals()['SST_allmonths']=my_shelf_t['SST_allmonths']      

        SST_allmonths_PAPA = SST_allmonths[:,135:140,210:220]

        if d_trend=='yes': # If it is requested to detrend the data before calculations
            SST_allmonths_PAPA= func_detrend_3d(SST_allmonths_PAPA)

        SST_monthlymean_PAPA = SST_allmonths_PAPA.reshape(( yrs_n ,12,SST_allmonths_PAPA.shape[1],SST_allmonths_PAPA.shape[2])) # Monthly means of SpCO2 over the time period
        SST_monthlymean_PAPA = np.nanmean(SST_monthlymean_PAPA,axis=0)
        SST_monthlymean_PAPA=np.squeeze(SST_monthlymean_PAPA)
        SST_monthlymean_PAPA_timeseries=np.nanmean( np.nanmean(SST_monthlymean_PAPA,axis=2) , axis=1); SST_monthlymean_PAPA_timeseries=np.squeeze(SST_monthlymean_PAPA_timeseries)

        ### Ploting ###
        plt.plot(P_Var_x, SST_monthlymean_PAPA_timeseries, c='gold',  label='SST', marker='D', markersize=10, markerfacecolor='W', linewidth=4.0)
        plt.yticks(fontsize = 14)
        #plt.ylabel('pCO2', fontsize=26)
        #plt.ylabel('Temp [C]', color='gold', fontsize=18)
        #plt.tick_params(axis='y', labelcolor='gold')
        plt.ylim(4,20)
        plt.xticks(P_Var_x, Time_months, rotation=45, fontsize=10) # Setting the X-ticks ## rotation='horizontal' # rotation='vertical'    
        plt.title(GCM, fontsize=12)

plt.subplots_adjust(left=0.08, bottom=0.07, right=0.95, top=0.9, hspace=0.27, wspace=0.2) # the amount of height/width reserved for space between subplots
fig.legend(shadow=True, loc='lower right', bbox_to_anchor=(0.96, 0.07), ncol=1 , prop={'size': 13}) # bbox_to_anchor adjusts the location of the legend
if d_trend=='yes':
    plt.suptitle('CMIP5 Sea Surface Temperature - North Pacific [45N-50N,140W-150W] - '+str(start_date_plt)+'-'+str(end_date_plt)+ ' - detrended' , fontsize=18)
else:
    plt.suptitle('CMIP5 Sea Surface Temperature - North Pacific [45N-50N,140W-150W] - '+str(start_date_plt)+'-'+str(end_date_plt) , fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full

if d_trend=='yes':
    fig.savefig(dir_figs+'AllModels_SST_detrended_PAPA_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
else:
    fig.savefig(dir_figs+'AllModels_SST_PAPA_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


####################
### NABE Plots #####
####################

fig=plt.figure()
for M_i in range(len(GCM_Names)): # M_i=0     # GCM=GCM_Names[0]
    
    ax = fig.add_subplot(n_r,n_c,M_i+1)  

    GCM=GCM_Names[M_i]
    print (GCM)     
    
    filename_in_s = (dir_data_in + 'AllResults_'+GCM+'_bio_1991_2010.out') # Directory to save processed data
    filename_in_t = (dir_data_in + 'AllResults_'+GCM+'_physics_1991_2010.out') # Directory to save processed data
    my_shelf_s = shelve.open(filename_in_s)
    if 'SpCO2_allmonths' in my_shelf_s: # If that variable exists in the saved data    
        globals()['SpCO2_allmonths']=my_shelf_s['SpCO2_allmonths']
            
        my_shelf_t = shelve.open(filename_in_t)
        globals()['SST_allmonths']=my_shelf_t['SST_allmonths']      

        SpCO2_allmonths_NABE = SpCO2_allmonths[:,135:140,325:335]
        SST_allmonths_NABE = SST_allmonths[:,135:140,325:335]

        if d_trend=='yes': # If it is requested to detrend the data before calculations
            SpCO2_allmonths_NABE= func_detrend_3d(SpCO2_allmonths_NABE)
            SST_allmonths_NABE= func_detrend_3d(SST_allmonths_NABE)
            
        #### In Takahashi et al formula for SpCO2 separation, average SpCO2 and SST are considerd as average of the whole study period #####
        SpCO2_allmonths_TEMP_eq1_NABE = empty((SpCO2_allmonths_NABE.shape[0], SpCO2_allmonths_NABE.shape[1],SpCO2_allmonths_NABE.shape[2])) * nan #SpCO2_allmonths_TEMP_eq1 = np.zeros([SpCO2_allmonths.shape[0], lat_n_regrid,lon_n_regrid])
        SpCO2_allmonths_NonTEMP_eq2_NABE = empty((SpCO2_allmonths_NABE.shape[0], SpCO2_allmonths_NABE.shape[1],SpCO2_allmonths_NABE.shape[2])) * nan #SpCO2_allmonths_NonTEMP_eq2 = np.zeros([SpCO2_allmonths.shape[0], lat_n_regrid,lon_n_regrid])
        for ii in range(SpCO2_allmonths_NABE.shape[0]):
            SpCO2_allmonths_TEMP_eq1_NABE[ii,:,:] = np.nanmean(SpCO2_allmonths_NABE,axis=0) * ( np.exp( 0.0423 * ( SST_allmonths_NABE[ii,:,:] - np.nanmean(SST_allmonths_NABE,axis=0) ) ) );
            SpCO2_allmonths_NonTEMP_eq2_NABE[ii,:,:] = SpCO2_allmonths_NABE[ii,:,:] * ( np.exp( 0.0423 * ( np.nanmean(SST_allmonths_NABE,axis=0) - SST_allmonths_NABE[ii,:,:] ) ) );

        SpCO2_monthlymean_NABE = SpCO2_allmonths_NABE.reshape(( yrs_n ,12,SpCO2_allmonths_NABE.shape[1],SpCO2_allmonths_NABE.shape[2])) # Monthly means of SpCO2 over the time period
        SpCO2_monthlymean_NABE = np.nanmean(SpCO2_monthlymean_NABE,axis=0)
        SpCO2_monthlymean_NABE=np.squeeze(SpCO2_monthlymean_NABE)
        SpCO2_monthlymean_NABE_timeseries=np.nanmean( np.nanmean(SpCO2_monthlymean_NABE,axis=2) , axis=1); SpCO2_monthlymean_NABE_timeseries=np.squeeze(SpCO2_monthlymean_NABE_timeseries)

        SpCO2_monthlymean_TEMP_eq1_NABE = SpCO2_allmonths_TEMP_eq1_NABE.reshape(( yrs_n ,12,SpCO2_allmonths_TEMP_eq1_NABE.shape[1],SpCO2_allmonths_TEMP_eq1_NABE.shape[2])) # Monthly means of SpCO2 over the time period
        SpCO2_monthlymean_TEMP_eq1_NABE = np.nanmean(SpCO2_monthlymean_TEMP_eq1_NABE,axis=0)
        SpCO2_monthlymean_TEMP_eq1_NABE=np.squeeze(SpCO2_monthlymean_TEMP_eq1_NABE)
        SpCO2_monthlymean_TEMP_eq1_NABE_timeseries=np.nanmean( np.nanmean(SpCO2_monthlymean_TEMP_eq1_NABE,axis=2) , axis=1); SpCO2_monthlymean_TEMP_eq1_NABE_timeseries=np.squeeze(SpCO2_monthlymean_TEMP_eq1_NABE_timeseries)

        SpCO2_monthlymean_NonTEMP_eq2_NABE = SpCO2_allmonths_NonTEMP_eq2_NABE.reshape(( yrs_n ,12,SpCO2_allmonths_NonTEMP_eq2_NABE.shape[1],SpCO2_allmonths_NonTEMP_eq2_NABE.shape[2])) # Monthly means of SpCO2 over the time period
        SpCO2_monthlymean_NonTEMP_eq2_NABE = np.nanmean(SpCO2_monthlymean_NonTEMP_eq2_NABE,axis=0)
        SpCO2_monthlymean_NonTEMP_eq2_NABE=np.squeeze(SpCO2_monthlymean_NonTEMP_eq2_NABE)
        SpCO2_monthlymean_NonTEMP_eq2_NABE_timeseries=np.nanmean( np.nanmean(SpCO2_monthlymean_NonTEMP_eq2_NABE,axis=2) , axis=1); SpCO2_monthlymean_NonTEMP_eq2_NABE_timeseries=np.squeeze(SpCO2_monthlymean_NonTEMP_eq2_NABE_timeseries)

        ### Ploting ###
        l = plt.axhline(y=  np.nanmean( SpCO2_monthlymean_NABE_timeseries, axis=0)    , color='k', linewidth=3.0)
        plt.plot(P_Var_x, SpCO2_monthlymean_NABE_timeseries, c='r',  label='pCO2', marker='o', markersize=10, markerfacecolor='W', linewidth=4.0)
        plt.plot(P_Var_x, SpCO2_monthlymean_TEMP_eq1_NABE_timeseries, c='b',  label='pCO2-T', marker='o', markersize=10, markerfacecolor='W', linewidth=4.0)
        plt.plot(P_Var_x, SpCO2_monthlymean_NonTEMP_eq2_NABE_timeseries, c='g',  label='pCO2-nonT', marker='o', markersize=10, markerfacecolor='W', linewidth=4.0)
        plt.yticks(fontsize = 14)
        #plt.ylabel('pCO2', fontsize=26)
        plt.ylim(260,460)
        plt.xticks(P_Var_x, Time_months, rotation=45, fontsize=10) # Setting the X-ticks ## rotation='horizontal' # rotation='vertical'    
        plt.title(GCM, fontsize=12)

plt.subplots_adjust(left=0.08, bottom=0.06, right=0.96, top=0.9, hspace=0.27, wspace=0.25) # the amount of height/width reserved for space between subplots
fig.legend(shadow=True, loc='lower right', bbox_to_anchor=(0.96, 0.07), ncol=1 , prop={'size': 13}) # bbox_to_anchor adjusts the location of the legend
if d_trend=='yes':
    plt.suptitle('CMIP5 SpCO2 sepration (ppm) - North Atlantic [45N-50N, 25W-35W] - '+str(start_date_plt)+'-'+str(end_date_plt)+ ' - detrended' , fontsize=18)
else:
    plt.suptitle('CMIP5 SpCO2 sepration (ppm) - North Atlantic [45N-50N, 25W-35W] - '+str(start_date_plt)+'-'+str(end_date_plt) , fontsize=18)
mng = plt.get_current_fig_manager()
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full

if d_trend=='yes':
    fig.savefig(dir_figs+'AllModels_SpCO2Separation_detrended_NABE_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
else:
    fig.savefig(dir_figs+'AllModels_SpCO2Separation_NABE_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')




fig=plt.figure()
for M_i in range(len(GCM_Names)): # M_i=0     # GCM=GCM_Names[0]
    
    ax = fig.add_subplot(n_r,n_c,M_i+1)  

    GCM=GCM_Names[M_i]
    print (GCM)     
    
    filename_in_s = (dir_data_in + 'AllResults_'+GCM+'_bio_1991_2010.out') # Directory to save processed data
    filename_in_t = (dir_data_in + 'AllResults_'+GCM+'_physics_1991_2010.out') # Directory to save processed data
    my_shelf_s = shelve.open(filename_in_s)
    if 'SpCO2_allmonths' in my_shelf_s: # If that variable exists in the saved data    
            
        my_shelf_t = shelve.open(filename_in_t)
        globals()['SST_allmonths']=my_shelf_t['SST_allmonths']      

        SST_allmonths_NABE = SST_allmonths[:,135:140,325:335]

        if d_trend=='yes': # If it is requested to detrend the data before calculations
            SST_allmonths_NABE= func_detrend_3d(SST_allmonths_NABE)

        SST_monthlymean_NABE = SST_allmonths_NABE.reshape(( yrs_n ,12,SST_allmonths_NABE.shape[1],SST_allmonths_NABE.shape[2])) # Monthly means of SpCO2 over the time period
        SST_monthlymean_NABE = np.nanmean(SST_monthlymean_NABE,axis=0)
        SST_monthlymean_NABE=np.squeeze(SST_monthlymean_NABE)
        SST_monthlymean_NABE_timeseries=np.nanmean( np.nanmean(SST_monthlymean_NABE,axis=2) , axis=1); SST_monthlymean_NABE_timeseries=np.squeeze(SST_monthlymean_NABE_timeseries)

        ### Ploting ###
        plt.plot(P_Var_x, SST_monthlymean_NABE_timeseries, c='gold',  label='SST', marker='D', markersize=10, markerfacecolor='W', linewidth=4.0)
        plt.yticks(fontsize = 14)
        #plt.ylabel('pCO2', fontsize=26)
        #plt.ylabel('Temp [C]', color='gold', fontsize=18)
        #plt.tick_params(axis='y', labelcolor='gold')
        plt.ylim(4,20)
        plt.xticks(P_Var_x, Time_months, rotation=45, fontsize=10) # Setting the X-ticks ## rotation='horizontal' # rotation='vertical'    
        plt.title(GCM, fontsize=12)

plt.subplots_adjust(left=0.08, bottom=0.07, right=0.95, top=0.9, hspace=0.27, wspace=0.2) # the amount of height/width reserved for space between subplots
fig.legend(shadow=True, loc='lower right', bbox_to_anchor=(0.96, 0.07), ncol=1 , prop={'size': 13}) # bbox_to_anchor adjusts the location of the legend
if d_trend=='yes':
    plt.suptitle('CMIP5 Sea Surface Temperature - North Atlantic [45N-50N, 25W-35W] - '+str(start_date_plt)+'-'+str(end_date_plt)+ ' - detrended' , fontsize=18)
else:
    plt.suptitle('CMIP5 Sea Surface Temperature - North Atlantic [45N-50N, 25W-35W] - '+str(start_date_plt)+'-'+str(end_date_plt) , fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full

if d_trend=='yes':
    fig.savefig(dir_figs+'AllModels_SST_detrended_NABE_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
else:
    fig.savefig(dir_figs+'AllModels_SST_NABE_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


#######################################################
###   Multi-GCM MultiVariable Plotting - Subplots  ####
#######################################################   
GCM_Names = ['CESM1-BGC' ,'GFDL-ESM2G','GFDL-ESM2M','IPSL-CM5A-MR','IPSL-CM5A-LR','MRI-ESM1','GISS-E2-H-CC','GISS-E2-R-CC','HadGEM2-ES','CanESM2','NorESM1-ME','MPI-ESM-MR','MPI-ESM-LR']
Colors    =  [[1, 1, 0],   [1, 0.8, 0],  [1, 0.6, 0], [1, 0.4, 0],   [1, 0, 0], [1, 0.8, 1],   [1, 0.5, 1],   [1, 0, 1],   [0, 0, 0.5], [0, 0, 1], [0, 0.5, 1],[0, 0.5, 0],  [0, 1, 0]]

Time_months = ('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec');
P_Var_x = np.linspace(1,12,12);

plot_ranges=[ [250,460], [250,460],  [250,460],  [-0.1e-6,5e-6],      [4,20],     [-0.001,0.02], [-0.1e-6,3.8e-6],[-0.001,0.018],  [-20,1000] ]
var_names=['SpCO2_allmonths',     '',       '','Chl_allmonths', 'SST_allmonths','Phyc_allmonths', 'Fe_allmonths','NO3_allmonths', 'MLD_allmonths']
plot_titles=[ 'pCO2 (ppm)', 'pCO2-T (ppm)', 'pCO2-nonT (ppm)', 'Chlorophyll (kg.m-3)', 'SST (C)', 'Biomass (mol/m3)', 'Iron - Fe (mol/m3)', 'Nitrate - NO3 (mol/m3)', 'MLD (m)'  ]  

d_trend='yes'

n_r=3 # Number of rows for subplot
n_c=3 # Number of columns for subplot
n_range=list(range(len(var_names)))

fig=plt.figure()
for V_i in range(len(var_names)): # V_i=0
    
    ax = fig.add_subplot(n_r,n_c,V_i+1)   

    for M_i in range(len(GCM_Names)): # M_i=0     # GCM=GCM_Names[0]
        
        GCM=GCM_Names[M_i]
        print (GCM)     
        
        if V_i==0 or V_i==3 or V_i==5 or V_i==6 or V_i==7:
            filename_in = (dir_data_in + 'AllResults_'+GCM+'_bio_1991_2010.out') # Directory to save processed data
        elif V_i==4 or V_i==8:
            filename_in = (dir_data_in + 'AllResults_'+GCM+'_physics_1991_2010.out') # Directory to save processed data
        else: # This condition would not happen
            filename_in = (dir_data_in + 'AllResults_'+GCM+'_bio_1991_2010.out') # Directory to save processed data

        if V_i==0 or V_i==3 or V_i==5 or V_i==6 or V_i==7 or V_i==4 or V_i==8:
          
            my_shelf = shelve.open(filename_in)
            if var_names[V_i] in my_shelf: # If that variable exists in the saved data    
                globals()['Variable_1']=my_shelf[var_names[V_i]]
    
                Var_1_allmonths=copy.deepcopy(Variable_1); 
                Var_1_allmonths_PAPA=Var_1_allmonths[:,135:140,210:220]
    
                if d_trend=='yes': # If it is requested to detrend the data before calculations
                    Var_1_allmonths_PAPA= func_detrend_3d(Var_1_allmonths_PAPA)
      
                Var_1_monthlymean_PAPA = Var_1_allmonths_PAPA.reshape(( np.int(Var_1_allmonths_PAPA.shape[0]/12) ,12,Var_1_allmonths_PAPA.shape[1],Var_1_allmonths_PAPA.shape[2])) # Monthly means of SST over the time period
                Var_1_monthlymean_PAPA = np.nanmean(Var_1_monthlymean_PAPA,axis=0)
                Var_1_monthlymean_PAPA=np.squeeze(Var_1_monthlymean_PAPA)           
                Var_1_monthlymean_PAPA_timeseries=np.nanmean( np.nanmean(Var_1_monthlymean_PAPA,axis=2) , axis=1); Var_1_monthlymean_PAPA_timeseries=np.squeeze(Var_1_monthlymean_PAPA_timeseries)
                
                plt.plot(P_Var_x, Var_1_monthlymean_PAPA_timeseries, c=Colors[M_i], label=GCM, marker='o', markersize=7, markerfacecolor='W', linewidth=2.0)
                #plt.xticks(fontsize = 26); plt.yticks(fontsize = 26) 
                plt.ylim(plot_ranges[V_i])

        else:

            filename_in_s = (dir_data_in + 'AllResults_'+GCM+'_bio_1991_2010.out') # Directory to save processed data
            filename_in_t = (dir_data_in + 'AllResults_'+GCM+'_physics_1991_2010.out') # Directory to save processed data
            my_shelf_s = shelve.open(filename_in_s)
            if 'SpCO2_allmonths' in my_shelf_s: # If that variable exists in the saved data    
                globals()['SpCO2_allmonths']=my_shelf_s['SpCO2_allmonths']
                    
                my_shelf_t = shelve.open(filename_in_t)
                globals()['SST_allmonths']=my_shelf_t['SST_allmonths']      
        
                SpCO2_allmonths_PAPA = SpCO2_allmonths[:,135:140,210:220]
                SST_allmonths_PAPA = SST_allmonths[:,135:140,210:220]

                if d_trend=='yes': # If it is requested to detrend the data before calculations
                    SpCO2_allmonths_PAPA= func_detrend_3d(SpCO2_allmonths_PAPA)
                    SST_allmonths_PAPA= func_detrend_3d(SST_allmonths_PAPA)
        
                #### In Takahashi et al formula for SpCO2 separation, average SpCO2 and SST are considerd as average of the whole study period #####
                SpCO2_allmonths_TEMP_eq1_PAPA = empty((SpCO2_allmonths_PAPA.shape[0], SpCO2_allmonths_PAPA.shape[1],SpCO2_allmonths_PAPA.shape[2])) * nan #SpCO2_allmonths_TEMP_eq1 = np.zeros([SpCO2_allmonths.shape[0], lat_n_regrid,lon_n_regrid])
                SpCO2_allmonths_NonTEMP_eq2_PAPA = empty((SpCO2_allmonths_PAPA.shape[0], SpCO2_allmonths_PAPA.shape[1],SpCO2_allmonths_PAPA.shape[2])) * nan #SpCO2_allmonths_NonTEMP_eq2 = np.zeros([SpCO2_allmonths.shape[0], lat_n_regrid,lon_n_regrid])
                for ii in range(SpCO2_allmonths_PAPA.shape[0]):
                    SpCO2_allmonths_TEMP_eq1_PAPA[ii,:,:] = np.nanmean(SpCO2_allmonths_PAPA,axis=0) * ( np.exp( 0.0423 * ( SST_allmonths_PAPA[ii,:,:] - np.nanmean(SST_allmonths_PAPA,axis=0) ) ) );
                    SpCO2_allmonths_NonTEMP_eq2_PAPA[ii,:,:] = SpCO2_allmonths_PAPA[ii,:,:] * ( np.exp( 0.0423 * ( np.nanmean(SST_allmonths_PAPA,axis=0) - SST_allmonths_PAPA[ii,:,:] ) ) );
        
                SpCO2_monthlymean_PAPA = SpCO2_allmonths_PAPA.reshape(( yrs_n ,12,SpCO2_allmonths_PAPA.shape[1],SpCO2_allmonths_PAPA.shape[2])) # Monthly means of SpCO2 over the time period
                SpCO2_monthlymean_PAPA = np.nanmean(SpCO2_monthlymean_PAPA,axis=0)
                SpCO2_monthlymean_PAPA=np.squeeze(SpCO2_monthlymean_PAPA)
                SpCO2_monthlymean_PAPA_timeseries=np.nanmean( np.nanmean(SpCO2_monthlymean_PAPA,axis=2) , axis=1); SpCO2_monthlymean_PAPA_timeseries=np.squeeze(SpCO2_monthlymean_PAPA_timeseries)
        
                SpCO2_monthlymean_TEMP_eq1_PAPA = SpCO2_allmonths_TEMP_eq1_PAPA.reshape(( yrs_n ,12,SpCO2_allmonths_TEMP_eq1_PAPA.shape[1],SpCO2_allmonths_TEMP_eq1_PAPA.shape[2])) # Monthly means of SpCO2 over the time period
                SpCO2_monthlymean_TEMP_eq1_PAPA = np.nanmean(SpCO2_monthlymean_TEMP_eq1_PAPA,axis=0)
                SpCO2_monthlymean_TEMP_eq1_PAPA=np.squeeze(SpCO2_monthlymean_TEMP_eq1_PAPA)
                SpCO2_monthlymean_TEMP_eq1_PAPA_timeseries=np.nanmean( np.nanmean(SpCO2_monthlymean_TEMP_eq1_PAPA,axis=2) , axis=1); SpCO2_monthlymean_TEMP_eq1_PAPA_timeseries=np.squeeze(SpCO2_monthlymean_TEMP_eq1_PAPA_timeseries)
        
                SpCO2_monthlymean_NonTEMP_eq2_PAPA = SpCO2_allmonths_NonTEMP_eq2_PAPA.reshape(( yrs_n ,12,SpCO2_allmonths_NonTEMP_eq2_PAPA.shape[1],SpCO2_allmonths_NonTEMP_eq2_PAPA.shape[2])) # Monthly means of SpCO2 over the time period
                SpCO2_monthlymean_NonTEMP_eq2_PAPA = np.nanmean(SpCO2_monthlymean_NonTEMP_eq2_PAPA,axis=0)
                SpCO2_monthlymean_NonTEMP_eq2_PAPA=np.squeeze(SpCO2_monthlymean_NonTEMP_eq2_PAPA)
                SpCO2_monthlymean_NonTEMP_eq2_PAPA_timeseries=np.nanmean( np.nanmean(SpCO2_monthlymean_NonTEMP_eq2_PAPA,axis=2) , axis=1); SpCO2_monthlymean_NonTEMP_eq2_PAPA_timeseries=np.squeeze(SpCO2_monthlymean_NonTEMP_eq2_PAPA_timeseries)                

                if V_i==1:
                    Var_1_monthlymean_PAPA_timeseries=copy.deepcopy(SpCO2_monthlymean_TEMP_eq1_PAPA_timeseries)
                elif V_i==2:
                    Var_1_monthlymean_PAPA_timeseries=copy.deepcopy(SpCO2_monthlymean_NonTEMP_eq2_PAPA_timeseries)
                else: # This condition would not happen
                    Var_1_monthlymean_PAPA_timeseries=0
            
                plt.plot(P_Var_x, Var_1_monthlymean_PAPA_timeseries, c=Colors[M_i], label=GCM, marker='o', markersize=7, markerfacecolor='W', linewidth=2.0)
                #plt.xticks(fontsize = 26); plt.yticks(fontsize = 26) 
                plt.ylim(plot_ranges[V_i])
            
        plt.xticks(P_Var_x, '', rotation=45) # Setting the X-ticks ## rotation='horizontal' # rotation='vertical'
        #plt.legend(prop={'size': 12}, loc='best', fancybox=True, framealpha=0.8)
        plt.xticks(P_Var_x, Time_months, rotation=45, fontsize=10) # Setting the X-ticks ## rotation='horizontal' # rotation='vertical'    
        plt.title(plot_titles[V_i], fontsize=12)

plt.subplots_adjust(left=0.08, bottom=0.07, right=0.85, top=0.9, hspace=0.27, wspace=0.23) # the amount of height/width reserved for space between subplots
fig.legend(shadow=True, loc='center right', bbox_to_anchor=(0.98, 0.5), ncol=1 , prop={'size': 13}) # bbox_to_anchor adjusts the location of the legend
if d_trend=='yes':
    plt.suptitle('CMIP5 SpCO2 separation and seasonal cylces - North Pacific [45N-50N, 140W-150W] - '+str(start_date_plt)+'-'+str(end_date_plt)+ ' - detrended' , fontsize=18)
else:
    plt.suptitle('CMIP5 SpCO2 separation and seasonal cylces - North Pacific [45N-50N, 140W-150W] - '+str(start_date_plt)+'-'+str(end_date_plt) , fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full

if d_trend=='yes':
    fig.savefig(dir_figs+'AllModels_SpCO2andSesonalCycle_detrended_PAPA_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
else:
    fig.savefig(dir_figs+'AllModels_SpCO2andSesonalCycle_PAPA_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


n_r=3 # Number of rows for subplot
n_c=3 # Number of columns for subplot
n_range=list(range(len(var_names)))

fig=plt.figure()
for V_i in range(len(var_names)): # V_i=0
    
    ax = fig.add_subplot(n_r,n_c,V_i+1)   

    for M_i in range(len(GCM_Names)): # M_i=0     # GCM=GCM_Names[0]
        
        GCM=GCM_Names[M_i]
        print (GCM)     
        
        if V_i==0 or V_i==3 or V_i==5 or V_i==6 or V_i==7:
            filename_in = (dir_data_in + 'AllResults_'+GCM+'_bio_1991_2010.out') # Directory to save processed data
        elif V_i==4 or V_i==8:
            filename_in = (dir_data_in + 'AllResults_'+GCM+'_physics_1991_2010.out') # Directory to save processed data
        else: # This condition would not happen
            filename_in = (dir_data_in + 'AllResults_'+GCM+'_bio_1991_2010.out') # Directory to save processed data

        if V_i==0 or V_i==3 or V_i==5 or V_i==6 or V_i==7 or V_i==4 or V_i==8:
          
            my_shelf = shelve.open(filename_in)
            if var_names[V_i] in my_shelf: # If that variable exists in the saved data    
                globals()['Variable_1']=my_shelf[var_names[V_i]]
    
                Var_1_allmonths=copy.deepcopy(Variable_1); 
                Var_1_allmonths_NABE=Var_1_allmonths[:,135:140,325:335]
    
                if d_trend=='yes': # If it is requested to detrend the data before calculations
                    Var_1_allmonths_NABE= func_detrend_3d(Var_1_allmonths_NABE)
      
                Var_1_monthlymean_NABE = Var_1_allmonths_NABE.reshape(( np.int(Var_1_allmonths_NABE.shape[0]/12) ,12,Var_1_allmonths_NABE.shape[1],Var_1_allmonths_NABE.shape[2])) # Monthly means of SST over the time period
                Var_1_monthlymean_NABE = np.nanmean(Var_1_monthlymean_NABE,axis=0)
                Var_1_monthlymean_NABE=np.squeeze(Var_1_monthlymean_NABE)           
                Var_1_monthlymean_NABE_timeseries=np.nanmean( np.nanmean(Var_1_monthlymean_NABE,axis=2) , axis=1); Var_1_monthlymean_NABE_timeseries=np.squeeze(Var_1_monthlymean_NABE_timeseries)
                
                plt.plot(P_Var_x, Var_1_monthlymean_NABE_timeseries, c=Colors[M_i], label=GCM, marker='o', markersize=7, markerfacecolor='W', linewidth=2.0)
                #plt.xticks(fontsize = 26); plt.yticks(fontsize = 26) 
                plt.ylim(plot_ranges[V_i])

        else:

            filename_in_s = (dir_data_in + 'AllResults_'+GCM+'_bio_1991_2010.out') # Directory to save processed data
            filename_in_t = (dir_data_in + 'AllResults_'+GCM+'_physics_1991_2010.out') # Directory to save processed data
            my_shelf_s = shelve.open(filename_in_s)
            if 'SpCO2_allmonths' in my_shelf_s: # If that variable exists in the saved data    
                globals()['SpCO2_allmonths']=my_shelf_s['SpCO2_allmonths']
                    
                my_shelf_t = shelve.open(filename_in_t)
                globals()['SST_allmonths']=my_shelf_t['SST_allmonths']      
        
                SpCO2_allmonths_NABE = SpCO2_allmonths[:,135:140,325:335]
                SST_allmonths_NABE = SST_allmonths[:,135:140,325:335]

                if d_trend=='yes': # If it is requested to detrend the data before calculations
                    SpCO2_allmonths_NABE= func_detrend_3d(SpCO2_allmonths_NABE)
                    SST_allmonths_NABE= func_detrend_3d(SST_allmonths_NABE)
            
                #### In Takahashi et al formula for SpCO2 separation, average SpCO2 and SST are considerd as average of the whole study period #####
                SpCO2_allmonths_TEMP_eq1_NABE = empty((SpCO2_allmonths_NABE.shape[0], SpCO2_allmonths_NABE.shape[1],SpCO2_allmonths_NABE.shape[2])) * nan #SpCO2_allmonths_TEMP_eq1 = np.zeros([SpCO2_allmonths.shape[0], lat_n_regrid,lon_n_regrid])
                SpCO2_allmonths_NonTEMP_eq2_NABE = empty((SpCO2_allmonths_NABE.shape[0], SpCO2_allmonths_NABE.shape[1],SpCO2_allmonths_NABE.shape[2])) * nan #SpCO2_allmonths_NonTEMP_eq2 = np.zeros([SpCO2_allmonths.shape[0], lat_n_regrid,lon_n_regrid])
                for ii in range(SpCO2_allmonths_NABE.shape[0]):
                    SpCO2_allmonths_TEMP_eq1_NABE[ii,:,:] = np.nanmean(SpCO2_allmonths_NABE,axis=0) * ( np.exp( 0.0423 * ( SST_allmonths_NABE[ii,:,:] - np.nanmean(SST_allmonths_NABE,axis=0) ) ) );
                    SpCO2_allmonths_NonTEMP_eq2_NABE[ii,:,:] = SpCO2_allmonths_NABE[ii,:,:] * ( np.exp( 0.0423 * ( np.nanmean(SST_allmonths_NABE,axis=0) - SST_allmonths_NABE[ii,:,:] ) ) );
        
                SpCO2_monthlymean_NABE = SpCO2_allmonths_NABE.reshape(( yrs_n ,12,SpCO2_allmonths_NABE.shape[1],SpCO2_allmonths_NABE.shape[2])) # Monthly means of SpCO2 over the time period
                SpCO2_monthlymean_NABE = np.nanmean(SpCO2_monthlymean_NABE,axis=0)
                SpCO2_monthlymean_NABE=np.squeeze(SpCO2_monthlymean_NABE)
                SpCO2_monthlymean_NABE_timeseries=np.nanmean( np.nanmean(SpCO2_monthlymean_NABE,axis=2) , axis=1); SpCO2_monthlymean_NABE_timeseries=np.squeeze(SpCO2_monthlymean_NABE_timeseries)
        
                SpCO2_monthlymean_TEMP_eq1_NABE = SpCO2_allmonths_TEMP_eq1_NABE.reshape(( yrs_n ,12,SpCO2_allmonths_TEMP_eq1_NABE.shape[1],SpCO2_allmonths_TEMP_eq1_NABE.shape[2])) # Monthly means of SpCO2 over the time period
                SpCO2_monthlymean_TEMP_eq1_NABE = np.nanmean(SpCO2_monthlymean_TEMP_eq1_NABE,axis=0)
                SpCO2_monthlymean_TEMP_eq1_NABE=np.squeeze(SpCO2_monthlymean_TEMP_eq1_NABE)
                SpCO2_monthlymean_TEMP_eq1_NABE_timeseries=np.nanmean( np.nanmean(SpCO2_monthlymean_TEMP_eq1_NABE,axis=2) , axis=1); SpCO2_monthlymean_TEMP_eq1_NABE_timeseries=np.squeeze(SpCO2_monthlymean_TEMP_eq1_NABE_timeseries)
        
                SpCO2_monthlymean_NonTEMP_eq2_NABE = SpCO2_allmonths_NonTEMP_eq2_NABE.reshape(( yrs_n ,12,SpCO2_allmonths_NonTEMP_eq2_NABE.shape[1],SpCO2_allmonths_NonTEMP_eq2_NABE.shape[2])) # Monthly means of SpCO2 over the time period
                SpCO2_monthlymean_NonTEMP_eq2_NABE = np.nanmean(SpCO2_monthlymean_NonTEMP_eq2_NABE,axis=0)
                SpCO2_monthlymean_NonTEMP_eq2_NABE=np.squeeze(SpCO2_monthlymean_NonTEMP_eq2_NABE)
                SpCO2_monthlymean_NonTEMP_eq2_NABE_timeseries=np.nanmean( np.nanmean(SpCO2_monthlymean_NonTEMP_eq2_NABE,axis=2) , axis=1); SpCO2_monthlymean_NonTEMP_eq2_NABE_timeseries=np.squeeze(SpCO2_monthlymean_NonTEMP_eq2_NABE_timeseries)                

                if V_i==1:
                    Var_1_monthlymean_NABE_timeseries=copy.deepcopy(SpCO2_monthlymean_TEMP_eq1_NABE_timeseries)
                elif V_i==2:
                    Var_1_monthlymean_NABE_timeseries=copy.deepcopy(SpCO2_monthlymean_NonTEMP_eq2_NABE_timeseries)
                else: # This condition would not happen
                    Var_1_monthlymean_NABE_timeseries=0
            
                plt.plot(P_Var_x, Var_1_monthlymean_NABE_timeseries, c=Colors[M_i], label=GCM, marker='o', markersize=7, markerfacecolor='W', linewidth=2.0)
                #plt.xticks(fontsize = 26); plt.yticks(fontsize = 26) 
                plt.ylim(plot_ranges[V_i])
            
        plt.xticks(P_Var_x, '', rotation=45) # Setting the X-ticks ## rotation='horizontal' # rotation='vertical'
        #plt.legend(prop={'size': 12}, loc='best', fancybox=True, framealpha=0.8)
        plt.xticks(P_Var_x, Time_months, rotation=45, fontsize=10) # Setting the X-ticks ## rotation='horizontal' # rotation='vertical'    
        plt.title(plot_titles[V_i], fontsize=12)

plt.subplots_adjust(left=0.08, bottom=0.07, right=0.85, top=0.9, hspace=0.27, wspace=0.23) # the amount of height/width reserved for space between subplots
fig.legend(shadow=True, loc='center right', bbox_to_anchor=(0.98, 0.5), ncol=1 , prop={'size': 13}) # bbox_to_anchor adjusts the location of the legend
if d_trend=='yes':
    plt.suptitle('CMIP5 SpCO2 separation and seasonal cylces - North Atlantic [45N-50N, 25W-35W] - '+str(start_date_plt)+'-'+str(end_date_plt)+ ' - detrended' , fontsize=18)
else:
    plt.suptitle('CMIP5 SpCO2 separation and seasonal cylces - North Atlantic [45N-50N, 25W-35W] - '+str(start_date_plt)+'-'+str(end_date_plt) , fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full

if d_trend=='yes':
    fig.savefig(dir_figs+'AllModels_SpCO2andSesonalCycle_detrended_NABE_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
else:
    fig.savefig(dir_figs+'AllModels_SpCO2andSesonalCycle_NABE_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


#####################################################
###   Correlation Multi-GCM Plotting - Subplots  ####
#####################################################
GCM_Names = ['CESM1-BGC' ,'GFDL-ESM2G','GFDL-ESM2M','IPSL-CM5A-MR','IPSL-CM5A-LR','MRI-ESM1','GISS-E2-H-CC','GISS-E2-R-CC','HadGEM2-ES','CanESM2','NorESM1-ME','MPI-ESM-MR','MPI-ESM-LR']
Colors    =  [[1, 1, 0],   [1, 0.8, 0],  [1, 0.6, 0], [1, 0.4, 0],   [1, 0, 0], [1, 0.8, 1],   [1, 0.5, 1],   [1, 0, 1],   [0, 0, 0.5], [0, 0, 1], [0, 0.5, 1],[0, 0.5, 0],  [0, 1, 0]]

#### Plots 1 #####
var_names_y=[ 'MLD_allmonths', 'MLD_allmonths', 'MLD_allmonths',  'Phyc_allmonths',      'SpCO2_allmonths',     'SpCO2_allmonths', 'SpCO2_T',                'SpCO2_nonT',                'SpCO2_T',           'Phyc_allmonths',      'Chl_allmonths',     'Phyc_allmonths',            'Phyc_allmonths',                'Phyc_allmonths',               'Phydiat_allmonths'       ]
var_names_x=[ 'Chl_allmonths', 'SST_allmonths', 'SpCO2_allmonths', 'MLD_allmonths',      'Phyc_allmonths',      'SST_allmonths',   'Phyc_allmonths',         'Phyc_allmonths',            'SST_allmonths',     'SST_allmonths',       'SpCO2_allmonths',  'Fe_allmonths',               'NO3_allmonths',                 'Si_allmonths',                 'Si_allmonths'         ]
plot_titles=[ 'MLD vs. log chl', 'MLD vs. SST',  'MLD vs. pCO2',  'log biomass vs. MLD', 'pCO2 vs. log biomass', 'pCO2 vs. SST',   'pCO2-T vs. log biomass', 'pCO2-nonT vs. log biomass', 'pCO2-T vs. SST',    'log biomass vs. SST', 'log chl vs. pCO2', 'log biomass vs. Iron (Fe)',  'log biomass vs. Nitrate (NO3)', 'log biomass vs. Silica (Si)',  'log diatoms vs. Silica (Si)'    ]
save_names= [ 'MLD_vs_Chl',      'MLD_vs_SST',   'MLD_vs_SpCO2',  'Phyc_vs_MLD',         'SpCO2_vs_Phyc',        'SpCO2_vs_SST',    'SpCO2Temp_vs_Phyc',     'SpCO2nponT_vs_Phyc',        'SpCO2Temp_vs_SST', 'Phyc_vs_SST',        'Chl_vs_SpCO2',       'Phyc_vs_Fe',                  'Phyc_vs_NO3',                   'Phyc_vs_Si',                  'Phydiat_vs_Si'       ]  
        
file_names_y=['physics',        'physics',       'physics',        'bio',                 'bio',                   'bio',            'bio',                    'bio',                       'bio',              'bio',                 'bio',               'bio',                        'bio',                          'bio',                          'bio'       ]
file_names_x=[ 'bio',           'physics',       'bio',            'physics',             'bio',                   'physics',        'bio',                    'bio',                       'physics',          'physics',             'bio',               'bio',                        'bio',                          'bio',                          'bio'      ]

#### Plots 2 #####
var_names_y=[ 'SpCO2_allmonths', 'SpCO2_allmonths',   'FgCO2_allmonths',          'FgCO2_allmonths',    'FgCO2_allmonths',       'SpCO2_nonT',          'SpCO2_nonT',             'PPint_allmonths',     'PPint_allmonths',         'PPint_allmonths',       'PPint_allmonths',  'EPC100_allmonths',       'EPC100_allmonths',           'EPC100_allmonths',         'EPC100_allmonths'                  ]
var_names_x=[ 'PPint_allmonths', 'EPC100_allmonths',  'Phyc_allmonths',           'PPint_allmonths',    'EPC100_allmonths',      'PPint_allmonths',     'EPC100_allmonths',       'Fe_allmonths',        'NO3_allmonths',           'Si_allmonths',          'RSNTDS_allmonths', 'Fe_allmonths',           'NO3_allmonths',              'Si_allmonths',             'RSNTDS_allmonths'                    ]
plot_titles=[ 'pCO2 vs. intPP',  'pCO2 vs. EPC 100m', 'Co2 flux vs. log biomass', 'Co2 flux vs. intPP', 'Co2 flux vs. EPC 100m', 'pCO2-nonT vs. intPP', 'pCO2-nonT vs. EPC 100m', 'intPP vs. Iron (Fe)', 'intPP vs. Nitrate (NO3)', 'intPP vs. Silica (Si)', 'intPP vs. Light',  'EPC 100m vs. Iron (Fe)', 'EPC 100m vs. Nitrate (NO3)', 'EPC 100m vs. Silica (Si)', 'EPC 100m vs. Light'                          ]
save_names= [ 'SpCO2_vs_PPint',  'SpCO2_vs_EPC100',   'FgCO2_vs_Phyc',            'FgCO2_vs_PPint',     'FgCO2_vs_EPC100',       'SpCO2nponT_vs_PPint', 'SpCO2nponT_vs_EPC100',   'PPint_vs_Fe',         'PPint_vs_NO3',            'PPint_vs_Si',           'PPint_vs_RSNTDS',  'EPC100_vs_Fe',           'EPC100_vs_NO3',              'EPC100_vs_Si',             'EPC100_vs_RSNTDS'        ]  
        
file_names_y=['bio',             'bio',               'bio',                      'bio',                 'bio',                   'bio',                 'bio',                    'bio',                 'bio',                     'bio',                   'bio',              'bio',                   'bio',                         'bio',                     'bio'          ]
file_names_x=['bio',             'bio',               'bio',                      'bio',                 'bio',                   'bio',                 'bio',                    'bio',                 'bio',                     'bio',                   'physics',          'bio',                   'bio',                         'bio',                     'physics'    ]

#### Plots 3 ##### 
var_names_y=[ 'Zooc_allmonths',        'Zooc_allmonths',            'Zooc_allmonths',          'Zooc_allmonths',    'Phyc_allmonths',            'Phyc_allmonths',                'Phyc_allmonths',              'Phyc_allmonths',        'Zoo_by_Phyto',            'Zoo_by_Phyto',                'Zoo_by_Phyto',              'Zoo_by_Phyto',        'Phyc_allmonths'                 ]
var_names_x=[ 'Fe_allmonths',          'NO3_allmonths',             'Si_allmonths',            'RSNTDS_allmonths',  'Fe_allmonths',              'NO3_allmonths',                 'Si_allmonths',                'RSNTDS_allmonths',      'Fe_allmonths',            'NO3_allmonths',               'Si_allmonths',              'RSNTDS_allmonths',    'Zooc_allmonths'                ]
plot_titles=[ 'log Zoo vs. Iron (Fe)', 'log Zoo vs. Nitrate (NO3)', 'log Zoo vs. Silica (Si)', 'log Zoo vs. Light', 'log biomass vs. Iron (Fe)', 'log biomass vs. Nitrate (NO3)', 'log biomass vs. Silica (Si)', 'log biomass vs. Light', 'Zoo/phyto vs. Iron (Fe)', 'Zoo/phyto vs. Nitrate (NO3)', 'Zoo/phyto vs. Silica (Si)', 'Zoo/phyto vs. Light', 'log biomass vs. log Zoo'   ]
save_names= [ 'Zooc_vs_Fe',            'Zooc_vs_NO3',               'Zooc_vs_Si',              'Zooc_vs_RSNTDS',    'Phyc_vs_Fe',                'Phyc_vs_NO3',                   'Phyc_vs_Si',                  'Phyc_vs_RSNTDS',        'ZoobyPhyc_vs_Fe',         'ZoobyPhyc_vs_NO3',            'ZoobyPhyc_vs_Si',           'ZoobyPhyc_vs_RSNTDS', 'Phyc_vs_Zoo'                      ]  
        
file_names_y=[ 'bio',                   'bio',                       'bio',                     'bio',               'bio',                        'bio',                          'bio',                         'bio',                  'bio',                     'bio',                         'bio',                       'bio',                  'bio'                       ]
file_names_x=[ 'bio',                   'bio',                       'bio',                     'physics',           'bio',                        'bio',                          'bio',                         'physics',              'bio',                     'bio',                         'bio',                       'physics',              'bio'                      ]
###################
 
P_cmap=plt.cm.bwr; P_c_fill='fill'; P_Lon=Lon_regrid_2D; P_Lat=Lat_regrid_2D; P_range=np.linspace(-1.,1.,51);
P_proj='cyl'; P_lon0=210.; P_latN=90.; P_latS=-90.; 

d_trend='NO'

n_r=4 ; n_c=4 ;
n_t=len(GCM_Names)

for V_i in range(len(var_names_y)): # V_i=0   V_i=1    V_i=2    V_i=3    V_i=4    V_i=5    V_i=6    V_i=7    V_i=8    V_i=9    V_i=10    V_i=11    V_i=12    V_i=13    V_i=14    V_i=15
    
    if d_trend=='yes':
        P_title='Correlations of '+plot_titles[V_i]+' - CMIP5 models - '+str(start_date_plt)+'-'+str(end_date_plt)+' - detrended'
    else:
        P_title='Correlations of '+plot_titles[V_i]+' - CMIP5 models - '+str(start_date_plt)+'-'+str(end_date_plt)

    fig=plt.figure()
    for M_i in range(len(GCM_Names)): # M_i=0     # GCM=GCM_Names[0]
        
        ax = fig.add_subplot(n_r,n_c,M_i+1)  
    
        GCM=GCM_Names[M_i]
        print (GCM) 
   
        if var_names_y[V_i]=='SpCO2_T' or var_names_y[V_i]=='SpCO2_nonT' or var_names_x[V_i]=='SpCO2_T' or var_names_x[V_i]=='SpCO2_nonT': # Only for SpCO2 sepration plots
                
            filename_in_spco2 = (dir_data_in + 'AllResults_'+GCM+'_'+'bio'+'_1991_2010.out') # Directory to save processed data
            filename_in_sst = (dir_data_in + 'AllResults_'+GCM+'_'+'physics'+'_1991_2010.out') # Directory to save processed data

            if var_names_y[V_i]=='SpCO2_T' or var_names_y[V_i]=='SpCO2_nonT':
                filename_in_xy = (dir_data_in + 'AllResults_'+GCM+'_'+file_names_x[V_i]+'_1991_2010.out') # Directory to save processed data
                var_names_xy=var_names_x[V_i]
            elif var_names_x[V_i]=='SpCO2_T' or var_names_x[V_i]=='SpCO2_nonT':
                filename_in_xy = (dir_data_in + 'AllResults_'+GCM+'_'+file_names_y[V_i]+'_1991_2010.out') # Directory to save processed data
                var_names_xy=var_names_y[V_i]

            my_shelf_spco2 = shelve.open(filename_in_spco2)
            my_shelf_xy = shelve.open(filename_in_xy)
            my_shelf_sst = shelve.open(filename_in_sst)
            
            if 'SpCO2_allmonths' in my_shelf_spco2  and 'SST_allmonths' in my_shelf_sst  and var_names_xy in my_shelf_xy: # If that variable exists in the saved data  
                
                globals()['Var_spco2']=my_shelf_spco2['SpCO2_allmonths']
                globals()['Var_sst']=my_shelf_sst['SST_allmonths']
                globals()['Variable_xy']=my_shelf_xy[var_names_xy]                     
                               
                Variable_xy[ np.abs(Variable_xy)>1e19 ]=nan; Var_spco2[ np.abs(Var_spco2)>1e19 ]=nan; Var_sst[ np.abs(Var_sst)>1e19 ]=nan  # Replacing missing values with NaN 

                if d_trend=='yes': # If it is requested to detrend the data before correlations
                    Variable_xy= func_detrend_3d(Variable_xy)
                    Var_spco2= func_detrend_3d(Var_spco2)
                    Var_sst= func_detrend_3d(Var_sst)

                if var_names_xy=='Chl_allmonths' or var_names_xy=='Phyc_allmonths' or var_names_xy=='Zooc_allmonths' or var_names_xy=='Phydiat_allmonths' or var_names_xy=='Phypico_allmonths': # For these variables, correlations are done on log10 of the variable
                    Variable_xy=np.log10(Variable_xy) 

                
                Var_spco2_TEMP_eq1 = empty((Var_spco2.shape[0], Var_spco2.shape[1],Var_spco2.shape[2])) * nan #SpCO2_allmonths_TEMP_eq1 = np.zeros([SpCO2_allmonths.shape[0], lat_n_regrid,lon_n_regrid])
                Var_spco2_NonTEMP_eq2 = empty((Var_spco2.shape[0], Var_spco2.shape[1],Var_spco2.shape[2])) * nan #SpCO2_allmonths_NonTEMP_eq2 = np.zeros([SpCO2_allmonths.shape[0], lat_n_regrid,lon_n_regrid])
                for ii in range(Var_spco2.shape[0]):
                    Var_spco2_TEMP_eq1[ii,:,:] = np.nanmean(Var_spco2,axis=0) * ( np.exp( 0.0423 * ( Var_sst[ii,:,:] - np.nanmean(Var_sst,axis=0) ) ) );
                    Var_spco2_NonTEMP_eq2[ii,:,:] = Var_spco2[ii,:,:] * ( np.exp( 0.0423 * ( np.nanmean(Var_sst,axis=0) - Var_sst[ii,:,:] ) ) );
                
                if var_names_y[V_i]=='SpCO2_T' or var_names_x[V_i]=='SpCO2_T': # Only for SpCO2 sepration plots                   
                    Corr_xy, Pvalue_xy = func_corr_3d(Var_spco2_TEMP_eq1, Variable_xy);  Corr_xy[Ocean_Index==10]=nan; Pvalue_xy[Ocean_Index==10]=nan;   #  Excluding the Arctic Ocean
                if var_names_y[V_i]=='SpCO2_nonT' or var_names_x[V_i]=='SpCO2_nonT': # Only for SpCO2 sepration plots       
                    Corr_xy, Pvalue_xy = func_corr_3d(Var_spco2_NonTEMP_eq2, Variable_xy);  Corr_xy[Ocean_Index==10]=nan; Pvalue_xy[Ocean_Index==10]=nan;   #  Excluding the Arctic Ocean
  
                Var_plot_ii=copy.deepcopy(Corr_xy)                
                color_land='k'
                
            else: # If we don't have both datasets for this model to do the correlation, just plot an empty map, to keep the order of GCMs constant in all models                
                Var_plot_ii=empty((lat_n_regrid,lon_n_regrid)) * nan
                color_land='w'


        elif var_names_y[V_i]=='Zoo_by_Phyto': # Only for SpCO2 sepration plots
            
            filename_in_y = (dir_data_in + 'AllResults_'+GCM+'_'+'bio'+'_1991_2010.out') # Directory to save processed data
            filename_in_x = (dir_data_in + 'AllResults_'+GCM+'_'+file_names_x[V_i]+'_1991_2010.out') # Directory to save processed data
            my_shelf_y = shelve.open(filename_in_y)
            my_shelf_x = shelve.open(filename_in_x)
     
            if 'Zooc_allmonths' in my_shelf_y and 'Phyc_allmonths' in my_shelf_y and var_names_x[V_i] in my_shelf_x: # If that variable exists in the saved data  
                
                globals()['Variable_1']=my_shelf_y['Zooc_allmonths'] 
                globals()['Variable_2']=my_shelf_y['Phyc_allmonths']
                globals()['Variable_x']=my_shelf_x[var_names_x[V_i]]        

                Variable_1[ np.abs(Variable_1)>1e19 ]=nan # Replacing missing values with NaN
                Variable_2[ np.abs(Variable_2)>1e19 ]=nan # Replacing missing values with NaN
                    
                if d_trend=='yes': # If it is requested to detrend the data before correlations
                    Variable_1= func_detrend_3d(Variable_1)
                    Variable_2= func_detrend_3d(Variable_2)
    
                Variable_y=Variable_1 / Variable_2

                Variable_y[ np.abs(Variable_y)>1e19 ]=nan; Variable_x[ np.abs(Variable_x)>1e19 ]=nan  # Replacing missing values with NaN 
                    
                if d_trend=='yes': # If it is requested to detrend the data before correlations
                    Variable_x= func_detrend_3d(Variable_x)            

                if var_names_x[V_i]=='Chl_allmonths' or var_names_x[V_i]=='Phyc_allmonths' or var_names_x[V_i]=='Zooc_allmonths' or var_names_x=='Phydiat_allmonths' or var_names_x=='Phypico_allmonths':
                    Variable_x=np.log10(Variable_x) 
                    
                Corr_xy, Pvalue_xy = func_corr_3d(Variable_y, Variable_x);  Corr_xy[Ocean_Index==10]=nan; Pvalue_xy[Ocean_Index==10]=nan;   #  Excluding the Arctic Ocean
    
                Var_plot_ii=copy.deepcopy(Corr_xy)
                
                color_land='k'
            else: # If we don't have both datasets for this model to do the correlation, just plot an empty map, to keep the order of GCMs constant in all models
                
                Var_plot_ii=empty((lat_n_regrid,lon_n_regrid)) * nan
                color_land='w'
                    
        else:

            filename_in_y = (dir_data_in + 'AllResults_'+GCM+'_'+file_names_y[V_i]+'_1991_2010.out') # Directory to save processed data
            filename_in_x = (dir_data_in + 'AllResults_'+GCM+'_'+file_names_x[V_i]+'_1991_2010.out') # Directory to save processed data
            my_shelf_y = shelve.open(filename_in_y)
            my_shelf_x = shelve.open(filename_in_x)
     
            if var_names_y[V_i] in my_shelf_y and var_names_x[V_i] in my_shelf_x: # If that variable exists in the saved data  
                
                globals()['Variable_y']=my_shelf_y[var_names_y[V_i]]
                globals()['Variable_x']=my_shelf_x[var_names_x[V_i]]        

                Variable_y[ np.abs(Variable_y)>1e19 ]=nan; Variable_x[ np.abs(Variable_x)>1e19 ]=nan  # Replacing missing values with NaN 
                    
                if d_trend=='yes': # If it is requested to detrend the data before correlations
                    Variable_y= func_detrend_3d(Variable_y)
                    Variable_x= func_detrend_3d(Variable_x)

                if var_names_y[V_i]=='Chl_allmonths' or var_names_y[V_i]=='Phyc_allmonths' or var_names_y[V_i]=='Zooc_allmonths' or var_names_y=='Phydiat_allmonths' or var_names_y=='Phypico_allmonths': # For these variables, correlations are done on log10 of the variable
                    Variable_y=np.log10(Variable_y)
                if var_names_x[V_i]=='Chl_allmonths' or var_names_x[V_i]=='Phyc_allmonths' or var_names_x[V_i]=='Zooc_allmonths' or var_names_x=='Phydiat_allmonths' or var_names_x=='Phypico_allmonths':
                    Variable_x=np.log10(Variable_x)                    
                    
                Corr_xy, Pvalue_xy = func_corr_3d(Variable_y, Variable_x);  Corr_xy[Ocean_Index==10]=nan; Pvalue_xy[Ocean_Index==10]=nan;   #  Excluding the Arctic Ocean
    
                Var_plot_ii=copy.deepcopy(Corr_xy)
                
                color_land='k'
            else: # If we don't have both datasets for this model to do the correlation, just plot an empty map, to keep the order of GCMs constant in all models
                
                Var_plot_ii=empty((lat_n_regrid,lon_n_regrid)) * nan
                color_land='w'
            
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
        
        if P_c_fill=='fill':
            m.fillcontinents(color='0')
        if color_land=='w':
            m.fillcontinents(color='1')
        
        m.drawcoastlines(linewidth=1.0, linestyle='solid', color=color_land, antialiased=1, ax=None, zorder=None)
        im=m.contourf(P_Lon, P_Lat, Var_plot_ii,P_range,latlon=True, cmap=P_cmap) #, extend='both')
        plt.title(GCM, fontsize=14)
            
    plt.subplots_adjust(left=0.05, bottom=0.05, right=0.93, top=0.92, hspace=0.3, wspace=0.02) # the amount of height/width reserved for space between subplots
    cbar_ax = fig.add_axes([0.94, 0.05, 0.015, 0.87]) # [right,bottom,width,height] 
    fig.colorbar(im, cax=cbar_ax, ticks=[-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1])
    cbar_ax.set_yticklabels(['- 1', '- 0.8', '- 0.6', '- 0.4', '- 0.2', '0', '0.2', '0.4', '0.6', '0.8', '1'])
    plt.suptitle(P_title, fontsize=20)    
    plt.show()
    mng = plt.get_current_fig_manager()
    mng.window.showMaximized() # Maximizes the plot window to save figures in full   

    if d_trend=='yes':
        fig.savefig(dir_figs+'AllModels_Correlations_Maps_detrended_'+save_names[V_i]+'_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
    else:
        fig.savefig(dir_figs+'AllModels_Correlations_Maps_'+save_names[V_i]+'_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


d_trend='yes'
n_r=4 ; n_c=4 ;
n_t=len(var_names_y)

fig=plt.figure()
for V_i in range(n_t): # V_i=0   V_i=1    V_i=2    V_i=3    V_i=4    V_i=5    V_i=6    V_i=7    V_i=8    V_i=9    V_i=10    V_i=11    V_i=12    V_i=13

    ax = fig.add_subplot(n_r,n_c,V_i+1) 

    for M_i in range(len(GCM_Names)): # M_i=0     # GCM=GCM_Names[0]

        GCM=GCM_Names[M_i]
        print (GCM) 
   
        if var_names_y[V_i]=='SpCO2_T' or var_names_y[V_i]=='SpCO2_nonT' or var_names_x[V_i]=='SpCO2_T' or var_names_x[V_i]=='SpCO2_nonT': # Only for SpCO2 sepration plots
                
            filename_in_spco2 = (dir_data_in + 'AllResults_'+GCM+'_'+'bio'+'_1991_2010.out') # Directory to save processed data
            filename_in_sst = (dir_data_in + 'AllResults_'+GCM+'_'+'physics'+'_1991_2010.out') # Directory to save processed data

            if var_names_y[V_i]=='SpCO2_T' or var_names_y[V_i]=='SpCO2_nonT':
                filename_in_xy = (dir_data_in + 'AllResults_'+GCM+'_'+file_names_x[V_i]+'_1991_2010.out') # Directory to save processed data
                var_names_xy=var_names_x[V_i]
            elif var_names_x[V_i]=='SpCO2_T' or var_names_x[V_i]=='SpCO2_nonT':
                filename_in_xy = (dir_data_in + 'AllResults_'+GCM+'_'+file_names_y[V_i]+'_1991_2010.out') # Directory to save processed data
                var_names_xy=var_names_y[V_i]

            my_shelf_spco2 = shelve.open(filename_in_spco2)
            my_shelf_xy = shelve.open(filename_in_xy)
            my_shelf_sst = shelve.open(filename_in_sst)
            
            if 'SpCO2_allmonths' in my_shelf_spco2  and 'SST_allmonths' in my_shelf_sst  and var_names_xy in my_shelf_xy: # If that variable exists in the saved data  
                
                globals()['Var_spco2']=my_shelf_spco2['SpCO2_allmonths']
                globals()['Var_sst']=my_shelf_sst['SST_allmonths']
                globals()['Variable_xy']=my_shelf_xy[var_names_xy]                     

                Variable_xy[ np.abs(Variable_xy)>1e19 ]=nan; Var_spco2[ np.abs(Var_spco2)>1e19 ]=nan; Var_sst[ np.abs(Var_sst)>1e19 ]=nan  # Replacing missing values with NaN 

                Var_spco2_PAPA = Var_spco2[:,135:140,210:220]
                Var_sst_PAPA = Var_sst[:,135:140,210:220]
                Variable_xy_PAPA = Variable_xy[:,135:140,210:220]
                
                Var_spco2_NABE = Var_spco2[:,135:140,325:335]
                Var_sst_NABE = Var_sst[:,135:140,325:335]
                Variable_xy_NABE = Variable_xy[:,135:140,325:335]               

                if d_trend=='yes': # If it is requested to detrend the data before correlations
                    Variable_xy_PAPA= func_detrend_3d(Variable_xy_PAPA)
                    Var_spco2_PAPA= func_detrend_3d(Var_spco2_PAPA)
                    Var_sst_PAPA= func_detrend_3d(Var_sst_PAPA)
                    Variable_xy_NABE= func_detrend_3d(Variable_xy_NABE)
                    Var_spco2_NABE= func_detrend_3d(Var_spco2_NABE)
                    Var_sst_NABE= func_detrend_3d(Var_sst_NABE)
 
                if var_names_xy=='Chl_allmonths' or var_names_xy=='Phyc_allmonths' or var_names_xy=='Zooc_allmonths' or var_names_xy=='Phydiat_allmonths' or var_names_xy=='Phypico_allmonths': # For these variables, correlations are done on log10 of the variable
                    Variable_xy_PAPA=np.log10(Variable_xy_PAPA)  
                    Variable_xy_NABE=np.log10(Variable_xy_NABE) 
              
                Var_spco2_TEMP_eq1_PAPA = empty((Var_spco2_PAPA.shape[0], Var_spco2_PAPA.shape[1],Var_spco2_PAPA.shape[2])) * nan #SpCO2_allmonths_TEMP_eq1 = np.zeros([SpCO2_allmonths.shape[0], lat_n_regrid,lon_n_regrid])
                Var_spco2_NonTEMP_eq2_PAPA = empty((Var_spco2_PAPA.shape[0], Var_spco2_PAPA.shape[1],Var_spco2_PAPA.shape[2])) * nan #SpCO2_allmonths_NonTEMP_eq2 = np.zeros([SpCO2_allmonths.shape[0], lat_n_regrid,lon_n_regrid])
                for ii in range(Var_spco2_PAPA.shape[0]):
                    Var_spco2_TEMP_eq1_PAPA[ii,:,:] = np.nanmean(Var_spco2_PAPA,axis=0) * ( np.exp( 0.0423 * ( Var_sst_PAPA[ii,:,:] - np.nanmean(Var_sst_PAPA,axis=0) ) ) );
                    Var_spco2_NonTEMP_eq2_PAPA[ii,:,:] = Var_spco2_PAPA[ii,:,:] * ( np.exp( 0.0423 * ( np.nanmean(Var_sst_PAPA,axis=0) - Var_sst_PAPA[ii,:,:] ) ) );

                Var_spco2_TEMP_eq1_NABE = empty((Var_spco2_NABE.shape[0], Var_spco2_NABE.shape[1],Var_spco2_NABE.shape[2])) * nan #SpCO2_allmonths_TEMP_eq1 = np.zeros([SpCO2_allmonths.shape[0], lat_n_regrid,lon_n_regrid])
                Var_spco2_NonTEMP_eq2_NABE = empty((Var_spco2_NABE.shape[0], Var_spco2_NABE.shape[1],Var_spco2_NABE.shape[2])) * nan #SpCO2_allmonths_NonTEMP_eq2 = np.zeros([SpCO2_allmonths.shape[0], lat_n_regrid,lon_n_regrid])
                for ii in range(Var_spco2_NABE.shape[0]):
                    Var_spco2_TEMP_eq1_NABE[ii,:,:] = np.nanmean(Var_spco2_NABE,axis=0) * ( np.exp( 0.0423 * ( Var_sst_NABE[ii,:,:] - np.nanmean(Var_sst_NABE,axis=0) ) ) );
                    Var_spco2_NonTEMP_eq2_NABE[ii,:,:] = Var_spco2_NABE[ii,:,:] * ( np.exp( 0.0423 * ( np.nanmean(Var_sst_NABE,axis=0) - Var_sst_NABE[ii,:,:] ) ) );

                
                if var_names_y[V_i]=='SpCO2_T' or var_names_x[V_i]=='SpCO2_T': # Only for SpCO2 sepration plots                                       
                    Corr_xy_PAPA, Pvalue_xy_PAPA = func_corr_3d(Var_spco2_TEMP_eq1_PAPA, Variable_xy_PAPA);
                    Corr_xy_NABE, Pvalue_xy_NABE = func_corr_3d(Var_spco2_TEMP_eq1_NABE, Variable_xy_NABE);                                       
                    
                if var_names_y[V_i]=='SpCO2_nonT' or var_names_x[V_i]=='SpCO2_nonT': # Only for SpCO2 sepration plots       

                    Corr_xy_PAPA, Pvalue_xy_PAPA = func_corr_3d(Var_spco2_NonTEMP_eq2_PAPA, Variable_xy_PAPA);
                    Corr_xy_NABE, Pvalue_xy_NABE = func_corr_3d(Var_spco2_NonTEMP_eq2_NABE, Variable_xy_NABE);  
                    
                Corr_xy_PAPA=np.nanmean(np.nanmean(Corr_xy_PAPA,axis=0) ,axis=0); Pvalue_xy_PAPA=np.nanmean(np.nanmean(Pvalue_xy_PAPA,axis=0) ,axis=0)
                Corr_xy_NABE=np.nanmean(np.nanmean(Corr_xy_NABE,axis=0) ,axis=0); Pvalue_xy_NABE=np.nanmean(np.nanmean(Pvalue_xy_NABE,axis=0) ,axis=0)
                
            else: # If we don't have both datasets for this model to do the correlation, just plot an empty map, to keep the order of GCMs constant in all models                
                Corr_xy_PAPA=nan; Pvalue_xy_PAPA=nan; Corr_xy_NABE=nan; Pvalue_xy_NABE=nan;

        elif var_names_y[V_i]=='Zoo_by_Phyto': # Only for SpCO2 sepration plots
            
            filename_in_y = (dir_data_in + 'AllResults_'+GCM+'_'+'bio'+'_1991_2010.out') # Directory to save processed data
            filename_in_x = (dir_data_in + 'AllResults_'+GCM+'_'+file_names_x[V_i]+'_1991_2010.out') # Directory to save processed data
            my_shelf_y = shelve.open(filename_in_y)
            my_shelf_x = shelve.open(filename_in_x)
     
            if 'Zooc_allmonths' in my_shelf_y and 'Phyc_allmonths' in my_shelf_y and var_names_x[V_i] in my_shelf_x: # If that variable exists in the saved data  
                
                globals()['Variable_1']=my_shelf_y['Zooc_allmonths'] 
                globals()['Variable_2']=my_shelf_y['Phyc_allmonths']
                globals()['Variable_x']=my_shelf_x[var_names_x[V_i]]        

                Variable_1[ np.abs(Variable_1)>1e19 ]=nan # Replacing missing values with NaN
                Variable_2[ np.abs(Variable_2)>1e19 ]=nan # Replacing missing values with NaN
                Variable_y=Variable_1 / Variable_2
                
                Variable_y[ np.abs(Variable_y)>1e19 ]=nan; Variable_x[ np.abs(Variable_x)>1e19 ]=nan  # Replacing missing values with NaN 

                Variable_y_PAPA = Variable_y[:,135:140,210:220]
                Variable_x_PAPA = Variable_x[:,135:140,210:220]                
                Variable_y_NABE = Variable_y[:,135:140,325:335]
                Variable_x_NABE = Variable_x[:,135:140,325:335]

                if d_trend=='yes': # If it is requested to detrend the data before correlations
                    Variable_y_PAPA= func_detrend_3d(Variable_y_PAPA)
                    Variable_x_PAPA= func_detrend_3d(Variable_x_PAPA)
                    Variable_y_NABE= func_detrend_3d(Variable_y_NABE)
                    Variable_x_NABE= func_detrend_3d(Variable_x_NABE)

                if var_names_y[V_i]=='Chl_allmonths' or var_names_y[V_i]=='Phyc_allmonths' or var_names_y[V_i]=='Zooc_allmonths' or var_names_y=='Phydiat_allmonths' or var_names_y=='Phypico_allmonths': # For these variables, correlations are done on log10 of the variable
                    Variable_y_PAPA=np.log10(Variable_y_PAPA)
                    Variable_y_NABE=np.log10(Variable_y_NABE)
                if var_names_x[V_i]=='Chl_allmonths' or var_names_x[V_i]=='Phyc_allmonths' or var_names_x[V_i]=='Zooc_allmonths' or var_names_x=='Phydiat_allmonths' or var_names_x=='Phypico_allmonths':
                    Variable_x_PAPA=np.log10(Variable_x_PAPA) 
                    Variable_x_NABE=np.log10(Variable_x_NABE) 
                     
                Corr_xy_PAPA, Pvalue_xy_PAPA = func_corr_3d(Variable_y_PAPA, Variable_x_PAPA);
                Corr_xy_NABE, Pvalue_xy_NABE = func_corr_3d(Variable_y_NABE, Variable_x_NABE);
                Corr_xy_PAPA=np.nanmean(np.nanmean(Corr_xy_PAPA,axis=0) ,axis=0); Pvalue_xy_PAPA=np.nanmean(np.nanmean(Pvalue_xy_PAPA,axis=0) ,axis=0)
                Corr_xy_NABE=np.nanmean(np.nanmean(Corr_xy_NABE,axis=0) ,axis=0); Pvalue_xy_NABE=np.nanmean(np.nanmean(Pvalue_xy_NABE,axis=0) ,axis=0)
                
            else: # If we don't have both datasets for this model to do the correlation, just plot an empty map, to keep the order of GCMs constant in all models
                
                Corr_xy_PAPA=nan; Pvalue_xy_PAPA=nan; Corr_xy_NABE=nan; Pvalue_xy_NABE=nan;

        else:

            filename_in_y = (dir_data_in + 'AllResults_'+GCM+'_'+file_names_y[V_i]+'_1991_2010.out') # Directory to save processed data
            filename_in_x = (dir_data_in + 'AllResults_'+GCM+'_'+file_names_x[V_i]+'_1991_2010.out') # Directory to save processed data
            my_shelf_y = shelve.open(filename_in_y)
            my_shelf_x = shelve.open(filename_in_x)
     
            if var_names_y[V_i] in my_shelf_y and var_names_x[V_i] in my_shelf_x: # If that variable exists in the saved data  
                
                globals()['Variable_y']=my_shelf_y[var_names_y[V_i]]
                globals()['Variable_x']=my_shelf_x[var_names_x[V_i]]        
                
                Variable_y[ np.abs(Variable_y)>1e19 ]=nan; Variable_x[ np.abs(Variable_x)>1e19 ]=nan  # Replacing missing values with NaN 
                
                Variable_y_PAPA = Variable_y[:,135:140,210:220]
                Variable_x_PAPA = Variable_x[:,135:140,210:220]                
                Variable_y_NABE = Variable_y[:,135:140,325:335]
                Variable_x_NABE = Variable_x[:,135:140,325:335]

                if d_trend=='yes': # If it is requested to detrend the data before correlations
                    Variable_y_PAPA= func_detrend_3d(Variable_y_PAPA)
                    Variable_x_PAPA= func_detrend_3d(Variable_x_PAPA)
                    Variable_y_NABE= func_detrend_3d(Variable_y_NABE)
                    Variable_x_NABE= func_detrend_3d(Variable_x_NABE)

                if var_names_y[V_i]=='Chl_allmonths' or var_names_y[V_i]=='Phyc_allmonths' or var_names_y[V_i]=='Zooc_allmonths': # For these variables, correlations are done on log10 of the variable
                    Variable_y_PAPA=np.log10(Variable_y_PAPA)
                    Variable_y_NABE=np.log10(Variable_y_NABE)
                if var_names_x[V_i]=='Chl_allmonths' or var_names_x[V_i]=='Phyc_allmonths' or var_names_x[V_i]=='Zooc_allmonths':
                    Variable_x_PAPA=np.log10(Variable_x_PAPA) 
                    Variable_x_NABE=np.log10(Variable_x_NABE) 
                     
                Corr_xy_PAPA, Pvalue_xy_PAPA = func_corr_3d(Variable_y_PAPA, Variable_x_PAPA);
                Corr_xy_NABE, Pvalue_xy_NABE = func_corr_3d(Variable_y_NABE, Variable_x_NABE);
                Corr_xy_PAPA=np.nanmean(np.nanmean(Corr_xy_PAPA,axis=0) ,axis=0); Pvalue_xy_PAPA=np.nanmean(np.nanmean(Pvalue_xy_PAPA,axis=0) ,axis=0)
                Corr_xy_NABE=np.nanmean(np.nanmean(Corr_xy_NABE,axis=0) ,axis=0); Pvalue_xy_NABE=np.nanmean(np.nanmean(Pvalue_xy_NABE,axis=0) ,axis=0)
                
            else: # If we don't have both datasets for this model to do the correlation, just plot an empty map, to keep the order of GCMs constant in all models
                
                Corr_xy_PAPA=nan; Pvalue_xy_PAPA=nan; Corr_xy_NABE=nan; Pvalue_xy_NABE=nan;
                
        plt.scatter(Corr_xy_PAPA, Corr_xy_NABE, s=100, marker='o', c=Colors[M_i], label=GCM) # Only the first pannel gets the lables so the legend location be alligned with that
        plt.scatter(Corr_xy_PAPA, Corr_xy_NABE, s=100, marker='o', facecolors='none', edgecolors='k')
        
    l = plt.axhline(y= 0 , color='k', linewidth=1.0); l = plt.axvline(x= 0 , color='k', linewidth=1.0)
    plt.plot([-2,2],[-2,2], 'k--', linewidth=0.75) ; plt.plot([-2,2],[2,-2], 'k--', linewidth=0.75)
    plt.xlim(-1.1,1.1) ; plt.ylim(-1.1,1.1)
    plt.xticks(fontsize=10); plt.yticks(fontsize = 10) # Setting the X-ticks ## rotation='horizontal' # rotation='vertical'   
    plt.title(plot_titles[V_i], fontsize=14)

plt.subplots_adjust(left=0.18, bottom=0.06, right=0.83, top=0.92, hspace=0.35, wspace=0.25) # the amount of height/width reserved for space between subplots
#plt.legend(shadow=True, loc='right', bbox_to_anchor=(3.3, 2.7), markerscale=1 , prop={'size': 14})
#fig.legend(shadow=True, loc='lower right', bbox_to_anchor=(0.96, 0.07), ncol=2 , prop={'size': 13}) # bbox_to_anchor adjusts the location of the legend
if d_trend=='yes':
    plt.suptitle('CMIP5 models - Average correlations of North Pacific (x-axis) and North Atlantic (y-axis) - '+str(start_date_plt)+'-'+str(end_date_plt)+' - detrended' , fontsize=18)
else:
    plt.suptitle('CMIP5 models - Average correlations of North Pacific (x-axis) and North Atlantic (y-axis) - '+str(start_date_plt)+'-'+str(end_date_plt) , fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full

if d_trend=='yes':
    fig.savefig(dir_figs+'AllModels_Correlations_Scatter_detrended_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
else:
    fig.savefig(dir_figs+'AllModels_Correlations_Scatter_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')



##########################################################
###   Seasonality scatter plots Multi-GCM - Subplots  ####
##########################################################
GCM_Names = ['CESM1-BGC' ,'GFDL-ESM2G','GFDL-ESM2M','IPSL-CM5A-MR','IPSL-CM5A-LR','MRI-ESM1','GISS-E2-H-CC','GISS-E2-R-CC','HadGEM2-ES','CanESM2','NorESM1-ME','MPI-ESM-MR','MPI-ESM-LR']
Colors    =  [[1, 1, 0],   [1, 0.8, 0],  [1, 0.6, 0], [1, 0.4, 0],   [1, 0, 0], [1, 0.8, 1],   [1, 0.5, 1],   [1, 0, 1],   [0, 0, 0.5], [0, 0, 1], [0, 0.5, 1],[0, 0.5, 0],  [0, 1, 0]]

plot_ranges=[[-0.01,9],       [-0.01,9],        [-0.01,9],          [-0.01,9],       [-0.01,4.7],   [-0.01,4.7],     [-0.01,1.3],     [-1,9],     [-1,25],         [-0.01,0.5],     [-0.01,0.5], [-0.01,0.5], [-0.01,2.3], [-0.01,4], [-0.01,2.3], [-0.01,2.3]]
var_names=['Phyc_allmonths','Chl_allmonths','Phydiat_allmonths','Zooc_allmonths','PPint_allmonths','EPC100_allmonths','SST_allmonths','MLD_allmonths','FgCO2_allmonths','SpCO2_allmonths','SpCO2_T','SpCO2_nonT','Fe_allmonths','NO3_allmonths','Si_allmonths','RSNTDS_allmonths']
plot_titles=['Biomass','Chlorophyll','Diatom', 'Zooplankton','Primary Production','Export Production at 100m',   'SST',   'MLD',  'Air-sea CO2 flux', 'SpCO2', 'SpCO2-T', 'SpCO2-nonT', 'Iron - Fe', 'Nitrate - NO3','Silica - Si', 'Light']  
file_names=[ 'bio',     'bio',      'bio',      'bio',        'bio',              'bio',                      'physics','physics', 'bio',              'bio',   'bio',      'bio',       'bio',      'bio',          'bio',       'physics'      ]

d_trend='yes'
n_r=4 ; n_c=4 ;
n_t=len(var_names)

fig=plt.figure()
for V_i in range(n_t): # V_i=0   V_i=1    V_i=2    V_i=3    V_i=4    V_i=5    V_i=6    V_i=7    V_i=8    V_i=9    V_i=10    V_i=11    V_i=12    V_i=13

    ax = fig.add_subplot(n_r,n_c,V_i+1) 

    for M_i in range(len(GCM_Names)): # M_i=0     # GCM=GCM_Names[0]

        GCM=GCM_Names[M_i]
        print (GCM) 
   
        if var_names[V_i]=='SpCO2_T' or var_names[V_i]=='SpCO2_nonT': # Only for SpCO2 sepration plots
                
            filename_in_spco2 = (dir_data_in + 'AllResults_'+GCM+'_'+'bio'+'_1991_2010.out') # Directory to save processed data
            filename_in_sst = (dir_data_in + 'AllResults_'+GCM+'_'+'physics'+'_1991_2010.out') # Directory to save processed data
            
            my_shelf_spco2 = shelve.open(filename_in_spco2)
            my_shelf_sst = shelve.open(filename_in_sst)
            
            if 'SpCO2_allmonths' in my_shelf_spco2  and 'SST_allmonths' in my_shelf_sst: # If that variable exists in the saved data  
                
                globals()['Var_spco2']=my_shelf_spco2['SpCO2_allmonths']
                globals()['Var_sst']=my_shelf_sst['SST_allmonths']                

                Var_spco2[ np.abs(Var_spco2)>1e19 ]=nan; Var_sst[ np.abs(Var_sst)>1e19 ]=nan  # Replacing missing values with NaN 

                Var_spco2_PAPA = Var_spco2[:,135:140,210:220]
                Var_sst_PAPA = Var_sst[:,135:140,210:220]                
                Var_spco2_NABE = Var_spco2[:,135:140,325:335]
                Var_sst_NABE = Var_sst[:,135:140,325:335]              

                if d_trend=='yes': # If it is requested to detrend the data before correlations
                    Var_spco2_PAPA= func_detrend_3d(Var_spco2_PAPA)
                    Var_sst_PAPA= func_detrend_3d(Var_sst_PAPA)
                    Var_spco2_NABE= func_detrend_3d(Var_spco2_NABE)
                    Var_sst_NABE= func_detrend_3d(Var_sst_NABE)
              
                Var_spco2_TEMP_eq1_PAPA = empty((Var_spco2_PAPA.shape[0], Var_spco2_PAPA.shape[1],Var_spco2_PAPA.shape[2])) * nan #SpCO2_allmonths_TEMP_eq1 = np.zeros([SpCO2_allmonths.shape[0], lat_n_regrid,lon_n_regrid])
                Var_spco2_NonTEMP_eq2_PAPA = empty((Var_spco2_PAPA.shape[0], Var_spco2_PAPA.shape[1],Var_spco2_PAPA.shape[2])) * nan #SpCO2_allmonths_NonTEMP_eq2 = np.zeros([SpCO2_allmonths.shape[0], lat_n_regrid,lon_n_regrid])
                for ii in range(Var_spco2_PAPA.shape[0]):
                    Var_spco2_TEMP_eq1_PAPA[ii,:,:] = np.nanmean(Var_spco2_PAPA,axis=0) * ( np.exp( 0.0423 * ( Var_sst_PAPA[ii,:,:] - np.nanmean(Var_sst_PAPA,axis=0) ) ) );
                    Var_spco2_NonTEMP_eq2_PAPA[ii,:,:] = Var_spco2_PAPA[ii,:,:] * ( np.exp( 0.0423 * ( np.nanmean(Var_sst_PAPA,axis=0) - Var_sst_PAPA[ii,:,:] ) ) );

                Var_spco2_TEMP_eq1_NABE = empty((Var_spco2_NABE.shape[0], Var_spco2_NABE.shape[1],Var_spco2_NABE.shape[2])) * nan #SpCO2_allmonths_TEMP_eq1 = np.zeros([SpCO2_allmonths.shape[0], lat_n_regrid,lon_n_regrid])
                Var_spco2_NonTEMP_eq2_NABE = empty((Var_spco2_NABE.shape[0], Var_spco2_NABE.shape[1],Var_spco2_NABE.shape[2])) * nan #SpCO2_allmonths_NonTEMP_eq2 = np.zeros([SpCO2_allmonths.shape[0], lat_n_regrid,lon_n_regrid])
                for ii in range(Var_spco2_NABE.shape[0]):
                    Var_spco2_TEMP_eq1_NABE[ii,:,:] = np.nanmean(Var_spco2_NABE,axis=0) * ( np.exp( 0.0423 * ( Var_sst_NABE[ii,:,:] - np.nanmean(Var_sst_NABE,axis=0) ) ) );
                    Var_spco2_NonTEMP_eq2_NABE[ii,:,:] = Var_spco2_NABE[ii,:,:] * ( np.exp( 0.0423 * ( np.nanmean(Var_sst_NABE,axis=0) - Var_sst_NABE[ii,:,:] ) ) );
                
                if var_names[V_i]=='SpCO2_T': # Only for SpCO2 sepration plots                                       
                    Variable_1_PAPA=copy.deepcopy(Var_spco2_TEMP_eq1_PAPA)
                    Variable_1_NABE=copy.deepcopy(Var_spco2_TEMP_eq1_NABE)                                      
                    
                if var_names[V_i]=='SpCO2_nonT': # Only for SpCO2 sepration plots       
                    Variable_1_PAPA=copy.deepcopy(Var_spco2_NonTEMP_eq2_PAPA)
                    Variable_1_NABE=copy.deepcopy(Var_spco2_NonTEMP_eq2_NABE) 
                    
                Variable_1_PAPA_monthlymean = Variable_1_PAPA.reshape(( np.int(Variable_1_PAPA.shape[0]/12) ,12,Variable_1_PAPA.shape[1],Variable_1_PAPA.shape[2])) # Monthly means of SST over the time period
                Variable_1_PAPA_monthlymean = np.nanmean(Variable_1_PAPA_monthlymean,axis=0)            
                Var_plot_PAPA= ( np.nanmax(Variable_1_PAPA_monthlymean, axis=0) - np.nanmin(Variable_1_PAPA_monthlymean, axis=0) ) / np.nanmean(Variable_1_PAPA_monthlymean, axis=0)
                Var_plot_PAPA=np.nanmean(np.nanmean(Var_plot_PAPA,axis=0) ,axis=0)

                Variable_1_NABE_monthlymean = Variable_1_NABE.reshape(( np.int(Variable_1_NABE.shape[0]/12) ,12,Variable_1_NABE.shape[1],Variable_1_NABE.shape[2])) # Monthly means of SST over the time period
                Variable_1_NABE_monthlymean = np.nanmean(Variable_1_NABE_monthlymean,axis=0)            
                Var_plot_NABE= ( np.nanmax(Variable_1_NABE_monthlymean, axis=0) - np.nanmin(Variable_1_NABE_monthlymean, axis=0) ) / np.nanmean(Variable_1_NABE_monthlymean, axis=0)
                Var_plot_NABE=np.nanmean(np.nanmean(Var_plot_NABE,axis=0) ,axis=0)
                
            else: # If we don't have both datasets for this model to do the correlation, just plot an empty map, to keep the order of GCMs constant in all models                
                Var_plot_PAPA=nan; Var_plot_NABE; 

        else:

            filename_in = (dir_data_in + 'AllResults_'+GCM+'_'+file_names[V_i]+'_1991_2010.out') # Directory to save processed data
            my_shelf = shelve.open(filename_in)
     
            if var_names[V_i] in my_shelf: # If that variable exists in the saved data  
                
                globals()['Variable_1']=my_shelf[var_names[V_i]]       
                
                Variable_1[ np.abs(Variable_1)>1e19 ]=nan; # Replacing missing values with NaN 
                
                Variable_1_PAPA = Variable_1[:,135:140,210:220]           
                Variable_1_NABE = Variable_1[:,135:140,325:335]

                if d_trend=='yes': # If it is requested to detrend the data before correlations
                    Variable_1_PAPA= func_detrend_3d(Variable_1_PAPA)
                    Variable_1_NABE= func_detrend_3d(Variable_1_NABE)

                Variable_1_PAPA_monthlymean = Variable_1_PAPA.reshape(( np.int(Variable_1_PAPA.shape[0]/12) ,12,Variable_1_PAPA.shape[1],Variable_1_PAPA.shape[2])) # Monthly means of SST over the time period
                Variable_1_PAPA_monthlymean = np.nanmean(Variable_1_PAPA_monthlymean,axis=0)            
                Var_plot_PAPA= ( np.nanmax(Variable_1_PAPA_monthlymean, axis=0) - np.nanmin(Variable_1_PAPA_monthlymean, axis=0) ) / np.nanmean(Variable_1_PAPA_monthlymean, axis=0)
                Var_plot_PAPA=np.nanmean(np.nanmean(Var_plot_PAPA,axis=0) ,axis=0)

                Variable_1_NABE_monthlymean = Variable_1_NABE.reshape(( np.int(Variable_1_NABE.shape[0]/12) ,12,Variable_1_NABE.shape[1],Variable_1_NABE.shape[2])) # Monthly means of SST over the time period
                Variable_1_NABE_monthlymean = np.nanmean(Variable_1_NABE_monthlymean,axis=0)            
                Var_plot_NABE= ( np.nanmax(Variable_1_NABE_monthlymean, axis=0) - np.nanmin(Variable_1_NABE_monthlymean, axis=0) ) / np.nanmean(Variable_1_NABE_monthlymean, axis=0)
                Var_plot_NABE=np.nanmean(np.nanmean(Var_plot_NABE,axis=0) ,axis=0)
                
            else: # If we don't have both datasets for this model to do the correlation, just plot an empty map, to keep the order of GCMs constant in all models
                
                Var_plot_PAPA=nan; Var_plot_NABE; 
                
        plt.scatter(Var_plot_PAPA, Var_plot_NABE, s=100, marker='o', c=Colors[M_i], label=GCM) # Only the first pannel gets the lables so the legend location be alligned with that
        plt.scatter(Var_plot_PAPA, Var_plot_NABE, s=100, marker='o', facecolors='none', edgecolors='k')
        
    plt.plot([-1000,1000],[-1000,1000], 'k--', linewidth=0.75)
    plt.xlim(plot_ranges[V_i][0],plot_ranges[V_i][1]) ; plt.ylim(plot_ranges[V_i][0],plot_ranges[V_i][1])
    plt.xticks(fontsize=10); plt.yticks(fontsize = 10) # Setting the X-ticks ## rotation='horizontal' # rotation='vertical'   
    plt.title(plot_titles[V_i], fontsize=14)

plt.subplots_adjust(left=0.18, bottom=0.06, right=0.83, top=0.92, hspace=0.35, wspace=0.25) # the amount of height/width reserved for space between subplots
plt.legend(shadow=True, loc='right', bbox_to_anchor=(2.1, 2.5), markerscale=1 , prop={'size': 14})
if d_trend=='yes':
    plt.suptitle('CMIP5 models - Seasonality [(max-min)/ave] of North Pacific (x-axis) and North Atlantic (y-axis) - '+str(start_date_plt)+'-'+str(end_date_plt)+' - detrended' , fontsize=18)
else:
    plt.suptitle('CMIP5 models - Seasonality [(max-min)/ave] of North Pacific (x-axis) and North Atlantic (y-axis) - '+str(start_date_plt)+'-'+str(end_date_plt) , fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full

if d_trend=='yes':
    fig.savefig(dir_figs+'AllModels_Seasonality_Scatter_detrended_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
else:
    fig.savefig(dir_figs+'AllModels_Seasonality_Scatter_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


############################################################
###   Average value scatter plots Multi-GCM - Subplots  ####
############################################################
GCM_Names = ['CESM1-BGC' ,'GFDL-ESM2G','GFDL-ESM2M','IPSL-CM5A-MR','IPSL-CM5A-LR','MRI-ESM1','GISS-E2-H-CC','GISS-E2-R-CC','HadGEM2-ES','CanESM2','NorESM1-ME','MPI-ESM-MR','MPI-ESM-LR']
Colors    =  [[1, 1, 0],   [1, 0.8, 0],  [1, 0.6, 0], [1, 0.4, 0],   [1, 0, 0], [1, 0.8, 1],   [1, 0.5, 1],   [1, 0, 1],   [0, 0, 0.5], [0, 0, 1], [0, 0.5, 1],[0, 0.5, 0],  [0, 1, 0]]

plot_ranges=[[-0.0001,0.008], [-0.1e-6,2e-6],  [-0.00005, 0.0025],  [-0.0001,0.008], [-0.05e-6,0.8e-6], [-0.1e-7,1.6e-7],   [6,18],     [-1,450],        [-1,29],         [300,370],       [300,370], [300,370], [-0.1e-6,3.8e-6], [-0.001,0.015], [-0.005,0.08], [80,140]]
var_names=['Phyc_allmonths','Chl_allmonths','Phydiat_allmonths','Zooc_allmonths','PPint_allmonths','EPC100_allmonths','SST_allmonths','MLD_allmonths','FgCO2_allmonths','SpCO2_allmonths','SpCO2_T','SpCO2_nonT','Fe_allmonths','NO3_allmonths','Si_allmonths','RSNTDS_allmonths']
plot_titles=['Biomass (mol/m3)','Chlorophyll (kg.m-3)','Diatom (mol/m3)', 'Zooplankton (mol/m3)','Primary Prod. (mol.m-2.s-1)','Export Prod. (mol.m-2.s-1)',   'SST (C)',   'MLD (m)',  'Air-sea CO2 flux(gram C /m2/year)', 'SpCO2 (ppm)', 'SpCO2-T (ppm)', 'SpCO2-nonT (ppm)', 'Iron - Fe (mol/m3)', 'Nitrate - NO3 (mol/m3)','Silica - Si (mol/m3)', 'Light (W/m2)']  
file_names=[ 'bio',     'bio',      'bio',      'bio',        'bio',              'bio',                      'physics','physics', 'bio',              'bio',   'bio',      'bio',       'bio',      'bio',          'bio',       'physics'      ]

d_trend='yes'
n_r=4 ; n_c=4 ;
n_t=len(var_names)

fig=plt.figure()
for V_i in range(n_t): # V_i=0   V_i=1    V_i=2    V_i=3    V_i=4    V_i=5    V_i=6    V_i=7    V_i=8    V_i=9    V_i=10    V_i=11    V_i=12    V_i=13

    ax = fig.add_subplot(n_r,n_c,V_i+1) 

    for M_i in range(len(GCM_Names)): # M_i=0     # GCM=GCM_Names[0]

        GCM=GCM_Names[M_i]
        print (GCM) 
   
        if var_names[V_i]=='SpCO2_T' or var_names[V_i]=='SpCO2_nonT': # Only for SpCO2 sepration plots
                
            filename_in_spco2 = (dir_data_in + 'AllResults_'+GCM+'_'+'bio'+'_1991_2010.out') # Directory to save processed data
            filename_in_sst = (dir_data_in + 'AllResults_'+GCM+'_'+'physics'+'_1991_2010.out') # Directory to save processed data
            
            my_shelf_spco2 = shelve.open(filename_in_spco2)
            my_shelf_sst = shelve.open(filename_in_sst)
            
            if 'SpCO2_allmonths' in my_shelf_spco2  and 'SST_allmonths' in my_shelf_sst: # If that variable exists in the saved data  
                
                globals()['Var_spco2']=my_shelf_spco2['SpCO2_allmonths']
                globals()['Var_sst']=my_shelf_sst['SST_allmonths']                

                Var_spco2[ np.abs(Var_spco2)>1e19 ]=nan; Var_sst[ np.abs(Var_sst)>1e19 ]=nan  # Replacing missing values with NaN 

                Var_spco2_PAPA = Var_spco2[:,135:140,210:220]
                Var_sst_PAPA = Var_sst[:,135:140,210:220]                
                Var_spco2_NABE = Var_spco2[:,135:140,325:335]
                Var_sst_NABE = Var_sst[:,135:140,325:335]              

                if d_trend=='yes': # If it is requested to detrend the data before correlations
                    Var_spco2_PAPA= func_detrend_3d(Var_spco2_PAPA)
                    Var_sst_PAPA= func_detrend_3d(Var_sst_PAPA)
                    Var_spco2_NABE= func_detrend_3d(Var_spco2_NABE)
                    Var_sst_NABE= func_detrend_3d(Var_sst_NABE)
              
                Var_spco2_TEMP_eq1_PAPA = empty((Var_spco2_PAPA.shape[0], Var_spco2_PAPA.shape[1],Var_spco2_PAPA.shape[2])) * nan #SpCO2_allmonths_TEMP_eq1 = np.zeros([SpCO2_allmonths.shape[0], lat_n_regrid,lon_n_regrid])
                Var_spco2_NonTEMP_eq2_PAPA = empty((Var_spco2_PAPA.shape[0], Var_spco2_PAPA.shape[1],Var_spco2_PAPA.shape[2])) * nan #SpCO2_allmonths_NonTEMP_eq2 = np.zeros([SpCO2_allmonths.shape[0], lat_n_regrid,lon_n_regrid])
                for ii in range(Var_spco2_PAPA.shape[0]):
                    Var_spco2_TEMP_eq1_PAPA[ii,:,:] = np.nanmean(Var_spco2_PAPA,axis=0) * ( np.exp( 0.0423 * ( Var_sst_PAPA[ii,:,:] - np.nanmean(Var_sst_PAPA,axis=0) ) ) );
                    Var_spco2_NonTEMP_eq2_PAPA[ii,:,:] = Var_spco2_PAPA[ii,:,:] * ( np.exp( 0.0423 * ( np.nanmean(Var_sst_PAPA,axis=0) - Var_sst_PAPA[ii,:,:] ) ) );

                Var_spco2_TEMP_eq1_NABE = empty((Var_spco2_NABE.shape[0], Var_spco2_NABE.shape[1],Var_spco2_NABE.shape[2])) * nan #SpCO2_allmonths_TEMP_eq1 = np.zeros([SpCO2_allmonths.shape[0], lat_n_regrid,lon_n_regrid])
                Var_spco2_NonTEMP_eq2_NABE = empty((Var_spco2_NABE.shape[0], Var_spco2_NABE.shape[1],Var_spco2_NABE.shape[2])) * nan #SpCO2_allmonths_NonTEMP_eq2 = np.zeros([SpCO2_allmonths.shape[0], lat_n_regrid,lon_n_regrid])
                for ii in range(Var_spco2_NABE.shape[0]):
                    Var_spco2_TEMP_eq1_NABE[ii,:,:] = np.nanmean(Var_spco2_NABE,axis=0) * ( np.exp( 0.0423 * ( Var_sst_NABE[ii,:,:] - np.nanmean(Var_sst_NABE,axis=0) ) ) );
                    Var_spco2_NonTEMP_eq2_NABE[ii,:,:] = Var_spco2_NABE[ii,:,:] * ( np.exp( 0.0423 * ( np.nanmean(Var_sst_NABE,axis=0) - Var_sst_NABE[ii,:,:] ) ) );
                
                if var_names[V_i]=='SpCO2_T': # Only for SpCO2 sepration plots                                       
                    Variable_1_PAPA=copy.deepcopy(Var_spco2_TEMP_eq1_PAPA)
                    Variable_1_NABE=copy.deepcopy(Var_spco2_TEMP_eq1_NABE)                                      
                    
                if var_names[V_i]=='SpCO2_nonT': # Only for SpCO2 sepration plots       
                    Variable_1_PAPA=copy.deepcopy(Var_spco2_NonTEMP_eq2_PAPA)
                    Variable_1_NABE=copy.deepcopy(Var_spco2_NonTEMP_eq2_NABE) 
                    
                Var_plot_PAPA=np.nanmean(np.nanmean(np.nanmean(Variable_1_PAPA,axis=0) ,axis=0) ,axis=0)
                Var_plot_NABE=np.nanmean(np.nanmean(np.nanmean(Variable_1_NABE,axis=0) ,axis=0) ,axis=0)
                    
            else: # If we don't have both datasets for this model to do the correlation, just plot an empty map, to keep the order of GCMs constant in all models                
                Var_plot_PAPA=nan; Var_plot_NABE=nan; 

        else:

            filename_in = (dir_data_in + 'AllResults_'+GCM+'_'+file_names[V_i]+'_1991_2010.out') # Directory to save processed data
            my_shelf = shelve.open(filename_in)
     
            if var_names[V_i] in my_shelf: # If that variable exists in the saved data  
                
                globals()['Variable_1']=my_shelf[var_names[V_i]]       
                
                Variable_1[ np.abs(Variable_1)>1e19 ]=nan; # Replacing missing values with NaN 
                
                Variable_1_PAPA = Variable_1[:,135:140,210:220]           
                Variable_1_NABE = Variable_1[:,135:140,325:335]

                if d_trend=='yes': # If it is requested to detrend the data before correlations
                    Variable_1_PAPA= func_detrend_3d(Variable_1_PAPA)
                    Variable_1_NABE= func_detrend_3d(Variable_1_NABE)

                Var_plot_PAPA=np.nanmean(np.nanmean(np.nanmean(Variable_1_PAPA,axis=0) ,axis=0) ,axis=0)
                Var_plot_NABE=np.nanmean(np.nanmean(np.nanmean(Variable_1_NABE,axis=0) ,axis=0) ,axis=0)
                
            else: # If we don't have both datasets for this model to do the correlation, just plot an empty map, to keep the order of GCMs constant in all models
                
                Var_plot_PAPA=nan; Var_plot_NABE=nan; 
                
        plt.scatter(Var_plot_PAPA, Var_plot_NABE, s=100, marker='o', c=Colors[M_i], label=GCM) # Only the first pannel gets the lables so the legend location be alligned with that
        plt.scatter(Var_plot_PAPA, Var_plot_NABE, s=100, marker='o', facecolors='none', edgecolors='k')
        
    plt.plot([-1000,1000],[-1000,1000], 'k--', linewidth=0.75)
    plt.xlim(plot_ranges[V_i][0],plot_ranges[V_i][1]) ; plt.ylim(plot_ranges[V_i][0],plot_ranges[V_i][1])
    plt.xticks(fontsize=10); plt.yticks(fontsize = 10) # Setting the X-ticks ## rotation='horizontal' # rotation='vertical'   
    plt.title(plot_titles[V_i], fontsize=14)

plt.subplots_adjust(left=0.18, bottom=0.06, right=0.83, top=0.92, hspace=0.35, wspace=0.25) # the amount of height/width reserved for space between subplots
plt.legend(shadow=True, loc='right', bbox_to_anchor=(2.1, 2.5), markerscale=1 , prop={'size': 14})
if d_trend=='yes':
    plt.suptitle('CMIP5 models - Average value for North Pacific (x-axis) and North Atlantic (y-axis) - '+str(start_date_plt)+'-'+str(end_date_plt)+' - detrended' , fontsize=18)
else:
    plt.suptitle('CMIP5 models - Average value for North Pacific (x-axis) and North Atlantic (y-axis) - '+str(start_date_plt)+'-'+str(end_date_plt) , fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full

if d_trend=='yes':
    fig.savefig(dir_figs+'AllModels_Average_Scatter_detrended_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
else:
    fig.savefig(dir_figs+'AllModels_Average_Scatter_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

##############################################################
###   Average variable map Multi-GCM Plotting - Subplots  ####
##############################################################
GCM_Names = ['CESM1-BGC' ,'GFDL-ESM2G','GFDL-ESM2M','IPSL-CM5A-MR','IPSL-CM5A-LR','MRI-ESM1','GISS-E2-H-CC','GISS-E2-R-CC','HadGEM2-ES','CanESM2','NorESM1-ME','MPI-ESM-MR','MPI-ESM-LR']
Colors    =  [[1, 1, 0],   [1, 0.8, 0],  [1, 0.6, 0], [1, 0.4, 0],   [1, 0, 0], [1, 0.8, 1],   [1, 0.5, 1],   [1, 0, 1],   [0, 0, 0.5], [0, 0, 1], [0, 0.5, 1],[0, 0.5, 0],  [0, 1, 0]]
 
var_names=['Phyc_allmonths','Chl_allmonths','Phydiat_allmonths','Zooc_allmonths','Zoo_by_Phyto','PPint_allmonths','EPC100_allmonths','FgCO2_allmonths','SpCO2_allmonths','Fe_allmonths','NO3_allmonths','Si_allmonths','SST_allmonths', 'MLD_allmonths', 'HFDS_allmonths','SpCO2_T', 'SpCO2_nonT']
plot_titles=['Biomass (mol/m3)','Chlorophyll (kg.m-3)','Diatom (mol/m3)', 'Zooplankton (mol/m3)', 'Zoop/Phyto (-)','Primary Production (mol.m-2.s-1)','Export Production at 100m (mol.m-2.s-1)','Air-sea CO2 flux (donward +) (gram C /m2/year)','Surface pCO2 (ppm)', 'Iron - Fe (mol/m3)', 'Nitrate - NO3 (mol/m3)','Silica - Si (mol/m3)', 'Sea Surface Temperature (C)', 'Mixed Layer Depth (m)', 'Downward Heat Flux at Sea Water Surface (W.m-2)', 'pCO2-T (ppm)', 'pCO2-nonT (ppm)'                                   ]  
file_names=['bio',  'bio',   'bio',    'bio',     'bio',   'bio',  'bio',   'bio',   'bio',   'bio',   'bio',  'bio', 'physics', 'physics', 'physics',   'bio',  'bio'  ]
save_names= [ 'Phyc','Chl', 'Phydiat', 'Zooc', 'ZoobyPhyto', 'PPint','EPC100', 'FgCO2', 'SpCO2', 'Fe', 'NO3', 'Si', 'SST', 'MLD', 'HFDS' , 'SpCO2_Temp', 'SpCO2_nonT'   ]  
plot_ranges=[[0,0.006], [0,8e-7],  [0, 0.0008],  [0,0.006],  [0,2], [0,1e-6], [0,1.6e-7],   [-12,12],   [300,450],    [0,3.8e-6], [0,0.015], [0,0.06], [0,32], [0,300], [-120,120],[300,450],[300,450]]

P_cmap=plt.cm.bwr; P_c_fill='fill'; P_Lon=Lon_regrid_2D; P_Lat=Lat_regrid_2D; P_range=51;
P_proj='cyl'; P_lon0=210.; P_latN=90.; P_latS=-90.; 

d_trend='NO'

n_r=4 ; n_c=4 ;
n_t=len(GCM_Names)

for V_i in range(len(var_names_y)): # V_i=0   V_i=1    V_i=2    V_i=3    V_i=4    V_i=5    V_i=6    V_i=7    V_i=8    V_i=9    V_i=10    V_i=11    V_i=12    V_i=13    V_i=14    V_i=15    V_i=16
    
    if d_trend=='yes':
        P_title='Average map of '+plot_titles[V_i]+' - CMIP5 models - '+str(start_date_plt)+'-'+str(end_date_plt)+' - detrended'
    else:
        P_title='Average map of '+plot_titles[V_i]+' - CMIP5 models - '+str(start_date_plt)+'-'+str(end_date_plt)

    fig=plt.figure()
    for M_i in range(len(GCM_Names)): # M_i=0     # GCM=GCM_Names[0]
        
        ax = fig.add_subplot(n_r,n_c,M_i+1)  
    
        GCM=GCM_Names[M_i]
        print (GCM) 
   
        if var_names[V_i]=='SpCO2_T' or var_names[V_i]=='SpCO2_nonT': # Only for SpCO2 sepration plots
                
            filename_in_spco2 = (dir_data_in + 'AllResults_'+GCM+'_'+'bio'+'_1991_2010.out') # Directory to save processed data
            filename_in_sst = (dir_data_in + 'AllResults_'+GCM+'_'+'physics'+'_1991_2010.out') # Directory to save processed data

            my_shelf_spco2 = shelve.open(filename_in_spco2)
            my_shelf_sst = shelve.open(filename_in_sst)
            
            if 'SpCO2_allmonths' in my_shelf_spco2  and 'SST_allmonths' in my_shelf_sst: # If that variable exists in the saved data  
                
                globals()['Var_spco2']=my_shelf_spco2['SpCO2_allmonths']
                globals()['Var_sst']=my_shelf_sst['SST_allmonths']                 
                               
                Var_spco2[ np.abs(Var_spco2)>1e19 ]=nan; Var_sst[ np.abs(Var_sst)>1e19 ]=nan  # Replacing missing values with NaN 

                if d_trend=='yes': # If it is requested to detrend the data before correlations
                    Var_spco2= func_detrend_3d(Var_spco2)
                    Var_sst= func_detrend_3d(Var_sst)
                
                Var_spco2_TEMP_eq1 = empty((Var_spco2.shape[0], Var_spco2.shape[1],Var_spco2.shape[2])) * nan #SpCO2_allmonths_TEMP_eq1 = np.zeros([SpCO2_allmonths.shape[0], lat_n_regrid,lon_n_regrid])
                Var_spco2_NonTEMP_eq2 = empty((Var_spco2.shape[0], Var_spco2.shape[1],Var_spco2.shape[2])) * nan #SpCO2_allmonths_NonTEMP_eq2 = np.zeros([SpCO2_allmonths.shape[0], lat_n_regrid,lon_n_regrid])
                for ii in range(Var_spco2.shape[0]):
                    Var_spco2_TEMP_eq1[ii,:,:] = np.nanmean(Var_spco2,axis=0) * ( np.exp( 0.0423 * ( Var_sst[ii,:,:] - np.nanmean(Var_sst,axis=0) ) ) );
                    Var_spco2_NonTEMP_eq2[ii,:,:] = Var_spco2[ii,:,:] * ( np.exp( 0.0423 * ( np.nanmean(Var_sst,axis=0) - Var_sst[ii,:,:] ) ) );
                
                if var_names[V_i]=='SpCO2_T' : # Only for SpCO2 sepration plots   
                    Var_plot_ii=np.nanmean(Var_spco2_TEMP_eq1, axis=0)
                if var_names[V_i]=='SpCO2_nonT': # Only for SpCO2 sepration plots     
                    Var_plot_ii=np.nanmean(Var_spco2_NonTEMP_eq2, axis=0)
  
                Var_plot_ii=copy.deepcopy(Corr_xy)                
                color_land='k'
                
            else: # If we don't have both datasets for this model to do the correlation, just plot an empty map, to keep the order of GCMs constant in all models                
                Var_plot_ii=empty((lat_n_regrid,lon_n_regrid)) * nan
                color_land='w'


        elif var_names[V_i]=='Zoo_by_Phyto' :

            filename_in = (dir_data_in + 'AllResults_'+GCM+'_'+'bio'+'_1991_2010.out') # Directory to save processed data
            my_shelf = shelve.open(filename_in)
     
            if 'Zooc_allmonths' in my_shelf and 'Phyc_allmonths' in my_shelf: # If that variable exists in the saved data  
                
                globals()['Variable_1']=my_shelf['Zooc_allmonths'] 
                globals()['Variable_2']=my_shelf['Phyc_allmonths'] 

                Variable_1[ np.abs(Variable_1)>1e19 ]=nan # Replacing missing values with NaN
                Variable_2[ np.abs(Variable_2)>1e19 ]=nan # Replacing missing values with NaN
                    
                if d_trend=='yes': # If it is requested to detrend the data before correlations
                    Variable_1= func_detrend_3d(Variable_1)
                    Variable_2= func_detrend_3d(Variable_2)
    
                Var_plot_ii=np.nanmean( Variable_1 / Variable_2, axis=0)
                
                color_land='k'
            else: # If we don't have both datasets for this model to do the correlation, just plot an empty map, to keep the order of GCMs constant in all models
                
                Var_plot_ii=empty((lat_n_regrid,lon_n_regrid)) * nan
                color_land='w'        

        else:

            filename_in = (dir_data_in + 'AllResults_'+GCM+'_'+file_names[V_i]+'_1991_2010.out') # Directory to save processed data
            my_shelf = shelve.open(filename_in)
     
            if var_names[V_i] in my_shelf: # If that variable exists in the saved data  
                
                globals()['Variable_1']=my_shelf[var_names[V_i]]      

                Variable_1[ np.abs(Variable_1)>1e19 ]=nan # Replacing missing values with NaN 
                    
                if d_trend=='yes': # If it is requested to detrend the data before correlations
                    Variable_1= func_detrend_3d(Variable_1)            
    
                Var_plot_ii=np.nanmean(Variable_1, axis=0)
                
                color_land='k'
            else: # If we don't have both datasets for this model to do the correlation, just plot an empty map, to keep the order of GCMs constant in all models
                
                Var_plot_ii=empty((lat_n_regrid,lon_n_regrid)) * nan
                color_land='w'


        if var_names[V_i]=='HFDS_allmonths' or var_names[V_i]=='FgCO2_allmonths':
            P_cmap=plt.cm.bwr
        else:
            P_cmap=plt.cm.jet
            
        #P_range=51;
        P_range=np.linspace(plot_ranges[V_i][0],plot_ranges[V_i][1],51);
        
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
        
        if P_c_fill=='fill':
            m.fillcontinents(color='0')
        if color_land=='w':
            m.fillcontinents(color='1')
        
        m.drawcoastlines(linewidth=1.0, linestyle='solid', color=color_land, antialiased=1, ax=None, zorder=None)
        im=m.contourf(P_Lon, P_Lat, Var_plot_ii,P_range,latlon=True, cmap=P_cmap, extend='both')
        plt.title(GCM, fontsize=14)
            
    plt.subplots_adjust(left=0.05, bottom=0.05, right=0.93, top=0.92, hspace=0.3, wspace=0.02) # the amount of height/width reserved for space between subplots
    cbar_ax = fig.add_axes([0.94, 0.05, 0.015, 0.87]) # [right,bottom,width,height] 
    fig.colorbar(im, cax=cbar_ax)
    #fig.colorbar(im, cax=cbar_ax, ticks=[-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1])
    #cbar_ax.set_yticklabels(['- 1', '- 0.8', '- 0.6', '- 0.4', '- 0.2', '0', '0.2', '0.4', '0.6', '0.8', '1'])
    plt.suptitle(P_title, fontsize=20)    
    plt.show()
    mng = plt.get_current_fig_manager()
    mng.window.showMaximized() # Maximizes the plot window to save figures in full   

    if d_trend=='yes':
        fig.savefig(dir_figs+'AllModels_Average_Maps_detrended_'+save_names[V_i]+'_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
    else:
        fig.savefig(dir_figs+'AllModels_Average_Maps_'+save_names[V_i]+'_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


###################################################################
###   Variable seasonality maps Multi-GCM Plotting - Subplots  ####
###################################################################
GCM_Names = ['CESM1-BGC' ,'GFDL-ESM2G','GFDL-ESM2M','IPSL-CM5A-MR','IPSL-CM5A-LR','MRI-ESM1','GISS-E2-H-CC','GISS-E2-R-CC','HadGEM2-ES','CanESM2','NorESM1-ME','MPI-ESM-MR','MPI-ESM-LR']
Colors    =  [[1, 1, 0],   [1, 0.8, 0],  [1, 0.6, 0], [1, 0.4, 0],   [1, 0, 0], [1, 0.8, 1],   [1, 0.5, 1],   [1, 0, 1],   [0, 0, 0.5], [0, 0, 1], [0, 0.5, 1],[0, 0.5, 0],  [0, 1, 0]]
 
var_names=['Phyc_allmonths','Chl_allmonths','Phydiat_allmonths','Zooc_allmonths','Zoo_by_Phyto','PPint_allmonths','EPC100_allmonths','FgCO2_allmonths','SpCO2_allmonths','Fe_allmonths','NO3_allmonths','Si_allmonths','SST_allmonths', 'MLD_allmonths', 'HFDS_allmonths','SpCO2_T', 'SpCO2_nonT']
plot_titles=['Biomass','Chlorophyll','Diatom', 'Zooplankton', 'Zoop/Phyto','Primary Production','Export Production at 100m','Air-sea CO2 flux','Surface pCO2', 'Iron - Fe', 'Nitrate - NO3','Silica - Si', 'Sea Surface Temperature', 'Mixed Layer Depth', 'Downward Heat Flux at Sea Water Surface', 'pCO2-T', 'pCO2-nonT'                                   ]  
file_names=['bio',  'bio',   'bio',    'bio',     'bio',   'bio',  'bio',   'bio',   'bio',   'bio',   'bio',  'bio', 'physics', 'physics', 'physics',   'bio',  'bio'  ]
save_names= [ 'Phyc','Chl', 'Phydiat', 'Zooc', 'ZoobyPhyto', 'PPint','EPC100', 'FgCO2', 'SpCO2', 'Fe', 'NO3', 'Si', 'SST', 'MLD', 'HFDS' , 'SpCO2_Temp', 'SpCO2_nonT'   ]  
plot_ranges=[[-3,3], [0,8e-7],  [0, 0.0008],  [0,0.006],  [0,2], [0,1e-6], [0,1.6e-7],   [0,1000],   [300,450],    [0,3.8e-6], [0,0.015], [0,0.06], [0,32], [0,300], [-120,120],[300,450],[300,450]]

P_cmap=plt.cm.bwr; P_c_fill='fill'; P_Lon=Lon_regrid_2D; P_Lat=Lat_regrid_2D; P_range=51;
P_proj='cyl'; P_lon0=210.; P_latN=90.; P_latS=-90.; 

d_trend='no'

n_r=4 ; n_c=4 ;
n_t=len(GCM_Names)

for V_i in range(len(var_names_y)): # V_i=0   V_i=1    V_i=2    V_i=3    V_i=4    V_i=5    V_i=6    V_i=7    V_i=8    V_i=9    V_i=10    V_i=11    V_i=12    V_i=13    V_i=14    V_i=15    V_i=16
    
    if d_trend=='yes':
        P_title='Seasonality [(max-min)/ave] map of '+plot_titles[V_i]+' - CMIP5 models - '+str(start_date_plt)+'-'+str(end_date_plt)+' - detrended'
    else:
        P_title='Seasonality [(max-min)/ave] map of '+plot_titles[V_i]+' - CMIP5 models - '+str(start_date_plt)+'-'+str(end_date_plt)

    fig=plt.figure()
    for M_i in range(len(GCM_Names)): # M_i=0     # GCM=GCM_Names[0]
        
        ax = fig.add_subplot(n_r,n_c,M_i+1)  
    
        GCM=GCM_Names[M_i]
        print (GCM) 
   
        if var_names[V_i]=='SpCO2_T' or var_names[V_i]=='SpCO2_nonT': # Only for SpCO2 sepration plots
            
            Var_plot_ii=empty((lat_n_regrid,lon_n_regrid)) * nan
            filename_in_spco2 = (dir_data_in + 'AllResults_'+GCM+'_'+'bio'+'_1991_2010.out') # Directory to save processed data
            filename_in_sst = (dir_data_in + 'AllResults_'+GCM+'_'+'physics'+'_1991_2010.out') # Directory to save processed data

            my_shelf_spco2 = shelve.open(filename_in_spco2)
            my_shelf_sst = shelve.open(filename_in_sst)
            
            if 'SpCO2_allmonths' in my_shelf_spco2  and 'SST_allmonths' in my_shelf_sst: # If that variable exists in the saved data  
                
                globals()['Var_spco2']=my_shelf_spco2['SpCO2_allmonths']
                globals()['Var_sst']=my_shelf_sst['SST_allmonths']                 
                               
                Var_spco2[ np.abs(Var_spco2)>1e19 ]=nan; Var_sst[ np.abs(Var_sst)>1e19 ]=nan  # Replacing missing values with NaN 

                if d_trend=='yes': # If it is requested to detrend the data before correlations
                    Var_spco2= func_detrend_3d(Var_spco2)
                    Var_sst= func_detrend_3d(Var_sst)
                
                Var_spco2_TEMP_eq1 = empty((Var_spco2.shape[0], Var_spco2.shape[1],Var_spco2.shape[2])) * nan #SpCO2_allmonths_TEMP_eq1 = np.zeros([SpCO2_allmonths.shape[0], lat_n_regrid,lon_n_regrid])
                Var_spco2_NonTEMP_eq2 = empty((Var_spco2.shape[0], Var_spco2.shape[1],Var_spco2.shape[2])) * nan #SpCO2_allmonths_NonTEMP_eq2 = np.zeros([SpCO2_allmonths.shape[0], lat_n_regrid,lon_n_regrid])
                for ii in range(Var_spco2.shape[0]):
                    Var_spco2_TEMP_eq1[ii,:,:] = np.nanmean(Var_spco2,axis=0) * ( np.exp( 0.0423 * ( Var_sst[ii,:,:] - np.nanmean(Var_sst,axis=0) ) ) );
                    Var_spco2_NonTEMP_eq2[ii,:,:] = Var_spco2[ii,:,:] * ( np.exp( 0.0423 * ( np.nanmean(Var_sst,axis=0) - Var_sst[ii,:,:] ) ) );
                
                if var_names[V_i]=='SpCO2_T' : # Only for SpCO2 sepration plots   
                    Variable_1=Var_spco2_TEMP_eq1
                if var_names[V_i]=='SpCO2_nonT': # Only for SpCO2 sepration plots     
                    Variable_1=Var_spco2_NonTEMP_eq2
  
                Variable_1_monthlymean = Variable_1.reshape(( np.int(Variable_1.shape[0]/12) ,12,Variable_1.shape[1],Variable_1.shape[2])) # Monthly means of SST over the time period
                Variable_1_monthlymean = np.nanmean(Variable_1_monthlymean,axis=0)
             
                Var_plot_ii= ( np.nanmax(Variable_1_monthlymean, axis=0) - np.nanmin(Variable_1_monthlymean, axis=0) ) / np.nanmean(Variable_1_monthlymean, axis=0)
                Var_plot_ii[Ocean_Index==10]=nan;   #  Excluding the Arctic Ocean   
                   
                color_land='k'
                
            else: # If we don't have both datasets for this model to do the correlation, just plot an empty map, to keep the order of GCMs constant in all models                
                Var_plot_ii=empty((lat_n_regrid,lon_n_regrid)) * nan
                color_land='w'

        elif var_names[V_i]=='Zoo_by_Phyto' :
            
            Var_plot_ii=empty((lat_n_regrid,lon_n_regrid)) * nan
            filename_in = (dir_data_in + 'AllResults_'+GCM+'_'+'bio'+'_1991_2010.out') # Directory to save processed data
            my_shelf = shelve.open(filename_in)
     
            if 'Zooc_allmonths' in my_shelf and 'Phyc_allmonths' in my_shelf: # If that variable exists in the saved data  
                
                globals()['Variable_1']=my_shelf['Zooc_allmonths'] 
                globals()['Variable_2']=my_shelf['Phyc_allmonths'] 

                Variable_1[ np.abs(Variable_1)>1e19 ]=nan # Replacing missing values with NaN
                Variable_2[ np.abs(Variable_2)>1e19 ]=nan # Replacing missing values with NaN
                    
                if d_trend=='yes': # If it is requested to detrend the data before correlations
                    Variable_1= func_detrend_3d(Variable_1)
                    Variable_2= func_detrend_3d(Variable_2)
    
                Variable_1= Variable_1 / Variable_2

                Variable_1_monthlymean = Variable_1.reshape(( np.int(Variable_1.shape[0]/12) ,12,Variable_1.shape[1],Variable_1.shape[2])) # Monthly means of SST over the time period
                Variable_1_monthlymean = np.nanmean(Variable_1_monthlymean,axis=0)
             
                Var_plot_ii= ( np.nanmax(Variable_1_monthlymean, axis=0) - np.nanmin(Variable_1_monthlymean, axis=0) ) / np.nanmean(Variable_1_monthlymean, axis=0)
                Var_plot_ii[Ocean_Index==10]=nan;   #  Excluding the Arctic Ocean
                
                color_land='k'
            else: # If we don't have both datasets for this model to do the correlation, just plot an empty map, to keep the order of GCMs constant in all models
                
                Var_plot_ii=empty((lat_n_regrid,lon_n_regrid)) * nan
                color_land='w'        

        else:
            
            Var_plot_ii=empty((lat_n_regrid,lon_n_regrid)) * nan
            filename_in = (dir_data_in + 'AllResults_'+GCM+'_'+file_names[V_i]+'_1991_2010.out') # Directory to save processed data
            my_shelf = shelve.open(filename_in)
     
            if var_names[V_i] in my_shelf: # If that variable exists in the saved data  
                
                globals()['Variable_1']=my_shelf[var_names[V_i]]      

                Variable_1[ np.abs(Variable_1)>1e19 ]=nan # Replacing missing values with NaN 
                    
                if d_trend=='yes': # If it is requested to detrend the data before correlations
                    Variable_1= func_detrend_3d(Variable_1)   
                
                Variable_1_monthlymean = Variable_1.reshape(( np.int(Variable_1.shape[0]/12) ,12,Variable_1.shape[1],Variable_1.shape[2])) # Monthly means of SST over the time period
                Variable_1_monthlymean = np.nanmean(Variable_1_monthlymean,axis=0)
             
                Var_plot_ii= ( np.nanmax(Variable_1_monthlymean, axis=0) - np.nanmin(Variable_1_monthlymean, axis=0) ) / np.nanmean(Variable_1_monthlymean, axis=0)
                Var_plot_ii[Ocean_Index==10]=nan;   #  Excluding the Arctic Ocean
                
#                Var_plot_ii=empty((lat_n_regrid,lon_n_regrid)) * nan
#                Var_plot_ii[0:90,:]= ( np.nanmean(Variable_1_monthlymean[0:2,0:90,:], axis=0) - np.nanmean(Variable_1_monthlymean[6:8,0:90,:], axis=0) ) / np.nanmean(Variable_1_monthlymean[:,0:90,:], axis=0) 
#                Var_plot_ii[90:180,:]= ( np.nanmean(Variable_1_monthlymean[6:8,90:180,:], axis=0) - np.nanmean(Variable_1_monthlymean[0:2,90:180,:], axis=0) ) / np.nanmean(Variable_1_monthlymean[:,90:180,:], axis=0)                  
#                Var_plot_ii[Ocean_Index==10]=nan;   #  Excluding the Arctic Ocean
                
                color_land='k'
            else: # If we don't have both datasets for this model to do the correlation, just plot an empty map, to keep the order of GCMs constant in all models
                
                Var_plot_ii=empty((lat_n_regrid,lon_n_regrid)) * nan
                color_land='w'

        P_cmap=plt.cm.jet
        P_range=51;
        #P_range=np.linspace(plot_ranges[V_i][0],plot_ranges[V_i][1],51);
        
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
        
        if P_c_fill=='fill':
            m.fillcontinents(color='0')
        if color_land=='w':
            m.fillcontinents(color='1')
        
        m.drawcoastlines(linewidth=1.0, linestyle='solid', color=color_land, antialiased=1, ax=None, zorder=None)
        im=m.contourf(P_Lon, P_Lat, Var_plot_ii,P_range,latlon=True, cmap=P_cmap, extend='both')
        plt.title(GCM, fontsize=14)
            
    plt.subplots_adjust(left=0.05, bottom=0.05, right=0.93, top=0.92, hspace=0.3, wspace=0.02) # the amount of height/width reserved for space between subplots
    cbar_ax = fig.add_axes([0.94, 0.05, 0.015, 0.87]) # [right,bottom,width,height] 
    fig.colorbar(im, cax=cbar_ax)
    #fig.colorbar(im, cax=cbar_ax, ticks=[-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1])
    #cbar_ax.set_yticklabels(['- 1', '- 0.8', '- 0.6', '- 0.4', '- 0.2', '0', '0.2', '0.4', '0.6', '0.8', '1'])
    plt.suptitle(P_title, fontsize=20)    
    plt.show()
    mng = plt.get_current_fig_manager()
    mng.window.showMaximized() # Maximizes the plot window to save figures in full   

    if d_trend=='yes':
        fig.savefig(dir_figs+'AllModels_Seasonality_Maps_detrended_'+save_names[V_i]+'_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
    else:
        fig.savefig(dir_figs+'AllModels_Seasonality_Maps_'+save_names[V_i]+'_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')



##################################
###   2D Plotting - Subplots  ####
##################################
GCM_Names = ['CESM1-BGC' ,'GFDL-ESM2G','GFDL-ESM2M','IPSL-CM5A-MR','IPSL-CM5A-LR','MRI-ESM1','GISS-E2-H-CC','GISS-E2-R-CC','HadGEM2-ES','CanESM2','NorESM1-ME','MPI-ESM-MR','MPI-ESM-LR']
Colors    =  [[1, 1, 0],   [1, 0.8, 0],  [1, 0.6, 0], [1, 0.4, 0],   [1, 0, 0], [1, 0.8, 1],   [1, 0.5, 1],   [1, 0, 1],   [0, 0, 0.5], [0, 0, 1], [0, 0.5, 1],[0, 0.5, 0],  [0, 1, 0]]

Time_months = ('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec');     
colors_months = ['navy', 'blue', 'royalblue', 'lime', 'limegreen', 'green', 'yellow', 'gold', 'orange', 'red', 'magenta','darkviolet']

#### Plots 1 #####
var_names_y=[ 'Phyc_allmonths',                'Phyc_allmonths',                    'Phyc_allmonths',        'Phyc_allmonths',      'Zooc_allmonths',            'Zooc_allmonths',                'Zooc_allmonths',    'Zooc_allmonths'       ]
var_names_x=[ 'Fe_allmonths',                  'NO3_allmonths',                     'RSNTDS_allmonths',      'MLD_allmonths',       'Fe_allmonths',              'NO3_allmonths',                 'RSNTDS_allmonths',  'MLD_allmonths'     ]
plot_titles=[ 'log biomass vs. log Iron (Fe)', 'log biomass vs. log Nitrate (NO3)', 'log biomass vs. light', 'log biomass vs. MLD', 'log Zoo vs. log Iron (Fe)', 'log Zoo vs. log Nitrate (NO3)', 'log Zoo vs. light', 'log Zoo vs. MLD'     ]
save_names= [ 'Phyc_vs_Fe',                    'Phyc_vs_NO3',                       'Phyc_vs_RSNTDS',        'Phyc_vs_MLD',         'Zooc_vs_Fe',                'Zooc_vs_NO3',                   'Zooc_vs_RSNTDS',    'Zooc_vs_MLD'       ]  

file_names_y=[ 'bio',                          'bio',                               'bio',                   'bio',                  'bio',                        'bio',                          'bio' ,             'bio'            ]
file_names_x=[ 'bio',                          'bio',                               'physics',               'physics',              'bio',                        'bio',                          'physics',          'physics'           ]

#### Plots 2 #####
var_names_y=[ 'PPint_allmonths',             'PPint_allmonths',                 'PPint_allmonths',     'PPint_allmonths',           'PPint_allmonths',   'EPC100_allmonths',               'EPC100_allmonths',                   'EPC100_allmonths',       'EPC100_allmonths'         ]
var_names_x=[ 'Fe_allmonths',                'NO3_allmonths',                   'RSNTDS_allmonths',    'Phyc_allmonths',            'MLD_allmonths',     'Fe_allmonths',                   'NO3_allmonths',                      'RSNTDS_allmonths',       'MLD_allmonths'               ]
plot_titles=[ 'log intPP vs. log Iron (Fe)', 'log intPP vs. log Nitrate (NO3)', 'log intPP vs. light', 'log intPP vs. log biomass', 'log intPP vs. MLD', 'log EPC 100m vs. log Iron (Fe)', 'log EPC 100m vs. log Nitrate (NO3)', 'log EPC 100m vs. light', 'log EPC 100m vs. MLD'  ]
save_names= [ 'PPint_vs_Fe',                 'PPint_vs_NO3',                    'PPint_vs_RSNTDS',     'PPint_vs_Phyc',             'PPint_vs_MLD',      'EPC100_vs_Fe',                   'EPC100_vs_NO3' ,                     'EPC100_vs_RSNTDS',       'EPC100_vs_MLD'                   ]  

file_names_y=[ 'bio',                         'bio',                             'bio',                'bio',                       'bio',                'bio',                            'bio',                                'bio',                   'bio'            ]
file_names_x=[ 'bio',                         'bio',                             'physics',            'bio',                       'physics',            'bio',                            'bio',                                'physics',               'physics'            ]

plot_ranges_y=[[0,0.006],  [-2, -7]]
plot_ranges_x=[[0,0.006],  [-1, -7]]


###   N PACIFIC and N ATLANTIC plotted separate - log scale
d_trend='yes'
n_r=4 ; n_c=4 ;
n_t=len(GCM_Names)
box_n='PAPA'   ;   box_n='NABE'


for V_i in range(n_t): # V_i=0   V_i=1    V_i=2    V_i=3    V_i=4    V_i=5    V_i=6    V_i=7    V_i=8    V_i=9    V_i=10    V_i=11    V_i=12    V_i=13

    if d_trend=='yes':
        if box_n=='PAPA':
            P_title=plot_titles[V_i]+' - North Pacific [45N-50N, 140W-150W] - '+str(start_date_plt)+'-'+str(end_date_plt)+' - detrended'
        elif box_n=='NABE':
            P_title=plot_titles[V_i]+' - North Atlantic [45N-50N, 25W-35W] - '+str(start_date_plt)+'-'+str(end_date_plt)+' - detrended'
    else:
        if box_n=='PAPA':
            P_title=plot_titles[V_i]+' - North Pacific [45N-50N, 140W-150W] - '+str(start_date_plt)+'-'+str(end_date_plt)
        elif box_n=='NABE':
            P_title=plot_titles[V_i]+' - North Atlantic [45N-50N, 25W-35W] - '+str(start_date_plt)+'-'+str(end_date_plt)
    
    fig=plt.figure()
    for M_i in range(len(GCM_Names)): # M_i=1     # GCM=GCM_Names[0]

        GCM=GCM_Names[M_i]
        print (GCM) 
        
        ax = fig.add_subplot(n_r,n_c,M_i+1) 

        filename_in_y = (dir_data_in + 'AllResults_'+GCM+'_'+file_names_y[V_i]+'_1991_2010.out') # Directory to save processed data
        filename_in_x = (dir_data_in + 'AllResults_'+GCM+'_'+file_names_x[V_i]+'_1991_2010.out') # Directory to save processed data
        my_shelf_y = shelve.open(filename_in_y)
        my_shelf_x = shelve.open(filename_in_x)
 
        if var_names_y[V_i] in my_shelf_y and var_names_x[V_i] in my_shelf_x: # If that variable exists in the saved data  
            
            globals()['Variable_y']=my_shelf_y[var_names_y[V_i]]
            globals()['Variable_x']=my_shelf_x[var_names_x[V_i]]        
            
            Variable_y[ np.abs(Variable_y)>1e19 ]=nan; Variable_x[ np.abs(Variable_x)>1e19 ]=nan  # Replacing missing values with NaN 

            if box_n=='PAPA':
                Variable_y_BOX = Variable_y[:,135:140,210:220]
                Variable_x_BOX = Variable_x[:,135:140,210:220]
            elif box_n=='NABE':
                Variable_y_BOX = Variable_y[:,135:140,325:335]
                Variable_x_BOX = Variable_x[:,135:140,325:335]
                    
            if d_trend=='yes': # If it is requested to detrend the data before correlations
                Variable_y_BOX= func_detrend_3d(Variable_y_BOX)
                Variable_x_BOX= func_detrend_3d(Variable_x_BOX)
            
            if var_names_y[V_i] !='RSNTDS_allmonths' and var_names_y[V_i] !='MLD_allmonths' : # For all variables except these ones, scatterplots are done on log10 of the variable
                Variable_y_BOX=np.log10(Variable_y_BOX)
            if var_names_x[V_i] !='RSNTDS_allmonths' and var_names_x[V_i] !='MLD_allmonths' : # For all variables except these ones, scatterplots are done on log10 of the variable
                Variable_x_BOX=np.log10(Variable_x_BOX)

            Variable_y_BOX_all_timeseries = Variable_y_BOX.reshape(( np.int(Variable_y_BOX.shape[0]*Variable_y_BOX.shape[1]*Variable_y_BOX.shape[2]) ))
            Variable_x_BOX_all_timeseries = Variable_x_BOX.reshape(( np.int(Variable_x_BOX.shape[0]*Variable_x_BOX.shape[1]*Variable_x_BOX.shape[2]) ))       

            Variable_y_BOX_reshaped = Variable_y_BOX.reshape(( yrs_n ,12,Variable_y_BOX.shape[1],Variable_y_BOX.shape[2])) # Monthly means of SpCO2 over the time period
            Variable_x_BOX_reshaped = Variable_x_BOX.reshape(( yrs_n ,12,Variable_x_BOX.shape[1],Variable_x_BOX.shape[2])) # Monthly means of SpCO2 over the time period        
    
            for ii in range(Variable_y_BOX_reshaped.shape[1]):
                plt.scatter(Variable_x_BOX_reshaped[:,ii,:,:], Variable_y_BOX_reshaped[:,ii,:,:], c= colors_months[ii] , s=1, marker='x', label=Time_months[ii] )
    #            ax.set_yscale('log', nonposy='clip')
    #            ax.set_xscale('log', nonposx='clip')
    #        plt.ylim(plot_ranges_y[V_i])
    #        plt.xlim(plot_ranges_x[V_i])
            plt.title(GCM, fontsize=14)
    
            if M_i == (n_c*(n_r-1)): # Adds longitude ranges only to the last subplots that appear at the bottom of plot
                plt.xlabel(var_names_x[V_i][:-10],fontsize=14)
                plt.ylabel(var_names_y[V_i][:-10],fontsize=14)
            elif M_i==0 or M_i==n_c or M_i==n_c*2 or M_i==n_c*3 or M_i==n_c*4 or M_i==n_c*5 or M_i==n_c*6 or M_i==n_c*7 or M_i==n_c*8:
                plt.ylabel(var_names_y[V_i][:-10],fontsize=14)
            elif M_i >= n_t-n_c and M_i != (n_c*(n_r-1)): # Adds longitude ranges only to the last subplots that appear at the bottom of plot            
                plt.xlabel(var_names_x[V_i][:-10],fontsize=14)

    plt.subplots_adjust(left=0.18, bottom=0.06, right=0.83, top=0.92, hspace=0.35, wspace=0.25) # the amount of height/width reserved for space between subplots
    plt.legend(shadow=True, loc='right', bbox_to_anchor=(5.5, 3), markerscale=16 , prop={'size': 16})
    plt.suptitle(P_title, fontsize=20)    
    plt.show()
    mng = plt.get_current_fig_manager()
    mng.window.showMaximized() # Maximizes the plot window to save figures in full
    
    if d_trend=='yes':
        if box_n=='PAPA':
            fig.savefig(dir_figs+'AllModels_ScatterPlots_PAPA_detrended_'+save_names[V_i]+'_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
        elif box_n=='NABE':
            fig.savefig(dir_figs+'AllModels_ScatterPlots_NABE_detrended_'+save_names[V_i]+'_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
    else:
        if box_n=='PAPA':
            fig.savefig(dir_figs+'AllModels_ScatterPlots_PAPA_'+save_names[V_i]+'_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
        elif box_n=='NABE':
            fig.savefig(dir_figs+'AllModels_ScatterPlots_NABE_'+save_names[V_i]+'_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


###   N PACIFIC and N ATLANTIC plotted together - log scale
            
for V_i in range(n_t): # V_i=0   V_i=1    V_i=2    V_i=3    V_i=4    V_i=5    V_i=6    V_i=7    V_i=8    V_i=9    V_i=10    V_i=11    V_i=12    V_i=13

    if d_trend=='yes':
        P_title=plot_titles[V_i]+' - North Pacific (circles) and North Atlantic (diamonds) - '+str(start_date_plt)+'-'+str(end_date_plt)+' - detrended'
    else:
        P_title=plot_titles[V_i]+' - North Pacific (circles) and North Atlantic (diamonds) - '+str(start_date_plt)+'-'+str(end_date_plt)
    
    fig=plt.figure()
    for M_i in range(len(GCM_Names)): # M_i=1     # GCM=GCM_Names[0]

        GCM=GCM_Names[M_i]
        print (GCM) 
        
        ax = fig.add_subplot(n_r,n_c,M_i+1) 

        filename_in_y = (dir_data_in + 'AllResults_'+GCM+'_'+file_names_y[V_i]+'_1991_2010.out') # Directory to save processed data
        filename_in_x = (dir_data_in + 'AllResults_'+GCM+'_'+file_names_x[V_i]+'_1991_2010.out') # Directory to save processed data
        my_shelf_y = shelve.open(filename_in_y)
        my_shelf_x = shelve.open(filename_in_x)
 
        if var_names_y[V_i] in my_shelf_y and var_names_x[V_i] in my_shelf_x: # If that variable exists in the saved data  
            
            globals()['Variable_y']=my_shelf_y[var_names_y[V_i]]
            globals()['Variable_x']=my_shelf_x[var_names_x[V_i]]        
            
            Variable_y[ np.abs(Variable_y)>1e19 ]=nan; Variable_x[ np.abs(Variable_x)>1e19 ]=nan  # Replacing missing values with NaN 

            Variable_x_PAPA = Variable_x[:,135:140,210:220]
            Variable_y_PAPA = Variable_y[:,135:140,210:220]

            Variable_x_NABE = Variable_x[:,135:140,325:335]
            Variable_y_NABE = Variable_y[:,135:140,325:335]
                    
            if d_trend=='yes': # If it is requested to detrend the data before correlations
                Variable_x_PAPA= func_detrend_3d(Variable_x_PAPA)
                Variable_y_PAPA= func_detrend_3d(Variable_y_PAPA)                
                Variable_x_NABE= func_detrend_3d(Variable_x_NABE)
                Variable_y_NABE= func_detrend_3d(Variable_y_NABE)
            
            if var_names_x[V_i] !='RSNTDS_allmonths' and var_names_x[V_i] !='MLD_allmonths' : # For all variables except these ones, scatterplots are done on log10 of the variable
                Variable_x_PAPA=np.log10(Variable_x_PAPA) 
                Variable_x_NABE=np.log10(Variable_x_NABE) 
            if var_names_y[V_i] !='RSNTDS_allmonths' and var_names_y[V_i] !='MLD_allmonths' : # For all variables except these ones, scatterplots are done on log10 of the variable
                Variable_y_PAPA=np.log10(Variable_y_PAPA)  
                Variable_y_NABE=np.log10(Variable_y_NABE)  

            Variable_x_PAPA_reshaped = Variable_x_PAPA.reshape(( yrs_n ,12,Variable_x_PAPA.shape[1],Variable_x_PAPA.shape[2])) # Monthly means of SpCO2 over the time period  
            Variable_y_PAPA_reshaped = Variable_y_PAPA.reshape(( yrs_n ,12,Variable_y_PAPA.shape[1],Variable_y_PAPA.shape[2])) # Monthly means of SpCO2 over the time period 
  
            Variable_x_NABE_reshaped = Variable_x_NABE.reshape(( yrs_n ,12,Variable_x_NABE.shape[1],Variable_x_NABE.shape[2])) # Monthly means of SpCO2 over the time period  
            Variable_y_NABE_reshaped = Variable_y_NABE.reshape(( yrs_n ,12,Variable_y_NABE.shape[1],Variable_y_NABE.shape[2])) # Monthly means of SpCO2 over the time period 
                
            Variable_y_PAPA_reshaped=np.nanmean(np.nanmean(np.nanmean(Variable_y_PAPA_reshaped,axis=0),axis=1) ,axis=1)
            Variable_x_PAPA_reshaped=np.nanmean(np.nanmean(np.nanmean(Variable_x_PAPA_reshaped,axis=0),axis=1) ,axis=1)
            
            Variable_y_NABE_reshaped=np.nanmean(np.nanmean(np.nanmean(Variable_y_NABE_reshaped,axis=0),axis=1) ,axis=1)
            Variable_x_NABE_reshaped=np.nanmean(np.nanmean(np.nanmean(Variable_x_NABE_reshaped,axis=0),axis=1) ,axis=1) 
            
            ax.plot(Variable_x_PAPA_reshaped, Variable_y_PAPA_reshaped, color='k', lw=1)    
            for ii in range(Variable_y_PAPA_reshaped.shape[0]):
                plt.scatter(Variable_x_PAPA_reshaped[ii], Variable_y_PAPA_reshaped[ii], c= colors_months[ii] , s=50, marker='o', label=Time_months[ii] )
            #ax.text(Variable_x_PAPA_reshaped[4], Variable_y_PAPA_reshaped[4], "N.P.", color='k', fontsize=14)

            ax.plot(Variable_x_NABE_reshaped, Variable_y_NABE_reshaped, color='k', lw=1, linestyle='--')    
            for ii in range(Variable_y_NABE_reshaped.shape[0]):
                plt.scatter(Variable_x_NABE_reshaped[ii], Variable_y_NABE_reshaped[ii], c= colors_months[ii] , s=50, marker='d', label=Time_months[ii] )
            #ax.text(Variable_x_NABE_reshaped[4], Variable_y_NABE_reshaped[4], "N.A.", color='k', fontsize=14) 
            
#            ax.set_xscale('log')
#            ax.set_yscale('log')
            
            plt.title(GCM, fontsize=14)
    
            if M_i == (n_c*(n_r-1)): # Adds longitude ranges only to the last subplots that appear at the bottom of plot
                plt.xlabel(var_names_x[V_i][:-10],fontsize=14)
                plt.ylabel(var_names_y[V_i][:-10],fontsize=14)
            elif M_i==0 or M_i==n_c or M_i==n_c*2 or M_i==n_c*3 or M_i==n_c*4 or M_i==n_c*5 or M_i==n_c*6 or M_i==n_c*7 or M_i==n_c*8:
                plt.ylabel(var_names_y[V_i][:-10],fontsize=14)
            elif M_i >= n_t-n_c and M_i != (n_c*(n_r-1)): # Adds longitude ranges only to the last subplots that appear at the bottom of plot            
                plt.xlabel(var_names_x[V_i][:-10],fontsize=14)

    plt.subplots_adjust(left=0.18, bottom=0.06, right=0.83, top=0.92, hspace=0.35, wspace=0.25) # the amount of height/width reserved for space between subplots
    #plt.legend(shadow=True, loc='right', bbox_to_anchor=(5.5, 3), markerscale=16 , prop={'size': 16})
    plt.suptitle(P_title, fontsize=20)    
    plt.show()
    mng = plt.get_current_fig_manager()
    mng.window.showMaximized() # Maximizes the plot window to save figures in full
    
    if d_trend=='yes':
        fig.savefig(dir_figs+'AllModels_ScatterPlots_ave_PAPA_NABE_detrended_'+save_names[V_i]+'_2.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
    else:
        fig.savefig(dir_figs+'AllModels_ScatterPlots_ave_PAPA_NABE_detrended_'+save_names[V_i]+'_2.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


###   N PACIFIC and N ATLANTIC plotted together - normal scale
#### Plots 1 #####
var_names_y=[ 'Phyc_allmonths',         'Phyc_allmonths',            'Phyc_allmonths',     'Phyc_allmonths',   'Zooc_allmonths',    'Zooc_allmonths',        'Zooc_allmonths',    'Zooc_allmonths',  'Phyc_allmonths',   'Phyc_allmonths',     'Zooc_allmonths',  'Zoo_by_Phyto_allmonths',      'Zoo_by_Phyto_allmonths',  'Zoo_by_Phyto_allmonths', 'Fe_allmonths'       ]
var_names_x=[ 'Fe_allmonths',           'NO3_allmonths',             'RSNTDS_allmonths',   'MLD_allmonths',    'Fe_allmonths',      'NO3_allmonths',         'RSNTDS_allmonths',  'MLD_allmonths',   'Zooc_allmonths',   'PPint_allmonths',    'Phyc_allmonths',  'NO3_allmonths',               'Fe_allmonths',            'RSNTDS_allmonths',       'NO3_allmonths'     ]
plot_titles=[ 'biomass vs. Iron (Fe)',  'biomass vs. Nitrate (NO3)', 'biomass vs. light',  'biomass vs. MLD',  'Zoo vs. Iron (Fe)', 'Zoo vs. Nitrate (NO3)', 'Zoo vs. light',     'Zoo vs. MLD',     'biomass vs. Zoo',  'biomass vs. intPP',  'Zoo vs. biomass', 'Zoo/phyto vs. Nitrate (NO3)', 'Zoo/phyto vs. Iron (Fe)', 'Zoo/phyto vs. light',    'Iron (Fe) vs. Nitrate (NO3)'      ]
save_names= [ 'Phyc_vs_Fe',             'Phyc_vs_NO3',               'Phyc_vs_RSNTDS',     'Phyc_vs_MLD',      'Zooc_vs_Fe',        'Zooc_vs_NO3',           'Zooc_vs_RSNTDS',    'Zooc_vs_MLD',     'Phyc_vs_Zooc',     'Phyc_vs_PPint',      'Zooc_vs_Phyc',    'ZoobyPhyc_vs_NO3',            'ZoobyPhyc_vs_Fe',         'ZoobyPhyc_vs_RSNTDS',    'Fe_vs_NO3'       ]  

file_names_y=[ 'bio',                   'bio',                       'bio',                'bio',              'bio',               'bio',                   'bio' ,              'bio',              'bio',             'bio',                'bio',             'bio',                         'bio',                      'bio',                   'bio'            ]
file_names_x=[ 'bio',                   'bio',                       'physics',            'physics',          'bio',               'bio',                   'physics',           'physics',          'bio',             'bio',                'bio',             'bio',                         'bio',                      'physics',               'bio'            ]

plot_ranges_y=[[-0.001,0.02],  [-0.001,0.02]]
plot_ranges_x=[[0,4E-6],  [-0.001, 1.6E-2]]

#### Plots 2 #####
var_names_y=[ 'PPint_allmonths',     'PPint_allmonths',         'PPint_allmonths',  'PPint_allmonths', 'EPC100_allmonths',        'EPC100_allmonths',          'EPC100_allmonths',    'EPC100_allmonths',  'PPint_allmonths', 'EPC100_allmonths', 'PP_by_Phyto_allmonths',       'PP_by_Phyto_allmonths',  'PP_by_Phyto_allmonths',  'PP_by_Phyto_allmonths',  'LnPP_by_LnPhyto_allmonths'         ]
var_names_x=[ 'Fe_allmonths',        'NO3_allmonths',           'RSNTDS_allmonths', 'MLD_allmonths',   'Fe_allmonths',            'NO3_allmonths',             'RSNTDS_allmonths',    'MLD_allmonths',     'Zooc_allmonths',  'Zooc_allmonths',   'NO3_allmonths',               'Fe_allmonths',           'RSNTDS_allmonths',       'Zooc_allmonths',         'Zooc_allmonths'               ]
plot_titles=[ 'intPP vs. Iron (Fe)', 'intPP vs. Nitrate (NO3)', 'intPP vs. light',  'intPP vs. MLD',   'EPC 100m vs. Iron (Fe)',  'EPC 100m vs. Nitrate (NO3)','EPC 100m vs. light',  'EPC 100m vs. MLD',  'intPP vs. Zoo',   'EPC 100m vs. Zoo', 'PP/phyto vs. Nitrate (NO3)',  'PP/phyto vs. Iron (Fe)', 'PP/phyto vs. light',     'PP/phyto vs. Zoo',       'Ln(PP)/Ln(phyto) vs. Zoo'   ]
save_names= [ 'PPint_vs_Fe',         'PPint_vs_NO3',            'PPint_vs_RSNTDS',  'PPint_vs_MLD',    'EPC100_vs_Fe',            'EPC100_vs_NO3' ,            'EPC100_vs_RSNTDS',    'EPC100_vs_MLD',     'PPint_vs_Zooc',   'EPC100_vs_Zooc',   'PPbyPhyc_vs_NO3',             'PPbyPhyc_vs_Fe',         'PPbyPhyc_vs_RSNTDS',     'PPbyPhyc_vs_Zooc',       'PP_by_Phyto_Ln_vs_Zooc'                   ]  

file_names_y=[ 'bio',                'bio',                     'bio',              'bio',             'bio',                     'bio',                       'bio',                 'bio',               'bio',             'bio',              'bio',                         'bio',                    'bio',                    'bio',                    'bio'            ]
file_names_x=[ 'bio',                'bio',                     'physics',          'physics',         'bio',                     'bio',                       'physics',             'physics',           'bio',             'bio',              'bio',                         'bio',                    'physics',                'bio',                    'bio'            ]

#### Plots 3 #####
var_names_y=[ 'Large Phyto_allmonths',       'Large Phyto_allmonths',        '% Large Phyto_allmonths',        'Phydiat_allmonths',        'LnPP_by_LnPhyto_allmonths',           'LnPP_by_LnPhyto_allmonths',           'Chl_allmonths',           'Chl_allmonths',        'Chl_allmonths',           'Chl_allmonths'                           ]
var_names_x=[ 'Phyc_allmonths',              'Phypico_allmonths',            'Phyc_allmonths',                 'Phyc_allmonths',           'NO3_allmonths',                       '% Large Phyto_allmonths',             'Phyc_allmonths',          'Fe_allmonths',         'NO3_allmonths',           'RSNTDS_allmonths'                                           ]
plot_titles=[ 'Large phyto vs. Total phyto', 'Large phyto vs. Small phyto',  '% Large phyto vs. Total phyto',  'Diatoms vs. Total phyto',  'Ln(PP)/Ln(phyto) vs. Nitrate (NO3)',  'Ln(PP)/Ln(phyto) vs. % Large phyto',  'Chlorophyll vs. biomass', 'Chlorophyll vs. Iron', 'Chlorophyll vs. Nitrate', 'Chlorophyll vs. light'                      ]
save_names= [ 'Phyto_Large_vs_Phyc',         'Phyto_Large_vs_Phypico',       'Phyto_Large_prcnt_vs_Phyc',      'Phydiat_vs_Phyc',          'PP_by_Phyto_Ln_vs_NO3' ,              'PP_by_Phyto_Ln_vs_Phyto_Large_prcnt', 'Chl_vs_Phyc',             'Chl_vs_Fe',            'Chl_vs_NO3',              'Chl_vs_RSNTDS'                              ]  

file_names_y=[ 'bio',                        'bio',                          'bio',                            'bio',                      'bio',                                 'bio',                                 'bio',                     'bio',                  'bio',                     'bio'                                              ]
file_names_x=[ 'bio',                        'bio',                          'bio',                            'bio',                      'bio',                                 'bio',                                 'bio',                     'bio',                  'bio',                     'physics'                                               ]


for V_i in range(n_t): # V_i=0   V_i=1    V_i=2    V_i=3    V_i=4    V_i=5    V_i=6    V_i=7    V_i=8    V_i=9    V_i=10    V_i=11    V_i=12    V_i=13   V_i=14

    if d_trend=='yes':
        P_title=plot_titles[V_i]+' - North Pacific (circles) and North Atlantic (diamonds) - '+str(start_date_plt)+'-'+str(end_date_plt)+' - detrended'
    else:
        P_title=plot_titles[V_i]+' - North Pacific (circles) and North Atlantic (diamonds) - '+str(start_date_plt)+'-'+str(end_date_plt)
    
    fig=plt.figure()
    for M_i in range(len(GCM_Names)): # M_i=0     # GCM=GCM_Names[0]

        GCM=GCM_Names[M_i]
        print (GCM) 
        
        ax = fig.add_subplot(n_r,n_c,M_i+1)      
        
        if var_names_y[V_i]=='Zoo_by_Phyto_allmonths' :
            filename_in_y = (dir_data_in + 'AllResults_'+GCM+'_'+'bio'+'_1991_2010.out') # Directory to save processed data
            filename_in_x = (dir_data_in + 'AllResults_'+GCM+'_'+file_names_x[V_i]+'_1991_2010.out') # Directory to save processed data
            my_shelf_y = shelve.open(filename_in_y)
            my_shelf_x = shelve.open(filename_in_x)            
            
            if 'Zooc_allmonths' in my_shelf_y and 'Phyc_allmonths' in my_shelf_y and var_names_x[V_i] in my_shelf_x: # If that variable exists in the saved data                       
                globals()['Variable_1']=my_shelf_y['Zooc_allmonths'] 
                globals()['Variable_2']=my_shelf_y['Phyc_allmonths']
                globals()['Variable_x']=my_shelf_x[var_names_x[V_i]]
                Variable_y=Variable_1/Variable_2
                Variable_y[ np.abs(Variable_y)>1e19 ]=nan; Variable_x[ np.abs(Variable_x)>1e19 ]=nan  # Replacing missing values with NaN 
            else:
                Variable_y=empty((yrs_n*12,lat_n_regrid,lon_n_regrid)) * nan
                Variable_x=empty((yrs_n*12,lat_n_regrid,lon_n_regrid)) * nan

        elif var_names_y[V_i]=='PP_by_Phyto_allmonths' or  var_names_y[V_i]=='LnPP_by_LnPhyto_allmonths' :
            
            if var_names_x[V_i]=='% Large Phyto_allmonths':
                filename_in_y = (dir_data_in + 'AllResults_'+GCM+'_'+'bio'+'_1991_2010.out') # Directory to save processed data
                my_shelf_y = shelve.open(filename_in_y)         
                
                if 'PPint_allmonths' in my_shelf_y and 'Phyc_allmonths' in my_shelf_y and 'Phypico_allmonths' in my_shelf_y: # If that variable exists in the saved data                       
                    globals()['Variable_1']=my_shelf_y['PPint_allmonths'] 
                    globals()['Variable_2']=my_shelf_y['Phyc_allmonths']
                    Variable_y=np.log(Variable_1)/np.log(Variable_2)
                    globals()['Variable_1']=my_shelf_y['Phyc_allmonths'] 
                    globals()['Variable_2']=my_shelf_y['Phypico_allmonths']                
                    Variable_x= ((Variable_1 - Variable_2) / Variable_1) * 100
                    
                    Variable_y[ np.abs(Variable_y)>1e19 ]=nan; Variable_x[ np.abs(Variable_x)>1e19 ]=nan  # Replacing missing values with NaN 
                    
                elif 'PPint_allmonths' in my_shelf_y and 'Phyc_allmonths' in my_shelf_y and 'Phydiat_allmonths' in my_shelf_y: # If that variable exists in the saved data 
                    globals()['Variable_1']=my_shelf_y['PPint_allmonths'] 
                    globals()['Variable_2']=my_shelf_y['Phyc_allmonths']
                    Variable_y=np.log(Variable_1)/np.log(Variable_2)
                    globals()['Variable_1']=my_shelf_y['Phyc_allmonths'] 
                    globals()['Variable_2']=my_shelf_y['Phydiat_allmonths']                
                    Variable_x= (Variable_2 / Variable_1) * 100
                    
                    Variable_y[ np.abs(Variable_y)>1e19 ]=nan; Variable_x[ np.abs(Variable_x)>1e19 ]=nan  # Replacing missing values with NaN                     
                    
                else:
                    Variable_y=empty((yrs_n*12,lat_n_regrid,lon_n_regrid)) * nan
                    Variable_x=empty((yrs_n*12,lat_n_regrid,lon_n_regrid)) * nan                

            else:
                filename_in_y = (dir_data_in + 'AllResults_'+GCM+'_'+'bio'+'_1991_2010.out') # Directory to save processed data
                filename_in_x = (dir_data_in + 'AllResults_'+GCM+'_'+file_names_x[V_i]+'_1991_2010.out') # Directory to save processed data
                my_shelf_y = shelve.open(filename_in_y)
                my_shelf_x = shelve.open(filename_in_x)            
                
                if 'PPint_allmonths' in my_shelf_y and 'Phyc_allmonths' in my_shelf_y and var_names_x[V_i] in my_shelf_x: # If that variable exists in the saved data                       
                    globals()['Variable_1']=my_shelf_y['PPint_allmonths'] 
                    globals()['Variable_2']=my_shelf_y['Phyc_allmonths']
                    globals()['Variable_x']=my_shelf_x[var_names_x[V_i]]
                    
                    if var_names_y[V_i]=='LnPP_by_LnPhyto_allmonths':
                        Variable_y=np.log(Variable_1)/np.log(Variable_2)
                    else:
                        Variable_y=Variable_1/Variable_2
                    
                    Variable_y[ np.abs(Variable_y)>1e19 ]=nan; Variable_x[ np.abs(Variable_x)>1e19 ]=nan  # Replacing missing values with NaN 
                else:
                    Variable_y=empty((yrs_n*12,lat_n_regrid,lon_n_regrid)) * nan
                    Variable_x=empty((yrs_n*12,lat_n_regrid,lon_n_regrid)) * nan

        elif var_names_y[V_i]=='Large Phyto_allmonths' or  var_names_y[V_i]=='% Large Phyto_allmonths' :
            filename_in_y = (dir_data_in + 'AllResults_'+GCM+'_'+'bio'+'_1991_2010.out') # Directory to save processed data
            filename_in_x = (dir_data_in + 'AllResults_'+GCM+'_'+file_names_x[V_i]+'_1991_2010.out') # Directory to save processed data
            my_shelf_y = shelve.open(filename_in_y)
            my_shelf_x = shelve.open(filename_in_x)            
            
            if 'Phyc_allmonths' in my_shelf_y and 'Phypico_allmonths' in my_shelf_y and var_names_x[V_i] in my_shelf_x: # If that variable exists in the saved data                       
                globals()['Variable_1']=my_shelf_y['Phyc_allmonths'] 
                globals()['Variable_2']=my_shelf_y['Phypico_allmonths']
                globals()['Variable_x']=my_shelf_x[var_names_x[V_i]]
                
                Variable_1[ np.abs(Variable_1)>1e19 ]=nan; Variable_2[ np.abs(Variable_2)>1e19 ]=nan; 
                
                if var_names_y[V_i]=='Large Phyto_allmonths':
                    Variable_y=Variable_1 - Variable_2
                elif var_names_y[V_i]=='% Large Phyto_allmonths':
                    Variable_y= ((Variable_1 - Variable_2) / Variable_1) * 100
                
                Variable_y[ np.abs(Variable_y)>1e19 ]=nan; Variable_x[ np.abs(Variable_x)>1e19 ]=nan  # Replacing missing values with NaN 

            elif 'Phyc_allmonths' in my_shelf_y and 'Phydiat_allmonths' in my_shelf_y and var_names_x[V_i] in my_shelf_x: # If that variable exists in the saved data                       
                globals()['Variable_1']=my_shelf_y['Phyc_allmonths'] 
                globals()['Variable_2']=my_shelf_y['Phydiat_allmonths']
                globals()['Variable_x']=my_shelf_x[var_names_x[V_i]]
                
                Variable_1[ np.abs(Variable_1)>1e19 ]=nan; Variable_2[ np.abs(Variable_2)>1e19 ]=nan; 
                
                if var_names_y[V_i]=='Large Phyto_allmonths':
                    Variable_y=Variable_2
                elif var_names_y[V_i]=='% Large Phyto_allmonths':
                    Variable_y= (Variable_2 / Variable_1) * 100
                
                Variable_y[ np.abs(Variable_y)>1e19 ]=nan; Variable_x[ np.abs(Variable_x)>1e19 ]=nan  # Replacing missing values with NaN 

            else:
                Variable_y=empty((yrs_n*12,lat_n_regrid,lon_n_regrid)) * nan
                Variable_x=empty((yrs_n*12,lat_n_regrid,lon_n_regrid)) * nan

        else:

            filename_in_y = (dir_data_in + 'AllResults_'+GCM+'_'+file_names_y[V_i]+'_1991_2010.out') # Directory to save processed data
            filename_in_x = (dir_data_in + 'AllResults_'+GCM+'_'+file_names_x[V_i]+'_1991_2010.out') # Directory to save processed data
            my_shelf_y = shelve.open(filename_in_y)
            my_shelf_x = shelve.open(filename_in_x)
     
            if var_names_y[V_i] in my_shelf_y and var_names_x[V_i] in my_shelf_x: # If that variable exists in the saved data  
                
                globals()['Variable_y']=my_shelf_y[var_names_y[V_i]]
                globals()['Variable_x']=my_shelf_x[var_names_x[V_i]]                        
                Variable_y[ np.abs(Variable_y)>1e19 ]=nan; Variable_x[ np.abs(Variable_x)>1e19 ]=nan  # Replacing missing values with NaN 
            else:
                Variable_y=empty((yrs_n*12,lat_n_regrid,lon_n_regrid)) * nan
                Variable_x=empty((yrs_n*12,lat_n_regrid,lon_n_regrid)) * nan                

    
        Variable_x_PAPA = Variable_x[:,135:140,210:220]
        Variable_y_PAPA = Variable_y[:,135:140,210:220]

        Variable_x_NABE = Variable_x[:,135:140,325:335]
        Variable_y_NABE = Variable_y[:,135:140,325:335]
                
        if d_trend=='yes': # If it is requested to detrend the data before correlations
            Variable_x_PAPA= func_detrend_3d(Variable_x_PAPA)
            Variable_y_PAPA= func_detrend_3d(Variable_y_PAPA)                
            Variable_x_NABE= func_detrend_3d(Variable_x_NABE)
            Variable_y_NABE= func_detrend_3d(Variable_y_NABE)

        Variable_x_PAPA_reshaped = Variable_x_PAPA.reshape(( yrs_n ,12,Variable_x_PAPA.shape[1],Variable_x_PAPA.shape[2])) # Monthly means of SpCO2 over the time period  
        Variable_y_PAPA_reshaped = Variable_y_PAPA.reshape(( yrs_n ,12,Variable_y_PAPA.shape[1],Variable_y_PAPA.shape[2])) # Monthly means of SpCO2 over the time period 
  
        Variable_x_NABE_reshaped = Variable_x_NABE.reshape(( yrs_n ,12,Variable_x_NABE.shape[1],Variable_x_NABE.shape[2])) # Monthly means of SpCO2 over the time period  
        Variable_y_NABE_reshaped = Variable_y_NABE.reshape(( yrs_n ,12,Variable_y_NABE.shape[1],Variable_y_NABE.shape[2])) # Monthly means of SpCO2 over the time period 
            
        Variable_y_PAPA_reshaped=np.nanmean(np.nanmean(np.nanmean(Variable_y_PAPA_reshaped,axis=0),axis=1) ,axis=1)
        Variable_x_PAPA_reshaped=np.nanmean(np.nanmean(np.nanmean(Variable_x_PAPA_reshaped,axis=0),axis=1) ,axis=1)
        
        Variable_y_NABE_reshaped=np.nanmean(np.nanmean(np.nanmean(Variable_y_NABE_reshaped,axis=0),axis=1) ,axis=1)
        Variable_x_NABE_reshaped=np.nanmean(np.nanmean(np.nanmean(Variable_x_NABE_reshaped,axis=0),axis=1) ,axis=1) 
     
        ax.plot(Variable_x_PAPA_reshaped, Variable_y_PAPA_reshaped, color='k', lw=1)    
        for ii in range(Variable_y_PAPA_reshaped.shape[0]):
            plt.scatter(Variable_x_PAPA_reshaped[ii], Variable_y_PAPA_reshaped[ii], c= colors_months[ii] , s=50, marker='o', label=Time_months[ii] )
        #ax.text(Variable_x_PAPA_reshaped[4], Variable_y_PAPA_reshaped[4], "N.P.", color='k', fontsize=14)

        ax.plot(Variable_x_NABE_reshaped, Variable_y_NABE_reshaped, color='k', lw=1, linestyle='--')    
        for ii in range(Variable_y_NABE_reshaped.shape[0]):
            plt.scatter(Variable_x_NABE_reshaped[ii], Variable_y_NABE_reshaped[ii], c= colors_months[ii] , s=50, marker='d', label=Time_months[ii] )
        #ax.text(Variable_x_NABE_reshaped[4], Variable_y_NABE_reshaped[4], "N.A.", color='k', fontsize=14) 
        
        plt.title(GCM, fontsize=14)
        #plt.xlim(plot_ranges_x[V_i][0],plot_ranges_x[V_i][1]) ; plt.ylim(plot_ranges_y[V_i][0],plot_ranges_y[V_i][1])
        ##plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        ##plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        if  np.nanmax( np.abs(Variable_x_NABE_reshaped) ) < 0.001 :
            plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))  
        if  np.nanmax( np.abs(Variable_y_NABE_reshaped) ) < 0.001 :
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))  

        if M_i == (n_c*(n_r-1)): # Adds longitude ranges only to the last subplots that appear at the bottom of plot
            plt.xlabel(var_names_x[V_i][:-10],fontsize=14)
            plt.ylabel(var_names_y[V_i][:-10],fontsize=14)
        elif M_i==0 or M_i==n_c or M_i==n_c*2 or M_i==n_c*3 or M_i==n_c*4 or M_i==n_c*5 or M_i==n_c*6 or M_i==n_c*7 or M_i==n_c*8:
            plt.ylabel(var_names_y[V_i][:-10],fontsize=14)
        elif M_i >= n_t-n_c and M_i != (n_c*(n_r-1)): # Adds longitude ranges only to the last subplots that appear at the bottom of plot            
            plt.xlabel(var_names_x[V_i][:-10],fontsize=14)

    plt.subplots_adjust(left=0.18, bottom=0.06, right=0.83, top=0.92, hspace=0.35, wspace=0.25) # the amount of height/width reserved for space between subplots
    #plt.legend(shadow=True, loc='right', bbox_to_anchor=(5.5, 3), markerscale=16 , prop={'size': 16})
    plt.suptitle(P_title, fontsize=20)    
    plt.show()
    mng = plt.get_current_fig_manager()
    mng.window.showMaximized() # Maximizes the plot window to save figures in full
    
    if d_trend=='yes':
        fig.savefig(dir_figs+'AllModels_ScatterPlots_ave_PAPA_NABE_detrended_'+save_names[V_i]+'_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
    else:
        fig.savefig(dir_figs+'AllModels_ScatterPlots_ave_PAPA_NABE_detrended_'+save_names[V_i]+'_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')



######################################################################
###   Max/Min of VarX vs VarY scatter plots Multi-GCM - Subplots  ####
######################################################################
GCM_Names = ['CESM1-BGC' ,'GFDL-ESM2G','GFDL-ESM2M','IPSL-CM5A-MR','IPSL-CM5A-LR','MRI-ESM1','GISS-E2-H-CC','GISS-E2-R-CC','HadGEM2-ES','CanESM2','NorESM1-ME','MPI-ESM-MR','MPI-ESM-LR']
Colors    =  [[1, 1, 0],   [1, 0.8, 0],  [1, 0.6, 0], [1, 0.4, 0],   [1, 0, 0], [1, 0.8, 1],   [1, 0.5, 1],   [1, 0, 1],   [0, 0, 0.5], [0, 0, 1], [0, 0.5, 1],[0, 0.5, 0],  [0, 1, 0]]

plot_ranges_y_max=[[-0.001,0.02],     [-0.001,0.02],    [-0.001,0.02],        [-0.001,0.02],      [-0.001,0.02],     [-0.001,0.02],       [-0.001,0.02],      [-0.001,0.02],      [-0.01,2.5],              [-0.01,2.5],               [-0.01,2.5],              [-1E-7,0.0000015],  [-1E-7,0.0000015],   [-0.00001,0.0005],        [-0.00001,0.0005],         [-0.00001,0.0005]           ]
plot_ranges_x_max=[[-1E-7,0.000004],  [-0.001,0.02],    [100,230],            [1,1000],           [-0.001,0.02],     [100,230],           [-0.001,0.02],      [-1E-7,0.0000015],  [-1E-7,0.000004],         [-0.001,0.02],             [100,230],                [-0.001,0.02],      [100,230],           [-0.001,0.02],            [100,230],                 [-0.001,0.02]      ]
var_names_y=[ 'Phyc_allmonths',  'Phyc_allmonths',  'Phyc_allmonths',     'Phyc_allmonths',   'Zooc_allmonths',  'Zooc_allmonths',    'Phyc_allmonths',   'Phyc_allmonths',  'Zoo_by_Phyto_allmonths',  'Zoo_by_Phyto_allmonths',  'Zoo_by_Phyto_allmonths', 'PPint_allmonths', 'PPint_allmonths',   'PP_by_Phyto_allmonths',   'PP_by_Phyto_allmonths',   'PP_by_Phyto_allmonths'       ]
var_names_x=[ 'Fe_allmonths',    'NO3_allmonths',   'RSNTDS_allmonths',   'MLD_allmonths',    'NO3_allmonths',   'RSNTDS_allmonths',  'Zooc_allmonths',   'PPint_allmonths', 'Fe_allmonths',            'NO3_allmonths',           'RSNTDS_allmonths',       'NO3_allmonths',   'RSNTDS_allmonths',  'NO3_allmonths',           'RSNTDS_allmonths',        'Zooc_allmonths'     ]
 
file_names_y=[ 'bio',            'bio',             'bio',                'bio',              'bio',             'bio' ,              'bio',              'bio',             'bio',                     'bio',                      'bio',                   'bio',             'bio',               'bio',                     'bio',                     'bio'            ]
file_names_x=[ 'bio',            'bio',             'physics',            'physics',          'bio',             'physics',           'bio',              'bio',             'bio',                     'bio',                      'physics',               'bio',             'physics',           'bio',                     'physics',                 'bio'          ]

plot_ranges_y_min=[[-0.0001,0.0025],  [-0.0001,0.0025],  [-0.0001,0.0025],     [-0.0001,0.0025],   [-0.0001,0.0025],  [-0.0001,0.0025],    [-0.0001,0.0025],   [-0.0001,0.0025],   [-0.01,1.25],           [-0.01,1.25],              [-0.01,1.25],            [-1E-8,2E-7],       [-1E-8,2E-7],      [-0.00001,0.0002],        [-0.00001,0.0002],         [-0.00001,0.0002]           ]
plot_ranges_x_min=[[-1E-7,0.000004],  [-0.0001,0.006],   [20,50],              [1,100],            [-0.0001,0.006],   [20,50],             [-0.0001,0.0025],   [-1E-8,2E-7],  [-1E-7,0.000004],            [-0.0001,0.006],           [20,50],                 [-0.0001,0.006],    [20,50],           [-0.0001,0.006],          [20,50],                    [-0.0001,0.0025]      ]


d_trend='yes'
n_r=4 ; n_c=4 ;
n_t=len(var_names_y)
max_min_ind='min'

fig=plt.figure()
for V_i in range(n_t): # V_i=0   V_i=1    V_i=2    V_i=3    V_i=4    V_i=5    V_i=6    V_i=7    V_i=8    V_i=9    V_i=10    V_i=11    V_i=12    V_i=13

    ax = fig.add_subplot(n_r,n_c,V_i+1) 
    
    if max_min_ind=='max':
        plot_ranges_y=plot_ranges_y_max; plot_ranges_x=plot_ranges_x_max
    elif max_min_ind=='min':
        plot_ranges_y=plot_ranges_y_min; plot_ranges_x=plot_ranges_x_min

    for M_i in range(len(GCM_Names)): # M_i=0     # GCM=GCM_Names[0]

        GCM=GCM_Names[M_i]
        print (GCM) 

        if var_names_y[V_i]=='Zoo_by_Phyto_allmonths' :
            filename_in_y = (dir_data_in + 'AllResults_'+GCM+'_'+'bio'+'_1991_2010.out') # Directory to save processed data
            filename_in_x = (dir_data_in + 'AllResults_'+GCM+'_'+file_names_x[V_i]+'_1991_2010.out') # Directory to save processed data
            my_shelf_y = shelve.open(filename_in_y)
            my_shelf_x = shelve.open(filename_in_x)            
            
            if 'Zooc_allmonths' in my_shelf_y and 'Phyc_allmonths' in my_shelf_y and var_names_x[V_i] in my_shelf_x: # If that variable exists in the saved data                       
                globals()['Variable_1']=my_shelf_y['Zooc_allmonths'] 
                globals()['Variable_2']=my_shelf_y['Phyc_allmonths']
                globals()['Variable_x']=my_shelf_x[var_names_x[V_i]]
                Variable_y=Variable_1/Variable_2
                Variable_y[ np.abs(Variable_y)>1e19 ]=nan; Variable_x[ np.abs(Variable_x)>1e19 ]=nan  # Replacing missing values with NaN 
            else:
                Variable_y=empty((yrs_n*12,lat_n_regrid,lon_n_regrid)) * nan
                Variable_x=empty((yrs_n*12,lat_n_regrid,lon_n_regrid)) * nan

        elif var_names_y[V_i]=='PP_by_Phyto_allmonths' or  var_names_y[V_i]=='LnPP_by_LnPhyto_allmonths' :
            filename_in_y = (dir_data_in + 'AllResults_'+GCM+'_'+'bio'+'_1991_2010.out') # Directory to save processed data
            filename_in_x = (dir_data_in + 'AllResults_'+GCM+'_'+file_names_x[V_i]+'_1991_2010.out') # Directory to save processed data
            my_shelf_y = shelve.open(filename_in_y)
            my_shelf_x = shelve.open(filename_in_x)            
            
            if 'PPint_allmonths' in my_shelf_y and 'Phyc_allmonths' in my_shelf_y and var_names_x[V_i] in my_shelf_x: # If that variable exists in the saved data                       
                globals()['Variable_1']=my_shelf_y['PPint_allmonths'] 
                globals()['Variable_2']=my_shelf_y['Phyc_allmonths']
                globals()['Variable_x']=my_shelf_x[var_names_x[V_i]]
                
                if var_names_y[V_i]=='LnPP_by_LnPhyto_allmonths':
                    Variable_y=np.log(Variable_1)/np.log(Variable_2)
                else:
                    Variable_y=Variable_1/Variable_2
                
                Variable_y[ np.abs(Variable_y)>1e19 ]=nan; Variable_x[ np.abs(Variable_x)>1e19 ]=nan  # Replacing missing values with NaN 
            else:
                Variable_y=empty((yrs_n*12,lat_n_regrid,lon_n_regrid)) * nan
                Variable_x=empty((yrs_n*12,lat_n_regrid,lon_n_regrid)) * nan

        else:

            filename_in_y = (dir_data_in + 'AllResults_'+GCM+'_'+file_names_y[V_i]+'_1991_2010.out') # Directory to save processed data
            filename_in_x = (dir_data_in + 'AllResults_'+GCM+'_'+file_names_x[V_i]+'_1991_2010.out') # Directory to save processed data
            my_shelf_y = shelve.open(filename_in_y)
            my_shelf_x = shelve.open(filename_in_x)
     
            if var_names_y[V_i] in my_shelf_y and var_names_x[V_i] in my_shelf_x: # If that variable exists in the saved data                  
                globals()['Variable_y']=my_shelf_y[var_names_y[V_i]]
                globals()['Variable_x']=my_shelf_x[var_names_x[V_i]]                        
                Variable_y[ np.abs(Variable_y)>1e19 ]=nan; Variable_x[ np.abs(Variable_x)>1e19 ]=nan  # Replacing missing values with NaN 
            else:
                Variable_y=empty((yrs_n*12,lat_n_regrid,lon_n_regrid)) * nan
                Variable_x=empty((yrs_n*12,lat_n_regrid,lon_n_regrid)) * nan  

        Variable_x_PAPA = Variable_x[:,135:140,210:220]
        Variable_y_PAPA = Variable_y[:,135:140,210:220]

        Variable_x_NABE = Variable_x[:,135:140,325:335]
        Variable_y_NABE = Variable_y[:,135:140,325:335]
                
        if d_trend=='yes': # If it is requested to detrend the data before correlations
            Variable_x_PAPA= func_detrend_3d(Variable_x_PAPA)
            Variable_y_PAPA= func_detrend_3d(Variable_y_PAPA)                
            Variable_x_NABE= func_detrend_3d(Variable_x_NABE)
            Variable_y_NABE= func_detrend_3d(Variable_y_NABE)

        Variable_x_PAPA_reshaped = Variable_x_PAPA.reshape(( yrs_n ,12,Variable_x_PAPA.shape[1],Variable_x_PAPA.shape[2])) # Monthly means of SpCO2 over the time period  
        Variable_y_PAPA_reshaped = Variable_y_PAPA.reshape(( yrs_n ,12,Variable_y_PAPA.shape[1],Variable_y_PAPA.shape[2])) # Monthly means of SpCO2 over the time period 
  
        Variable_x_NABE_reshaped = Variable_x_NABE.reshape(( yrs_n ,12,Variable_x_NABE.shape[1],Variable_x_NABE.shape[2])) # Monthly means of SpCO2 over the time period  
        Variable_y_NABE_reshaped = Variable_y_NABE.reshape(( yrs_n ,12,Variable_y_NABE.shape[1],Variable_y_NABE.shape[2])) # Monthly means of SpCO2 over the time period 
        
        Variable_y_PAPA_reshaped=np.nanmean(np.nanmean(np.nanmean(Variable_y_PAPA_reshaped,axis=0),axis=1) ,axis=1)
        Variable_x_PAPA_reshaped=np.nanmean(np.nanmean(np.nanmean(Variable_x_PAPA_reshaped,axis=0),axis=1) ,axis=1)
        
        Variable_y_NABE_reshaped=np.nanmean(np.nanmean(np.nanmean(Variable_y_NABE_reshaped,axis=0),axis=1) ,axis=1)
        Variable_x_NABE_reshaped=np.nanmean(np.nanmean(np.nanmean(Variable_x_NABE_reshaped,axis=0),axis=1) ,axis=1) 
        
        if max_min_ind=='max':
 
            vec_new=np.concatenate((Variable_x_PAPA_reshaped, Variable_x_PAPA_reshaped*1.000000000001,Variable_x_PAPA_reshaped), axis=0) # put vector 3 times in a row and nultiply the middle one by a small value so the maximum value alwasy found in the middle matrix
            ii_maxmin=np.argmax(vec_new); Var_plot_x_PAPA=( ( vec_new[ii_maxmin-1]+vec_new[ii_maxmin]+vec_new[ii_maxmin+1]) / 3 ) / 1.000000000001
            if var_names_x[V_i]=='RSNTDS_allmonths': # In case of Var_x==Light, take ii_maxmin of the other variable same as that of the Light variable
                vec_new=np.concatenate((Variable_y_PAPA_reshaped, Variable_y_PAPA_reshaped*1.000000000001,Variable_y_PAPA_reshaped), axis=0) # put vector 3 times in a row and nultiply the middle one by a small value so the maximum value alwasy found in the middle matrix
                Var_plot_y_PAPA=( ( vec_new[ii_maxmin-1]+vec_new[ii_maxmin]+vec_new[ii_maxmin+1]) / 3 ) / 1.000000000001
            else:
                vec_new=np.concatenate((Variable_y_PAPA_reshaped, Variable_y_PAPA_reshaped*1.000000000001,Variable_y_PAPA_reshaped), axis=0) # put vector 3 times in a row and nultiply the middle one by a small value so the maximum value alwasy found in the middle matrix
                ii_maxmin=np.argmax(vec_new); Var_plot_y_PAPA=( ( vec_new[ii_maxmin-1]+vec_new[ii_maxmin]+vec_new[ii_maxmin+1]) / 3 ) / 1.000000000001               

            vec_new=np.concatenate((Variable_x_NABE_reshaped, Variable_x_NABE_reshaped*1.000000000001,Variable_x_NABE_reshaped), axis=0) # put vector 3 times in a row and nultiply the middle one by a small value so the maximum value alwasy found in the middle matrix
            ii_maxmin=np.argmax(vec_new); Var_plot_x_NABE=( ( vec_new[ii_maxmin-1]+vec_new[ii_maxmin]+vec_new[ii_maxmin+1]) / 3 ) / 1.000000000001
            if var_names_x[V_i]=='RSNTDS_allmonths': # In case of Var_x==Light, take ii_maxmin of the other variable same as that of the Light variable
                vec_new=np.concatenate((Variable_y_NABE_reshaped, Variable_y_NABE_reshaped*1.000000000001,Variable_y_PAPA_reshaped), axis=0) # put vector 3 times in a row and nultiply the middle one by a small value so the maximum value alwasy found in the middle matrix
                Var_plot_y_NABE=( ( vec_new[ii_maxmin-1]+vec_new[ii_maxmin]+vec_new[ii_maxmin+1]) / 3 ) / 1.000000000001
            else:
                vec_new=np.concatenate((Variable_y_NABE_reshaped, Variable_y_NABE_reshaped*1.000000000001,Variable_y_NABE_reshaped), axis=0) # put vector 3 times in a row and nultiply the middle one by a small value so the maximum value alwasy found in the middle matrix
                ii_maxmin=np.argmax(vec_new); Var_plot_y_NABE=( ( vec_new[ii_maxmin-1]+vec_new[ii_maxmin]+vec_new[ii_maxmin+1]) / 3 ) / 1.000000000001               

        if max_min_ind=='min':
            
            vec_new=np.concatenate((Variable_x_PAPA_reshaped, Variable_x_PAPA_reshaped*1.000000000001,Variable_x_PAPA_reshaped), axis=0) # put vector 3 times in a row and nultiply the middle one by a small value so the maximum value alwasy found in the middle matrix
            ii_maxmin=np.argmin(vec_new); Var_plot_x_PAPA=( ( vec_new[ii_maxmin-1]+vec_new[ii_maxmin]+vec_new[ii_maxmin+1]) / 3 ) / 1.000000000001
            if var_names_x[V_i]=='RSNTDS_allmonths': # In case of Var_x==Light, take ii_maxmin of the other variable same as that of the Light variable
                vec_new=np.concatenate((Variable_y_PAPA_reshaped, Variable_y_PAPA_reshaped*1.000000000001,Variable_y_PAPA_reshaped), axis=0) # put vector 3 times in a row and nultiply the middle one by a small value so the maximum value alwasy found in the middle matrix
                Var_plot_y_PAPA=( ( vec_new[ii_maxmin-1]+vec_new[ii_maxmin]+vec_new[ii_maxmin+1]) / 3 ) / 1.000000000001
            else:
                vec_new=np.concatenate((Variable_y_PAPA_reshaped, Variable_y_PAPA_reshaped*1.000000000001,Variable_y_PAPA_reshaped), axis=0) # put vector 3 times in a row and nultiply the middle one by a small value so the maximum value alwasy found in the middle matrix
                ii_maxmin=np.argmin(vec_new); Var_plot_y_PAPA=( ( vec_new[ii_maxmin-1]+vec_new[ii_maxmin]+vec_new[ii_maxmin+1]) / 3 ) / 1.000000000001               

            vec_new=np.concatenate((Variable_x_NABE_reshaped, Variable_x_NABE_reshaped*1.000000000001,Variable_x_NABE_reshaped), axis=0) # put vector 3 times in a row and nultiply the middle one by a small value so the maximum value alwasy found in the middle matrix
            ii_maxmin=np.argmin(vec_new); Var_plot_x_NABE=( ( vec_new[ii_maxmin-1]+vec_new[ii_maxmin]+vec_new[ii_maxmin+1]) / 3 ) / 1.000000000001
            if var_names_x[V_i]=='RSNTDS_allmonths': # In case of Var_x==Light, take ii_maxmin of the other variable same as that of the Light variable
                vec_new=np.concatenate((Variable_y_NABE_reshaped, Variable_y_NABE_reshaped*1.000000000001,Variable_y_PAPA_reshaped), axis=0) # put vector 3 times in a row and nultiply the middle one by a small value so the maximum value alwasy found in the middle matrix
                Var_plot_y_NABE=( ( vec_new[ii_maxmin-1]+vec_new[ii_maxmin]+vec_new[ii_maxmin+1]) / 3 ) / 1.000000000001
            else:
                vec_new=np.concatenate((Variable_y_NABE_reshaped, Variable_y_NABE_reshaped*1.000000000001,Variable_y_NABE_reshaped), axis=0) # put vector 3 times in a row and nultiply the middle one by a small value so the maximum value alwasy found in the middle matrix
                ii_maxmin=np.argmin(vec_new); Var_plot_y_NABE=( ( vec_new[ii_maxmin-1]+vec_new[ii_maxmin]+vec_new[ii_maxmin+1]) / 3 ) / 1.000000000001               
                
        plt.scatter(Var_plot_x_PAPA, Var_plot_y_PAPA, s=100, marker='o', c=Colors[M_i], label=GCM) # Only the first pannel gets the lables so the legend location be alligned with that
        plt.scatter(Var_plot_x_PAPA, Var_plot_y_PAPA, s=100, marker='o', facecolors='none', edgecolors='k')

        plt.scatter(Var_plot_x_NABE, Var_plot_y_NABE, s=100, marker='d', c=Colors[M_i], label=GCM) # Only the first pannel gets the lables so the legend location be alligned with that
        plt.scatter(Var_plot_x_NABE, Var_plot_y_NABE, s=100, marker='d', facecolors='none', edgecolors='k')
        
        plt.xlabel(var_names_x[V_i][:-10],fontsize=14)
        plt.ylabel(var_names_y[V_i][:-10],fontsize=14)  

    plt.xlim(plot_ranges_x[V_i][0],plot_ranges_x[V_i][1])
    plt.ylim(plot_ranges_y[V_i][0],plot_ranges_y[V_i][1])
    plt.xticks(fontsize=10); plt.yticks(fontsize = 10) # Setting the X-ticks ## rotation='horizontal' # rotation='vertical'   
    if  np.abs( plot_ranges_x[V_i][1] ) < 0.001 :
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))  
    if  np.abs( plot_ranges_y[V_i][1] ) < 0.001 :
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))  
    
plt.subplots_adjust(left=0.09, bottom=0.06, right=0.83, top=0.92, hspace=0.35, wspace=0.35) # the amount of height/width reserved for space between subplots
#plt.legend(shadow=True, loc='right', bbox_to_anchor=(2.1, 2.5), markerscale=1 , prop={'size': 14})
if d_trend=='yes':
    if max_min_ind=='max':
        plt.suptitle('CMIP5 models - Max 3-month value for North Pacific (circles) and North Atlantic (diamonds) - '+str(start_date_plt)+'-'+str(end_date_plt)+' - detrended' , fontsize=18)
    else:
        plt.suptitle('CMIP5 models - Min 3-month value for North Pacific (circles) and North Atlantic (diamonds) - '+str(start_date_plt)+'-'+str(end_date_plt)+' - detrended' , fontsize=18)
else:
    if max_min_ind=='max':
        plt.suptitle('CMIP5 models - Max 3-month value for North Pacific (circles) and North Atlantic (diamonds) - '+str(start_date_plt)+'-'+str(end_date_plt) , fontsize=18)
    else:
        plt.suptitle('CMIP5 models - Min 3-month value for North Pacific (circles) and North Atlantic (diamonds) - '+str(start_date_plt)+'-'+str(end_date_plt) , fontsize=18)

mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full

if d_trend=='yes':
    if max_min_ind=='max':
        fig.savefig(dir_figs+'AllModels_Max_3m_Scatter_detrended_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
    else:
        fig.savefig(dir_figs+'AllModels_Min_3m_Scatter_detrended_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
else:
    if max_min_ind=='max':
        fig.savefig(dir_figs+'AllModels_Max_3m_Scatter_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
    else:
        fig.savefig(dir_figs+'AllModels_Min_3m_Scatter_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')





##################################
###   3D Plotting - Subplots  ####
##################################
GCM_Names = ['CESM1-BGC' ,'GFDL-ESM2G','GFDL-ESM2M','IPSL-CM5A-MR','IPSL-CM5A-LR','MRI-ESM1','GISS-E2-H-CC','GISS-E2-R-CC','HadGEM2-ES','CanESM2','NorESM1-ME','MPI-ESM-MR','MPI-ESM-LR']

Time_months = ('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec');     
colors_months = ['navy', 'blue', 'royalblue', 'lime', 'limegreen', 'green', 'yellow', 'gold', 'orange', 'red', 'magenta','darkviolet']
#### Plots 1 #####
var_names_z=[ 'Phyc_allmonths',                                         'PPint_allmonths',                                       'Zooc_allmonths',                                    'EPC100_allmonths'             ]
var_names_x=[ 'Fe_allmonths',                                           'Fe_allmonths',                                          'Fe_allmonths',                                      'Fe_allmonths'                ]
var_names_y=[ 'NO3_allmonths',                                          'NO3_allmonths',                                         'NO3_allmonths',                                     'NO3_allmonths'              ]
plot_titles=[ 'log [biomass vs. Iron (Fe) vs. Nitrate (NO3)] [mol/m3]', 'log [intPP vs. Iron (Fe) vs. Nitrate (NO3)] [mol/m3]', 'log [Zoo vs. Iron (Fe) vs. Nitrate (NO3)] [mol/m3]', 'log [EPC 100m vs. Iron (Fe) vs. Nitrate (NO3)] [mol/m3]'   ]
save_names= [ 'Phyc_vs_Fe_vs_NO3',                                      'PPint_vs_Fe_vs_NO3',                                    'Zooc_vs_Fe_vs_NO3',                                 'EPC100_vs_Fe_vs_NO3'                 ]  

file_names_z=[ 'bio',                                                    'bio',                                                   'bio',                                              'bio'                 ]
file_names_x=[ 'bio',                                                    'bio',                                                   'bio',                                              'bio'                  ]
file_names_y=[ 'bio',                                                    'bio',                                                   'bio',                                              'bio'                  ]
#### Plots 1 #####
var_names_z=[ 'Phyc_allmonths',                                     'PPint_allmonths',                                  'Zooc_allmonths',                                 'EPC100_allmonths'             ]
var_names_x=[ 'RSNTDS_allmonths',                                   'RSNTDS_allmonths',                                 'RSNTDS_allmonths',                               'RSNTDS_allmonths'             ]
var_names_y=[ 'MLD_allmonths',                                      'MLD_allmonths',                                    'MLD_allmonths',                                  'MLD_allmonths'                ]
plot_titles=[ 'log (biomass)[mol/m3] vs. light [w/m2] vs. MLD [m]', 'log (intPP)[mol/m3] vs. light [w/m2] vs. MLD [m]', 'log (Zoo)[mol/m3] vs. light [w/m2] vs. MLD [m]', 'log (EPC 100m)[mol/m3] vs. light [w/m2] vs. MLD [m]'   ]
save_names= [ 'Phyc_vs_RSNTDS_vs_MLD',                              'PPint_vs_RSNTDS_vs_MLD',                           'Zooc_vs_RSNTDS_vs_MLD',                          'EPC100_vs_RSNTDS_vs_MLD'                 ]  

file_names_z=[ 'bio',                                               'bio',                                              'bio',                                            'bio'                      ]
file_names_x=[ 'physics',                                           'physics',                                          'physics',                                        'physics'                  ]
file_names_y=[ 'physics',                                           'physics',                                          'physics',                                        'physics'                  ]

#plot_ranges_y=[[0,0.006],  [-2, -7]]
#plot_ranges_x=[[0,0.006],  [-1, -7]]

d_trend='yes'
n_r=4 ; n_c=4 ;
n_t=len(GCM_Names)
box_n='PAPA'   ;   box_n='NABE'

for V_i in range(n_t): # V_i=0   V_i=1    V_i=2    V_i=3    V_i=4    V_i=5    V_i=6    V_i=7    V_i=8    V_i=9    V_i=10    V_i=11    V_i=12    V_i=13

    if d_trend=='yes':
        if box_n=='PAPA':
            P_title=plot_titles[V_i]+' - North Pacific [45N-50N, 140W-150W] - '+str(start_date_plt)+'-'+str(end_date_plt)+' - detrended'
        elif box_n=='NABE':
            P_title=plot_titles[V_i]+' - North Atlantic [45N-50N, 25W-35W] - '+str(start_date_plt)+'-'+str(end_date_plt)+' - detrended'
    else:
        if box_n=='PAPA':
            P_title=plot_titles[V_i]+' - North Pacific [45N-50N, 140W-150W] - '+str(start_date_plt)+'-'+str(end_date_plt)
        elif box_n=='NABE':
            P_title=plot_titles[V_i]+' - North Atlantic [45N-50N, 25W-35W] - '+str(start_date_plt)+'-'+str(end_date_plt)
    
    fig=plt.figure()
    for M_i in range(len(GCM_Names)): # M_i=1     # GCM=GCM_Names[0]

        GCM=GCM_Names[M_i]
        print (GCM) 
        
        ax = fig.add_subplot(n_r,n_c,M_i+1,projection='3d') 

        filename_in_z = (dir_data_in + 'AllResults_'+GCM+'_'+file_names_z[V_i]+'_1991_2010.out') # Directory to save processed data
        filename_in_x = (dir_data_in + 'AllResults_'+GCM+'_'+file_names_x[V_i]+'_1991_2010.out') # Directory to save processed data
        filename_in_y = (dir_data_in + 'AllResults_'+GCM+'_'+file_names_y[V_i]+'_1991_2010.out') # Directory to save processed data

        my_shelf_z = shelve.open(filename_in_z)
        my_shelf_x = shelve.open(filename_in_x)
        my_shelf_y = shelve.open(filename_in_y)
 
        if var_names_x[V_i] in my_shelf_x and var_names_y[V_i] in my_shelf_y and var_names_z[V_i] in my_shelf_z: # If that variable exists in the saved data  
 
            globals()['Variable_z']=my_shelf_z[var_names_z[V_i]]            
            globals()['Variable_x']=my_shelf_x[var_names_x[V_i]] 
            globals()['Variable_y']=my_shelf_y[var_names_y[V_i]]
                
            Variable_x[ np.abs(Variable_x)>1e19 ]=nan; Variable_y[ np.abs(Variable_y)>1e19 ]=nan; Variable_z[ np.abs(Variable_z)>1e19 ]=nan  # Replacing missing values with NaN 

            if box_n=='PAPA':
                Variable_z_BOX = Variable_z[:,135:140,210:220]
                Variable_x_BOX = Variable_x[:,135:140,210:220]
                Variable_y_BOX = Variable_y[:,135:140,210:220]

            elif box_n=='NABE':
                Variable_z_BOX = Variable_z[:,135:140,325:335]
                Variable_x_BOX = Variable_x[:,135:140,325:335]
                Variable_y_BOX = Variable_y[:,135:140,325:335]
                    
            if d_trend=='yes': # If it is requested to detrend the data before correlations
                Variable_z_BOX= func_detrend_3d(Variable_z_BOX)
                Variable_x_BOX= func_detrend_3d(Variable_x_BOX)
                Variable_y_BOX= func_detrend_3d(Variable_y_BOX)

            if var_names_z[V_i] !='RSNTDS_allmonths' and var_names_z[V_i] !='MLD_allmonths' : # For all variables except these ones, scatterplots are done on log10 of the variable
                Variable_z_BOX=np.log10(Variable_z_BOX)
            if var_names_x[V_i] !='RSNTDS_allmonths' and var_names_x[V_i] !='MLD_allmonths' : # For all variables except these ones, scatterplots are done on log10 of the variable
                Variable_x_BOX=np.log10(Variable_x_BOX)            
            if var_names_y[V_i] !='RSNTDS_allmonths' and var_names_y[V_i] !='MLD_allmonths' : # For all variables except these ones, scatterplots are done on log10 of the variable
                Variable_y_BOX=np.log10(Variable_y_BOX)               
            
            Variable_z_BOX_reshaped = Variable_z_BOX.reshape(( yrs_n ,12,Variable_z_BOX.shape[1],Variable_z_BOX.shape[2])) # Monthly means of SpCO2 over the time period  
            Variable_x_BOX_reshaped = Variable_x_BOX.reshape(( yrs_n ,12,Variable_x_BOX.shape[1],Variable_x_BOX.shape[2])) # Monthly means of SpCO2 over the time period  
            Variable_y_BOX_reshaped = Variable_y_BOX.reshape(( yrs_n ,12,Variable_y_BOX.shape[1],Variable_y_BOX.shape[2])) # Monthly means of SpCO2 over the time period 
    
#            for ii in range(Variable_z_BOX_reshaped.shape[1]):
#                ax.scatter(Variable_x_BOX_reshaped[:,ii,:,:], Variable_y_BOX_reshaped[:,ii,:,:], Variable_z_BOX_reshaped[:,ii,:,:], c= colors_months[ii] , s=0.25, marker='x', label=Time_months[ii] )
            
            Variable_y_BOX_reshaped=np.nanmean(Variable_y_BOX_reshaped,axis=0)
            Variable_x_BOX_reshaped=np.nanmean(Variable_x_BOX_reshaped,axis=0)
            Variable_z_BOX_reshaped=np.nanmean(Variable_z_BOX_reshaped,axis=0)
    
            for ii in range(Variable_y_BOX_reshaped.shape[0]):
                ax.scatter(Variable_x_BOX_reshaped[ii,:,:], Variable_y_BOX_reshaped[ii,:,:], Variable_z_BOX_reshaped[ii,:,:], c= colors_months[ii] , s=3 , marker='x', label=Time_months[ii] )
    
            ax.set_xlabel(var_names_x[V_i][:-10],fontsize=10, labelpad=14)
            ax.set_ylabel(var_names_y[V_i][:-10],fontsize=10, labelpad=20)
            ax.set_zlabel(var_names_z[V_i][:-10],fontsize=10, labelpad=3, rotation=90)
            plt.xticks(rotation=45, fontsize=10); plt.yticks(rotation=-45, fontsize=10)

        plt.title(GCM, fontsize=14)
        plt.subplots_adjust(left=0.06, bottom=0.04, right=0.83, top=0.92, hspace=0.25, wspace=0.1) # the amount of height/width reserved for space between subplots
    plt.legend(shadow=True, loc='right', bbox_to_anchor=(5, 3), markerscale=40 , prop={'size': 16})
    plt.suptitle(P_title, fontsize=20)    
    plt.show()
    mng = plt.get_current_fig_manager()
    mng.window.showMaximized() # Maximizes the plot window to save figures in full
    
    if d_trend=='yes':
        if box_n=='PAPA':
            fig.savefig(dir_figs+'AllModels_ScatterPlots_3D_PAPA_detrended_'+save_names[V_i]+'_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
        elif box_n=='NABE':
            fig.savefig(dir_figs+'AllModels_ScatterPlots_3D_NABE_detrended_'+save_names[V_i]+'_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
    else:
        if box_n=='PAPA':
            fig.savefig(dir_figs+'AllModels_ScatterPlots_3D_PAPA_'+save_names[V_i]+'_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
        elif box_n=='NABE':
            fig.savefig(dir_figs+'AllModels_ScatterPlots_3D_NABE_'+save_names[V_i]+'_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')



for V_i in range(n_t): # V_i=0   V_i=1    V_i=2    V_i=3    V_i=4    V_i=5    V_i=6    V_i=7    V_i=8    V_i=9    V_i=10    V_i=11    V_i=12    V_i=13

    if d_trend=='yes':
        P_title=plot_titles[V_i]+' - North Pacific (circles) and North Atlantic (diamonds) - '+str(start_date_plt)+'-'+str(end_date_plt)+' - detrended'
    else:
        P_title=plot_titles[V_i]+' - North Pacific (circles) and North Atlantic (diamonds) - '+str(start_date_plt)+'-'+str(end_date_plt)
    
    fig=plt.figure()
    for M_i in range(len(GCM_Names)): # M_i=1     # GCM=GCM_Names[0]

        GCM=GCM_Names[M_i]
        print (GCM) 
        
        ax = fig.add_subplot(n_r,n_c,M_i+1,projection='3d') 

        filename_in_z = (dir_data_in + 'AllResults_'+GCM+'_'+file_names_z[V_i]+'_1991_2010.out') # Directory to save processed data
        filename_in_x = (dir_data_in + 'AllResults_'+GCM+'_'+file_names_x[V_i]+'_1991_2010.out') # Directory to save processed data
        filename_in_y = (dir_data_in + 'AllResults_'+GCM+'_'+file_names_y[V_i]+'_1991_2010.out') # Directory to save processed data

        my_shelf_z = shelve.open(filename_in_z)
        my_shelf_x = shelve.open(filename_in_x)
        my_shelf_y = shelve.open(filename_in_y)
 
        if var_names_x[V_i] in my_shelf_x and var_names_y[V_i] in my_shelf_y and var_names_z[V_i] in my_shelf_z: # If that variable exists in the saved data  
 
            globals()['Variable_z']=my_shelf_z[var_names_z[V_i]]            
            globals()['Variable_x']=my_shelf_x[var_names_x[V_i]] 
            globals()['Variable_y']=my_shelf_y[var_names_y[V_i]]
                
            Variable_x[ np.abs(Variable_x)>1e19 ]=nan; Variable_y[ np.abs(Variable_y)>1e19 ]=nan; Variable_z[ np.abs(Variable_z)>1e19 ]=nan  # Replacing missing values with NaN 

            Variable_z_PAPA = Variable_z[:,135:140,210:220]
            Variable_x_PAPA = Variable_x[:,135:140,210:220]
            Variable_y_PAPA = Variable_y[:,135:140,210:220]

            Variable_z_NABE = Variable_z[:,135:140,325:335]
            Variable_x_NABE = Variable_x[:,135:140,325:335]
            Variable_y_NABE = Variable_y[:,135:140,325:335]
                    
            if d_trend=='yes': # If it is requested to detrend the data before correlations
                Variable_z_PAPA= func_detrend_3d(Variable_z_PAPA)
                Variable_x_PAPA= func_detrend_3d(Variable_x_PAPA)
                Variable_y_PAPA= func_detrend_3d(Variable_y_PAPA)                
                Variable_z_NABE= func_detrend_3d(Variable_z_NABE)
                Variable_x_NABE= func_detrend_3d(Variable_x_NABE)
                Variable_y_NABE= func_detrend_3d(Variable_y_NABE)

            if var_names_z[V_i] !='RSNTDS_allmonths' and var_names_z[V_i] !='MLD_allmonths' : # For all variables except these ones, scatterplots are done on log10 of the variable
                Variable_z_PAPA=np.log10(Variable_z_PAPA)
                Variable_z_NABE=np.log10(Variable_z_NABE)
            if var_names_x[V_i] !='RSNTDS_allmonths' and var_names_x[V_i] !='MLD_allmonths' : # For all variables except these ones, scatterplots are done on log10 of the variable
                Variable_x_PAPA=np.log10(Variable_x_PAPA) 
                Variable_x_NABE=np.log10(Variable_x_NABE) 
            if var_names_y[V_i] !='RSNTDS_allmonths' and var_names_y[V_i] !='MLD_allmonths' : # For all variables except these ones, scatterplots are done on log10 of the variable
                Variable_y_PAPA=np.log10(Variable_y_PAPA)  
                Variable_y_NABE=np.log10(Variable_y_NABE)  
            
            Variable_z_PAPA_reshaped = Variable_z_PAPA.reshape(( yrs_n ,12,Variable_z_PAPA.shape[1],Variable_z_PAPA.shape[2])) # Monthly means of SpCO2 over the time period  
            Variable_x_PAPA_reshaped = Variable_x_PAPA.reshape(( yrs_n ,12,Variable_x_PAPA.shape[1],Variable_x_PAPA.shape[2])) # Monthly means of SpCO2 over the time period  
            Variable_y_PAPA_reshaped = Variable_y_PAPA.reshape(( yrs_n ,12,Variable_y_PAPA.shape[1],Variable_y_PAPA.shape[2])) # Monthly means of SpCO2 over the time period 

            Variable_z_NABE_reshaped = Variable_z_NABE.reshape(( yrs_n ,12,Variable_z_NABE.shape[1],Variable_z_NABE.shape[2])) # Monthly means of SpCO2 over the time period  
            Variable_x_NABE_reshaped = Variable_x_NABE.reshape(( yrs_n ,12,Variable_x_NABE.shape[1],Variable_x_NABE.shape[2])) # Monthly means of SpCO2 over the time period  
            Variable_y_NABE_reshaped = Variable_y_NABE.reshape(( yrs_n ,12,Variable_y_NABE.shape[1],Variable_y_NABE.shape[2])) # Monthly means of SpCO2 over the time period 
                
            Variable_y_PAPA_reshaped=np.nanmean(np.nanmean(np.nanmean(Variable_y_PAPA_reshaped,axis=0),axis=1) ,axis=1)
            Variable_x_PAPA_reshaped=np.nanmean(np.nanmean(np.nanmean(Variable_x_PAPA_reshaped,axis=0),axis=1) ,axis=1)
            Variable_z_PAPA_reshaped=np.nanmean(np.nanmean(np.nanmean(Variable_z_PAPA_reshaped,axis=0),axis=1) ,axis=1)
            
            Variable_y_NABE_reshaped=np.nanmean(np.nanmean(np.nanmean(Variable_y_NABE_reshaped,axis=0),axis=1) ,axis=1)
            Variable_x_NABE_reshaped=np.nanmean(np.nanmean(np.nanmean(Variable_x_NABE_reshaped,axis=0),axis=1) ,axis=1)
            Variable_z_NABE_reshaped=np.nanmean(np.nanmean(np.nanmean(Variable_z_NABE_reshaped,axis=0),axis=1) ,axis=1)
            
            ax.plot(Variable_x_PAPA_reshaped, Variable_y_PAPA_reshaped, Variable_z_PAPA_reshaped, color='k', lw=1)
            for ii in range(Variable_y_PAPA_reshaped.shape[0]):
                ax.scatter(Variable_x_PAPA_reshaped[ii], Variable_y_PAPA_reshaped[ii], Variable_z_PAPA_reshaped[ii], c= colors_months[ii] , s=50 , marker='o',zorder=2, label=Time_months[ii] )
            #ax.text(Variable_x_PAPA_reshaped[4], Variable_y_PAPA_reshaped[4], Variable_z_PAPA_reshaped[4], "N.P.", color='k', fontsize=14)

            ax.plot(Variable_x_NABE_reshaped, Variable_y_NABE_reshaped, Variable_z_NABE_reshaped, color='k', lw=1, linestyle='--')
            for ii in range(Variable_y_NABE_reshaped.shape[0]):
                ax.scatter(Variable_x_NABE_reshaped[ii], Variable_y_NABE_reshaped[ii], Variable_z_NABE_reshaped[ii], c= colors_months[ii] , s=50 , marker='d',zorder=2, label=Time_months[ii] )
            #ax.text(Variable_x_NABE_reshaped[4], Variable_y_NABE_reshaped[4], Variable_z_NABE_reshaped[4], "N.A.", color='k', fontsize=14)
             
            ax.set_xlabel(var_names_x[V_i][:-10],fontsize=10, labelpad=14)
            ax.set_ylabel(var_names_y[V_i][:-10],fontsize=10, labelpad=20)
            ax.set_zlabel(var_names_z[V_i][:-10],fontsize=10, labelpad=3, rotation=90)
            plt.xticks(rotation=45, fontsize=10); plt.yticks(rotation=-45, fontsize=10)

        plt.title(GCM, fontsize=14)
        plt.subplots_adjust(left=0.06, bottom=0.04, right=0.83, top=0.92, hspace=0.25, wspace=0.1) # the amount of height/width reserved for space between subplots
    #plt.legend(shadow=True, loc='right', bbox_to_anchor=(5, 3), markerscale=2 , prop={'size': 16})
    plt.suptitle(P_title, fontsize=20)    
    plt.show()
    mng = plt.get_current_fig_manager()
    mng.window.showMaximized() # Maximizes the plot window to save figures in full

    if d_trend=='yes':
        fig.savefig(dir_figs+'AllModels_ScatterPlots_3D_ave_PAPA_NABE_detrended_'+save_names[V_i]+'_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
    else:
        fig.savefig(dir_figs+'AllModels_ScatterPlots_3D_ave_PAPA_NABE_'+save_names[V_i]+'_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')




#######################################################
###   Multi-GCM MultiVariable Plotting - Subplots  ####
####################################################### 
GCM_Names = ['CESM1-BGC' ,'GFDL-ESM2G','GFDL-ESM2M','IPSL-CM5A-MR','IPSL-CM5A-LR','MRI-ESM1','GISS-E2-H-CC','GISS-E2-R-CC','HadGEM2-ES','CanESM2','NorESM1-ME','MPI-ESM-MR','MPI-ESM-LR']
Colors    =  [[1, 1, 0],   [1, 0.8, 0],  [1, 0.6, 0], [1, 0.4, 0],   [1, 0, 0], [1, 0.8, 1],   [1, 0.5, 1],   [1, 0, 1],   [0, 0, 0.5], [0, 0, 1], [0, 0.5, 1],[0, 0.5, 0],  [0, 1, 0]]

Time_months = ('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec');
P_Var_x = np.linspace(1,12,12);

plot_ranges=  [ [-1,20],              [-0.5,13],                [-0.05,1.4],                              [-0.1e-1,2.7e-1],                              [-20,40],                                        [270,460],          [-0.1,3.8],                [-1,18],                    [0,250]]
var_names=    [ 'Phyc_allmonths',     'Zooc_allmonths',         'PPint_allmonths',                       'EPC100_allmonths',                             'FgCO2_allmonths',                                'SpCO2_allmonths', 'Fe_allmonths',             'NO3_allmonths',           'RSNTDS_allmonths']
plot_titles=  [ 'Biomass (m mol/m3)', 'Zooplankton (m mol/m3)', 'Primary Production (micromol.m-2.s-1)', 'Export Production at 100m (micromol.m-2.s-1)', 'Air-sea CO2 flux (donward +) (gram C /m2/year)', 'SpCO2 (ppm)',     'Iron - Fe (micro mol/m3)', 'Nitrate - NO3 (m mol/m3)', 'Light - Surface Downward SW Radiation (W/m2)']  
file_names= [   'bio',                 'bio',                    'bio',                                  'bio',                                          'bio',                                            'bio',             'bio',                      'bio',                      'physics'                ]
var_multiplier= [ 1000,                1000,                      1E6,                                    1E6,                                            1,                                                1,                 1E6,                        1000,                       1     ]
colors_vars = ['navy', 'blue', 'royalblue', 'lime', 'limegreen', 'green', 'gold', 'orange', 'red', 'magenta','darkviolet']

colors_vars_1 = [ 'blue', 'lime', 'limegreen', 'green', 'gold', 'orange',  'magenta','darkviolet']
colors_vars_2 = ['navy', 'blue', 'royalblue', 'lime', 'limegreen', 'green', 'gold', 'orange', 'red', 'magenta','darkviolet']

plot_ranges_1=  [ [-1,20],                    [-1,20],                   [0,2.5],                 [-0.01,0.5],                              [-0.05,1.4],                              [-0.1e-1,2.7e-1],                              [-0.1,1.5],                 [-1,18],                    [0,240]]
var_names_1=    [ 'Phyc_allmonths',           'Zooc_allmonths',         'Zoo_by_Phyto_allmonths', 'PP_by_Phyto_allmonths',                  'PPint_allmonths',                       'EPC100_allmonths',                             'Fe_allmonths',             'NO3_allmonths',            'RSNTDS_allmonths']
plot_titles_1=  [ 'Total Biomass (m mol/m3)', 'Zooplankton (m mol/m3)', 'Zoo/phyto',              'Specific Growth rate (PP/phyto * 1E-3)', 'Primary Production (micromol.m-2.s-1)', 'Export Production at 100m (micromol.m-2.s-1)', 'Iron - Fe (micro mol/m3)', 'Nitrate - NO3 (m mol/m3)', 'Light - Surface Downward SW Radiation (W/m2)']  
file_names_1= [   'bio',                      'bio',                    'bio',                    'bio',                                    'bio',                                   'bio',                                           'bio',                     'bio',                      'physics'                ]
var_multiplier_1=[ 1000,                      1000,                      1,                        1000,                                     1E6,                                     1E6,                                            1E6,                        1000,                        1         ]
colors_vars_1 = [ 'red',                     'navy',                    'deepskyblue',            'grey',                                  'lime',                                  'green',                                        'gold',                     'orange',                     'magenta'                            ]


plot_ranges_2=  [ [-1,10],                  [-1,110],                      [-1,5],                   [4,20],          [-110,215],                             [0,1000],                 [-20,40],                                        [270,460],          ]
var_names_2=    [ 'Large Phyto_allmonths',  '% Large Phyto_allmonths',    'Small Phyto_allmonths',  'SST_allmonths', 'HFDS_allmonths',                       'MLD_allmonths',         'FgCO2_allmonths',                                'SpCO2_allmonths', ]
plot_titles_2=  [ 'Large Phyto (m mol/m3)', 'Large Phyto percentage (%)', 'Small Phyto (m mol/m3)', 'SST (C)',       'Downward Heat Flux at Surface (W/m2)', 'Mixed Layer Depth (m)', 'Air-sea CO2 flux (donward +) (gram C /m2/year)', 'SpCO2 (ppm)',    ]  
file_names_2=   [ 'bio',                    'bio',                        'bio',                    'physics' ,      'physics' ,                             'physics',               'bio',                                            'bio',              ]
var_multiplier_2=[1000,                     1,                            1000,                      1,               1,                                      1,                       1,                                                 1,               ]
colors_vars_2 = ['navy',                   'blue',                        'deepskyblue',            'orange',        'crimson',                              'darkviolet',            'lime',                                            'green']


d_trend='yes'

for M_i in range(len(GCM_Names)): # M_i=8     # GCM=GCM_Names[0]
    
    #M_i=0
    
    GCM=GCM_Names[M_i]
    print (GCM)  

    fig=plt.figure() #; fig, ax = plt.subplots()     
    for AX_i in range(4):    
        ax = fig.add_subplot(2,2,AX_i+1)  
        axes = [ax, ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx()] 
        
        if AX_i+1==1 or AX_i+1==2:  # Axes 1 and 2 (first row) are for the first set of variables
            var_names=var_names_1 ; plot_titles=plot_titles_1 ; file_names=file_names_1 ; var_multiplier=var_multiplier_1 ; colors_vars=colors_vars_1 ; plot_ranges=plot_ranges_1
        elif AX_i+1==3 or AX_i+1==4:   # Axes 3 and 4 (second row) are for the second set of variables
            var_names=var_names_2 ; plot_titles=plot_titles_2 ; file_names=file_names_2 ; var_multiplier=var_multiplier_2 ; colors_vars=colors_vars_2 ; plot_ranges=plot_ranges_2
        
        V_i=-1
        for ax in axes:            
            V_i += 1
            if V_i < len(var_names):
                if AX_i+1==1 or AX_i+1==3:
                    if V_i>0:
                        axes[V_i].spines['right'].set_position(('axes', -0.1-V_i*0.13)); axes[V_i].set_frame_on(True); axes[V_i].patch.set_visible(False)
                else:
                    axes[V_i].spines['right'].set_position(('axes', -10-V_i*0.14)); axes[V_i].set_frame_on(True); axes[V_i].patch.set_visible(False)
                
                if var_names[V_i]=='Zoo_by_Phyto_allmonths':
                    filename_in = (dir_data_in + 'AllResults_'+GCM+'_'+'bio'+'_1991_2010.out') # Directory to save processed data
                    my_shelf = shelve.open(filename_in)          
                    
                    if 'Zooc_allmonths' in my_shelf and 'Phyc_allmonths' in my_shelf: # If that variable exists in the saved data                       
                        globals()['Variable_1']=my_shelf['Zooc_allmonths'] 
                        globals()['Variable_2']=my_shelf['Phyc_allmonths']
                        Variable_P=Variable_1/Variable_2
                    else:
                        Variable_P=empty((yrs_n*12,lat_n_regrid,lon_n_regrid)) * nan              

                elif var_names[V_i]=='PP_by_Phyto_allmonths' or  var_names[V_i]=='LnPP_by_LnPhyto_allmonths' :
                    
                    filename_in = (dir_data_in + 'AllResults_'+GCM+'_'+'bio'+'_1991_2010.out') # Directory to save processed data
                    my_shelf = shelve.open(filename_in)         
                    if 'PPint_allmonths' in my_shelf and 'Phyc_allmonths' in my_shelf: # If that variable exists in the saved data                       
                        globals()['Variable_1']=my_shelf['PPint_allmonths'] 
                        globals()['Variable_2']=my_shelf['Phyc_allmonths']
                        
                        if var_names[V_i]=='LnPP_by_LnPhyto_allmonths':
                            Variable_P=np.log(Variable_1)/np.log(Variable_2)
                        else:
                            Variable_P=Variable_1/Variable_2
                    else:
                        Variable_P=empty((yrs_n*12,lat_n_regrid,lon_n_regrid)) * nan

                elif var_names[V_i]=='Large Phyto_allmonths' or  var_names[V_i]=='% Large Phyto_allmonths' :
                    filename_in = (dir_data_in + 'AllResults_'+GCM+'_'+'bio'+'_1991_2010.out') # Directory to save processed data
                    my_shelf = shelve.open(filename_in)

                    if 'Phyc_allmonths' in my_shelf and 'Phypico_allmonths' in my_shelf: # If that variable exists in the saved data                       
                        globals()['Variable_1']=my_shelf['Phyc_allmonths'] 
                        globals()['Variable_2']=my_shelf['Phypico_allmonths']                       
                        
                        if var_names[V_i]=='Large Phyto_allmonths':
                            Variable_P=Variable_1 - Variable_2
                        elif var_names[V_i]=='% Large Phyto_allmonths':
                            Variable_P= ((Variable_1 - Variable_2) / Variable_1) * 100

                    elif 'Phyc_allmonths' in my_shelf and 'Phydiat_allmonths' in my_shelf: # If that variable exists in the saved data                       
                        globals()['Variable_1']=my_shelf['Phyc_allmonths'] 
                        globals()['Variable_2']=my_shelf['Phydiat_allmonths']                       
                        
                        if var_names[V_i]=='Large Phyto_allmonths':
                            Variable_P=Variable_2
                        elif var_names[V_i]=='% Large Phyto_allmonths':
                            Variable_P= (Variable_2 / Variable_1) * 100

                    else:
                        Variable_P=empty((yrs_n*12,lat_n_regrid,lon_n_regrid)) * nan

                elif var_names[V_i]=='Small Phyto_allmonths':
                    filename_in = (dir_data_in + 'AllResults_'+GCM+'_'+'bio'+'_1991_2010.out') # Directory to save processed data
                    my_shelf = shelve.open(filename_in)

                    if 'Phypico_allmonths' in my_shelf: # If that variable exists in the saved data                       
                        globals()['Variable_2']=my_shelf['Phypico_allmonths']  
                        Variable_P=Variable_2

                    elif 'Phyc_allmonths' in my_shelf and 'Phydiat_allmonths' in my_shelf: # If that variable exists in the saved data                       
                        globals()['Variable_1']=my_shelf['Phyc_allmonths'] 
                        globals()['Variable_2']=my_shelf['Phydiat_allmonths']
                        Variable_P = Variable_1 - Variable_2

                    else:
                        Variable_P=empty((yrs_n*12,lat_n_regrid,lon_n_regrid)) * nan

                else:                   

                    filename_in = (dir_data_in + 'AllResults_'+GCM+'_'+file_names[V_i]+'_1991_2010.out') # Directory to save processed data
                    my_shelf = shelve.open(filename_in) 
                    if var_names[V_i] in my_shelf: # If that variable exists in the saved data    
                        globals()['Variable_1']=my_shelf[var_names[V_i]]
                        Variable_P=Variable_1
                    else:
                        Variable_P=empty((yrs_n*12,lat_n_regrid,lon_n_regrid)) * nan
                
                Variable_P=Variable_P * var_multiplier[V_i]
                Variable_P[ np.abs(Variable_P)>1e19 ]=nan
                Var_P_allmonths=copy.deepcopy(Variable_P); 
                if AX_i+1==1 or AX_i+1==3: # Plot PAPA results in these axes
                    Var_P_allmonths_BOX=Var_P_allmonths[:,135:140,210:220]
                elif AX_i+1==2 or AX_i+1==4: # Plot NABEresults in these axes
                    Var_P_allmonths_BOX=Var_P_allmonths[:,135:140,325:335]
        
                if d_trend=='yes': # If it is requested to detrend the data before calculations
                    Var_P_allmonths_BOX= func_detrend_3d(Var_P_allmonths_BOX)
          
                Var_P_monthlymean_BOX = Var_P_allmonths_BOX.reshape(( np.int(Var_P_allmonths_BOX.shape[0]/12) ,12,Var_P_allmonths_BOX.shape[1],Var_P_allmonths_BOX.shape[2])) # Monthly means of SST over the time period
                Var_P_monthlymean_BOX = np.nanmean(Var_P_monthlymean_BOX,axis=0)
                Var_P_monthlymean_BOX=np.squeeze(Var_P_monthlymean_BOX)           
                Var_P_monthlymean_BOX_timeseries=np.nanmean( np.nanmean(Var_P_monthlymean_BOX,axis=2) , axis=1); Var_P_monthlymean_BOX_timeseries=np.squeeze(Var_P_monthlymean_BOX_timeseries)
           
                ax.plot(P_Var_x, Var_P_monthlymean_BOX_timeseries, c=colors_vars[V_i], label=plot_titles[V_i], marker='o', markersize=7, markerfacecolor='W', linewidth=2.0)    
                ax.set_ylim(plot_ranges[V_i])     
                ax.set_ylabel(plot_titles[V_i], color=colors_vars[V_i])#, labelpad = -5-V_i*5)
                ax.tick_params(axis='y', colors=colors_vars[V_i])
            
        plt.xticks(P_Var_x, Time_months, rotation=45, fontsize=12)
        #axes[0].set_xticks(Time_months)
        if AX_i+1==1 or AX_i+1==3: # Plot PAPA results in these axes
            ax.set_title('North Pacific [45N-50N, 140W-150W]')
        elif AX_i+1==2 or AX_i+1==4: # Plot NABEresults in these axes
            ax.set_title('North Atlantic [45N-50N, 25W-35W]')  
            axes[0].set_ylabel(' ')

    plt.subplots_adjust(left=0.36, bottom=0.05, right=0.96, top=0.9, hspace=0.2, wspace=0.1) # the amount of height/width reserved for space between subplots  
    plt.show()
    mng = plt.get_current_fig_manager()
    mng.window.showMaximized() # Maximizes the plot window to save figures in full        
    if d_trend=='yes':
        plt.suptitle( GCM+' - Variable seasonal cylce - '+str(start_date_plt)+'-'+str(end_date_plt)+ ' - detrended' , fontsize=18)
    else:
        plt.suptitle( GCM+' - Variable seasonal cylce - '+str(start_date_plt)+'-'+str(end_date_plt) , fontsize=18)

    if d_trend=='yes':
        fig.savefig(dir_figs+GCM+'_SesonalCycle_AllVar_detrended_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
    else:
        fig.savefig(dir_figs+GCM+'_SesonalCycle_AllVar_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
































