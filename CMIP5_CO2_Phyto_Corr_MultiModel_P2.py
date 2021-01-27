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
import warnings; warnings.simplefilter('ignore')
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


#####################################################
#####################################################
GCM_Names = ['CESM1-BGC' ,'GFDL-ESM2G','GFDL-ESM2M','IPSL-CM5A-MR','IPSL-CM5A-LR','MRI-ESM1','GISS-E2-H-CC','GISS-E2-R-CC','HadGEM2-ES','CanESM2','NorESM1-ME','MPI-ESM-MR','MPI-ESM-LR']
Colors    =  [[1, 1, 0],   [1, 0.8, 0],  [1, 0.6, 0], [1, 0.4, 0],   [1, 0, 0], [1, 0.8, 1],   [1, 0.5, 1],   [1, 0, 1],   [0, 0, 0.5], [0, 0, 1], [0, 0.5, 1],[0, 0.5, 0],  [0, 1, 0]]

plot_ranges_y=[[270, 370],    [0, 160],       [0, 160],            [1,4000]     ]
plot_ranges_x=[[-0.0001,0.008], [-0.001,0.02],  [-1E-7,0.0000015],   [0,6.5]       ]
var_names_y=[ 'SpCO2_nonT (JJA ave)', 'SpCO2_nonT (max - min)', 'SpCO2_nonT (max - min)',  '(SST(max)-SST(min))/(NO3(max)-NO3(min))'  ]
var_names_x=[ 'Nitrate (JJA ave)',    'Nitrate (max - min)',    'intPP (max - min)',       '(pCO2-T(max)-pCO2-T(min))/(pCO2-nonT(max)-pCO2-nonT(min))'  ]


d_trend='yes'
n_r=2 ; n_c=2 ;
n_t=4

fig=plt.figure()
for V_i in range(n_t): # V_i=0   V_i=1    V_i=2    V_i=3    V_i=4    V_i=5    V_i=6    V_i=7    V_i=8    V_i=9    V_i=10    V_i=11    V_i=12    V_i=13

    ax = fig.add_subplot(n_r,n_c,V_i+1) 

    for M_i in range(len(GCM_Names)): # M_i=0     # GCM=GCM_Names[0]

        GCM=GCM_Names[M_i]
        print (GCM) 

        Variable_y=empty((lat_n_regrid,lon_n_regrid)) * nan
        filename_in_spco2 = (dir_data_in + 'AllResults_'+GCM+'_'+'bio'+'_1991_2010.out') # Directory to save processed data
        filename_in_sst = (dir_data_in + 'AllResults_'+GCM+'_'+'physics'+'_1991_2010.out') # Directory to save processed data

        my_shelf_spco2 = shelve.open(filename_in_spco2)
        my_shelf_sst = shelve.open(filename_in_sst)
        
        if 'SpCO2_allmonths' in my_shelf_spco2  and 'SST_allmonths' in my_shelf_sst: # If that variable exists in the saved data  
            
            globals()['Var_spco2']=my_shelf_spco2['SpCO2_allmonths']
            globals()['Var_sst']=my_shelf_sst['SST_allmonths']                 
                           
            Var_spco2[ np.abs(Var_spco2)>1e19 ]=nan; Var_sst[ np.abs(Var_sst)>1e19 ]=nan  # Replacing missing values with NaN 
            
            Var_spco2_TEMP_eq1 = empty((Var_spco2.shape[0], Var_spco2.shape[1],Var_spco2.shape[2])) * nan #SpCO2_allmonths_TEMP_eq1 = np.zeros([SpCO2_allmonths.shape[0], lat_n_regrid,lon_n_regrid])
            Var_spco2_NonTEMP_eq2 = empty((Var_spco2.shape[0], Var_spco2.shape[1],Var_spco2.shape[2])) * nan #SpCO2_allmonths_NonTEMP_eq2 = np.zeros([SpCO2_allmonths.shape[0], lat_n_regrid,lon_n_regrid])
            for ii in range(Var_spco2.shape[0]):
                Var_spco2_TEMP_eq1[ii,:,:] = np.nanmean(Var_spco2,axis=0) * ( np.exp( 0.0423 * ( Var_sst[ii,:,:] - np.nanmean(Var_sst,axis=0) ) ) );
                Var_spco2_NonTEMP_eq2[ii,:,:] = Var_spco2[ii,:,:] * ( np.exp( 0.0423 * ( np.nanmean(Var_sst,axis=0) - Var_sst[ii,:,:] ) ) );
            
        else: # If we don't have both datasets for this model to do the correlation, just plot an empty map, to keep the order of GCMs constant in all models  
            Var_spco2_TEMP_eq1=empty((yrs_n*12,lat_n_regrid,lon_n_regrid)) * nan             
            Var_spco2_NonTEMP_eq2=empty((yrs_n*12,lat_n_regrid,lon_n_regrid)) * nan


        if 'SST_allmonths' in my_shelf_sst: # If that variable exists in the saved data  
            globals()['Var_sst']=my_shelf_sst['SST_allmonths']                    
            Var_sst[ np.abs(Var_sst)>1e19 ]=nan  # Replacing missing values with NaN 
            
        else: # If we don't have both datasets for this model to do the correlation, just plot an empty map, to keep the order of GCMs constant in all models  
            Var_sst=empty((yrs_n*12,lat_n_regrid,lon_n_regrid)) * nan             


        if V_i == 0 or V_i == 1 or V_i == 3:
            
            filename_in_x = (dir_data_in + 'AllResults_'+GCM+'_'+'bio'+'_1991_2010.out') # Directory to save processed data
            my_shelf_x = shelve.open(filename_in_x)
     
            if 'NO3_allmonths' in my_shelf_x: # If that variable exists in the saved data                  
                globals()['Variable_x']=my_shelf_x['NO3_allmonths']                        
                Variable_x[ np.abs(Variable_x)>1e19 ]=nan  # Replacing missing values with NaN 
                
            else:
                Variable_x=empty((yrs_n*12,lat_n_regrid,lon_n_regrid)) * nan  
            
        else:
            
            filename_in_x = (dir_data_in + 'AllResults_'+GCM+'_'+'bio'+'_1991_2010.out') # Directory to save processed data
            my_shelf_x = shelve.open(filename_in_x)
     
            if 'PPint_allmonths' in my_shelf_x: # If that variable exists in the saved data                  
                globals()['Variable_x']=my_shelf_x['PPint_allmonths']                        
                Variable_x[ np.abs(Variable_x)>1e19 ]=nan  # Replacing missing values with NaN 
                
            else:
                Variable_x=empty((yrs_n*12,lat_n_regrid,lon_n_regrid)) * nan 


        Variable_x_PAPA = Variable_x[:,135:140,210:220]
        Var_spco2_TEMP_eq1_PAPA = Var_spco2_TEMP_eq1[:,135:140,210:220]
        Var_spco2_NonTEMP_eq2_PAPA = Var_spco2_NonTEMP_eq2[:,135:140,210:220]
        Var_sst_PAPA = Var_sst[:,135:140,210:220]

        Variable_x_NABE = Variable_x[:,135:140,325:335]
        Var_spco2_TEMP_eq1_NABE = Var_spco2_TEMP_eq1[:,135:140,325:335]
        Var_spco2_NonTEMP_eq2_NABE = Var_spco2_NonTEMP_eq2[:,135:140,325:335]
        Var_sst_NABE = Var_sst[:,135:140,325:335]
                
        if d_trend=='yes': # If it is requested to detrend the data before correlations
            Variable_x_PAPA= func_detrend_3d(Variable_x_PAPA)
            Var_spco2_TEMP_eq1_PAPA= func_detrend_3d(Var_spco2_TEMP_eq1_PAPA) 
            Var_spco2_NonTEMP_eq2_PAPA= func_detrend_3d(Var_spco2_NonTEMP_eq2_PAPA) 
            Var_sst_PAPA= func_detrend_3d(Var_sst_PAPA) 
            Variable_x_NABE= func_detrend_3d(Variable_x_NABE)
            Var_spco2_TEMP_eq1_NABE= func_detrend_3d(Var_spco2_TEMP_eq1_NABE)
            Var_spco2_NonTEMP_eq2_NABE= func_detrend_3d(Var_spco2_NonTEMP_eq2_NABE)
            Var_sst_NABE= func_detrend_3d(Var_sst_NABE) 

        Variable_x_PAPA_reshaped = Variable_x_PAPA.reshape(( yrs_n ,12,Variable_x_PAPA.shape[1],Variable_x_PAPA.shape[2])) # Monthly means of SpCO2 over the time period  
        Var_spco2_TEMP_eq1_PAPA_reshaped = Var_spco2_TEMP_eq1_PAPA.reshape(( yrs_n ,12,Var_spco2_TEMP_eq1_PAPA.shape[1],Var_spco2_TEMP_eq1_PAPA.shape[2])) # Monthly means of SpCO2 over the time period 
        Var_spco2_NonTEMP_eq2_PAPA_reshaped = Var_spco2_NonTEMP_eq2_PAPA.reshape(( yrs_n ,12,Var_spco2_NonTEMP_eq2_PAPA.shape[1],Var_spco2_NonTEMP_eq2_PAPA.shape[2])) # Monthly means of SpCO2 over the time period 
        Var_sst_PAPA_reshaped = Var_sst_PAPA.reshape(( yrs_n ,12,Var_sst_PAPA.shape[1],Var_sst_PAPA.shape[2])) # Monthly means of SpCO2 over the time period 
  
        Variable_x_NABE_reshaped = Variable_x_NABE.reshape(( yrs_n ,12,Variable_x_NABE.shape[1],Variable_x_NABE.shape[2])) # Monthly means of SpCO2 over the time period  
        Var_spco2_TEMP_eq1_NABE_reshaped = Var_spco2_TEMP_eq1_NABE.reshape(( yrs_n ,12,Var_spco2_TEMP_eq1_NABE.shape[1],Var_spco2_TEMP_eq1_NABE.shape[2])) # Monthly means of SpCO2 over the time period 
        Var_spco2_NonTEMP_eq2_NABE_reshaped = Var_spco2_NonTEMP_eq2_NABE.reshape(( yrs_n ,12,Var_spco2_NonTEMP_eq2_NABE.shape[1],Var_spco2_NonTEMP_eq2_NABE.shape[2])) # Monthly means of SpCO2 over the time period 
        Var_sst_NABE_reshaped = Var_sst_NABE.reshape(( yrs_n ,12,Var_sst_NABE.shape[1],Var_sst_NABE.shape[2])) # Monthly means of SpCO2 over the time period 
        
        Var_spco2_TEMP_eq1_PAPA_reshaped= np.nanmean(Var_spco2_TEMP_eq1_PAPA_reshaped,axis=0)
        Var_spco2_NonTEMP_eq2_PAPA_reshaped= np.nanmean(Var_spco2_NonTEMP_eq2_PAPA_reshaped,axis=0)
        Variable_x_PAPA_reshaped= np.nanmean(Variable_x_PAPA_reshaped,axis=0)
        Var_sst_PAPA_reshaped= np.nanmean(Var_sst_PAPA_reshaped,axis=0)
        
        Var_spco2_TEMP_eq1_NABE_reshaped= np.nanmean(Var_spco2_TEMP_eq1_NABE_reshaped,axis=0)
        Var_spco2_NonTEMP_eq2_NABE_reshaped= np.nanmean(Var_spco2_NonTEMP_eq2_NABE_reshaped,axis=0)
        Variable_x_NABE_reshaped= np.nanmean(Variable_x_NABE_reshaped,axis=0)
        Var_sst_NABE_reshaped= np.nanmean(Var_sst_NABE_reshaped,axis=0)
        
        
        if V_i == 0:
            
            Var_plot_x_PAPA = np.nanmean(np.nanmean(  np.nanmean(Variable_x_PAPA_reshaped[5:8,:,:], axis=0)  , axis=0) , axis=0) #NO3 (JJA average)
            Var_plot_x_NABE = np.nanmean(np.nanmean(  np.nanmean(Variable_x_NABE_reshaped[5:8,:,:], axis=0)  , axis=0) , axis=0)
            
            Var_plot_y_PAPA = np.nanmean(np.nanmean(  np.nanmean(Var_spco2_NonTEMP_eq2_PAPA_reshaped[5:8,:,:], axis=0)  , axis=0) , axis=0) #pco2-NonT (JJA average)
            Var_plot_y_NABE = np.nanmean(np.nanmean(  np.nanmean(Var_spco2_NonTEMP_eq2_NABE_reshaped[5:8,:,:], axis=0)  , axis=0) , axis=0)
        
        if V_i == 1:
            
            Var_plot_x_PAPA = np.nanmean(np.nanmean(  np.nanmax(Variable_x_PAPA_reshaped,axis=0) - np.nanmin(Variable_x_PAPA_reshaped,axis=0)  , axis=0) , axis=0) #NO3 (max - min)
            Var_plot_x_NABE = np.nanmean(np.nanmean(  np.nanmax(Variable_x_NABE_reshaped,axis=0) - np.nanmin(Variable_x_NABE_reshaped,axis=0)  , axis=0) , axis=0)
            
            Var_plot_y_PAPA = np.nanmean(np.nanmean(  np.nanmax(Var_spco2_NonTEMP_eq2_PAPA_reshaped,axis=0) - np.nanmin(Var_spco2_NonTEMP_eq2_PAPA_reshaped,axis=0)  , axis=0) , axis=0) #pco2-NonT (max - min)
            Var_plot_y_NABE = np.nanmean(np.nanmean(  np.nanmax(Var_spco2_NonTEMP_eq2_NABE_reshaped,axis=0) - np.nanmin(Var_spco2_NonTEMP_eq2_NABE_reshaped,axis=0)  , axis=0) , axis=0)
        
        if V_i == 2:
            
            Var_plot_x_PAPA = np.nanmean(np.nanmean(  np.nanmax(Variable_x_PAPA_reshaped,axis=0) - np.nanmin(Variable_x_PAPA_reshaped,axis=0)  , axis=0) , axis=0) #PP_int (max - min)
            Var_plot_x_NABE = np.nanmean(np.nanmean(  np.nanmax(Variable_x_NABE_reshaped,axis=0) - np.nanmin(Variable_x_NABE_reshaped,axis=0)  , axis=0) , axis=0)
            
            Var_plot_y_PAPA = np.nanmean(np.nanmean(  np.nanmax(Var_spco2_NonTEMP_eq2_PAPA_reshaped,axis=0) - np.nanmin(Var_spco2_NonTEMP_eq2_PAPA_reshaped,axis=0)  , axis=0) , axis=0) #pco2-NonT (max - min)
            Var_plot_y_NABE = np.nanmean(np.nanmean(  np.nanmax(Var_spco2_NonTEMP_eq2_NABE_reshaped,axis=0) - np.nanmin(Var_spco2_NonTEMP_eq2_NABE_reshaped,axis=0)  , axis=0) , axis=0)
        
        if V_i == 3:
            # (pCO2-T(max)-pCO2-T(min))/(pCO2-nonT(max)-pCO2-nonT(min))            
            Var_plot_x_PAPA = np.nanmean(np.nanmean(  ((np.nanmax(Var_spco2_TEMP_eq1_PAPA_reshaped,axis=0) - np.nanmin(Var_spco2_TEMP_eq1_PAPA_reshaped,axis=0))/(np.nanmax(Var_spco2_NonTEMP_eq2_PAPA_reshaped,axis=0) - np.nanmin(Var_spco2_NonTEMP_eq2_PAPA_reshaped,axis=0)))  , axis=0) , axis=0)
            Var_plot_x_NABE = np.nanmean(np.nanmean(  ((np.nanmax(Var_spco2_TEMP_eq1_NABE_reshaped,axis=0) - np.nanmin(Var_spco2_TEMP_eq1_NABE_reshaped,axis=0))/(np.nanmax(Var_spco2_NonTEMP_eq2_NABE_reshaped,axis=0) - np.nanmin(Var_spco2_NonTEMP_eq2_NABE_reshaped,axis=0)))  , axis=0) , axis=0)

            # (SST(max)-SST(min))/(NO3(max)-NO3(min)) 
            Var_plot_y_PAPA = np.nanmean(np.nanmean(  ((np.nanmax(Var_sst_PAPA_reshaped,axis=0) - np.nanmin(Var_sst_PAPA_reshaped,axis=0))/(np.nanmax(Variable_x_PAPA_reshaped,axis=0) - np.nanmin(Variable_x_PAPA_reshaped,axis=0)))  , axis=0) , axis=0)
            Var_plot_y_NABE = np.nanmean(np.nanmean(  ((np.nanmax(Var_sst_NABE_reshaped,axis=0) - np.nanmin(Var_sst_NABE_reshaped,axis=0))/(np.nanmax(Variable_x_NABE_reshaped,axis=0) - np.nanmin(Variable_x_NABE_reshaped,axis=0)))  , axis=0) , axis=0)
              
        
        plt.scatter(Var_plot_x_PAPA, Var_plot_y_PAPA, s=100, marker='o', c=Colors[M_i], label=GCM) # Only the first pannel gets the lables so the legend location be alligned with that
        plt.scatter(Var_plot_x_PAPA, Var_plot_y_PAPA, s=100, marker='o', facecolors='none', edgecolors='k')

        plt.scatter(Var_plot_x_NABE, Var_plot_y_NABE, s=100, marker='d', c=Colors[M_i], label=GCM) # Only the first pannel gets the lables so the legend location be alligned with that
        plt.scatter(Var_plot_x_NABE, Var_plot_y_NABE, s=100, marker='d', facecolors='none', edgecolors='k')
        
        plt.xlabel(var_names_x[V_i],fontsize=14)
        plt.ylabel(var_names_y[V_i],fontsize=14)  

        plt.xlim(plot_ranges_x[V_i][0],plot_ranges_x[V_i][1])
        plt.ylim(plot_ranges_y[V_i][0],plot_ranges_y[V_i][1])
        plt.xticks(fontsize=10); plt.yticks(fontsize = 10) # Setting the X-ticks ## rotation='horizontal' # rotation='vertical'   
        if  np.abs( plot_ranges_x[V_i][1] ) < 0.001 :
            plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))  
        if  np.abs( plot_ranges_y[V_i][1] ) < 0.001 :
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))  
    
plt.subplots_adjust(left=0.1, bottom=0.1, right=0.8, top=0.92, hspace=0.35, wspace=0.35) # the amount of height/width reserved for space between subplots
#plt.legend(shadow=True, loc='right', bbox_to_anchor=(2.1, 2.5), markerscale=1 , prop={'size': 14})
if d_trend=='yes':
    plt.suptitle('CMIP5 models - North Pacific (circles) and North Atlantic (diamonds) - '+str(start_date_plt)+'-'+str(end_date_plt)+' - detrended' , fontsize=18)
else:
    plt.suptitle('CMIP5 models - North Pacific (circles) and North Atlantic (diamonds) - '+str(start_date_plt)+'-'+str(end_date_plt) , fontsize=18)

mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full

#if d_trend=='yes':
#    fig.savefig(dir_figs+'AllModels_Scatter_detrended_3.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
#else:
#    fig.savefig(dir_figs+'AllModels_Scatter_3.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
#


#####################################################
#####################################################
dir_data_in_1 = (dir_pwd + '/CMIP6_TXT/') # Directory to save results

GCM_Names = ['CESM2', 'GFDL-CM4', 'GFDL-ESM4', 'CESM2-WACCM', 'IPSL-CM6A-LR', 'MRI-ESM2-0','GISS-E2-1-G', 'GISS-E2-1-G-CC', 'UKESM1-0-LL', 'CanESM5', 'NorESM2-LM', 'MPI-ESM1-2-HR', 'NorCPM1', 'MIROC-ES2L', 'CNRM-ESM2-1']
Colors = [[1, 1, 0], [1, 0.8, 0], [1, 0.6, 0], [1, 0.4, 0], [1, 0, 0], [1, 0.8, 1], [1, 0.5, 1], [1, 0, 1], [0, 0, 0.5], [0, 0, 1], [0, 0.5, 1], [0, 0.5, 0], [0, 1, 0], [0, 0.6, 0], [0.5, 0.2, 0]]

#plot_ranges_x=[[0,12], [0,20], [0,20],  [0,4],   [0,600],  [0,6],  [0,6],  [0,8],  [0,8],  [0,8],  [0,8],  [0,8],  [0,8],  [0,8] ]
#plot_ranges_y=[[0,15], [0,20], [0,600], [0,600], [0,20],   [0,8],  [0,20], [0,8],  [0,8],  [0,8],  [0,8],  [0,8],  [0,8],  [0,8] ]

plot_ranges_x=[[-0.001,0.008],  [-0.001,0.011], [-0.001,0.011], [-1E-6,2E-5],  [-1E-7,2E-6],   [-0.0001,0.0014], [-0.0001,0.0008], [0,5], [0,5], [0,5],  [-0.001,0.018],  [-0.001,0.018],  [-0.001,0.018],  [-0.001,0.011] ]
plot_ranges_y=[[-0.0001,0.004], [-0.001,0.018], [-1E-7,2E-6],   [-1E-7,2E-6],  [-0.001,0.018], [-0.0001,0.004],  [-0.001,0.018],   [0,5], [0,5], [0,5],  [-0.0001,0.0008],  [-0.001,0.005],  [-1E-7,2E-6],  [-0.001,0.018] ]

var_label_x = ['(ave)nitrate','(max)nitrate','(max)nitrate','(max)fe','(max)pp','(ave)zoo','(max)chlo','(max-min/ave)phyto','(max-min/ave)phyto','(max-min/ave)phyto', '(max)phyto','(max)phyto','(max)phyto','(max)nitrate']
var_label_y = ['(ave)phyto','(max)phyto','(max)pp','(max)pp','(max)phyto','(ave)phyto','(max)phyto','(max-min/ave)chlo','(max-min/ave)zoo','(max-min/ave)pp','(max)chlo','(max)zoo','(max)pp','(max)phyto']

var_names_x = ['no3','no3','no3','dfe','intpp','zooc','chl','phyc','phyc','phyc', 'phyc','phyc','phyc','no3']
var_names_y = ['phyc','phyc','intpp','intpp','phyc','phyc','phyc','chl','zooc','intpp','chl','zooc','intpp','phyc']

var_operator_x = ['ave','max','max','max','max','ave','max','norm','norm','norm', 'max','max','max','max']
var_operator_y = ['ave','max','max','max','max','ave','max','norm','norm','norm','max','max','max','max']

############################################################
############################################################

#var_names_lst=['phyc','chl','phydiat','zooc','intpp','epc100','fgco2','spco2','dfe','no3','si','rsntds']

n_r=4 ; n_c=5 ;
n_t=len(var_names_x)

fig=plt.figure()
for V_i in range(n_t): # V_i=0   V_i=1    V_i=2    V_i=3    V_i=4    V_i=5    V_i=6    V_i=7    V_i=8    V_i=9    V_i=10    V_i=11    V_i=12    V_i=13

    ax = fig.add_subplot(n_r,n_c,V_i+1) 

    for M_i in range(len(GCM_Names)): # M_i=0   M_i=6    # GCM=GCM_Names[0]

        GCM=GCM_Names[M_i]
        print (GCM) 
        
        dir_x_PAPA = (dir_data_in_1 + GCM + '_' + var_names_x[V_i] + '_np.txt' )
        dir_y_PAPA = (dir_data_in_1 + GCM + '_' + var_names_y[V_i] + '_np.txt' )
        if os.path.isfile(dir_x_PAPA) and os.path.isfile(dir_y_PAPA):
            Var_x_PAPA = np.array( open(dir_x_PAPA, "r").readlines() , dtype='float64')
            Var_y_PAPA = np.array( open(dir_y_PAPA, "r").readlines() , dtype='float64')
            if np.nanmean(Var_x_PAPA) < 0:
                Var_x_PAPA = Var_x_PAPA * -1
            if np.nanmean(Var_y_PAPA) < 0:
                Var_y_PAPA = Var_y_PAPA * -1            
            
            if var_operator_x[V_i] == 'ave':
                Var_plot_x_PAPA = np.nanmean(Var_x_PAPA)
            elif var_operator_x[V_i] == 'max':
                Var_plot_x_PAPA = np.nanmax(Var_x_PAPA)
            elif var_operator_x[V_i] == 'min':
                Var_plot_x_PAPA = np.nanmin(Var_x_PAPA)
            elif var_operator_x[V_i] == 'norm':  
                Var_plot_x_PAPA = (np.nanmax(Var_x_PAPA) - np.nanmin(Var_x_PAPA)) / np.nanmean(Var_x_PAPA)
            else:
                Var_plot_x_PAPA = np.nan

            if var_operator_y[V_i] == 'ave':
                Var_plot_y_PAPA = np.nanmean(Var_y_PAPA)
            elif var_operator_y[V_i] == 'max':
                Var_plot_y_PAPA = np.nanmax(Var_y_PAPA)
            elif var_operator_y[V_i] == 'min':
                Var_plot_y_PAPA = np.nanmin(Var_y_PAPA)
            elif var_operator_x[V_i] == 'norm':  
                Var_plot_y_PAPA = (np.nanmax(Var_y_PAPA) - np.nanmin(Var_y_PAPA)) / np.nanmean(Var_y_PAPA)
            else:
                Var_plot_y_PAPA = np.nan
            
        else:
            Var_plot_x_PAPA = np.nan
            Var_plot_y_PAPA = np.nan


        dir_x_NABE = (dir_data_in_1 + GCM + '_' + var_names_x[V_i] + '_na.txt' )
        dir_y_NABE = (dir_data_in_1 + GCM + '_' + var_names_y[V_i] + '_na.txt' ) 
        if os.path.isfile(dir_x_NABE) and os.path.isfile(dir_y_NABE):
            Var_x_NABE = np.array( open(dir_x_NABE, "r").readlines() , dtype='float64')
            Var_y_NABE = np.array( open(dir_y_NABE, "r").readlines() , dtype='float64')
            if np.nanmean(Var_x_NABE) < 0:
                Var_x_NABE = Var_x_NABE * -1
            if np.nanmean(Var_y_NABE) < 0:
                Var_y_NABE = Var_y_NABE * -1  
                
            if var_operator_x[V_i] == 'ave':
                Var_plot_x_NABE = np.nanmean(Var_x_NABE)
            elif var_operator_x[V_i] == 'max':
                Var_plot_x_NABE = np.nanmax(Var_x_NABE)
            elif var_operator_x[V_i] == 'min':
                Var_plot_x_NABE = np.nanmin(Var_x_NABE)
            elif var_operator_x[V_i] == 'norm':  
                Var_plot_x_NABE = (np.nanmax(Var_x_NABE) - np.nanmin(Var_x_NABE)) / np.nanmean(Var_x_NABE)
            else:
                Var_plot_x_NABE = np.nan

            if var_operator_y[V_i] == 'ave':
                Var_plot_y_NABE = np.nanmean(Var_y_NABE)
            elif var_operator_y[V_i] == 'max':
                Var_plot_y_NABE = np.nanmax(Var_y_NABE)
            elif var_operator_y[V_i] == 'min':
                Var_plot_y_NABE = np.nanmin(Var_y_NABE)
            elif var_operator_x[V_i] == 'norm':  
                Var_plot_y_NABE = (np.nanmax(Var_y_NABE) - np.nanmin(Var_y_NABE)) / np.nanmean(Var_y_NABE)
            else:
                Var_plot_y_NABE = np.nan
            
        else:
            Var_plot_x_NABE = np.nan
            Var_plot_y_NABE = np.nan            

        plt.scatter(Var_plot_x_PAPA, Var_plot_y_PAPA, s=100, marker='o', c=Colors[M_i], label=GCM) # Only the first pannel gets the lables so the legend location be alligned with that
        plt.scatter(Var_plot_x_PAPA, Var_plot_y_PAPA, s=100, marker='o', facecolors='none', edgecolors='k')

        plt.scatter(Var_plot_x_NABE, Var_plot_y_NABE, s=100, marker='d', c=Colors[M_i], label=GCM) # Only the first pannel gets the lables so the legend location be alligned with that
        plt.scatter(Var_plot_x_NABE, Var_plot_y_NABE, s=100, marker='d', facecolors='none', edgecolors='k')
        
#        if V_i >= 7:
#            plt.plot([-1000,1000],[-1000,1000], 'k--', linewidth=0.75)
        
        plt.xlabel(var_label_x[V_i],fontsize=14)
        plt.ylabel(var_label_y[V_i],fontsize=14)  

        plt.xlim(plot_ranges_x[V_i][0],plot_ranges_x[V_i][1])
        plt.ylim(plot_ranges_y[V_i][0],plot_ranges_y[V_i][1])
        plt.xticks(fontsize=10); plt.yticks(fontsize = 10) # Setting the X-ticks ## rotation='horizontal' # rotation='vertical'   
        if  np.abs( plot_ranges_x[V_i][1] ) < 0.001 :
            plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))  
        if  np.abs( plot_ranges_y[V_i][1] ) < 0.001 :
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))  
    
plt.subplots_adjust(left=0.05, bottom=0.1, right=0.85, top=0.92, hspace=0.35, wspace=0.42) # the amount of height/width reserved for space between subplots
#plt.legend(shadow=True, loc='right', bbox_to_anchor=(2.1, 2.5), markerscale=1 , prop={'size': 14})
if d_trend=='yes':
    plt.suptitle('CMIP5 models - North Pacific (circles) and North Atlantic (diamonds) - '+str(start_date_plt)+'-'+str(end_date_plt)+' - detrended' , fontsize=18)
else:
    plt.suptitle('CMIP5 models - North Pacific (circles) and North Atlantic (diamonds) - '+str(start_date_plt)+'-'+str(end_date_plt) , fontsize=18)

mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full

#if d_trend=='yes':
#    fig.savefig(dir_figs+'AllModels_Scatter_detrended_4.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
#else:
#    fig.savefig(dir_figs+'AllModels_Scatter_4.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


#####################################################
#####################################################
dir_data_in_1 = (dir_pwd + '/CMIP6_TXT/') # Directory to save results

GCM_Names = ['CESM2', 'GFDL-CM4', 'GFDL-ESM4', 'CESM2-WACCM', 'IPSL-CM6A-LR', 'MRI-ESM2-0','GISS-E2-1-G', 'GISS-E2-1-G-CC', 'UKESM1-0-LL', 'CanESM5', 'NorESM2-LM', 'MPI-ESM1-2-HR', 'NorCPM1', 'MIROC-ES2L', 'CNRM-ESM2-1']
Colors = [[1, 1, 0], [1, 0.8, 0], [1, 0.6, 0], [1, 0.4, 0], [1, 0, 0], [1, 0.8, 1], [1, 0.5, 1], [1, 0, 1], [0, 0, 0.5], [0, 0, 1], [0, 0.5, 1], [0, 0.5, 0], [0, 1, 0], [0, 0.6, 0], [0.5, 0.2, 0]]

plot_ranges_x=[[-0.001,0.008],  [-0.001,0.011], [-0.001,0.011], [-1E-6,2E-5],  [-1E-7,2E-6],   [-0.0001,0.0014], [-0.0001,0.0008], [0,5], [0,5], [0,5],  [-0.001,0.018],  [-0.001,0.018],  [-0.001,0.018],  [-0.001,0.011],  [-0.001,0.005],  [-0.00005,0.0008], [-5E-8,2E-7],    [-5E-8,2E-7],      [-0.003,0.011], [-1E-8,3E-7] ]
plot_ranges_y=[[-0.0001,0.004], [-0.001,0.018], [-1E-7,2E-6],   [-1E-7,2E-6],  [-0.001,0.018], [-0.0001,0.004],  [-0.001,0.018],   [0,5], [0,5], [0,5],  [-0.0001,0.0008],  [-0.001,0.005],  [-1E-7,2E-6],  [-0.001,0.018],  [-0.001,0.018],  [-0.0001,0.002],   [-0.0001,0.002], [-0.00001,0.0008], [-0.003,0.02],  [-0.03,0.8]  ]

var_label_x = ['(ave)nitrate','(max)nitrate','(max)nitrate','(max)fe','(max)pp','(ave)zoo','(max)chlo','(max-min/ave)phyto','(max-min/ave)phyto','(max-min/ave)phyto', '(max)phyto','(max)phyto','(max)phyto','(max)nitrate','(max)zoo','(min)zoo','(min)pp','(min)pp','(max)nitrate','(min)fe']
var_label_y = ['(ave)phyto','(max)phyto','(max)pp','(max)pp','(max)phyto','(ave)phyto','(max)phyto','(max-min/ave)chlo','(max-min/ave)zoo','(max-min/ave)pp','(max)chlo','(max)zoo','(max)pp','(max)phyto','(max)phyto','(min)phyto','(min)phyto','(min)zoo','(max)(phyto+zoo)s','min_phyto / min_zoo']

var_names_x = ['no3','no3','no3','dfe','intpp','zooc','chl','phyc','phyc','phyc', 'phyc','phyc','phyc','no3','zooc','zooc','intpp','intpp','no3','dfe']
var_names_y = ['phyc','phyc','intpp','intpp','phyc','phyc','phyc','chl','zooc','intpp','chl','zooc','intpp','phyc','phyc','phyc','phyc','zooc','phyc','phyc']

var_operator_x = ['ave','max','max','max','max','ave','max','norm','norm','norm', 'max','max','max','max','max','min','min','min','max','min']
var_operator_y = ['ave','max','max','max','max','ave','max','norm','norm','norm','max','max','max','max','max','min','min','min','max','min']


n_r=4 ; n_c=5 ;
n_t=len(var_names_x)

fig=plt.figure()
for V_i in range(n_t): # V_i=0   V_i=1    V_i=2    V_i=3    V_i=4    V_i=5    V_i=6    V_i=7    V_i=8    V_i=9    V_i=10    V_i=11    V_i=12    V_i=13

    ax = fig.add_subplot(n_r,n_c,V_i+1) 

    for M_i in range(len(GCM_Names)): # M_i=0   M_i=6    # GCM=GCM_Names[0]

        GCM=GCM_Names[M_i]
        print (GCM) 
        
        if V_i == 18:
            dir_x_PAPA = (dir_data_in_1 + GCM + '_' + 'no3' + '_np.txt' )
            dir_y1_PAPA = (dir_data_in_1 + GCM + '_' + 'phyc' + '_np.txt' )  
            dir_y2_PAPA = (dir_data_in_1 + GCM + '_' + 'zooc' + '_np.txt' )
            if os.path.isfile(dir_x_PAPA) and os.path.isfile(dir_y1_PAPA) and os.path.isfile(dir_y2_PAPA):
                Var_x_PAPA = np.array( open(dir_x_PAPA, "r").readlines() , dtype='float64')
                Var_y1_PAPA = np.array( open(dir_y1_PAPA, "r").readlines() , dtype='float64')
                Var_y2_PAPA = np.array( open(dir_y2_PAPA, "r").readlines() , dtype='float64')
                
                Var_plot_x_PAPA = np.nanmax(Var_x_PAPA)
                Var_plot_y_PAPA = np.nanmax(Var_y1_PAPA+Var_y2_PAPA)

            else:
                Var_plot_x_PAPA = np.nan
                Var_plot_y_PAPA = np.nan        
        
            dir_x_NABE = (dir_data_in_1 + GCM + '_' + 'no3' + '_na.txt' )
            dir_y1_NABE = (dir_data_in_1 + GCM + '_' + 'phyc' + '_na.txt' )  
            dir_y2_NABE = (dir_data_in_1 + GCM + '_' + 'zooc' + '_na.txt' )
            if os.path.isfile(dir_x_NABE) and os.path.isfile(dir_y1_NABE)and os.path.isfile(dir_y2_NABE):
                Var_x_NABE = np.array( open(dir_x_NABE, "r").readlines() , dtype='float64')
                Var_y1_NABE = np.array( open(dir_y1_NABE, "r").readlines() , dtype='float64')
                Var_y2_NABE = np.array( open(dir_y2_NABE, "r").readlines() , dtype='float64')
                
                Var_plot_x_NABE = np.nanmax(Var_x_NABE)
                Var_plot_y_NABE = np.nanmax(Var_y1_NABE+Var_y2_NABE)

            else:
                Var_plot_x_NABE = np.nan
                Var_plot_y_NABE = np.nan         


        elif V_i == 19:
            dir_x_PAPA = (dir_data_in_1 + GCM + '_' + 'dfe' + '_np.txt' )
            dir_y1_PAPA = (dir_data_in_1 + GCM + '_' + 'zooc' + '_np.txt' )  
            dir_y2_PAPA = (dir_data_in_1 + GCM + '_' + 'phyc' + '_np.txt' )
            if os.path.isfile(dir_x_PAPA) and os.path.isfile(dir_y1_PAPA) and os.path.isfile(dir_y2_PAPA):
                Var_x_PAPA = np.array( open(dir_x_PAPA, "r").readlines() , dtype='float64')
                Var_y1_PAPA = np.array( open(dir_y1_PAPA, "r").readlines() , dtype='float64')
                Var_y2_PAPA = np.array( open(dir_y2_PAPA, "r").readlines() , dtype='float64')
                
                Var_plot_x_PAPA = np.nanmin(Var_x_PAPA)
                Var_plot_y_PAPA = np.nanmin(Var_y1_PAPA)/np.nanmin(Var_y2_PAPA)

            else:
                Var_plot_x_PAPA = np.nan
                Var_plot_y_PAPA = np.nan

            dir_x_NABE = (dir_data_in_1 + GCM + '_' + 'dfe' + '_na.txt' )
            dir_y1_NABE = (dir_data_in_1 + GCM + '_' + 'zooc' + '_na.txt' )  
            dir_y2_NABE = (dir_data_in_1 + GCM + '_' + 'phyc' + '_na.txt' )
            if os.path.isfile(dir_x_NABE) and os.path.isfile(dir_y1_NABE) and os.path.isfile(dir_y2_NABE):
                Var_x_NABE = np.array( open(dir_x_NABE, "r").readlines() , dtype='float64')
                Var_y1_NABE = np.array( open(dir_y1_NABE, "r").readlines() , dtype='float64')
                Var_y2_NABE = np.array( open(dir_y2_NABE, "r").readlines() , dtype='float64')
                
                Var_plot_x_NABE = np.nanmin(Var_x_NABE)
                Var_plot_y_NABE = np.nanmin(Var_y1_NABE)/np.nanmin(Var_y2_NABE)

            else:
                Var_plot_x_NABE = np.nan
                Var_plot_y_NABE = np.nan
        
        else:
            dir_x_PAPA = (dir_data_in_1 + GCM + '_' + var_names_x[V_i] + '_np.txt' )
            dir_y_PAPA = (dir_data_in_1 + GCM + '_' + var_names_y[V_i] + '_np.txt' )
            if os.path.isfile(dir_x_PAPA) and os.path.isfile(dir_y_PAPA):
                Var_x_PAPA = np.array( open(dir_x_PAPA, "r").readlines() , dtype='float64')
                Var_y_PAPA = np.array( open(dir_y_PAPA, "r").readlines() , dtype='float64')
                if np.nanmean(Var_x_PAPA) < 0:
                    Var_x_PAPA = Var_x_PAPA * -1
                if np.nanmean(Var_y_PAPA) < 0:
                    Var_y_PAPA = Var_y_PAPA * -1            
                
                if var_operator_x[V_i] == 'ave':
                    Var_plot_x_PAPA = np.nanmean(Var_x_PAPA)
                elif var_operator_x[V_i] == 'max':
                    Var_plot_x_PAPA = np.nanmax(Var_x_PAPA)
                elif var_operator_x[V_i] == 'min':
                    Var_plot_x_PAPA = np.nanmin(Var_x_PAPA)
                elif var_operator_x[V_i] == 'norm':  
                    Var_plot_x_PAPA = (np.nanmax(Var_x_PAPA) - np.nanmin(Var_x_PAPA)) / np.nanmean(Var_x_PAPA)
                else:
                    Var_plot_x_PAPA = np.nan
    
                if var_operator_y[V_i] == 'ave':
                    Var_plot_y_PAPA = np.nanmean(Var_y_PAPA)
                elif var_operator_y[V_i] == 'max':
                    Var_plot_y_PAPA = np.nanmax(Var_y_PAPA)
                elif var_operator_y[V_i] == 'min':
                    Var_plot_y_PAPA = np.nanmin(Var_y_PAPA)
                elif var_operator_x[V_i] == 'norm':  
                    Var_plot_y_PAPA = (np.nanmax(Var_y_PAPA) - np.nanmin(Var_y_PAPA)) / np.nanmean(Var_y_PAPA)
                else:
                    Var_plot_y_PAPA = np.nan
                
            else:
                Var_plot_x_PAPA = np.nan
                Var_plot_y_PAPA = np.nan
    
    
            dir_x_NABE = (dir_data_in_1 + GCM + '_' + var_names_x[V_i] + '_na.txt' )
            dir_y_NABE = (dir_data_in_1 + GCM + '_' + var_names_y[V_i] + '_na.txt' ) 
            if os.path.isfile(dir_x_NABE) and os.path.isfile(dir_y_NABE):
                Var_x_NABE = np.array( open(dir_x_NABE, "r").readlines() , dtype='float64')
                Var_y_NABE = np.array( open(dir_y_NABE, "r").readlines() , dtype='float64')
                if np.nanmean(Var_x_NABE) < 0:
                    Var_x_NABE = Var_x_NABE * -1
                if np.nanmean(Var_y_NABE) < 0:
                    Var_y_NABE = Var_y_NABE * -1  
                    
                if var_operator_x[V_i] == 'ave':
                    Var_plot_x_NABE = np.nanmean(Var_x_NABE)
                elif var_operator_x[V_i] == 'max':
                    Var_plot_x_NABE = np.nanmax(Var_x_NABE)
                elif var_operator_x[V_i] == 'min':
                    Var_plot_x_NABE = np.nanmin(Var_x_NABE)
                elif var_operator_x[V_i] == 'norm':  
                    Var_plot_x_NABE = (np.nanmax(Var_x_NABE) - np.nanmin(Var_x_NABE)) / np.nanmean(Var_x_NABE)
                else:
                    Var_plot_x_NABE = np.nan
    
                if var_operator_y[V_i] == 'ave':
                    Var_plot_y_NABE = np.nanmean(Var_y_NABE)
                elif var_operator_y[V_i] == 'max':
                    Var_plot_y_NABE = np.nanmax(Var_y_NABE)
                elif var_operator_y[V_i] == 'min':
                    Var_plot_y_NABE = np.nanmin(Var_y_NABE)
                elif var_operator_x[V_i] == 'norm':  
                    Var_plot_y_NABE = (np.nanmax(Var_y_NABE) - np.nanmin(Var_y_NABE)) / np.nanmean(Var_y_NABE)
                else:
                    Var_plot_y_NABE = np.nan
                
            else:
                Var_plot_x_NABE = np.nan
                Var_plot_y_NABE = np.nan            

        plt.scatter(Var_plot_x_PAPA, Var_plot_y_PAPA, s=100, marker='o', c=Colors[M_i], label=GCM) # Only the first pannel gets the lables so the legend location be alligned with that
        plt.scatter(Var_plot_x_PAPA, Var_plot_y_PAPA, s=100, marker='o', facecolors='none', edgecolors='k')

        plt.scatter(Var_plot_x_NABE, Var_plot_y_NABE, s=100, marker='d', c=Colors[M_i], label=GCM) # Only the first pannel gets the lables so the legend location be alligned with that
        plt.scatter(Var_plot_x_NABE, Var_plot_y_NABE, s=100, marker='d', facecolors='none', edgecolors='k')
        
#        if V_i >= 7:
#            plt.plot([-1000,1000],[-1000,1000], 'k--', linewidth=0.75)
        
        plt.xlabel(var_label_x[V_i],fontsize=14)
        plt.ylabel(var_label_y[V_i],fontsize=14)  

        plt.xlim(plot_ranges_x[V_i][0],plot_ranges_x[V_i][1])
        plt.ylim(plot_ranges_y[V_i][0],plot_ranges_y[V_i][1])
        plt.xticks(fontsize=10); plt.yticks(fontsize = 10) # Setting the X-ticks ## rotation='horizontal' # rotation='vertical'   
        if  np.abs( plot_ranges_x[V_i][1] ) < 0.001 :
            plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))  
        if  np.abs( plot_ranges_y[V_i][1] ) < 0.001 :
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))  
    
plt.subplots_adjust(left=0.05, bottom=0.1, right=0.85, top=0.92, hspace=0.35, wspace=0.43) # the amount of height/width reserved for space between subplots
#plt.legend(shadow=True, loc='right', bbox_to_anchor=(2.1, 2.5), markerscale=1 , prop={'size': 14})
if d_trend=='yes':
    plt.suptitle('CMIP6 models - North Pacific (circles) and North Atlantic (diamonds) - '+str(start_date_plt)+'-'+str(end_date_plt)+' - detrended' , fontsize=18)
else:
    plt.suptitle('CMIP6 models - North Pacific (circles) and North Atlantic (diamonds) - '+str(start_date_plt)+'-'+str(end_date_plt) , fontsize=18)

mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full

#if d_trend=='yes':
#    fig.savefig(dir_figs+'AllModels_CMIP6_Scatter_detrended_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
#else:
#    fig.savefig(dir_figs+'AllModels_CMIP6_Scatter_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')














