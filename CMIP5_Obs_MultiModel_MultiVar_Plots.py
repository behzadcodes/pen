

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



###########################################################
###   Observational MultiVariable Plotting - Subplots  ####
###########################################################
        
from pathlib import Path    
dir_data_in_o=('/data5/scratch/Behzad/ocean/ocean_biogeochemistry/processed_data/python_data_Obs_bio/')
GCM_Names = ['Obs']

Time_months = ('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec');
P_Var_x = np.linspace(1,12,12);

plot_ranges_1=  [ [-0.01,0.2],                       [-0.01,20],                      [-0.01,20],                        [-0.01,0.5],                              [-0.05,1.4],                              [-1,25],                              [-0.1,1.5],                 [-1,18],                    [0,50]]
var_names_1=    [ 'Biomass_BRF_allmonths',           'Biomass_STK_allmonths',         'Biomass_THK_allmonths',           'PP_by_Phyto_allmonths',                  'PPint_allmonths',                       'Si_allmonths',                                 'Fe_allmonths',             'NO3_allmonths',            'PAR_allmonths']
plot_titles_1=  [ 'Biomass - Behrenfeld (m mol/m3)', 'Biomass - Stramski (m mol/m3)', 'Biomass - Kostadinov (m mol/m3)', 'Specific Growth rate (PP/phyto * 1E-3)', 'Primary Production (micromol.m-2.s-1)', 'Silicate - Si (m mol/m3)',                     'Iron - Fe (micro mol/m3)', 'Nitrate - NO3 (m mol/m3)', 'Light - (PAR)']  
var_multiplier_1=[ 1000,                              1000,                             1000,                             1000,                                     1E6,                                     1000,                                          1E6,                        1000,                        1         ]
colors_vars_1 = [ 'red',                             'navy',                           'green',                          'grey',                                   'lime',                                  'deepskyblue',                                  'gold',                     'orange',                     'magenta'                            ]

plot_ranges_2=  [ [-1,10],                  [-1,110],                      [-1,5],                   [4,20],          [-110,215],                             [0,1000],                 [-20,40],                                        [270,460],          ]
var_names_2=    [ 'Large Phyto_allmonths',  '% Large Phyto_allmonths',    'Small Phyto_allmonths',  'SST_allmonths', 'HFDS_allmonths',                       'MLD_allmonths',         'FgCO2_allmonths',                                'SpCO2_allmonths', ]
plot_titles_2=  [ 'Large Phyto (m mol/m3)', 'Large Phyto percentage (%)', 'Small Phyto (m mol/m3)', 'SST (C)',       'Downward Heat Flux at Surface (W/m2)', 'Mixed Layer Depth (m)', 'Air-sea CO2 flux (donward +) (gram C /m2/year)', 'SpCO2 (ppm)',    ]  
var_multiplier_2=[1000,                     1,                            1000,                      1,               1,                                      1,                       1,                                                 1,               ]
colors_vars_2 = ['navy',                   'blue',                        'deepskyblue',            'orange',        'crimson',                              'darkviolet',            'lime',                                            'green']


d_trend='yes'

for M_i in range(len(GCM_Names)): # M_i=0     # GCM=GCM_Names[0]

    #M_i=0

    fig=plt.figure() #; fig, ax = plt.subplots()     
    for AX_i in range(4):    
        ax = fig.add_subplot(2,2,AX_i+1)  
        axes = [ax, ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx()] 
        
        if AX_i+1==1 or AX_i+1==2:  # Axes 1 and 2 (first row) are for the first set of variables
            var_names=var_names_1 ; plot_titles=plot_titles_1 ; var_multiplier=var_multiplier_1 ; colors_vars=colors_vars_1 ; plot_ranges=plot_ranges_1
        elif AX_i+1==3 or AX_i+1==4:   # Axes 3 and 4 (second row) are for the second set of variables
            var_names=var_names_2 ; plot_titles=plot_titles_2 ; var_multiplier=var_multiplier_2 ; colors_vars=colors_vars_2 ; plot_ranges=plot_ranges_2
        
        V_i=-1
        for ax in axes:            
            V_i += 1
            if V_i < len(var_names):
                if AX_i+1==1 or AX_i+1==3:
                    if V_i>0:
                        axes[V_i].spines['right'].set_position(('axes', -0.13-V_i*0.13)); axes[V_i].set_frame_on(True); axes[V_i].patch.set_visible(False)
                else:
                    axes[V_i].spines['right'].set_position(('axes', -10-V_i*0.14)); axes[V_i].set_frame_on(True); axes[V_i].patch.set_visible(False)
                
                if var_names[V_i]=='Zoo_by_Phyto_allmonths':
                    filename_in = (dir_data_in_o + var_names[V_i]+'_Obs.npy') # Directory to save processed data # Directory to save processed data
                    my_file = Path(filename_in)         
            

                else:                   

                    filename_in = (dir_data_in_o + var_names[V_i]+'_Obs.npy') # Directory to save processed data # Directory to save processed data
                    my_file = Path(filename_in)
                    if my_file.is_file(): # If that variable exists in the saved data    
                        Variable_1 = np.load(filename_in)
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
        plt.suptitle( 'Observations - Variable seasonal cylce - '+str(start_date_plt)+'-'+str(end_date_plt)+ ' - detrended' , fontsize=18)
    else:
        plt.suptitle( 'Observations - Variable seasonal cylce - '+str(start_date_plt)+'-'+str(end_date_plt) , fontsize=18)


    if d_trend=='yes':
        fig.savefig(dir_figs+'Observations_SesonalCycle_AllVar_detrended_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
    else:
        fig.savefig(dir_figs+'Observations__SesonalCycle_AllVar_1.png', format='png', dpi=300, transparent=True, bbox_inches='tight')











for M_i in range(len(GCM_Names)): # M_i=5     # GCM=GCM_Names[0]
    
    GCM=GCM_Names[M_i]
    print (GCM)  
       

    fig=plt.figure() #; fig, ax = plt.subplots() 
    ax = fig.add_subplot(2,2,1)  
    axes = [ax, ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx()]       
    
    V_i=-1
    for ax in axes:
        
        V_i += 1
        if V_i < len(var_names):
            if V_i>0:
                axes[V_i].spines['right'].set_position(('axes', -0.1-V_i*0.14)); axes[V_i].set_frame_on(True); axes[V_i].patch.set_visible(False)
                #axes[1].spines['right'].set_position(('axes', -0.2)); axes[1].set_frame_on(True); axes[1].patch.set_visible(False)
            filename_in = (dir_data_in + 'AllResults_'+GCM+'_'+file_names[V_i]+'_1991_2010.out') # Directory to save processed data
            my_shelf = shelve.open(filename_in)
            if var_names[V_i] in my_shelf: # If that variable exists in the saved data    
                globals()['Variable_1']=my_shelf[var_names[V_i]]
                Variable_1=Variable_1 * var_multiplier[V_i]
        
                Var_1_allmonths=copy.deepcopy(Variable_1); 
                Var_1_allmonths_PAPA=Var_1_allmonths[:,135:140,210:220]
        
                if d_trend=='yes': # If it is requested to detrend the data before calculations
                    Var_1_allmonths_PAPA= func_detrend_3d(Var_1_allmonths_PAPA)
          
                Var_1_monthlymean_PAPA = Var_1_allmonths_PAPA.reshape(( np.int(Var_1_allmonths_PAPA.shape[0]/12) ,12,Var_1_allmonths_PAPA.shape[1],Var_1_allmonths_PAPA.shape[2])) # Monthly means of SST over the time period
                Var_1_monthlymean_PAPA = np.nanmean(Var_1_monthlymean_PAPA,axis=0)
                Var_1_monthlymean_PAPA=np.squeeze(Var_1_monthlymean_PAPA)           
                Var_1_monthlymean_PAPA_timeseries=np.nanmean( np.nanmean(Var_1_monthlymean_PAPA,axis=2) , axis=1); Var_1_monthlymean_PAPA_timeseries=np.squeeze(Var_1_monthlymean_PAPA_timeseries)
           
                ax.plot(P_Var_x, Var_1_monthlymean_PAPA_timeseries, c=colors_vars[V_i], label=plot_titles[V_i], marker='o', markersize=7, markerfacecolor='W', linewidth=2.0)    
            ax.set_ylim(plot_ranges[V_i])   
            ax.set_ylabel(plot_titles[V_i], color=colors_vars[V_i])#, labelpad = -5-V_i*5)
            ax.tick_params(axis='y', colors=colors_vars[V_i])
    plt.xticks(P_Var_x, Time_months, rotation=45, fontsize=10)
    #axes[0].set_xticks(Time_months)
    ax.set_title('North Pacific [45N-50N, 140W-150W]')

    ax = fig.add_subplot(2,2,2)   
    axes = [ax, ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx()]   

    V_i=-1
    for ax in axes:
        
        V_i += 1
        if V_i < len(var_names):
            if V_i>0:
                axes[V_i].spines['right'].set_position(('axes', -5-V_i*0.1)); axes[V_i].set_frame_on(True); axes[V_i].patch.set_visible(False)
            filename_in = (dir_data_in + 'AllResults_'+GCM+'_'+file_names[V_i]+'_1991_2010.out') # Directory to save processed data
            my_shelf = shelve.open(filename_in)
            if var_names[V_i] in my_shelf: # If that variable exists in the saved data    
                globals()['Variable_1']=my_shelf[var_names[V_i]]
                Variable_1=Variable_1 * var_multiplier[V_i]
        
                Var_1_allmonths=copy.deepcopy(Variable_1); 
                Var_1_allmonths_NABE=Var_1_allmonths[:,135:140,325:335]
        
                if d_trend=='yes': # If it is requested to detrend the data before calculations
                    Var_1_allmonths_NABE= func_detrend_3d(Var_1_allmonths_NABE)
          
                Var_1_monthlymean_NABE = Var_1_allmonths_NABE.reshape(( np.int(Var_1_allmonths_NABE.shape[0]/12) ,12,Var_1_allmonths_NABE.shape[1],Var_1_allmonths_NABE.shape[2])) # Monthly means of SST over the time period
                Var_1_monthlymean_NABE = np.nanmean(Var_1_monthlymean_NABE,axis=0)
                Var_1_monthlymean_NABE=np.squeeze(Var_1_monthlymean_NABE)           
                Var_1_monthlymean_NABE_timeseries=np.nanmean( np.nanmean(Var_1_monthlymean_NABE,axis=2) , axis=1); Var_1_monthlymean_NABE_timeseries=np.squeeze(Var_1_monthlymean_NABE_timeseries)
           
                ax.plot(P_Var_x, Var_1_monthlymean_NABE_timeseries, c=colors_vars[V_i], label=plot_titles[V_i], marker='o', markersize=7, markerfacecolor='W', linewidth=2.0)    
            ax.set_ylim(plot_ranges[V_i])      
            ax.set_ylabel(plot_titles[V_i], color=colors_vars[V_i])#, labelpad = -5-V_i*5)
            ax.tick_params(axis='y', colors=colors_vars[V_i])
    plt.xticks(P_Var_x, Time_months, rotation=45, fontsize=10)
    axes[0].set_ylabel(' ')
    #axes[0].set_xticks(Time_months)
    ax.set_title('North Atlantic [45N-50N, 25W-35W]')
    
    plt.subplots_adjust(left=0.36, bottom=0.05, right=0.96, top=0.9, hspace=0.2, wspace=0.1) # the amount of height/width reserved for space between subplots
    if d_trend=='yes':
        plt.suptitle( GCM+' - Variable seasonal cylce - '+str(start_date_plt)+'-'+str(end_date_plt)+ ' - detrended' , fontsize=18)
    else:
        plt.suptitle( GCM+' - Variable seasonal cylce - '+str(start_date_plt)+'-'+str(end_date_plt) , fontsize=18)






















import matplotlib.pyplot as plt
import numpy as np


fig=plt.figure()
ax = fig.add_subplot(n_r,n_c,1)  
#fig, ax = plt.subplots()

axes = [ax, ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx(), ax.twinx()]

fig.subplots_adjust(left=0.5)

#axes[1].spines['right'].set_position(('axes', 1.1)); axes[1].set_frame_on(True); axes[1].patch.set_visible(False)
#axes[2].spines['right'].set_position(('axes', 1.2)); axes[2].set_frame_on(True); axes[2].patch.set_visible(False)

axes[1].spines['right'].set_position(('axes', -0.2)); axes[1].set_frame_on(True); axes[1].patch.set_visible(False)
axes[2].spines['right'].set_position(('axes', -0.3)); axes[2].set_frame_on(True); axes[2].patch.set_visible(False)
axes[3].spines['right'].set_position(('axes', -0.4)); axes[3].set_frame_on(True); axes[3].patch.set_visible(False)
axes[4].spines['right'].set_position(('axes', -0.5)); axes[4].set_frame_on(True); axes[4].patch.set_visible(False)
axes[5].spines['right'].set_position(('axes', -0.6)); axes[5].set_frame_on(True); axes[5].patch.set_visible(False)
axes[6].spines['right'].set_position(('axes', -0.7)); axes[6].set_frame_on(True); axes[6].patch.set_visible(False)
axes[7].spines['right'].set_position(('axes', -0.8)); axes[7].set_frame_on(True); axes[7].patch.set_visible(False)
axes[8].spines['right'].set_position(('axes', -0.9)); axes[8].set_frame_on(True); axes[8].patch.set_visible(False)

V_i=-1
for ax in axes:
    
    V_i += 1
    if V_i < len(var_names):
        filename_in = (dir_data_in + 'AllResults_'+GCM+'_'+file_names[V_i]+'_1991_2010.out') # Directory to save processed data
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
       
            ax.plot(P_Var_x, Var_1_monthlymean_PAPA_timeseries, c=colors_vars[V_i], label=plot_titles[V_i], marker='o', markersize=7, markerfacecolor='W', linewidth=2.0)    
        
        #ax.xticks(P_Var_x, '', rotation=45) # Setting the X-ticks ## rotation='horizontal' # rotation='vertical'
        #ax.xticks(P_Var_x, Time_months, rotation=45, fontsize=10) # Setting the X-ticks ## rotation='horizontal' # rotation='vertical'      
        ax.set_ylabel(plot_titles[V_i], color=colors_vars[V_i])#, labelpad = -5-V_i*5)
        ax.tick_params(axis='y', colors=colors_vars[V_i])
axes[0].set_xlabel('Months')
#axes[0].set_xticks(Time_months)

plt.show()












