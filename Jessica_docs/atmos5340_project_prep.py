#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 11:57:29 2024

@author: u6044586
"""

import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
import xarray as xr 
import yaml 
import re
import os

met_dir='/uufs/chpc.utah.edu/common/home/haskins-group1/data/ExtData/GEOS_0.5x0.625_NA/MERRA2/2016/06/'

# Regex pattern to match dif Met files
met0_pth_str="MERRA2.201606\d{2}.A1.05x0625.NA.nc4"
met1_pth_str="MERRA2.201606\d{2}.I3.05x0625.NA.nc4"

# Get a list of met files:
met0_files = [f"{met_dir}/{x}" for x in os.listdir(met_dir) if re.match(met0_pth_str, x)]
met1_files = [f"{met_dir}/{x}" for x in os.listdir(met_dir) if re.match(met1_pth_str, x)]

file=met0_files[0]
met0= xr.open_mfdataset(file)
for var in met0.data_vars: 
    print(var.ljust(15), ': ', met0[var].attrs['long_name'], met0[var].coords.keys()) 
    
    
    
    
met0_keep=['CLDTOT','GRN','GWETTOP','LAI','PLBH','PRECTOT','T2M','U10M','V10M','T10M']
met1_keep=['T']







# infopath='/uufs/chpc.utah.edu/common/home/u6044586/ATMOS_5340/final_project/var_info_5340.csv'
# varz0=pd.read_csv(infopath)

# # vocs=[item for item in varz.iloc[np.where(varz['Is_VOC']==True)]['Species']]
# # nitrates=[item for item in varz.iloc[np.where(varz['Is_Nitrate']==True)]['Species']]
# # oxidants=[item for item in varz.iloc[np.where(varz['IS_Oxidant']==True)]['Species']]+['OH']
# # halogens=[item for item in varz.iloc[np.where(varz['Is_Halogen']==True)]['Species']]
# # cfcs=[item for item in varz.iloc[np.where(varz['Is_CFC']==True)]['Species']]
# # NOx_res=[item for item in varz.iloc[np.where(varz['Is_NOx_reservoir']==True)]['Species']]

# # Open the geoschem_config.yml file and load contents as dict: 
# spdb_yaml='/uufs/chpc.utah.edu/common/home/haskins-group1/users/amayhew/GEOS-Chem_Models/runs/H2O2_runs/gc_05x0625_NA_47L_merra2_fullchem/species_database.yml'      
# config_yaml='/uufs/chpc.utah.edu/common/home/haskins-group1/users/amayhew/GEOS-Chem_Models/runs/H2O2_runs/gc_05x0625_NA_47L_merra2_fullchem/geoschem_config.yml'      

# with open(spdb_yaml, 'r') as f:
#     spdb=yaml.load(f, Loader=yaml.FullLoader)
# with open(config_yaml, 'r') as f:
#     config=yaml.load(f, Loader=yaml.FullLoader)
    
# conc_path='/uufs/chpc.utah.edu/common/home/u6044586/ATMOS_5340/final_project/conc.nc4'
# conc = xr.open_mfdataset(conc_path)

# varz=varz0.copy()
# for sp in list(spdb.keys()):
#     if sp not in varz0['Species'] and sp in conc.data_vars: 
#         new_row=pd.DataFrame({'Species':[sp]})
#         for key in varz.columns: 
#             if key in list(spdb[sp].keys()):
#                 new_row[key]=spdb[sp][key]
#             elif key != 'Species': 
#                 new_row[key]=np.nan
                
#         varz = pd.concat([varz,new_row], ignore_index=True)
        
# varz.to_csv('/uufs/chpc.utah.edu/common/home/u6044586/ATMOS_5340/final_project/var_info.csv')
                
    
    