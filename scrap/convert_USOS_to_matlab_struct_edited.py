#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 16:55:49 2025

@author: u6044586
"""
import os 
import xarray as xr
import pandas as pd
from scipy.io import savemat
from collections import OrderedDict


def find_files(directory, search_string, ext):
    matches = []
    for root, dirs, files in os.walk(directory):
        for filename in files:
            if search_string in filename and filename.endswith(ext):
                matches.append(os.path.join(root, filename))
    return matches

def dataframe_to_nested_dict(df):
    nested_dict = {}
    for column in df.columns:
        nested_dict[column] = df[column][0:-1].to_numpy()
        
    # Add the index as a key-value pair
    nested_dict['time_UTC'] = df.index[0:-1].to_numpy()
    return nested_dict

# %% BEGIN MAIN 

# Define the directory and file pattern (dirpath has no / after it!)
dirpath = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_merges/R0/CSL_MobileLab_Parked/rev_30min'
file_pattern = 'USOS_rev30min_Merged_Parked'

# Call the function and store the results
files = find_files(dirpath, file_pattern, '.nc')

#for file in matching_files[0:1]: 
for file_num in range(0,len(files)):
    file=files[file_num]
    print('loading file: ' + file)

    # Load in the USOS 30min parked data on 1 day: 
    ds=xr.open_dataset(file)

    # Convert the Dataset to a DataFrame
    df = ds.to_dataframe()

    # Convert everything to F0AM units: 
    #df['CO2_ppb']= df['CO2_ppm_1min']*1e3

    # Convert units of all vars to F0AM expected units & return dict with units of each var.
    #df1,units_dict=convert_to_F0AM_units(df0, show=False)   

    # Convert the dataframe to a nested dictionary (so scipy can output to a matalb structure!) 
    ddict=dataframe_to_nested_dict(df)

    # Sort alphabetically so not annoying in MATLAB...  
    ddict= OrderedDict(sorted(ddict.items())) 

    # Save the USOS data in an output .mat file: 
    outpath = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_merges/matlab_merges/30min_Parked/'
    matfilename = os.path.basename(file).replace('.nc','.mat')
    savemat(outpath+matfilename,{'USOS': ddict})
    print('Finished merging: ' + outpath+matfilename)


