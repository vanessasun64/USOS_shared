#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 16:55:49 2025

@author: u6044586
"""
import os 
import xarray as xr
import pandas as pd
#from scipy.io import savemat
from collections import OrderedDict
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as ticker

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
# Load in data about all MCM variables: 
# mcm_info='/uufs/chpc.utah.edu/common/home/u6044586/python_scripts/modules/pyMCM/data/v331_AllRxn_scrape_all.xlsx'
# mcm=pd.read_excel(mcm_info, index_col=0)

# meta_info='/uufs/chpc.utah.edu/common/home/haskins-group1/data/Campaign_Data/Raw_Data/USOS_2024/RA/parked/merged/revised_USOS_parked_var_info_spreadsheet.xlsx'
# meta=pd.read_excel(meta_info, index_col=0)

# #  Figure out which variables in the MCM mechanism match the values in the data based on InChI matches. 
# merged_df = pd.merge(meta, mcm, on='InChI', how='inner')

# # Make a dict with MCM varname as key and data value varname as value. 
# grouped = merged_df.groupby('MCM_Name')['NEW_VARNAME'].apply(list)
# result_dict = grouped.to_dict()

# # Define column widths for alignment
# key_width = max(len(key) for key in result_dict.keys()) + 2  # Add space for quotes
# value_width = max(len(value) for values in result_dict.values() for value in values) + 10  # Add length for 'USOS.'


# # Open the file in write mode
# for_foam='/uufs/chpc.utah.edu/common/home/u6044586/mcm_usos_matches.txt'
# with open(for_foam, 'w') as file:
#     # Iterate over each key-value pair in the dictionary
#     for key, values in result_dict.items():
#         # Iterate over each value in the list associated with the key
#         for index, value in enumerate(values):
#             # Format the line as specified
#             line = f"'{key}'".ljust(key_width) + f"USOS.{value}".ljust(value_width) + "hold;\n"
#             # Write the formatted line to the file
#             file.write(line)

using=['O3_ppbv','CO_Piccaro','NO_LIF','NO2_LIF','PAN_CIMS','PPN_CIMS','iPropONO2_WAS',
       'nPropONO2_WAS','HNO3_CIMS','N2O5_CIMS','Isoprene_PTR','Alpha_Pinene_WAS',
       'Beta_Pinene_WAS','Limonene_WAS','CH4_Piccaro','Ethyne_WAS','Ethene_WAS',
       'Ethane_WAS','Propene_WAS','Propane_WAS','iButane_WAS','iPentane_WAS',
       'nButane_WAS','nPentane_WAS','nHexane_WAS','nHeptane_WAS','nOctane_WAS',
       'nNonane_WAS','Acrolein_PTR','Benzene_PTR','Toluene_PTR','Ethyl_WAS',
       'x135_TriMethylBenzene_WAS','o_Xylene_WAS','m_p_Xylene_WAS','Benzaldehyde_PTR', 
       'Styrene_PTR ','Acetone_WAS','Acetaldehyde_PTR','HCHO_CRDS','HCOOH_CIMS',
       'Methanol_PTR','Ethanol_PTR','MACR_WAS','Octanal_PTR','Nonanal_PTR', 
       'Acetone_Propanal_PTR','C9Aromatics_PTR','Monoterpenes_PTR','MVK_MACR_PTR','HONO_CIMS']
        
# Define the directory and file pattern (dirpath has no / after it!)
#for file in matching_files[0:1]: 
all_parked='/uufs/chpc.utah.edu/common/home/haskins-group1/data/Campaign_Data/Raw_Data/USOS_2024/RA/parked/merged/rev_30min/all_parked_rev30min.nc'
# Load in the USOS 30min parked data on 1 day: 
ds=xr.open_dataset(all_parked)


# Convert the Dataset to a DataFrame
df = ds.to_dataframe()

no_have=[var for var in using if var not in df.columns] 
for var in no_have: 
    using.remove(var)

# Setting 'time' as the index for more effective interpolation
df.set_index('time_local', inplace=True)

df_interp=df.copy()
for i,col in enumerate( using):
    if i<10: 
        # Non-linear interpolation using cubic method
        df_interp[col] = df[col].interpolate(method='linear')
        fig, ax = plt.subplots()
        plt.plot(df.index, df[col], color='k', marker='o',label='Original')
        plt.plot(df_interp.index, df_interp[col], color='r', marker='x', label='Interpolated')
        plt.legend()
    
        ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
        ax.xaxis.set_major_locator(mdates.HourLocator(interval=3))
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        ax.set_title(col)
        plt.show()
#we might want to save the interpolation output to show how close it looks



# # Convert the dataframe to a nested dictionary (so scipy can output to a matalb structure!) 
# ddict=dataframe_to_nested_dict(df)

# # Sort alphabetically so not annoying in MATLAB...  
# ddict= OrderedDict(sorted(ddict.items())) 

# # Save the USOS data in an output .mat file: 
# outpath = '/uufs/chpc.utah.edu/common/home/u1545774/Documents/F0AM/F0AM-4.3.0.1/Setups/'
# matfilename = os.path.basename(file).replace('.nc','.mat')
# savemat(outpath+matfilename,{'dat': ddict})


