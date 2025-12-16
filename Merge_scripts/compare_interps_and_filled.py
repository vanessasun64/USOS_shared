import os 
import sys 
import re 
import yaml
import inspect 
import numpy as np 
import pandas as pd
import xarray as xr
from collections import defaultdict
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.colors import ListedColormap

from scipy.io import savemat
from collections import OrderedDict

#Plot formatting
mpl.rcParams['xtick.labelsize'] = 15
mpl.rcParams['ytick.labelsize'] = 15
mpl.rcParams['legend.fontsize'] = 16
mpl.rcParams['axes.labelsize'] = 18
mpl.rcParams['axes.titlesize'] = 28
mpl.rcParams['axes.xmargin'] = 0

# Set the font family to 'serif'
mpl.rcParams['font.family'] = 'serif'
# Specify preferred serif font (Computer Modern Roman is 'cmr10')
mpl.rcParams['font.serif'] = 'Lato' 
# Optionally, configure mathtext to use Computer Modern fonts as well
mpl.rcParams['mathtext.fontset'] = 'cm'
# Ensure minus signs are rendered correctly with CM fonts
mpl.rcParams['axes.unicode_minus'] = False

# Main datapaths
dirpath = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/'
#define path for Hawthorne data directory
savepath= dirpath + 'Merge_scripts/'
plots_savepath= dirpath + 'Merge_scripts/plots/'

ml_nc_filepath = dirpath + 'CampaignData_and_Merges/R0/CSL_MobileLab_Parked/merged/rev_15min/all_CSL_MobileLab_Parked_rev1hr_iWASupdated.nc'
ml_data = xr.open_dataset(ml_nc_filepath, engine='netcdf4')
print(ml_data.data_vars)
# df_ml_data = ml_data.to_dataframe()
ml_data.close()
print('Finished!')

# df_ml_data.reset_index(inplace=True)
# df_ml_data.set_index('time_local', inplace=True)
# print('hello world')

# #We want to see if we can reasonably interpolate values for the missing hours. Plot the time series comparison.
# df_interp = df_iwas.copy()

# for species in df_iwas.columns[3:60]:
#     n_baddies= len([item for item in df_iwas[species] if item <0 or np.isnan(item)])
#     if n_baddies > 0:
#         df_interp[species] = df_interp[species].interpolate(method='linear')
#         fig, ax = plt.subplots(figsize = (20,4), constrained_layout=True)
#         plt.plot(df_iwas.index, df_iwas[species], color = 'k', marker = 'o', label=f'Original (baddies={n_baddies})')
#         plt.plot(df_iwas.index, df_interp[species], color = 'r', marker = 'x', label = 'Interpolated')
#         plt.xlabel('Date/Time (UTC)')
#         plt.margins(x=0)
#         midnight_vals = []
#         for midnight_idx in range(6,len(df_iwas.index),24):
#             midnight_vals.append(df_iwas.index[midnight_idx])
#         for day_pos in midnight_vals:
#             ax.axvline(day_pos, color = 'black', linestyle = 'dotted', alpha = 0.7)
#         plt.title(species)
#         plt.legend()
