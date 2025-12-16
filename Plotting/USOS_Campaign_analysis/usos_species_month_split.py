import numpy as np 
import os 
import xarray as xr
import pandas as pd

#from scipy.io import savemat
from collections import OrderedDict
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
import matplotlib as mpl
from matplotlib import colors
from matplotlib import colormaps
import matplotlib.colors as mcolors
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import cm
import matplotlib

import mplcursors

#load file for all dates
all_days_filepath = '../../CampaignData_and_Merges/R0/CSL_MobileLab_Parked/merged/rev_30min/all_CSL_MobileLab_Parked_rev30minv4.nc'
all_days_filepath_load = xr.open_dataset(all_days_filepath)
df_alldays = all_days_filepath_load.to_dataframe()
df_alldays.reset_index(inplace=True)
df_alldays.set_index('time_local', inplace=True, drop=False)

#adds padding of NaNs for first day and last day of campaign, since they need to have the same length of time as the other days in order to plot
new_start_time = pd.Timestamp('2024-07-14 00:00:00')
new_end_time = pd.Timestamp('2024-08-18 23:30:00')

# Create a new datetime index from new_start to the end of existing index with same frequency
new_index = pd.date_range(start=new_start_time, end=new_end_time, freq='30T')

# Reindex the dataframe to include new rows
df_alldays = df_alldays.reindex(new_index)

fig, ax = plt.subplots(1,2, figsize = (15,10))
hour_range = np.arange(0,24,0.5)

cmap1 = mpl.colormaps['YlOrRd']
cmap2 = mpl.colormaps['YlGnBu']

# Create 20 color values for gradient in the colormap
color_values = np.linspace(0.2, 1, 20) 

# july_ozone_list = []
# aug_ozone_list = []
overall_ozone_list = []

#Loop through each day for July and plot the ozone values on first subplot
for day_val in range(0,18):
    days_indexing = df_alldays.index.date[day_val*48]
    one_day_df = df_alldays.sort_index().loc[str(days_indexing)]
    species_vals = one_day_df['O3_ppbv'].values
    #july_ozone_list.append(species_vals)
    overall_ozone_list.append(species_vals)
    toplot, = ax[0].plot(hour_range, species_vals, color=cmap1(color_values[day_val]), label = days_indexing)
    ax[0].set_xlabel('Hour')
    ax[0].set_ylabel('Ozone Concentration (ppb)')
    ax[0].set_xticks(np.arange(0,24))
    ax[0].set_ylim([0,100])
    ax[0].set_yticks(np.arange(0,105,5))
    ax[0].set_title('July Ozone')

# avg_july_ozone = np.nanmean(july_ozone_list, axis = 0)
# ax[0].plot(hour_range,avg_july_ozone, color='k', label = 'Average July 14-31')

# std_avg_july_ozone = np.std(avg_july_ozone, axis=0)
# ax[0].fill_between(hour_range, avg_july_ozone - std_avg_july_ozone, avg_july_ozone + std_avg_july_ozone, alpha = 0.2)

#Loop through each day for August and plot the ozone values on second subplot
for day_val in range(0,18):
    days_indexing = df_alldays.index.date[(day_val*48)+864]
    one_day_df = df_alldays.sort_index().loc[str(days_indexing)]
    species_vals = one_day_df['O3_ppbv'].values
    #aug_ozone_list.append(species_vals)
    overall_ozone_list.append(species_vals)
    ax[1].plot(hour_range, species_vals, color=cmap2(color_values[day_val]), label = days_indexing)
    ax[1].set_xlabel('Hour')
    ax[1].set_ylabel('Ozone Concentration (ppb)')
    ax[1].set_xticks(np.arange(0,24))
    ax[1].set_ylim([0,100])
    ax[1].set_yticks(np.arange(0,105,5))
    ax[1].set_title('August Ozone')

# avg_aug_ozone = np.nanmean(aug_ozone_list, axis = 0)
# ax[1].plot(hour_range,avg_aug_ozone, color='k', label = 'Average Aug. 14-31')

# std_avg_aug_ozone = np.std(avg_aug_ozone, axis=0)
# ax[1].fill_between(hour_range, avg_aug_ozone - std_avg_aug_ozone, avg_aug_ozone + std_avg_aug_ozone, alpha = 0.2)

avg_overall_ozone = np.nanmean(overall_ozone_list, axis = 0)
std_avg_overall_ozone = np.std(avg_overall_ozone, axis=0)
ax[0].plot(hour_range,avg_overall_ozone, color='k', label = 'Average for campaign', linewidth = 3)
ax[0].fill_between(hour_range, avg_overall_ozone - std_avg_overall_ozone, avg_overall_ozone + std_avg_overall_ozone, alpha = 0.2, color = 'tab:gray')
ax[1].plot(hour_range,avg_overall_ozone, color='k', label = 'Average for campaign', linewidth = 3)
ax[1].fill_between(hour_range, avg_overall_ozone - std_avg_overall_ozone, avg_overall_ozone + std_avg_overall_ozone, alpha = 0.2, color = 'tab:gray')

#get labels for every 3 days to show on legend
july_handles, july_labels = ax[0].get_legend_handles_labels()

skip=3
july_selected_handles = july_handles[::skip]
july_selected_labels = july_labels[::skip]
ax[0].legend(july_selected_handles, july_selected_labels)
ax[0].margins(x=0)
ax[0].grid()

#get labels for every 3 days to show on legend
aug_handles, aug_labels = ax[1].get_legend_handles_labels()

skip=3
aug_selected_handles = aug_handles[::skip]
aug_selected_labels = aug_labels[::skip]
ax[1].legend(aug_selected_handles, aug_selected_labels)
ax[1].margins(x=0)
ax[1].grid()

ax[0].hlines(y=70, xmin = 0, xmax = 23.5, color = 'tab:red')
ax[1].hlines(y=70, xmin = 0, xmax = 23.5, color = 'tab:red')

plt.show()