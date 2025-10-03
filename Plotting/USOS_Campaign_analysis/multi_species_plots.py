import numpy as np 
import os 
import xarray as xr
import pandas as pd

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

import os 
import sys 
import re 
import yaml
import inspect 
import numpy as np 
import pandas as pd
import xarray as xr
from collections import defaultdict
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt

from sklearn.metrics import r2_score

#get path for USOS_shared directory for correct laptop
#input parameter should be 'CHPC', 'Mac', or 'Windows'
# current_dir = os.path.dirname(os.path.abspath(__file__))
# global_scripts_path = os.path.abspath(os.path.join(current_dir, "..", "global_scripts"))
# sys.path.insert(0,global_scripts_path)
# from dirpath import filepath_source
# dirpath = filepath_source('CHPC')
dirpath = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/'

#load file for all dates
#all_days_filepath = '../../CampaignData_and_Merges/R0/CSL_MobileLab_Parked/merged/rev_30min/all_CSL_MobileLab_Parked_rev30minv4.nc'
# all_days_filepath = '../../CampaignData_and_Merges/R0/CSL_MobileLab_Parked/merged/rev_30min/all_CSL_MobileLab_Parked_rev30minv4.nc'
# all_days_filepath_load = xr.open_dataset(all_days_filepath)
# df_alldays = all_days_filepath_load.to_dataframe()
# df_alldays.reset_index(inplace=True)
# df_alldays.set_index('time_local', inplace=True, drop=False)

# #adds padding of NaNs for first day and last day of campaign, since they need to have the same length of time as the other days in order to plot
# new_start_time = pd.Timestamp('2024-07-14 00:00:00')
# new_end_time = pd.Timestamp('2024-08-18 23:30:00')

# # Create a new datetime index from new_start to the end of existing index with same frequency
# new_index = pd.date_range(start=new_start_time, end=new_end_time, freq='30T')

# # Reindex the dataframe to include new rows
# df_alldays = df_alldays.reindex(new_index)
# # Now df_expanded has rows starting at 00:00 with NaNs for the new rows

# #add NO and NO2 to get NOx values
# df_alldays['NOx_ppbv'] = df_alldays['NO_LIF'] + df_alldays['NO2_LIF']

# month_index = df_alldays.index.month
# date_index = df_alldays.index.date
# day_names_index = df_alldays.index.day_name()
# hr_of_day_index = df_alldays.index.time

days_of_week = ['Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Saturday', 'Sunday']

# hour_range = np.arange(0,24,0.5)
hour_range = np.arange(0,24,1)

def species_campaign_total_monthly_split_plot(species_var, species_name, y_limit, y_ticks, color_map1, color_map2, plot_name):
    df_species = df_alldays[species_var]
    df_july_species = df_species.iloc[df_species.index.month == 7]
    july_dates = pd.date_range(start='2024-07-14', end='2024-07-31')
    july_dates_list = july_dates.strftime('%Y-%m-%d')

    df_august_species = df_species.iloc[df_species.index.month == 8]
    august_dates = pd.date_range(start='2024-08-01', end='2024-08-18')
    august_dates_list = august_dates.strftime('%Y-%m-%d')

    overall_species_avg = df_species.groupby([hr_of_day_index]).mean()
    overall_species_std = np.std(overall_species_avg)
    overall_dates = pd.date_range(start='2024-07-14', end='2024-08-18')
    overall_dates_list = overall_dates.strftime('%Y-%m-%d')

    fig, ax = plt.subplots(1,2, figsize = (20,10))

    cmap1 = mpl.colormaps[color_map1]
    cmap2 = mpl.colormaps[color_map2]

    # Create 20 color values for gradient in the colormap
    color_values = np.linspace(0.15, 0.95, 20) 

    #Loop through each day for July and plot the ozone values on first subplot
    for day_val in range(0,18):
        toplot, = ax[0].plot(hour_range, df_july_species[july_dates_list[day_val]], color=cmap1(color_values[day_val]), label = july_dates_list[day_val])
        ax[0].set_xlabel('Hour')
        ax[0].set_ylabel(species_name + ' Concentration (ppb)')
        ax[0].set_xticks(np.arange(0,24))
        ax[0].set_ylim(y_limit)
        ax[0].set_yticks(y_ticks)
        ax[0].set_title('July ' + species_name)
    ax[0].plot(hour_range, overall_species_avg, color = 'k', label = 'Campaign Average')
    ax[0].fill_between(hour_range, overall_species_avg - overall_species_std, overall_species_avg + overall_species_std, color='k',alpha = 0.1)

    #get labels for legend
    handles, labels = ax[0].get_legend_handles_labels()
    skip=3
    selected_handles = handles[::skip]
    selected_labels = labels[::skip]
    ax[0].legend(selected_handles, selected_labels)
    ax[0].margins(x=0)
    ax[0].grid()

    #Loop through each day for August and plot the ozone values on second subplot
    for day_val in range(0,18):
        toplot, = ax[1].plot(hour_range, df_august_species[august_dates_list[day_val]], color=cmap2(color_values[day_val]), label = august_dates_list[day_val])
        ax[1].set_xlabel('Hour')
        ax[1].set_ylabel(species_name + ' Concentration (ppb)')
        ax[1].set_xticks(np.arange(0,24))
        ax[1].set_ylim(y_limit)
        ax[1].set_yticks(y_ticks)
        ax[1].set_title('August ' + species_name)
    ax[1].plot(hour_range, overall_species_avg, color = 'k', label = 'Campaign Average')
    ax[1].fill_between(hour_range, overall_species_avg - overall_species_std, overall_species_avg + overall_species_std, color='k', alpha = 0.1)

    #get labels for legend
    handles, labels = ax[1].get_legend_handles_labels()
    skip=3
    selected_handles = handles[::skip]
    selected_labels = labels[::skip]
    ax[1].legend(selected_handles, selected_labels)
    ax[1].margins(x=0)
    ax[1].grid()

    plt.tight_layout()
    #plt.savefig('../haskins-group1/users/vsun/USOS_shared//Plotting/USOS_Campaign_analysis/plots/' + plot_name)
    plt.show()

def ozone_campaign_total_monthly_split_plot(species_var, species_name, y_limit, y_ticks, color_map1, color_map2, plot_name):
    df_ozone = df_alldays[species_var]
    df_july_ozone = df_ozone.iloc[df_ozone.index.month == 7]
    july_dates = pd.date_range(start='2024-07-14', end='2024-07-31')
    july_dates_list = july_dates.strftime('%Y-%m-%d')

    df_august_ozone = df_ozone.iloc[df_ozone.index.month == 8]
    august_dates = pd.date_range(start='2024-08-01', end='2024-08-18')
    august_dates_list = august_dates.strftime('%Y-%m-%d')

    overall_ozone_avg = df_ozone.groupby([hr_of_day_index]).mean()
    overall_ozone_std = np.std(overall_ozone_avg)
    overall_dates = pd.date_range(start='2024-07-14', end='2024-08-18')
    overall_dates_list = overall_dates.strftime('%Y-%m-%d')

    fig, ax = plt.subplots(1,2, figsize = (20,10))

    cmap1 = mpl.colormaps[color_map1]
    cmap2 = mpl.colormaps[color_map2]

    # Create 20 color values for gradient in the colormap
    color_values = np.linspace(0.15, 0.95, 20) 

    #Loop through each day for July and plot the ozone values on first subplot
    for day_val in range(0,18):
        toplot, = ax[0].plot(hour_range, df_july_ozone[july_dates_list[day_val]], color=cmap1(color_values[day_val]), label = july_dates_list[day_val])
        ax[0].set_xlabel('Hour')
        ax[0].set_ylabel(species_name + ' Concentration (ppb)')
        ax[0].set_xticks(np.arange(0,24))
        ax[0].set_ylim(y_limit)
        ax[0].set_yticks(y_ticks)
        ax[0].set_title('July ' + species_name)
    ax[0].plot(hour_range, overall_ozone_avg, color = 'k', label = 'Campaign Average')
    ax[0].fill_between(hour_range, overall_ozone_avg - overall_ozone_std, overall_ozone_avg + overall_ozone_std, color='k',alpha = 0.1)
    ax[0].hlines(y=70, xmin = 0, xmax = 23.5, color = 'tab:red')

    #get labels for legend
    handles, labels = ax[0].get_legend_handles_labels()
    skip=3
    selected_handles = handles[::skip]
    selected_labels = labels[::skip]
    ax[0].legend(selected_handles, selected_labels)
    ax[0].margins(x=0)
    ax[0].grid()

    #Loop through each day for August and plot the ozone values on second subplot
    for day_val in range(0,18):
        toplot, = ax[1].plot(hour_range, df_august_ozone[august_dates_list[day_val]], color=cmap2(color_values[day_val]), label = august_dates_list[day_val])
        ax[1].set_xlabel('Hour')
        ax[1].set_ylabel(species_name + ' Concentration (ppb)')
        ax[1].set_xticks(np.arange(0,24))
        ax[1].set_ylim(y_limit)
        ax[1].set_yticks(y_ticks)
        ax[1].set_title('August ' + species_name)
    ax[1].plot(hour_range, overall_ozone_avg, color = 'k', label = 'Campaign Average')
    ax[1].fill_between(hour_range, overall_ozone_avg - overall_ozone_std, overall_ozone_avg + overall_ozone_std, color='k', alpha = 0.1)
    ax[1].hlines(y=70, xmin = 0, xmax = 23.5, color = 'tab:red')

    #get labels for legend
    handles, labels = ax[1].get_legend_handles_labels()
    skip=3
    selected_handles = handles[::skip]
    selected_labels = labels[::skip]
    ax[1].legend(selected_handles, selected_labels)
    ax[1].margins(x=0)
    ax[1].grid()

    plt.tight_layout()
    #plt.savefig('../haskins-group1/users/vsun/USOS_shared/Plotting/USOS_Campaign_analysis/plots/' + plot_name)
    plt.show()

def species_by_day_of_week_plot(species_var, species_name, y_limit, y_ticks, color_map, plot_name):
    day_of_week_species_mean = df_alldays.groupby([day_names_index, hr_of_day_index])[species_var].mean()

    fig, ax = plt.subplots(figsize = (15,10))
    cmap = mpl.colormaps[color_map] #set colormap
    color_values = np.linspace(0.2, 1, 7) 

    #Loop through each day for July and plot the ozone values on first subplot
    for day_val in range(0,7):
        plt.plot(hour_range, day_of_week_species_mean[days_of_week[day_val]], color=cmap(color_values[day_val]))
    plt.legend(days_of_week)
    plt.xlabel('Hour')
    plt.ylabel(species_name + ' Concentration (ppb)')
    plt.xticks(np.arange(0,24))
    plt.ylim(y_limit)
    plt.yticks(y_ticks)
    plt.margins(x=0)
    plt.title('Average ' + species_name + ' Concentration by Day of Week')

    plt.grid()
    plt.tight_layout()
    #plt.savefig('../haskins-group1/users/vsun/USOS_shared/Plotting/USOS_Campaign_analysis/plots/' + plot_name)
    plt.show()

def species_by_day_of_week_monthly_split_plot(species_var, species_name, y_limit, y_ticks, color_map1, color_map2, plot_name):
    day_of_week_species_mean_split_monthly = df_alldays.groupby([month_index, day_names_index, hr_of_day_index])[species_var].mean()

    fig, ax = plt.subplots(1,2, figsize = (20,10))

    cmap1 = mpl.colormaps[color_map1]
    cmap2 = mpl.colormaps[color_map2]

    # Create 20 color values for gradient in the colormap
    color_values = np.linspace(0.2, 1, 7) 

    #Loop through each day for July and plot the ozone values on first subplot
    for day_val in range(0,7):
        toplot, = ax[0].plot(hour_range, day_of_week_species_mean_split_monthly[7][days_of_week[day_val]], color=cmap1(color_values[day_val]), label = days_of_week[day_val])
        ax[0].set_xlabel('Hour')
        ax[0].set_ylabel(species_name + ' Concentration (ppb)')
        ax[0].set_xticks(np.arange(0,24))
        ax[0].set_ylim(y_limit)
        ax[0].set_yticks(y_ticks)
        ax[0].set_title('July ' + species_name)

    #get labels for legend
    handles, labels = ax[0].get_legend_handles_labels()

    ax[0].legend(handles, labels)
    ax[0].margins(x=0)
    ax[0].grid()

    #Loop through each day for August and plot the ozone values on second subplot
    for day_val in range(0,7):
        toplot, = ax[1].plot(hour_range, day_of_week_species_mean_split_monthly[8][days_of_week[day_val]], color=cmap2(color_values[day_val]), label = days_of_week[day_val])
        ax[1].set_xlabel('Hour')
        ax[1].set_ylabel(species_name + ' Concentration (ppb)')
        ax[1].set_xticks(np.arange(0,24))
        ax[1].set_ylim(y_limit)
        ax[1].set_yticks(y_ticks)
        ax[1].set_title('August ' + species_name)

    #get labels for legend
    handles, labels = ax[1].get_legend_handles_labels()

    ax[1].legend(handles, labels)
    ax[1].margins(x=0)
    ax[1].grid()

    plt.tight_layout()
    #plt.savefig('../haskins-group1/users/vsun/USOS_shared/Plotting/USOS_Campaign_analysis/plots/' + plot_name)
    plt.show()

def species_campaign_total_plot(species_var1, species_var2, species_name1, species_name2,  y_limit1, y_ticks1,y_limit2, y_ticks2, color_map1, color_map2,plot_name):
    #shows half hour intervals 

    df_nox = df_alldays[species_var1]
    overall_nox_avg = df_nox.groupby([hr_of_day_index]).mean()
    overall_nox_std = np.std(overall_nox_avg)

    #shows half hour intervals of ozone
    df_ozone = df_alldays[species_var2]

    overall_ozone_avg = df_ozone.groupby([hr_of_day_index]).mean()
    overall_ozone_std = np.std(overall_ozone_avg)

    overall_dates = pd.date_range(start='2024-07-16', end='2024-08-18')
    overall_dates_list = overall_dates.strftime('%Y-%m-%d')
    
    fig, ax = plt.subplots(1,2, figsize = (20,10))
    cmap1 = mpl.colormaps[color_map1]
    cmap2 = mpl.colormaps[color_map2]

    # Create 20 color values for gradient in the colormap
    color_values = np.linspace(0.15, 0.95, 34) 

    #Loop through each day and plot the ozone values on first subplot
    for day_val in range(0, 34):
        ax[0].plot(hour_range, df_nox[overall_dates_list[day_val]], color=cmap1(color_values[day_val]), label = overall_dates_list[day_val])
        ax[0].set_xlabel('Hour')
        ax[0].set_ylabel(species_name1 + ' Concentration (ppb)')
        ax[0].set_xticks(np.arange(0,24))
        ax[0].set_ylim(y_limit1)
        ax[0].set_yticks(y_ticks1)
        ax[0].set_title('USOS ' + species_name1)

        ax[1].plot(hour_range, df_ozone[overall_dates_list[day_val]], color=cmap2(color_values[day_val]), label = overall_dates_list[day_val])
        ax[1].set_xlabel('Hour')
        ax[1].set_ylabel(species_name2 + ' Concentration (ppb)')
        ax[1].set_xticks(np.arange(0,24))
        ax[1].set_ylim(y_limit2)
        ax[1].set_yticks(y_ticks2)
        ax[1].set_title('USOS ' + species_name2)
    ax[0].plot(hour_range, overall_nox_avg, color = 'k', label = 'Campaign Average')
    ax[0].fill_between(hour_range, overall_nox_avg - overall_nox_std, overall_nox_avg + overall_nox_std, color='k',alpha = 0.1)
    
    ax[1].plot(hour_range, overall_ozone_avg, color = 'k', label = 'Campaign Average')
    ax[1].fill_between(hour_range, overall_ozone_avg - overall_ozone_std, overall_ozone_avg + overall_ozone_std, color='k',alpha = 0.1)
    ax[1].hlines(y=70, xmin = 0, xmax = 23.5, color = 'tab:red')
    
def single_species_plot(ozone_filepath):
    df_ozone = pd.read_csv(ozone_filepath)
    date_index = pd.date_range(start='2024-07-16', periods=len(df_ozone), freq='1h')
    df_ozone.index = date_index
    df_merged_ozone =  df_ozone['merged_ozone']
    month_index = df_ozone.index.month
    date_index = df_ozone.index.date
    day_names_index =df_ozone.index.day_name()
    hr_of_day_index = df_ozone.index.time

    # mda8_exceedance = np.where(df_ozone['Exceedance_day']==True, df_ozone['MDA8_O3'].values , np.nan)
    # print(mda8_exceedance)

    overall_ozone_avg = df_ozone['merged_ozone'].groupby([hr_of_day_index]).mean()
    overall_ozone_std = np.std(overall_ozone_avg)

    overall_dates = pd.date_range(start='2024-07-16', end='2024-08-18')
    overall_dates_list = overall_dates.strftime('%Y-%m-%d')

    color_values = np.linspace(0.2, 1, 34) 
    cmap1 = mpl.colormaps['YlGnBu']

    # Get unique days
    unique_days = df_ozone.index.normalize().unique()
    num_days = len(unique_days)
    print(unique_days)

    fig = plt.figure(figsize=(10, 10))
    for i, day in enumerate(unique_days):
        mask = df_ozone.index.normalize() == day
        day_df = df_ozone[mask]
        # print(df_ozone.index.date)
        # print(day)
        # day_df = df_ozone[df_ozone.index.date == day]
        # print(day_df)
        #x = day_df.index.hour
        linestyle = '-' if day_df['Exceedance_day'].any() else '--'
        label = str(day.date()) if i % 3 == 0 else None  # label every 3rd day, else no label
        plt.plot(hour_range, day_df['merged_ozone'], linestyle=linestyle, color=cmap1(color_values[i]), label=label)
    plt.plot(hour_range,overall_ozone_avg, linestyle='solid', linewidth=2, color='k')
    plt.fill_between(hour_range, overall_ozone_avg - overall_ozone_std, overall_ozone_avg + overall_ozone_std, color='k',alpha = 0.1)
    plt.hlines(y=70, xmin = 0, xmax = 23, color = 'tab:red')
    plt.xticks(np.arange(0,24,1))
    plt.xlabel('Hour')
    plt.ylabel('Ozone (ppbv)')
    plt.title('Ozone at Hawthorne during USOS')
    plt.legend()
    plt.margins(x=0)
    plt.grid()

    plt.tight_layout()

    true_count = df_ozone['Exceedance_day'].sum()
    print(true_count)
    print(len(df_ozone.index))
    # for day, day_df in df_ozone.groupby(df_ozone.index.date):
    #     print(day_df.index)
    # # Decide line style for the whole day:
    # # e.g. if ANY condition True -> solid, else dashed
    # # You can change this logic as needed
    #     if day_df['Exceedance_day'].any():
    #         linestyle = '-'
    #     else:
    #         linestyle = '--'

    #     plt.plot(hour_range, day_df['merged_ozone'], linestyle=linestyle, label=str(day), color=cmap1(color_values[day]))


    # Plot each segment with appropriate style
    # plt.figure(figsize=(10, 5))
    # for _, segment_df in df_ozone.groupby('segment'):
    #     linestyle = '-' if segment_df['Exceedance_day'].iloc[0] else '--'
    #     plt.plot(segment_df.index, segment_df['merged_O3'], linestyle=linestyle, color='blue')

    
    # for day_val in range(0, 34):
    #     plt.plot(hour_range, df_ozone['merged_ozone'][overall_dates_list[day_val]], color=cmap1(color_values[day_val]), label = overall_dates_list[day_val], linestyle = 'dashed')
    #     # plt.plot(df_ozone.index, mda8_exceedance, color='r', linestyle = 'solid')
    #     plt.xlabel('Hour')
    #     plt.ylabel('Ozone Concentration (ppb)')
    #     #plt.xticks(np.arange(0,24))
    #     #plt.ylim(y_limit1)
    #     #plt.yticks(y_ticks1)
    #     plt.title('USOS Ozone')


#####
#Function calls

# Ozone
# ozone_campaign_total_monthly_split_plot(
#     species_var = 'O3_ppbv',
#     species_name = 'Ozone',
#     y_limit = [0,100],
#     y_ticks = np.arange(0,100,5),
#     color_map1 = 'YlOrRd',
#     color_map2 = 'YlGnBu',
#     plot_name = 'ozone_monthly_split.png' 
# )

# species_by_day_of_week_plot(
#     species_var = 'O3_ppbv',
#     species_name = 'Ozone',
#     y_limit = [0,85],
#     y_ticks = np.arange(0,85,5),
#     color_map = 'Greens',
#     plot_name = 'ozone_day_of_week.png'
# )

# species_by_day_of_week_monthly_split_plot(
#     species_var = 'O3_ppbv',
#     species_name = 'Ozone',
#     y_limit = [0,85],
#     y_ticks = np.arange(0,85,5),
#     color_map1 = 'YlOrRd',
#     color_map2 = 'YlGnBu',
#     plot_name = 'ozone_day_of_week_monthly_split.png'
# )

# # NOx
# species_by_day_of_week_plot(
#     species_var = 'NOx_ppbv',
#     species_name = 'NOx',
#     y_limit = [0,50],
#     y_ticks = np.arange(0,55,5),
#     color_map = 'Greens',
#     plot_name = 'ozone_day_of_week.png'
# )

# species_campaign_total_monthly_split_plot(
#     species_var = 'NOx_ppbv',
#     species_name = 'NOx',
#     y_limit = [0,50],
#     y_ticks = np.arange(0,50,5),
#     color_map1 = 'YlOrRd',
#     color_map2 = 'YlGnBu',
#     plot_name = 'nox_monthly_split.png'
# )

# #Acetonitrile
# species_campaign_total_monthly_split_plot(
#     species_var = 'Acetonitrile_PTR',
#     species_name = 'Acetonitrile',
#     y_limit = [0,1],
#     y_ticks = np.arange(0,1,0.05),
#     color_map1 = 'YlOrRd',
#     color_map2 = 'YlGnBu',
#     plot_name = 'acetonitrile_monthly_split.png'  
# )

# species_campaign_total_plot(
#     species_var1 = 'NOx_ppbv',
#     species_var2 = 'O3_ppbv',
#     species_name1 = 'NOx',
#     species_name2 = 'Ozone',
#     y_limit1 = [0,45],
#     y_limit2 = [10,100],
#     y_ticks1 = np.arange(0,45,5),
#     y_ticks2 = np.arange(10,100,5),
#     color_map1 = 'YlGnBu',
#     color_map2 = 'YlOrRd',
#     plot_name = 'nox_ozone_usos.png'
# )
#     
single_species_plot(
    ozone_filepath = dirpath + '/Plotting/USOS_Campaign_analysis/all_ozone_during_usos_with_exceedances.csv'
)
#     
#     
#    