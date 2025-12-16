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
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerBase
import matplotlib.patches as mpatches
from mpl_axes_aligner import align
import cmasher as cmr

import sys 
import re 
import yaml
import inspect 
from collections import defaultdict

dirpath = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/'

##### Global variables #####
days_of_week = ['Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Saturday', 'Sunday']
weekday_type_list = ['Weekday', 'Weekend']
weekday_type_list_sunday = ['Weekday','Sunday']
hour_range = np.arange(0,24,1)

##### Load files #####
#Combined file with USOS Ozone, UDAQ Ozone, Merge of USOS & UDAQ
#Includes columns for MDA8 and Exceedance Days for all 3 measurement types
ozone_filepath = dirpath + '/Plotting/USOS_Campaign_analysis/all_ozone_during_usos_with_exceedances.csv'
df_ozone = pd.read_csv(ozone_filepath, index_col='dt', parse_dates=True)

#All species from USOS campaign
all_days_filepath = dirpath + 'CampaignData_and_Merges/R0/CSL_MobileLab_Parked/merged/rev_30min/all_CSL_MobileLab_Parked_rev30minv4.nc'
all_days_filepath_load = xr.open_dataset(all_days_filepath)
df_alldays = all_days_filepath_load.to_dataframe()
df_alldays.reset_index(inplace=True)
df_alldays.set_index('time_local', inplace=True, drop=False)
df_alldays_resample = df_alldays.resample('1h').mean()
#adds padding of NaNs for first day and last day of campaign, since they need to have the same length of time as the other days in order to plot
new_start_time = pd.Timestamp('2024-07-15 00:00:00')
new_end_time = pd.Timestamp('2024-08-18 23:00:00')
# Create a new datetime index from new_start to the end of existing index with same frequency
new_index = pd.date_range(start=new_start_time, end=new_end_time, freq='1h')
# Reindex the dataframe to include new rows
df_alldays = df_alldays_resample.reindex(new_index)
#add NO and NO2 to get NOx column
df_alldays['NOx_ppbv'] = df_alldays['NO_LIF'] + df_alldays['NO2_LIF']

##### Datetime and Date Index information
month_index_allspecies = df_alldays.index.month
date_index_allspecies = df_alldays.index.date
day_names_index_allspecies = df_alldays.index.day_name()
hr_of_day_index_allspecies = df_alldays.index.time
weekday_index_allspecies = df_alldays.index.weekday
df_alldays['Weekday'] = np.where(weekday_index_allspecies < 5, 'Weekday', 'Weekend')
df_alldays['Sunday'] = np.where(weekday_index_allspecies == 6, 'Sunday', 'Not-Sunday')

month_index_ozone = df_ozone.index.month
date_index_ozone = df_ozone.index.date
day_names_index_ozone =df_ozone.index.day_name()
hr_of_day_index_ozone = df_ozone.index.time
weekday_index_ozone = df_ozone.index.weekday
df_ozone['Weekday'] = np.where(weekday_index_ozone < 5, 'Weekday', 'Weekend')
df_ozone['Sunday'] = np.where(weekday_index_ozone == 6, 'Sunday', 'Not-Sunday')

overall_dates = pd.date_range(start='2024-07-15', end='2024-08-18')
overall_dates_list = overall_dates.strftime('%Y-%m-%d')

#DatetimeIndex of days of the campaign from 07-15 to 08-18
unique_days = df_ozone.index.normalize().unique()

##### Averages used for various plots #####
#USOS Ozone only
overall_ozone_avg_usos = df_ozone['Mobile Lab O3'].groupby([hr_of_day_index_ozone]).mean()
overall_ozone_std_usos = np.std(overall_ozone_avg_usos)
day_of_week_ozone_mean_usos = df_ozone.groupby([day_names_index_ozone, hr_of_day_index_ozone])['Mobile Lab O3'].mean()
weekday_ozone_mean_usos = df_ozone.groupby(['Weekday', hr_of_day_index_ozone])['Mobile Lab O3'].mean()

#Merged Ozone
overall_ozone_avg_merged = df_ozone['merged_ozone'].groupby([hr_of_day_index_ozone]).mean()
overall_ozone_std_merged = np.std(overall_ozone_avg_merged)
day_of_week_ozone_mean_merged = df_ozone.groupby([day_names_index_ozone, hr_of_day_index_ozone])['merged_ozone'].mean()
weekday_ozone_mean_merged = df_ozone.groupby(['Weekday', hr_of_day_index_ozone])['merged_ozone'].mean()

#NOx
df_nox = df_alldays['NOx_ppbv']
overall_nox_avg = df_nox.groupby([hr_of_day_index_allspecies]).mean()
overall_nox_std = np.std(overall_nox_avg)
day_of_week_allspecies_mean = df_alldays.groupby([day_names_index_allspecies, hr_of_day_index_allspecies])['NOx_ppbv'].mean()
weekday_allspecies_mean = df_alldays.groupby(['Weekday', hr_of_day_index_allspecies])['NOx_ppbv'].mean()

##### Plotting Formatting #####
mpl.rcParams['xtick.labelsize'] = 15
mpl.rcParams['ytick.labelsize'] = 15
mpl.rcParams['legend.fontsize'] = 16
mpl.rcParams['axes.labelsize'] = 18
mpl.rcParams['axes.titlesize'] = 28
mpl.rcParams['axes.xmargin'] = 0

# Set the font family to 'serif'
mpl.rcParams['font.family'] = 'serif'
# Specify 'cmr10' as the preferred serif font (Computer Modern Roman)
mpl.rcParams['font.serif'] = 'cmr10' 
# Optionally, configure mathtext to use Computer Modern fonts as well
mpl.rcParams['mathtext.fontset'] = 'cm'
# Ensure minus signs are rendered correctly with CM fonts
mpl.rcParams['axes.unicode_minus'] = False

def ozone_daily(colormap, o3_measured_from, o3_avg_from, o3_std_from, merged_or_usos):
    color_values = np.linspace(0.3, 1.0, 39) 
    cmap1 = mpl.colormaps[colormap]
    fig, ax = plt.subplots(figsize=(10, 10), layout = 'constrained')
    for i, day in enumerate(unique_days):
        mask = df_ozone.index.normalize() == day
        day_df = df_ozone[mask]
        exceedance_string = 'Exceedance_day_' + o3_measured_from
        dashed_linestyle = '-' if day_df[exceedance_string].any() else '--'
        alpha_amount = 1 if day_df[exceedance_string].any() else 0.8
        exceedance_labels = str(day.date()) if day_df[exceedance_string].any() else '_'
        #label = str(day.date()) if i % 3 == 0 else None  # label every 3rd day, else no label
        plt.plot(hour_range, day_df[o3_measured_from], linestyle=dashed_linestyle, color=cmap1(color_values[i]), label=exceedance_labels, alpha = alpha_amount)

    plt.plot(hour_range, o3_avg_from, linestyle='solid', linewidth=2, color='k')
    plt.fill_between(hour_range, o3_avg_from - o3_std_from, o3_avg_from + o3_std_from, color='k',alpha = 0.1)
    plt.hlines(y=70, xmin = 0, xmax = 23, color = 'tab:red')

    plt.xlabel('Hour')
    ax.xaxis.set_minor_locator(MultipleLocator(3))
    ax.xaxis.set_major_locator(MultipleLocator(6))
    plt.xticks()
    plt.margins(x=0)

    plt.ylabel('Ozone Concentration (ppbv)')
    plt.title('Ozone at Hawthorne during USOS (' + merged_or_usos + ')')

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    #add grid
    ax.grid(which='both', linestyle='--', linewidth=0.5)

    plt.savefig(dirpath + '/Plotting/USOS_Campaign_analysis/plots/hawthorne_ozone_alldays_' + merged_or_usos + '.png', dpi =300)
    plt.show()

    true_count = df_ozone[exceedance_string].sum()
    print('Exceedance days ', int(true_count/24), ' out of ', int(len(df_ozone.index)/24)-2, ' days of ozone measurements (2024-07-15 to 2024-08-18)')

def nox_ozone_daily(colormap1, colormap2, o3_measured_from, o3_avg_from, o3_std_from, merged_or_usos):
    fig, ax = plt.subplots(1,2, figsize = (20,10))

    # Create color values for gradient in the colormap
    color_values = np.linspace(0.15, 0.95, 35) 
    cmap1 = mpl.colormaps[colormap1]
    cmap2 = mpl.colormaps[colormap2]

    for i, day in enumerate(unique_days):
        mask = df_ozone.index.normalize() == day
        day_df_o3 = df_ozone[mask]
        day_df_nox = df_nox[mask]
        exceedance_string = 'Exceedance_day_' + o3_measured_from
        dashed_linestyle = '-' if day_df_o3[exceedance_string].any() else '--'
        alpha_amount = 1 if day_df_o3[exceedance_string].any() else 0.8
        exceedance_labels = str(day.date()) if day_df_o3[exceedance_string].any() else '_'
        #label = str(day.date()) if i % 3 == 0 else None  # label every 3rd day, else no label

        ax[0].plot(hour_range, day_df_nox, linestyle=dashed_linestyle, color=cmap1(color_values[i]), label=exceedance_labels, alpha = alpha_amount)
        ax[0].set_xlabel('Hour')
        ax[0].set_ylabel('NOx Concentration (ppb)')
        ax[0].set_ylim([0,40])
        ax[0].set_title('USOS NOx (' + merged_or_usos  + ')')
        
        ax[1].plot(hour_range, day_df_o3[o3_measured_from], linestyle=dashed_linestyle, color=cmap2(color_values[i]), label=exceedance_labels, alpha = alpha_amount)
        ax[1].set_xlabel('Hour')
        ax[1].set_ylabel('Ozone Concentration (ppb)')
        ax[1].set_ylim([10,100])
        ax[1].set_title('USOS Ozone (' + merged_or_usos + ')')
    ax[0].plot(hour_range, overall_nox_avg, color = 'k', label = 'Campaign Average')
    ax[0].fill_between(hour_range, overall_nox_avg - overall_nox_std, overall_nox_avg + overall_nox_std, color='k',alpha = 0.1)
    ax[0].xaxis.set_minor_locator(MultipleLocator(3))
    ax[0].xaxis.set_major_locator(MultipleLocator(6))
    ax[0].grid(which='both', linestyle='--', linewidth=0.5)

    ax[1].plot(hour_range, o3_avg_from, color = 'k', label = 'Campaign Average')
    ax[1].fill_between(hour_range, o3_avg_from - o3_std_from, o3_avg_from + o3_std_from, color='k',alpha = 0.1)
    ax[1].xaxis.set_minor_locator(MultipleLocator(3))
    ax[1].xaxis.set_major_locator(MultipleLocator(6))
    ax[1].grid(which='both', linestyle='--', linewidth=0.5)
    ax[1].hlines(y=70, xmin = 0, xmax = 23, color = 'tab:red')

    plt.savefig(dirpath+'/Plotting/USOS_Campaign_analysis/plots/hawthorne_nox_ozone_comparison_'+ merged_or_usos +'.png', dpi =300)
    plt.show()

def leighton_ratio_daily():
    df_alldays['NO_LIF'] = df_alldays['NO_LIF'].mask(df_alldays['NO_LIF'] < 0)
    pressure_mb = df_alldays['Press_mb_POPS']
    temp_celsius = df_alldays['Temp_C_POPS']
    photolysis_no2_meas = df_alldays['jNO2_meas']
    no_mobilelab = df_alldays['NO_LIF']
    no2_mobilelab = df_alldays['NO2_LIF']
    o3_mobilelab = df_alldays['O3_ppbv']
    df_alldays['NOx_ppbv'] = df_alldays['NO_LIF'] + df_alldays['NO2_LIF']

    #ppbv to molecs/cm^3 conversion
    pressure_pa = pressure_mb * 10**2
    temp_kelvin = temp_celsius + 273.15
    gas_constant = 8.31446261815324 * 10**6 #cm^3 * Pa * K^−1 * mol^−1

    #70 ppbv= 70x10^-9 mol/mol 
    volume_mixing_ratio_no = no_mobilelab/10**9
    volume_mixing_ratio_no2 = no2_mobilelab/10**9
    volume_mixing_ratio_o3 = o3_mobilelab/10**9
    num_density_air = ((6.022*10**23) * pressure_pa) / (gas_constant * temp_kelvin)
    num_density_no =  num_density_air * volume_mixing_ratio_no
    num_density_no2 = num_density_air * volume_mixing_ratio_no2
    num_density_o3 = num_density_air * volume_mixing_ratio_o3

    df_alldays['num_density_no'] = num_density_no
    df_alldays['num_density_no2'] = num_density_no2
    df_alldays['num_density_nox'] = num_density_no + num_density_no2
    df_alldays['num_density_o3'] = num_density_o3

    df_alldays['exceedances mobile lab'] = df_ozone['Exceedance_day_Mobile Lab O3']

    #based off JPL kinetics publication
    k_rate_o3_no = (3.0*10**-12) * np.exp(-1500*(1/temp_kelvin))

    prod_o3_num = photolysis_no2_meas * num_density_no2
    loss_o3_den = k_rate_o3_no * num_density_no * num_density_o3
    leighton_ratio = prod_o3_num / loss_o3_den
    no_no2_ratio = (df_alldays['num_density_no'])/(df_alldays['num_density_no2'])

    df_alldays['prod_o3_num'] = prod_o3_num
    df_alldays['loss_o3_den'] = loss_o3_den
    df_alldays['leighton ratio'] = leighton_ratio
    df_alldays['no_no2_ratio'] = no_no2_ratio

    overall_leighton_ratio_avg_mobile_lab = df_alldays['leighton ratio'].groupby([hr_of_day_index_allspecies]).mean()
    overall_leighton_ratio_std_mobile_lab = np.std(overall_leighton_ratio_avg_mobile_lab)

    fig, ax = plt.subplots(figsize=(10, 10), layout = 'tight')
    colors = cmr.take_cmap_colors('cmr.freeze_r', 35, cmap_range=(0.1, 0.8), return_fmt='hex')

    for i, day in enumerate(unique_days):
        mask = df_alldays.index.normalize() == day
        day_df = df_alldays[mask]
        exceedance_string = 'exceedances mobile lab'
        dashed_linestyle = '-' if day_df[exceedance_string].any() else '--'
        alpha_amount = 1 if day_df[exceedance_string].any() else 0.8
        exceedance_labels = str(day.date()) if day_df[exceedance_string].any() else '_'
        plt.plot(hour_range, day_df['leighton ratio'], linestyle=dashed_linestyle, 
                 color=colors[i], label=exceedance_labels, alpha = alpha_amount)

    plt.plot(hour_range, overall_leighton_ratio_avg_mobile_lab, linestyle='solid', linewidth=2, color='k')
    plt.fill_between(hour_range, overall_leighton_ratio_avg_mobile_lab - overall_leighton_ratio_std_mobile_lab, 
                     overall_leighton_ratio_avg_mobile_lab + overall_leighton_ratio_std_mobile_lab, color='k',alpha = 0.1)

    plt.xlabel('Hour')
    ax.xaxis.set_minor_locator(MultipleLocator(3))
    ax.xaxis.set_major_locator(MultipleLocator(6))
    plt.xticks()
    plt.margins(x=0)

    plt.ylabel('Leighton Ratio')
    plt.title('Leighton Ratio at Hawthorne during USOS')

    # box = ax.get_position()
    # ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='upper left')

    #add grid
    ax.grid(which='both', linestyle='--', linewidth=0.5)

    plt.show()

def leighton_ratio_daily_exceedances_only():
    df_alldays['NO_LIF'] = df_alldays['NO_LIF'].mask(df_alldays['NO_LIF'] < 0)
    pressure_mb = df_alldays['Press_mb_POPS']
    temp_celsius = df_alldays['Temp_C_POPS']
    photolysis_no2_meas = df_alldays['jNO2_meas']
    no_mobilelab = df_alldays['NO_LIF']
    no2_mobilelab = df_alldays['NO2_LIF']
    o3_mobilelab = df_alldays['O3_ppbv']
    df_alldays['NOx_ppbv'] = df_alldays['NO_LIF'] + df_alldays['NO2_LIF']

    #ppbv to molecs/cm^3 conversion
    pressure_pa = pressure_mb * 10**2
    temp_kelvin = temp_celsius + 273.15
    gas_constant = 8.31446261815324 * 10**6 #cm^3 * Pa * K^−1 * mol^−1

    #70 ppbv= 70x10^-9 mol/mol 
    volume_mixing_ratio_no = no_mobilelab/10**9
    volume_mixing_ratio_no2 = no2_mobilelab/10**9
    volume_mixing_ratio_o3 = o3_mobilelab/10**9
    num_density_air = ((6.022*10**23) * pressure_pa) / (gas_constant * temp_kelvin)
    num_density_no =  num_density_air * volume_mixing_ratio_no
    num_density_no2 = num_density_air * volume_mixing_ratio_no2
    num_density_o3 = num_density_air * volume_mixing_ratio_o3

    df_alldays['num_density_no'] = num_density_no
    df_alldays['num_density_no2'] = num_density_no2
    df_alldays['num_density_nox'] = num_density_no + num_density_no2
    df_alldays['num_density_o3'] = num_density_o3

    df_alldays['exceedances mobile lab'] = df_ozone['Exceedance_day_Mobile Lab O3']

    print(df_alldays['exceedances mobile lab'])
    #based off JPL kinetics publication
    k_rate_o3_no = (3.0*10**-12) * np.exp(-1500*(1/temp_kelvin))

    prod_o3_num = photolysis_no2_meas * num_density_no2
    loss_o3_den = k_rate_o3_no * num_density_no * num_density_o3
    leighton_ratio = prod_o3_num / loss_o3_den
    no_no2_ratio = (df_alldays['num_density_no'])/(df_alldays['num_density_no2'])

    df_alldays['prod_o3_num'] = prod_o3_num
    df_alldays['loss_o3_den'] = loss_o3_den
    df_alldays['leighton ratio'] = leighton_ratio
    df_alldays['no_no2_ratio'] = no_no2_ratio

    overall_leighton_ratio_avg_mobile_lab = df_alldays['leighton ratio'].groupby([hr_of_day_index_allspecies]).mean()
    overall_leighton_ratio_std_mobile_lab = np.std(overall_leighton_ratio_avg_mobile_lab)

    # Extract date and hour for convenience
    df_alldays['date'] = df_alldays.index.date
    df_alldays['hour'] = df_alldays.index.hour
    # Group by each day
    grouped = df_alldays.groupby('date')

    fig, ax = plt.subplots(figsize=(10, 10), layout = 'tight')
    colors = cmr.take_cmap_colors('cmr.freeze_r', 8, cmap_range=(0.1, 0.8), return_fmt='hex')

    # Counter for assigning colors
    color_counter = 0
    for date, group in grouped:
        # Check if *any* hour of this day has 'include' == True
        if group['exceedances mobile lab'].any():
            color_chosen = colors[color_counter % len(colors)]  # cycle if more days than colors
            # Plot hour vs value
            plt.plot(group['hour'], group['leighton ratio'], color = color_chosen,
                     label=str(date))
            color_counter += 1
            
            ###### Leighton Ratio > 1 counter #####
            # threshold = 1.0
            # count = (group['leighton ratio'] > threshold).sum()
            # print('Hours with Leighton Ratio > 1 ', date,', ', count)

    plt.plot(hour_range, overall_leighton_ratio_avg_mobile_lab, linestyle='solid', linewidth=2, color='k', label = 'Overall Avg')
    plt.fill_between(hour_range, overall_leighton_ratio_avg_mobile_lab - overall_leighton_ratio_std_mobile_lab, 
                     overall_leighton_ratio_avg_mobile_lab + overall_leighton_ratio_std_mobile_lab, color='k',alpha = 0.1)

    plt.xlabel('Hour')
    ax.xaxis.set_minor_locator(MultipleLocator(3))
    ax.xaxis.set_major_locator(MultipleLocator(6))
    plt.xticks()
    plt.margins(x=0)

    plt.ylabel('Leighton Ratio')
    plt.title('Leighton Ratio at Hawthorne \n During Exceedance Days')

    # box = ax.get_position()
    # ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='upper left')

    #add grid
    ax.grid(which='both', linestyle='--', linewidth=0.5)

    plt.savefig(dirpath + '/Plotting/USOS_Campaign_analysis/plots/leighton_ratio_during_exceedances.png', dpi =150)
    plt.show()



def leighton_ratio_daily_exceedance_comparison():
    df_alldays['NO_LIF'] = df_alldays['NO_LIF'].mask(df_alldays['NO_LIF'] < 0)
    pressure_mb = df_alldays['Press_mb_POPS']
    temp_celsius = df_alldays['Temp_C_POPS']
    photolysis_no2_meas = df_alldays['jNO2_meas']
    no_mobilelab = df_alldays['NO_LIF']
    no2_mobilelab = df_alldays['NO2_LIF']
    o3_mobilelab = df_alldays['O3_ppbv']
    df_alldays['NOx_ppbv'] = df_alldays['NO_LIF'] + df_alldays['NO2_LIF']

    #ppbv to molecs/cm^3 conversion
    pressure_pa = pressure_mb * 10**2
    temp_kelvin = temp_celsius + 273.15
    gas_constant = 8.31446261815324 * 10**6 #cm^3 * Pa * K^−1 * mol^−1

    #70 ppbv= 70x10^-9 mol/mol 
    volume_mixing_ratio_no = no_mobilelab/10**9
    volume_mixing_ratio_no2 = no2_mobilelab/10**9
    volume_mixing_ratio_o3 = o3_mobilelab/10**9
    num_density_air = ((6.022*10**23) * pressure_pa) / (gas_constant * temp_kelvin)
    num_density_no =  num_density_air * volume_mixing_ratio_no
    num_density_no2 = num_density_air * volume_mixing_ratio_no2
    num_density_o3 = num_density_air * volume_mixing_ratio_o3

    df_alldays['num_density_no'] = num_density_no
    df_alldays['num_density_no2'] = num_density_no2
    df_alldays['num_density_nox'] = num_density_no + num_density_no2
    df_alldays['num_density_o3'] = num_density_o3

    df_alldays['exceedances mobile lab'] = df_ozone['Exceedance_day_Mobile Lab O3']

    #based off JPL kinetics publication
    k_rate_o3_no = (3.0*10**-12) * np.exp(-1500*(1/temp_kelvin))

    prod_o3_num = photolysis_no2_meas * num_density_no2
    loss_o3_den = k_rate_o3_no * num_density_no * num_density_o3
    leighton_ratio = prod_o3_num / loss_o3_den
    no_no2_ratio = (df_alldays['num_density_no'])/(df_alldays['num_density_no2'])

    df_alldays['prod_o3_num'] = prod_o3_num
    df_alldays['loss_o3_den'] = loss_o3_den
    df_alldays['leighton ratio'] = leighton_ratio
    df_alldays['no_no2_ratio'] = no_no2_ratio

    overall_leighton_ratio_avg_mobile_lab = df_alldays['leighton ratio'].groupby([hr_of_day_index_allspecies]).mean()
    overall_leighton_ratio_std_mobile_lab = np.std(overall_leighton_ratio_avg_mobile_lab)

    # Extract date and hour for convenience
    df_alldays['date'] = df_alldays.index.date
    df_alldays['hour'] = df_alldays.index.hour
    # Group by each day
    grouped = df_alldays.groupby('date')

    fig, ax = plt.subplots(1,2, figsize = (20,10), layout = 'tight')
    colors1 = cmr.take_cmap_colors('cmr.freeze_r', 8, cmap_range=(0.1, 0.8), return_fmt='hex')
    colors2 = cmr.take_cmap_colors('cmr.freeze_r', 26, cmap_range=(0.1, 0.8), return_fmt='hex')

    # Counter for assigning colors
    color_counter1 = 0
    color_counter2 = 0
    for date, group in grouped:
        # Check if *any* hour of this day has 'include' == True
        if group['exceedances mobile lab'].any():
            color_chosen1 = colors1[color_counter1 % len(colors1)]  # cycle if more days than colors
            # Plot hour vs value
            ax[0].plot(group['hour'], group['leighton ratio'], color = color_chosen1,
                     label=str(date))
            color_counter1 += 1

            ###### Leighton Ratio > 1 counter #####
            threshold = 1.0
            count = (group['leighton ratio'] > threshold).sum()
            print('Hours with Leighton Ratio > 1, exceedance day, ', date,', ', count)

        else:
            color_chosen2 = colors2[color_counter2 % len(colors2)]  # cycle if more days than colors
            ax[1].plot(group['hour'], group['leighton ratio'], color = color_chosen2,
                     label=str(date))
            color_counter2 += 1
            ###### Leighton Ratio > 1 counter #####
            threshold = 1.0
            count = (group['leighton ratio'] > threshold).sum()
            print('Hours with Leighton Ratio > 1, non-exceedance day, ', date,', ', count)

    ax[0].plot(hour_range, overall_leighton_ratio_avg_mobile_lab, linestyle='solid', linewidth=2, color='k', label = 'Overall Avg')
    ax[0].fill_between(hour_range, overall_leighton_ratio_avg_mobile_lab - overall_leighton_ratio_std_mobile_lab, 
                     overall_leighton_ratio_avg_mobile_lab + overall_leighton_ratio_std_mobile_lab, color='k',alpha = 0.1)
    
    ax[1].plot(hour_range, overall_leighton_ratio_avg_mobile_lab, linestyle='solid', linewidth=2, color='k', label = 'Overall Avg')
    ax[1].fill_between(hour_range, overall_leighton_ratio_avg_mobile_lab - overall_leighton_ratio_std_mobile_lab, 
                     overall_leighton_ratio_avg_mobile_lab + overall_leighton_ratio_std_mobile_lab, color='k',alpha = 0.1)

    ax[0].set_xlabel('Hour')
    ax[0].xaxis.set_minor_locator(MultipleLocator(3))
    ax[0].xaxis.set_major_locator(MultipleLocator(6))
    ax[0].margins(x=0)
    ax[0].set_ylabel('Leighton Ratio')
    ax[0].set_title('Leighton Ratio at Hawthorne \n During Exceedance Days')
    ax[0].legend(loc='upper left')
    #add grid
    ax[0].grid(which='both', linestyle='--', linewidth=0.5)

    ax[1].set_xlabel('Hour')
    ax[1].xaxis.set_minor_locator(MultipleLocator(3))
    ax[1].xaxis.set_major_locator(MultipleLocator(6))
    ax[1].margins(x=0)
    ax[1].set_ylabel('Leighton Ratio')
    ax[1].set_title('Leighton Ratio at Hawthorne \n During Non-Exceedance Days')
    ax[1].legend(loc='upper left')
    #add grid
    ax[1].grid(which='both', linestyle='--', linewidth=0.5)

    plt.savefig(dirpath + '/Plotting/USOS_Campaign_analysis/plots/leighton_ratio_during_exceedances.png', dpi =150)
    plt.show()


def day_of_week_nox_ozone(colormap1, colormap2, day_of_week_ozone_from, o3_avg_from, o3_std_from, merged_or_usos):
    fig, ax = plt.subplots(1,2, figsize = (20,10))
    # Create 20 color values for gradient in the colormap
    color_values_week = np.linspace(0.2, 1, 8)

    cmap1 = mpl.colormaps[colormap1]
    cmap2 = mpl.colormaps[colormap2]

    class TwoColorHandler(HandlerBase):
        def create_artists(self, legend, orig_handle,
                        xdescent, ydescent, width, height, fontsize, trans):
            patch_width = width / 2.5  # width of each square
            patch_height = height * 0.6  # height of each square
            spacing = 4
            y = ydescent + (height - patch_height) / 2

            # Two colored rectangles side by side
            patch1 = mpatches.Rectangle((xdescent, y), patch_width, patch_height,
                                        facecolor=orig_handle[0], edgecolor=orig_handle[0], transform=trans)
            patch2 = mpatches.Rectangle((xdescent + patch_width + spacing, y), patch_width, patch_height,
                                        facecolor=orig_handle[1], edgecolor=orig_handle[1], transform=trans)
            return [patch1, patch2]
    legend_handles = []
    legend_labels = []

    for day_val in range(0,7):
        ax[0].plot(hour_range, day_of_week_allspecies_mean[days_of_week[day_val]], color=cmap1(color_values_week[day_val]), label = days_of_week[day_val])
        ax[0].set_xlabel('Hour')
        ax[0].set_ylabel('NOx Concentration (ppb)')
        ax[0].set_ylim([0,30])
        ax[0].set_title('USOS NOx ' + '(' + merged_or_usos + ')')

        ax[1].plot(hour_range, day_of_week_ozone_from[days_of_week[day_val]], color=cmap2(color_values_week[day_val]), label = days_of_week[day_val])
        ax[1].set_xlabel('Hour')
        ax[1].set_ylabel('Ozone Concentration (ppb)')
        ax[1].set_ylim([20,80])
        ax[1].set_title('USOS Ozone ' + '(' + merged_or_usos + ')' )
        
        legend_handles.append((cmap1(color_values_week[day_val]), cmap2(color_values_week[day_val])))
        legend_labels.append(days_of_week[day_val])

        # Add to legend
        ax[0].legend(legend_handles, legend_labels, handler_map={tuple: TwoColorHandler()})

    ax[0].plot(hour_range, overall_nox_avg, color = 'k', label = 'Campaign Average')
    ax[0].fill_between(hour_range, overall_nox_avg - overall_nox_std, overall_nox_avg + overall_nox_std, color='k',alpha = 0.1)
    ax[0].xaxis.set_minor_locator(MultipleLocator(3))
    ax[0].xaxis.set_major_locator(MultipleLocator(6))
    ax[0].grid(which='both', linestyle='--', linewidth=0.5)

    # box = ax[0].get_position()
    # ax[0].set_position([box.x0, box.y0, box.width * 0.8, box.height])
    # ax[0].legend(loc='center left', bbox_to_anchor=(1, 0.5))

    ax[1].plot(hour_range, o3_avg_from, color = 'k', label = 'Campaign Average')
    ax[1].fill_between(hour_range, o3_avg_from - o3_std_from, o3_avg_from + o3_std_from, color='k',alpha = 0.1)
    ax[1].xaxis.set_minor_locator(MultipleLocator(3))
    ax[1].xaxis.set_major_locator(MultipleLocator(6))
    ax[1].grid(which='both', linestyle='--', linewidth=0.5)
    ax[1].hlines(y=70, xmin = 0, xmax = 23, color = 'tab:red')
    plt.savefig(dirpath + '/Plotting/USOS_Campaign_analysis/plots/hawthorne_nox_ozone_day_of_week_comparison_' + merged_or_usos + '.png', dpi =300)
    plt.show()

def day_of_week_nox_ozone_no_avg_std(colormap1, colormap2, day_of_week_ozone_from, merged_or_usos):
    fig, ax = plt.subplots(1,2, figsize = (20,10))
    # Create 20 color values for gradient in the colormap
    color_values_week = np.linspace(0.2, 1, 8)

    cmap1 = mpl.colormaps[colormap1]
    cmap2 = mpl.colormaps[colormap2]

    class TwoColorHandler(HandlerBase):
        def create_artists(self, legend, orig_handle,
                        xdescent, ydescent, width, height, fontsize, trans):
            patch_width = width / 2.5  # width of each square
            patch_height = height * 0.6  # height of each square
            spacing = 4
            y = ydescent + (height - patch_height) / 2

            # Two colored rectangles side by side
            patch1 = mpatches.Rectangle((xdescent, y), patch_width, patch_height,
                                        facecolor=orig_handle[0], edgecolor=orig_handle[0], transform=trans)
            patch2 = mpatches.Rectangle((xdescent + patch_width + spacing, y), patch_width, patch_height,
                                        facecolor=orig_handle[1], edgecolor=orig_handle[1], transform=trans)
            return [patch1, patch2]
    legend_handles = []
    legend_labels = []

    for day_val in range(0,7):
        ax[0].plot(hour_range, day_of_week_allspecies_mean[days_of_week[day_val]], color=cmap1(color_values_week[day_val]), label = days_of_week[day_val])
        ax[0].set_xlabel('Hour')
        ax[0].set_ylabel('NOx Concentration (ppb)')
        ax[0].set_ylim([0,30])
        ax[0].set_title('USOS NOx ' + '(' + merged_or_usos + ')')
        
        ax[1].plot(hour_range, day_of_week_ozone_from[days_of_week[day_val]], color=cmap2(color_values_week[day_val]), label = days_of_week[day_val])
        ax[1].set_xlabel('Hour')
        ax[1].set_ylabel('Ozone Concentration (ppb)')
        ax[1].set_ylim([20,80])
        ax[1].set_title('USOS Ozone ' + '(' + merged_or_usos + ')' )

        legend_handles.append((cmap1(color_values_week[day_val]), cmap2(color_values_week[day_val])))
        legend_labels.append(days_of_week[day_val])

        # Add to legend
        ax[0].legend(legend_handles, legend_labels, handler_map={tuple: TwoColorHandler()}, loc='center left', bbox_to_anchor=(1, 0.5))

    ax[0].xaxis.set_minor_locator(MultipleLocator(3))
    ax[0].xaxis.set_major_locator(MultipleLocator(6))
    ax[0].grid(which='both', linestyle='--', linewidth=0.5)

    box = ax[0].get_position()
    ax[0].set_position([box.x0, box.y0, box.width * 0.8, box.height])


    ax[1].xaxis.set_minor_locator(MultipleLocator(3))
    ax[1].xaxis.set_major_locator(MultipleLocator(6))
    ax[1].grid(which='both', linestyle='--', linewidth=0.5)
    ax[1].hlines(y=70, xmin = 0, xmax = 23, color = 'tab:red')
    plt.savefig(dirpath + '/Plotting/USOS_Campaign_analysis/plots/hawthorne_nox_ozone_day_of_week_comparison_' + merged_or_usos + '_no_avg_std.png', dpi =300)
    plt.show()

def weekday_weekend_nox_ozone(colormap1, colormap2, weekday_ozone_from, merged_or_usos):
    fig, ax = plt.subplots(1,2, figsize = (20,10))
    # Create 20 color values for gradient in the colormap
    color_values_week = np.linspace(0.2, 1, 8)

    cmap1 = mpl.colormaps[colormap1]
    cmap2 = mpl.colormaps[colormap2]

    class TwoColorHandler(HandlerBase):
        def create_artists(self, legend, orig_handle,
                        xdescent, ydescent, width, height, fontsize, trans):
            patch_width = width / 2.5  # width of each square
            patch_height = height * 0.6  # height of each square
            spacing = 4
            y = ydescent + (height - patch_height) / 2

            # Two colored rectangles side by side
            patch1 = mpatches.Rectangle((xdescent, y), patch_width, patch_height,
                                        facecolor=orig_handle[0], edgecolor=orig_handle[0], transform=trans)
            patch2 = mpatches.Rectangle((xdescent + patch_width + spacing, y), patch_width, patch_height,
                                        facecolor=orig_handle[1], edgecolor=orig_handle[1], transform=trans)
            return [patch1, patch2]
    legend_handles = []
    legend_labels = []

    for day_val in range(0,2):
        ax[0].plot(hour_range, weekday_allspecies_mean[weekday_type_list[day_val]], color=cmap1(color_values_week[day_val]), label = days_of_week[day_val])
        ax[0].set_xlabel('Hour')
        ax[0].set_ylabel('NOx Concentration (ppb)')
        ax[0].set_ylim([0,30])
        ax[0].set_title('USOS NOx ' + '(' + merged_or_usos + ')')

        ax[1].plot(hour_range, weekday_ozone_from[weekday_type_list[day_val]], color=cmap2(color_values_week[day_val]), label = days_of_week[day_val])
        ax[1].set_xlabel('Hour')
        ax[1].set_ylabel('Ozone Concentration (ppb)')
        ax[1].set_ylim([20,80])
        ax[1].set_title('USOS Ozone ' + '(' + merged_or_usos + ')' )

        legend_handles.append((cmap1(color_values_week[day_val]), cmap2(color_values_week[day_val])))
        legend_labels.append(weekday_type_list[day_val])

        # Add to legend
        ax[0].legend(legend_handles, legend_labels, handler_map={tuple: TwoColorHandler()}, loc='center left', bbox_to_anchor=(1, 0.5))
    ax[0].xaxis.set_minor_locator(MultipleLocator(3))
    ax[0].xaxis.set_major_locator(MultipleLocator(6))
    ax[0].grid(which='both', linestyle='--', linewidth=0.5)

    ax[1].xaxis.set_minor_locator(MultipleLocator(3))
    ax[1].xaxis.set_major_locator(MultipleLocator(6))
    ax[1].grid(which='both', linestyle='--', linewidth=0.5)
    #ax[1].hlines(y=70, xmin = 0, xmax = 23, color = 'tab:red')

    plt.savefig(dirpath + '/Plotting/USOS_Campaign_analysis/plots/hawthorne_nox_ozone_weekday_weekend_comparison_' + merged_or_usos + '.png', dpi =300)
    plt.show()

def weekday_saturday_sunday_nox_ozone(colormap1, colormap2, weekday_ozone_from, day_of_week_ozone_from, merged_or_usos):
    fig, ax = plt.subplots(1,2, figsize = (20,10))
    color_values_week = np.linspace(0.2, 1.2, 3)
    cmap1 = mpl.colormaps[colormap1]
    cmap2 = mpl.colormaps[colormap2]

    class TwoColorHandler(HandlerBase):
        def create_artists(self, legend, orig_handle,
                        xdescent, ydescent, width, height, fontsize, trans):
            patch_width = width / 2.5  # width of each square
            patch_height = height * 0.6  # height of each square
            spacing = 4
            y = ydescent + (height - patch_height) / 2

            # Two colored rectangles side by side
            patch1 = mpatches.Rectangle((xdescent, y), patch_width, patch_height,
                                        facecolor=orig_handle[0], edgecolor=orig_handle[0], transform=trans)
            patch2 = mpatches.Rectangle((xdescent + patch_width + spacing, y), patch_width, patch_height,
                                        facecolor=orig_handle[1], edgecolor=orig_handle[1], transform=trans)

            return [patch1, patch2]
    legend_handles = []
    legend_labels = []

    # for day_val in range(0,2):
    ax[0].plot(hour_range, weekday_allspecies_mean['Weekday'], color=cmap1(color_values_week[0]), label = 'Weekday')
    ax[0].plot(hour_range, day_of_week_allspecies_mean['Saturday'], color=cmap1(color_values_week[1]), label = 'Saturday')
    ax[0].plot(hour_range, day_of_week_allspecies_mean['Sunday'], color=cmap1(color_values_week[2]), label = 'Sunday')
    ax[0].set_xlabel('Hour')
    ax[0].set_ylabel('NOx Concentration (ppb)')
    ax[0].set_xticks(np.arange(0,24))
    #ax[0].set_ylim(y_limit1_weekend)
    ax[0].set_title('USOS NOx ' + '(' + merged_or_usos + ')')

    ax[1].plot(hour_range, weekday_ozone_from['Weekday'], color=cmap2(color_values_week[0]), label = 'Weekday')
    ax[1].plot(hour_range, day_of_week_ozone_from['Saturday'], color=cmap2(color_values_week[1]), label = 'Saturday')
    ax[1].plot(hour_range, day_of_week_ozone_from['Sunday'], color=cmap2(color_values_week[2]), label = 'Sunday')
    ax[1].set_xlabel('Hour')
    ax[1].set_ylabel('Ozone Concentration (ppb)')
    #ax[1].set_ylim(y_limit2_weekend)
    ax[1].set_title('USOS Ozone ' + '(' + merged_or_usos + ')')

    legend_handles.append((cmap1(color_values_week[0]), cmap2(color_values_week[0])))
    legend_handles.append((cmap1(color_values_week[1]), cmap2(color_values_week[1])))
    legend_handles.append((cmap1(color_values_week[2]), cmap2(color_values_week[2])))
    legend_labels.append('Weekday')
    legend_labels.append('Saturday')
    legend_labels.append('Sunday')

    # Add to legend
    ax[0].legend(legend_handles, legend_labels, handler_map={tuple: TwoColorHandler()}, loc='center left', bbox_to_anchor=(1, 0.5))

            # Create two small colored markers for one label

    # ax[0].plot(hour_range, overall_nox_avg, color = 'k', label = 'Campaign Average')
    # ax[0].fill_between(hour_range, overall_nox_avg - overall_nox_std, overall_nox_avg + overall_nox_std, color='k',alpha = 0.1)
    ax[0].xaxis.set_minor_locator(MultipleLocator(3))
    ax[0].xaxis.set_major_locator(MultipleLocator(6))
    ax[0].set_grid(which='both', linestyle='--', linewidth=0.5)

    # box = ax[0].get_position()
    # ax[0].set_position([box.x0, box.y0, box.width * 0.8, box.height])
    #ax[0].legend(loc='center left', bbox_to_anchor=(1, 0.5))

    # ax[1].plot(hour_range, overall_ozone_avg, color = 'k', label = 'Campaign Average')
    # ax[1].fill_between(hour_range, overall_ozone_avg_merged - overall_ozone_std_merged, overall_ozone_avg_merged + overall_ozone_std_merged, color='k',alpha = 0.1)
    ax[1].xaxis.set_minor_locator(MultipleLocator(3))
    ax[1].xaxis.set_major_locator(MultipleLocator(6))
    ax[1].grid(which='both', linestyle='--', linewidth=0.5)
    #ax[1].hlines(y=70, xmin = 0, xmax = 23, color = 'tab:red')

    plt.savefig(dirpath + '/Plotting/USOS_Campaign_analysis/plots/hawthorne_nox_ozone_weekday_saturday_sunday_comparison_' + merged_or_usos + '.png', dpi =300)
    plt.show()

def exceedance_count():
    count_exceedances_merged = df_ozone[df_ozone['Exceedance_day_merged_ozone'] == True].shape[0]
    print('Exceedances (merged data): ', int(count_exceedances_merged/24))
    count_weekday_exceedances_merged = df_ozone[(df_ozone['Weekday'] == 'Weekday') & (df_ozone['Exceedance_day_merged_ozone'] == True)].shape[0]
    print('Exceedances on weekdays (merged data): ', int(count_weekday_exceedances_merged/24))

    count_exceedances_usos = df_ozone[df_ozone['Exceedance_day_Mobile Lab O3'] == True].shape[0]
    print('Exceedances (Mobile Lab only): ', int(count_exceedances_usos/24))
    count_weekday_exceedances_usos = df_ozone[(df_ozone['Weekday'] == 'Weekday') & (df_ozone['Exceedance_day_Mobile Lab O3'] == True)].shape[0]
    print('Exceedances on weekdays (Mobile Lab only): ', int(count_weekday_exceedances_usos/24))

    #Group by each date and check which ones have all values in the Exceedance column as True
    group_separation_check_true_merged = df_ozone['Exceedance_day_merged_ozone'].groupby(df_ozone.index.normalize()).all()

    #Filter which dates have all values in the Exceedance column as True
    true_days_merged = group_separation_check_true_merged[group_separation_check_true_merged]

    weekdays_merged = true_days_merged.index.day_name() #Get days of week for each exceedance date
    day_of_week_counts_merged = weekdays_merged.value_counts().sort_index() #count how many days of the week per day
    print('Exceedance by day of week count (merged):', day_of_week_counts_merged)

    #Group by each date and check which ones have all values in the Exceedance column as True
    group_separation_check_true_usos = df_ozone['Exceedance_day_Mobile Lab O3'].groupby(df_ozone.index.normalize()).all()

    #Filter which dates have all values in the Exceedance column as True
    true_days_usos = group_separation_check_true_usos[group_separation_check_true_usos]

    weekdays_usos = true_days_usos.index.day_name() #Get days of week for each exceedance date
    day_of_week_counts_usos = weekdays_usos.value_counts().sort_index() #count how many days of the week per day
    print('Exceedance by day of week count (Mobile Lab only):', day_of_week_counts_usos)

##### Call Functions #####
# ozone_daily(
#     colormap = 'viridis',
#     o3_measured_from = 'merged_ozone',
#     o3_avg_from = overall_ozone_avg_merged,
#     o3_std_from = overall_ozone_std_merged,
#     merged_or_usos = 'merged'
# )

# ozone_daily(
#     colormap = 'viridis',
#     o3_measured_from = 'USOS O3',
#     o3_avg_from = overall_ozone_avg_usos,
#     o3_std_from = overall_ozone_std_usos,
#     merged_or_usos = 'usos'
# )

# nox_ozone_daily(
#     colormap1 = 'YlGnBu_r',
#     colormap2 = 'YlOrRd_r',
#     o3_measured_from = 'merged_ozone',
#     o3_avg_from = overall_ozone_avg_merged,
#     o3_std_from = overall_ozone_std_merged,
#     merged_or_usos = 'merged'    
# )
# nox_ozone_daily(
#     colormap1 = 'YlGnBu_r',
#     colormap2 = 'YlOrRd_r',
#     o3_measured_from = 'USOS O3',
#     o3_avg_from = overall_ozone_avg_usos,
#     o3_std_from = overall_ozone_std_usos,
#     merged_or_usos = 'usos'   
# )

# day_of_week_nox_ozone(
#     colormap1 = 'YlGnBu',
#     colormap2 = 'YlOrRd',
#     day_of_week_ozone_from = day_of_week_ozone_mean_merged,
#     o3_avg_from = overall_ozone_avg_merged,
#     o3_std_from = overall_ozone_std_merged,
#     merged_or_usos = 'merged'
# )

# day_of_week_nox_ozone(
#     colormap1 = 'YlGnBu',
#     colormap2 = 'YlOrRd',
#     day_of_week_ozone_from = day_of_week_ozone_mean_usos,
#     o3_avg_from = overall_ozone_avg_usos,
#     o3_std_from = overall_ozone_std_usos,
#     merged_or_usos = 'usos'
# )

# day_of_week_nox_ozone_no_avg_std(
#     colormap1 = 'YlGnBu',
#     colormap2 = 'YlOrRd',
#     day_of_week_ozone_from = day_of_week_ozone_mean_merged,
#     merged_or_usos = 'merged'
# )

# day_of_week_nox_ozone_no_avg_std(
#     colormap1 = 'YlGnBu',
#     colormap2 = 'YlOrRd',
#     day_of_week_ozone_from = day_of_week_ozone_mean_usos,
#     merged_or_usos = 'usos'
# )

# weekday_weekend_nox_ozone(
#     colormap1 = 'YlGnBu',
#     colormap2 = 'YlOrRd',
#     weekday_ozone_from = weekday_ozone_mean_merged,
#     merged_or_usos = 'merged'
# )

# weekday_weekend_nox_ozone(
#     colormap1 = 'YlGnBu',
#     colormap2 = 'YlOrRd',
#     weekday_ozone_from = weekday_ozone_mean_usos,
#     merged_or_usos = 'usos'
# )
# weekday_saturday_sunday_nox_ozone(
#     colormap1 = 'YlGnBu',
#     colormap2 = 'YlOrRd',
#     weekday_ozone_from = weekday_ozone_mean_merged,
#     day_of_week_ozone_from = day_of_week_ozone_mean_merged,
#     merged_or_usos = 'merged'
# )

# weekday_saturday_sunday_nox_ozone(
#     colormap1 = 'YlGnBu',
#     colormap2 = 'YlOrRd',
#     weekday_ozone_from = weekday_ozone_mean_usos,
#     day_of_week_ozone_from = day_of_week_ozone_mean_usos,
#     merged_or_usos = 'usos'
# )

exceedance_count()

leighton_ratio_daily()
leighton_ratio_daily_exceedances_only()
leighton_ratio_daily_exceedance_comparison()