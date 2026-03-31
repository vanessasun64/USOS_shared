#Working on adding columns to calibration_adjustment_odr function, under #Evaluation of adjustment validity
#to evaluate improvement of RMSE, slope

import pandas as pd
import xarray as xr
import numpy as np

import pytz

import statsmodels.api as sm
from scipy import stats
from scipy.stats import pearsonr, spearmanr
# from scipy import odr
from odrpack import odr_fit

from sklearn.model_selection import KFold
from sklearn.metrics import mean_absolute_error

import re
from collections import defaultdict

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator
#import cmasher as cmr
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
from datetime import datetime

import os


dirpath = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/'
hawthorne_data_dir = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/Hawthorne_data/data/'

#region: Plot formatting
mpl.rcParams['xtick.labelsize'] = 15
mpl.rcParams['ytick.labelsize'] = 15
mpl.rcParams['legend.fontsize'] = 16
mpl.rcParams['axes.labelsize'] = 18
mpl.rcParams['axes.titlesize'] = 22
mpl.rcParams['axes.xmargin'] = 0

# Set the font family to 'serif'
mpl.rcParams['font.family'] = 'serif'
# Specify preferred serif font (Computer Modern Roman is 'cmr10')
mpl.rcParams['font.serif'] = 'Lato' 
# Optionally, configure mathtext to use Computer Modern fonts as well
mpl.rcParams['mathtext.fontset'] = 'cm'
# Ensure minus signs are rendered correctly with CM fonts
mpl.rcParams['axes.unicode_minus'] = False

#endregion

#region: old defunct functions
# def plot_both(noaa, udaq, noaa_species, udaq_species, date):
#     start = pd.to_datetime(date)
#     stop = start + pd.Timedelta(days=1)

#     noaa_sel = noaa[udaq.index.to_series().between(start, stop, inclusive='left')]
#     udaq_sel = udaq[udaq.index.to_series().between(start, stop, inclusive='left')]

#     noaa_name= noaa[noaa_species]
#     udaq_name= udaq[udaq_species]

#     #Local time as index
#     noaa_var=noaa_sel[noaa_name] 
#     udaq_var=udaq_sel[udaq_name]

#     mx=np.nanmax([np.nanmax(noaa_var), np.nanmax(udaq_var)])

#     fig, ax1=plt.subplots(nrows=1, ncols=1)
#     ax1.scatter(noaa_sel.index, noaa_var, color = 'r', label='NOAA')
#     ax1.scatter(udaq_sel.index, udaq_var, color = 'b', label='UDAQ')
#     fig.suptitle(udaq_species)
#     ax1.set_ylim([0,1.05*mx])
#     ax1.legend()
#     plt.tight_layout()
#     plt.show()
#     return

# def which_is_better(file, date, VOC, preferred):

#     line=f"{date}, {VOC}, {preferred}\n"

#     with open(file, 'a') as f :
#         f.write(line)
#     return
#endregion

#Load UDAQ and NOAA ML VOC files respectively
#Time index is every 15 min from 07/14/2024 18:00:00 to 8/18/2024 17:45:00 Local Time
noaa_f= dirpath + '/CampaignData_and_Merges/R0/CSL_MobileLab_Parked/merged/rev_15min/all_CSL_MobileLab_Parked_rev15min_iWASupdated.nc'

ds = xr.open_dataset(noaa_f)
noaa = ds.to_dataframe()
noaa = noaa.set_index(['time_local'])
#localize time zone to MDT
noaa.index = noaa.index.tz_localize('America/Denver')
#Set index to only span from 2024-07-15 00:00:00 to 2024-08-18 17:45:00
noaa = noaa.loc['2024-07-15 00:00:00':'2024-08-18 17:45:00']

display(noaa['Styrene_PTR'])

# Set any inf, neg. inf, and negative values to NaN
for col in noaa.columns:
    noaa[col] = noaa[col].replace([np.inf, -np.inf], np.nan).mask(noaa[col] < 0, np.nan)
#Give explicit variable to CO that will be used in functions. Turn any negatives into NaNs.
ml_co_raw = noaa['CO_Piccaro']
#There seems to be a data issue with the CO Picarro during 07/22 15:00 MDT to 07/24 00:00 MDT so mask that range into NaNs
ml_co_raw.loc['2024-07-22 11:00:00':'2024-07-22 18:00:00'] = np.nan
ml_co_raw.loc['2024-07-23 10:30:00':'2024-07-23 12:00:00'] = np.nan

udaq_f = dirpath + '/Hawthorne_data/data/script_output/hawthorne_udaq_all_vocs_15min_timezone_carbon_number_updated.csv'
udaq = pd.read_csv(udaq_f, index_col='time_local', parse_dates=True)
#For some reason, pandas is reading the time_local as UTC so we change it back to reading as the UTC-6 time zone
udaq.index = udaq.index.tz_localize(None)
udaq.index = udaq.index.tz_localize('America/Denver')
udaq = udaq.loc['2024-07-15 00:00:00':'2024-08-18 17:45:00']

# Set any inf, neg. inf, and negative values to NaN
for col in udaq.columns[1:len(udaq.columns)+1]:
    udaq[col] = udaq[col].replace([np.inf, -np.inf], np.nan).mask(udaq[col] < 0, np.nan)

#get mappings to match ML names to UDAQ names
mappings_filepath = dirpath + '/Hawthorne_data/mappings/manually_edited/UDAQ_Hawthorne_CRACMM_GEOSCHEM_CB6r5h_mapped_updated_01292026.csv'
df_mapping_parameters = pd.read_csv(mappings_filepath)
df_mapping_parameters = df_mapping_parameters.drop([0]) #drop Total NMVOCs

mapping = {
    (c1 if pd.notna(c1) else f'ML_NoVar_{i}'):(c2 if pd.notna(c2) else f'UDAQ_NoVar_{i}')
    for i, (c1, c2) in enumerate(zip(df_mapping_parameters['USOS Mapping'], df_mapping_parameters['UDAQ_Variable']))
}

#Sum mEthylToluene and pEthylToluene and change mapping so that 'x3_x4_EthylToluene_WAS' maps to their sum
udaq['mpEthyltoluene'] = udaq['mEthyltoluene'] + udaq['pEthyltoluene']
udaq = udaq.drop(['mEthyltoluene','pEthyltoluene'], axis = 1)
mapping['x3_x4_EthylToluene_WAS'] = 'mpEthyltoluene'
print(mapping)

#List of all the WAS species that aren't in the mapping (We can't use any UDAQ data for filling holes but we'll need to interpolate the gaps values in order to use in F0AM)
#Takes all the WAS species in the original Mobile Lab data and removes all those in the mapping, so that we know which ones we won't be cross-calibrating UDAQ data with
#We can still see how well the CO tracer works to fill in holes
was_species_list = []
for colname in noaa.columns:
    if 'WAS' in colname:
        was_species_list.append(colname)
    else:
        pass
was_species_not_in_mapping = list(set(was_species_list) - set(mapping.keys()))
print(was_species_not_in_mapping)

#Create log file
logfile='tracerlog.txt'

hour_range = np.arange(0,24,1)

def plot_both(noaa, udaq, noaa_species, udaq_species, dates):
    '''
    This function plots the concentration of a VOC species for one day, for the NOAA Mobile Lab and UDAQ data
    '''

    only_date = dates.date()
    start_time = pd.to_datetime(dates)
    stop_time = start_time + pd.Timedelta(days=1)

    noaa_sel = noaa[noaa.index.to_series().between(start_time, stop_time, inclusive='left')]
    udaq_sel = udaq[udaq.index.to_series().between(start_time, stop_time, inclusive='left')]
        
    if 'ML_NoVar' in noaa_species:
        nan_series = pd.Series(index=noaa_sel.index)
        noaa_var = nan_series
        udaq_var=udaq_sel[udaq_species]
    elif 'ML_NoVar' in udaq_species:
        noaa_var=noaa_sel[noaa_species] 
        nan_series = pd.Series(index=noaa_sel.index)
        udaq_var = nan_series
    else:
        noaa_var=noaa_sel[noaa_species]
        udaq_var=udaq_sel[udaq_species]

    # mx=np.nanmax([np.nanmax(noaa_var), np.nanmax(udaq_var)])
    fig, ax1=plt.subplots(nrows=1, ncols=1)
    ax1.scatter(noaa_sel.index, noaa_var, color = 'r', label='NOAA')
    ax1.scatter(udaq_sel.index, udaq_var, color = 'b', label='UDAQ')
    ax1.xaxis.set_major_locator(mdates.DayLocator())
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%H'))
    #Major ticks: every 6 hours
    ax1.xaxis.set_major_locator(mdates.HourLocator(byhour=[0, 6, 12, 18]))
    # Minor ticks: every 3 hours
    ax1.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21]))
    # Rotate and format tick labels
    ax1.tick_params(axis='x', which='major')
    ax1.tick_params(axis='x', which='minor', length=3, color='gray')
    fig.suptitle(f'{udaq_species}, {only_date}')
    ax1.set_xlabel('Time (MDT)')
    ax1.set_ylabel(udaq_species)
    # ax1.set_ylim([0,1.05*mx])
    ax1.set_xlim([start_time, stop_time])
    ax1.legend()
    plt.tight_layout()
    plt.show()
    return

def which_is_better(file, dates, VOC, preferred):
	line=f"{date}, {VOC}, {preferred}\n"

	with open(file, 'a') as f :
			f.write(line)
	return
		
def gap_measuring_daily(noaa_species, dates):
    '''
    This function measures the gaps between measurements for one day for the NOAA Mobile Lab and UDAQ data.
    Output:
        gap_df: a dataframe of all the gaps for the day of data
    '''
    start = pd.to_datetime(dates)
    stop = start + pd.Timedelta(days=1)
    noaa_sel = noaa[noaa.index.to_series().between(start, stop, inclusive='left')]

    if 'ML_NoVar' in noaa_species:
        nan_series = pd.Series(index=noaa_sel.index)
        noaa_var = nan_series
    else:
        #Identify where real measurements exist
        is_measured = noaa_sel[noaa_species].notna()
        #print('is_measured:', is_measured)

        #Compute time gaps only between valid measurements
        measurement_times = noaa_sel.index[is_measured]
        gaps = pd.Series(measurement_times, index=measurement_times).diff()
        #print('gaps:', gaps)

        #Store dataframe of gaps
        gap_df = pd.DataFrame({
            'gap_start': pd.Series(measurement_times, index = measurement_times).shift(1),
            'gap_end': measurement_times,
            'gap_size': gaps
        }).reset_index(drop=True)
        #Create column that gives boolean on if the gap is a large gap (> 3 hours)
        gap_df['is_large_gap'] = gap_df['gap_size'] > pd.Timedelta(hours=3)
        #print('gap_df:', gap_df)
        return gap_df
def gap_measuring_total(noaa_species):
    gap_records_fulltime = []

    if 'ML_NoVar' in noaa_species:
        nan_series = pd.Series(index=noaa.index)
        noaa_var = nan_series
    else:
        #Identify where real measurements exist
        is_measured_fulltime = noaa[noaa_species].notna()
        #print('is_measured:', is_measured)

        #Compute time gaps only between valid measurements
        measurement_times_fulltime = noaa.index[is_measured_fulltime]
        #gaps_fulltime =  pd.Series(measurement_times_fulltime) - pd.Series(measurement_times_fulltime).shift(1)
        #print('gaps:', gaps)

        # Internal gaps
        if len(measurement_times_fulltime) > 1:
            for start_of_gap, end_of_gap in zip(measurement_times_fulltime[:-1], measurement_times_fulltime[1:]):
                gap_size = end_of_gap - start_of_gap
                gap_records_fulltime.append({
                    'gap_start': start_of_gap,
                    'gap_end': end_of_gap,
                    'gap_size': gap_size,
                    'is_large_gap': gap_size >  pd.Timedelta(hours=3)
                })

        # Leading gap
        if len(measurement_times_fulltime) >= 1:
            leading_gap = measurement_times_fulltime[0] - noaa.index.min()
            if leading_gap > pd.Timedelta(0):
                gap_records_fulltime.insert(0, {
                    'gap_start': noaa.index.min(),
                    'gap_end': measurement_times_fulltime[0],
                    'gap_size': leading_gap,
                    'is_large_gap': leading_gap >  pd.Timedelta(hours=3)
                })

        # Trailing gap
        if len(measurement_times_fulltime) >= 1:
            trailing_gap = noaa.index.max() - measurement_times_fulltime[-1]
            if trailing_gap > pd.Timedelta(0):
                gap_records_fulltime.append({
                    'gap_start': measurement_times_fulltime.max(),
                    'gap_end': noaa.index.max(),
                    'gap_size': trailing_gap,
                    'is_large_gap': trailing_gap >  pd.Timedelta(hours=3)
                })

        # #Store dataframe of gaps
        # gap_df_fulltime = pd.DataFrame({
        #     'gap_start': pd.Series(measurement_times_fulltime, index = measurement_times_fulltime).shift(1),
        #     'gap_end': measurement_times_fulltime,
        #     'gap_size': gaps_fulltime
        # }).reset_index(drop=True)
        #Create column that gives boolean on if the gap is a large gap (> 3 hours)
        gap_df_fulltime = pd.DataFrame(gap_records_fulltime)
        display(gap_df_fulltime)
        # gap_df_fulltime['is_large_gap'] = gap_records_fulltime['gap_size'] > pd.Timedelta(hours=3)
        # print('gap_df_fulltime: \n', gap_df_fulltime)
        # return gap_df_fulltime

def formaldehyde_only_calibration_adjustment_ols(time_interval_formaldehyde):
    # read UDAQ Formaldehyde data, provided by and Bart (UDAQ), which is in UTC
    udaq_formaldehyde_load = hawthorne_data_dir + 'hw_zero_corrected_data_formaldehyde.csv'
    df_udaq_formaldehyde = pd.read_csv(udaq_formaldehyde_load, index_col='dt', parse_dates=True)
    df_udaq_formaldehyde.index = df_udaq_formaldehyde.index.tz_localize('UTC')

    df_udaq_formaldehyde_usos_only = df_udaq_formaldehyde.sort_index().loc['2024-07-15 00:00:00':'2024-08-18 23:59:00']
    df_udaq_formaldehyde_usos_only.index = df_udaq_formaldehyde_usos_only.index.rename('time_UTC')
    keep_colnames = ['H2CO_Corrected']
    df_udaq_formaldehyde_usos_only = df_udaq_formaldehyde_usos_only[keep_colnames]

    #This averaging in this section is based off the same code used for averaging the icartt data
    #as in the function align2master_timeline in time_utils.py of the icartt_read_and_merge

    # Get the average native sampling frequency in total seconds:
    tseries = df_udaq_formaldehyde_usos_only.index.to_series()
    #Avg native sampling frequency
    min_sep = int(np.round(tseries.diff().median().total_seconds()))
    #Intended interval
    step_S = time_interval_formaldehyde

    #Reindex to avg native sampling frequency (in this case, 1 minute)
    new_start_time = pd.Timestamp('2024-07-15 00:00:00').tz_localize('UTC')
    new_end_time = pd.Timestamp('2024-08-18 23:59:00').tz_localize('UTC')
    dts = pd.date_range(new_start_time, new_end_time, freq=str(min_sep) + 's')
    dfn = df_udaq_formaldehyde_usos_only.reindex(dts, method='nearest', fill_value=np.nan)

    # Take a centered boxcar average around the 900s (15 min) avg. (for numerical columns only)
    #NOTE: .mean() handles Nans like np.nanmean() in this context!!! 
    df_nums=dfn.select_dtypes(exclude=['datetime64'])
    df_nums_new = df_nums.rolling(str(int(step_S)) + 's').mean().resample(str(step_S) + 's').mean()

    df_avg_formaldehyde = df_nums_new
    #Convert to local time (MDT)
    df_avg_formaldehyde.index = df_avg_formaldehyde.index.tz_convert('America/Denver')
    df_avg_formaldehyde.index = df_avg_formaldehyde.index.rename('time_local')
    #print(df_new_formaldehyde)
    df_avg_formaldehyde = df_avg_formaldehyde.rename(columns = {'H2CO_Corrected': 'H2CO_UDAQ_Corrected'})
    df_avg_formaldehyde = df_avg_formaldehyde.loc['2024-07-15 00:00:00':'2024-08-18 17:45:00']
    
    df_avg_formaldehyde['ML Data']= noaa['HCHO_CRDS']

    # Set any inf, neg. inf, and negative values to NaN
    # Take the overlap in the ML and UDAQ data only
    df_drop_formaldehyde = df_avg_formaldehyde.replace([np.inf, -np.inf], np.nan).mask(udaq[col] < 0, np.nan)
    df_clean_formaldehyde = df_drop_formaldehyde.dropna()

    # Add constant term for intercept
    M_with_const = sm.add_constant(df_clean_formaldehyde['H2CO_UDAQ_Corrected'])

    # Fit OLS regression: O = a + b*M
    ols_calibration = sm.OLS(df_clean_formaldehyde['ML Data'], M_with_const).fit()

    # Get calibration coefficients
    a = ols_calibration.params['const']  # additive offset
    b = ols_calibration.params[df_clean_formaldehyde['H2CO_UDAQ_Corrected'].name]   # multiplicative factor

    calibration_eq_formaldehyde = f"Calibration equation: O_hat = {a:.4f} + {b:.4f} * M"
    print(calibration_eq_formaldehyde)
    #print(ols_calibration.summary())
    udaq_formaldehyde_adjusted = a+b*df_drop_formaldehyde['H2CO_UDAQ_Corrected']

    #To determine if calibration adjustment aligns closer with ML data, use Root Mean Squared Error
    # Original RMSE
    rmse_original_formaldehyde = np.sqrt(np.mean((df_clean_formaldehyde['H2CO_UDAQ_Corrected'] - df_clean_formaldehyde['ML Data'])**2))

    # Calibrated RMSE
    #Get series of overlap between odr-fitted udaq data with the initial clean data (no NaNs between ML Data and UDAQ Data)
    udaq_overlap_adjusted_with_initial = udaq_formaldehyde_adjusted.loc[udaq_formaldehyde_adjusted.index.intersection(df_clean_formaldehyde.index)]
    rmse_calibrated_formaldehyde = np.sqrt(np.mean((udaq_overlap_adjusted_with_initial - df_clean_formaldehyde['ML Data'])**2))

    print(f"Original RMSE: {rmse_original_formaldehyde:.4f}")
    print(f"Calibrated RMSE: {rmse_calibrated_formaldehyde:.4f}")

    #Calculate percent improvement
    rmse_percent_improvement_formaldehyde = 100*((rmse_original_formaldehyde - rmse_calibrated_formaldehyde) / (rmse_original_formaldehyde))
    print(f'Percent improvement of calibration adjusted RMSE: {rmse_percent_improvement_formaldehyde:.4f}')

    #Merge the old formaldehyde ML data with the UDAQ data with calibration adjustment applied to fill holes
    merged_formaldehyde = noaa['HCHO_CRDS'].fillna(udaq_formaldehyde_adjusted)

    #change variable names to label as Ordinary Least Squares applied
    ml_initial_formaldehyde_ols = df_drop_formaldehyde['ML Data']
    ml_initial_formaldehyde_overlap_ols = df_clean_formaldehyde['ML Data']

    udaq_initial_formaldehyde_ols = df_drop_formaldehyde['H2CO_UDAQ_Corrected']
    udaq_initial_formaldehyde_overlap_ols = df_clean_formaldehyde['H2CO_UDAQ_Corrected']

    udaq_adjusted_formaldehyde_ols = udaq_formaldehyde_adjusted
    udaq_overlap_adjusted_formaldehyde_ols = udaq_overlap_adjusted_with_initial

    calibration_eq_formaldehyde_ols = calibration_eq_formaldehyde
    rmse_original_formaldehyde_ols = rmse_original_formaldehyde
    rmse_calibrated_formaldehyde_ols = rmse_calibrated_formaldehyde
    rmse_percent_improvement_formaldehyde_ols = rmse_percent_improvement_formaldehyde

    merged_formaldehyde_ols = merged_formaldehyde

    return (ml_initial_formaldehyde_ols, ml_initial_formaldehyde_overlap_ols, udaq_initial_formaldehyde_ols, 
    udaq_initial_formaldehyde_overlap_ols, udaq_adjusted_formaldehyde_ols, udaq_overlap_adjusted_formaldehyde_ols,
    calibration_eq_formaldehyde_ols, rmse_original_formaldehyde_ols, 
    rmse_calibrated_formaldehyde_ols, rmse_percent_improvement_formaldehyde_ols, merged_formaldehyde_ols)
def plot_formaldehyde_only_calibration_adjustment_ols(ml_initial_formaldehyde_ols, ml_initial_formaldehyde_overlap_ols, udaq_initial_formaldehyde_ols, 
    udaq_initial_formaldehyde_overlap_ols, udaq_adjusted_formaldehyde_ols, udaq_overlap_adjusted_formaldehyde_ols,
    calibration_eq_formaldehyde_ols, merged_formaldehyde_ols):

    fig, (ax1, ax2) = plt.subplots(2,1, figsize = (16,8), tight_layout=True)
    xlim_start_jul = pd.to_datetime('2024-07-15 00:00:00').tz_localize('America/Denver')
    xlim_end_jul = pd.to_datetime('2024-07-31 23:00:00').tz_localize('America/Denver')
    xlim_start_aug = pd.to_datetime('2024-08-01 00:00:00').tz_localize('America/Denver')
    xlim_end_aug = pd.to_datetime('2024-08-18 23:45:00').tz_localize('America/Denver')

    #ax1 is the first row of subplot, for July only
    ax1.plot(ml_initial_formaldehyde_overlap_ols.index, ml_initial_formaldehyde_overlap_ols, linestyle = 'solid', color='m', marker='x',label='Mobile Lab', alpha = 0.7)
    ax1.plot(udaq_initial_formaldehyde_overlap_ols.index, udaq_initial_formaldehyde_overlap_ols,linestyle = 'solid', color = 'b', marker = '+', label = 'Initial UDAQ')
    ax1.plot(udaq_overlap_adjusted_formaldehyde_ols.index, udaq_overlap_adjusted_formaldehyde_ols, linestyle='solid', color = 'y', marker = '.', label = 'Calibrated UDAQ')

    #Set x ticks
    tz_mdt = noaa.index.tz #this time zone should be in Mountain Daylight Time
    ax1.xaxis.set_major_locator(mdates.DayLocator(tz=tz_mdt))
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d', tz=tz_mdt))
    # Minor ticks: every 3 hours
    ax1.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21]))
    # Rotate and format tick labels
    ax1.tick_params(axis='x', which='major')
    ax1.tick_params(axis='x', which='minor', length=3, color='gray')
    #ax.grid(True, which='both')
    
    ax1.set_ylabel('Formaldehyde (ppb)')
    #ax1.set_xlabel('Date')
    ax1.margins(x=0)
    ax1.set_xlim([xlim_start_jul, xlim_end_jul])

    ax1.legend(loc = 'upper right')
    
    #ax2 is the second row of subplot, for August only
    ax2.plot(ml_initial_formaldehyde_overlap_ols.index, ml_initial_formaldehyde_overlap_ols, linestyle = 'solid', color='m', marker='x',label='Mobile Lab', alpha = 0.7)
    ax2.plot(udaq_initial_formaldehyde_overlap_ols.index, udaq_initial_formaldehyde_overlap_ols,linestyle = 'solid', color = 'b', marker = '+', label = 'Initial UDAQ')
    ax2.plot(udaq_overlap_adjusted_formaldehyde_ols.index, udaq_overlap_adjusted_formaldehyde_ols, linestyle='solid', color = 'y', marker = '.', label = 'Calibrated UDAQ')

    #Set x ticks
    ax2.xaxis.set_major_locator(mdates.DayLocator(tz=tz_mdt))
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d', tz=tz_mdt))
    # Minor ticks: every 3 hours
    ax2.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21]))
    # Rotate and format tick labels
    ax2.tick_params(axis='x', which='major')
    ax2.tick_params(axis='x', which='minor', length=3, color='gray')

    ax2.set_ylabel('Formaldehyde (ppb)')
    ax2.set_xlabel('Time (MDT)')
    ax2.margins(x=0)

    ax2.set_xlim([xlim_start_aug, xlim_end_aug])
    ax2.legend(loc = 'upper right')

    #Mark midnight for every day
    midnight_vals = []
    for midnight_idx in range(0,len(noaa.index),96):
        midnight_vals.append(noaa.index[midnight_idx])
    for day_pos in midnight_vals:
        ax1.axvline(day_pos, color = 'black', linestyle = 'dotted')
        ax2.axvline(day_pos, color = 'black', linestyle = 'dotted')

    plt.suptitle(calibration_eq_formaldehyde_ols)
    plt.savefig(dirpath + '/Merge_scripts/calibration_adjustments/calibration_adjustment_plots/ols_fitting/overlap/ols_cal_adj_hawthorne_udaq_ml_formaldehyde_overlap_comparison_july_aug_timeseries.png', dpi =300)
    plt.show()

    #PLOT 2
    fig, (ax1, ax2) = plt.subplots(2,1, figsize = (16,8), tight_layout=True)
    xlim_start_jul = pd.to_datetime('2024-07-15 00:00:00').tz_localize('America/Denver')
    xlim_end_jul = pd.to_datetime('2024-07-31 23:00:00').tz_localize('America/Denver')
    xlim_start_aug = pd.to_datetime('2024-08-01 00:00:00').tz_localize('America/Denver')
    xlim_end_aug = pd.to_datetime('2024-08-18 23:45:00').tz_localize('America/Denver')

    #ax1 is the first row of subplot, for July only
    ax1.plot(ml_initial_formaldehyde_overlap_ols.dropna().index, ml_initial_formaldehyde_overlap_ols.dropna(), linestyle = 'solid', color='m', marker='x',label='Mobile Lab', alpha = 0.7)
    ax1.plot(udaq_initial_formaldehyde_overlap_ols.dropna().index, udaq_initial_formaldehyde_overlap_ols.dropna(),linestyle = 'solid', color = 'b', marker = '+', label = 'Initial UDAQ')
    ax1.plot(udaq_overlap_adjusted_formaldehyde_ols.dropna().index, udaq_overlap_adjusted_formaldehyde_ols.dropna(), linestyle='solid', color = 'y', marker = '.', label = 'Calibrated UDAQ')

    #Set x ticks
    tz_mdt = noaa.index.tz #this time zone should be in Mountain Daylight Time
    ax1.xaxis.set_major_locator(mdates.DayLocator(tz=tz_mdt))
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d', tz=tz_mdt))
    # Minor ticks: every 3 hours
    ax1.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21]))
    # Rotate and format tick labels
    ax1.tick_params(axis='x', which='major')
    ax1.tick_params(axis='x', which='minor', length=3, color='gray')
    #ax.grid(True, which='both')
    
    ax1.set_ylabel('Formaldehyde (ppb)')
    #ax1.set_xlabel('Date')
    ax1.margins(x=0)
    ax1.set_xlim([xlim_start_jul, xlim_end_jul])

    ax1.legend(loc = 'upper right')
    
    #ax2 is the second row of subplot, for August only
    ax2.plot(ml_initial_formaldehyde_ols.dropna().index, ml_initial_formaldehyde_ols.dropna(), linestyle = 'solid', color='m', marker='x',label='Mobile Lab', alpha = 0.7)
    ax2.plot(udaq_initial_formaldehyde_ols.dropna().index, udaq_initial_formaldehyde_ols.dropna(),linestyle = 'solid', color = 'b', marker = '+', label = 'Initial UDAQ')
    ax2.plot(udaq_adjusted_formaldehyde_ols.dropna().index, udaq_adjusted_formaldehyde_ols.dropna(), linestyle='solid', color = 'y', marker = '.', label = 'Calibrated UDAQ')

    #Set x ticks
    ax2.xaxis.set_major_locator(mdates.DayLocator(tz=tz_mdt))
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d', tz=tz_mdt))
    # Minor ticks: every 3 hours
    ax2.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21]))
    # Rotate and format tick labels
    ax2.tick_params(axis='x', which='major')
    ax2.tick_params(axis='x', which='minor', length=3, color='gray')

    ax2.set_ylabel('Formaldehyde (ppb)')
    ax2.set_xlabel('Time (MDT)')
    ax2.margins(x=0)

    ax2.set_xlim([xlim_start_aug, xlim_end_aug])
    ax2.legend(loc = 'upper right')

    #Mark midnight for every day
    midnight_vals = []
    for midnight_idx in range(0,len(noaa.index),96):
        midnight_vals.append(noaa.index[midnight_idx])
    for day_pos in midnight_vals:
        ax1.axvline(day_pos, color = 'black', linestyle = 'dotted')
        ax2.axvline(day_pos, color = 'black', linestyle = 'dotted')

    plt.suptitle(calibration_eq_formaldehyde_ols)
    plt.savefig(dirpath + '/Merge_scripts/calibration_adjustments/calibration_adjustment_plots/ols_fitting/full_time_series/ols_cal_adj_hawthorne_udaq_ml_Formaldehyde_overlap_comparison_july_aug_timeseries.png', dpi =300)
    plt.show()

def formaldehyde_only_calibration_adjustment_odr(time_interval_formaldehyde):
    # read UDAQ Formaldehyde data, provided by and Bart (UDAQ), which is in UTC
    udaq_formaldehyde_load = hawthorne_data_dir + '/from_Bart/hw_zero_corrected_data_formaldehyde.csv'
    df_udaq_formaldehyde = pd.read_csv(udaq_formaldehyde_load, index_col='dt', parse_dates=True)
    df_udaq_formaldehyde.index = df_udaq_formaldehyde.index.tz_localize('UTC')

    df_udaq_formaldehyde_usos_only = df_udaq_formaldehyde.sort_index().loc['2024-07-15 00:00:00':'2024-08-18 23:59:00']
    df_udaq_formaldehyde_usos_only.index = df_udaq_formaldehyde_usos_only.index.rename('time_UTC')
    keep_colnames = ['H2CO_Corrected']
    df_udaq_formaldehyde_usos_only = df_udaq_formaldehyde_usos_only[keep_colnames]

    #This averaging in this section is based off the same code used for averaging the icartt data
    #as in the function align2master_timeline in time_utils.py of the icartt_read_and_merge

    # Get the average native sampling frequency in total seconds:
    tseries = df_udaq_formaldehyde_usos_only.index.to_series()
    #Avg native sampling frequency
    min_sep = int(np.round(tseries.diff().median().total_seconds()))
    #Intended interval
    step_S = time_interval_formaldehyde

    #Reindex to avg native sampling frequency (in this case, 1 minute)
    new_start_time = pd.Timestamp('2024-07-15 00:00:00').tz_localize('UTC')
    new_end_time = pd.Timestamp('2024-08-18 23:59:00').tz_localize('UTC')
    dts = pd.date_range(new_start_time, new_end_time, freq=str(min_sep) + 's')
    dfn = df_udaq_formaldehyde_usos_only.reindex(dts, method='nearest', fill_value=np.nan)

    # Take a centered boxcar average around the 900s (15 min) avg. (for numerical columns only)
    #NOTE: .mean() handles Nans like np.nanmean() in this context!!! 
    df_nums=dfn.select_dtypes(exclude=['datetime64'])
    df_nums_new = df_nums.rolling(str(int(step_S)) + 's').mean().resample(str(step_S) + 's').mean()

    df_avg_formaldehyde = df_nums_new
    #Convert to local time (MDT)
    df_avg_formaldehyde.index = df_avg_formaldehyde.index.tz_convert('America/Denver')
    df_avg_formaldehyde.index = df_avg_formaldehyde.index.rename('time_local')
    #print(df_new_formaldehyde)
    df_avg_formaldehyde = df_avg_formaldehyde.rename(columns = {'H2CO_Corrected': 'Formaldehyde'})
    df_avg_formaldehyde = df_avg_formaldehyde.loc['2024-07-15 00:00:00':'2024-08-18 17:45:00']

    # Set any inf, neg. inf, and negative values to NaN
    # Take the overlap in the ML and UDAQ data only
    series_nanconvert_formaldehyde = df_avg_formaldehyde['Formaldehyde'].replace([np.inf, -np.inf], np.nan).mask(df_avg_formaldehyde['Formaldehyde'] < 0, np.nan)
    
    mask_nans_formaldehyde = (~noaa['HCHO_CRDS'].isna()) & (~series_nanconvert_formaldehyde.isna())
    ml_clean_formaldehyde = noaa['HCHO_CRDS'][mask_nans_formaldehyde]
    udaq_clean_formaldehyde = series_nanconvert_formaldehyde[mask_nans_formaldehyde]

    udaq_overlap_formaldehyde_np_arr = udaq_clean_formaldehyde.to_numpy()
    #xdata_2d = np.array([udaq_overlap_formaldehyde_np_arr])
    ml_overlap_formaldehyde_np_arr = ml_clean_formaldehyde.to_numpy()
    #ydata_2d = np.array([ml_overlap_formaldehyde_np_arr])

    #print('xdata_2d shape: ', np.shape(udaq_overlap_formaldehyde_np_arr))
    #print('ydata shape: ', np.shape(ml_overlap_formaldehyde_np_arr))

    #From the SciPy and ODRpack documentation
    # Define the function we want to fit against
    def linear_model(x, beta):
        return beta[0] + beta[1] * x

    #beta0 is an initial guess
    odr_fit_result = odr_fit(f = linear_model, xdata = udaq_overlap_formaldehyde_np_arr, ydata = ml_overlap_formaldehyde_np_arr, beta0 = [0.0, 1.0])

    odr_intercept, odr_slope = odr_fit_result.beta

    # print('odr_fit_result: ', odr_fit_result)
    # print('odr_intercept: ', odr_intercept)
    # print('odr_slope: ', odr_slope)

    calibration_eq_formaldehyde = f"Calibration equation: O_hat = {odr_intercept:.4f} + {odr_slope:.4f} * M"
    #print(calibration_eq_formaldehyde)

    udaq_formaldehyde_adjusted = odr_intercept + odr_slope * series_nanconvert_formaldehyde

    #To determine if calibration adjustment aligns closer with ML data, use Root Mean Squared Error
    # Original RMSE
    rmse_original_formaldehyde = np.sqrt(np.mean((udaq_clean_formaldehyde - ml_clean_formaldehyde)**2))

    # Calibrated RMSE
    udaq_overlap_adjusted_with_initial = udaq_formaldehyde_adjusted.loc[udaq_formaldehyde_adjusted.index.intersection(udaq_clean_formaldehyde.index)]
    rmse_adjusted_formaldehyde = np.sqrt(np.mean((udaq_overlap_adjusted_with_initial - ml_clean_formaldehyde)**2))

    #print(f"Original RMSE: {rmse_original_formaldehyde:.4f}")
    #print(f"Calibrated RMSE: {rmse_adjusted_formaldehyde:.4f}")

    #Calculate percent improvement
    rmse_percent_improvement_formaldehyde = 100*((rmse_original_formaldehyde - rmse_adjusted_formaldehyde) / (rmse_original_formaldehyde))
    #print(f'Percent improvement of calibration adjusted RMSE: {rmse_percent_improvement_formaldehyde:.4f}')

    #Merge the old formaldehyde ML data with the UDAQ data with calibration adjustment applied to fill holes
    merged_formaldehyde = noaa['HCHO_CRDS'].fillna(udaq_formaldehyde_adjusted)

    #Compare stats before and after applying the adjustment to the UDAQ data
    slope_initial, intercept_initial = np.polyfit(udaq_clean_formaldehyde, ml_clean_formaldehyde, 1)
    slope_adjusted, intercept_adjusted = np.polyfit(udaq_overlap_adjusted_with_initial, ml_clean_formaldehyde, 1)

    #r squared calculation
    r_initial = np.corrcoef(udaq_clean_formaldehyde, ml_clean_formaldehyde)[0, 1]
    r2_initial = r_initial**2
    r_adjusted = np.corrcoef(udaq_overlap_adjusted_with_initial, ml_clean_formaldehyde)[0, 1]
    r2_adjusted = r_adjusted**2

    #Evaluation of adjustment validity
    did_rmse_improve = bool(rmse_adjusted_formaldehyde > rmse_original_formaldehyde)
    initial_slope_distance = abs(1-slope_initial)
    adjusted_slope_distance = abs(1-slope_adjusted)
    slope_distance_difference = adjusted_slope_distance - initial_slope_distance
    did_slope_improve = bool(adjusted_slope_distance < initial_slope_distance)

    #change variable names to label as Orthogonal Distance Regression applied
    ml_initial_formaldehyde_odr = noaa['HCHO_CRDS']
    ml_initial_formaldehyde_overlap_odr = ml_clean_formaldehyde

    udaq_initial_formaldehyde_odr = series_nanconvert_formaldehyde
    udaq_initial_formaldehyde_overlap_odr = udaq_clean_formaldehyde

    udaq_adjusted_formaldehyde_odr = udaq_formaldehyde_adjusted
    udaq_overlap_adjusted_formaldehyde_odr = udaq_overlap_adjusted_with_initial

    calibration_eq_formaldehyde_odr = calibration_eq_formaldehyde

    return (ml_initial_formaldehyde_odr, ml_initial_formaldehyde_overlap_odr, udaq_initial_formaldehyde_odr, 
    udaq_initial_formaldehyde_overlap_odr, udaq_adjusted_formaldehyde_odr, udaq_overlap_adjusted_formaldehyde_odr,
    calibration_eq_formaldehyde_odr, rmse_original_formaldehyde, rmse_adjusted_formaldehyde, did_rmse_improve, rmse_percent_improvement_formaldehyde, 
    slope_initial, initial_slope_distance, slope_adjusted, adjusted_slope_distance, slope_distance_difference, did_slope_improve, 
    intercept_initial, intercept_adjusted, r_initial, r2_initial , r_adjusted , r2_adjusted, merged_formaldehyde)
def plot_formaldehyde_only_calibration_adjustment_odr(ml_initial_formaldehyde_odr, ml_initial_formaldehyde_overlap_odr, udaq_initial_formaldehyde_odr, 
    udaq_initial_formaldehyde_overlap_odr, udaq_adjusted_formaldehyde_odr, udaq_overlap_adjusted_formaldehyde_odr,
    calibration_eq_formaldehyde_odr, merged_formaldehyde_odr):

    #First plot shows the overlap between initial ML, initial UDAQ, and fitted-UDAQ data
    #Second plot shows the initial ML, initial UDAQ, and fitted-UDAQ data regardless of data availability

    fig, (ax1, ax2) = plt.subplots(2,1, figsize = (16,8), tight_layout=True)
    xlim_start_jul = pd.to_datetime('2024-07-15 00:00:00').tz_localize('America/Denver')
    xlim_end_jul = pd.to_datetime('2024-07-31 23:00:00').tz_localize('America/Denver')
    xlim_start_aug = pd.to_datetime('2024-08-01 00:00:00').tz_localize('America/Denver')
    xlim_end_aug = pd.to_datetime('2024-08-18 23:45:00').tz_localize('America/Denver')

    #ax1 is the first row of subplot, for July only
    ax1.plot(ml_initial_formaldehyde_overlap_odr.index, ml_initial_formaldehyde_overlap_odr, linestyle = 'solid', color='m', marker='x',label='Mobile Lab', alpha = 0.7)
    ax1.plot(udaq_initial_formaldehyde_overlap_odr.index, udaq_initial_formaldehyde_overlap_odr,linestyle = 'solid', color = 'b', marker = '+', label = 'Initial UDAQ')
    ax1.plot(udaq_overlap_adjusted_formaldehyde_odr.index, udaq_overlap_adjusted_formaldehyde_odr, linestyle='solid', color = 'y', marker = '.', label = 'Calibrated UDAQ')

    #Set x ticks
    tz_mdt = noaa.index.tz #this time zone should be in Mountain Daylight Time
    ax1.xaxis.set_major_locator(mdates.DayLocator(tz=tz_mdt))
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d', tz=tz_mdt))
    # Minor ticks: every 3 hours
    ax1.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21]))
    # Rotate and format tick labels
    ax1.tick_params(axis='x', which='major')
    ax1.tick_params(axis='x', which='minor', length=3, color='gray')
    #ax.grid(True, which='both')
    
    ax1.set_ylabel('Formaldehyde (ppb)')
    #ax1.set_xlabel('Date')
    ax1.margins(x=0)
    ax1.set_xlim([xlim_start_jul, xlim_end_jul])

    ax1.legend(loc = 'upper right')
    
    #ax2 is the second row of subplot, for August only
    ax2.plot(ml_initial_formaldehyde_overlap_odr.index, ml_initial_formaldehyde_overlap_odr, linestyle = 'solid', color='m', marker='x',label='Mobile Lab', alpha = 0.7)
    ax2.plot(udaq_initial_formaldehyde_overlap_odr.index, udaq_initial_formaldehyde_overlap_odr,linestyle = 'solid', color = 'b', marker = '+', label = 'Initial UDAQ')
    ax2.plot(udaq_overlap_adjusted_formaldehyde_odr.index, udaq_overlap_adjusted_formaldehyde_odr, linestyle='solid', color = 'y', marker = '.', label = 'Calibrated UDAQ')

    #Set x ticks
    ax2.xaxis.set_major_locator(mdates.DayLocator(tz=tz_mdt))
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d', tz=tz_mdt))
    # Minor ticks: every 3 hours
    ax2.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21]))
    # Rotate and format tick labels
    ax2.tick_params(axis='x', which='major')
    ax2.tick_params(axis='x', which='minor', length=3, color='gray')

    ax2.set_ylabel('Formaldehyde (ppb)')
    ax2.set_xlabel('Time (MDT)')
    ax2.margins(x=0)

    ax2.set_xlim([xlim_start_aug, xlim_end_aug])
    ax2.legend(loc = 'upper right')

    #Mark midnight for every day
    midnight_vals = []
    for midnight_idx in range(0,len(noaa.index),96):
        midnight_vals.append(noaa.index[midnight_idx])
    for day_pos in midnight_vals:
        ax1.axvline(day_pos, color = 'black', linestyle = 'dotted')
        ax2.axvline(day_pos, color = 'black', linestyle = 'dotted')

    plt.suptitle(calibration_eq_formaldehyde_odr)
    plt.savefig(dirpath + '/Merge_scripts/calibration_adjustments/calibration_adjustment_plots/odr_fitting/overlap/odr_cal_adj_hawthorne_udaq_ml_Formaldehyde_overlap_comparison_july_aug_timeseries.png', dpi =300)
    plt.show()

    #PLOT 2
    fig, (ax1, ax2) = plt.subplots(2,1, figsize = (16,8), tight_layout=True)
    xlim_start_jul = pd.to_datetime('2024-07-15 00:00:00').tz_localize('America/Denver')
    xlim_end_jul = pd.to_datetime('2024-07-31 23:00:00').tz_localize('America/Denver')
    xlim_start_aug = pd.to_datetime('2024-08-01 00:00:00').tz_localize('America/Denver')
    xlim_end_aug = pd.to_datetime('2024-08-18 23:45:00').tz_localize('America/Denver')

    #ax1 is the first row of subplot, for July only
    ax1.plot(ml_initial_formaldehyde_overlap_odr.dropna().index, ml_initial_formaldehyde_overlap_odr.dropna(), linestyle = 'solid', color='m', marker='x',label='Mobile Lab', alpha = 0.7)
    ax1.plot(udaq_initial_formaldehyde_overlap_odr.dropna().index, udaq_initial_formaldehyde_overlap_odr.dropna(),linestyle = 'solid', color = 'b', marker = '+', label = 'Initial UDAQ')
    ax1.plot(udaq_overlap_adjusted_formaldehyde_odr.dropna().index, udaq_overlap_adjusted_formaldehyde_odr.dropna(), linestyle='solid', color = 'y', marker = '.', label = 'Calibrated UDAQ')

    #Set x ticks
    tz_mdt = noaa.index.tz #this time zone should be in Mountain Daylight Time
    ax1.xaxis.set_major_locator(mdates.DayLocator(tz=tz_mdt))
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d', tz=tz_mdt))
    # Minor ticks: every 3 hours
    ax1.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21]))
    # Rotate and format tick labels
    ax1.tick_params(axis='x', which='major')
    ax1.tick_params(axis='x', which='minor', length=3, color='gray')
    #ax.grid(True, which='both')
    
    ax1.set_ylabel('Formaldehyde (ppb)')
    #ax1.set_xlabel('Date')
    ax1.margins(x=0)
    ax1.set_xlim([xlim_start_jul, xlim_end_jul])

    ax1.legend(loc = 'upper right')
    
    #ax2 is the second row of subplot, for August only
    ax2.plot(ml_initial_formaldehyde_odr.dropna().index, ml_initial_formaldehyde_odr.dropna(), linestyle = 'solid', color='m', marker='x',label='Mobile Lab', alpha = 0.7)
    ax2.plot(udaq_initial_formaldehyde_odr.dropna().index, udaq_initial_formaldehyde_odr.dropna(),linestyle = 'solid', color = 'b', marker = '+', label = 'Initial UDAQ')
    ax2.plot(udaq_adjusted_formaldehyde_odr.dropna().index, udaq_adjusted_formaldehyde_odr.dropna(), linestyle='solid', color = 'y', marker = '.', label = 'Calibrated UDAQ')

    #Set x ticks
    ax2.xaxis.set_major_locator(mdates.DayLocator(tz=tz_mdt))
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d', tz=tz_mdt))
    # Minor ticks: every 3 hours
    ax2.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21]))
    # Rotate and format tick labels
    ax2.tick_params(axis='x', which='major')
    ax2.tick_params(axis='x', which='minor', length=3, color='gray')

    ax2.set_ylabel('Formaldehyde (ppb)')
    ax2.set_xlabel('Time (MDT)')
    ax2.margins(x=0)

    ax2.set_xlim([xlim_start_aug, xlim_end_aug])
    ax2.legend(loc = 'upper right')

    #Mark midnight for every day
    midnight_vals = []
    for midnight_idx in range(0,len(noaa.index),96):
        midnight_vals.append(noaa.index[midnight_idx])
    for day_pos in midnight_vals:
        ax1.axvline(day_pos, color = 'black', linestyle = 'dotted')
        ax2.axvline(day_pos, color = 'black', linestyle = 'dotted')

    plt.suptitle(calibration_eq_formaldehyde_odr)
    plt.savefig(dirpath + '/Merge_scripts/calibration_adjustments/calibration_adjustment_plots/odr_fitting/full_time_series/odr_cal_adj_hawthorne_udaq_ml_Formaldehyde_overlap_comparison_july_aug_timeseries.png', dpi =300)
    plt.show()

def calibration_adjustment_ols(noaa, udaq, noaa_species, udaq_species, time_interval_formaldehyde_used):
    if 'ML_NoVar' in noaa_species:
        print('Skipping ', udaq_species, ' due to not being a var in ML Data')
        noaa_species_name = noaa_species
        udaq_species_name = udaq_species
        ml_initial_ols = np.nan
        ml_initial_overlap_ols = np.nan
        udaq_initial_ols = udaq[udaq_species]
        udaq_initial_overlap_ols = udaq[udaq_species]
        udaq_adjusted_ols = udaq[udaq_species]
        udaq_overlap_adjusted_ols = udaq[udaq_species]
        calibration_eq_ols = np.nan
        rmse_original_ols = np.nan
        rmse_calibrated_ols = np.nan
        rmse_percent_improvement_ols = np.nan
        #udaq_vals_adjusted = udaq[udaq_species]
        merged_voc_ols = udaq[udaq_species]
    
    elif 'UDAQ_NoVar' in udaq_species:
        print('Skipping ', noaa_species, ' due to not being a var in UDAQ Data')
        noaa_species_name = noaa_species
        udaq_species_name = udaq_species
        ml_initial_ols = noaa[noaa_species_name]
        ml_initial_overlap_ols = noaa[noaa_species_name]
        udaq_initial_ols = np.nan
        udaq_initial_overlap_ols = np.nan
        udaq_adjusted_ols = np.nan
        udaq_overlap_adjusted_ols = np.nan
        calibration_eq_ols = np.nan
        rmse_original_ols = np.nan
        rmse_calibrated_ols = np.nan
        rmse_percent_improvement_ols = np.nan
        #udaq_vals_adjusted = noaa[noaa_species_name]
        merged_voc_ols = noaa[noaa_species_name]

    elif 'Formaldehyde' in udaq_species:
        print('Calculating Formaldehyde calibration adjustment.')
        (ml_initial_formaldehyde_ols, ml_initial_formaldehyde_overlap_ols, 
         udaq_initial_formaldehyde_ols, udaq_initial_formaldehyde_overlap_ols, 
         udaq_adjusted_formaldehyde_ols, udaq_overlap_adjusted_formaldehyde_ols,
         calibration_eq_formaldehyde_ols, rmse_original_formaldehyde_ols, 
         rmse_calibrated_formaldehyde_ols, rmse_percent_improvement_formaldehyde_ols, merged_formaldehyde_ols) = formaldehyde_only_calibration_adjustment_ols(time_interval_formaldehyde = time_interval_formaldehyde_used)
        
        noaa_species_name = 'HCHO_CRDS'
        udaq_species_name = udaq_species
        ml_initial_ols = ml_initial_formaldehyde_ols
        ml_initial_overlap_ols = ml_initial_formaldehyde_overlap_ols
        udaq_initial_ols = udaq_initial_formaldehyde_ols
        udaq_initial_overlap_ols = udaq_initial_formaldehyde_overlap_ols
        udaq_adjusted_ols = udaq_adjusted_formaldehyde_ols
        udaq_overlap_adjusted_ols = udaq_overlap_adjusted_formaldehyde_ols
        calibration_eq_ols = calibration_eq_formaldehyde_ols
        rmse_original_ols = rmse_original_formaldehyde_ols
        rmse_calibrated_ols = rmse_calibrated_formaldehyde_ols
        rmse_percent_improvement_ols = rmse_percent_improvement_formaldehyde_ols
        #udaq_vals_adjusted = udaq_formaldehyde_adjusted
        merged_voc_ols = merged_formaldehyde_ols

    else:
        print('Applying calibration adjustment for ' + str(udaq_species))
        noaa_species_name = noaa_species
        udaq_species_name = udaq_species
        df_vars = pd.DataFrame({'ML Data': noaa[noaa_species_name], 'UDAQ Data': udaq[udaq_species_name]})
        df_vars_clean = df_vars.dropna() #Only include overlap with no NaNs

        if len(df_vars_clean) == 0:
            print(udaq_species, ' has no UDAQ Data (all NaNs), use ML Data')
            ml_initial_ols = noaa[noaa_species_name]
            ml_initial_overlap_ols = noaa[noaa_species_name]
            udaq_initial_ols = np.nan
            udaq_initial_overlap_ols = np.nan
            udaq_adjusted_ols = np.nan
            udaq_overlap_adjusted_ols = np.nan
            calibration_eq_ols = np.nan
            rmse_original_ols = np.nan
            rmse_calibrated_ols = np.nan
            rmse_percent_improvement_ols = np.nan
            #udaq_vals_adjusted = noaa[noaa_species_name]
            merged_voc_ols = noaa[noaa_species_name]

        else:
            M_with_const = sm.add_constant(df_vars_clean['UDAQ Data'])

            # Fit OLS regression: O = a + b*M
            ols_calibration = sm.OLS(df_vars_clean['ML Data'], M_with_const).fit()
            # Get calibration coefficients
            a = ols_calibration.params['const']  # additive offset
            b = ols_calibration.params[df_vars_clean['UDAQ Data'].name]   # multiplicative factor
            calibration_eq = f"Calibration equation: O_hat = {a:.4f} + {b:.4f} * M"
            #print(calibration_eq)
            #print(ols_calibration.summary())
            udaq_vals_adjusted = a+b*df_vars['UDAQ Data']

            #To determine if calibration adjustment aligns closer with ML data, use Root Mean Squared Error
            # Original RMSE
            rmse_original = np.sqrt(np.mean((df_vars_clean['UDAQ Data'] - df_vars_clean['ML Data'])**2))

            # Calibrated RMSE
            udaq_overlap_adjusted_with_initial = udaq_vals_adjusted.loc[udaq_vals_adjusted.index.intersection(df_vars_clean.index)]
            rmse_calibrated = np.sqrt(np.mean((udaq_overlap_adjusted_with_initial - df_vars_clean['ML Data'])**2))

            #print(f"Original RMSE: {rmse_original:.4f}")
            print(f"Calibration adjusted RMSE: {rmse_calibrated:.4f}")

            #Calculate percent improvement
            rmse_percent_improvement = 100*((rmse_original - rmse_calibrated) / (rmse_original))
            print(f'Percent improvement of calibration adjusted RMSE: {rmse_percent_improvement:.4f}')

            merged_voc = noaa[noaa_species_name].fillna(udaq_vals_adjusted)
            print('Completed filling the holes in ML Data for '+ str(noaa_species_name))

            #Compare stats before and after applying the adjustment to the UDAQ data
            slope_initial, intercept_initial, r_initial, p_initial, stderr_initial = stats.linregress(df_vars_clean['UDAQ Data'], df_vars_clean['ML Data'])
            slope_adjusted, intercept_adjusted, r_adjusted, p_adjusted, stderr_adjusted = stats.linregress(udaq_overlap_adjusted_with_initial, df_vars_clean['ML Data'])

            #change variable names to label as Ordinary Least Squares applied
            ml_initial_ols = df_vars['ML Data']
            ml_initial_overlap_ols = df_vars_clean['ML Data']
            udaq_initial_ols = df_vars['UDAQ Data']
            udaq_initial_overlap_ols = df_vars_clean['UDAQ Data']
            udaq_adjusted_ols = udaq_vals_adjusted
            udaq_overlap_adjusted_ols = udaq_overlap_adjusted_with_initial
            calibration_eq_ols = calibration_eq
            rmse_original_ols = rmse_original
            rmse_calibrated_ols = rmse_calibrated
            rmse_percent_improvement_ols = rmse_percent_improvement
            merged_voc_ols = merged_voc

    return (noaa_species_name, udaq_species_name, ml_initial_ols, ml_initial_overlap_ols, udaq_initial_ols, 
    udaq_initial_overlap_ols, udaq_adjusted_ols, udaq_overlap_adjusted_ols,
    calibration_eq_ols, rmse_original_ols, rmse_calibrated_ols, rmse_percent_improvement_ols, merged_voc_ols)
def plot_calibration_adjustment_ols(udaq_species_name, ml_initial_ols, ml_initial_overlap_ols, udaq_initial_ols, udaq_initial_overlap_ols, udaq_adjusted_ols, udaq_overlap_adjusted_ols, calibration_eq_ols, merged_voc_ols):
    fig, (ax1, ax2) = plt.subplots(2,1, figsize = (16,8), tight_layout=True)
    xlim_start_jul = pd.to_datetime('2024-07-15 00:00:00').tz_localize('America/Denver')
    xlim_end_jul = pd.to_datetime('2024-07-31 23:00:00').tz_localize('America/Denver')
    xlim_start_aug = pd.to_datetime('2024-08-01 00:00:00').tz_localize('America/Denver')
    xlim_end_aug = pd.to_datetime('2024-08-18 23:45:00').tz_localize('America/Denver')

    #ax1 is the first row of subplot, for July only
    ax1.plot(ml_initial_overlap_ols.index, ml_initial_overlap_ols, linestyle = 'solid', color='m', marker='x',label='Mobile Lab', alpha = 0.7)
    ax1.plot(udaq_initial_overlap_ols.index, udaq_initial_overlap_ols, linestyle = 'solid', color = 'b', marker = '+', label = 'Initial UDAQ')
    ax1.plot(udaq_overlap_adjusted_ols.index, udaq_overlap_adjusted_ols, linestyle='solid', color = 'y', marker = '.', label = 'Calibrated UDAQ')

    #Set x ticks
    tz_mdt = noaa.index.tz #this time zone should be in Mountain Daylight Time
    ax1.xaxis.set_major_locator(mdates.DayLocator(tz=tz_mdt))
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d', tz=tz_mdt))
    # Minor ticks: every 3 hours
    ax1.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21]))
    # Rotate and format tick labels
    ax1.tick_params(axis='x', which='major')
    ax1.tick_params(axis='x', which='minor', length=3, color='gray')
    #ax.grid(True, which='both')
    
    ax1.set_ylabel(udaq_species_name + ' (ppb)')
    #ax1.set_xlabel('Date')
    ax1.margins(x=0)
    ax1.set_xlim([xlim_start_jul, xlim_end_jul])

    ax1.legend(loc = 'upper right')
    
    #ax2 is the second row of subplot, for August only
    ax2.plot(ml_initial_overlap_ols.index, ml_initial_overlap_ols, linestyle = 'solid', color='m', marker='x',label='Mobile Lab', alpha = 0.7)
    ax2.plot(udaq_initial_overlap_ols.index, udaq_initial_overlap_ols, linestyle = 'solid', color = 'b', marker = '+', label = 'Initial UDAQ')
    ax2.plot(udaq_overlap_adjusted_ols.index, udaq_overlap_adjusted_ols, linestyle='solid', color = 'y', marker = '.', label = 'Calibrated UDAQ')
    #Set x ticks
    ax2.xaxis.set_major_locator(mdates.DayLocator(tz=tz_mdt))
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d', tz=tz_mdt))
    # Minor ticks: every 3 hours
    ax2.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21]))
    # Rotate and format tick labels
    ax2.tick_params(axis='x', which='major')
    ax2.tick_params(axis='x', which='minor', length=3, color='gray')

    ax2.set_ylabel(udaq_species_name + ' (ppb)')
    ax2.set_xlabel('Time (MDT)')
    ax2.margins(x=0)

    ax2.set_xlim([xlim_start_aug, xlim_end_aug])
    ax2.legend(loc = 'upper right')

    #Mark midnight for every day
    midnight_vals = []
    for midnight_idx in range(0,len(noaa.index),96):
        midnight_vals.append(noaa.index[midnight_idx])
    for day_pos in midnight_vals:
        ax1.axvline(day_pos, color = 'black', linestyle = 'dotted')
        ax2.axvline(day_pos, color = 'black', linestyle = 'dotted')

    plt.suptitle(calibration_eq_ols)
    plt.savefig(dirpath + '/Merge_scripts/calibration_adjustments/calibration_adjustment_plots/ols_fitting/ols_cal_adj_hawthorne_udaq_ml_'+ str(udaq_species_name) + '_overlap_comparison_july_aug_timeseries.png', dpi =300)
    plt.show()

    #PLOT 2
    fig, (ax1, ax2) = plt.subplots(2,1, figsize = (16,8), tight_layout=True)
    xlim_start_jul = pd.to_datetime('2024-07-15 00:00:00').tz_localize('America/Denver')
    xlim_end_jul = pd.to_datetime('2024-07-31 23:00:00').tz_localize('America/Denver')
    xlim_start_aug = pd.to_datetime('2024-08-01 00:00:00').tz_localize('America/Denver')
    xlim_end_aug = pd.to_datetime('2024-08-18 23:45:00').tz_localize('America/Denver')

    #ax1 is the first row of subplot, for July only
    ax1.plot(ml_initial_ols.dropna().index, ml_initial_ols.dropna(), linestyle = 'solid', color='m', marker='x',label='Mobile Lab', alpha = 0.7)
    ax1.plot(udaq_initial_ols.dropna().index, udaq_initial_ols.dropna(), linestyle = 'solid', color = 'b', marker = '+', label = 'Initial UDAQ')
    ax1.plot(udaq_adjusted_ols.dropna().index, udaq_adjusted_ols.dropna(), linestyle='solid', color = 'y', marker = '.', label = 'Calibrated UDAQ')

    #Set x ticks
    tz_mdt = noaa.index.tz #this time zone should be in Mountain Daylight Time
    ax1.xaxis.set_major_locator(mdates.DayLocator(tz=tz_mdt))
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d', tz=tz_mdt))
    # Minor ticks: every 3 hours
    ax1.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21]))
    # Rotate and format tick labels
    ax1.tick_params(axis='x', which='major')
    ax1.tick_params(axis='x', which='minor', length=3, color='gray')
    #ax.grid(True, which='both')
    
    ax1.set_ylabel(udaq_species_name + ' (ppb)')
    #ax1.set_xlabel('Date')
    ax1.margins(x=0)
    ax1.set_xlim([xlim_start_jul, xlim_end_jul])

    ax1.legend(loc = 'upper right')
    
    #ax2 is the second row of subplot, for August only
    ax2.plot(ml_initial_ols.dropna().index, ml_initial_ols.dropna(), linestyle = 'solid', color='m', marker='x',label='Mobile Lab', alpha = 0.7)
    ax2.plot(udaq_initial_ols.dropna().index, udaq_initial_ols.dropna(), linestyle = 'solid', color = 'b', marker = '+', label = 'Initial UDAQ')
    ax2.plot(udaq_adjusted_ols.dropna().index, udaq_adjusted_ols.dropna(), linestyle='solid', color = 'y', marker = '.', label = 'Calibrated UDAQ')
    #Set x ticks
    ax2.xaxis.set_major_locator(mdates.DayLocator(tz=tz_mdt))
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d', tz=tz_mdt))
    # Minor ticks: every 3 hours
    ax2.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21]))
    # Rotate and format tick labels
    ax2.tick_params(axis='x', which='major')
    ax2.tick_params(axis='x', which='minor', length=3, color='gray')

    ax2.set_ylabel(udaq_species_name + ' (ppb)')
    ax2.set_xlabel('Time (MDT)')
    ax2.margins(x=0)

    ax2.set_xlim([xlim_start_aug, xlim_end_aug])
    ax2.legend(loc = 'upper right')

    #Mark midnight for every day
    midnight_vals = []
    for midnight_idx in range(0,len(noaa.index),96):
        midnight_vals.append(noaa.index[midnight_idx])
    for day_pos in midnight_vals:
        ax1.axvline(day_pos, color = 'black', linestyle = 'dotted')
        ax2.axvline(day_pos, color = 'black', linestyle = 'dotted')

    plt.suptitle(calibration_eq_ols)
    plt.savefig(dirpath + '/Merge_scripts/calibration_adjustments/calibration_adjustment_plots/ols_fitting/full_time_series/ols_cal_adj_hawthorne_udaq_ml_'+ str(udaq_species_name) + '_overlap_comparison_july_aug_timeseries.png', dpi =300)
    plt.show()

def calibration_adjustment_odr(noaa, udaq, noaa_species, udaq_species, time_interval_formaldehyde_used):
    if 'ML_NoVar' in noaa_species:
        print('Skipping ', udaq_species, ' calibration due to not being a var in ML Data')
        noaa_species_name = noaa_species
        udaq_species_name = udaq_species
        ml_initial_odr = np.nan
        ml_initial_overlap_odr = np.nan
        udaq_initial_odr = udaq[udaq_species]
        udaq_initial_overlap_odr = udaq[udaq_species]
        udaq_adjusted_odr = udaq[udaq_species]
        udaq_overlap_adjusted_odr = udaq[udaq_species]
        calibration_eq_odr = np.nan
        rmse_original = np.nan
        rmse_adjusted = np.nan
        did_rmse_improve = np.nan
        rmse_percent_improvement = np.nan

        slope_initial = np.nan
        initial_slope_distance = np.nan
        slope_adjusted = np.nan
        adjusted_slope_distance = np.nan
        slope_distance_difference = np.nan
        did_slope_improve = np.nan
        
        intercept_initial = np.nan
        intercept_adjusted = np.nan
        
        r_initial = np.nan
        r2_initial = np.nan
        r_adjusted = np.nan
        r2_adjusted = np.nan
        merged_voc = udaq[udaq_species]
    
    elif 'UDAQ_NoVar' in udaq_species:
        print('Skipping ', noaa_species, 'calibration due to not being a var in UDAQ Data')
        noaa_species_name = noaa_species
        udaq_species_name = udaq_species
        ml_initial_odr = noaa[noaa_species_name]
        ml_initial_overlap_odr = noaa[noaa_species_name]
        udaq_initial_odr = np.nan
        udaq_initial_overlap_odr = np.nan
        udaq_adjusted_odr = np.nan
        udaq_overlap_adjusted_odr = np.nan
        calibration_eq_odr = np.nan

        rmse_original = np.nan
        rmse_adjusted = np.nan
        did_rmse_improve = np.nan
        rmse_percent_improvement = np.nan

        slope_initial = np.nan
        initial_slope_distance = np.nan
        slope_adjusted = np.nan
        adjusted_slope_distance = np.nan
        slope_distance_difference = np.nan
        did_slope_improve = np.nan
        
        intercept_initial = np.nan
        intercept_adjusted = np.nan
        
        r_initial = np.nan
        r2_initial = np.nan
        r_adjusted = np.nan
        r2_adjusted = np.nan

        merged_voc = noaa[noaa_species_name]

    elif 'Formaldehyde' in udaq_species:
        print('Calculating Formaldehyde calibration.')
        (ml_initial_formaldehyde_odr, ml_initial_formaldehyde_overlap_odr, udaq_initial_formaldehyde_odr, 
        udaq_initial_formaldehyde_overlap_odr, udaq_adjusted_formaldehyde_odr, udaq_overlap_adjusted_formaldehyde_odr,
        calibration_eq_formaldehyde_odr, rmse_original_formaldehyde, 
        rmse_adjusted_formaldehyde, did_rmse_improve, rmse_percent_improvement_formaldehyde, 
        slope_initial, initial_slope_distance, slope_adjusted, adjusted_slope_distance, slope_distance_difference, did_slope_improve, 
        intercept_initial, intercept_adjusted, r_initial, r2_initial , r_adjusted , r2_adjusted, merged_formaldehyde)= formaldehyde_only_calibration_adjustment_odr(time_interval_formaldehyde = time_interval_formaldehyde_used)
        
        noaa_species_name = 'HCHO_CRDS'
        udaq_species_name = udaq_species
        ml_initial_odr = ml_initial_formaldehyde_odr
        ml_initial_overlap_odr = ml_initial_formaldehyde_overlap_odr
        udaq_initial_odr = udaq_initial_formaldehyde_odr
        udaq_initial_overlap_odr = udaq_initial_formaldehyde_overlap_odr
        udaq_adjusted_odr = udaq_adjusted_formaldehyde_odr
        udaq_overlap_adjusted_odr = udaq_overlap_adjusted_formaldehyde_odr
        calibration_eq_odr = calibration_eq_formaldehyde_odr
        rmse_original = rmse_original_formaldehyde
        rmse_adjusted = rmse_adjusted_formaldehyde

        rmse_percent_improvement = rmse_percent_improvement_formaldehyde
        merged_voc = merged_formaldehyde

    else:
        print('Applying calibration adjustment for ' + str(udaq_species))
        noaa_species_name = noaa_species
        udaq_species_name = udaq_species
        df_vars = pd.DataFrame({'ML Data': noaa[noaa_species_name], 'UDAQ Data': udaq[udaq_species_name]})
        df_vars_clean = df_vars.dropna() #Only include overlap with no NaNs
    
        if len(df_vars_clean) == 0:
            print(udaq_species, ' has no UDAQ Data (all NaNs), use ML Data')
            ml_initial_odr = noaa[noaa_species_name]
            ml_initial_overlap_odr = noaa[noaa_species_name]
            udaq_initial_odr = np.nan
            udaq_initial_overlap_odr = np.nan
            udaq_adjusted_odr = np.nan
            udaq_overlap_adjusted_odr = np.nan
            calibration_eq_odr = np.nan
            rmse_original = np.nan
            rmse_adjusted = np.nan
            did_rmse_improve = np.nan
            rmse_percent_improvement = np.nan

            slope_initial = np.nan
            initial_slope_distance = np.nan
            slope_adjusted = np.nan
            adjusted_slope_distance = np.nan
            slope_distance_difference = np.nan
            did_slope_improve = np.nan
            
            intercept_initial = np.nan
            intercept_adjusted = np.nan
            
            r_initial = np.nan
            r2_initial = np.nan
            r_adjusted = np.nan
            r2_adjusted = np.nan

            merged_voc = noaa[noaa_species_name]

        else:
            #From the ODRpack documentation
            # Define the function we want to fit against
            udaq_overlap_voc_np_arr = df_vars_clean['UDAQ Data'].to_numpy()
            #xdata_2d = np.array([udaq_overlap_voc_np_arr])
            ml_overlap_voc_np_arr = df_vars_clean['ML Data'].to_numpy()
            #ydata_2d = np.array([ml_overlap_voc_np_arr])

            #print('xdata_2d shape: ', np.shape(udaq_overlap_voc_np_arr))
            #print('ydata shape: ', np.shape(ml_overlap_voc_np_arr))

            def linear_model(x, beta):
                return beta[0] + beta[1] * x
            
            #beta0 is an initial guess
            odr_fit_result = odr_fit(linear_model, udaq_overlap_voc_np_arr, ml_overlap_voc_np_arr, [0.0, 1.0])

            odr_intercept, odr_slope = odr_fit_result.beta

            calibration_eq = f"Calibration equation: O_hat = {odr_intercept:.4f} + {odr_slope:.4f} * M"
            #print(calibration_eq)

            udaq_vals_adjusted = odr_intercept + odr_slope * df_vars['UDAQ Data']
            #overlap between initial and fitted data
            udaq_overlap_adjusted_with_initial = udaq_vals_adjusted.loc[udaq_vals_adjusted.index.intersection(df_vars_clean.index)]
        
            #To determine if calibration adjustment aligns closer with ML data, use Root Mean Squared Error
            # Original RMSE
            rmse_original = np.sqrt(np.mean((df_vars_clean['UDAQ Data'] - df_vars_clean['ML Data'])**2))

            # Calibrated RMSE
            rmse_adjusted = np.sqrt(np.mean((udaq_overlap_adjusted_with_initial - df_vars_clean['ML Data'])**2))

            #print(f"Original RMSE: {rmse_original:.4f}")
            #print(f"Calibration adjusted RMSE: {rmse_adjusted:.4f}")

            #Calculate percent improvement
            rmse_percent_improvement = 100*((rmse_original - rmse_adjusted) / (rmse_original))
            #print(f'Percent improvement of calibration adjusted RMSE: {rmse_percent_improvement:.4f}')

            merged_voc = noaa[noaa_species_name].fillna(udaq_vals_adjusted)
            print('Completed filling the holes in ML Data for '+ str(noaa_species_name))

            #Compare stats before and after applying the adjustment to the UDAQ data
            slope_initial, intercept_initial = np.polyfit(df_vars_clean['UDAQ Data'], df_vars_clean['ML Data'], 1)
            slope_adjusted, intercept_adjusted = np.polyfit(udaq_overlap_adjusted_with_initial, df_vars_clean['ML Data'], 1)
            #r squared calculation
            r_initial = np.corrcoef(df_vars_clean['UDAQ Data'], df_vars_clean['ML Data'])[0, 1]
            r2_initial = r_initial**2
            r_adjusted = np.corrcoef(udaq_overlap_adjusted_with_initial, df_vars_clean['ML Data'])[0, 1]
            r2_adjusted = r_adjusted**2

            #Evaluation of adjustment validity
            did_rmse_improve = bool(rmse_adjusted > rmse_original)
            initial_slope_distance = abs(1-slope_initial)
            adjusted_slope_distance = abs(1-slope_adjusted)
            slope_distance_difference = adjusted_slope_distance - initial_slope_distance
            did_slope_improve = bool(adjusted_slope_distance < initial_slope_distance)
    
            #change variable names to label as ODR
            ml_initial_odr = df_vars['ML Data']
            ml_initial_overlap_odr = df_vars_clean['ML Data']
            udaq_initial_odr = df_vars['UDAQ Data']
            udaq_initial_overlap_odr = df_vars_clean['UDAQ Data']
            udaq_adjusted_odr = udaq_vals_adjusted
            udaq_overlap_adjusted_odr = udaq_overlap_adjusted_with_initial
            calibration_eq_odr = calibration_eq

    return (noaa_species_name, udaq_species_name, ml_initial_odr, ml_initial_overlap_odr, udaq_initial_odr, 
    udaq_initial_overlap_odr, udaq_adjusted_odr, udaq_overlap_adjusted_odr, calibration_eq_odr,
    rmse_original, rmse_adjusted, did_rmse_improve, rmse_percent_improvement, 
    slope_initial, initial_slope_distance, slope_adjusted, adjusted_slope_distance, slope_distance_difference, did_slope_improve, 
    intercept_initial, intercept_adjusted, r_initial, r2_initial, r_adjusted, r2_adjusted, merged_voc)

def plot_calibration_adjustment_odr(udaq_species_name, ml_initial_odr, ml_initial_overlap_odr, udaq_initial_odr, udaq_initial_overlap_odr, udaq_adjusted_odr, udaq_overlap_adjusted_odr, calibration_eq_odr, merged_voc_odr):
    #Plot 1: Time Series
    fig, (ax1, ax2) = plt.subplots(2,1, figsize = (16,8), tight_layout=True)
    xlim_start_jul = pd.to_datetime('2024-07-15 00:00:00').tz_localize('America/Denver')
    xlim_end_jul = pd.to_datetime('2024-07-31 23:00:00').tz_localize('America/Denver')
    xlim_start_aug = pd.to_datetime('2024-08-01 00:00:00').tz_localize('America/Denver')
    xlim_end_aug = pd.to_datetime('2024-08-18 23:45:00').tz_localize('America/Denver')

    #ax1 is the first row of subplot, for July only
    ax1.plot(ml_initial_overlap_odr.index, ml_initial_overlap_odr, linestyle = 'solid', color='m', marker='x',label='Mobile Lab', alpha = 0.7)
    ax1.plot(udaq_initial_overlap_odr.index, udaq_initial_overlap_odr, linestyle = 'solid', color = 'b', marker = '+', label = 'Initial UDAQ')
    ax1.plot(udaq_overlap_adjusted_odr.index, udaq_overlap_adjusted_odr, linestyle='solid', color = 'y', marker = '.', label = 'Calibrated UDAQ')

    #Set x ticks
    tz_mdt = noaa.index.tz #this time zone should be in Mountain Daylight Time
    ax1.xaxis.set_major_locator(mdates.DayLocator(tz=tz_mdt))
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d', tz=tz_mdt))
    # Minor ticks: every 3 hours
    ax1.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21]))
    # Rotate and format tick labels
    ax1.tick_params(axis='x', which='major')
    ax1.tick_params(axis='x', which='minor', length=3, color='gray')
    #ax.grid(True, which='both')
    
    ax1.set_ylabel(udaq_species_name + ' (ppb)')
    #ax1.set_xlabel('Date')
    ax1.margins(x=0)
    ax1.set_xlim([xlim_start_jul, xlim_end_jul])

    ax1.legend(loc = 'upper right')
    
    #ax2 is the second row of subplot, for August only
    ax2.plot(ml_initial_overlap_odr.index, ml_initial_overlap_odr, linestyle = 'solid', color='m', marker='x',label='Mobile Lab', alpha = 0.7)
    ax2.plot(udaq_initial_overlap_odr.index, udaq_initial_overlap_odr, linestyle = 'solid', color = 'b', marker = '+', label = 'Initial UDAQ')
    ax2.plot(udaq_overlap_adjusted_odr.index, udaq_overlap_adjusted_odr, linestyle='solid', color = 'y', marker = '.', label = 'Calibrated UDAQ')

    #Set x ticks
    ax2.xaxis.set_major_locator(mdates.DayLocator(tz=tz_mdt))
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d', tz=tz_mdt))
    # Minor ticks: every 3 hours
    ax2.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21]))
    # Rotate and format tick labels
    ax2.tick_params(axis='x', which='major')
    ax2.tick_params(axis='x', which='minor', length=3, color='gray')

    ax2.set_ylabel(udaq_species_name + ' (ppb)')
    ax2.set_xlabel('Time (MDT)')
    ax2.margins(x=0)

    ax2.set_xlim([xlim_start_aug, xlim_end_aug])
    ax2.legend(loc = 'upper right')

    #Mark midnight for every day
    midnight_vals = []
    for midnight_idx in range(0,len(noaa.index),96):
        midnight_vals.append(noaa.index[midnight_idx])
    for day_pos in midnight_vals:
        ax1.axvline(day_pos, color = 'black', linestyle = 'dotted')
        ax2.axvline(day_pos, color = 'black', linestyle = 'dotted')

    plt.suptitle(calibration_eq_odr)
    plt.savefig(dirpath + '/Merge_scripts/calibration_adjustments/calibration_adjustment_plots/odr_fitting/overlap/odr_cal_adj_hawthorne_udaq_ml_'+ str(udaq_species_name) + '_overlap_comparison_july_aug_timeseries.png', dpi =300)
    plt.show()

    #PLOT 2
    fig, (ax1, ax2) = plt.subplots(2,1, figsize = (16,8), tight_layout=True)
    xlim_start_jul = pd.to_datetime('2024-07-15 00:00:00').tz_localize('America/Denver')
    xlim_end_jul = pd.to_datetime('2024-07-31 23:00:00').tz_localize('America/Denver')
    xlim_start_aug = pd.to_datetime('2024-08-01 00:00:00').tz_localize('America/Denver')
    xlim_end_aug = pd.to_datetime('2024-08-18 23:45:00').tz_localize('America/Denver')

    #ax1 is the first row of subplot, for July only
    ax1.plot(ml_initial_odr.dropna().index, ml_initial_odr.dropna(), linestyle = 'solid', color='m', marker='x',label='Mobile Lab', alpha = 0.7)
    ax1.plot(udaq_initial_odr.dropna().index, udaq_initial_odr.dropna(), linestyle = 'solid', color = 'b', marker = '+', label = 'Initial UDAQ')
    ax1.plot(udaq_adjusted_odr.dropna().index, udaq_adjusted_odr.dropna(), linestyle='solid', color = 'y', marker = '.', label = 'Calibrated UDAQ')

    #Set x ticks
    tz_mdt = noaa.index.tz #this time zone should be in Mountain Daylight Time
    ax1.xaxis.set_major_locator(mdates.DayLocator(tz=tz_mdt))
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d', tz=tz_mdt))
    # Minor ticks: every 3 hours
    ax1.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21]))
    # Rotate and format tick labels
    ax1.tick_params(axis='x', which='major')
    ax1.tick_params(axis='x', which='minor', length=3, color='gray')
    #ax.grid(True, which='both')
    
    ax1.set_ylabel(udaq_species_name + ' (ppb)')
    #ax1.set_xlabel('Date')
    ax1.margins(x=0)
    ax1.set_xlim([xlim_start_jul, xlim_end_jul])

    ax1.legend(loc = 'upper right')
    
    #ax2 is the second row of subplot, for August only
    ax2.plot(ml_initial_odr.dropna().index, ml_initial_odr.dropna(), linestyle = 'solid', color='m', marker='x',label='Mobile Lab', alpha = 0.7)
    ax2.plot(udaq_initial_odr.dropna().index, udaq_initial_odr.dropna(), linestyle = 'solid', color = 'b', marker = '+', label = 'Initial UDAQ')
    ax2.plot(udaq_adjusted_odr.dropna().index, udaq_adjusted_odr.dropna(), linestyle='solid', color = 'y', marker = '.', label = 'Calibrated UDAQ')
    #Set x ticks
    ax2.xaxis.set_major_locator(mdates.DayLocator(tz=tz_mdt))
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d', tz=tz_mdt))
    # Minor ticks: every 3 hours
    ax2.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21]))
    # Rotate and format tick labels
    ax2.tick_params(axis='x', which='major')
    ax2.tick_params(axis='x', which='minor', length=3, color='gray')

    ax2.set_ylabel(udaq_species_name + ' (ppb)')
    ax2.set_xlabel('Time (MDT)')
    ax2.margins(x=0)

    ax2.set_xlim([xlim_start_aug, xlim_end_aug])
    ax2.legend(loc = 'upper right')

    #Mark midnight for every day
    midnight_vals = []
    for midnight_idx in range(0,len(noaa.index),96):
        midnight_vals.append(noaa.index[midnight_idx])
    for day_pos in midnight_vals:
        ax1.axvline(day_pos, color = 'black', linestyle = 'dotted')
        ax2.axvline(day_pos, color = 'black', linestyle = 'dotted')

    plt.suptitle(calibration_eq_odr)
    plt.savefig(dirpath + '/Merge_scripts/calibration_adjustments/calibration_adjustment_plots/odr_fitting/full_time_series/odr_cal_adj_hawthorne_udaq_ml_'+ str(udaq_species_name) + '_overlap_comparison_july_aug_timeseries.png', dpi =300)
    plt.show()

#region: older code for applying calibration adjustment
list_calibration_adj_total_ols = []
list_calibration_adj_total_odr = []
data_merged_dict_ols = {}
data_merged_dict_odr = {}
data_merged_overlap_dict_odr = {}
for voc_spec in mapping.keys():
    (noaa_species_name, udaq_species_name, ml_initial_odr, ml_initial_overlap_odr, udaq_initial_odr, 
    udaq_initial_overlap_odr, udaq_adjusted_odr, udaq_overlap_adjusted_odr, calibration_eq_odr,
    rmse_original, rmse_adjusted, did_rmse_improve, rmse_percent_improvement, 
    slope_initial, initial_slope_distance, slope_adjusted, adjusted_slope_distance, slope_distance_difference, did_slope_improve, 
    intercept_initial, intercept_adjusted, r_initial, r2_initial, r_adjusted, r2_adjusted, merged_voc) = calibration_adjustment_odr(noaa, udaq, 
                                                                                                        noaa_species = voc_spec, 
                                                                                                        udaq_species = mapping[voc_spec],
                                                                                                        time_interval_formaldehyde_used = 15*60)

    # (noaa_species_name, udaq_species_name, 
    # ml_initial_ols, ml_initial_overlap_ols, 
    # udaq_initial_ols, udaq_initial_overlap_ols, 
    # udaq_adjusted_ols, udaq_overlap_adjusted_ols,
    # calibration_eq_ols, rmse_original_ols, 
    # rmse_calibrated_ols, rmse_percent_improvement_ols, merged_voc_ols) = calibration_adjustment_ols(noaa, udaq, 
    #                                                                     noaa_species = voc_spec, 
    #                                                                     udaq_species = mapping[voc_spec],
    #                                                                     time_interval_formaldehyde_used = 15*60)

    list_calibration_adj_total_odr.append({
        'NOAA_Species_Name':noaa_species_name, 
        'UDAQ_Species_Name':udaq_species_name,
        'Calibration_Equation': calibration_eq_odr,
        'RMSE_Initial':rmse_original, 
        'RMSE_Calibrated': rmse_adjusted,
        'Did_RMSE_Improve': did_rmse_improve,
        'RMSE_Percent_Improvement': rmse_percent_improvement,

        'Slope_ML_to_Initial_UDAQ': slope_initial, 
        'Distance_from_1_Slope_ML_to_Initial_UDAQ': initial_slope_distance, 
        'Slope_ML_to_Calibrated_UDAQ': slope_adjusted, 
        'Distance_from_1_Slope_ML_to_Calibrated_UDAQ':adjusted_slope_distance, 
        'Distance_from_1_Slopes_Comparison': slope_distance_difference, 
        'Did_Slope_Improve': did_slope_improve,

        'Intercept_ML_to_Initial_UDAQ': intercept_initial, 
        'Intercept_ML_to_Calibrated_UDAQ':intercept_adjusted,

        'R_ML_to_Initial_UDAQ': r_initial, 
        'R_Squared_ML_to_Initial_UDAQ': r2_initial, 
        'R_ML_to_Calibrated_UDAQ': r_adjusted, 
        'R_Squared_ML_to_Calibrated_UDAQ': r2_adjusted
    })

    # list_calibration_adj_total_ols.append({
    #     'NOAA_Species_Name':noaa_species_name, 
    #     'UDAQ_Species_Name':udaq_species_name,
    #     'Calibration_Equation': calibration_eq_ols,
    #     'RMSE_Original':rmse_original_ols, 
    #     'RMSE_Calibrated': rmse_calibrated_ols,
    #     'RMSE_Improvement': rmse_percent_improvement_ols
    # })

    if 'ML_NoVar' in voc_spec: 
        data_merged_dict_odr['NOAA_'+str(noaa_species_name)+'_Initial'] = pd.Series(np.nan, index=udaq[udaq_species_name].index)
        data_merged_dict_odr['UDAQ_'+str(udaq_species_name)+'_Initial'] = udaq[udaq_species_name]
        data_merged_dict_odr['UDAQ_'+str(udaq_species_name)+'_Adjusted'] = udaq_adjusted_odr
        data_merged_dict_odr['Filled_'+str(noaa_species_name)+'_ML_with_UDAQ_Adjusted'] = udaq[udaq_species_name]

        data_merged_overlap_dict_odr['NOAA_'+str(noaa_species_name)+'_Initial'] = pd.Series(np.nan, index=udaq_initial_overlap_odr.index)
        data_merged_overlap_dict_odr['UDAQ_'+str(udaq_species_name)+'_Initial'] = udaq_initial_overlap_odr
        data_merged_overlap_dict_odr['UDAQ_'+str(udaq_species_name)+'_Adjusted'] = udaq_overlap_adjusted_odr
    
        # data_merged_dict_ols['NOAA_'+str(noaa_species_name)+'_Initial'] = pd.Series(np.nan, index=udaq[udaq_species_name].index)
        # data_merged_dict_ols['UDAQ_'+str(udaq_species_name)+'_Initial'] = udaq[udaq_species_name]
        # data_merged_dict_ols['UDAQ_'+str(udaq_species_name)+'_Adjusted'] = udaq_adjusted_ols
        # data_merged_dict_ols['Filled_'+str(udaq_species_name)+'_ML_with_UDAQ_Adjusted'] = udaq[udaq_species_name]
        
    elif 'UDAQ_NoVar' in mapping[voc_spec]: 
        data_merged_dict_odr['NOAA_'+str(noaa_species_name)+'_Initial'] = noaa[noaa_species_name]
        data_merged_dict_odr['UDAQ_'+str(udaq_species_name)+'_Initial'] = pd.Series(np.nan, index=noaa[noaa_species_name].index)
        data_merged_dict_odr['UDAQ_'+str(udaq_species_name)+'_Adjusted'] = udaq_adjusted_odr
        data_merged_dict_odr['Filled_'+str(noaa_species_name)+'_ML_with_UDAQ_Adjusted'] = noaa[noaa_species_name]

        data_merged_overlap_dict_odr['NOAA_'+str(noaa_species_name)+'_Initial'] = ml_initial_overlap_odr
        data_merged_overlap_dict_odr['UDAQ_'+str(udaq_species_name)+'_Initial'] = pd.Series(np.nan, index=noaa[noaa_species_name].index)
        data_merged_overlap_dict_odr['UDAQ_'+str(udaq_species_name)+'_Adjusted'] = udaq_overlap_adjusted_odr
        data_merged_overlap_dict_odr['Filled_'+str(noaa_species_name)+'_ML_with_UDAQ_Adjusted'] = ml_initial_overlap_odr

        # data_merged_dict_ols['NOAA_'+str(noaa_species_name)+'_Initial'] = noaa[noaa_species_name]
        # data_merged_dict_ols['UDAQ_'+str(udaq_species_name)+'_Initial'] = pd.Series(np.nan, index=noaa[noaa_species_name].index)
        # data_merged_dict_ols['UDAQ_'+str(udaq_species_name)+'_Adjusted'] = udaq_adjusted_ols
        # data_merged_dict_ols['Filled_'+str(noaa_species_name)+'_ML_with_UDAQ_Adjusted'] = noaa[noaa_species_name]

    elif 'Formaldehyde' in mapping[voc_spec]: 
        # (ml_initial_formaldehyde_odr, ml_initial_formaldehyde_overlap_odr, udaq_initial_formaldehyde_odr, 
        # udaq_initial_formaldehyde_overlap_odr, udaq_adjusted_formaldehyde_odr, udaq_overlap_adjusted_formaldehyde_odr,
        # calibration_eq_formaldehyde_odr, rmse_original_formaldehyde_odr, 
        # rmse_calibrated_formaldehyde_odr, rmse_percent_improvement_formaldehyde_odr, merged_formaldehyde_odr) = formaldehyde_only_calibration_adjustment_odr(time_interval_formaldehyde =  15*60)

    
        # data_merged_dict_odr['NOAA_Formaldehyde_Initial'] = ml_initial_formaldehyde_odr
        # data_merged_dict_odr['UDAQ_Formaldehyde_Initial'] = udaq_initial_formaldehyde_odr
        # data_merged_dict_odr['Filled_Formaldehyde_ML_with_UDAQ_Adjusted'] = merged_voc_odr

        data_merged_dict_odr['NOAA_Formaldehyde_Initial'] = ml_initial_odr
        data_merged_dict_odr['UDAQ_Formaldehyde_Initial'] = udaq_initial_odr
        data_merged_dict_odr['UDAQ_Formaldehyde_Adjusted'] = udaq_adjusted_odr
        data_merged_dict_odr['Filled_Formaldehyde_ML_with_UDAQ_Adjusted'] = merged_voc

        data_merged_overlap_dict_odr['NOAA_Formaldehyde_Initial'] = ml_initial_overlap_odr
        data_merged_overlap_dict_odr['UDAQ_Formaldehyde_Initial'] = udaq_initial_overlap_odr
        data_merged_overlap_dict_odr['UDAQ_Formaldehyde_Adjusted'] = udaq_overlap_adjusted_odr

        
        #Comment out if you don't want plots
        # plot_formaldehyde_only_calibration_adjustment_odr(ml_initial_formaldehyde_odr = ml_initial_odr, 
        #                                                 ml_initial_formaldehyde_overlap_odr = ml_initial_overlap_odr, 
        #                                                 udaq_initial_formaldehyde_odr = udaq_initial_odr, 
        #                                                 udaq_initial_formaldehyde_overlap_odr = udaq_initial_overlap_odr, 
        #                                                 udaq_adjusted_formaldehyde_odr = udaq_adjusted_odr, 
        #                                                 udaq_overlap_adjusted_formaldehyde_odr = udaq_overlap_adjusted_odr,
        #                                                 calibration_eq_formaldehyde_odr = calibration_eq_odr, 
        #                                                 merged_formaldehyde_odr = merged_voc)
    
        # (ml_initial_formaldehyde_ols, ml_initial_formaldehyde_overlap_ols, udaq_initial_formaldehyde_ols, 
        # udaq_initial_formaldehyde_overlap_ols, udaq_adjusted_formaldehyde_ols, udaq_overlap_adjusted_formaldehyde_ols,
        # calibration_eq_formaldehyde_ols, rmse_original_formaldehyde_ols, 
        # rmse_calibrated_formaldehyde_ols, rmse_percent_improvement_formaldehyde_ols, merged_formaldehyde_ols) = formaldehyde_only_calibration_adjustment_ols(time_interval_formaldehyde =  15*60)
        
        # data_merged_dict_ols['NOAA_Formaldehyde_Initial'] = ml_initial_ols
        # data_merged_dict_ols['UDAQ_Formaldehyde_Initial'] = udaq_initial_ols
        # data_merged_dict_odr['UDAQ_Formaldehyde_Adjusted'] = udaq_adjusted_ols
        # data_merged_dict_ols['Filled_Formaldehyde_ML_with_UDAQ_Adjusted'] = merged_voc_ols

        # #Comment out if you don't want plots
        # plot_formaldehyde_only_calibration_adjustment_ols(ml_initial_formaldehyde_ols = ml_initial_ols, 
        #                                                   ml_initial_formaldehyde_overlap_ols = ml_initial_overlap_ols, 
        #                                                   udaq_initial_formaldehyde_ols = udaq_initial_ols, 
        #                                                   udaq_initial_formaldehyde_overlap_ols = udaq_initial_overlap_ols, 
        #                                                   udaq_adjusted_formaldehyde_ols = udaq_adjusted_ols, 
        #                                                   udaq_overlap_adjusted_formaldehyde_ols = udaq_overlap_adjusted_ols,
        #                                                   calibration_eq_formaldehyde_ols = calibration_eq_ols, 
        #                                                   merged_formaldehyde_ols = merged_voc_ols)

    else: 
        data_merged_dict_odr['NOAA_'+str(noaa_species_name)+'_Initial'] = ml_initial_odr
        data_merged_dict_odr['UDAQ_'+str(udaq_species_name)+'_Initial'] = udaq_initial_odr
        data_merged_dict_odr['UDAQ_'+str(udaq_species_name)+'_Adjusted'] = udaq_adjusted_odr
        data_merged_dict_odr['Filled_'+str(noaa_species_name)+'_ML_with_UDAQ_Adjusted'] = merged_voc

        data_merged_overlap_dict_odr['NOAA_'+str(noaa_species_name)+'_Initial'] = ml_initial_overlap_odr
        data_merged_overlap_dict_odr['UDAQ_'+str(udaq_species_name)+'_Initial'] = udaq_initial_overlap_odr
        data_merged_overlap_dict_odr['UDAQ_'+str(udaq_species_name)+'_Adjusted'] = udaq_overlap_adjusted_odr
    
        # data_merged_dict_ols['NOAA_'+str(noaa_species_name)+'_Initial'] = noaa[noaa_species_name]
        # data_merged_dict_ols['UDAQ_'+str(udaq_species_name)+'_Initial'] = udaq[udaq_species_name]
        # data_merged_dict_ols['UDAQ_'+str(udaq_species_name)+'_Adjusted'] = udaq_adjusted_ols
        # data_merged_dict_ols['Filled_'+str(udaq_species_name)+'_ML_with_UDAQ_Adjusted'] = merged_voc_ols
        
        #For now, the if statement should only catch Styrene
        if np.isscalar(udaq_initial_odr):
            print(udaq_species_name, ' has no UDAQ Data (all NaNs), using ML Data only')

        # else:
            # #Comment out if you don't want plots
            # plot_calibration_adjustment_odr(udaq_species_name, ml_initial_odr, ml_initial_overlap_odr, 
            #                         udaq_initial_odr, udaq_initial_overlap_odr, 
            #                         udaq_adjusted_odr, udaq_overlap_adjusted_odr, 
            #                         calibration_eq_odr, merged_voc)
            
            # #Comment out if you don't want plots
            # plot_calibration_adjustment_ols(udaq_species_name, ml_initial_ols, ml_initial_overlap_ols, 
            #                         udaq_initial_ols, udaq_initial_overlap_ols, 
            #                         udaq_adjusted_ols, udaq_overlap_adjusted_ols, 
            #                         calibration_eq_ols, merged_voc_ols)
    

df_calibration_adj_total_odr = pd.DataFrame(list_calibration_adj_total_odr)
savepath_odr_info = dirpath + '/Merge_scripts/ml_filling_with_udaq_data/calibration_data/odr_cal_adj_voc_info.csv'
df_calibration_adj_total_odr.to_csv(savepath_odr_info)

# df_calibration_adj_total_ols = pd.DataFrame(list_calibration_adj_total_ols)
# savepath_ols_info = dirpath + '/Merge_scripts/calibration_adjustments/calibration_data/ols_cal_adj_voc_info.csv'
# df_calibration_adj_total_ols.to_csv(savepath_ols_info)

df_filled_ML_with_UDAQ_Adjusted_odr = pd.DataFrame(data_merged_dict_odr)
#print('df_filled_ML_with_UDAQ_Adjusted_odr index time zone: ', df_filled_ML_with_UDAQ_Adjusted_odr.index.tz)
savepath_odr_vals = dirpath + '/Merge_scripts/ml_filling_with_udaq_data/calibration_data/odr_ml_udaq_initial_and_filled_ml_with_udaq_adjusted.csv'
df_filled_ML_with_UDAQ_Adjusted_odr.to_csv(savepath_odr_vals)
print('Saved data for ODR fitting applied UDAQ data, then filled gaps in ML data: ', savepath_odr_vals)

df_filled_ML_with_UDAQ_Adjusted_odr_overlap = pd.DataFrame(data_merged_overlap_dict_odr)
#print('df_filled_ML_with_UDAQ_Adjusted_odr index time zone: ', df_filled_ML_with_UDAQ_Adjusted_odr.index.tz)
savepath_odr_vals_overlap = dirpath + '/Merge_scripts/ml_filling_with_udaq_data/calibration_data/odr_ml_udaq_initial_and_filled_ml_with_udaq_adjusted_overlap.csv'
df_filled_ML_with_UDAQ_Adjusted_odr_overlap.to_csv(savepath_odr_vals_overlap)
print('Saved data for ODR fitting applied UDAQ data, then filled gaps in ML data: ', savepath_odr_vals_overlap)

# df_filled_ML_with_UDAQ_Adjusted_ols = pd.DataFrame(data_merged_dict_ols)
# print('df_filled_ML_with_UDAQ_Adjusted_ols index time zone: ', df_filled_ML_with_UDAQ_Adjusted_ols.index.tz)
# savepath_ols_vals = dirpath + '/Merge_scripts/calibration_adjustments/calibration_data/ols_ml_udaq_initial_and_filled_ml_with_udaq_adjusted.csv'
# df_filled_ML_with_UDAQ_Adjusted_ols.to_csv(savepath_ols_vals)
# print('Saved data for OLS fitting applied UDAQ data, then filled gaps in ML data: ', savepath_ols_vals)
#endregion
#####################



#####################
ml_udaq_initial_and_filled_ml_load = dirpath + '/Merge_scripts/ml_filling_with_udaq_data/calibration_data/odr_ml_udaq_initial_and_filled_ml_with_udaq_adjusted.csv'
df_ml_udaq_initial_and_filled_ml = pd.read_csv(ml_udaq_initial_and_filled_ml_load, index_col='time_local', parse_dates=True)

#For some reason, pandas is reading the time_local as UTC so we convert it back to the UTC-6 time zone
df_ml_udaq_initial_and_filled_ml.index = df_ml_udaq_initial_and_filled_ml.index.tz_localize(None)
df_ml_udaq_initial_and_filled_ml.index = df_ml_udaq_initial_and_filled_ml.index.tz_localize('America/Denver')

def plot_merged_data_time_series(df_index, init_col, filled_col, vocname_udaq):
    fig, (ax1, ax2) = plt.subplots(2,1, figsize = (16,8), tight_layout=True)
    xlim_start_jul = pd.to_datetime('2024-07-15 00:00:00').tz_localize('America/Denver')
    xlim_end_jul = pd.to_datetime('2024-07-31 23:45:00').tz_localize('America/Denver')
    xlim_start_aug = pd.to_datetime('2024-08-01 00:00:00').tz_localize('America/Denver')
    xlim_end_aug = pd.to_datetime('2024-08-18 23:45:00').tz_localize('America/Denver')

    #ax1 is the first row of subplot, for July only
    valid_ml_points_initial = ~np.isnan(init_col)
    ax1.plot(df_index[valid_ml_points_initial], init_col[valid_ml_points_initial], linestyle='solid', color = 'm', marker = 'x', label = 'ML Obs', alpha = 0.7)
    valid_points_filled = ~np.isnan(filled_col)
    ax1.plot(df_index[valid_points_filled], filled_col[valid_points_filled], linestyle='solid', color = 'y', marker = '.', label = 'ML Obs filled with UDAQ Adj Obs')

    #Set x ticks
    tz_mdt = noaa.index.tz #this time zone should be in Mountain Daylight Time
    print(tz_mdt)
    ax1.xaxis.set_major_locator(mdates.DayLocator(tz=tz_mdt))
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d', tz=tz_mdt))
    # Minor ticks: every 3 hours
    ax1.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21], tz=tz_mdt))
    # Rotate and format tick labels
    ax1.tick_params(axis='x', which='major')
    ax1.tick_params(axis='x', which='minor', length=3, color='gray')
    
    ax1.set_ylabel(vocname_udaq + ' (ppb)')
    #ax1.set_xlabel('Date')
    ax1.margins(x=0)
    ax1.set_xlim([xlim_start_jul, xlim_end_jul])

    ax1.legend(loc = 'upper right')

    #ax2 is the second row of subplot, for August only
    valid_ml_points_initial = ~np.isnan(init_col)
    ax2.plot(df_index[valid_ml_points_initial], init_col[valid_ml_points_initial], linestyle='solid', color = 'm', marker = 'x', label = 'ML Obs', alpha = 0.7)
    valid_points_filled = ~np.isnan(filled_col)
    ax2.plot(df_index[valid_points_filled], filled_col[valid_points_filled], linestyle='solid', color = 'y', marker = '.', label = 'ML Obs filled with UDAQ Adj Obs')

    #Set x ticks
    ax2.xaxis.set_major_locator(mdates.DayLocator(tz=tz_mdt))
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d', tz=tz_mdt))
    # Minor ticks: every 3 hours
    ax2.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21], tz=tz_mdt))
    # Rotate and format tick labels
    ax2.tick_params(axis='x', which='major')
    ax2.tick_params(axis='x', which='minor', length=3, color='gray')

    ax2.set_ylabel(vocname_udaq + ' (ppb)')
    ax2.set_xlabel('Time (MDT)')
    ax2.margins(x=0)

    ax2.set_xlim([xlim_start_aug, xlim_end_aug])
    ax2.legend(loc = 'upper right')

    #Mark midnight for every day
    midnight_vals = []
    for midnight_idx in range(0,len(noaa.index),96):
        midnight_vals.append(noaa.index[midnight_idx])
    for day_pos in midnight_vals:
        ax1.axvline(day_pos, color = 'black', linestyle = 'dotted')
        ax2.axvline(day_pos, color = 'black', linestyle = 'dotted')

    plt.savefig(dirpath + '/Merge_scripts/calibration_adjustments/filling_ml_with_udaq_plots/filling_ml_with_udaq_hawthorne_'+ str(vocname_udaq) + '_comparison_july_aug_timeseries.png', dpi =300)
    plt.show()
def scatterplot_odr():
    #Plot 1: ML to UDAQ relationship initial
    #Plot scatterplot relationship between VOC and CO
    #Looking for: Clear positive trend, not a vertical cloud, no obvious curvature or thresholds
    co_align, voc_align = ml_co_raw.align(filled_col, join='inner') #include only overlap between VOC and CO
    mask = co_align.notna() & voc_align.notna()
    co_nonans  = filled_col[mask]
    voc_nonans = voc_align[mask]

    plt.figure(figsize=(5,5))
    plt.scatter(co_nonans, voc_nonans, s=10, alpha=0.3)
    plt.xlabel('CO (ppb)')
    plt.ylabel(vocname_udaq + ' (ppb)')
    plt.tight_layout()
    #plt.savefig(dirpath + '/Merge_scripts/calibration_adjustments/co_tracer_plots/scatterplot_co_vs_voc/ml_co_with_'+ str(vocname_udaq) + '_scatterplot.png', dpi =300)
    plt.show()

    #Plot 2: ML to UDAQ relationship after applying ODR fit

def apply_interpolation_to_small_gaps(df_index, init_col, filled_col, vocname):
    gap_records_fulltime = []

    #Identify where real measurements exist
    is_measured_fulltime = filled_col.notna()
    #print('is_measured:', is_measured)

    #Compute time gaps only between valid measurements
    measurement_times_fulltime = df_index[is_measured_fulltime]
    gaps_fulltime =  pd.Series(measurement_times_fulltime) - pd.Series(measurement_times_fulltime).shift(1)

    # Internal gaps
    if len(measurement_times_fulltime) > 1:
        for start_of_gap, end_of_gap in zip(measurement_times_fulltime[:-1], measurement_times_fulltime[1:]):
            gap_size = end_of_gap - start_of_gap
            gap_records_fulltime.append({
                'gap_start': start_of_gap,
                'gap_end': end_of_gap,
                'gap_size': gap_size,
                'is_large_gap': gap_size >  pd.Timedelta(hours=3)
            })

    # Leading gap
    if len(measurement_times_fulltime) >= 1:
        leading_gap = measurement_times_fulltime[0] - df_index.min()
        if leading_gap > pd.Timedelta(0):
            gap_records_fulltime.insert(0, {
                'gap_start': df_index.min(),
                'gap_end': measurement_times_fulltime[0],
                'gap_size': leading_gap,
                'is_large_gap': leading_gap >  pd.Timedelta(hours=3)
            })

    # Trailing gap
    if len(measurement_times_fulltime) >= 1:
        trailing_gap = df_index.max() - measurement_times_fulltime[-1]
        if trailing_gap > pd.Timedelta(0):
            gap_records_fulltime.append({
                'gap_start': measurement_times_fulltime.max(),
                'gap_end': df_index.max(),
                'gap_size': trailing_gap,
                'is_large_gap': trailing_gap >  pd.Timedelta(hours=3)
            })

    gap_df_fulltime = pd.DataFrame(gap_records_fulltime)
    display(gap_df_fulltime)
    
    ds_filled_data = filled_col.copy()

    if gaps_fulltime.empty:
        print(vocname, ' has no measurements.')
        nan_fill_series = pd.Series(index=df_index)
        filled_var = nan_fill_series
        ds_filled_interp_final = filled_var

    else:
        interp_mask = pd.Series(False, index=df_index)

        #Find all small gaps where the gaps are <= 1 hour
        small_gaps = gap_df_fulltime[gap_df_fulltime['gap_size'] <= pd.Timedelta(hours=1)]
        gaps_larger_than_one_hour = gap_df_fulltime[gap_df_fulltime['gap_size'] > pd.Timedelta(hours=1)]

        for _, row in small_gaps.iterrows():
            interp_mask |= (df_index >= row['gap_start']) & (df_index <= row['gap_end'])

        ds_filled_interp = ds_filled_data.interpolate(method="time")

        ds_filled_interp_final = ds_filled_data.copy()
        ds_filled_interp_final.loc[interp_mask] = ds_filled_interp.loc[interp_mask]

        ml_initial_gaps = init_col.isna().sum()
        adjusted_filled_gaps = filled_col.isna().sum()
        interp_gaps = ds_filled_interp_final.isna().sum()

        number_gaps_filled_by_adjusted_udaq_data = ml_initial_gaps - adjusted_filled_gaps
        percent_of_gaps_filled_by_udaq = ((ml_initial_gaps - adjusted_filled_gaps) / ml_initial_gaps) * 100

        number_gaps_filled_by_interp = adjusted_filled_gaps - interp_gaps
        percent_of_gaps_filled_by_interp = ((adjusted_filled_gaps - interp_gaps) / adjusted_filled_gaps) * 100

        percent_of_initial_gaps_filled_by_udaq_then_interp = (((ml_initial_gaps - adjusted_filled_gaps) + (adjusted_filled_gaps - interp_gaps)) / ml_initial_gaps) * 100

        print('# Gaps in initial ML data: ', ml_initial_gaps)
        print('# Gaps in ML Data after filling with calibration adjusted UDAQ data: ', adjusted_filled_gaps)
        print('# Gaps filled by calibration adjusted UDAQ data: ', number_gaps_filled_by_adjusted_udaq_data)
        print('Percentage of gaps filled from initial ML data after filling with calibration adjusted UDAQ data: ', percent_of_gaps_filled_by_udaq)
        print('# Gaps filled by interpolation: ', number_gaps_filled_by_interp)
        print('Percentage of gaps filled from remaining gaps after filling with calibration adjusted UDAQ data from interpolation: ', percent_of_gaps_filled_by_interp)
        print('Percentage of gaps filled from initial ML data by filling with calibration + interpolation: ', percent_of_initial_gaps_filled_by_udaq_then_interp)
    
    print(ds_filled_interp_final)
    return ds_filled_interp_final
    

    # for gap_row in range(0,len(gap_df_fulltime)):
    #     if gap_df_fulltime['gap_size'] <= pd.Timedelta(hours=1):
    #         df_filled_data.sort_index().loc[gap_df_fulltime['gap_start']:gap_df_fulltime['gap_end']].interpolate(method='linear')
    #         print('Interpolated gap from ', gap_df_fulltime['gap_start'], ' to ', gap_df_fulltime['gap_end'])
    #     else:
def apply_interpolation_to_small_gaps_2hrs(df_index, init_col, filled_col, vocname):
    gap_records_fulltime = []

    #Identify where real measurements exist
    is_measured_fulltime = filled_col.notna()
    #print('is_measured:', is_measured)

    #Compute time gaps only between valid measurements
    measurement_times_fulltime = df_index[is_measured_fulltime]
    gaps_fulltime =  pd.Series(measurement_times_fulltime) - pd.Series(measurement_times_fulltime).shift(1)

    # Internal gaps
    if len(measurement_times_fulltime) > 1:
        for start_of_gap, end_of_gap in zip(measurement_times_fulltime[:-1], measurement_times_fulltime[1:]):
            gap_size = end_of_gap - start_of_gap
            gap_records_fulltime.append({
                'gap_start': start_of_gap,
                'gap_end': end_of_gap,
                'gap_size': gap_size,
                'is_large_gap': gap_size >  pd.Timedelta(hours=3)
            })

    # Leading gap
    if len(measurement_times_fulltime) >= 1:
        leading_gap = measurement_times_fulltime[0] - df_index.min()
        if leading_gap > pd.Timedelta(0):
            gap_records_fulltime.insert(0, {
                'gap_start': df_index.min(),
                'gap_end': measurement_times_fulltime[0],
                'gap_size': leading_gap,
                'is_large_gap': leading_gap >  pd.Timedelta(hours=3)
            })

    # Trailing gap
    if len(measurement_times_fulltime) >= 1:
        trailing_gap = df_index.max() - measurement_times_fulltime[-1]
        if trailing_gap > pd.Timedelta(0):
            gap_records_fulltime.append({
                'gap_start': measurement_times_fulltime.max(),
                'gap_end': df_index.max(),
                'gap_size': trailing_gap,
                'is_large_gap': trailing_gap >  pd.Timedelta(hours=3)
            })

    gap_df_fulltime = pd.DataFrame(gap_records_fulltime)
    display(gap_df_fulltime)
    
    ds_filled_data = filled_col.copy()

    if gaps_fulltime.empty:
        print(vocname, ' has no measurements.')
        nan_fill_series = pd.Series(index=df_index)
        filled_var = nan_fill_series
        ds_filled_interp_final = filled_var

    else:
        interp_mask = pd.Series(False, index=df_index)

        #Find all small gaps where the gaps are <= 1 hour
        small_gaps = gap_df_fulltime[gap_df_fulltime['gap_size'] <= pd.Timedelta(hours=2)]
        #gaps_larger_than_one_hour = gap_df_fulltime[gap_df_fulltime['gap_size'] > pd.Timedelta(hours=1)]

        for _, row in small_gaps.iterrows():
            interp_mask |= (df_index >= row['gap_start']) & (df_index <= row['gap_end'])

        ds_filled_interp = ds_filled_data.interpolate(method="time")

        ds_filled_interp_final = ds_filled_data.copy()
        ds_filled_interp_final.loc[interp_mask] = ds_filled_interp.loc[interp_mask]

        ml_initial_gaps = init_col.isna().sum()
        adjusted_filled_gaps = filled_col.isna().sum()
        interp_gaps = ds_filled_interp_final.isna().sum()

        number_gaps_filled_by_adjusted_udaq_data = ml_initial_gaps - adjusted_filled_gaps
        percent_of_gaps_filled_by_udaq = ((ml_initial_gaps - adjusted_filled_gaps) / ml_initial_gaps) * 100

        number_gaps_filled_by_interp = adjusted_filled_gaps - interp_gaps
        percent_of_gaps_filled_by_interp = ((adjusted_filled_gaps - interp_gaps) / adjusted_filled_gaps) * 100

        percent_of_initial_gaps_filled_by_udaq_then_interp = (((ml_initial_gaps - adjusted_filled_gaps) + (adjusted_filled_gaps - interp_gaps)) / ml_initial_gaps) * 100

        print('# Gaps in initial ML data: ', ml_initial_gaps)
        print('# Gaps in ML Data after filling with calibration adjusted UDAQ data: ', adjusted_filled_gaps)
        print('# Gaps filled by calibration adjusted UDAQ data: ', number_gaps_filled_by_adjusted_udaq_data)
        print('Percentage of gaps filled from initial ML data after filling with calibration adjusted UDAQ data: ', percent_of_gaps_filled_by_udaq)
        print('# Gaps filled by interpolation: ', number_gaps_filled_by_interp)
        print('Percentage of gaps filled from remaining gaps after filling with calibration adjusted UDAQ data from interpolation: ', percent_of_gaps_filled_by_interp)
        print('Percentage of gaps filled from initial ML data by filling with calibration + interpolation: ', percent_of_initial_gaps_filled_by_udaq_then_interp)
    
    print(ds_filled_interp_final)
    return ds_filled_interp_final
def plot_interpolation_time_series(df_index, init_col, filled_col, interpolated_index, interpolated_col, vocname):
    fig, (ax1, ax2) = plt.subplots(2,1, figsize = (16,8), tight_layout=True)
    xlim_start_jul = pd.to_datetime('2024-07-15 00:00:00').tz_localize('America/Denver')
    xlim_end_jul = pd.to_datetime('2024-07-31 23:45:00').tz_localize('America/Denver')
    xlim_start_aug = pd.to_datetime('2024-08-01 00:00:00').tz_localize('America/Denver')
    xlim_end_aug = pd.to_datetime('2024-08-18 23:45:00').tz_localize('America/Denver')

    #ax1 is the first row of subplot, for July only
    valid_ml_points_initial = ~np.isnan(init_col)
    ax1.plot(df_index[valid_ml_points_initial], init_col[valid_ml_points_initial], linestyle='solid', color = 'm', marker = 's', label = 'ML Obs', alpha = 0.7)
    valid_points_filled = ~np.isnan(filled_col)
    ax1.plot(df_index[valid_points_filled], filled_col[valid_points_filled], linestyle='solid', color = 'y', marker = 'o', label = 'ML Obs filled with UDAQ Adj Obs')
    valid_points_interpolated = ~np.isnan(interpolated_col)
    ax1.plot(interpolated_index[valid_points_filled], interpolated_col[valid_points_filled], linestyle='solid', color = 'g', marker = '+', label = 'ML Obs filled with UDAQ Adj Obs and Interpolation', alpha = 0.5)

    #Set x ticks
    tz_mdt = df_index.tz #this time zone should be in Mountain Daylight Time
    ax1.xaxis.set_major_locator(mdates.DayLocator(tz=tz_mdt))
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d', tz=tz_mdt))
    # Minor ticks: every 3 hours
    ax1.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21], tz=tz_mdt))
    # Rotate and format tick labels
    ax1.tick_params(axis='x', which='major')
    ax1.tick_params(axis='x', which='minor', length=3, color='gray')
    #ax.grid(True, which='both')
    
    ax1.set_ylabel(vocname + ' (ppb)')
    #ax1.set_xlabel('Date')
    ax1.margins(x=0)
    ax1.set_xlim([xlim_start_jul, xlim_end_jul])

    ax1.legend(loc = 'upper right')

    #ax2 is the second row of subplot, for August only
    
    valid_ml_points_initial = ~np.isnan(init_col)
    ax2.plot(df_index[valid_ml_points_initial], init_col[valid_ml_points_initial], linestyle='solid', color = 'm', marker = 's', label = 'ML Obs', alpha = 0.7)
    valid_points_filled = ~np.isnan(filled_col)
    ax2.plot(df_index[valid_points_filled], filled_col[valid_points_filled], linestyle='solid', color = 'y', marker = 'o', label = 'ML Obs filled with UDAQ Adj Obs')
    valid_points_interpolated = ~np.isnan(interpolated_col)
    ax2.plot(interpolated_index[valid_points_interpolated], interpolated_col[valid_points_interpolated], linestyle='solid', color = 'g', marker = '+', label = 'ML Obs filled with UDAQ Adj Obs and Interpolation', alpha = 0.5)

    #Set x ticks
    ax2.xaxis.set_major_locator(mdates.DayLocator())
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
    # Minor ticks: every 3 hours
    ax2.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21]))
    # Rotate and format tick labels
    ax2.tick_params(axis='x', which='major')
    ax2.tick_params(axis='x', which='minor', length=3, color='gray')

    ax2.set_ylabel(vocname + ' (ppb)')
    ax2.set_xlabel('Time (MDT)')
    ax2.margins(x=0)

    ax2.set_xlim([pd.to_datetime('2024-08-01 00:00:00'), pd.to_datetime('2024-08-18 23:00:00')])
    ax2.legend(loc = 'upper right')

    #Mark midnight for every day
    midnight_vals = []
    for midnight_idx in range(0,len(noaa.index),96):
        midnight_vals.append(noaa.index[midnight_idx])
    for day_pos in midnight_vals:
        ax1.axvline(day_pos, color = 'black', linestyle = 'dotted')
        ax2.axvline(day_pos, color = 'black', linestyle = 'dotted')

    plt.savefig(dirpath + '/Merge_scripts/calibration_adjustments/interpolation_plots/interpolation_after_filling_ml_data_with_udaq_hawthorne_'+ str(vocname) + '_comparison_july_aug_timeseries.png', dpi =300)
    plt.show()
def plot_interpolation_2hr_time_series(df_index, init_col, filled_col, interpolated_index, interpolated_col, vocname):
    fig, (ax1, ax2) = plt.subplots(2,1, figsize = (16,8), tight_layout=True)
    xlim_start_jul = pd.to_datetime('2024-07-15 00:00:00').tz_localize('America/Denver')
    xlim_end_jul = pd.to_datetime('2024-07-31 23:45:00').tz_localize('America/Denver')
    xlim_start_aug = pd.to_datetime('2024-08-01 00:00:00').tz_localize('America/Denver')
    xlim_end_aug = pd.to_datetime('2024-08-18 23:45:00').tz_localize('America/Denver')

    #ax1 is the first row of subplot, for July only
    valid_ml_points_initial = ~np.isnan(init_col)
    ax1.plot(df_index[valid_ml_points_initial], init_col[valid_ml_points_initial], linestyle='solid', color = 'm', marker = 's', label = 'ML Obs', alpha = 0.7)
    valid_points_filled = ~np.isnan(filled_col)
    ax1.plot(df_index[valid_points_filled], filled_col[valid_points_filled], linestyle='solid', color = 'y', marker = 'o', label = 'ML Obs filled with UDAQ Adj Obs')
    valid_points_interpolated = ~np.isnan(interpolated_col)
    ax1.plot(interpolated_index[valid_points_filled], interpolated_col[valid_points_filled], linestyle='solid', color = 'g', marker = '+', label = 'ML Obs filled with UDAQ Adj Obs and Interpolation', alpha = 0.5)

    #Set x ticks
    tz_mdt = df_index.tz #this time zone should be in Mountain Daylight Time
    print(tz_mdt)
    ax1.xaxis.set_major_locator(mdates.DayLocator(tz=tz_mdt))
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d', tz=tz_mdt))
    # Minor ticks: every 3 hours
    ax1.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21], tz=tz_mdt))
    # Rotate and format tick labels
    ax1.tick_params(axis='x', which='major')
    ax1.tick_params(axis='x', which='minor', length=3, color='gray')
    #ax.grid(True, which='both')
    
    ax1.set_ylabel(vocname + ' (ppb)')
    #ax1.set_xlabel('Date')
    ax1.margins(x=0)
    ax1.set_xlim([xlim_start_jul, xlim_end_jul])

    ax1.legend(loc = 'upper right')

    #ax2 is the second row of subplot, for August only
    
    valid_ml_points_initial = ~np.isnan(init_col)
    ax2.plot(df_index[valid_ml_points_initial], init_col[valid_ml_points_initial], linestyle='solid', color = 'm', marker = 's', label = 'ML Obs', alpha = 0.7)
    valid_points_filled = ~np.isnan(filled_col)
    ax2.plot(df_index[valid_points_filled], filled_col[valid_points_filled], linestyle='solid', color = 'y', marker = 'o', label = 'ML Obs filled with UDAQ Adj Obs')
    valid_points_interpolated = ~np.isnan(interpolated_col)
    ax2.plot(interpolated_index[valid_points_interpolated], interpolated_col[valid_points_interpolated], linestyle='solid', color = 'g', marker = '+', label = 'ML Obs filled with UDAQ Adj Obs and Interpolation', alpha = 0.5)

    #Set x ticks
    ax2.xaxis.set_major_locator(mdates.DayLocator(tz=tz_mdt))
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d', tz=tz_mdt))
    # Minor ticks: every 3 hours
    ax2.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21], tz=tz_mdt))
    # Rotate and format tick labels
    ax2.tick_params(axis='x', which='major')
    ax2.tick_params(axis='x', which='minor', length=3, color='gray')

    ax2.set_ylabel(vocname + ' (ppb)')
    ax2.set_xlabel('Time (MDT)')
    ax2.margins(x=0)

    ax2.set_xlim([xlim_start_aug, xlim_end_aug])
    ax2.legend(loc = 'upper right')

    #Mark midnight for every day
    midnight_vals = []
    for midnight_idx in range(0,len(noaa.index),96):
        midnight_vals.append(noaa.index[midnight_idx])
    for day_pos in midnight_vals:
        ax1.axvline(day_pos, color = 'black', linestyle = 'dotted')
        ax2.axvline(day_pos, color = 'black', linestyle = 'dotted')

    plt.savefig(dirpath + '/Merge_scripts/calibration_adjustments/interpolation_2hr_plots/interpolation_2hrs_after_filling_ml_data_with_udaq_hawthorne_'+ str(vocname) + '_comparison_july_aug_timeseries.png', dpi =300)
    plt.show()


def plot_carbon_monoxide_tracer_timeseries(df_index, filled_col, vocname_udaq):
    if filled_col.isna().all():
        print(vocname_udaq, ' has no measurements.')
    else: 
        fig, (ax1, ax2) = plt.subplots(2,1, figsize = (16,8), tight_layout=True)
        xlim_start_jul = pd.to_datetime('2024-07-15 00:00:00').tz_localize('America/Denver')
        xlim_end_jul = pd.to_datetime('2024-07-31 23:45:00').tz_localize('America/Denver')
        xlim_start_aug = pd.to_datetime('2024-08-01 00:00:00').tz_localize('America/Denver')
        xlim_end_aug = pd.to_datetime('2024-08-18 23:45:00').tz_localize('America/Denver')

        #ax1 is the first row of subplot, for July only
        valid_ml_points_after_filling = ~np.isnan(filled_col)
        ax1.plot(df_index[valid_ml_points_after_filling], filled_col[valid_ml_points_after_filling], linestyle='solid', color = 'm', marker = 'x', label = f'ML Filled with Adjusted UDAQ {vocname_udaq}', alpha = 0.7)
        ax1_co = ax1.twinx()
        valid_points_co = ~np.isnan(ml_co_raw)
        ax1_co.plot(df_index[valid_points_co], ml_co_raw[valid_points_co], linestyle='solid', color = 'y', marker = '.', label = 'ML CO Obs')

        #Set x ticks
        tz_mdt = df_index.tz #this time zone should be in Mountain Daylight Time
        print(tz_mdt)
        ax1.xaxis.set_major_locator(mdates.DayLocator(tz=tz_mdt))
        ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d', tz=tz_mdt))
        # Minor ticks: every 3 hours
        ax1.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21], tz=tz_mdt))
        # Rotate and format tick labels
        ax1.tick_params(axis='x', which='major')
        ax1.tick_params(axis='x', which='minor', length=3, color='gray')
        #ax.grid(True, which='both')

        ax1.set_ylabel(vocname_udaq + ' (ppb)')
        ax1_co.set_ylabel('CO (ppb)')
        ax1.margins(x=0)
        ax1.set_xlim([xlim_start_jul, xlim_end_jul])

        lines1, labels1 = ax1.get_legend_handles_labels()
        lines1_co, labels1_co = ax1_co.get_legend_handles_labels()
        ax1.legend(lines1 + lines1_co, labels1 + labels1_co, loc='upper right')

        #ax2 is the second row of subplot, for August only
        ax2.plot(df_index[valid_ml_points_after_filling], filled_col[valid_ml_points_after_filling], linestyle='solid', color = 'm', marker = 'x', label = f'ML Filled with Adjusted UDAQ {vocname_udaq}', alpha = 0.7)
        ax2_co = ax2.twinx()
        ax2_co.plot(df_index[valid_points_co], ml_co_raw[valid_points_co], linestyle='solid', color = 'y', marker = '.', label = 'ML CO Obs')
       
        #Set x ticks
        ax2.xaxis.set_major_locator(mdates.DayLocator(tz=tz_mdt))
        ax2.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d', tz=tz_mdt))
        # Minor ticks: every 3 hours
        ax2.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21], tz=tz_mdt))
        # Rotate and format tick labels
        ax2.tick_params(axis='x', which='major')
        ax2.tick_params(axis='x', which='minor', length=3, color='gray')

        ax2.set_ylabel(vocname_udaq + ' (ppb)')
        ax2_co.set_ylabel('CO (ppb)')
        ax2.set_xlabel('Time (MDT)')
        ax2.margins(x=0)

        ax2.set_xlim([xlim_start_aug, xlim_end_aug])

        lines2, labels2 = ax2.get_legend_handles_labels()
        lines2_co, labels2_co = ax2_co.get_legend_handles_labels()
        ax2.legend(lines2 + lines2_co, labels2 + labels2_co, loc='upper right')

        #Mark midnight for every day
        midnight_vals = []
       
        for midnight_idx in range(0,len(noaa.index),96):
            midnight_vals.append(noaa.index[midnight_idx])
        for day_pos in midnight_vals:
            ax1.axvline(day_pos, color = 'black', linestyle = 'dotted')
            ax2.axvline(day_pos, color = 'black', linestyle = 'dotted')

        plt.savefig(dirpath + '/Merge_scripts/calibration_adjustments/co_tracer_plots/time_series_ml_co_vs_voc_filled/compare_ml_co_with_'+ str(vocname_udaq) + '_comparison_july_aug_timeseries.png', dpi =300)
        plt.show()
def scatterplot_co_and_voc(filled_col, vocname_udaq):
    #Plot scatterplot relationship between VOC and CO
    #Looking for: Clear positive trend, not a vertical cloud, no obvious curvature or thresholds
    co_align, voc_align = ml_co_raw.align(filled_col, join='inner') #include only overlap between VOC and CO
    mask = co_align.notna() & voc_align.notna()
    co_nonans  = co_align[mask]
    voc_nonans = voc_align[mask]

    plt.figure(figsize=(5,5))
    plt.scatter(co_nonans, voc_nonans, s=10, alpha=0.3)
    plt.xlabel('CO (ppb)')
    plt.ylabel(vocname_udaq + ' (ppb)')
    plt.tight_layout()
    #plt.savefig(dirpath + '/Merge_scripts/calibration_adjustments/co_tracer_plots/scatterplot_co_vs_voc/ml_co_with_'+ str(vocname_udaq) + '_scatterplot.png', dpi =300)
    plt.show()
    return voc_nonans, co_nonans
def co_tracer_voc_regression_ols(idx, filled_col, vocname_udaq):
    #Get clean version of only overlap between CO and VOC that are not nans
    co_align, voc_align = ml_co_raw.align(filled_col, join='inner') 
    mask =  co_align.notna() & voc_align.notna()
    co_nonans  = co_align[mask]
    voc_nonans = voc_align[mask]

    #Make a DataFrame of the clean data
    df_tracer_clean = pd.DataFrame({'CO': co_nonans, 'VOC': voc_nonans})
    #Perform fitting over the clean data
    x_co_nonans = sm.add_constant(df_tracer_clean['CO'])
    y_voc_nonans = df_tracer_clean['VOC']
    tracer_model_fit = sm.OLS(y_voc_nonans, x_co_nonans).fit()

    #Save the slope, y-intercept, and R squared value of our linear fit
    print(f"Slope = {tracer_model_fit.params['const']:.3f}, " f"Intercept = {tracer_model_fit.params['CO']:.3f}, " f"R$^2$ = {tracer_model_fit.rsquared:.3f}")
    tracer_model_fit_slope = tracer_model_fit.params['const']
    tracer_model_fit_intercept = tracer_model_fit.params['CO']
    tracer_model_fit_rsquared = tracer_model_fit.rsquared

    #Apply the regression to all the points of the dataset
    co_15min_constant = sm.add_constant(ml_co_raw)
    voc_reconstructed = tracer_model_fit.predict(co_15min_constant)
    
    #See if we already have data saved. If so, we don't need to add anything new.
    outputs = ['Slope_tracer_model_fit', 'Intercept_tracer_model_fit', 'Rsquared_tracer_model_fit']
    # if df_co_and_voc_results.loc[idx, outputs].notna().all():
    #     print('We already have some outputs for this species.')
    #     return df_co_and_voc_results, voc_reconstructed
    
    #Save the slope, y-intercept, and R squared value of our linear fit so that we can write that into the CSV file when we call the function
    df_co_and_voc_results.loc[idx, 'Slope_tracer_model_fit'] = tracer_model_fit_slope
    df_co_and_voc_results.loc[idx, 'Intercept_tracer_model_fit'] = tracer_model_fit_intercept
    df_co_and_voc_results.loc[idx, 'Rsquared_tracer_model_fit'] = tracer_model_fit_rsquared

    print('voc_reconstructed: ', voc_reconstructed)
    return df_co_and_voc_results, voc_reconstructed
def co_tracer_voc_regression_odr(idx, filled_col, vocname_udaq):
    #Get clean version of only overlap between CO and VOC that are not nans
    co_align, voc_align = ml_co_raw.align(filled_col, join='inner') 
    mask =  co_align.notna() & voc_align.notna()
    co_nonans  = co_align[mask]
    voc_nonans = voc_align[mask]

    #Make a DataFrame of the clean data
    df_tracer_clean = pd.DataFrame({'CO': co_nonans, 'VOC': voc_nonans})
   
   #From the ODRpack documentation
    # Define the function we want to fit against
    voc_overlap_np_arr = df_tracer_clean['VOC'].to_numpy()
    #xdata_2d = np.array([udaq_overlap_voc_np_arr])
    co_overlap_np_arr = df_tracer_clean['CO'].to_numpy()
    #ydata_2d = np.array([ml_overlap_voc_np_arr])

    print('xdata_2d shape: ', np.shape(voc_overlap_np_arr))
    print('ydata shape: ', np.shape(co_overlap_np_arr))

    def linear_model(x, beta):
        return beta[0] + beta[1] * x
    
    #beta0 is an initial guess
    odr_fit_result = odr_fit(linear_model, voc_overlap_np_arr, co_overlap_np_arr, [0.0, 1.0])

    odr_intercept, odr_slope = odr_fit_result.beta

    print('odr_fit_result: ', odr_fit_result)
    print('odr_intercept: ', odr_intercept)
    print('odr_slope: ', odr_slope)

    calibration_eq = f"Calibration equation: O_hat = {odr_intercept:.4f} + {odr_slope:.4f} * M"
    print(calibration_eq)

    voc_vals_adjusted = odr_intercept + odr_slope * df_tracer_clean['VOC']
    #overlap between initial and fitted data
    voc_overlap_adjusted_with_initial = voc_vals_adjusted.loc[voc_vals_adjusted.index.intersection(df_tracer_clean['VOC'].index)]

    #Save the slope, y-intercept, and R squared value of our linear fit
    print(f"Slope = {odr_slope:.3f}, " f"Intercept = {odr_intercept:.3f}, " f"R$^2$ = {tracer_model_fit.rsquared:.3f}")
    model_fit_slope = odr_slope
    model_fit_intercept = odr_intercept
    model_fit_rsquared = 

    
    #See if we already have data saved. If so, we don't need to add anything new.
    outputs = ['Slope_tracer_model_fit', 'Intercept_tracer_model_fit', 'Rsquared_tracer_model_fit']
    # if df_co_and_voc_results.loc[idx, outputs].notna().all():
    #     print('We already have some outputs for this species.')
    #     return df_co_and_voc_results, voc_reconstructed
    
    #Save the slope, y-intercept, and R squared value of our linear fit so that we can write that into the CSV file when we call the function
    df_co_and_voc_results.loc[idx, 'Slope_tracer_model_fit'] = model_fit_slope
    df_co_and_voc_results.loc[idx, 'Intercept_tracer_model_fit'] = model_fit_intercept
    df_co_and_voc_results.loc[idx, 'Rsquared_tracer_model_fit'] = model_fit_rsquared

    print('voc_reconstructed: ', voc_reconstructed)
    return df_co_and_voc_results, voc_reconstructed
def plot_voc_reconstructed(voc_reconstructed, df_index, filled_col, vocname_udaq):
    fig, (ax1, ax2) = plt.subplots(2,1, figsize = (16,8), tight_layout=True)
    xlim_start_jul = pd.to_datetime('2024-07-15 00:00:00').tz_localize('America/Denver')
    xlim_end_jul = pd.to_datetime('2024-07-31 23:45:00').tz_localize('America/Denver')
    xlim_start_aug = pd.to_datetime('2024-08-01 00:00:00').tz_localize('America/Denver')
    xlim_end_aug = pd.to_datetime('2024-08-18 23:45:00').tz_localize('America/Denver')

    #ax1 is the first row of subplot, for July only
    valid_ml_points_after_filling = ~np.isnan(filled_col)
    ax1.plot(df_index[valid_ml_points_after_filling], filled_col[valid_ml_points_after_filling], linestyle='solid', color = 'm', marker = 'x', label = f'ML Filled with Adjusted UDAQ', alpha = 0.7)
    ax1.plot(voc_reconstructed.index, voc_reconstructed, linestyle='solid', color = 'y', marker = '.', label = 'Fitted with CO')

    #Set x ticks
    tz_mdt = df_index.tz #this time zone should be in Mountain Daylight Time
    print(tz_mdt)
    ax1.xaxis.set_major_locator(mdates.DayLocator(tz=tz_mdt))
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d', tz=tz_mdt))
    # Minor ticks: every 3 hours
    ax1.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21], tz=tz_mdt))
    # Rotate and format tick labels
    ax1.tick_params(axis='x', which='major')
    ax1.tick_params(axis='x', which='minor', length=3, color='gray')
    #ax.grid(True, which='both')
    
    ax1.set_ylabel(vocname_udaq + ' (ppb)')
    #ax1.set_xlabel('Date')
    ax1.margins(x=0)
    ax1.set_xlim([xlim_start_jul, xlim_end_jul])

    lines1, labels1 = ax1.get_legend_handles_labels()
    ax1.legend(lines1, labels1, loc='upper right')

    #ax2 is the second row of subplot, for August only
    ax2.plot(df_index[valid_ml_points_after_filling], filled_col[valid_ml_points_after_filling], linestyle='solid', color = 'm', marker = 'x', label = f'ML Filled with Adjusted UDAQ', alpha = 0.7)
    ax2.plot(voc_reconstructed.index, voc_reconstructed, linestyle='solid', color = 'y', marker = '.', label = 'Fitted with CO')

    #Set x ticks
    ax2.xaxis.set_major_locator(mdates.DayLocator(tz=tz_mdt))
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d', tz=tz_mdt))
    # Minor ticks: every 3 hours
    ax2.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21], tz=tz_mdt))
    # Rotate and format tick labels
    ax2.tick_params(axis='x', which='major')
    ax2.tick_params(axis='x', which='minor', length=3, color='gray')

    ax2.set_ylabel(vocname_udaq + ' (ppb)')
    ax2.set_xlabel('Time (MDT)')
    ax2.margins(x=0)

    ax2.set_xlim([xlim_start_aug, xlim_end_aug])
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines2, labels2, loc='upper right')

    #Mark midnight for every day
    midnight_vals = []
    for midnight_idx in range(0,len(noaa.index),96):
        midnight_vals.append(noaa.index[midnight_idx])
    for day_pos in midnight_vals:
        ax1.axvline(day_pos, color = 'black', linestyle = 'dotted')
        ax2.axvline(day_pos, color = 'black', linestyle = 'dotted')

    plt.savefig(dirpath + '/Merge_scripts/calibration_adjustments/co_tracer_plots/time_series_voc_adjusted_vs_fitted/compare_adjusted_voc_with_'+ str(vocname_udaq) + '_comparison_july_aug_timeseries.png', dpi =300)
    plt.show()


reverse_mapping = {v: k for k, v in mapping.items()}
co_and_voc_results = []
calculated_fields = ['Slope_tracer_model_fit', 'Intercept_tracer_model_fit', 'Rsquared_tracer_model_fit']
input_required_fields = ['InputFields_placeholder']
all_fields = ['ML_VOC_Name'] + ['UDAQ_VOC_Name'] + calculated_fields + input_required_fields
tracer_csv_savepath = dirpath + '/Merge_scripts/calibration_adjustments/co_tracer_info.csv'


if os.path.exists(tracer_csv_savepath):
    df_co_and_voc_results = pd.read_csv(tracer_csv_savepath)
else:
    df_co_and_voc_results = pd.DataFrame(columns=all_fields)

# # #############################

def get_row_index(df_co_and_voc_results, vocname_ml, vocname_udaq):
    print('df_co_and_voc_results: ', df_co_and_voc_results)
    print("df_co_and_voc_results['ML_VOC_Name']: ", df_co_and_voc_results['ML_VOC_Name'])
    print('vocname_ml: ', vocname_ml)

    matches = df_co_and_voc_results.index[df_co_and_voc_results['ML_VOC_Name'] == vocname_ml]
    if len(matches) > 0:
        return matches[0]
    
    # create new row
    row = {k: pd.NA for k in all_fields}
    row['ML_VOC_Name'] = vocname_ml
    row['UDAQ_VOC_Name'] = vocname_udaq

    #df_co_and_voc_results.loc[len(df_co_and_voc_results)] = row
    return len(df_co_and_voc_results) - 1

for col in df_ml_udaq_initial_and_filled_ml.columns:
    #NOAA_ML_NoVar_0_Initial,UDAQ_nDodecane_Initial,UDAQ_nDodecane_Adjusted,Filled_nDodecane_ML_with_UDAQ_Adjusted
# 'ML_NoVar' noaa udaq udaq udaq
# UDAQ_NoVar noaa udaq udaq noaa
# Formaldehyde
# else noaa udaq udaq udaq

    if col.startswith('Filled_') and col.endswith('_ML_with_UDAQ_Adjusted'):
        filled_vocname = col[len('Filled_'):-len('_ML_with_UDAQ_Adjusted')]
        print('filled_voc_name: ', filled_vocname)
        filled_colname = f'Filled_{filled_vocname}_ML_with_UDAQ_Adjusted'

        if filled_vocname not in reverse_mapping:
            continue
        
        ml_name_init = reverse_mapping[filled_vocname]
        init_colname = f'NOAA_{ml_name_init}_Initial'
        print('name_init: ', ml_name_init)

        if init_colname not in df_ml_udaq_initial_and_filled_ml.columns:
            continue

        idx = get_row_index(df_co_and_voc_results, vocname_ml = ml_name_init, vocname_udaq = filled_vocname)
        
        if df_ml_udaq_initial_and_filled_ml[filled_colname].isna().all():
            print(filled_vocname, ' has no measurements.')
            # # explicitly mark all auto + annotation fields as NaN
            # for field in calculated_fields + input_required_fields:
            #     df_co_and_voc_results.loc[idx, field] = np.nan
            # df_co_and_voc_results.to_csv(tracer_csv_savepath, index=False)
            # continue
        
        else:

            # plot_carbon_monoxide_tracer_timeseries(df_index = df_ml_udaq_initial_and_filled_ml.index, 
            #                                        filled_col = df_ml_udaq_initial_and_filled_ml[filled_colname], 
            #                                        vocname_udaq = filled_vocname)
            # (voc_nonans, co_nonans) = scatterplot_co_and_voc(
            #     filled_col = df_ml_udaq_initial_and_filled_ml[filled_colname],
            #     vocname_udaq = filled_vocname)
            # (df_co_and_voc_results, voc_reconstructed) = co_tracer_voc_regression(
            #     idx,
            #     filled_col = df_ml_udaq_initial_and_filled_ml[filled_colname],
            #     vocname_udaq = filled_vocname)
            
            # plot_voc_reconstructed(voc_reconstructed, 
            #                        df_index =  df_ml_udaq_initial_and_filled_ml.index, 
            #                        filled_col = df_ml_udaq_initial_and_filled_ml[filled_colname],
            #                        vocname_udaq = filled_vocname)
            # #save VOC data
            # df_co_and_voc_results = df_co_and_voc_results.to_csv(tracer_csv_savepath, index=False)