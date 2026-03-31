#region: packages
import pandas as pd
import numpy as np
from scipy import odr

import xarray as xr

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
import pickle
#endregion

#region: necessary filepaths
dirpath = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/'
hawthorne_data_dir = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/Hawthorne_data/data/'
hawthorne_script_output_data_dir = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/Hawthorne_data/data/script_output/'
gap_filling_data_dir = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/Merge_scripts/gap_filling/'
stage_data_dirs = {
    'stage1': gap_filling_data_dir + 'stage1_odr_fitting/',
    'stage2': gap_filling_data_dir + 'stage2_filling/',
    'stage3': gap_filling_data_dir + 'stage3_co_tracer/',
    'stage3_5': gap_filling_data_dir + 'stage3_5_co_tracer/',
}
#Make directory for cache storage
CACHE_DIR = '../cache'
os.makedirs(CACHE_DIR, exist_ok=True)
#endregion
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
###############################################################

def merge_all_data_sources():
    #Load UDAQ and NOAA ML VOC files respectively
    udaq_f = hawthorne_script_output_data_dir + 'hawthorne_udaq_all_vocs_15min_timezone_carbon_number_updated.csv'
    udaq = pd.read_csv(udaq_f, index_col='time_local', parse_dates=True)
    #For some reason, pandas is reading the time_local as UTC so we change it back to reading as the UTC-6 time zone
    udaq.index = udaq.index.tz_localize(None)
    udaq.index = udaq.index.tz_localize('America/Denver')
    #Datetime index has 15 min intervals, and from from 07/14/2024 18:00:00 to 8/18/2024 17:45:00 MDT
    #Set index to only include 2024-07-15 00:00:00 to 2024-08-18 17:45:00 MDT
    udaq = udaq.loc['2024-07-15 00:00:00':'2024-08-18 17:45:00']
    # Set any inf, neg. inf, and negative values to NaN
    for col in udaq.columns[1:len(udaq.columns)+1]:
        udaq[col] = udaq[col].replace([np.inf, -np.inf], np.nan).mask(udaq[col] < 0, np.nan)
    udaq = udaq.drop(['Formaldehyde'], axis=1)

    #Sum mEthylToluene and pEthylToluene into 'mpEthyltoluene' and change mapping so that 'mpEthyltoluene' maps to 'x3_x4_EthylToluene_WAS'
    udaq['mpEthyltoluene'] = udaq['mEthyltoluene'] + udaq['pEthyltoluene']
    udaq = udaq.drop(['mEthyltoluene','pEthyltoluene'], axis = 1)
    # voc_mapping['mpEthyltoluene'] = 'x3_x4_EthylToluene_WAS'
    # del voc_mapping['mEthyltoluene']
    # del voc_mapping['pEthyltoluene']
    for col in udaq.columns:
        print('udaq col: ', col)
    
    formaldehyde_f = hawthorne_script_output_data_dir + 'hawthorne_udaq_Formaldehyde_15min_reindexed_timezone_updated.csv'
    df_formaldehyde = pd.read_csv(formaldehyde_f, index_col='time_local', parse_dates=True)
    df_formaldehyde.index = df_formaldehyde.index.tz_localize(None)
    df_formaldehyde.index = df_formaldehyde.index.tz_localize('America/Denver')
    df_formaldehyde['Formaldehyde'] = df_formaldehyde['Formaldehyde'].replace([np.inf, -np.inf], np.nan).mask(df_formaldehyde['Formaldehyde'] < 0, np.nan)
    df_formaldehyde = df_formaldehyde.loc['2024-07-15 00:00:00':'2024-08-18 17:45:00']
    
    ozone_f = hawthorne_script_output_data_dir +'hawthorne_udaq_o3_2024_15min_reindexed_timezone_updated.csv'
    df_o3 = pd.read_csv(ozone_f, index_col='time_local', parse_dates=True)
    df_o3.index = df_o3.index.tz_localize(None)
    df_o3.index = df_o3.index.tz_localize('America/Denver')
    df_o3 = df_o3.loc['2024-07-15 00:00:00':'2024-08-18 17:45:00']
    
    terpenes_f = hawthorne_script_output_data_dir + 'hawthorne_udaq_isoprene_alpha_beta_pinene_07152024_08012024_15min_reindexed_timezone_updated.csv'
    df_terpenes =pd.read_csv(terpenes_f, index_col='time_local', parse_dates=True)
    df_terpenes.index = df_terpenes.index.tz_localize(None)
    df_terpenes.index = df_terpenes.index.tz_localize('America/Denver')
    #Set index to only include 2024-07-15 00:00:00 to 2024-08-18 17:45:00 MDT
    df_terpenes = df_terpenes.loc['2024-07-15 00:00:00':'2024-08-18 17:45:00']
    df_terpenes = df_terpenes.drop(['Isoprene'], axis=1)
    for col in df_terpenes.columns:
        df_terpenes[col] = df_terpenes[col].replace([np.inf, -np.inf], np.nan).mask(df_terpenes[col] < 0, np.nan)

    no_f = hawthorne_script_output_data_dir + 'hawthorne_udaq_no_07152024_08012024_15min_reindexed_timezone_updated.csv'
    df_no =pd.read_csv(no_f, index_col='time_local', parse_dates=True)
    df_no.index = df_no.index.tz_localize(None)
    df_no.index = df_no.index.tz_localize('America/Denver')
    #Set index to only include 2024-07-15 00:00:00 to 2024-08-18 17:45:00 MDT
    df_no = df_no.loc['2024-07-15 00:00:00':'2024-08-18 17:45:00']
    
    df_no['NO'] = df_no['NO'].replace([np.inf, -np.inf], np.nan).mask(df_no['NO'] < 0, np.nan)

    noy_f = hawthorne_script_output_data_dir + 'hawthorne_udaq_noy_07152024_08012024_15min_reindexed_timezone_updated.csv'
    df_noy =pd.read_csv(noy_f, index_col='time_local', parse_dates=True)
    df_noy.index = df_noy.index.tz_localize(None)
    df_noy.index = df_noy.index.tz_localize('America/Denver')
    #Set index to only include 2024-07-15 00:00:00 to 2024-08-18 17:45:00 MDT
    df_noy = df_noy.loc['2024-07-15 00:00:00':'2024-08-18 17:45:00']
    df_noy['NOy'] = df_noy['NOy'].replace([np.inf, -np.inf], np.nan).mask(df_noy['NOy'] < 0, np.nan)

    no2_f = hawthorne_script_output_data_dir + 'hawthorne_udaq_no2_07152024_08012024_15min_reindexed_timezone_updated.csv'
    df_no2 =pd.read_csv(no2_f, index_col='time_local', parse_dates=True)
    df_no2.index = df_no2.index.tz_localize(None)
    df_no2.index = df_no2.index.tz_localize('America/Denver')
    #Set index to only include 2024-07-15 00:00:00 to 2024-08-18 17:45:00 MDT
    df_no2 = df_no2.loc['2024-07-15 00:00:00':'2024-08-18 17:45:00']
    df_no2['NO2'] = df_no2['NO2'].replace([np.inf, -np.inf], np.nan).mask(df_no2['NO2'] < 0, np.nan)

    noaa_f= dirpath + '/CampaignData_and_Merges/R0/CSL_MobileLab_Parked/merged/rev_15min/all_CSL_MobileLab_Parked_rev15min_iWASupdated.nc'

    ds = xr.open_dataset(noaa_f)
    noaa = ds.to_dataframe()
    noaa = noaa.set_index(['time_local'])
    #localize time zone to MDT
    noaa.index = noaa.index.tz_localize('America/Denver')
    #Set index to only include 2024-07-15 00:00:00 to 2024-08-18 17:45:00 MDT
    noaa = noaa.loc['2024-07-15 00:00:00':'2024-08-18 17:45:00']

    # Set any inf, neg. inf, and negative values to NaN
    for col in noaa.columns:
        noaa[col] = noaa[col].replace([np.inf, -np.inf], np.nan).mask(noaa[col] < 0, np.nan)

    #There seems to be a data issue with the CO Picarro during 07/22 15:00 MDT to 07/24 00:00 MDT so mask that range into NaNs
    noaa.loc['2024-07-22 11:00:00':'2024-07-22 18:00:00', 'CO_Piccaro'] = np.nan
    noaa.loc['2024-07-23 10:30:00':'2024-07-23 12:00:00', 'CO_Piccaro'] = np.nan


    df_all_measured = noaa.join([udaq, df_formaldehyde, df_o3, df_terpenes, df_no, df_noy, df_no2])
    #save as new csv file
    savepath = gap_filling_data_dir + 'all_measured_species.csv'
    df_all_measured.to_csv(savepath)
    print('Saved to:' + savepath)

#Given series as inputs, return tuple that masks overlap
#vars_overlap, original = mask_overlap(s1, s2, s3)

def mask_overlap(*vars_list):
    mask = vars_list[0].notna()
    for s in vars_list[1:]:
        mask &= s.notna()

    vars_overlap = [s.loc[mask] for s in vars_list]
    return vars_overlap, list(vars_list)
def run_odr(x_data_odr, y_data_odr):
    (vars_overlap, vars_original) = mask_overlap(x_data_odr, y_data_odr)
    xvar_overlap = vars_overlap[0]
    yvar_overlap = vars_overlap[1]
    xvar_init_data = vars_original[0]
    yvar_init_data = vars_original[1]
    print('len(vars_overlap[0]): ', len(vars_overlap[0]))
    #print('vars_overlap[0]: ', vars_overlap[0].to_numpy())
    print('len(vars_overlap[1]): ', len(vars_overlap[1]))
    #print('vars_overlap[1]: ', vars_overlap[1].to_numpy())
    print('len(vars_original[0]): ', len(vars_original[0]))
    print('len(vars_original[1]): ', len(vars_original[1]))

    #Get count of how many points are overlapping, and used in the ODR fitting
    points_considered_in_odr = len(vars_overlap[0])

    #From the ODRpack documentation
    # Define the function we want to fit against
    x_data_overlap_np_arr = vars_overlap[0].to_numpy()
    #xdata_2d = np.array([udaq_overlap_voc_np_arr])
    y_data_overlap_np_arr = vars_overlap[1].to_numpy()
    #ydata_2d = np.array([ml_overlap_voc_np_arr])

    # print('xdata_2d shape: ', np.shape(udaq_overlap_voc_np_arr))
    # print('ydata shape: ', np.shape(ml_overlap_voc_np_arr))

    def linear_model(x, beta):
        return beta[0] + beta[1] * x
    
    #beta0 is an initial guess
    odr_fit_result = odr_fit(linear_model, x_data_overlap_np_arr, y_data_overlap_np_arr, [0.0, 1.0])

    odr_intercept, odr_slope = odr_fit_result.beta

    correction_eq = f"Correction applied: O_hat = {odr_intercept:.4f} + {odr_slope:.4f} * M"

    xvals_corrected = odr_intercept + odr_slope * x_data_odr

    #overlap between initial and fitted data, no NaNs
    xcorrected_overlap = xvals_corrected.loc[xvals_corrected.index.intersection(vars_overlap[0].index)]
    return xvar_overlap, yvar_overlap, xvar_init_data, yvar_init_data, correction_eq, xvals_corrected, xcorrected_overlap, points_considered_in_odr
def fitting_metrics(x_overlap, y_overlap):
    slope, intercept = np.polyfit(x_overlap, y_overlap, 1)
    rmse = np.sqrt(np.mean((x_overlap - y_overlap)**2))
    
    #R squared calculation
    r = np.corrcoef(x_overlap, y_overlap)[0, 1]
    r2 = r**2

    return slope, intercept, rmse, r2
def odr_improves(metrics_init_input, metrics_corr_input, y_overlap):
    slope_i, intercept_i, rmse_i, r2_i = metrics_init_input
    slope_c, intercept_c, rmse_c, r2_c = metrics_corr_input

    did_slope_improve = bool(abs(1 - slope_c) < abs(1 - slope_i))
    slope_distance_i = abs(1 - slope_i)
    slope_distance_c = abs(1 - slope_c)
    slope_distance_improvement_val = slope_distance_i - slope_distance_c

    intercept_err_i = abs(intercept_i) / np.nanmean(y_overlap)
    intercept_err_c = abs(intercept_c) / np.nanmean(y_overlap)

    did_rmse_improve = bool(rmse_c < rmse_i)
    rmse_percent_improvement = 100*((rmse_i - rmse_c) / (rmse_i))
    rmse_norm_i = rmse_i/np.nanmean(y_overlap)
    rmse_norm_c = rmse_c/np.nanmean(y_overlap)

    r2_round_init = round(r2_i, 2)
    r2_at_least_half = bool(r2_round_init >= 0.5)

    score_i = slope_distance_i + intercept_err_i + rmse_norm_i
    score_c = slope_distance_c + intercept_err_c + rmse_norm_c
    #Apply correction if initial score is higher than corrected score
    apply_correction_based_off_score = bool(score_c < score_i)

    #points_considered_in_odr_init should be same as points_considered_in_odr_corr

    return did_slope_improve, slope_distance_i, slope_distance_c, slope_distance_improvement_val, intercept_err_i, intercept_err_c, did_rmse_improve, rmse_percent_improvement, rmse_norm_i, rmse_norm_c, r2_at_least_half, score_i, score_c, apply_correction_based_off_score
def process_one_voc_stage1(var_species_name, duplicate_species_name, var_species_data, duplicate_species_data):
    cache_file_stage1 = f'{CACHE_DIR}/{var_species_name}_stage1.pkl'
    
    # ---- load if already processed ----
    # if os.path.exists(cache_file_stage1):
    #     with open(cache_file_stage1, 'rb') as f:
    #         return pickle.load(f)

        # If any of the results need to be modified, comment out the return pickle.load(f) line and replace with code below with alterations. 
        # This one changes a key previously named 'metrics_ml_udaq_corrected' into 'metrics_ml_udaq_corr'
        #     results = pickle.load(f)
        # print('results: ', results)
        # if 'metrics_ml_udaq_corrected' in results:
        #     results['metrics_ml_udaq_corr'] = results.pop('metrics_ml_udaq_corrected')
        # with open(cache_file, 'wb') as f:
        #     pickle.dump(results, f)

    # -----------------------------
    # CASE: All NaNs in var
    # -----------------------------
    if var_species_data.isna().all():

        print('Processing ', var_species_name, ': All NaNs in var')

        stage1_results = {
            'x_species_name': duplicate_species_name,
            'y_species_name': var_species_name,
            'case': 'All NaNs in var',
            'points_considered_in_odr': np.nan,
            'odr_eq_adj': np.nan,
            'metrics_init': (np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan),
            'metrics_corr': (np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan),
            'did_slope_improve': np.nan, 
            'slope_distance_from_1_init': np.nan,
            'slope_distance_from_1_corr': np.nan,
            'slope_distance_improvement_val': np.nan,
            'intercept_error_init': np.nan,
            'intercept_error_corr': np.nan,
            'did_rmse_improve': np.nan, 
            'rmse_percent_improvement': np.nan,
            'rmse_normalized_init': np.nan,
            'rmse_normalized_corr': np.nan,
            'r2_at_least_half': np.nan,
            'total_slope_intercept_rmse_score_init': np.nan, 
            'total_slope_intercept_rmse_score_corr': np.nan,
            'should_correction_be_applied_from_score_eval': np.nan,
            'var_data_init_full': var_species_data,
            'var_data_init_overlap': pd.Series(np.nan, index=df_all_measured_species.index),
            'dup_data_init_full': duplicate_species_data,
            'dup_data_init_overlap': pd.Series(np.nan, index=df_all_measured_species.index),
            'dup_data_corr_full': pd.Series(np.nan, index=df_all_measured_species.index),
            'dup_data_corr_overlap': pd.Series(np.nan, index=df_all_measured_species.index)
        }
    elif isinstance(duplicate_species_name, list):
        print('Processing ', var_species_name, ': Has duplicate list')

        run_odr_outputs = [run_odr(x_data_odr = data, y_data_odr= var_species_data) for data in duplicate_species_data]
        for i, out in enumerate(run_odr_outputs):
            print(i, len(out))
        #Metrics/Evaluation of ODR fitting       
        #next_outputs = [next_function(out, other_input) for xvar_overlap, yvar_overlap in (xvar_overlap, yvar_overlap, xvar_init_data, yvar_init_data, correction_eq, xvals_corrected, xcorrected_overlap, points_considered_in_odr)]
        metrics_init = [fitting_metrics(xvar_overlap, yvar_overlap) for xvar_overlap, yvar_overlap, xvar_init_data, yvar_init_data, correction_eq, xvals_corrected, xcorrected_overlap, points_considered_in_odr in run_odr_outputs]
        metrics_corr  = [fitting_metrics(xcorrected_overlap, yvar_overlap) for xvar_overlap, yvar_overlap, xvar_init_data, yvar_init_data, correction_eq, xvals_corrected, xcorrected_overlap, points_considered_in_odr in run_odr_outputs]
        
        # #Quality control decision based off of if improvements are made by applying correction to UDAQ data
        improvement_stats = [odr_improves(metrics_init_input = metrics_init_tuple, metrics_corr_input = metrics_corr_tuple, y_overlap = var_overlap) for metrics_init_tuple, metrics_corr_tuple, (_, var_overlap, *_) in zip(metrics_init, metrics_corr, run_odr_outputs)]
       
        stage1_results = {
            'x_species_name': duplicate_species_name,
            'y_species_name': var_species_name,
            'case': 'Has duplicate list, has duplicate data not all NaNs, has vardata not all NaNs',
            'points_considered_in_odr': [out[7] for out in run_odr_outputs],
            'odr_eq_adj': [out[4] for out in run_odr_outputs], #equation used to correct the UDAQ data
            'metrics_init': metrics_init,
            'metrics_corr': metrics_corr,
            'did_slope_improve': [out[0] for out in improvement_stats], 
            'slope_distance_from_1_init': [out[1] for out in improvement_stats], 
            'slope_distance_from_1_corr': [out[2] for out in improvement_stats], 
            'slope_distance_improvement_val': [out[3] for out in improvement_stats], 
            'intercept_error_init': [out[4] for out in improvement_stats], 
            'intercept_error_corr': [out[5] for out in improvement_stats], 
            'did_rmse_improve': [out[6] for out in improvement_stats], 
            'rmse_percent_improvement': [out[7] for out in improvement_stats], 
            'rmse_normalized_init': [out[8] for out in improvement_stats], 
            'rmse_normalized_corr': [out[9] for out in improvement_stats], 
            'r2_at_least_half': [out[10] for out in improvement_stats], 
            'total_slope_intercept_rmse_score_init': [out[11] for out in improvement_stats], 
            'total_slope_intercept_rmse_score_corr': [out[12] for out in improvement_stats], 
            'should_correction_be_applied_from_score_eval': [out[13] for out in improvement_stats], 
            'var_data_init_full': var_species_data,
            'var_data_init_overlap': [out[1] for out in run_odr_outputs],
            'dup_data_init_full': [out[2] for out in run_odr_outputs],
            'dup_data_init_overlap': [out[0] for out in run_odr_outputs],
            'dup_data_corr_full': [out[5] for out in run_odr_outputs],
            'dup_data_corr_overlap': [out[6] for out in run_odr_outputs]
        }
        
    # -----------------------------
    # CASE: All NaNs in Duplicate
    # -----------------------------
    elif duplicate_species_data.isna().all():

        print('Processing ', var_species_name, ': All NaNs in duplicate')

        stage1_results = {
            'x_species_name': duplicate_species_name,
            'y_species_name': var_species_name,
            'case': 'All NaNs in duplicate',
            'points_considered_in_odr': np.nan,
            'odr_eq_adj': np.nan,
            'metrics_init': (np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan),
            'metrics_corr': (np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan),
            'did_slope_improve': np.nan, 
            'slope_distance_from_1_init': np.nan,
            'slope_distance_from_1_corr': np.nan,
            'slope_distance_improvement_val': np.nan,
            'intercept_error_init': np.nan,
            'intercept_error_corr': np.nan,
            'did_rmse_improve': np.nan, 
            'rmse_percent_improvement': np.nan,
            'rmse_normalized_init': np.nan,
            'rmse_normalized_corr': np.nan,
            'r2_at_least_half': np.nan,
            'total_slope_intercept_rmse_score_init': np.nan, 
            'total_slope_intercept_rmse_score_corr': np.nan,
            'should_correction_be_applied_from_score_eval': np.nan,
            'var_data_init_full': var_species_data,
            'var_data_init_overlap': pd.Series(np.nan, index=df_all_measured_species.index),
            'dup_data_init_full': duplicate_species_data,
            'dup_data_init_overlap': pd.Series(np.nan, index=df_all_measured_species.index),
            'dup_data_corr_full': pd.Series(np.nan, index=df_all_measured_species.index),
            'dup_data_corr_overlap': pd.Series(np.nan, index=df_all_measured_species.index)
        }

    else:

        print('Processing ', var_species_name, ': Has duplicate, has duplicate data not all NaNs, has vardata not all NaNs')
        
        #ODR fitting on UDAQ data
        (xvar_overlap, yvar_overlap, xvar_init_data, yvar_init_data, correction_eq, xvals_corrected, xcorrected_overlap, points_considered_in_odr) = run_odr(x_data_odr = duplicate_species_data, 
                                                                                                y_data_odr = var_species_data)

        #Metrics/Evaluation of ODR fitting
        (vars_overlap, vars_original) = mask_overlap(duplicate_species_data, var_species_data)
        
        metrics_init = fitting_metrics(vars_overlap[0], vars_overlap[1])

        metrics_corr  = fitting_metrics(xcorrected_overlap, vars_overlap[1])

        #Quality control decision based off of if improvements are made by applying correction to UDAQ data
        (did_slope_improve, slope_distance_i, 
        slope_distance_c, slope_distance_improvement_val, 
        intercept_err_i, intercept_err_c,
        did_rmse_improve, rmse_percent_improvement,
        rmse_norm_i, rmse_norm_c,
        r2_at_least_half, score_i, score_c,
        apply_correction_based_off_score) = odr_improves(metrics_init, metrics_corr, vars_overlap[1])

        stage1_results = {
            'x_species_name': duplicate_species_name,
            'y_species_name': var_species_name,
            'case': 'Has duplicate, has duplicate data not all NaNs, has vardata not all NaNs',
            'points_considered_in_odr': points_considered_in_odr,
            'odr_eq_adj': correction_eq, #equation used to correct the UDAQ data
            'metrics_init': metrics_init,
            'metrics_corr': metrics_corr,
            'did_slope_improve': did_slope_improve, 
            'slope_distance_from_1_init': slope_distance_i,
            'slope_distance_from_1_corr': slope_distance_c,
            'slope_distance_improvement_val': slope_distance_improvement_val,
            'intercept_error_init': intercept_err_i,
            'intercept_error_corr': intercept_err_c,
            'did_rmse_improve': did_rmse_improve, 
            'rmse_percent_improvement': rmse_percent_improvement,
            'rmse_normalized_init': rmse_norm_i, 
            'rmse_normalized_corr': rmse_norm_c,
            'r2_at_least_half': r2_at_least_half, 
            'total_slope_intercept_rmse_score_init': score_i, 
            'total_slope_intercept_rmse_score_corr': score_c,
            'should_correction_be_applied_from_score_eval': apply_correction_based_off_score,
            'var_data_init_full': var_species_data,
            'var_data_init_overlap': vars_overlap[1],
            'dup_data_init_full': duplicate_species_data,
            'dup_data_init_overlap': vars_overlap[0],
            'dup_data_corr_full': xvals_corrected,
            'dup_data_corr_overlap': xcorrected_overlap
        }
        print('Resolved ' + var_species_name + ' Stage1, one VOC')
    
    # save so ODR never recomputes
    with open(cache_file_stage1, 'wb') as f:
        pickle.dump(stage1_results, f)

    return stage1_results
def process_all_vocs_stage1():
    stage1_results_all = []
    
    for row_idx in range(0, len(df_duplicates_and_tracers.index)):
        spec_name = df_duplicates_and_tracers['Varname'].iloc[row_idx]
        print('Processing Stage 1 for: ', spec_name)
        dup_name = df_duplicates_and_tracers['Duplicate_name'].iloc[row_idx]

        print('type(dup_name): ', type(dup_name))
        if pd.isna(dup_name):
            dup_species_data = pd.Series(np.nan, index=df_all_measured_species.index)
        elif pd.notna(dup_name) and ';' in dup_name:
            dup_name_list = dup_name.split('; ')
            print('Has more than one dup:', dup_name_list)
            dup_name = dup_name_list
            dup_species_data = [df_all_measured_species[dup_name[0]], df_all_measured_species[dup_name[0]]]
            print('dup_species_data type should be a list: ', type(dup_species_data))

        else:
            print('Processing Stage 1 with ', dup_name)
            dup_species_data = df_all_measured_species[dup_name]

        stage1_results = process_one_voc_stage1(var_species_name = spec_name, 
                                                duplicate_species_name = dup_name, 
                                                var_species_data = df_all_measured_species[spec_name], 
                                                duplicate_species_data = dup_species_data)
        
        stage1_results_all.append(stage1_results)
    return stage1_results_all
def build_summary_stage1(stage1_results_all):
    rows = []
    for r in stage1_results_all:
        print('Saving for :', r['y_species_name'])
        print(r['x_species_name'])
        print(type(r['x_species_name']))

        if isinstance(r['x_species_name'], list):
            print('Benzene Catch')
            rows.append({    
            'varname': r['y_species_name'],
            'Duplicate_name': r['x_species_name'],
            'points_considered_in_odr': r['points_considered_in_odr'],
            'case': r['case'],
            'odr_eq_adj': r['odr_eq_adj'],
            'metrics_init_slope': [out[0] for out in r['metrics_init']], 
            'metrics_init_intercept': [out[1] for out in r['metrics_init']],
            'metrics_init_rmse': [out[2] for out in r['metrics_init']],
            'metrics_init_r2': [out[3] for out in r['metrics_init']],
            'metrics_corr_slope': [out[0] for out in r['metrics_corr']], 
            'metrics_corr_intercept': [out[1] for out in r['metrics_corr']], 
            'metrics_corr_rmse': [out[2] for out in r['metrics_corr']], 
            'metrics_corr_r2': [out[3] for out in r['metrics_corr']], 
            'did_slope_improve': r['did_slope_improve'],
            'slope_distance_from_1_init': r['slope_distance_from_1_init'],
            'slope_distance_from_1_corr': r[ 'slope_distance_from_1_corr'],
            'slope_distance_improvement_val': r['slope_distance_improvement_val'],
            'intercept_error_init': r['intercept_error_init'],
            'intercept_error_corr': r['intercept_error_corr'],
            'did_rmse_improve': r['did_rmse_improve'], 
            'rmse_percent_improvement': r['rmse_percent_improvement'],
            'rmse_normalized_init': r['rmse_normalized_init'],
            'rmse_normalized_corr': r['rmse_normalized_corr'],
            'r2_at_least_half': r['r2_at_least_half'],
            'total_slope_intercept_rmse_score_init': r['total_slope_intercept_rmse_score_init'],
            'total_slope_intercept_rmse_score_corr': r['total_slope_intercept_rmse_score_corr'],
            'should_correction_be_applied_from_score_eval': r['should_correction_be_applied_from_score_eval']
            })

        else:
            rows.append({    
            'varname': r['y_species_name'],
            'Duplicate_name': r['x_species_name'],
            'points_considered_in_odr': r['points_considered_in_odr'],
            'case': r['case'],
            'odr_eq_adj': r['odr_eq_adj'],
            'metrics_init_slope': r['metrics_init'][0],
            'metrics_init_intercept': r['metrics_init'][1],
            'metrics_init_rmse': r['metrics_init'][2],
            'metrics_init_r2': r['metrics_init'][3],
            'metrics_corr_slope': r['metrics_corr'][0],
            'metrics_corr_intercept': r['metrics_corr'][1],
            'metrics_corr_rmse': r['metrics_corr'][2],
            'metrics_corr_r2': r['metrics_corr'][3],
            'did_slope_improve': r['did_slope_improve'],
            'slope_distance_from_1_init': r['slope_distance_from_1_init'],
            'slope_distance_from_1_corr': r[ 'slope_distance_from_1_corr'],
            'slope_distance_improvement_val': r['slope_distance_improvement_val'],
            'intercept_error_init': r['intercept_error_init'],
            'intercept_error_corr': r['intercept_error_corr'],
            'did_rmse_improve': r['did_rmse_improve'], 
            'rmse_percent_improvement': r['rmse_percent_improvement'],
            'rmse_normalized_init': r['rmse_normalized_init'],
            'rmse_normalized_corr': r['rmse_normalized_corr'],
            'r2_at_least_half': r['r2_at_least_half'],
            'total_slope_intercept_rmse_score_init': r['total_slope_intercept_rmse_score_init'],
            'total_slope_intercept_rmse_score_corr': r['total_slope_intercept_rmse_score_corr'],
            'should_correction_be_applied_from_score_eval': r['should_correction_be_applied_from_score_eval']
            })
    stage1_summary_savepath = stage_data_dirs['stage1'] + 'csv_storage/stage1_summary_new.csv'
    pd.DataFrame(rows).to_csv(stage1_summary_savepath, index=True)
    print('Saved summary to: ' + str(stage1_summary_savepath))

def odr_linear_func(B, x):
    """Linear function for ODR: B0*x0 + B1*x1 + ... + Bn*xn + intercept"""
    #return np.dot(x, B[:-1]) + B[-1]
    print('length B: ', len(B[:-1]))
    print('x shape: ', x.shape)
    #return B[0]*x[:,0] + B[1]   # x must be 2D
    #return np.dot(x, B[:-1]) + B[-1]
    return np.dot(B[:-1], x) + B[-1]

def fit_odr(X_fit, y_fit):
    mask = ~np.isnan(y_fit) & ~np.isnan(X_fit).any(axis=1)
    X_fit_val = X_fit[mask].astype(float)
    y_fit_val = y_fit[mask].astype(float)

    print(np.nanstd(X_fit_val))

    # Convert to numpy
    X_np = X_fit_val.to_numpy()
    y_np = y_fit_val.to_numpy()
    # IMPORTANT: transpose predictors for ODR
    X_odr = X_np.T   # shape (n_tracers, n_obs)

    beta0 = np.zeros(X_odr.shape[0] + 1)
    """Fit ODR linear model"""
    # print('Running ODR model')
    model = odr.Model(odr_linear_func)
    print(model)

    data = odr.Data(X_odr, y_np)
    # print('X.values shape: ', X_fit_val.shape, 'X.values type: ', type(X_fit_val))
    # print('y.values shape: ', y_fit_val.shape, 'y.values type: ', type(y_fit_val))
    #print('X_fit_val: ', X_fit_val['BrCl_CIMS'].tolist())
    #print(X_fit_val.shape)
    #print('y_fit_val: ', y_fit_val.tolist())
 
    odr_run = odr.ODR(data, model, beta0=beta0)
    out = odr_run.run()
    return out

def predict_odr(out, X):
    """Predict y using fitted ODR model"""
    print('Predict y using fitted ODR model')
    return odr_linear_func(out.beta, X.values)

def fill_species_gaps_conditional_odr(df, species_tracer_dict, regime_col='day_night',
                                      coverage_threshold=0.2, corr_threshold=0.9):
    """
    Conditional-regime gap filling for multiple species using ODR regression.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with species and tracers.
    species_tracer_dict : dict
        Keys: target species
        Values: list of tracer columns
    regime_col : str or None
        Column in df that indicates regime (e.g., day/night). If None, all data treated as one regime.
    coverage_threshold : float
        Minimum fraction of full overlap for multiple regression to be considered.
    corr_threshold : float
        Maximum allowed correlation between tracers for multiple regression.
        
    Returns
    -------
    filled_df : pd.DataFrame
        Filled species columns.
    fill_source : pd.DataFrame
        Source of each filled value.
    """
    filled_df = df.copy()
    fill_source = pd.DataFrame(index=df.index)
    
    # Determine regimes
    if regime_col is not None:
        regimes = df[regime_col].unique()
    else:
        regimes = [None]
    
    for regime in regimes:
        if regime is not None:
            mask_regime = df[regime_col] == regime
        else:
            mask_regime = pd.Series(True, index=df.index)
        
        for target, tracers in species_tracer_dict.items():
            print('Target species: ', target)
            print('Tracer species: ', tracers)

            y = filled_df.loc[mask_regime, target]
            
            filled_col = y.copy()
            source_col = pd.Series(['measured']*len(y), index=y.index)
            
            print('Now processing Step 1: Multiple regression feasibility')
            # --- Step 1: Multiple regression feasibility ---
            predictor_overlap = filled_df.loc[mask_regime, tracers].notna().all(axis=1).mean()
            
            if predictor_overlap < coverage_threshold:
                use_multi = False
            else:
                corr_matrix = filled_df.loc[mask_regime, tracers].corr().abs()
                np.fill_diagonal(corr_matrix.values, 0)
                max_corr = corr_matrix.max().max()
                use_multi = max_corr < corr_threshold
            
            # --- Step 2: Fit multiple regression using ODR ---
            print('Now processing Step 2: Fit multiple regression using ODR')
            if use_multi:
                mask_fit = y.notna() & filled_df.loc[mask_regime, tracers].notna().all(axis=1)
                X_fit = filled_df.loc[mask_regime, tracers].loc[mask_fit]
                y_fit = y.loc[mask_fit]
                if len(y_fit) > 1:
                    odr_out = fit_odr(X_fit, y_fit)
                    
                    # Apply ODR to missing points where all tracers exist
                    mask_pred = y.isna() & filled_df.loc[mask_regime, tracers].notna().all(axis=1)
                    if mask_pred.any():
                        filled_col.loc[mask_pred] = predict_odr(odr_out, filled_df.loc[mask_regime, tracers].loc[mask_pred])
                        source_col.loc[mask_pred] = 'multi_odr'

            print('Now processing Step 3: Sequential single-tracer ODR filling')
            # --- Step 3: Sequential single-tracer ODR filling ---
            remaining_mask = filled_col.isna()
            for tracer in tracers:
                print('Tracer: ', tracer)
                mask_tracer = remaining_mask & filled_df.loc[mask_regime, tracer].notna()
                if mask_tracer.any():
                    mask_fit = y.notna() & filled_df.loc[mask_regime, tracer].notna()
                    X_fit = filled_df.loc[mask_regime, [tracer]].loc[mask_fit]

                    print('X_fit Shape: ', X_fit.shape)
                    y_fit = y.loc[mask_fit]
                    print('y_fit Shape: ', y_fit.shape)
                    print('length y_fit: ', len(y_fit))
                    if len(y_fit) > 1:
                        print('Entering odr_out_single')
                        odr_out_single = fit_odr(X_fit, y_fit)
                        print('Entering predict_odr')
                        filled_col.loc[mask_tracer] = predict_odr(odr_out_single, filled_df.loc[mask_regime, [tracer]].loc[mask_tracer])
                        source_col.loc[mask_tracer] = f'{tracer}_odr_fill'
                        remaining_mask = filled_col.isna()
                        if not remaining_mask.any():
                            break

            print('Now processing Step 4: Save results')
            # --- Step 4: Save results ---
            filled_df.loc[mask_regime, target + '_filled'] = filled_col
            fill_source.loc[mask_regime, target + '_source'] = source_col
    
    return filled_df, fill_source

#TO DO
#Save formaldehyde data
#Merge UDAQ and NOAA ML data sources
if __name__ == "__main__":
    all_measured_species_f = gap_filling_data_dir + 'all_measured_species.csv'
    df_all_measured_species = pd.read_csv(all_measured_species_f, index_col='time_local', parse_dates=True)
    df_all_measured_species.index = df_all_measured_species.index.tz_localize(None)
    df_all_measured_species.index = df_all_measured_species.index.tz_localize('America/Denver')

    duplicates_and_tracers_f = gap_filling_data_dir + 'species_duplicates_tracers.csv'
    df_duplicates_and_tracers = pd.read_csv(duplicates_and_tracers_f, index_col='Index')

    # print('Starting Stage 1')

    # stage1_results_all = process_all_vocs_stage1()
    # build_summary_stage1(stage1_results_all)


    df_duplicates_and_tracers['Expected_tracers'] = df_duplicates_and_tracers['Expected_tracers'].apply(lambda x: x.split('; '))
    species_tracer_dict = {}
    #print(len(df_duplicates_and_tracers.index))
    
    for idx in range(0, len(df_duplicates_and_tracers.index)):
        species_tracer_dict.update({df_duplicates_and_tracers['Varname'][idx] : df_duplicates_and_tracers['Expected_tracers'][idx]})
    

    df_all_measured_species['hour'] = df_all_measured_species.index.hour
    df_all_measured_species['day_night'] = (df_all_measured_species['hour'] >= 6) & (df_all_measured_species['hour'] < 18)

    #print(df_all_measured_species)
    print(species_tracer_dict)

    filled_df, fill_source = fill_species_gaps_conditional_odr(df_all_measured_species, species_tracer_dict,
                                                        regime_col='day_night', coverage_threshold=0.2, corr_threshold=0.9)

    print(filled_df)
    print(fill_source)
