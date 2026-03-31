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
from sklearn.metrics import root_mean_squared_error, r2_score

import re
from collections import defaultdict

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator
#import cmasher as cmr
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
from matplotlib.lines import Line2D
from datetime import datetime

import os
import pickle

from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from statsmodels.stats.outliers_influence import variance_inflation_factor
from itertools import combinations

from scipy.io import savemat
from collections import OrderedDict
#endregion

#region: necessary filepaths
dirpath = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/'
hawthorne_data_dir = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/Hawthorne_data/data/'
hawthorne_script_output_data_dir = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/Hawthorne_data/data/script_output/'
gap_filling_data_dir = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/Merge_scripts/gap_filling/'
gap_filling_csv_data_dir = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/Merge_scripts/gap_filling/csv/'
stage_data_dirs = {
    'stage1': gap_filling_data_dir + 'stage1_odr_fitting/',
    'stage2': gap_filling_data_dir + 'stage2_tracers/',
    'stage3': gap_filling_data_dir + 'stage3_filling/',
    'stage4': gap_filling_data_dir + 'stage4_interpolation/'
}
#Make directory for cache storage
CACHE_DIR = '../cache'
os.makedirs(CACHE_DIR, exist_ok=True)
#endregion
#region: Plot formatting
mpl.rcParams['xtick.labelsize'] = 15
mpl.rcParams['ytick.labelsize'] = 15
mpl.rcParams['legend.fontsize'] = 13
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
    # udaq = udaq.loc['2024-07-15 00:00:00':'2024-08-18 17:45:00']
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
    df_formaldehyde = df_formaldehyde.loc['2024-07-14 18:00:00':'2024-08-18 17:45:00']
    
    ozone_f = hawthorne_script_output_data_dir +'hawthorne_udaq_o3_2024_15min_reindexed_timezone_updated.csv'
    df_o3 = pd.read_csv(ozone_f, index_col='time_local', parse_dates=True)
    df_o3.index = df_o3.index.tz_localize(None)
    df_o3.index = df_o3.index.tz_localize('America/Denver')
    df_o3 = df_o3.loc['2024-07-14 18:00:00':'2024-08-18 17:45:00']
    
    terpenes_f = hawthorne_script_output_data_dir + 'hawthorne_udaq_isoprene_alpha_beta_pinene_07152024_08012024_15min_reindexed_timezone_updated.csv'
    df_terpenes =pd.read_csv(terpenes_f, index_col='time_local', parse_dates=True)
    df_terpenes.index = df_terpenes.index.tz_localize(None)
    df_terpenes.index = df_terpenes.index.tz_localize('America/Denver')
    #Set index to only include 2024-07-15 00:00:00 to 2024-08-18 17:45:00 MDT
    df_terpenes = df_terpenes.loc['2024-07-14 18:00:00':'2024-08-18 17:45:00']
    df_terpenes = df_terpenes.drop(['Isoprene'], axis=1)
    for col in df_terpenes.columns:
        df_terpenes[col] = df_terpenes[col].replace([np.inf, -np.inf], np.nan).mask(df_terpenes[col] < 0, np.nan)

    no_f = hawthorne_script_output_data_dir + 'hawthorne_udaq_no_07152024_08012024_15min_reindexed_timezone_updated.csv'
    df_no =pd.read_csv(no_f, index_col='time_local', parse_dates=True)
    df_no.index = df_no.index.tz_localize(None)
    df_no.index = df_no.index.tz_localize('America/Denver')
    #Set index to only include 2024-07-15 00:00:00 to 2024-08-18 17:45:00 MDT
    df_no = df_no.loc['2024-07-14 18:00:00':'2024-08-18 17:45:00']
    
    df_no['NO'] = df_no['NO'].replace([np.inf, -np.inf], np.nan).mask(df_no['NO'] < 0, np.nan)

    noy_f = hawthorne_script_output_data_dir + 'hawthorne_udaq_noy_07152024_08012024_15min_reindexed_timezone_updated.csv'
    df_noy =pd.read_csv(noy_f, index_col='time_local', parse_dates=True)
    df_noy.index = df_noy.index.tz_localize(None)
    df_noy.index = df_noy.index.tz_localize('America/Denver')
    #Set index to only include 2024-07-15 00:00:00 to 2024-08-18 17:45:00 MDT
    df_noy = df_noy.loc['2024-07-14 18:00:00':'2024-08-18 17:45:00']
    df_noy['NOy'] = df_noy['NOy'].replace([np.inf, -np.inf], np.nan).mask(df_noy['NOy'] < 0, np.nan)

    no2_f = hawthorne_script_output_data_dir + 'hawthorne_udaq_no2_07152024_08012024_15min_reindexed_timezone_updated.csv'
    df_no2 =pd.read_csv(no2_f, index_col='time_local', parse_dates=True)
    df_no2.index = df_no2.index.tz_localize(None)
    df_no2.index = df_no2.index.tz_localize('America/Denver')
    #Set index to only include 2024-07-15 00:00:00 to 2024-08-18 17:45:00 MDT
    df_no2 = df_no2.loc['2024-07-14 18:00:00':'2024-08-18 17:45:00']
    df_no2['NO2'] = df_no2['NO2'].replace([np.inf, -np.inf], np.nan).mask(df_no2['NO2'] < 0, np.nan)

    noaa_f= dirpath + '/CampaignData_and_Merges/R0/CSL_MobileLab_Parked/merged/rev_15min/all_CSL_MobileLab_Parked_rev15min_iWASupdated.nc'

    ds = xr.open_dataset(noaa_f)
    noaa = ds.to_dataframe()
    noaa = noaa.set_index(['time_local'])
    #localize time zone to MDT
    noaa.index = noaa.index.tz_localize('America/Denver')
    #Set index to only include 2024-07-15 00:00:00 to 2024-08-18 17:45:00 MDT
    noaa = noaa.loc['2024-07-14 18:00:00':'2024-08-18 17:45:00']

    # Set any inf, neg. inf, and negative values to NaN
    for col in noaa.columns:
        noaa[col] = noaa[col].replace([np.inf, -np.inf], np.nan).mask(noaa[col] < 0, np.nan)

    #There seems to be a data issue with the CO Picarro during 07/22 15:00 MDT to 07/24 00:00 MDT so mask that range into NaNs
    noaa.loc['2024-07-22 11:00:00':'2024-07-22 18:00:00', 'CO_Piccaro'] = np.nan
    noaa.loc['2024-07-23 10:30:00':'2024-07-23 12:00:00', 'CO_Piccaro'] = np.nan

    df_all_measured = noaa.join([udaq, df_formaldehyde, df_o3, df_terpenes, df_no, df_noy, df_no2])
    #save as new csv file
    savepath = gap_filling_csv_data_dir + 'all_measured_species.csv'
    df_all_measured.to_csv(savepath)
    print('Saved to:' + savepath)

#Given series as inputs, return tuple that masks overlap
#vars_overlap, original = mask_overlap(s1, s2, s3)
def mask_overlap_revised1(*vars_list, extra_mask = None):
    mask = vars_list[0].notna()
    for s in vars_list[1:]:
        mask &= s.notna()
    #add extra mask for night/day regime
    if extra_mask is not None:
        mask &= extra_mask
    vars_overlap = [s.loc[mask] for s in vars_list]

    return vars_overlap, list(vars_list)
def mask_overlap_revised(*vars_list, extra_mask = None):
    mask = vars_list[0].notna()
    for s in vars_list[1:]:
        mask &= s.notna()
    #add extra mask for night/day regime
    if extra_mask is not None:
        mask &= extra_mask
        init_data_with_regime = [v.loc[extra_mask] for v in vars_list]
    vars_overlap = [s.loc[mask] for s in vars_list]

    return vars_overlap, list(vars_list), init_data_with_regime
def gap_fill_metrics(tracer_data, target_data, extra_mask):
    gap_fill_potential = ((target_data.isna() & tracer_data.notna()).loc[extra_mask]).sum()
    gaps_needed_to_fill_for_regime = ((target_data.loc[extra_mask]).isna()).sum()
    return gap_fill_potential, gaps_needed_to_fill_for_regime
def run_odr(x_data_odr, y_data_odr):
    (vars_overlap, vars_original) = mask_overlap_revised1(x_data_odr, y_data_odr)

    xvar_overlap = vars_overlap[0]
    yvar_overlap = vars_overlap[1]
    xvar_init_data = vars_original[0]
    yvar_init_data = vars_original[1]
    #print('len(vars_overlap[0]): ', len(vars_overlap[0]))
    #print('vars_overlap[0]: ', vars_overlap[0].to_numpy())
    #print('len(vars_overlap[1]): ', len(vars_overlap[1]))
    #print('vars_overlap[1]: ', vars_overlap[1].to_numpy())
    #print('len(vars_original[0]): ', len(vars_original[0]))
    #print('len(vars_original[1]): ', len(vars_original[1]))

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

    xvals_corrected_overlap = odr_intercept + odr_slope * xvar_overlap
    xvals_corrected = odr_intercept + odr_slope * xvar_init_data

    #Metrics
    rmse= root_mean_squared_error(yvar_overlap, xvar_overlap)
    norm_rmse = rmse / np.std(yvar_overlap)
    r2 = r2_score(yvar_overlap, xvar_overlap)

    # return {'slope': odr_slope,
    #         'intercept': odr_intercept,
    #         'rmse': rmse,
    #         'norm_rmse': norm_rmse,
    #         'r2': r2,
    #         'overlap_points_counted': overlap_points_count,
    #         'score':score_estimate}
    return {'xvar_overlap': xvar_overlap, 
            'yvar_overlap': yvar_overlap, 
            'xvar_init_data': xvar_init_data, 
            'yvar_init_data': yvar_init_data, 
            'odr_slope': odr_slope, 
            'odr_intercept': odr_intercept, 
            'correction_eq': correction_eq, 
            'rmse': rmse, 
            'norm_rmse': norm_rmse, 
            'r2': r2, 
            'xvals_corrected': xvals_corrected,
            'xvals_corrected_overlap': xvals_corrected_overlap,
            'points_considered_in_odr': points_considered_in_odr}

def run_odr_revised(x_data_odr, y_data_odr, mask_type):
    (vars_overlap, vars_original, init_data_with_regime) = mask_overlap_revised(x_data_odr, y_data_odr, extra_mask=mask_type)

    xvar_overlap = vars_overlap[0]
    yvar_overlap = vars_overlap[1]
    xvar_init_data = vars_original[0]
    yvar_init_data = vars_original[1]
    xvar_init_data_with_regime = init_data_with_regime[0]
    yvar_init_data_with_regime = init_data_with_regime[1]
    # print('len(vars_overlap[0]): ', len(vars_overlap[0]))
    # #print('vars_overlap[0]: ', vars_overlap[0].to_numpy())
    # print('len(vars_overlap[1]): ', len(vars_overlap[1]))
    # #print('vars_overlap[1]: ', vars_overlap[1].to_numpy())
    # print('len(vars_original[0]): ', len(vars_original[0]))
    # print('len(vars_original[1]): ', len(vars_original[1]))

    #Get count of how many points are overlapping, and used in the ODR fitting
    overlap_points_count = len(vars_overlap[0])

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

    xvals_corrected_overlap = odr_intercept + odr_slope * xvar_overlap
    #gives our predicted y value for applying tracer to only regime
    xvals_corrected_with_regime = odr_intercept + odr_slope * xvar_init_data_with_regime

    #Because the ODR fitting might produce negative values, we have to exclude them from our metrics (unusable)
    valid_mask = (~np.isnan(xvals_corrected_overlap)) & (xvals_corrected_overlap >= 0) & (~yvar_overlap.isna())
    xvals_corrected_with_regime_no_negs = xvals_corrected_with_regime.mask(xvals_corrected_with_regime < 0, np.nan)


    # Apply mask
    y_pred_valid = xvals_corrected_overlap[valid_mask]
    y_overlap_valid = yvar_overlap[valid_mask]

    #Metrics
    rmse = np.sqrt(np.mean((y_overlap_valid - y_pred_valid)**2))
    #rmse= root_mean_squared_error(yvar_overlap, xvals_corrected_overlap)
    norm_rmse = rmse / np.std(y_overlap_valid)

    ss_res = np.sum((y_overlap_valid - (y_pred_valid))**2)
    ss_tot = np.sum((y_overlap_valid - np.mean(y_overlap_valid))**2)
    r2 = 1 - ss_res/ss_tot
    #r2 = r2_score(yvar_overlap, xvals_corrected_overlap)

    gap_fill_potential, gaps_needed_to_fill_for_regime = gap_fill_metrics(xvals_corrected_with_regime, yvar_init_data, extra_mask=mask_type)

    return {'xvar_overlap': xvar_overlap,
            'yvar_overlap': yvar_overlap,
            'slope': odr_slope,
            'intercept': odr_intercept,
            'rmse': rmse,
            'norm_rmse': norm_rmse,
            'r2': r2,
            'overlap_points_counted': overlap_points_count,
            'xvals_corrected_all_points': xvals_corrected_with_regime,
            'xvals_corrected_all_points_no_negs': xvals_corrected_with_regime_no_negs,
            'xvals_corrected_overlap': xvals_corrected_overlap,
            'xvals_corrected_overlap_no_negs': y_pred_valid,
            'y_overlap_valid_no_negs': y_overlap_valid,
            'gap_fill_potential': gap_fill_potential,
            'gaps_needed_to_fill_for_regime': gaps_needed_to_fill_for_regime}
def odr_improves(odr_init_outputs, odr_corr_outputs): 
    (xvar_overlap_i, yvar_overlap_i, xvar_init_data_i, 
    yvar_init_data_i, slope_i, intercept_i, correction_eq_i, 
    rmse_i, norm_rmse_i, r2_i, xvals_corrected_i, xvals_corrected_overlap_i,
    points_considered_in_odr_i) = odr_init_outputs.values()

    (xvar_overlap_c, yvar_overlap_c, xvar_init_data_c, 
    yvar_init_data_c, slope_c, intercept_c, 
    correction_eq_c, rmse_c, norm_rmse_c, r2_c, 
    xvals_corrected_c, xvals_corrected_overlap_c,
    points_considered_in_odr_c) = odr_corr_outputs.values()

    did_slope_improve = bool(abs(1 - slope_c) < abs(1 - slope_i))
    slope_distance_i = abs(1 - slope_i)
    slope_distance_c = abs(1 - slope_c)
    slope_distance_improvement_val = slope_distance_i - slope_distance_c

    intercept_err_i = abs(intercept_i) / np.nanmean(yvar_overlap_i)
    intercept_err_c = abs(intercept_c) / np.nanmean(yvar_overlap_c)

    did_rmse_improve = bool(rmse_c < rmse_i)
    rmse_percent_improvement = 100*((rmse_i - rmse_c) / (rmse_i))

    r2_round_init = round(r2_i, 2)
    r2_at_least_half_init = bool(r2_round_init >= 0.5)
    r2_round_corr = round(r2_i, 2)
    r2_at_least_half_corr = bool(r2_round_corr >= 0.5)

    #points_considered_in_odr_init should be same as points_considered_in_odr_corr

    return {'did_slope_improve': did_slope_improve, 
            'slope_distance_i': slope_distance_i, 
            'slope_distance_c': slope_distance_c, 
            'slope_distance_improvement_val': slope_distance_improvement_val, 
            'intercept_err_i': intercept_err_i, 
            'intercept_err_c': intercept_err_c, 
            'did_rmse_improve': did_rmse_improve, 
            'rmse_percent_improvement': rmse_percent_improvement, 
            }


    # slope_i, intercept_i, rmse_i, r2_i = metrics_init_input
    # slope_c, intercept_c, rmse_c, r2_c = metrics_corr_input

    # did_slope_improve = bool(abs(1 - slope_c) < abs(1 - slope_i))
    # slope_distance_i = abs(1 - slope_i)
    # slope_distance_c = abs(1 - slope_c)
    # slope_distance_improvement_val = slope_distance_i - slope_distance_c

    # intercept_err_i = abs(intercept_i) / np.nanmean(y_overlap)
    # intercept_err_c = abs(intercept_c) / np.nanmean(y_overlap)

    # did_rmse_improve = bool(rmse_c < rmse_i)
    # rmse_percent_improvement = 100*((rmse_i - rmse_c) / (rmse_i))
    # rmse_norm_i = rmse_i/np.nanmean(y_overlap)
    # rmse_norm_c = rmse_c/np.nanmean(y_overlap)

    # r2_round_init = round(r2_i, 2)
    # r2_at_least_half = bool(r2_round_init >= 0.5)

    # score_i = slope_distance_i + intercept_err_i + rmse_norm_i
    # score_c = slope_distance_c + intercept_err_c + rmse_norm_c
    # #Apply correction if initial score is higher than corrected score
    # apply_correction_based_off_score = bool(score_c < score_i)

    # #points_considered_in_odr_init should be same as points_considered_in_odr_corr

    # return did_slope_improve, slope_distance_i, slope_distance_c, slope_distance_improvement_val, intercept_err_i, intercept_err_c, did_rmse_improve, rmse_percent_improvement, rmse_norm_i, rmse_norm_c, r2_at_least_half, score_i, score_c, apply_correction_based_off_score
def process_one_voc_stage1(var_species_name, duplicate_species_name, var_species_data, duplicate_species_data):
    cache_file_stage1 = f'{CACHE_DIR}/stage1/{var_species_name}_duplicate_odr_fit.pkl'
    
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
            'metrics_init_slope': np.nan,
            'metrics_init_intercept': np.nan,
            'metrics_init_rmse': np.nan,
            'metrics_init_norm_rmse': np.nan,
            'metrics_init_r2': np.nan,
            'metrics_corr_slope': np.nan,
            'metrics_corr_intercept': np.nan,
            'metrics_corr_rmse': np.nan,
            'metrics_corr_norm_rmse': np.nan,
            'metrics_corr_r2': np.nan,
            'did_slope_improve': np.nan, 
            'slope_distance_from_1_init': np.nan,
            'slope_distance_from_1_corr': np.nan,
            'slope_distance_improvement_val': np.nan,
            'intercept_error_init': np.nan,
            'intercept_error_corr': np.nan,
            'did_rmse_improve': np.nan, 
            'rmse_percent_improvement': np.nan,
            'var_data_init_full': var_species_data,
            'var_data_init_overlap': pd.Series(np.nan, index=df_all_measured_species.index),
            'dup_data_init_full': duplicate_species_data,
            'dup_data_init_overlap': pd.Series(np.nan, index=df_all_measured_species.index),
            'dup_data_corr_full': pd.Series(np.nan, index=df_all_measured_species.index),
            'dup_data_corr_overlap': pd.Series(np.nan, index=df_all_measured_species.index)
        }
    elif isinstance(duplicate_species_name, list):
        print('Processing ', var_species_name, ': Has duplicate list')
        #print('duplicate_species_data: ', duplicate_species_data)

        #Run ODR first time
        run_odr_outputs_init = [run_odr(x_data_odr = data, 
                                  y_data_odr = var_species_data) for data in duplicate_species_data]
        # print('run_odr_outputs_init:', type(run_odr_outputs_init[0]))
        # for i, out_init in enumerate(run_odr_outputs_init):
        #     print(i, len(out_init))
    
        #Gives us initial: xvar_overlap, yvar_overlap, xvar_init_data, 
        # yvar_init_data, odr_slope, odr_intercept, correction_eq, 
        # rmse, norm_rmse, r2, xvals_corrected, 
        # points_considered_in_odr
        run_odr_outputs_corr = [run_odr(out_init['xvals_corrected'], out_init['yvar_overlap']) for out_init
                                in run_odr_outputs_init]
        #run_odr_outputs_corr = [run_odr(xvals_corrected, yvar_overlap) for xvar_overlap, yvar_overlap, 
                                # xvar_init_data, yvar_init_data, odr_slope, odr_intercept, correction_eq, 
                                # rmse, norm_rmse, r2, xvals_corrected, points_considered_in_odr
                                # in run_odr_outputs_init]
            #run_odr_outputs_init[10], 
                                  #y_data_odr = run_odr_outputs_init[ord][1])

        # #Quality control decision based off of if improvements are made by applying correction to UDAQ data
        odr_improves_out = [odr_improves(odr_init_outputs = out_init, 
                                     odr_corr_outputs = out_corr) for (out_init, out_corr)
                                     in (run_odr_outputs_init, run_odr_outputs_corr)]

        # print('type(run_odr_outputs_init):', type(run_odr_outputs_init))
        # print('run_odr_outputs_init:', run_odr_outputs_init)
        # #Quality control decision based off of if improvements are made by applying correction to UDAQ data

        stage1_results = {
            'x_species_name': duplicate_species_name,
            'y_species_name': var_species_name,
            'case': 'Has duplicate list, has duplicate data not all NaNs, has vardata not all NaNs',
            'points_considered_in_odr': [odr_results['points_considered_in_odr'] for odr_results in run_odr_outputs_init],
            'odr_eq_adj': [odr_results['correction_eq'] for odr_results in run_odr_outputs_init], #equation used to correct the UDAQ data
            'metrics_init_slope': [odr_results['odr_slope'] for odr_results in run_odr_outputs_init],
            'metrics_init_intercept': [odr_results['odr_intercept'] for odr_results in run_odr_outputs_init],
            'metrics_init_rmse': [odr_results['rmse'] for odr_results in run_odr_outputs_init],
            'metrics_init_norm_rmse': [odr_results['norm_rmse'] for odr_results in run_odr_outputs_init],
            'metrics_init_r2': [odr_results['r2'] for odr_results in run_odr_outputs_init],
            'metrics_corr_slope': [odr_results['odr_slope'] for odr_results in run_odr_outputs_corr],
            'metrics_corr_intercept': [odr_results['odr_intercept'] for odr_results in run_odr_outputs_corr],
            'metrics_corr_rmse': [odr_results['rmse'] for odr_results in run_odr_outputs_corr],
            'metrics_corr_norm_rmse': [odr_results['norm_rmse'] for odr_results in run_odr_outputs_corr],
            'metrics_corr_r2': [odr_results['r2'] for odr_results in run_odr_outputs_corr],
            'did_slope_improve': [odr_results['did_slope_improve'] for odr_results in odr_improves_out],
            'slope_distance_from_1_init': [odr_results['slope_distance_i'] for odr_results in odr_improves_out],
            'slope_distance_from_1_corr': [odr_results['slope_distance_c'] for odr_results in odr_improves_out],
            'slope_distance_improvement_val': [odr_results['slope_distance_improvement_val'] for odr_results in odr_improves_out],
            'intercept_error_init': [odr_results['intercept_err_i'] for odr_results in odr_improves_out],
            'intercept_error_corr': [odr_results['intercept_err_c'] for odr_results in odr_improves_out],
            'did_rmse_improve': [odr_results['did_rmse_improve'] for odr_results in odr_improves_out],
            'rmse_percent_improvement': [odr_results['rmse_percent_improvement'] for odr_results in odr_improves_out],
            'var_data_init_full': var_species_data,
            'var_data_init_overlap': [odr_results['yvar_overlap'] for odr_results in run_odr_outputs_init],
            'dup_data_init_full': duplicate_species_data,
            'dup_data_init_overlap': [odr_results['xvar_overlap'] for odr_results in run_odr_outputs_init],
            'dup_data_corr_full': [odr_results['xvals_corrected'] for odr_results in run_odr_outputs_init],
            'dup_data_corr_overlap': [odr_results['xvals_corrected_overlap'] for odr_results in run_odr_outputs_init]
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
            'metrics_init_slope': np.nan,
            'metrics_init_intercept': np.nan,
            'metrics_init_rmse': np.nan,
            'metrics_init_norm_rmse': np.nan,
            'metrics_init_r2': np.nan,
            'metrics_corr_slope': np.nan,
            'metrics_corr_intercept': np.nan,
            'metrics_corr_rmse': np.nan,
            'metrics_corr_norm_rmse': np.nan,
            'metrics_corr_r2': np.nan,
            'did_slope_improve': np.nan, 
            'slope_distance_from_1_init': np.nan,
            'slope_distance_from_1_corr': np.nan,
            'slope_distance_improvement_val': np.nan,
            'intercept_error_init': np.nan,
            'intercept_error_corr': np.nan,
            'did_rmse_improve': np.nan, 
            'rmse_percent_improvement': np.nan,
            'var_data_init_full': var_species_data,
            'var_data_init_overlap': pd.Series(np.nan, index=df_all_measured_species.index),
            'dup_data_init_full': duplicate_species_data,
            'dup_data_init_overlap': pd.Series(np.nan, index=df_all_measured_species.index),
            'dup_data_corr_full': pd.Series(np.nan, index=df_all_measured_species.index),
            'dup_data_corr_overlap': pd.Series(np.nan, index=df_all_measured_species.index)
        }
    else:
        print('Processing ', var_species_name, ': Has duplicate, has duplicate data not all NaNs, has vardata not all NaNs')

        #Run ODR first time
        run_odr_outputs_init = run_odr(x_data_odr = duplicate_species_data, 
                                  y_data_odr = var_species_data)
        
        #Gives us initial: xvar_overlap, yvar_overlap, xvar_init_data, 
        # yvar_init_data, odr_slope, odr_intercept, correction_eq, 
        # rmse, norm_rmse, r2, xvals_corrected,
        # points_considered_in_odr
        run_odr_outputs_corr = run_odr(x_data_odr = run_odr_outputs_init['xvals_corrected'], 
                                        y_data_odr = run_odr_outputs_init['yvar_overlap'])

        #Quality control decision based off of if improvements are made by applying correction to UDAQ data
        odr_improves_out = odr_improves(odr_init_outputs = run_odr_outputs_init, 
                                       odr_corr_outputs = run_odr_outputs_corr)
        # print(odr_improves_out)

        stage1_results = {
            'x_species_name': duplicate_species_name,
            'y_species_name': var_species_name,
            'case': 'Has duplicate, has duplicate data not all NaNs, has vardata not all NaNs',
            'points_considered_in_odr': run_odr_outputs_init['points_considered_in_odr'],
            'odr_eq_adj': run_odr_outputs_init['correction_eq'], #equation used to correct the UDAQ data
            'metrics_init_slope': run_odr_outputs_init['odr_slope'],
            'metrics_init_intercept': run_odr_outputs_init['odr_intercept'],
            'metrics_init_rmse': run_odr_outputs_init['rmse'],
            'metrics_init_norm_rmse': run_odr_outputs_init['norm_rmse'],
            'metrics_init_r2': run_odr_outputs_init['r2'],
            'metrics_corr_slope': run_odr_outputs_corr['odr_slope'],
            'metrics_corr_intercept': run_odr_outputs_corr['odr_intercept'],
            'metrics_corr_rmse': run_odr_outputs_corr['rmse'],
            'metrics_corr_norm_rmse': run_odr_outputs_corr['norm_rmse'],
            'metrics_corr_r2': run_odr_outputs_corr['r2'],
            'did_slope_improve': odr_improves_out['did_slope_improve'], 
            'slope_distance_from_1_init': odr_improves_out['slope_distance_i'],
            'slope_distance_from_1_corr': odr_improves_out['slope_distance_c'],
            'slope_distance_improvement_val': odr_improves_out['slope_distance_improvement_val'],
            'intercept_error_init': odr_improves_out['intercept_err_i'],
            'intercept_error_corr': odr_improves_out['intercept_err_c'],
            'did_rmse_improve': odr_improves_out['did_rmse_improve'],
            'rmse_percent_improvement': odr_improves_out['rmse_percent_improvement'],
            'var_data_init_full': var_species_data,
            'var_data_init_overlap': run_odr_outputs_init['yvar_overlap'],
            'dup_data_init_full': duplicate_species_data,
            'dup_data_init_overlap': run_odr_outputs_init['xvar_overlap'],
            'dup_data_corr_full': run_odr_outputs_init['xvals_corrected'],
            'dup_data_corr_overlap': run_odr_outputs_init['xvals_corrected_overlap']
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
        print('Processing Stage 1 for var: ', spec_name)
        dup_name = df_duplicates_and_tracers['Duplicate_name'].iloc[row_idx]

        #print('type(dup_name): ', type(dup_name))
        if pd.isna(dup_name):
            dup_species_data = pd.Series(np.nan, index=df_all_measured_species.index)
        elif pd.notna(dup_name) and ';' in dup_name:
            dup_name_list = dup_name.split('; ')
            print('Has more than one dup: ', dup_name_list)
            dup_name = dup_name_list
            dup_species_data = [df_all_measured_species[dup_name[0]], df_all_measured_species[dup_name[1]]]
            #print('dup_species_data type should be a list: ', type(dup_species_data))

        else:
            print('Processing Stage 1 with duplicate: ', dup_name)
            dup_species_data = df_all_measured_species[dup_name]

        stage1_results = process_one_voc_stage1(var_species_name = spec_name, 
                                                duplicate_species_name = dup_name, 
                                                var_species_data = df_all_measured_species[spec_name], 
                                                duplicate_species_data = dup_species_data)
        
        stage1_results_all.append(stage1_results)
    return stage1_results_all
def build_summary_stage1(stage1_results_all):
    summary_rows = []
    for r in stage1_results_all:
        print('Saving for :', r['y_species_name'])

        summary_rows.append({    
            'varname': r['y_species_name'],
            'Duplicate_name': r['x_species_name'],
            'case': r['case'],
            'points_considered_in_odr': r['points_considered_in_odr'],
            'odr_eq_adj': r['odr_eq_adj'],
            'metrics_init_slope': r['metrics_init_slope'],
            'metrics_init_intercept': r['metrics_init_intercept'],
            'metrics_init_rmse': r['metrics_init_rmse'],
            'metrics_init_norm_rmse': r['metrics_init_norm_rmse'],
            'metrics_init_r2': r['metrics_init_r2'],
            'metrics_corr_slope': r['metrics_corr_slope'],
            'metrics_corr_intercept': r['metrics_corr_intercept'],
            'metrics_corr_rmse': r['metrics_corr_rmse'],
            'metrics_corr_norm_rmse':r['metrics_corr_norm_rmse'],
            'metrics_corr_r2': r['metrics_corr_r2'],
            'did_slope_improve': r['did_slope_improve'],
            'slope_distance_from_1_init': r['slope_distance_from_1_init'],
            'slope_distance_from_1_corr': r[ 'slope_distance_from_1_corr'],
            'slope_distance_improvement_val': r['slope_distance_improvement_val'],
            'intercept_error_init': r['intercept_error_init'],
            'intercept_error_corr': r['intercept_error_corr'],
            'did_rmse_improve': r['did_rmse_improve'], 
            'rmse_percent_improvement': r['rmse_percent_improvement'],
        })
    stage1_summary_savepath = stage_data_dirs['stage1'] + 'csv/stage1_summary_new.csv'
    pd.DataFrame(summary_rows).to_csv(stage1_summary_savepath, index=True)
    print('Saved summary to: ' + str(stage1_summary_savepath))

###############################################################
def stage2_odr_fitting_for_tracer():
    corr_matrices = {}
    tracer_dfs = {}
    good_tracers_names = {}
    tracer_filtered_data = {}

    #Get overlap between tracers and the target species. 
    for target, tracer in species_tracer_dict.items():
        print('target: ', target)
        print('tracer length: ', len(tracer))
        print('all potential tracers: ', tracer)

        all_results_tracer = {}

        for tracer_length in range(0, len(tracer)):
            #Perform ODR fitting between the target species and the tracer
            all_results_tracer[tracer[tracer_length]] = {}
            #loop through all mask types (for regime)
            for regime, mask_applied in regime_mask.items():
                run_odr_outputs = run_odr_revised(x_data_odr = df_all_measured_species[tracer[tracer_length]],
                                                  y_data_odr = df_all_measured_species[target], 
                                                  mask_type = mask_applied)
                all_results_tracer[tracer[tracer_length]][regime] = run_odr_outputs
        #Makes a dictionary of dictionaries so you'd get odr fittings for 
        #all_results_tracer[tracer name][regime]
        #print('all_results_tracer: ', all_results_tracer)


        # Organize one dataframe for all the tracers per target
        rows = []
        for tracer, regimes in all_results_tracer.items():
            for regime, metrics in regimes.items():
                if metrics is None:
                    continue
                row = {"tracer": tracer,
                        "regime": regime,
                        **metrics}
                rows.append(row)
        odr_results_df = pd.DataFrame(rows)
        print('odr_results_df: ', odr_results_df)

        cache_file_stage2_dfs_all_tracers_fitted = f'{CACHE_DIR}/stage2/{target}_df_all_tracers_fitted.pkl'
        odr_results_df.to_pickle(cache_file_stage2_dfs_all_tracers_fitted)

def stage2_tracer_scoring(corr_threshold_tracers, norm_rmse_max, r2_min):
    for target in species_tracer_dict:
        print('target: ', target)

        cache_file_stage2_dfs_all_tracers_fitted = f'{CACHE_DIR}/stage2/{target}_df_all_tracers_fitted.pkl'
        with open(cache_file_stage2_dfs_all_tracers_fitted, 'rb') as f:
            odr_results_df = pickle.load(f)
        
        df_odr_scoring = odr_results_df.copy()
        #remove negative r squared values
        df_odr_scoring_filtered = df_odr_scoring[df_odr_scoring['r2'] >= 0]
        #remove any overlap less than 50
        df_odr_scoring_filtered = df_odr_scoring_filtered[df_odr_scoring_filtered['overlap_points_counted'] >= 25]
        df_odr_scoring_filtered = df_odr_scoring_filtered[df_odr_scoring_filtered['norm_rmse'] <= norm_rmse_max]
        df_odr_scoring_filtered = df_odr_scoring_filtered[df_odr_scoring_filtered['r2'] >= r2_min]
        df_odr_scoring_filtered = df_odr_scoring_filtered[df_odr_scoring_filtered['gap_fill_potential'] > 0]
        # top_10_norm_rmse = df_odr_scoring_filtered.nsmallest(10, 'norm_rmse')
        # top_10_overlap = df_odr_scoring_filtered.nlargest(10, 'overlap_points_counted')
        # #print('df_odr_scoring_filtered: ', df_odr_scoring_filtered)
        # #print('top_10_rmse: ', top_10_norm_rmse)
        # #print('top_10_overlap: ', top_10_overlap)
        for regime, mask_applied in regime_mask.items():
            print('regime: ', regime)
            #Sort tracers based on regime
            df_reg = df_odr_scoring_filtered[df_odr_scoring_filtered['regime'] == regime]
            print('df_reg: ', df_reg)
            #print(df_reg['norm_rmse'])
            # Top 10 by rmse
            top10_by_norm_rmse = df_reg.nsmallest(10, 'norm_rmse')

            top_tracer_names = top10_by_norm_rmse['tracer'].tolist()
            #print('top_tracer_names: ', top_tracer_names)
            #print('top10_by_norm_rmse: ', top10_by_norm_rmse[['tracer', 'regime', 'norm_rmse', 'r2', 'gap_fill_potential', 'gaps_needed_to_fill_for_regime', 'overlap_points_counted']])
            
            #Test for collinearity
            corr_matrix = pd.DataFrame(index=top_tracer_names, columns=top_tracer_names, dtype=float)

            for tracer1, tracer2 in combinations(top_tracer_names, 2):
                print('Mask overlap')
                (vars_overlap, vars_original, init_data_with_regime) = mask_overlap_revised(df_all_measured_species[tracer1], 
                                                                     df_all_measured_species[tracer2], 
                                                                     extra_mask = mask_applied)
                print('Computing correlation')
                # compute correlation
                corr = vars_overlap[0].corr(vars_overlap[1])
                #print(f"{tracer1} vs {tracer2}: correlation = {corr:.3f}")
            
                corr_matrix.loc[tracer1, tracer2] = corr
                corr_matrix.loc[tracer2, tracer1] = corr  # symmetric

            # Fill diagonal with 1
            np.fill_diagonal(corr_matrix.values, 1.0)
            #print('corr_matrix: ', corr_matrix)

            selected_tracers = []
            for top_tracer in top_tracer_names:
                # Check if this tracer is collinear with any already selected
                if all(abs(corr_matrix.loc[top_tracer, kept]) <= corr_threshold_tracers for kept in selected_tracers):
                    selected_tracers.append(top_tracer)
            print("Selected tracers after removing collinear ones:", selected_tracers)
            
            # Filter dataframe
            df_filtered_collinearity = top10_by_norm_rmse[top10_by_norm_rmse['tracer'].isin(selected_tracers)]

            print('Selected tracers: ', df_filtered_collinearity)

            #save pickle for top tracers per regime
            cache_file_name_top_10_per_regime = f'{CACHE_DIR}/stage2/{target}_top10_tracers_{regime}_stage2.pkl'
            df_filtered_collinearity.to_pickle(cache_file_name_top_10_per_regime)

            # #######################################

            # for selected_spec in selected_tracers:
            #     print('selected_spec: ', selected_spec)
            #     print('regime: ', regime)
            #     selected_tracer_position = selected_tracers.index(selected_spec)
            #     print('selected tracer data:')
            #     print(df_reg['xvar_overlap'].values[selected_tracer_position])

            #     # search_term = (df_odr_scoring_filtered['tracer'] == selected_spec) & (df_odr_scoring_filtered['regime'] == regime)
            #     # if search_term.any():
            #     #     print(selected_spec, 'overlap: ', df_odr_scoring_filtered.loc[search_term,'overlap_points_counted'].values)


                
                # ######################################

def stage2_savedata_to_csv():
    regime_list = ['day', 'night']
    for reg in regime_list:
        target_tracer_match_out_list = []
        for row_idx in range(0, len(df_duplicates_and_tracers.index)):
            target = df_duplicates_and_tracers['Varname'].iloc[row_idx]
            cache_file_name_top_10_per_regime = f'{CACHE_DIR}/stage2/{target}_top10_tracers_{reg}_stage2.pkl'
            df_top_species_per_regime = pd.read_pickle(cache_file_name_top_10_per_regime)
            #print('pickle load df_top_species_per_regime: ', df_top_species_per_regime)
            top_tracers = df_top_species_per_regime['tracer'].tolist()
            norm_rmse_vals =  df_top_species_per_regime['norm_rmse'].tolist()
            r2_vals = df_top_species_per_regime['r2'].tolist()
            gap_fill_potential_vals = df_top_species_per_regime['gap_fill_potential'].tolist()
            gaps_needed_to_fill_for_regime_vals = df_top_species_per_regime['gaps_needed_to_fill_for_regime'].tolist()
            overlap_points_counted_vals = df_top_species_per_regime['overlap_points_counted'].tolist()

            tracer_match = {'target': target, 'tracers': top_tracers,
                            'norm_rmse': norm_rmse_vals, 'r2': r2_vals,
                            'gap_fill_potential': gap_fill_potential_vals,
                            'gaps_needed_to_fill_for_regime': gaps_needed_to_fill_for_regime_vals,
                            'overlap_points_counted': overlap_points_counted_vals}
            target_tracer_match_out_list.append(tracer_match)
        df_target_tracer_match=pd.DataFrame(target_tracer_match_out_list)
        df_target_tracer_match.to_csv(stage_data_dirs['stage3'] + f"csv/{reg}_target_tracer_match.csv", index=True)

#######################################
#Plots to perform quality control check: determine if ODR fit was helpful to correct UDAQ data
def scatterplots_for_comparing_init_and_corr_odr_fits(odr_fit_results, stage_type, flag_type):
    if flag_type == 'Duplicate list':
        for dup in range(0, 1):
            x_voc_name = odr_fit_results['x_species_name'][dup]
            y_voc_name = odr_fit_results['y_species_name']

            x_init_overlap = odr_fit_results['dup_data_init_overlap'][dup]
            x_corr_overlap = odr_fit_results['dup_data_corr_overlap'][dup]
            y_overlap = odr_fit_results['var_data_init_overlap'][dup]

            slope_init = odr_fit_results['metrics_init_slope'][dup]
            slope_corr = odr_fit_results['metrics_corr_slope'][dup]

            intercept_init = odr_fit_results['metrics_init_intercept'][dup]
            intercept_corr = odr_fit_results['metrics_corr_intercept'][dup]

            rmse_init = odr_fit_results['metrics_init_rmse'][dup]
            rmse_corr = odr_fit_results['metrics_corr_rmse'][dup]

            r2_init = odr_fit_results['metrics_init_r2'][dup]
            r2_corr = odr_fit_results['metrics_corr_r2'][dup]

            #scatterplot of ML vs UDAQ initial
            fig, ax = plt.subplots(1, 2, figsize=(10,10), tight_layout=True)
            ax[0].scatter(x_init_overlap, y_overlap, s=10, alpha=0.5)

            #To draw regression line, we need a continuous line. Since some species are in ppt, we need
            #to select an appropriate scale to the step
            if np.nanmax(x_init_overlap) < 0.1:
                step = 0.001
            else:
                step = 0.1

            xrange_init = np.arange(0,np.nanmax(x_init_overlap), step)

            ax[0].plot(xrange_init, (slope_init * xrange_init + intercept_init))
            ax[0].set_xlabel(x_voc_name + ' (ppb)')
            ax[0].set_title('Initial')
            # if dup == 0:
            #     label0 = odr_fit_results['x_species_name'][0]
            # elif dup == 1:
            #     label0 = odr_fit_results['x_species_name'][1]

            ax[0].set_ylabel(str(y_voc_name) + ' (ppb)')
            ax[0].set_ylim([0, np.nanmax(y_overlap)*1.1])
            ax[0].set_xlim([0, np.nanmax(x_init_overlap)*1.1])


            ax[0].text(0.05, 0.96, "Slope = " + str(round(slope_init, 3)), transform=ax[0].transAxes)
            ax[0].text(0.05, 0.94, "Intercept = " + str(round(intercept_init, 3)), transform=ax[0].transAxes)
            ax[0].text(0.05, 0.92, "R$^2$= " + str(round(r2_init, 3)), transform=ax[0].transAxes) 
            ax[0].text(0.05, 0.90, "RMSE:  " + str(round(rmse_init, 3)), transform=ax[0].transAxes)


            ax[1].scatter(x_corr_overlap, y_overlap, s=10, alpha=0.5)
            ax[1].set_xlabel(x_voc_name + ' (ppb)')
            if np.nanmax(x_corr_overlap) < 0.1:
                step = 0.001
            else:
                step = 0.1

            xrange_corr = np.arange(0,np.nanmax(x_corr_overlap),step)
            
            ax[1].plot(xrange_corr, (slope_corr * xrange_corr + intercept_corr))
            ax[1].set_title('Corrected')
            ax[1].set_ylabel(y_voc_name + ' (ppb)')
            ax[1].set_ylim([0, np.nanmax(y_overlap)*1.1])
            ax[1].set_xlim([0, np.nanmax(x_init_overlap)*1.1])

            ax[1].text(0.05, 0.96, "Slope = " + str(round(slope_corr, 3)), transform=ax[1].transAxes)
            ax[1].text(0.05, 0.94, "Intercept = " + str(round(intercept_corr, 3)), transform=ax[1].transAxes)
            ax[1].text(0.05, 0.92, "R$^2$= " + str(round(r2_corr, 3)), transform=ax[1].transAxes) 
            ax[1].text(0.05, 0.90, "RMSE:  " + str(round(rmse_corr, 3)), transform=ax[1].transAxes)

            plt.savefig(stage_data_dirs[stage_type] + f'plots/init_corr_scatterplot_comparison_{y_voc_name}_{x_voc_name} + .png', dpi = 150)
            plt.show()
            plt.close()
    else:    
        x_voc_name = odr_fit_results['x_species_name']
        y_voc_name = odr_fit_results['y_species_name']

        x_init_overlap = odr_fit_results['dup_data_init_overlap']
        x_corr_overlap = odr_fit_results['dup_data_corr_overlap']
        y_overlap = odr_fit_results['var_data_init_overlap']

        slope_init = odr_fit_results['metrics_init_slope']
        slope_corr = odr_fit_results['metrics_corr_slope']

        intercept_init = odr_fit_results['metrics_init_intercept']
        intercept_corr = odr_fit_results['metrics_corr_intercept']

        rmse_init = odr_fit_results['metrics_init_rmse']
        rmse_corr = odr_fit_results['metrics_corr_rmse']

        r2_init = odr_fit_results['metrics_init_r2']
        r2_corr = odr_fit_results['metrics_corr_r2']

        fig, ax = plt.subplots(1, 2, figsize=(10,10), tight_layout=True)
        ax[0].scatter(x_init_overlap, y_overlap, s=10, alpha=0.5)

        #To draw regression line, we need a continuous line. Since some species are in ppt, we need
        #to select an appropriate scale to the step
        if np.nanmax(x_init_overlap) < 0.1:
            step = 0.001
        elif x_voc_name == 'H2O_Piccaro' or y_voc_name == 'H2O_Piccaro':
            step = 1000
        else:
            step = 0.1

        xrange_init = np.arange(0,np.nanmax(x_init_overlap), step)
        

        ax[0].plot(xrange_init, (slope_init * xrange_init + intercept_init))
        ax[0].set_title('Initial')
        ax[0].set_xlabel(x_voc_name + ' (ppb)')
        ax[0].set_ylabel(y_voc_name + ' (ppb)')
        ax[0].set_ylim([0, np.nanmax(y_overlap)])
        ax[0].set_xlim([0, np.nanmax(x_init_overlap)])

        ax[0].text(0.05, 0.96, "Slope = " + str(round(slope_init, 3)), transform=ax[0].transAxes)
        ax[0].text(0.05, 0.94, "Intercept = " + str(round(intercept_init, 3)), transform=ax[0].transAxes)
        ax[0].text(0.05, 0.92, "R$^2$= " + str(round(r2_init, 3)), transform=ax[0].transAxes) 
        ax[0].text(0.05, 0.90, "RMSE:  " + str(round(rmse_init, 3)), transform=ax[0].transAxes)


        ax[1].scatter(x_corr_overlap, y_overlap, s=10, alpha=0.5)

        if np.nanmax(x_corr_overlap) < 0.1:
            step = 0.001
        elif x_voc_name == 'H2O_Piccaro' or y_voc_name == 'H2O_Piccaro':
            step = 1000
        else:
            step = 0.1

        xrange_corr = np.arange(0,np.nanmax(x_corr_overlap),step)
        
        ax[1].plot(xrange_corr, (slope_corr * xrange_corr + intercept_corr))
        ax[1].set_title('Corrected')
        ax[1].set_xlabel(str(x_voc_name) + ' (ppb)')
        ax[1].set_ylabel(str(y_voc_name) + ' (ppb)')
        ax[1].set_ylim([0, np.nanmax(y_overlap)])
        ax[1].set_xlim([0, np.nanmax(x_init_overlap)])

        ax[1].text(0.05, 0.96, "Slope = " + str(round(slope_corr, 3)), transform=ax[1].transAxes)
        ax[1].text(0.05, 0.94, "Intercept = " + str(round(intercept_corr, 3)), transform=ax[1].transAxes)
        ax[1].text(0.05, 0.92, "R$^2$= " + str(round(r2_corr, 3)), transform=ax[1].transAxes) 
        ax[1].text(0.05, 0.90, "RMSE:  " + str(round(rmse_corr, 3)), transform=ax[1].transAxes)

        plt.savefig(stage_data_dirs[stage_type] + 'plots/init_corr_scatterplot_comparison_' + str(y_voc_name) + '_' +str(x_voc_name) + '.png', dpi = 150)
        plt.show()
        plt.close()
def tracer_scatterplots(regime, xdata, ydata, ypred, ypred_no_negs, yvar_no_negs, xname, yname, slope_plot, intercept_plot, r2_plot, norm_rmse_plot, gap_fill_potential_plot, gap_needed_to_fill_plot):
    yline = np.linspace(min(ydata), max(ydata), 100)
    
    mask_negs = ypred <= 0

    ypred_negs = ypred[mask_negs]
    ydata_negs = ydata[mask_negs]
    
    plt.figure(figsize=(7,7), tight_layout = True)
    plt.scatter(ypred_no_negs, yvar_no_negs, color = 'tab:blue', alpha=0.5)
    plt.scatter(ypred_negs, ydata_negs, color = 'tab:gray', alpha=0.5)
    one_to_one_line = plt.plot(yline, yline, color = 'tab:green', linestyle = '--', label = '1:1')
    
    plt.text(0.03, 0.88,"R$^2$= " + str(round(r2_plot, 3)) + '\n' +  "Norm RMSE = " + str(round(norm_rmse_plot, 3)) + '\n' + f'Tracer fills {gap_fill_potential_plot}/{gap_needed_to_fill_plot} gaps', transform=plt.gca().transAxes, verticalalignment='top')
    plt.legend(loc = 'upper left')
    plt.xlabel(f'Predicted {yname} \n based on '+ xname + ' tracer (ppb)')
    plt.ylabel(yname + ' (ppb)')
    plt.title("Tracer relationship (" +  regime + ")")
    plt.savefig(stage_data_dirs['stage2'] + 'plots/scatterplots/tracer_scatterplot_comparison_' + yname + '_' + xname + '_' + regime + '.png', dpi = 150)
    plt.show()
    
#######################################
def gap_counter(data_before_change, data_after_change):
    #gap counter
    before_change_gap_count = data_before_change.isna().sum()
    after_change_gap_count = data_after_change.isna().sum()
    
    number_gaps_filled_by_change = before_change_gap_count - after_change_gap_count
    percent_of_gaps_filled_by_change = ((before_change_gap_count - after_change_gap_count) / before_change_gap_count) * 100

    return before_change_gap_count, after_change_gap_count, number_gaps_filled_by_change, percent_of_gaps_filled_by_change
def gap_filling(x_data_full, y_data_full):
    #Fill NaNs in ML Data with corrected UDAQ data
    filled_voc = y_data_full.fillna(x_data_full)
    return filled_voc
#def gaps_filled_data_only():
def stage3_dupfilling_one_voc(var_species_name):
    cache_file_stage1 = f'{CACHE_DIR}/stage1/{var_species_name}_duplicate_odr_fit.pkl'
    cache_file_stage3 = f'{CACHE_DIR}/stage3/{var_species_name}_duplicate_filled.pkl'

    # ---- load if already processed ----
    # if os.path.exists(cache_file_stage3):
    #     with open(cache_file_stage3, 'rb') as f:
    #         return pickle.load(f)
    # else:
    with open(cache_file_stage1, 'rb') as f:
        r = pickle.load(f)
    
    var_name = r['y_species_name']
    print('Filling ' + var_name + ' in Stage 3')

    dup_name = r['x_species_name']
    #print(dup_name)
    dup_init_data_full = r['dup_data_init_full']
    dup_corr_data_full = r['dup_data_corr_full']
    var_init_data_full = r['var_data_init_full']
    #print('var_init_data_full: ', var_init_data_full)
    r_squared_init = r['metrics_init_r2']
    #print('r_squared_init: ', r_squared_init)
    r_squared_corr = r['metrics_corr_r2']
    rmse_init = r['metrics_init_rmse']
    rmse_corr = r['metrics_corr_rmse']
    did_rmse_improve = r['did_rmse_improve']
    odr_eq_adj = r['odr_eq_adj']

    if isinstance(dup_name, list):
        threshold = 0.45
        # store decisions
        y_data_input = None
        filled_voc = []
        filled_with_name = []
        fill_case1 = []
        odr_eq_adj_used = []
        rmse_used = []
        r2_used = []
        before_change_gap_count1 = []
        after_change_gap_count1 = []
        number_gaps_filled_by_change1 = []
        percent_of_gaps_filled_by_change1 = []

        # Determine evaluation order without modifying the lists
        if r_squared_init[0] >= r_squared_init[1]:
            order = [0, 1]
        else:
            order = [1, 0]

        for pos in order:
            init_val = r_squared_init[pos]
            corr_val = r_squared_corr[pos]
            # choose x data
            if init_val >= threshold and corr_val >= init_val:
                print(f"apply corrected value from corrected[{pos}] ({dup_name[pos]})")
                filled_with_name.append(dup_name[pos])
                fill_case1.append('Filled with duplicate corrected')
                odr_eq_adj_used.append(odr_eq_adj[pos])
                rmse_used.append(rmse_corr[pos])
                r2_used.append(r_squared_corr[pos])
                x_data_full = (dup_corr_data_full[pos]).mask(dup_corr_data_full[pos] < 0, np.nan)
                #print('x_data_full: ', x_data_full)
                
            elif init_val >= threshold and corr_val <= init_val:
                print(f"apply initial value from initial[{pos}] ({dup_name[pos]})")
                filled_with_name.append(dup_name[pos])
                fill_case1.append('Filled with duplicate initial')
                odr_eq_adj_used.append(np.nan)
                rmse_used.append(rmse_init[pos])
                r2_used.append(r_squared_init[pos])
                x_data_full = dup_init_data_full[pos]
                #print('x_data_full: ', x_data_full)

            elif init_val <= threshold and corr_val >= threshold:
                print(f"apply corrected value from corrected[{pos}] ({dup_name[pos]})")
                filled_with_name.append(dup_name[pos])
                fill_case1.append('Filled with duplicate corrected')
                odr_eq_adj_used.append(odr_eq_adj[pos])
                rmse_used.append(rmse_corr[pos])
                r2_used.append(r_squared_corr[pos])
                x_data_full = (dup_corr_data_full[pos]).mask(dup_corr_data_full[pos] < 0, np.nan)


            elif init_val <= threshold and corr_val <= threshold:
                print(f"neither initial or corrected value from initial[{pos}] or corrected[{pos}] applied ({dup_name[pos]})")
                filled_with_name.append(np.nan)
                fill_case1.append('No fill from dup due to low r squared')
                odr_eq_adj_used.append(np.nan)
                rmse_used.append(np.nan)
                r2_used.append(np.nan)
                x_data_full = var_init_data_full
                print('x_data_full: ', x_data_full)

            # choose y input
            print('y_data_input: ', y_data_input)
            if y_data_input is None:
                y_data_input = var_init_data_full   #first pass uses original y
                #print('y_data_input2: ', y_data_input)
            else:
                y_data_input = filled_voc_partial  #carry forward
                #print('y_data_input2: ', y_data_input)

            #print('type x_data_full: ', type(x_data_full), '\n type y_data_input: ', type(y_data_input))

            # apply function
            filled_voc_partial = gap_filling(x_data_full, y_data_full = y_data_input)
            #print('filled_voc_partial: ', filled_voc_partial)

            (before_change_gap_count1_partial, after_change_gap_count1_partial, 
            number_gaps_filled_by_change1_partial, percent_of_gaps_filled_by_change1_partial) = gap_counter(data_before_change = y_data_input,
                                                                                            data_after_change = filled_voc_partial)
            #store results
            filled_voc.append(filled_voc_partial)

            before_change_gap_count1.append(before_change_gap_count1_partial)
            after_change_gap_count1.append(after_change_gap_count1_partial)
            number_gaps_filled_by_change1.append(number_gaps_filled_by_change1_partial)
            percent_of_gaps_filled_by_change1.append(percent_of_gaps_filled_by_change1_partial)
    
    #Case: Rsquared is good for initial, RMSE does not improve with correction
    elif (round(r_squared_init, 2) >= 0.45) & (rmse_init < rmse_corr):
        fill_case1 = 'Filled with duplicate initial'
        filled_with_name = dup_name + ' initial'
        print(fill_case1)
        odr_eq_adj_used = np.nan
        rmse_used = rmse_init
        r2_used = r_squared_init
        filled_voc = gap_filling(x_data_full = dup_init_data_full, 
                                y_data_full = var_init_data_full)

        (before_change_gap_count1, after_change_gap_count1, 
        number_gaps_filled_by_change1, percent_of_gaps_filled_by_change1) = gap_counter(data_before_change = var_init_data_full, 
                                                                                        data_after_change = filled_voc)
        
    #Case: Rsquared is good for corrected, RMSE improves with correction
    elif (round(r_squared_corr, 2) >= 0.45) & (rmse_init > rmse_corr):
        fill_case1 = 'Filled with duplicate corrected'
        filled_with_name = dup_name + ' corrected'
        print(fill_case1)
        odr_eq_adj_used = odr_eq_adj
        rmse_used = rmse_corr
        r2_used = r_squared_corr
        filled_voc = gap_filling(x_data_full = (dup_corr_data_full.mask(dup_corr_data_full < 0, np.nan)), 
                                y_data_full = var_init_data_full)

        (before_change_gap_count1, after_change_gap_count1, 
        number_gaps_filled_by_change1, percent_of_gaps_filled_by_change1) = gap_counter(data_before_change = var_init_data_full, 
                                                                                        data_after_change = filled_voc)
    # #Case: Rsquared is good for corrected, RMSE improves with correction, R squared improves
    # elif (r_squared_corr >= 0.45) & (did_rmse_improve == 'True') & (r_squared_init < r_squared_corr):
    #     fill_case = 'Filled with duplicate corrected'
    #     print(fill_case)
    #     merged_voc = gap_filling(x_data_full = dup_corr_data_full, 
    #                             y_data_full = var_init_data_full)

    #     (before_change_gap_count, after_change_gap_count, 
    #     number_gaps_filled_by_change, percent_of_gaps_filled_by_change) = gap_counter(data_before_change = var_init_data_full, 
    #                                                                                     data_after_change = merged_voc)
    elif (r_squared_init <= 0.45):
        fill_case1 = 'No fill from dup due to low r squared'
        filled_with_name = np.nan
        odr_eq_adj_used = np.nan
        rmse_used = np.nan
        r2_used = np.nan
        filled_voc = gap_filling(x_data_full = var_init_data_full, 
                                y_data_full = var_init_data_full)
        (before_change_gap_count1, after_change_gap_count1, 
        number_gaps_filled_by_change1, percent_of_gaps_filled_by_change1) = gap_counter(data_before_change = var_init_data_full, 
                                                                                    data_after_change = filled_voc)
    elif dup_init_data_full.isna().all():
        fill_case1 = 'No dup so stays vars'
        filled_with_name = np.nan
        #print(fill_case)
        odr_eq_adj_used = np.nan
        rmse_used = np.nan
        r2_used = np.nan
        filled_voc = gap_filling(x_data_full = var_init_data_full, 
                                y_data_full = var_init_data_full)
        (before_change_gap_count1, after_change_gap_count1, 
        number_gaps_filled_by_change1, percent_of_gaps_filled_by_change1) = gap_counter(data_before_change = var_init_data_full, 
                                                                                    data_after_change = filled_voc)
    elif var_init_data_full.isna().all():
        fill_case1 = 'No vardata so filled with init dup'
        filled_with_name = np.nan
        #print(fill_case)
        odr_eq_adj_used = np.nan
        rmse_used = np.nan
        r2_used = np.nan
        filled_voc = gap_filling(x_data_full = dup_init_data_full, 
                                y_data_full = var_init_data_full)
        (before_change_gap_count1, after_change_gap_count1, 
        number_gaps_filled_by_change1, percent_of_gaps_filled_by_change1) = gap_counter(data_before_change = var_init_data_full, 
                                                                                    data_after_change = filled_voc)
    else:
        print('Failed case')

    stage3_dupfilling_results = {
        'varname': var_name,
        'potential_duplicates': dup_name,
        'filled_with_name': filled_with_name,
        'fill_case1': fill_case1,
        'odr_eq_adj_used': odr_eq_adj_used,
        'rmse_used': rmse_used,
        'r2_used': r2_used,
        'before_change_gap_count1': before_change_gap_count1, 
        'after_change_gap_count1': after_change_gap_count1, 
        'number_gaps_filled_by_change1': number_gaps_filled_by_change1, 
        'percent_of_gaps_filled_by_change1': percent_of_gaps_filled_by_change1,
        'filled_voc': filled_voc,
        'var_init_data_full': var_init_data_full
    }

    cache_file_stage3 = f'{CACHE_DIR}/stage3/{var_name}_duplicate_filled.pkl'
    with open(cache_file_stage3, 'wb') as f:
        pickle.dump(stage3_dupfilling_results, f)
    return stage3_dupfilling_results 
def stage3_dupfilling_all_voc():
    stage3_dupfilling_results_all = []
    for row_idx in range(0, len(df_duplicates_and_tracers.index)):
        spec_name = df_duplicates_and_tracers['Varname'].iloc[row_idx]
        stage3_dupfilling_results = stage3_dupfilling_one_voc(var_species_name = spec_name)
        stage3_dupfilling_results_all.append(stage3_dupfilling_results)
    return stage3_dupfilling_results_all
def build_summary_stage3(stage3_dupfilling_results_all):
    summary_rows = []
    for r in stage3_dupfilling_results_all:
        summary_rows.append({    
            'Varname': r['varname'],
            'Duplicate_name': r['potential_duplicates'],
            'filled_with_name': r['filled_with_name'],
            'fillcase1': r['fill_case1'],
            'before_change_gap_count1': r['before_change_gap_count1'],
            'after_change_gap_count1': r['after_change_gap_count1'],
            'number_gaps_filled_by_change1': r['number_gaps_filled_by_change1'],
            'percent_of_gaps_filled_by_change1': r['percent_of_gaps_filled_by_change1'],
            'odr_eq_adj_used_for_correction': r['odr_eq_adj_used'],
            'rmse_used': r['rmse_used'],
            'r2_used': r['r2_used']
        })
    stage3_summary_savepath = stage_data_dirs['stage3'] + 'csv/stage3_summary_dup_fill.csv'
    pd.DataFrame(summary_rows).to_csv(stage3_summary_savepath, index=True)
    print('Saved summary to: ' + str(stage3_summary_savepath))
def duplicate_fill_timeseries(dupfilling_results, stage_type, flag_type):
    if flag_type == 'Duplicate list':
        var_name = dupfilling_results['varname']
        var_initial_data = dupfilling_results['var_init_data_full']
        #print('var_initial_data: ', var_initial_data)
        filled_name_0 = dupfilling_results['filled_with_name'][0]
        filled_data_0 = dupfilling_results['filled_voc'][0]
        #print(filled_name_0)
        filled_name_1 = dupfilling_results['filled_with_name'][1]
        filled_data_1 = dupfilling_results['filled_voc'][1]
        #print(filled_name_1)

        if pd.isna(filled_name_0) and pd.isna(filled_name_1):
            print('No dups filled')

        elif pd.isna(filled_name_0) and not pd.isna(filled_name_1):
            fig, (ax1, ax2) = plt.subplots(2,1, figsize = (16,8), tight_layout=True)
            xlim_start_jul = pd.to_datetime('2024-07-15 00:00:00').tz_localize('America/Denver')
            xlim_end_jul = pd.to_datetime('2024-07-31 23:45:00').tz_localize('America/Denver')
            xlim_start_aug = pd.to_datetime('2024-08-01 00:00:00').tz_localize('America/Denver')
            xlim_end_aug = pd.to_datetime('2024-08-18 23:45:00').tz_localize('America/Denver')

            #ax1 is the first row of subplot, for July only
            valid_points_initial = ~np.isnan(var_initial_data)
            ax1.plot(var_initial_data.index[valid_points_initial], var_initial_data[valid_points_initial], linestyle='solid', color = 'm', marker = 'x', label = f'Obs (gaps = {var_initial_data.isna().sum()})')
            # valid_points_filled_0 = ~np.isnan(filled_data_0)
            # ax1.plot(filled_data_0.index[valid_points_filled_0], filled_data_0[valid_points_filled_0], linestyle='solid', color = 'y', marker = '.', label = f'Filled with {filled_name_0} (gaps = {filled_name_0.isna().sum()})', alpha = 0.7)
            valid_points_filled_1 = ~np.isnan(filled_data_1)
            ax1.plot(filled_data_1.index[valid_points_filled_1], filled_data_1[valid_points_filled_1], linestyle='solid', color = 'y', marker = '.', label = f'Filled with {filled_name_1} (gaps = {filled_data_1.isna().sum()})', alpha = 0.7)
            
            #Set x ticks
            tz_mdt = df_all_measured_species.index.tz #this time zone should be in Mountain Daylight Time
            ax1.xaxis.set_major_locator(mdates.DayLocator(tz=tz_mdt))
            ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d', tz=tz_mdt))
            # Minor ticks: every 3 hours
            ax1.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21], tz=tz_mdt))
            # Rotate and format tick labels
            ax1.tick_params(axis='x', which='major')
            ax1.tick_params(axis='x', which='minor', length=3, color='gray')
            

            ax1.set_ylabel(var_name + ' (ppb)')
            #ax1.set_xlabel('Date')
            ax1.margins(x=0)
            ax1.set_xlim([xlim_start_jul, xlim_end_jul])

            ax1.legend(loc = 'upper right')

            #ax2 is the second row of subplot, for August only
            ax2.plot(var_initial_data.index[valid_points_initial], var_initial_data[valid_points_initial], linestyle='solid', color = 'm', marker = 'x', label = f'Obs (gaps = {var_initial_data.isna().sum()})')
            # ax2.plot(filled_data_0.index[valid_points_filled_0], filled_data_0[valid_points_filled_0], linestyle='solid', color = 'y', marker = '.', label = f'Filled with {filled_name_0} (gaps = {filled_data_0.isna().sum()})', alpha = 0.7)
            ax2.plot(filled_data_1.index[valid_points_filled_1], filled_data_1[valid_points_filled_1], linestyle='solid', color = 'y', marker = '.', label = f'Filled with {filled_name_1} (gaps = {filled_data_1.isna().sum()})', alpha = 0.7)
            
            #Set x ticks
            ax2.xaxis.set_major_locator(mdates.DayLocator(tz=tz_mdt))
            ax2.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d', tz=tz_mdt))
            # Minor ticks: every 3 hours
            ax2.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21], tz=tz_mdt))
            # Rotate and format tick labels
            ax2.tick_params(axis='x', which='major')
            ax2.tick_params(axis='x', which='minor', length=3, color='gray')

            ax2.set_ylabel(var_name + ' (ppb)')
            ax2.set_xlabel('Time (MDT)')
            ax2.margins(x=0)

            ax2.set_xlim([xlim_start_aug, xlim_end_aug])
            ax2.legend(loc = 'upper right')

            #Mark midnight for every day
            midnight_vals = []
            for midnight_idx in range(0,len(df_all_measured_species.index),96):
                midnight_vals.append(df_all_measured_species.index[midnight_idx])
            for day_pos in midnight_vals:
                ax1.axvline(day_pos, color = 'black', linestyle = 'dotted')
                ax2.axvline(day_pos, color = 'black', linestyle = 'dotted')

            plt.savefig(stage_data_dirs['stage3'] + 'plots/timeseries_full/filled_' + str(var_name) + '_from_dup_' + str(filled_name_0) + '_' + str(filled_name_1) + '_july_aug_timeseries.png', dpi =300)
            plt.show()
        
        elif not pd.isna(filled_name_0) and pd.isna(filled_name_1):
            fig, (ax1, ax2) = plt.subplots(2,1, figsize = (16,8), tight_layout=True)
            xlim_start_jul = pd.to_datetime('2024-07-15 00:00:00').tz_localize('America/Denver')
            xlim_end_jul = pd.to_datetime('2024-07-31 23:45:00').tz_localize('America/Denver')
            xlim_start_aug = pd.to_datetime('2024-08-01 00:00:00').tz_localize('America/Denver')
            xlim_end_aug = pd.to_datetime('2024-08-18 23:45:00').tz_localize('America/Denver')

            #ax1 is the first row of subplot, for July only
            valid_points_initial = ~np.isnan(var_initial_data)
            ax1.plot(var_initial_data.index[valid_points_initial], var_initial_data[valid_points_initial], linestyle='solid', color = 'm', marker = 'x', label = f'Obs (gaps = {var_initial_data.isna().sum()})')
            valid_points_filled_0 = ~np.isnan(filled_data_0)
            ax1.plot(filled_data_0.index[valid_points_filled_0], filled_data_0[valid_points_filled_0], linestyle='solid', color = 'y', marker = '.', label = f'Filled with {filled_name_0} (gaps = {filled_data_0.isna().sum()})', alpha = 0.7)
            
            #Set x ticks
            tz_mdt = df_all_measured_species.index.tz #this time zone should be in Mountain Daylight Time
            ax1.xaxis.set_major_locator(mdates.DayLocator(tz=tz_mdt))
            ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d', tz=tz_mdt))
            # Minor ticks: every 3 hours
            ax1.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21], tz=tz_mdt))
            # Rotate and format tick labels
            ax1.tick_params(axis='x', which='major')
            ax1.tick_params(axis='x', which='minor', length=3, color='gray')
            

            ax1.set_ylabel(var_name + ' (ppb)')
            #ax1.set_xlabel('Date')
            ax1.margins(x=0)
            ax1.set_xlim([xlim_start_jul, xlim_end_jul])

            ax1.legend(loc = 'upper right')

            #ax2 is the second row of subplot, for August only
            ax2.plot(var_initial_data.index[valid_points_initial], var_initial_data[valid_points_initial], linestyle='solid', color = 'm', marker = 'x', label = f'Obs (gaps = {var_initial_data.isna().sum()})')
            ax2.plot(filled_data_0.index[valid_points_filled_0], filled_data_0[valid_points_filled_0], linestyle='solid', color = 'y', marker = '.', label = f'Filled with {filled_name_0} (gaps = {filled_data_0.isna().sum()})', alpha = 0.7)
            
            #Set x ticks
            ax2.xaxis.set_major_locator(mdates.DayLocator(tz=tz_mdt))
            ax2.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d', tz=tz_mdt))
            # Minor ticks: every 3 hours
            ax2.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21], tz=tz_mdt))
            # Rotate and format tick labels
            ax2.tick_params(axis='x', which='major')
            ax2.tick_params(axis='x', which='minor', length=3, color='gray')

            ax2.set_ylabel(var_name + ' (ppb)')
            ax2.set_xlabel('Time (MDT)')
            ax2.margins(x=0)

            ax2.set_xlim([xlim_start_aug, xlim_end_aug])
            ax2.legend(loc = 'upper right')

            #Mark midnight for every day
            midnight_vals = []
            for midnight_idx in range(0,len(var_initial_data.index),96):
                midnight_vals.append(var_initial_data.index[midnight_idx])
            for day_pos in midnight_vals:
                ax1.axvline(day_pos, color = 'black', linestyle = 'dotted')
                ax2.axvline(day_pos, color = 'black', linestyle = 'dotted')

            plt.savefig(stage_data_dirs['stage3'] + 'plots/timeseries_full/filled_' + str(var_name) + '_from_dup_' + str(filled_name_0) + '_july_aug_timeseries.png', dpi =300)
            plt.show()

        else:
            fig, (ax1, ax2) = plt.subplots(2,1, figsize = (16,8), tight_layout=True)
            xlim_start_jul = pd.to_datetime('2024-07-15 00:00:00').tz_localize('America/Denver')
            xlim_end_jul = pd.to_datetime('2024-07-31 23:45:00').tz_localize('America/Denver')
            xlim_start_aug = pd.to_datetime('2024-08-01 00:00:00').tz_localize('America/Denver')
            xlim_end_aug = pd.to_datetime('2024-08-18 23:45:00').tz_localize('America/Denver')

            #ax1 is the first row of subplot, for July only
            valid_points_initial = ~np.isnan(var_initial_data)
            ax1.plot(var_initial_data.index[valid_points_initial], var_initial_data[valid_points_initial], linestyle='solid', color = 'm', marker = 'x', label = f'Obs (gaps = {var_initial_data.isna().sum()})')
            valid_points_filled_0 = ~np.isnan(filled_data_0)
            ax1.plot(filled_data_0.index[valid_points_filled_0], filled_data_0[valid_points_filled_0], linestyle='solid', color = 'y', marker = '.', label = f'Filled with {filled_name_0} (gaps = {filled_data_0.isna().sum()})', alpha = 0.7)
            valid_points_filled_1 = ~np.isnan(filled_data_1)
            ax1.plot(filled_data_1.index[valid_points_filled_1], filled_data_1[valid_points_filled_1], linestyle='solid', color = 'g', marker = '.', label = f'Filled with {filled_name_1} (gaps = {filled_data_1.isna().sum()})', alpha = 0.7)
            
            #Set x ticks
            tz_mdt = df_all_measured_species.index.tz #this time zone should be in Mountain Daylight Time
            ax1.xaxis.set_major_locator(mdates.DayLocator(tz=tz_mdt))
            ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d', tz=tz_mdt))
            # Minor ticks: every 3 hours
            ax1.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21], tz=tz_mdt))
            # Rotate and format tick labels
            ax1.tick_params(axis='x', which='major')
            ax1.tick_params(axis='x', which='minor', length=3, color='gray')
            

            ax1.set_ylabel(var_name + ' (ppb)')
            #ax1.set_xlabel('Date')
            ax1.margins(x=0)
            ax1.set_xlim([xlim_start_jul, xlim_end_jul])

            ax1.legend(loc = 'upper right')

            #ax2 is the second row of subplot, for August only
            ax2.plot(var_initial_data.index[valid_points_initial], var_initial_data[valid_points_initial], linestyle='solid', color = 'm', marker = 'x', label = f'Obs (gaps = {var_initial_data.isna().sum()})')
            ax2.plot(filled_data_0.index[valid_points_filled_0], filled_data_0[valid_points_filled_0], linestyle='solid', color = 'y', marker = '.', label = f'Filled with {filled_name_0} (gaps = {filled_data_0.isna().sum()})', alpha = 0.7)
            ax2.plot(filled_data_1.index[valid_points_filled_1], filled_data_1[valid_points_filled_1], linestyle='solid', color = 'g', marker = '.', label = f'Filled with {filled_name_1} (gaps = {filled_data_1.isna().sum()})', alpha = 0.7)
            
            #Set x ticks
            ax2.xaxis.set_major_locator(mdates.DayLocator(tz=tz_mdt))
            ax2.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d', tz=tz_mdt))
            # Minor ticks: every 3 hours
            ax2.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21], tz=tz_mdt))
            # Rotate and format tick labels
            ax2.tick_params(axis='x', which='major')
            ax2.tick_params(axis='x', which='minor', length=3, color='gray')

            ax2.set_ylabel(var_name + ' (ppb)')
            ax2.set_xlabel('Time (MDT)')
            ax2.margins(x=0)

            ax2.set_xlim([xlim_start_aug, xlim_end_aug])
            ax2.legend(loc = 'upper right')

            #Mark midnight for every day
            midnight_vals = []
            for midnight_idx in range(0,len(var_initial_data.index),96):
                midnight_vals.append(var_initial_data.index[midnight_idx])
            for day_pos in midnight_vals:
                ax1.axvline(day_pos, color = 'black', linestyle = 'dotted')
                ax2.axvline(day_pos, color = 'black', linestyle = 'dotted')

            plt.savefig(stage_data_dirs['stage3'] + 'plots/timeseries_full/filled_' + str(var_name) + '_from_dup_' + str(filled_name_1) + '_july_aug_timeseries.png', dpi =300)
            plt.show()

    else:        
        var_name = dupfilling_results['varname']
        var_initial_data = dupfilling_results['var_init_data_full']
        filled_name = dupfilling_results['filled_with_name']
        filled_data = dupfilling_results['filled_voc']
        
        fig, (ax1, ax2) = plt.subplots(2,1, figsize = (16,8), tight_layout=True)
        xlim_start_jul = pd.to_datetime('2024-07-15 00:00:00').tz_localize('America/Denver')
        xlim_end_jul = pd.to_datetime('2024-07-31 23:45:00').tz_localize('America/Denver')
        xlim_start_aug = pd.to_datetime('2024-08-01 00:00:00').tz_localize('America/Denver')
        xlim_end_aug = pd.to_datetime('2024-08-18 23:45:00').tz_localize('America/Denver')

        #ax1 is the first row of subplot, for July only
        valid_points_initial = ~np.isnan(var_initial_data)
        ax1.plot(var_initial_data.index[valid_points_initial], var_initial_data[valid_points_initial], linestyle='solid', color = 'm', marker = 'x', label = f'Obs (gaps = {var_initial_data.isna().sum()})')
        valid_points_filled = ~np.isnan(filled_data)
        ax1.plot(filled_data.index[valid_points_filled], filled_data[valid_points_filled], linestyle='solid', color = 'y', marker = '.', label = f'Filled with {filled_name} (gaps = {filled_data.isna().sum()})', alpha = 0.7)
        
        #Set x ticks
        tz_mdt = df_all_measured_species.index.tz #this time zone should be in Mountain Daylight Time
        ax1.xaxis.set_major_locator(mdates.DayLocator(tz=tz_mdt))
        ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d', tz=tz_mdt))
        # Minor ticks: every 3 hours
        ax1.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21], tz=tz_mdt))
        # Rotate and format tick labels
        ax1.tick_params(axis='x', which='major')
        ax1.tick_params(axis='x', which='minor', length=3, color='gray')
        

        ax1.set_ylabel(var_name + ' (ppb)')
        #ax1.set_xlabel('Date')
        ax1.margins(x=0)
        ax1.set_xlim([xlim_start_jul, xlim_end_jul])

        ax1.legend(loc = 'upper right')

        #ax2 is the second row of subplot, for August only
        ax2.plot(var_initial_data.index[valid_points_initial], var_initial_data[valid_points_initial], linestyle='solid', color = 'm', marker = 'x', label = f'Obs (gaps = {var_initial_data.isna().sum()})')
        ax2.plot(filled_data.index[valid_points_filled], filled_data[valid_points_filled], linestyle='solid', color = 'y', marker = '.', label = f'Filled with {filled_name} (gaps = {filled_data.isna().sum()})', alpha = 0.7)
        
        #Set x ticks
        ax2.xaxis.set_major_locator(mdates.DayLocator(tz=tz_mdt))
        ax2.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d', tz=tz_mdt))
        # Minor ticks: every 3 hours
        ax2.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21], tz=tz_mdt))
        # Rotate and format tick labels
        ax2.tick_params(axis='x', which='major')
        ax2.tick_params(axis='x', which='minor', length=3, color='gray')

        ax2.set_ylabel(var_name + ' (ppb)')
        ax2.set_xlabel('Time (MDT)')
        ax2.margins(x=0)

        ax2.set_xlim([xlim_start_aug, xlim_end_aug])
        ax2.legend(loc = 'upper right')

        #Mark midnight for every day
        midnight_vals = []
        for midnight_idx in range(0,len(var_initial_data.index),96):
            midnight_vals.append(var_initial_data.index[midnight_idx])
        for day_pos in midnight_vals:
            ax1.axvline(day_pos, color = 'black', linestyle = 'dotted')
            ax2.axvline(day_pos, color = 'black', linestyle = 'dotted')

        plt.savefig(stage_data_dirs['stage3'] + 'plots/timeseries_full/filled_' + str(var_name) + '_from_dup_' + str(filled_name) + '_july_aug_timeseries.png', dpi =300)
        plt.show()
def tracer_timeseries():

    fig, (ax1, ax2) = plt.subplots(2,1, figsize = (16,8), tight_layout=True)
    xlim_start_jul = pd.to_datetime('2024-07-15 00:00:00').tz_localize('America/Denver')
    xlim_end_jul = pd.to_datetime('2024-07-31 23:45:00').tz_localize('America/Denver')
    xlim_start_aug = pd.to_datetime('2024-08-01 00:00:00').tz_localize('America/Denver')
    xlim_end_aug = pd.to_datetime('2024-08-18 23:45:00').tz_localize('America/Denver')

    #ax1 is the first row of subplot, for July only
    valid_points_initial = ~np.isnan(filled_init)
    ax1.plot(filled_init.index[valid_points_initial], filled_init[valid_points_initial], linestyle='solid', color = 'm', marker = 'x', label = f'Obs (gaps = {filled_init.isna().sum()})')
    valid_points_filled = ~np.isnan(fill_first_tracer)
    ax1.plot(fill_first_tracer.index[valid_points_filled], fill_first_tracer[valid_points_filled], linestyle='solid', color = 'y', marker = '.', label = f'Filled (gaps = {fill_first_tracer.isna().sum()})', alpha = 0.7)

    #Set x ticks
    tz_mdt = df_all_measured_species.index.tz #this time zone should be in Mountain Daylight Time
    ax1.xaxis.set_major_locator(mdates.DayLocator(tz=tz_mdt))
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d', tz=tz_mdt))
    # Minor ticks: every 3 hours
    ax1.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21], tz=tz_mdt))
    # Rotate and format tick labels
    ax1.tick_params(axis='x', which='major')
    ax1.tick_params(axis='x', which='minor', length=3, color='gray')
    

    ax1.set_ylabel(target + ' (ppb)')
    #ax1.set_xlabel('Date')
    ax1.margins(x=0)
    ax1.set_xlim([xlim_start_jul, xlim_end_jul])

    ax1.legend(loc = 'upper right')

    #ax2 is the second row of subplot, for August only
    ax2.plot(filled_init.index[valid_points_initial], filled_init[valid_points_initial], linestyle='solid', color = 'm', marker = 'x', label = f'Obs (gaps = {filled_init.isna().sum()})')
    ax2.plot(fill_first_tracer.index[valid_points_filled], fill_first_tracer[valid_points_filled], linestyle='solid', color = 'y', marker = '.', label = f'Filled with {tracer_spec} (gaps = {fill_first_tracer.isna().sum()})', alpha = 0.7)
    #Set x ticks
    ax2.xaxis.set_major_locator(mdates.DayLocator(tz=tz_mdt))
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d', tz=tz_mdt))
    # Minor ticks: every 3 hours
    ax2.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21], tz=tz_mdt))
    # Rotate and format tick labels
    ax2.tick_params(axis='x', which='major')
    ax2.tick_params(axis='x', which='minor', length=3, color='gray')

    ax2.set_ylabel(target + ' (ppb)')
    ax2.set_xlabel('Time (MDT)')
    ax2.margins(x=0)

    ax2.set_xlim([xlim_start_aug, xlim_end_aug])
    ax2.legend(loc = 'upper right')

    #Mark midnight for every day
    midnight_vals = []
    # for midnight_idx in range(0,len(noaa.index),96):
    #     midnight_vals.append(noaa.index[midnight_idx])
    # for day_pos in midnight_vals:
    #     ax1.axvline(day_pos, color = 'black', linestyle = 'dotted')
    #     ax2.axvline(day_pos, color = 'black', linestyle = 'dotted')

    #plt.savefig(stage_data_dirs['stage2'] + 'plots/timeseries_full/filled_' + str(target) + '_day_from_tracer_' + tracer_spec + '_comparison_july_aug_timeseries.png', dpi =300)
    plt.show()

    # fig, (ax1, ax2) = plt.subplots(2,1, figsize = (16,8), tight_layout=True)
    # xlim_start_jul = pd.to_datetime('2024-07-15 00:00:00').tz_localize('America/Denver')
    # xlim_end_jul = pd.to_datetime('2024-07-31 23:45:00').tz_localize('America/Denver')
    # xlim_start_aug = pd.to_datetime('2024-08-01 00:00:00').tz_localize('America/Denver')
    # xlim_end_aug = pd.to_datetime('2024-08-18 23:45:00').tz_localize('America/Denver')

    # #ax1 is the first row of subplot, for July only
    # valid_points_initial = ~np.isnan(filled_init)
    # ax1.plot(filled_init.index[valid_points_initial], filled_init[valid_points_initial], linestyle='solid', color = 'm', marker = 'x', label = f'Obs (gaps = {filled_init.isna().sum()})')
    # valid_points_filled = ~np.isnan(fill_first_tracer)
    # ax1.plot(fill_first_tracer.index[valid_points_filled], fill_first_tracer[valid_points_filled], linestyle='solid', color = 'y', marker = '.', label = f'Filled (gaps = {fill_first_tracer.isna().sum()})', alpha = 0.7)

    # #Set x ticks
    # tz_mdt = df_all_measured_species.index.tz #this time zone should be in Mountain Daylight Time
    # ax1.xaxis.set_major_locator(mdates.DayLocator(tz=tz_mdt))
    # ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d', tz=tz_mdt))
    # # Minor ticks: every 3 hours
    # ax1.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21], tz=tz_mdt))
    # # Rotate and format tick labels
    # ax1.tick_params(axis='x', which='major')
    # ax1.tick_params(axis='x', which='minor', length=3, color='gray')
    

    # ax1.set_ylabel(target + ' (ppb)')
    # #ax1.set_xlabel('Date')
    # ax1.margins(x=0)
    # ax1.set_xlim([xlim_start_jul, xlim_end_jul])

    # ax1.legend(loc = 'upper right')

    # #ax2 is the second row of subplot, for August only
    # ax2.plot(filled_init.index[valid_points_initial], filled_init[valid_points_initial], linestyle='solid', color = 'm', marker = 'x', label = f'Obs (gaps = {filled_init.isna().sum()})')
    # ax2.plot(fill_first_tracer.index[valid_points_filled], fill_first_tracer[valid_points_filled], linestyle='solid', color = 'y', marker = '.', label = f'Filled with {tracer_spec} (gaps = {fill_first_tracer.isna().sum()})', alpha = 0.7)
    # #Set x ticks
    # ax2.xaxis.set_major_locator(mdates.DayLocator(tz=tz_mdt))
    # ax2.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d', tz=tz_mdt))
    # # Minor ticks: every 3 hours
    # ax2.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21], tz=tz_mdt))
    # # Rotate and format tick labels
    # ax2.tick_params(axis='x', which='major')
    # ax2.tick_params(axis='x', which='minor', length=3, color='gray')

    # ax2.set_ylabel(target + ' (ppb)')
    # ax2.set_xlabel('Time (MDT)')
    # ax2.margins(x=0)

    # ax2.set_xlim([xlim_start_aug, xlim_end_aug])
    # ax2.legend(loc = 'upper right')

    # #Mark midnight for every day
    # midnight_vals = []
    # # for midnight_idx in range(0,len(noaa.index),96):
    # #     midnight_vals.append(noaa.index[midnight_idx])
    # # for day_pos in midnight_vals:
    # #     ax1.axvline(day_pos, color = 'black', linestyle = 'dotted')
    # #     ax2.axvline(day_pos, color = 'black', linestyle = 'dotted')

    # plt.savefig(stage_data_dirs['stage2'] + 'plots/timeseries_full/filled_' + str(target) + '_night_from_tracer_' + tracer_spec + '_comparison_july_aug_timeseries.png', dpi =300)
    # plt.show()
def stage3_tracer_filling_one_voc(var_species_name):
    target = var_species_name

    #get filled species from stage3 dupfilling
    cache_file_stage3 = f'{CACHE_DIR}/stage3/{target}_duplicate_filled.pkl'
    with open(cache_file_stage3, 'rb') as f:
        target_filled_load = pickle.load(f)

    if isinstance(target_filled_load['filled_voc'], list):
        target_filled_data = target_filled_load['filled_voc'][1]
    else:
        target_filled_data = target_filled_load['filled_voc']

    #get tracer data
    cache_file_name_top_10_day = f'{CACHE_DIR}/stage2/{target}_top10_tracers_day_stage2.pkl'
    df_top_species_day = pd.read_pickle(cache_file_name_top_10_day)
    #print('pickle load df_top_species_per_regime: ', df_top_species_per_regime)
    print(df_top_species_day.columns)
    top_tracers_day = df_top_species_day['tracer']
    print(f'Top tracers for {target} during the day are: ', top_tracers_day)

    filled_day_new = []
    before_change_gap_count1_day_new = []
    after_change_gap_count1_day_new = []
    number_gaps_filled_by_change1_day_new = []
    percent_of_gaps_filled_by_change1_day_new = []
    y_data_input = None
    y_data_init_list = []

    if df_top_species_day.empty:
        print(f'No tracer species for day fill for {target}')
        y_data_input = target_filled_data
        (before_change_gap_count1, after_change_gap_count1, 
            number_gaps_filled_by_change1, percent_of_gaps_filled_by_change1) = gap_counter(data_before_change = y_data_input, 
                                                                                        data_after_change = target_filled_data)
        stage3_tracerfilling_results_day = {
            'varname': target,
            'tracers_filled_name': np.nan,
            'odr_eq_adj_used': np.nan,
            'norm_rmse': np.nan,
            'r2_tracer': np.nan,
            'before_change_gap_count1': before_change_gap_count1, 
            'after_change_gap_count1': after_change_gap_count1,
            'number_gaps_filled_by_change1': number_gaps_filled_by_change1, 
            'percent_of_gaps_filled_by_change1': percent_of_gaps_filled_by_change1,
            'filled_voc': target_filled_data,
            'var_init_data_full': target_filled_data
        }

    else:
        for tracer_number in range(0, len(df_top_species_day)):
            filled_init = target_filled_data.copy()
            slope_tracer = df_top_species_day.iloc[tracer_number]['slope']
            intercept_tracer = df_top_species_day.iloc[tracer_number]['intercept']
            norm_rmse_tracer = df_top_species_day.iloc[tracer_number]['norm_rmse']
            rsquared_tracer = df_top_species_day.iloc[tracer_number]['r2']

            tracer_spec = df_top_species_day.iloc[tracer_number]['tracer']
            print('Filling ' + target +  'with ' + tracer_spec + ' tracer')

            xvals_corrected_overlap_tracer_species = df_top_species_day.iloc[tracer_number]['xvals_corrected_overlap_no_negs']

            full_index = pd.date_range(start="2024-07-15 00:00:00", end="2024-08-18 17:45:00", 
                                       freq="15min", 
                                       tz = df_all_measured_species.index.tz)
            #print('df_top_species_day.iloc[tracer_number]: ', df_top_species_day.iloc[tracer_number])

            xvals_corrected_tracer_species = df_top_species_day.iloc[tracer_number]['xvals_corrected_all_points_no_negs']
            xvals_corrected_tracer_species.reindex(full_index)

            day_mask = (hours >= 6) & (hours < 18)

            if y_data_input is None:
                y_data_input = filled_init   #first pass uses original y
                #print('y_data_input2: ', y_data_input)
            else:
                y_data_input = fill_first_tracer  #carry forward
                #print('y_data_input2: ', y_data_input)
            #print('xvals_corrected_tracer_species: ', xvals_corrected_tracer_species)
            #print('xvals_corrected_tracer_species[day_mask]: ', xvals_corrected_tracer_species)
            fill_first_tracer = y_data_input.fillna(xvals_corrected_tracer_species)
            
            (before_change_gap_count1, after_change_gap_count1, 
            number_gaps_filled_by_change1, percent_of_gaps_filled_by_change1) = gap_counter(data_before_change = y_data_input, 
                                                                                        data_after_change = fill_first_tracer)
            filled_day_new.append(fill_first_tracer)
            before_change_gap_count1_day_new.append(before_change_gap_count1)
            after_change_gap_count1_day_new.append(after_change_gap_count1)
            number_gaps_filled_by_change1_day_new.append(number_gaps_filled_by_change1)
            percent_of_gaps_filled_by_change1_day_new.append(percent_of_gaps_filled_by_change1)
            y_data_init_list.append(y_data_input)


        stage3_tracerfilling_results_day = {
            'varname': target,
            'tracers_filled_name': top_tracers_day,
            'odr_eq_adj_used': 'y = '+ str(intercept_tracer) + ' + ' + str(slope_tracer) + '* M',
            'norm_rmse': norm_rmse_tracer,
            'r2_tracer': rsquared_tracer,
            'before_change_gap_count1': before_change_gap_count1_day_new, 
            'after_change_gap_count1': after_change_gap_count1_day_new, 
            'number_gaps_filled_by_change1': number_gaps_filled_by_change1_day_new, 
            'percent_of_gaps_filled_by_change1': percent_of_gaps_filled_by_change1_day_new,
            'filled_voc': filled_day_new,
            'var_init_data_full': y_data_init_list
        }
    print("stage3_tracerfilling_results_day['filled_voc']: ", stage3_tracerfilling_results_day['filled_voc'])

    #get tracer data
    cache_file_name_top_10_night = f'{CACHE_DIR}/stage2/{target}_top10_tracers_night_stage2.pkl'
    df_top_species_night = pd.read_pickle(cache_file_name_top_10_night)
    #print('pickle load df_top_species_per_regime: ', df_top_species_per_regime)
    print(df_top_species_night.columns)
    top_tracers_night = df_top_species_night['tracer']
    print(f'Top tracers for {target} during the night are: ', top_tracers_night)

    filled_night_new = []
    before_change_gap_count1_night_new = []
    after_change_gap_count1_night_new = []
    number_gaps_filled_by_change1_night_new = []
    percent_of_gaps_filled_by_change1_night_new = []
    y_data_input = None
    y_data_init_list = []

    if df_top_species_night.empty:
        print(f'No tracer species for night fill for {target}')
        if isinstance(stage3_tracerfilling_results_day['filled_voc'], list):
            y_data_input = stage3_tracerfilling_results_day['filled_voc'][-1]
            print('y_data_input for empty night fill: ', y_data_input)
            (before_change_gap_count1, after_change_gap_count1, 
                number_gaps_filled_by_change1, percent_of_gaps_filled_by_change1) = gap_counter(data_before_change = y_data_input, 
                                                                                            data_after_change = y_data_input)
            print('before_change_gap_count1 empty night1: ', before_change_gap_count1)
        elif isinstance(stage3_tracerfilling_results_day['filled_voc'], pd.Series):
            y_data_input = stage3_tracerfilling_results_day['filled_voc']
            print('y_data_input for empty night fill: ', y_data_input)
            (before_change_gap_count1, after_change_gap_count1, 
                number_gaps_filled_by_change1, percent_of_gaps_filled_by_change1) = gap_counter(data_before_change = y_data_input, 
                                                                                            data_after_change = y_data_input)
            print('before_change_gap_count1 empty night2: ', before_change_gap_count1)
        stage3_tracerfilling_results_night = {
            'varname': target,
            'tracers_filled_name': np.nan,
            'odr_eq_adj_used': np.nan,
            'norm_rmse': np.nan,
            'r2_tracer': np.nan,
            'before_change_gap_count1': before_change_gap_count1, 
            'after_change_gap_count1': after_change_gap_count1,
            'number_gaps_filled_by_change1': number_gaps_filled_by_change1, 
            'percent_of_gaps_filled_by_change1': percent_of_gaps_filled_by_change1,
            'filled_voc': y_data_input,
            'var_init_data_full': y_data_init_list.append(y_data_input)
        }


    else:
        for tracer_number in range(0, len(df_top_species_night)):
            # filled_init = target_filled_data.copy()
            slope_tracer = df_top_species_night.iloc[tracer_number]['slope']
            intercept_tracer = df_top_species_night.iloc[tracer_number]['intercept']
            norm_rmse_tracer = df_top_species_night.iloc[tracer_number]['norm_rmse']
            rsquared_tracer = df_top_species_night.iloc[tracer_number]['r2']

            tracer_spec = df_top_species_night.iloc[tracer_number]['tracer']
            print('Filling ' + target +  'with ' + tracer_spec + ' tracer')

            xvals_corrected_overlap_tracer_species = df_top_species_night.iloc[tracer_number]['xvals_corrected_overlap_no_negs']

            full_index = pd.date_range(start="2024-07-15 00:00:00", end="2024-08-18 17:45:00", 
                                       freq="15min", 
                                       tz = df_all_measured_species.index.tz)
            #print('df_top_species_night.iloc[tracer_number]: ', df_top_species_night.iloc[tracer_number])

            xvals_corrected_tracer_species = df_top_species_night.iloc[tracer_number]['xvals_corrected_all_points_no_negs']
            xvals_corrected_tracer_species.reindex(full_index)
            print('xvals_corrected_tracer_species: ', xvals_corrected_tracer_species)

            night_mask = (hours < 6) & (hours >= 18)

            if isinstance(stage3_tracerfilling_results_day['filled_voc'], list):
                if y_data_input is None:
                    y_data_input = stage3_tracerfilling_results_day['filled_voc'][-1]   #first pass uses original y
                    print('y_data_input1: ', y_data_input)

                else:
                    y_data_input = fill_first_tracer  #carry forward
                    print('y_data_input2: ', y_data_input)

                fill_first_tracer = y_data_input.fillna(xvals_corrected_tracer_species)
            
                (before_change_gap_count1, after_change_gap_count1, 
                number_gaps_filled_by_change1, percent_of_gaps_filled_by_change1) = gap_counter(data_before_change = y_data_input, 
                                                                                            data_after_change = fill_first_tracer)
            elif isinstance(stage3_tracerfilling_results_day['filled_voc'], pd.Series):
                if y_data_input is None:
                    y_data_input = stage3_tracerfilling_results_day['filled_voc']   #first pass uses original y
                    print('y_data_input1: ', y_data_input)
                else:
                    y_data_input = fill_first_tracer  #carry forward
                    print('y_data_input2: ', y_data_input)
                fill_first_tracer = y_data_input.fillna(xvals_corrected_tracer_species)
                (before_change_gap_count1, after_change_gap_count1, 
                number_gaps_filled_by_change1, percent_of_gaps_filled_by_change1) = gap_counter(data_before_change = y_data_input, 
                                                                                            data_after_change = fill_first_tracer)

            #print('xvals_corrected_tracer_species: ', xvals_corrected_tracer_species)
            #print('xvals_corrected_tracer_species[night_mask]: ', xvals_corrected_tracer_species)
            print('xvals_corrected_tracer_species: ', xvals_corrected_tracer_species)

            filled_night_new.append(fill_first_tracer)
            before_change_gap_count1_night_new.append(before_change_gap_count1)
            after_change_gap_count1_night_new.append(after_change_gap_count1)
            number_gaps_filled_by_change1_night_new.append(number_gaps_filled_by_change1)
            percent_of_gaps_filled_by_change1_night_new.append(percent_of_gaps_filled_by_change1)
            y_data_init_list.append(y_data_input)
            

        stage3_tracerfilling_results_night = {
            'varname': target,
            'tracers_filled_name': top_tracers_night,
            'odr_eq_adj_used': 'y = '+ str(intercept_tracer) + ' + ' + str(slope_tracer) + '* M',
            'norm_rmse': norm_rmse_tracer,
            'r2_tracer': rsquared_tracer,
            'before_change_gap_count1': before_change_gap_count1_night_new, 
            'after_change_gap_count1': after_change_gap_count1_night_new, 
            'number_gaps_filled_by_change1': number_gaps_filled_by_change1_night_new, 
            'percent_of_gaps_filled_by_change1': percent_of_gaps_filled_by_change1_night_new,
            'filled_voc': filled_night_new,
            'var_init_data_full': y_data_init_list
        }
    cache_file_stage3_tracer_fill_day = f'{CACHE_DIR}/stage3/{target}_day_tracer_filled.pkl'
    with open(cache_file_stage3_tracer_fill_day, 'wb') as f:
        pickle.dump(stage3_tracerfilling_results_day, f)
    cache_file_stage3_tracer_fill_night = f'{CACHE_DIR}/stage3/{target}_night_tracer_filled.pkl'
    with open(cache_file_stage3_tracer_fill_night, 'wb') as f:
        pickle.dump(stage3_tracerfilling_results_night, f)

    return stage3_tracerfilling_results_day, stage3_tracerfilling_results_night
def stage3_tracer_filling_all_voc():
    # Load pickle files with target and tracers
    print('Loading tracers for sequential filling')

    stage3_tracerfilling_results_day_list = []
    stage3_tracerfilling_results_night_list = []
    for row_idx in range(0, len(df_duplicates_and_tracers.index)):
        spec_name = df_duplicates_and_tracers['Varname'].iloc[row_idx]
        stage3_tracerfilling_results_day, stage3_tracerfilling_results_night = stage3_tracer_filling_one_voc(var_species_name = spec_name)
        stage3_tracerfilling_results_day_list.append(stage3_tracerfilling_results_day)
        stage3_tracerfilling_results_night_list.append(stage3_tracerfilling_results_night)
    return stage3_tracerfilling_results_day_list, stage3_tracerfilling_results_night_list
def build_summary_stage3_tracer(results, regime_type):
    summary_rows = []
    if regime_type == 'day':
        for r in stage3_tracerfilling_results_day_list:
            summary_rows.append({    
                'Varname': r['varname'],
                'tracers_filled_name': r['tracers_filled_name'],
                'before_change_gap_count1': r['before_change_gap_count1'],
                'after_change_gap_count1': r['after_change_gap_count1'],
                'number_gaps_filled_by_change1': r['number_gaps_filled_by_change1'],
                'percent_of_gaps_filled_by_change1': r['percent_of_gaps_filled_by_change1'],
                'odr_eq_adj_used_for_correction': r['odr_eq_adj_used'],
                'norm_rmse_used': r['norm_rmse'],
                'r2_used': r['r2_tracer']
            })
        stage3_summary_savepath = stage_data_dirs['stage3'] + 'csv/stage3_summary_tracer_fill_day.csv'
        pd.DataFrame(summary_rows).to_csv(stage3_summary_savepath, index=True)
        print('Saved summary to: ' + str(stage3_summary_savepath))

    elif regime_type == 'night':
        for r in stage3_tracerfilling_results_night_list:
            summary_rows.append({    
                'Varname': r['varname'],
                'tracers_filled_name': r['tracers_filled_name'],
                'before_change_gap_count1': r['before_change_gap_count1'],
                'after_change_gap_count1': r['after_change_gap_count1'],
                'number_gaps_filled_by_change1': r['number_gaps_filled_by_change1'],
                'percent_of_gaps_filled_by_change1': r['percent_of_gaps_filled_by_change1'],
                'odr_eq_adj_used_for_correction': r['odr_eq_adj_used'],
                'norm_rmse_used': r['norm_rmse'],
                'r2_used': r['r2_tracer']
            })
        stage3_summary_savepath = stage_data_dirs['stage3'] + 'csv/stage3_summary_tracer_fill_night.csv'
        pd.DataFrame(summary_rows).to_csv(stage3_summary_savepath, index=True)
        print('Saved summary to: ' + str(stage3_summary_savepath))

def interpolation_fill(target):
    cache_file_stage3_tracer_fill_night = f'{CACHE_DIR}/stage3/{target}_night_tracer_filled.pkl'
    with open(cache_file_stage3_tracer_fill_night, 'rb') as f:
        r = pickle.load(f)

    if isinstance(r['filled_voc'], list):
        filled_data = r['filled_voc'][-1]
    else:
        filled_data = r['filled_voc']

    print(filled_data)


    # #convert datetime index to UTC

    # df_iwas_updated = df_iwas_updated.rename(columns=mapping)
    # #rename index to 'time_UTC' to match netcdf file index name
    # df_iwas_updated.index.rename('time_UTC', inplace=True)
    # #We need to use timezone naive datetime index for the netCDF file
    # df_iwas_updated.index = df_iwas_updated.index.tz_convert('UTC').tz_localize(None)
    
    # #Now we are going to replace the old netCDF data with the new data
    # ds_nc = xr.open_dataset('/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/CampaignData_and_Merges/R0/CSL_MobileLab_Parked/merged/rev_15min/all_CSL_MobileLab_Parked_rev15min.nc')
    # print(ds_nc)
    # #Convert the dataframe of updated iWAS data to xarray Dataset
    # new_iwas_data = xr.Dataset.from_dataframe(df_iwas_updated)
    # #Align time coordinates to NetCDF time
    # new_iwas_data = new_iwas_data.reindex(time_UTC = ds_nc.time_UTC)

    # for var in new_iwas_data.data_vars:
    #     # Replace old values with new ones for each relevant variable
    #     ds_nc[var] = new_iwas_data[var]
    # save_ncfilename = savefilename + '.nc'
    # full_savepath = merged_data_dir_15min + save_ncfilename
    # ds_nc.to_netcdf(full_savepath, format = 'NETCDF4',mode='w')
    # print('Saved netCDF file with updated iWAS measurements to: ', full_savepath)

    # #Check that the format looks consistent with the original ds_nc
    # ds_new = xr.open_dataset(full_savepath)
    # print(ds_new)
def all_species_interp(df_filled_data):
    df_interp=df_filled_data.copy()
    #print(df_interp.columns.tolist())

    # interpolation won't work with the time_UTC column, need to get a new dataframe that doesn't include it to loop over
    df_filled_data.loc[df_filled_data['Lon'] <= -111.87, 'Lon'] = -111.87211

    #Manually edit some issues with longitude
    df_interp.loc[:, 'Lon'] = df_interp['Lon'].interpolate(method='linear')
    df_interp.loc[:, 'Lon'] = df_interp['Lon'].fillna(-111.87211)

    # #Convert any negative except in Longitude into a NaN
    cols_lon_exclude = df_interp.columns.difference(['Lon', 'time_UTC', 'time_local'])
    df_interp.loc[:, cols_lon_exclude] = df_interp.loc[:, cols_lon_exclude].mask(df_interp.loc[:, cols_lon_exclude] < 0, np.nan)

    # Apply interpolation to every variable except Longitude
    df_interp.loc[:, cols_lon_exclude] = df_interp.loc[:, cols_lon_exclude].interpolate(method='linear')

    return df_filled_data, df_interp 
def plot_interps(df_filled_data, df_interp, col):
    """
    This function shows the interpolation plots for our 15 minute revised merges. Currently, any time there is data missing from the USOS Campaign, 
    there is a NaN as the value. These plots will show what it looks like to fill a NaN with an interpolated value.
    """
    #Plot interpolated columns so we can see how well the interpolation performed

    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4,1, figsize = (30,20), tight_layout=True)

    #Each subplot shows approximately 1 week
    xlim_start_w1 = pd.to_datetime('2024-07-15 00:00:00')
    xlim_end_w1 = pd.to_datetime('2024-07-23 23:00:00')
    xlim_start_w2 = pd.to_datetime('2024-07-24 00:00:00')
    xlim_end_w2 = pd.to_datetime('2024-07-31 23:00:00')
    xlim_start_w3 = pd.to_datetime('2024-08-01 00:00:00')
    xlim_end_w3 = pd.to_datetime('2024-08-08 23:00:00')
    xlim_start_w4 = pd.to_datetime('2024-08-09 00:00:00')
    xlim_end_w4 = pd.to_datetime('2024-08-18 23:00:00')

    n_baddies= len([item for item in df_filled_data[col] if item <0 or np.isnan(item)]) 
    ax1.plot(df_filled_data.index, df_filled_data[col], color='k', marker='o',label=f'Original (baddies={n_baddies})')
    ax1.plot(df_interp.index, df_interp[col], color='r', marker='x', label='Interpolated')

    #Set x ticks
    tz_mdt = df_all_measured_species.index.tz
    ax1.xaxis.set_major_locator(mdates.DayLocator(tz=tz_mdt))
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d', tz=tz_mdt))
    # Minor ticks: every 3 hours
    ax1.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21], tz=tz_mdt))
    # Rotate and format tick labels
    ax1.tick_params(axis='x', which='major')
    ax1.tick_params(axis='x', which='minor', length=3, color='gray')
    ax1.set_xlim([xlim_start_w1, xlim_end_w1])
    ax1.legend(loc = 'upper right')

    ax2.plot(df_filled_data.index, df_filled_data[col], color='k', marker='o',label=f'Original (baddies={n_baddies})')
    ax2.plot(df_interp.index, df_interp[col], color='r', marker='x', label='Interpolated')

    # #Set x ticks
    ax2.xaxis.set_major_locator(mdates.DayLocator(tz=tz_mdt))
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d', tz=tz_mdt))
    # Minor ticks: every 3 hours
    ax2.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21], tz=tz_mdt))
    # Rotate and format tick labels
    ax2.tick_params(axis='x', which='major')
    ax2.tick_params(axis='x', which='minor', length=3, color='gray')
    ax2.set_xlim([xlim_start_w2, xlim_end_w2])
    ax2.legend(loc = 'upper right')

    ax3.plot(df_filled_data.index, df_filled_data[col], color='k', marker='o',label=f'Original (baddies={n_baddies})')
    ax3.plot(df_interp.index, df_interp[col], color='r', marker='x', label='Interpolated')

    # #Set x ticks
    ax3.xaxis.set_major_locator(mdates.DayLocator(tz=tz_mdt))
    ax3.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d', tz=tz_mdt))
    # Minor ticks: every 3 hours
    ax3.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21], tz=tz_mdt))
    # Rotate and format tick labels
    ax3.tick_params(axis='x', which='major')
    ax3.tick_params(axis='x', which='minor', length=3, color='gray')
    ax3.set_xlim([xlim_start_w3, xlim_end_w3])
    ax3.legend(loc = 'upper right')

    ax4.plot(df_filled_data.index, df_filled_data[col], color='k', marker='o',label=f'Original (baddies={n_baddies})')
    ax4.plot(df_interp.index, df_interp[col], color='r', marker='x', label='Interpolated')

    # #Set x ticks
    ax4.xaxis.set_major_locator(mdates.DayLocator(tz=tz_mdt))
    ax4.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d', tz=tz_mdt))
    # Minor ticks: every 3 hours
    ax4.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21], tz=tz_mdt))
    # Rotate and format tick labels
    ax4.tick_params(axis='x', which='major')
    ax4.tick_params(axis='x', which='minor', length=3, color='gray')
    ax4.set_xlim([xlim_start_w4, xlim_end_w4])
    ax4.legend(loc = 'upper right')

    plt.suptitle(col)
    plt.savefig(stage_data_dirs['stage4'] + 'plots/timeseries_full/interpolation_' + str(col) + '_july_aug_timeseries.png', dpi =150)
    plt.show()
def dataframe_to_nested_dict(df):
    nested_dict = {}
    for column in df.columns:
        nested_dict[column] = df[column].to_numpy()
        
    # Add the index as a key-value pair
    #nested_dict['time_local'] = df[column].to_numpy()
    return nested_dict
def subset_day_interp_save_to_matlab(df_interp, var_name):

    # Put UTC column back into the dataframe
    df_interp['time_UTC'] = df_filled_data['time_UTC'].tz_localize('UTC')

    for date, df_day in df_interp.resample('D'):
        # if df_day.empty:
        #     continue  # skip days with no data
    
        start_str = date.strftime('%Y%m%d')
        end_str = (date + pd.Timedelta(days=1)).strftime('%Y%m%d')

        measurement_filename_str = f'{start_str}_{end_str}_15min_CSL_mobile_lab_parked_shortened'
        measurement_final_filename = f'/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/CampaignData_and_Merges/R0/CSL_MobileLab_Parked/F0AM_filled/{measurement_filename_str}.csv'
        
        df_day.to_csv(measurement_final_filename)
        print('Saved CSV to:' + measurement_final_filename)

        # Convert the dataframe to a nested dictionary (so scipy can output to a matlab structure!) 
        ddict=dataframe_to_nested_dict(df_day)

        # Sort alphabetically so not annoying in MATLAB...  
        ddict= OrderedDict(sorted(ddict.items())) 

        # Save the USOS data in an output .mat file: 
        outpath = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/F0AM-4.4.2/Campaign_Data/matlab_merge/parked/original/'
        matfilename = measurement_filename_str + '.mat'
        savemat(outpath+matfilename,{var_name: ddict})
        print('Saved MATLAB file to:' + outpath + matfilename)
def all_days_interp_save_to_matlab(df_interp, var_name):
    # Put UTC column back into the dataframe
    #df_interp['time_UTC'] = df_filled_data['time_UTC'].tz_localize('UTC')
    # df_interp['time_local'] = df_interp['time_local'].tz_localize('America/Denver')
    print(df_interp['time_local'])
    print(df_interp['time_UTC'])

    measurement_filename_str = f'filled_alldays_15min_CSL_mobile_lab_parked'
    measurement_final_filename = f'/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/CampaignData_and_Merges/R0/CSL_MobileLab_Parked/F0AM_filled/{measurement_filename_str}.csv'
        
    df_interp.to_csv(measurement_final_filename)
    print('Saved CSV to:' + measurement_final_filename)

    #MATLAB equates the time zones to the same UTC so we can get rid of the time_local column
    df_interp['time_UTC_unix'] = (df_interp['time_UTC'].astype('int64') / 1e9)
    df_interp = df_interp.drop(columns = ['time_local', 'time_UTC'])

    # Convert the dataframe to a nested dictionary (so scipy can output to a matlab structure!) 
    ddict=dataframe_to_nested_dict(df_interp)
    print(ddict)

    # Sort alphabetically so not annoying in MATLAB...  
    ddict= OrderedDict(sorted(ddict.items())) 
    print(ddict)
    
    # Save the USOS data in an output .mat file: 
    outpath = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/F0AM-4.4.2/Campaign_Data/matlab_merge/parked/original/'
    matfilename = measurement_filename_str + '.mat'
    savemat(outpath+matfilename,{var_name: ddict})
    print('Saved MATLAB file to:' + outpath + matfilename)

if __name__ == "__main__":
    # merge_all_data_sources()
    all_measured_species_f = gap_filling_csv_data_dir + 'all_measured_species.csv'
    df_all_measured_species = pd.read_csv(all_measured_species_f, index_col='time_local', parse_dates=True)
    df_all_measured_species.index = df_all_measured_species.index.tz_localize(None)
    df_all_measured_species.index = df_all_measured_species.index.tz_localize('America/Denver')

    #can't find code i used to make this csv file???
    duplicates_and_tracers_f = gap_filling_csv_data_dir + 'species_duplicates_tracers.csv'
    df_duplicates_and_tracers = pd.read_csv(duplicates_and_tracers_f, index_col='Index')
    #print(' df_duplicates_and_tracers: ', df_duplicates_and_tracers)
# #####################################################################
    # print('Starting Stage 1')

    # stage1_results_all = process_all_vocs_stage1()
    # build_summary_stage1(stage1_results_all)
    
    #stage1_results_all is a list of dictionaries, storing keys: 'x_species_name', 
    # 'y_species_name', 'case', 'points_considered_in_odr', 'odr_eq_adj', 
    # 'metrics_init_slope', 'metrics_init_intercept', 'metrics_init_rmse', 
    # 'metrics_init_norm_rmse', 'metrics_init_r2', 'metrics_corr_slope', 
    # 'metrics_corr_intercept', 'metrics_corr_rmse', 'metrics_corr_norm_rmse',
    #  'metrics_corr_r2', 'did_slope_improve' , 'slope_distance_from_1_init', 
    # 'slope_distance_from_1_corr', 'slope_distance_improvement_val', 'intercept_error_init', 
    # 'intercept_error_corr', 'did_rmse_improve', 'rmse_percent_improvement', 
    # 'var_data_init_full', 'var_data_init_overlap', 'dup_data_init_full',
    #  'dup_data_init_overlap', 'dup_data_corr_full', 'dup_data_corr_overlap'

    # # # Comment out if unnecessary to save or see plots
    # for item in stage1_results_all:
    #     if item['case'] == 'All NaNs in duplicate':
    #         print('No duplicate data for ' + str(item['y_species_name']) + ' in Stage 1.')
    #     elif item['case'] == 'All NaNs in var':
    #         print('No var data for ' + str(item['y_species_name']) + ' in Stage 1.')
    #     elif item['case'] == 'Has duplicate list, has duplicate data not all NaNs, has vardata not all NaNs':
    #         scatterplots_for_comparing_init_and_corr_odr_fits(odr_fit_results = item, stage_type = 'stage1', flag_type = 'Duplicate list')
    #     else:
    #         scatterplots_for_comparing_init_and_corr_odr_fits(odr_fit_results = item, stage_type = 'stage1', flag_type = None)

# #####################################################################
    # STAGE 2: Determine Tracers
    # hours = df_all_measured_species.index.hour
    # regime_mask = {"day": (hours >= 6) & (hours < 18),
    #                 "night": (hours < 6) | (hours >= 18),
    #                 "total": np.ones(len(df_all_measured_species), dtype=bool)}
    
    # df_duplicates_and_tracers['Expected_tracers'] = df_duplicates_and_tracers['Expected_tracers'].apply(lambda x: x.split('; '))
    # species_tracer_dict = {}
    # for idx in range(0, len(df_duplicates_and_tracers.index)):
    #     species_tracer_dict.update({df_duplicates_and_tracers['Varname'][idx] : df_duplicates_and_tracers['Expected_tracers'][idx]})
    #don't need to run stage2_odr_fitting_for_tracer every time if it's updated, stage2_tracer_scoring will use cache
    # stage2_odr_fitting_for_tracer()

    # set function parameters for scoring tracers:
    # currently set to have a maximum normalized RMSE of 0.6, meet a minimum of 0.05 for r squared value
    # stage2_tracer_scoring(corr_threshold_tracers = 0.75, norm_rmse_max = 0.6, r2_min = 0.05)

#####################################################################
    # stage2_savedata_to_csv()

    # # Plot scatterplots for each tracer with its predicted value for each target species
    # regime_list = ['day', 'night']
    
    # for target in species_tracer_dict:
    #     for reg in regime_list:
    #         cache_file_name_top_10_per_regime = f'{CACHE_DIR}/stage2/{target}_top10_tracers_{reg}_stage2.pkl'
    #         df_top_species_per_regime = pd.read_pickle(cache_file_name_top_10_per_regime)
    #         for row_idx in range(0, len(df_top_species_per_regime.index)):
    #             # target = df_duplicates_and_tracers['Varname'].iloc[row_idx]
    #             ypred_plot = df_top_species_per_regime['xvals_corrected_overlap'].values[row_idx]
    #             ypred_no_negs = df_top_species_per_regime['xvals_corrected_overlap_no_negs'].values[row_idx]
    #             xvars_plot = df_top_species_per_regime['xvar_overlap'].values[row_idx]
    #             yvars_plot = df_top_species_per_regime['yvar_overlap'].values[row_idx]
    #             yvar_no_negs = df_top_species_per_regime['y_overlap_valid_no_negs'].values[row_idx]
    #             tracer_name = df_top_species_per_regime['tracer'].values[row_idx]
    #             slope_plot = df_top_species_per_regime['slope'].values[row_idx]
    #             intercept_plot = df_top_species_per_regime['intercept'].values[row_idx]
    #             r2_plot =  df_top_species_per_regime['r2'].values[row_idx]
    #             norm_rmse_plot = df_top_species_per_regime['norm_rmse'].values[row_idx]
    #             gap_fill_potential_plot = df_top_species_per_regime['gap_fill_potential'].values[row_idx]
    #             gap_needed_to_fill_plot = df_top_species_per_regime['gaps_needed_to_fill_for_regime'].values[row_idx]
    #             tracer_scatterplots(regime = reg, 
    #                                 xdata = xvars_plot, 
    #                                 ydata = yvars_plot,
    #                                 ypred = ypred_plot,
    #                                 ypred_no_negs = ypred_no_negs,
    #                                 yvar_no_negs = yvar_no_negs,
    #                                 xname = tracer_name,
    #                                 yname = target, slope_plot = slope_plot, 
    #                                 intercept_plot = intercept_plot,
    #                                 r2_plot = r2_plot, norm_rmse_plot = norm_rmse_plot,
    #                                 gap_fill_potential_plot = gap_fill_potential_plot, 
    #                                 gap_needed_to_fill_plot = gap_needed_to_fill_plot)

# #####################################################################
    # #stage 3 filling
    # stage3_dupfilling_results_all = stage3_dupfilling_all_voc()
    # build_summary_stage3(stage3_dupfilling_results_all)

    # # for item in stage3_dupfilling_results_all:
    # #     if isinstance(item['fill_case1'], list):
    # #         duplicate_fill_timeseries(dupfilling_results = item, stage_type = 'stage3', 
    # #                                   flag_type = 'Duplicate list')
    # #     elif item['fill_case1'] == 'No dup so stays vars':
    # #         print(f'No duplicate for {str(item['varname'])} in Stage 3.')
    # #     elif item['fill_case1'] == 'No vardata so filled with init dup':
    # #         print(f'No var data for {str(item['varname'])}, fill with init dup {str(item['filled_with_name'])} in Stage 3.')
    # #     elif item['fill_case1'] == 'No fill from dup due to low r squared':
    # #         print(f'No duplicate filling {str(item['varname'])} in Stage 3 due to low r squared.')
    # #     else:
    # #         duplicate_fill_timeseries(dupfilling_results = item, stage_type = 'stage3', flag_type = None)
    
    # stage3_tracerfilling_results_day_list, stage3_tracerfilling_results_night_list = stage3_tracer_filling_all_voc()
    # build_summary_stage3_tracer(results = stage3_tracerfilling_results_day_list, regime_type = 'day')
    # build_summary_stage3_tracer(results = stage3_tracerfilling_results_night_list, regime_type = 'night')
    
    # #####################################################################
    # Interpolation
    # data_for_df = {}
    # for spec in df_duplicates_and_tracers['Varname'].values:
    #     cache_file_stage3_tracer_fill_night = f'{CACHE_DIR}/stage3/{spec}_night_tracer_filled.pkl'
    #     with open(cache_file_stage3_tracer_fill_night, 'rb') as f:
    #         r = pickle.load(f)

    #     if isinstance(r['filled_voc'], list):
    #         filled_data = r['filled_voc'][-1]
    #     else:
    #         filled_data = r['filled_voc']

    #     data_for_df[spec] = filled_data
    # df_filled_data = pd.DataFrame(data_for_df)
    # print(df_filled_data)
    # # print(df_filled_data.columns)
    # # print(df_filled_data)
    # #df_filled_data['time_local'] = df_filled_data.index

    # df_filled_data.index = df_filled_data.index.tz_convert('UTC').tz_localize(None)
    # print(df_filled_data.index)

    # ds_nc = xr.open_dataset('/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/CampaignData_and_Merges/R0/CSL_MobileLab_Parked/merged/rev_15min/all_CSL_MobileLab_Parked_rev15min_iWASupdated.nc').load()
    # # print(ds_nc.coords['time_UTC'].values)
    # # print(ds_nc.dims)

    # df_filled_data.index.rename('time_UTC', inplace=True)
    # #We need to use timezone naive datetime index for the netCDF file
    # print(df_filled_data.index)

    # #Convert the dataframe of updated iWAS data to xarray Dataset
    # new_iwas_data = xr.Dataset.from_dataframe(df_filled_data)
    # #Align time coordinates to NetCDF time
    # new_iwas_data = new_iwas_data.reindex(time_UTC = ds_nc.time_UTC)

    # for var in new_iwas_data.data_vars:
    #     # Replace old values with new ones for each relevant variable
    #     ds_nc[var] = new_iwas_data[var]

    # save_ncfile_filled = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/CampaignData_and_Merges/R0/CSL_MobileLab_Parked/merged/rev_15min/all_CSL_MobileLab_Parked_rev15min_iWASupdated_filled.nc'

    # ds_nc.to_netcdf(save_ncfile_filled, format = 'NETCDF4', mode='w')
    # print('Saved netCDF file with updated iWAS measurements to: ', save_ncfile_filled)
    # ds_nc.close()

    # #Check that the format looks consistent with the original ds_nc
    # # ds_new = xr.open_dataset(save_ncfilename)
    # # print(ds_new)
    saved_ncfile = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/CampaignData_and_Merges/R0/CSL_MobileLab_Parked/merged/rev_15min/all_CSL_MobileLab_Parked_rev15min_iWASupdated_filled.nc'
    ds_filled_data = xr.open_dataset(saved_ncfile).load()
    df_filled_data = ds_filled_data.to_dataframe()
    ds_filled_data.close()
    
    df_filled_data.index = df_filled_data.index.tz_localize('UTC')
    df_filled_data['time_UTC'] = df_filled_data.index

    df_filled_data['time_local'] = df_filled_data['time_local'].dt.tz_localize('America/Denver')
    df_filled_data.set_index(['time_local'], inplace = True, drop = True)

    # Boundary Layer Height
    blh_filepath = dirpath+'Boundary_layer_height/BLH_Nell.csv'
    df_blh = pd.read_csv(blh_filepath)
    #Convert from Igor Pro Time to local time
    df_blh['time_local'] = (pd.to_datetime("1904-01-01") + pd.to_timedelta(df_blh['Time_start_15min_local'], unit="s"))
    df_blh = df_blh.rename(columns={'UDAQ_mixing_height_m_15min_avg':'BLH_Nell_m', 'UDAQ_mixing_height_m_15min_avg_smoothed':'BLH_Nell_m_smoothed'})

    df_blh['time_local'] = df_blh['time_local'].dt.tz_localize('America/Denver')

    df_blh.set_index(['time_local'], inplace = True, drop = False)
    df_blh = df_blh.drop(columns=['Time_start_15min_local'])
    df_filled_data = df_filled_data.join(df_blh, how="inner")   # only matching index values

    #rename the index to time_local_idx so it doesn't get mixed up with column assignment
    df_filled_data.index.name = 'time_local_idx'

    #Call function to interpolate:
    df_filled_data, df_interp = all_species_interp(df_filled_data)

    # get a ratio for jNO2 measured to TUV
    df_interp['jNO2_ratio'] = df_interp['jNO2_meas']/df_interp['jNO2']
    #level out the inf values and values that are too high for the jNO2 ratio
    msk = ((df_interp['jNO2_ratio'] ==np.inf)  | (df_interp['jNO2_ratio'] >10) )
    df_interp.loc[msk,'jNO2_ratio'] = 1.0

    # Get a count of the number of NaNs per day after interpolation, store in csv
    daily_nan_count = df_interp.isna().resample('D').sum()
    nan_count_csv_savepath = stage_data_dirs['stage4'] + 'csv/nan_count_interpolation_per_day_per_species_15min.csv'
    daily_nan_count.to_csv(nan_count_csv_savepath)

    #Save interpolated values into a CSV file for easy reading
    interpolated_values_filename = 'interpolated_values_15min'
    interpolation_savepath = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/CampaignData_and_Merges/R0/CSL_MobileLab_Parked/' + 'F0AM_filled/' + interpolated_values_filename + '.csv'
    df_interp.to_csv(interpolation_savepath)
    print('Saved CSV to:' + interpolation_savepath)

    #make interpolation plots for each species
    # for col in df_filled_data:
    #     # Calc number of points that are negative or Nans: 
    #     n_baddies= len([item for item in df_filled_data[col] if item <0 or np.isnan(item)]) 
    #     if n_baddies > 0: 
    #         plot_interps(df_filled_data, df_interp, col)

    all_days_interp_save_to_matlab(df_interp, 
                                     var_name = 'USOS')

    # subset_day_interp_save_to_matlab(df_interp, 
    #                                  var_name = 'USOS')
