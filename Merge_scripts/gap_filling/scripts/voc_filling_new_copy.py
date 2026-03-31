#region: packages
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
import pickle
#endregion

#region: necessary filepaths
dirpath = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/'
hawthorne_data_dir = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/Hawthorne_data/data/from_Bart/'
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
#Mapping functions

def voc_mapping_setup(udaq):
    #get mappings to match UDAQ names to ML names
    mappings_filepath = dirpath + '/Hawthorne_data/mappings/manually_edited/UDAQ_Hawthorne_CRACMM_GEOSCHEM_CB6r5h_mapped_updated_01292026.csv'
    df_mapping_parameters = pd.read_csv(mappings_filepath)
    df_mapping_parameters = df_mapping_parameters.drop([0]) #drop Total NMVOCs

    voc_mapping = {
        (c1 if pd.notna(c1) else f'UDAQ_NoVar_{c2}'):(c2 if pd.notna(c2) else f'ML_NoVar_{c1}')
        for c1, c2 in zip(df_mapping_parameters['UDAQ_Variable'], df_mapping_parameters['USOS Mapping'])
    }

    #Sum mEthylToluene and pEthylToluene into 'mpEthyltoluene' and change mapping so that 'mpEthyltoluene' maps to 'x3_x4_EthylToluene_WAS'
    udaq['mpEthyltoluene'] = udaq['mEthyltoluene'] + udaq['pEthyltoluene']
    udaq = udaq.drop(['mEthyltoluene','pEthyltoluene'], axis = 1)
    voc_mapping['mpEthyltoluene'] = 'x3_x4_EthylToluene_WAS'
    del voc_mapping['mEthyltoluene']
    del voc_mapping['pEthyltoluene']
    print(voc_mapping)

    #List of all the WAS species that aren't in the mapping (We can't use any UDAQ data for filling holes but we'll need to interpolate the gaps values in order to use in F0AM)
    #Takes all the WAS species in the original Mobile Lab data and removes all those in the mapping, so that we know which ones we won't be cross-calibrating UDAQ data with
    #We can still see how well the CO tracer works to fill in holes
    was_species_list = []
    for colname in noaa.columns:
        if 'WAS' in colname:
            was_species_list.append(colname)
        else:
            pass
    was_species_not_in_mapping = list(set(was_species_list) - set(voc_mapping.values()))
    print(was_species_not_in_mapping)
    #['nPropylONO2_WAS', 'MACR_WAS', 'CF2Cl2_WAS', 'C2HCl3_WAS', 'Beta_Pinene_WAS', 'C2Cl4_WAS', 'iPropylONO2_WAS', 'CCl4_WAS', 'Acrolein_WAS', 'Acetone_WAS', 'CH3Br_WAS', 'Furan_WAS', 'CFCl3_WAS', 'Limonene_WAS', 'Benzene_WAS', 'Alpha_Pinene_WAS', 'CH2Cl2_WAS', 'Toluene_WAS']

    # Remove special cases:
    # For 'Alpha_Pinene_WAS', 'Beta_Pinene_WAS', and 'Limonene_WAS', we want to compare if Monoterpenes_PTR is approximately equal to 'Alpha_Pinene_WAS' + 'Beta_Pinene_WAS' + 'Limonene_WAS'. We will analyze this in its own function.
    # For Benzene and Toluene, we want to have the option to use either WAS or PTR data. 
    # For Acrolein and Acetone, we want to use PTR measurements to fill
    special_case_vocs = ['Alpha_Pinene_WAS', 'Beta_Pinene_WAS', 'Limonene_WAS', 'Benzene_WAS', 'Toluene_WAS', 'Acrolein_WAS', 'Acetone_WAS']
    for c in special_case_vocs:
        was_species_not_in_mapping.remove(c)
    #Add the other ML variables into our voc_mapping (excluding the special cases)
    for var in was_species_not_in_mapping:
        voc_mapping[f'UDAQ_NoVar_{var}'] = var

    return voc_mapping, special_case_vocs

############################################
#Helper functions
def mask_overlap(*series_list):

    mask = series_list[0].notna()

    for s in series_list[1:]:
        mask &= s.notna()

    return tuple(s.loc[mask] for s in series_list)

#Get overlap between two datasets
def get_overlap(x_data_overlap, y_data_overlap):
    df_vars = pd.DataFrame({'xdata': x_data_overlap, 'ydata': y_data_overlap})
    df_vars_overlap = df_vars.dropna() #Only include overlap with no NaNs
    return df_vars, df_vars_overlap

#Handle formaldehyde data
#If F0AM runs are changed to any interval other than 15 minutes, need to change the time_interval_formaldehyde value in Main run block's function call for process_all_vocs
def fill_formaldehyde_data(time_interval_formaldehyde):
    # read UDAQ Formaldehyde data, provided by and Bart (UDAQ), which is in UTC
    udaq_formaldehyde_load = hawthorne_data_dir + 'hw_zero_corrected_data_formaldehyde.csv'
    df_udaq_formaldehyde = pd.read_csv(udaq_formaldehyde_load, index_col='dt', parse_dates=True)
    df_udaq_formaldehyde.index = df_udaq_formaldehyde.index.tz_localize('UTC')

    df_udaq_formaldehyde_usos_only = df_udaq_formaldehyde.sort_index().loc['2024-07-15 00:00:00':'2024-08-18 23:59:00']
    df_udaq_formaldehyde_usos_only.index = df_udaq_formaldehyde_usos_only.index.rename('time_UTC')
    keep_colnames = ['H2CO_Corrected'] #This is the formaldehyde data that we want to use from UDAQ's 1 minute measurements
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
    df_avg_formaldehyde = df_avg_formaldehyde.rename(columns = {'H2CO_Corrected': 'Formaldehyde'}) #Change col name
    df_avg_formaldehyde = df_avg_formaldehyde.loc['2024-07-15 00:00:00':'2024-08-18 17:45:00']
    
    # #Now we have a dataframe df_avg_formaldehyde with columns 'ML Data' and 'UDAQ Data' with datetime index 2024-07-15 00:00:00 to 2024-08-18 17:45:00 (MDT)
    # df_avg_formaldehyde['ML Data']= noaa['HCHO_CRDS']

    # Set any inf, neg. inf, and negative values to NaN
    series_drop_formaldehyde = df_avg_formaldehyde['Formaldehyde'].replace([np.inf, -np.inf], np.nan).mask(df_avg_formaldehyde['Formaldehyde'] < 0, np.nan)
    
    return series_drop_formaldehyde

#Determine if the VOC species is:
#Non existent in the UDAQ data
#Non existent in the ML data
#Exists as a variable in the UDAQ data but is actually full of only NaNs
#Neither

def detect_voc_case(udaq_species_name, noaa_species_name, udaq):
    if 'UDAQ_NoVar' in udaq_species_name:
        case = 'UDAQ_NoVar'
        return case

    if 'ML_NoVar' in noaa_species_name:
        case = 'ML_NoVar'
        return case

    #UDAQ variable (column) exists but all NaN
    if isinstance(udaq, pd.DataFrame):
        if udaq[udaq_species_name].isna().all():
            case = 'UDAQ_All_NaN'
            return case
        else: 
            pass

    elif isinstance(udaq, pd.Series):
        if udaq.isna().all():
            case = 'UDAQ_All_NaN'
            return case
        else:
            pass
    else:
        case = 'Normal'
        return case

#Runs ODR fit for UDAQ and NOAA data
def run_odr(x_data_odr, y_data_odr):

    (df_vars, df_vars_overlap) = get_overlap(x_data_overlap = x_data_odr, 
                                             y_data_overlap = y_data_odr)

    #Get count of how many points are overlapping, and used in the ODR fitting
    points_considered_in_odr = len(df_vars_overlap['xdata'])

    #From the ODRpack documentation
    # Define the function we want to fit against
    x_data_overlap_np_arr = df_vars_overlap['xdata'].to_numpy()
    #xdata_2d = np.array([udaq_overlap_voc_np_arr])
    y_data_overlap_np_arr = df_vars_overlap['ydata'].to_numpy()
    #ydata_2d = np.array([ml_overlap_voc_np_arr])

    # print('xdata_2d shape: ', np.shape(udaq_overlap_voc_np_arr))
    # print('ydata shape: ', np.shape(ml_overlap_voc_np_arr))

    def linear_model(x, beta):
        return beta[0] + beta[1] * x
    
    #beta0 is an initial guess
    odr_fit_result = odr_fit(linear_model, x_data_overlap_np_arr, y_data_overlap_np_arr, [0.0, 1.0])

    odr_intercept, odr_slope = odr_fit_result.beta

    correction_eq = f"Correction applied: O_hat = {odr_intercept:.4f} + {odr_slope:.4f} * M"

    xvals_corrected = odr_intercept + odr_slope * df_vars['xdata']

    #overlap between initial and fitted data, no NaNs
    xcorrected_overlap = xvals_corrected.loc[xvals_corrected.index.intersection(df_vars_overlap.index)]
    return correction_eq, xvals_corrected, xcorrected_overlap, points_considered_in_odr

#Determine metrics of ODR fitting
def fitting_metrics(x_overlap, y_overlap):
    slope, intercept = np.polyfit(x_overlap, y_overlap, 1)
    rmse = np.sqrt(np.mean((x_overlap - y_overlap)**2))
    
    #R squared calculation
    r = np.corrcoef(x_overlap, y_overlap)[0, 1]
    r2 = r**2

    return slope, intercept, rmse, r2

#Determine if ODR fitting applied to correct the UDAQ data actually gets it closer to looking like the ML data
def odr_improves(metrics_init, metrics_corr, y_overlap):
    slope_i, intercept_i, rmse_i, r2_i = metrics_init
    slope_c, intercept_c, rmse_c, r2_c = metrics_corr

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

###############################################################
#Plots to perform quality control check: determine if ODR fit was helpful to correct UDAQ data
def scatterplots_for_comparing_init_and_corr_x_to_y(odr_fit_results, stage_type):
    x_voc_name = odr_fit_results['x_species_name']
    y_voc_name = odr_fit_results['y_species_name']

    x_init_overlap = odr_fit_results['x_data_init_overlap']
    x_corr_overlap = odr_fit_results['x_data_corr_overlap']
    y_overlap = odr_fit_results['y_data_init_overlap']

    slope_init = odr_fit_results['metrics_init'][0]
    slope_corr = odr_fit_results['metrics_corr'][0]

    intercept_init = odr_fit_results['metrics_init'][1]
    intercept_corr = odr_fit_results['metrics_corr'][1]

    rmse_init = odr_fit_results['metrics_init'][2]
    rmse_corr = odr_fit_results['metrics_corr'][2]

    r2_init = odr_fit_results['metrics_init'][3]
    r2_corr = odr_fit_results['metrics_corr'][3]

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
    ax[0].set_title('Initial')
    ax[0].set_xlabel(str(x_voc_name) + ' (ppb)')
    ax[0].set_ylabel(str(y_voc_name) + ' (ppb)')

    ax[0].text(0.05, 0.96, "Slope = " + str(round(slope_init, 3)), transform=ax[0].transAxes)
    ax[0].text(0.05, 0.94, "Intercept = " + str(round(intercept_init, 3)), transform=ax[0].transAxes)
    ax[0].text(0.05, 0.92, "R$^2$= " + str(round(r2_init, 3)), transform=ax[0].transAxes) 
    ax[0].text(0.05, 0.90, "RMSE:  " + str(round(rmse_init, 3)), transform=ax[0].transAxes)


    ax[1].scatter(x_corr_overlap, y_overlap, s=10, alpha=0.5)

    if np.nanmax(x_corr_overlap) < 0.1:
        step = 0.001
    else:
        step = 0.1

    xrange_corr = np.arange(0,np.nanmax(x_corr_overlap),step)
    
    ax[1].plot(xrange_corr, (slope_corr * xrange_corr + intercept_corr))
    ax[1].set_title('Corrected')
    ax[1].set_xlabel(str(x_voc_name) + ' (ppb)')
    ax[1].set_ylabel(str(y_voc_name) + ' (ppb)')

    ax[1].text(0.05, 0.96, "Slope = " + str(round(slope_corr, 3)), transform=ax[1].transAxes)
    ax[1].text(0.05, 0.94, "Intercept = " + str(round(intercept_corr, 3)), transform=ax[1].transAxes)
    ax[1].text(0.05, 0.92, "R$^2$= " + str(round(r2_corr, 3)), transform=ax[1].transAxes) 
    ax[1].text(0.05, 0.90, "RMSE:  " + str(round(rmse_corr, 3)), transform=ax[1].transAxes)

    if stage_type == 'stage3' or stage_type == 'stage3_5':
        plt.savefig(stage_data_dirs[stage_type] + 'plots/init_corr_scatterplot_comparison_' + str(y_voc_name) + '.png', dpi = 300)
        plt.show()
    else:
        plt.savefig(stage_data_dirs[stage_type] + 'plots/init_corr_scatterplot_comparison_' + str(x_voc_name) + '.png', dpi = 300)
        plt.show()

###############################################################
def gap_counter(data_before_change, data_after_change):
    #gap counter
    before_change_gap_count = data_before_change.isna().sum()
    after_change_gap_count = data_after_change.isna().sum()
    
    number_gaps_filled_by_change = before_change_gap_count - after_change_gap_count
    percent_of_gaps_filled_by_change = ((before_change_gap_count - after_change_gap_count) / before_change_gap_count) * 100

    return before_change_gap_count, after_change_gap_count, number_gaps_filled_by_change, percent_of_gaps_filled_by_change

def gap_filling(x_data_full, y_data_full):
    #Fill NaNs in ML Data with corrected UDAQ data
    merged_voc = y_data_full.fillna(x_data_full)
    return merged_voc

###############################################################
# Function to process each VOC
def process_one_voc_stage1(x_species_name, y_species_name, x_data, y_data):
    cache_file_stage1 = f'{CACHE_DIR}/{x_species_name}_stage1.pkl'
    
    # ---- load if already processed ----
    if os.path.exists(cache_file_stage1):
        with open(cache_file_stage1, 'rb') as f:
            return pickle.load(f)
        # If any of the results need to be modified, comment out the return pickle.load(f) line and replace with code below with alterations. 
        # This one changes a key previously named 'metrics_ml_udaq_corrected' into 'metrics_ml_udaq_corr'
        #     results = pickle.load(f)
        # print('results: ', results)
        # if 'metrics_ml_udaq_corrected' in results:
        #     results['metrics_ml_udaq_corr'] = results.pop('metrics_ml_udaq_corrected')
        # with open(cache_file, 'wb') as f:
        #     pickle.dump(results, f)


    case = detect_voc_case(x_species_name, y_species_name, x_data)
    print('Now processing VOC '+ str(x_species_name) + ' or ' + str(y_species_name) + ' in Stage 1.')
    print('Case: ', case)




    # -----------------------------
    # CASE: ML_NoVar
    # -----------------------------
    if case == 'ML_NoVar':
        x_init_data = x_data[x_species_name]
        y_init_data = pd.Series(np.nan, index = x_data.index)

        stage1_results = {
            'x_species_name': x_species_name,
            'y_species_name': y_species_name,
            'case': case,
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
            'y_data_init_full': y_init_data,
            'y_data_init_overlap': pd.Series(np.nan, index=udaq.index),
            'x_data_init_full': x_init_data,
            'x_data_init_overlap': pd.Series(np.nan, index=udaq.index),
            'x_data_corr_full': pd.Series(np.nan, index=udaq.index),
            'x_data_corr_overlap': pd.Series(np.nan, index=udaq.index)
        }
    # -----------------------------
    # CASE: UDAQ_NoVar
    # -----------------------------
    elif case == 'UDAQ_NoVar': 
        x_init_data = pd.Series(np.nan, index = y_data.index)
        y_init_data = y_data[y_species_name]

        stage1_results = {
            'x_species_name': x_species_name,
            'y_species_name': y_species_name,
            'case': case,
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
            'y_data_init_full': y_init_data,
            'y_data_init_overlap': pd.Series(np.nan, index=y_init_data.index),
            'x_data_init_full': x_init_data,
            'x_data_init_overlap': pd.Series(np.nan, index=y_init_data.index),
            'x_data_corr_full': pd.Series(np.nan, index=y_init_data.index),
            'x_data_corr_overlap': pd.Series(np.nan, index=y_init_data.index)
        }

    # -----------------------------
    # CASE: UDAQ_All_NaN
    # -----------------------------
    
    elif case == 'UDAQ_All_NaN':
        x_init_data = pd.Series(np.nan, index = y_data.index)
        y_init_data = y_data[y_species_name]

        stage1_results = {
            'x_species_name': x_species_name,
            'y_species_name': y_species_name,
            'case': case,
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
            'y_data_init_full': y_init_data,
            'y_data_init_overlap': pd.Series(np.nan, index=y_init_data.index),
            'x_data_init_full': x_init_data,
            'x_data_init_overlap': pd.Series(np.nan, index=y_init_data.index),
            'x_data_corr_full': pd.Series(np.nan, index=y_init_data.index),
            'x_data_corr_overlap': pd.Series(np.nan, index=y_init_data.index)
        }
    
    else:
        x_init_data = x_data[x_species_name]
        y_init_data = y_data[y_species_name]
        
        #Formaldehyde special case
        if x_species_name == 'Formaldehyde':
            series_drop_formaldehyde = fill_formaldehyde_data(time_interval_formaldehyde = 15*60)
            x_init_data = series_drop_formaldehyde
        
        #ODR fitting on UDAQ data
        (correction_eq, xvals_corrected, xcorrected_overlap, points_considered_in_odr) = run_odr(x_data_odr = x_init_data, 
                                                                                                    y_data_odr = y_init_data)

        #Metrics/Evaluation of ODR fitting
        (df_vars, df_vars_overlap) = get_overlap(x_data_overlap = x_init_data, 
                                                    y_data_overlap = y_init_data)
        
        metrics_init = fitting_metrics(df_vars_overlap['xdata'], df_vars_overlap['ydata'])
        metrics_corr  = fitting_metrics(xcorrected_overlap, df_vars_overlap['ydata'])

        #Quality control decision based off of if improvements are made by applying correction to UDAQ data
        (did_slope_improve, slope_distance_i, 
        slope_distance_c, slope_distance_improvement_val, 
        intercept_err_i, intercept_err_c,
        did_rmse_improve, rmse_percent_improvement,
        rmse_norm_i, rmse_norm_c,
        r2_at_least_half, score_i, score_c,
        apply_correction_based_off_score) = odr_improves(metrics_init, metrics_corr, df_vars_overlap['ydata'])

        stage1_results = {
            'x_species_name': x_species_name,
            'y_species_name': y_species_name,
            'case': case,
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
            'y_data_init_full': df_vars['ydata'],
            'y_data_init_overlap': df_vars_overlap['ydata'],
            'x_data_init_full': df_vars['xdata'],
            'x_data_init_overlap': df_vars_overlap['xdata'],
            'x_data_corr_full': xvals_corrected,
            'x_data_corr_overlap': xcorrected_overlap
        }
    
    # save so ODR never recomputes
    with open(cache_file_stage1, 'wb') as f:
        pickle.dump(stage1_results, f)

    return stage1_results
    
# Function to apply stage 2
def process_one_voc_stage2_udaq_ml(x_species_name, y_species_name, x_data):
    cache_file_stage1 = f'{CACHE_DIR}/{x_species_name}_stage1.pkl'
    cache_file_stage2 = f'{CACHE_DIR}/{x_species_name}_stage2.pkl'

    # ---- load if already processed ----
    if os.path.exists(cache_file_stage2):
        with open(cache_file_stage2, 'rb') as f:
            return pickle.load(f)
    else:
        with open(cache_file_stage1, 'rb') as f:
            r = pickle.load(f)

        udaq_voc_init_data_full = r['x_data_init_full']
        udaq_voc_corr_data_full = r['x_data_corr_full']
        noaa_voc_obs_data_full = r['y_data_init_full']

        case = detect_voc_case(x_species_name, y_species_name, udaq = x_data)

        if case == 'ML_NoVar':
            fill_case = 'Only UDAQ var so ml_filled_with_udaq is init UDAQ data'
            merged_voc = udaq_voc_init_data_full
            (before_change_gap_count, after_change_gap_count, 
            number_gaps_filled_by_change, percent_of_gaps_filled_by_change) = gap_counter(data_before_change = udaq_voc_init_data_full, 
                                                                                            data_after_change = udaq_voc_init_data_full)
            higher_error_params = np.nan
            

        elif case == 'UDAQ_NoVar':
            fill_case = 'Only ML var so ml_filled_with_udaq is ML data'
            merged_voc = noaa_voc_obs_data_full
            (before_change_gap_count, after_change_gap_count, 
            number_gaps_filled_by_change, percent_of_gaps_filled_by_change) = gap_counter(data_before_change = noaa_voc_obs_data_full, 
                                                                                            data_after_change = noaa_voc_obs_data_full)
            higher_error_params = np.nan
            
        elif case == 'UDAQ_All_NaN':
            fill_case = 'UDAQ data all NaNs so ml_filled_with_udaq is ML data'
            merged_voc = noaa_voc_obs_data_full
            (before_change_gap_count, after_change_gap_count, 
            number_gaps_filled_by_change, percent_of_gaps_filled_by_change) = gap_counter(data_before_change = noaa_voc_obs_data_full, 
                                                                                            data_after_change = noaa_voc_obs_data_full)
            higher_error_params = np.nan

        else:
            initial_udaq_params = {
            'slope': r['slope_distance_from_1_init'],
            'intercept': r['intercept_error_init'],
            'rmse': r['rmse_normalized_init']
            }

            corrected_udaq_params = {
            'slope': r['slope_distance_from_1_corr'],
            'intercept': r['intercept_error_corr'],
            'rmse': r['rmse_normalized_corr']
            }


            if r['r2_at_least_half'] is True and r['should_correction_be_applied_from_score_eval'] is True:
                fill_case = 'ML gaps filled by corr UDAQ data'

                # Set any negative values to NaN
                udaq_voc_corr_negs_to_nans = udaq_voc_corr_data_full.mask(udaq_voc_corr_data_full < 0, np.nan)
        
                merged_voc = gap_filling(x_data_full = udaq_voc_corr_negs_to_nans,
                                        y_data_full =  noaa_voc_obs_data_full)
                
                #identifies which parameters (slope, intercept, rmse) had higher error in the initial UDAQ data
                higher_error_params = [p for p in initial_udaq_params if initial_udaq_params[p] > corrected_udaq_params[p]]


                (before_change_gap_count, after_change_gap_count, 
                number_gaps_filled_by_change, percent_of_gaps_filled_by_change) = gap_counter(data_before_change = noaa_voc_obs_data_full, 
                                                                                            data_after_change = merged_voc)

            elif r['r2_at_least_half'] is True and r['should_correction_be_applied_from_score_eval'] is False:
                fill_case = 'ML gaps filled by init UDAQ data'

                #use the udaq initial data
                merged_voc = gap_filling(x_data_full = udaq_voc_init_data_full, 
                                        y_data_full = noaa_voc_obs_data_full)
            
                #identifies which parameters (slope, intercept, rmse) had higher error in the initial UDAQ data
                higher_error_params = [p for p in initial_udaq_params if initial_udaq_params[p] > corrected_udaq_params[p]]

                (before_change_gap_count, after_change_gap_count, 
                number_gaps_filled_by_change, percent_of_gaps_filled_by_change) = gap_counter(data_before_change = noaa_voc_obs_data_full, 
                                                                                            data_after_change = merged_voc)

            else:
                fill_case = 'R squared too low so no ML gaps filled by UDAQ data'
                merged_voc = noaa_voc_obs_data_full

                higher_error_params = np.nan
            
                #use the ML data for before and after (should have no gap changes)
                (before_change_gap_count, after_change_gap_count, 
                number_gaps_filled_by_change, percent_of_gaps_filled_by_change) = gap_counter(data_before_change = noaa_voc_obs_data_full, 
                                                                                            data_after_change = noaa_voc_obs_data_full)
                
        stage2_results = {
            'x_species_name': r['x_species_name'],
            'y_species_name': r['y_species_name'],
            'fill_case': fill_case,
            'higher_error_params': higher_error_params,
            'y_data_init_full': r['y_data_init_full'],
            'y_data_init_overlap': r['y_data_init_overlap'],
            'x_data_init_full': r['x_data_init_full'],
            'x_data_init_overlap': r['x_data_init_overlap'],
            'x_data_corr_full': r['x_data_corr_full'],
            'x_data_corr_overlap': r['x_data_corr_overlap'],
            'y_data_init_filled': merged_voc,
            'y_data_init_gap_count': before_change_gap_count, 
            'gap_count_left_after_y_data_init_filled': after_change_gap_count, 
            'gaps_filled_by_x': number_gaps_filled_by_change, 
            'percent_gaps_filled_by_x': percent_of_gaps_filled_by_change
        }
            
        with open(f'{CACHE_DIR}/{x_species_name}_stage2.pkl', 'wb') as f:
            pickle.dump(stage2_results, f)

        return stage2_results
    
def process_one_voc_stage3(x_species_name, y_species_name):
    #Compare CO (x variable) with VOC
    cache_file_stage2 = f'{CACHE_DIR}/{x_species_name}_stage2.pkl'
    cache_file_stage3 = f'{CACHE_DIR}/{x_species_name}_stage3.pkl'

    # # ---- load if already processed ----
    if os.path.exists(cache_file_stage3):
        with open(cache_file_stage3, 'rb') as f:
            return pickle.load(f)
    else:
        with open(cache_file_stage2, 'rb') as f:
            r = pickle.load(f)

    co_init_data = ml_co_raw
    voc_filled_data = r['y_data_init_filled']

    if voc_filled_data.isna().all():
        print('All vals are NaNs for ' + str(x_species_name))
        stage3_results = {
            'x_species_name': 'CO',
            'y_species_name': x_species_name,
            'voc_alt_name': y_species_name,
            'points_considered_in_odr': np.nan,
            'odr_eq_adj': np.nan, #equation used to correct the CO data
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
            'y_data_init_full': voc_filled_data,
            'y_data_init_overlap': pd.Series(np.nan, index=co_init_data.index),
            'x_data_init_full': co_init_data,
            'x_data_init_overlap': pd.Series(np.nan, index=co_init_data.index),
            'x_data_corr_full': pd.Series(np.nan, index=co_init_data.index),
            'x_data_corr_overlap': pd.Series(np.nan, index=co_init_data.index)
        }
        with open(f'{CACHE_DIR}/{x_species_name}_stage3.pkl', 'wb') as f:
            pickle.dump(stage3_results, f)

        return stage3_results
    
    #If the VOC values aren't all NaNs, then apply ODR fitting with the CO and VOC data
    (correction_eq, xvals_corrected, xcorrected_overlap, points_considered_in_odr) = run_odr(x_data_odr = co_init_data, 
                                                                                            y_data_odr = voc_filled_data)

    #Metrics/Evaluation of ODR fitting
    (df_vars, df_vars_overlap) = get_overlap(x_data_overlap = co_init_data, 
                                            y_data_overlap = voc_filled_data)
    
    metrics_init = fitting_metrics(df_vars_overlap['xdata'], df_vars_overlap['ydata'])
    metrics_corr  = fitting_metrics(xcorrected_overlap, df_vars_overlap['ydata'])

    #Quality control decision based off of if improvements are made by applying correction to UDAQ data
    (did_slope_improve, slope_distance_i, 
    slope_distance_c, slope_distance_improvement_val, 
    intercept_err_i, intercept_err_c,
    did_rmse_improve, rmse_percent_improvement,
    rmse_norm_i, rmse_norm_c,
    r2_at_least_half, score_i, score_c,
    apply_correction_based_off_score) = odr_improves(metrics_init, metrics_corr, df_vars_overlap['ydata'])


    stage3_results = {
        'x_species_name': 'CO',
        'y_species_name': x_species_name,
        'voc_alt_name': y_species_name,
        'points_considered_in_odr': points_considered_in_odr,
        'odr_eq_adj': correction_eq, #equation used to correct the CO data
        'metrics_corr': metrics_corr,
        'did_slope_improve': did_slope_improve, 
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
        'y_data_init_full': df_vars['ydata'],
        'y_data_init_overlap': df_vars_overlap['ydata'],
        'x_data_init_full': df_vars['xdata'],
        'x_data_init_overlap': df_vars_overlap['xdata'],
        'x_data_corr_full': xvals_corrected,
        'x_data_corr_overlap': xcorrected_overlap
    }

    with open(f'{CACHE_DIR}/{x_species_name}_stage3.pkl', 'wb') as f:
            pickle.dump(stage3_results, f)

    return stage3_results

def process_one_voc_stage3_5(x_species_name, y_species_name):
    #Compare CO (x variable) with VOC
    cache_file_stage2 = f'{CACHE_DIR}/{x_species_name}_stage2.pkl'
    #cache_file_stage3 = f'{CACHE_DIR}/{x_species_name}_stage3.pkl'
    cache_file_stage2_benzene = f'{CACHE_DIR}/Benzene_stage2.pkl'
    cache_file_stage2_toluene = f'{CACHE_DIR}/Toluene_stage2.pkl'

    # ---- load if already processed ----
    # if os.path.exists(cache_file_stage3):
    #     with open(cache_file_stage3, 'rb') as f:
    #         return pickle.load(f)
    #else:
    with open(cache_file_stage2, 'rb') as f:
        r = pickle.load(f)

    print('Loaded Stage2 Pickle file')

    co_init_data = ml_co_raw
    voc_filled_data = r['y_data_init_filled']

    mask_co = (
        ~ml_co_raw.isna() & 
        ~voc_filled_data.isna()
    )

    if voc_filled_data.isna().all():
        print('All vals are NaNs for ' + str(x_species_name))
        stage3_5_results = {
            'udaq_species_name': x_species_name,
            'ml_species_name': y_species_name,
            'points_considered_in_odr_CO_tracer': np.nan,
            'y_pred_CO_tracer': np.nan,
            'r_CO_tracer': np.nan,
            'r2_CO_tracer': np.nan,
            'rmse_CO_tracer': np.nan,
            'points_considered_in_odr_CO_Benzene_tracer': np.nan,
            'y_pred_CO_Benzene_tracer': np.nan,
            'r_CO_Benzene_tracer': np.nan,
            'r2_CO_Benzene_tracer': np.nan,
            'rmse_CO_Benzene_tracer': np.nan,
            'points_considered_in_odr_CO_Toluene_tracer': np.nan,
            'y_pred_CO_Toluene_tracer': np.nan,
            'r_CO_Toluene_tracer': np.nan,
            'r2_CO_Toluene_tracer': np.nan,
            'rmse_CO_Toluene_tracer': np.nan,
            'points_considered_in_odr_CO_Benzene_Toluene_tracer': np.nan,
            'y_pred_CO_Benzene_Toluene_tracer': np.nan,
            'r_CO_Benzene_Toluene_tracer': np.nan,
            'r2_CO_Benzene_Toluene_tracer': np.nan,
            'rmse_CO_Benzene_Toluene_tracer': np.nan,     
            'y_data_init_full': pd.Series(np.nan, index = co_init_data.index),
            'y_data_CO_overlap': pd.Series(np.nan, index = co_init_data.index),
            'y_data_CO_Benzene_overlap': pd.Series(np.nan, index = co_init_data.index),
            'y_data_CO_Toluene_overlap': pd.Series(np.nan, index = co_init_data.index),
            'y_data_CO_Benzene_Toluene_overlap': pd.Series(np.nan, index = co_init_data.index)
        }

        print('Dumping NaN VOC data into Stage 3.5 Pickle file')

        with open(f'{CACHE_DIR}/{x_species_name}_stage3_5.pkl', 'wb') as f:
            pickle.dump(stage3_5_results, f)

        return stage3_5_results
    
    print('Loading initial Benzene and Toluene Pickle files')
    with open(cache_file_stage2_benzene, 'rb') as f:
        b = pickle.load(f)
    with open(cache_file_stage2_toluene, 'rb') as f:
        t = pickle.load(f)

    benzene_series = b['y_data_init_filled']
    toluene_series = t['y_data_init_filled']
    voc_series = r['y_data_init_filled']

    mask_co_only = (~ml_co_raw.isna() & ~voc_series.isna())

    mask_benz_co = (~ml_co_raw.isna() & ~benzene_series.isna() & ~voc_series.isna())

    mask_tol_co = (~ml_co_raw.isna() &~toluene_series.isna() &~voc_series.isna())

    mask_total = (~ml_co_raw.isna() &~benzene_series.isna() & ~toluene_series.isna() & ~voc_series.isna())

    #make tracer stack CO:
    co_only = np.vstack([ml_co_raw[mask_co_only].values])
    print('co_only shape: ', co_only.shape)
    voc_series_co_overlap = voc_series[mask_co_only].values
    print('voc_series_co_overlap shape: ', voc_series_co_overlap.shape)

    #make tracer stack of CO + Benz:
    benz_co_stack = np.vstack([ml_co_raw[mask_benz_co].values, benzene_series[mask_benz_co].values])
    print('benz_co_stack shape: ', benz_co_stack.shape)
    voc_series_benz_co_overlap = voc_series[mask_benz_co].values
    print('voc_series_benz_co_overlap shape: ', voc_series_benz_co_overlap.shape)

    #make tracer stack of CO + Tol:
    tol_co_stack = np.vstack([ml_co_raw[mask_tol_co].values, toluene_series[mask_tol_co].values])
    print('tol_co_stack shape: ', tol_co_stack.shape)
    voc_series_tol_co_overlap = voc_series[mask_tol_co].values
    print('voc_series_tol_co_overlap shape: ', voc_series_tol_co_overlap.shape)

    #make tracer stack of CO, Benz, Tol:
    total_tracer_stack = np.vstack([ml_co_raw[mask_total].values, benzene_series[mask_total].values, toluene_series[mask_total].values])
    print('total_tracer_stack shape: ', total_tracer_stack.shape)
    voc_series_total_overlap = voc_series[mask_total].values
    print('voc_series_total_overlap shape: ', voc_series_total_overlap.shape)

    tracer_overlap_x = [co_only, benz_co_stack, tol_co_stack, total_tracer_stack]
    voc_overlap_y = [voc_series_co_overlap, voc_series_benz_co_overlap, voc_series_tol_co_overlap, voc_series_total_overlap]

    y_pred_store_all = []
    r_store_all = []
    r2_store_all = []
    rmse_store_all = []
    points_considered_store_all = []

    print('Starting ODR for different tracers (CO, Benzene, Toluene)')
    for xdata_tracer, ydata_voc in zip(tracer_overlap_x, voc_overlap_y):
        #Get count of how many points are overlapping, and used in the ODR fitting
        points_considered_in_odr = xdata_tracer.shape[1]
        print('points_considered_in_odr: ', points_considered_in_odr)
        points_considered_store_all.append(points_considered_in_odr)

        #From the ODRpack documentation
        # Define the function we want to fit against
        x_data_overlap_np_arr = xdata_tracer
        print('xdata_tracer type:', type(xdata_tracer))
        print('x_data_overlap_np_arr shape: ', x_data_overlap_np_arr.shape)
        #xdata_2d = np.array([udaq_overlap_voc_np_arr])
        y_data_overlap_np_arr = ydata_voc
        print('y_data_overlap_np_arr shape: ', y_data_overlap_np_arr.shape)
        #ydata_2d = np.array([ml_overlap_voc_np_arr])

        # print('xdata_2d shape: ', np.shape(udaq_overlap_voc_np_arr))
        # print('ydata shape: ', np.shape(ml_overlap_voc_np_arr))

        def linear_model(x, beta):
            x = np.atleast_2d(x)
            if x.shape[0] != len(beta) - 1:
                x = x.T
            return np.dot(beta[:-1], x) + beta[-1]
            #return beta[0] + beta[1] * x
        # def linear_model_two_factor(x, beta):
        #     return beta[0] + beta[1] * x[1] + beta[2] * x[2]
        # def linear_model_three_factor(x, beta):
        #     return beta[0] + beta[1] * x[1] + beta[2] * x[2] + beta[3] * x[3]
        beta0 = np.ones(x_data_overlap_np_arr.shape[0] + 1)

        test_shape = linear_model(x_data_overlap_np_arr, beta0)
        print('test_shape shape: ', test_shape.shape)
        #beta0 is an initial guess
        odr_fit_result = odr_fit(linear_model, x_data_overlap_np_arr, y_data_overlap_np_arr, beta0)

        #beta = odr_fit_result.beta
        y_pred = linear_model(x_data_overlap_np_arr, odr_fit_result.beta)
        print('y_pred shape: ', y_pred.shape)
        r = np.corrcoef(y_data_overlap_np_arr, y_pred)[0,1]
        r2 = r**2
        rmse = np.sqrt(np.mean((y_data_overlap_np_arr - y_pred)**2))

        y_pred_store_all.append(y_pred)
        r_store_all.append(r)
        r2_store_all.append(r2)
        rmse_store_all.append(rmse)

        ######################

        #if CO is the only tracer
        # if xdata_tracer.shape[0] == 1:
        #     test_shape = linear_model_one_factor(x_data_overlap_np_arr, [0.0, 1.0])
        #     print(test_shape.shape)
        #     #beta0 is an initial guess
        #     odr_fit_result = odr_fit(linear_model_one_factor, x_data_overlap_np_arr, y_data_overlap_np_arr, [0.0, 1.0])

        #     beta = odr_fit_result.beta
        #     y_pred = linear_model_one_factor(x_data_overlap_np_arr, beta)
        #     r = np.corrcoef(y_data_overlap_np_arr, y_pred)[0,1]
        #     r2 = r**2
        #     rmse = np.sqrt(np.mean((y_data_overlap_np_arr - y_pred)**2))

        #     y_pred_store_all.append(y_pred)
        #     r_store_all.append(r)
        #     r2_store_all.append(r2)
        #     rmse_store_all.append(rmse)

        # #if CO + another tracer
        # elif xdata_tracer.shape[0] == 2:
        #     odr_fit_result = odr_fit(linear_model_two_factor, x_data_overlap_np_arr, y_data_overlap_np_arr, [0.0, 0.1, 0.1])

        #     beta = odr_fit_result.beta
        #     y_pred = linear_model_two_factor(beta, x_data_overlap_np_arr)
        #     r = np.corrcoef(y_data_overlap_np_arr, y_pred)[0,1]
        #     r2 = r**2
        #     rmse = np.sqrt(np.mean((y_data_overlap_np_arr - y_pred)**2))

        #     y_pred_store_all.append(y_pred)
        #     r_store_all.append(r)
        #     r2_store_all.append(r2)
        #     rmse_store_all.append(rmse)

        # #if CO + benz + tol tracer
        # elif xdata_tracer.shape[0] == 3:
        #     odr_fit_result = odr_fit(linear_model_three_factor, x_data_overlap_np_arr, y_data_overlap_np_arr, [0.0, 0.5, 0.5, 0.001])

        #     beta = odr_fit_result.beta
        #     y_pred = linear_model_three_factor(beta, x_data_overlap_np_arr)
        #     r = np.corrcoef(y_data_overlap_np_arr, y_pred)[0,1]
        #     r2 = r**2
        #     rmse = np.sqrt(np.mean((y_data_overlap_np_arr - y_pred)**2))

        #     y_pred_store_all.append(y_pred)
        #     r_store_all.append(r)
        #     r2_store_all.append(r2)
        #     rmse_store_all.append(rmse)
###################################################
    stage3_5_results = {
        'udaq_species_name': x_species_name,
        'ml_species_name': y_species_name,
        'points_considered_in_odr_CO_tracer': points_considered_store_all[0],
        'y_pred_CO_tracer': y_pred_store_all[0],
        'r_CO_tracer': r_store_all[0],
        'r2_CO_tracer': r2_store_all[0],
        'rmse_CO_tracer': rmse_store_all[0],
        'points_considered_in_odr_CO_Benzene_tracer': points_considered_store_all[1],
        'y_pred_CO_Benzene_tracer': y_pred_store_all[1],
        'r_CO_Benzene_tracer': r_store_all[1],
        'r2_CO_Benzene_tracer': r2_store_all[1],
        'rmse_CO_Benzene_tracer': rmse_store_all[1],
        'points_considered_in_odr_CO_Toluene_tracer': points_considered_store_all[2],
        'y_pred_CO_Toluene_tracer': y_pred_store_all[2],
        'r_CO_Toluene_tracer': r_store_all[2],
        'r2_CO_Toluene_tracer': r2_store_all[2],
        'rmse_CO_Toluene_tracer': rmse_store_all[2],
        'points_considered_in_odr_CO_Benzene_Toluene_tracer': points_considered_store_all[3],
        'y_pred_CO_Benzene_Toluene_tracer': y_pred_store_all[3],
        'r_CO_Benzene_Toluene_tracer': r_store_all[3],
        'r2_CO_Benzene_Toluene_tracer': r2_store_all[3],
        'rmse_CO_Benzene_Toluene_tracer': rmse_store_all[3],        
        'y_data_init_full': voc_series,
        'y_data_CO_overlap': voc_series_co_overlap,
        'y_data_CO_Benzene_overlap': voc_series_benz_co_overlap,
        'y_data_CO_Toluene_overlap': voc_series_tol_co_overlap,
        'y_data_CO_Benzene_Toluene_overlap': voc_series_total_overlap
    }

    print('Dumping most VOC data into Stage 3.5 Pickle file')
    with open(f'{CACHE_DIR}/{x_species_name}_stage3_5.pkl', 'wb') as f:
        pickle.dump(stage3_5_results, f)

    return stage3_5_results
##########################################
#Loops through mappings to run for all VOCs
def process_all_vocs_stage1(x_data, y_data, mapping_species):
    stage1_results_all = []
    for x_species_name, y_species_name in mapping_species.items():
        stage1_results = process_one_voc_stage1(x_species_name, y_species_name, x_data, y_data)
        stage1_results_all.append(stage1_results)
    return stage1_results_all

#Build summary CSV file to save information needed
def build_summary_stage1(stage1_results_all):
    rows = []
    for r in stage1_results_all:
        rows.append({
        'voc_udaq_name': r['x_species_name'],
        'voc_noaa_name': r['y_species_name'],
        'points_considered_in_odr': r['points_considered_in_odr'],
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
    stage1_summary_savepath = stage_data_dirs['stage1'] + 'csv_storage/stage1_summary.csv'
    pd.DataFrame(rows).to_csv(stage1_summary_savepath, index=True)
    print('Saved summary to: ' + str(stage1_summary_savepath))
#Run Stage 2 for all VOCs
def process_all_vocs_stage2_udaq_ml(x_data, mapping_species):
    stage2_results_all = []
    for x_species_name, y_species_name in mapping_species.items():
        stage2_results = process_one_voc_stage2_udaq_ml(x_species_name, y_species_name, x_data)
        stage2_results_all.append(stage2_results)
    return stage2_results_all

def build_summary_stage2(stage2_results_all):
    rows = []
    for r in stage2_results_all:
        rows.append({
        'voc_udaq_name': r['x_species_name'],
        'voc_noaa_name': r['y_species_name'],
        'fill_case': r['fill_case'],
        'higher_error_params': r['higher_error_params'],
        'ml_initial_gap_count': r['y_data_init_gap_count'], 
        'gap_count_left_after_ml_filled': r['gap_count_left_after_y_data_init_filled'], 
        'gaps_filled_by_udaq': r['gaps_filled_by_x'], 
        'percent_gaps_filled_by_udaq': r['percent_gaps_filled_by_x']
        })
    stage2_summary_savepath = stage_data_dirs['stage2'] + 'csv_storage/stage2_summary.csv'
    pd.DataFrame(rows).to_csv(stage2_summary_savepath, index=True)
    print('Saved summary to: ' + str(stage2_summary_savepath))

def process_all_vocs_stage3(mapping_species):
    stage3_results_all = []
    for x_species_name, y_species_name in mapping_species.items():
        stage3_results = process_one_voc_stage3(x_species_name, y_species_name)
        stage3_results_all.append(stage3_results)
    # print(stage3_results_all)
    return stage3_results_all

def process_all_vocs_stage3_5(mapping_species):
    stage3_5_results_all = []
    for x_species_name, y_species_name in mapping_species.items():
        print('Processing Stage 3.5 for: ', x_species_name)
        stage3_5_results = process_one_voc_stage3_5(x_species_name, y_species_name)
        stage3_5_results_all.append(stage3_5_results)
    #print('stage3_5_results_all in process_all_vocs_stage3_5 function', stage3_5_results_all)
    return stage3_5_results_all

def build_summary_stage3(stage3_results_all):
    rows = []
    for r in stage3_results_all:
        rows.append({
        'voc_udaq_name': r['y_species_name'],
        'voc_noaa_name': r['voc_alt_name'],
        'points_considered_in_odr': r['points_considered_in_odr'],
        'odr_eq_adj': r['odr_eq_adj'],
        'metrics_corr_slope': r['metrics_corr'][0],
        'metrics_corr_intercept': r['metrics_corr'][1],
        'metrics_corr_rmse': r['metrics_corr'][2],
        'metrics_corr_r2': r['metrics_corr'][3],
        'slope_distance_from_1_corr': r[ 'slope_distance_from_1_corr'],
        'intercept_error_corr': r['intercept_error_corr'],
        'rmse_normalized_corr': r['rmse_normalized_corr'],
        'r2_at_least_half': r['r2_at_least_half'],
        'total_slope_intercept_rmse_score_corr': r['total_slope_intercept_rmse_score_corr']
        })
    stage3_summary_savepath = stage_data_dirs['stage3'] + 'csv_storage/stage3_summary.csv'
    pd.DataFrame(rows).to_csv(stage3_summary_savepath, index=True)
    print('Saved summary to: ' + str(stage3_summary_savepath))

def build_summary_stage3_5(stage3_5_results_all):
    rows = []
    for r in stage3_5_results_all:
        rows.append({
            'voc_udaq_name': r['udaq_species_name'],
            'voc_noaa_name': r['ml_species_name'],
            'points_considered_in_odr_CO_tracer': r['points_considered_in_odr_CO_tracer'],
            # 'y_pred_CO_tracer': r['y_pred_CO_tracer'],
            'r_CO_tracer': r['r_CO_tracer'], 
            'r2_CO_tracer': r['r2_CO_tracer'],
            'rmse_CO_tracer': r['rmse_CO_tracer'],
            'points_considered_in_odr_CO_Benzene_tracer': r['points_considered_in_odr_CO_Benzene_tracer'],
            # 'y_pred_CO_Benzene_tracer': r['y_pred_CO_Benzene_tracer'],
            'r_CO_Benzene_tracer': r['r_CO_Benzene_tracer'],
            'r2_CO_Benzene_tracer': r['r2_CO_Benzene_tracer'],
            'rmse_CO_Benzene_tracer': r['rmse_CO_Benzene_tracer'],
            'points_considered_in_odr_CO_Toluene_tracer': r['points_considered_in_odr_CO_Toluene_tracer'],
            # 'y_pred_CO_Toluene_tracer': r['y_pred_CO_Toluene_tracer'],
            'r_CO_Toluene_tracer': r['r_CO_Toluene_tracer'],
            'r2_CO_Toluene_tracer': r['r2_CO_Toluene_tracer'],
            'rmse_CO_Toluene_tracer': r['rmse_CO_Toluene_tracer'],
            'points_considered_in_odr_CO_Benzene_Toluene_tracer': r['points_considered_in_odr_CO_Benzene_Toluene_tracer'],
            # 'y_pred_CO_Benzene_Toluene_tracer': r['y_pred_CO_Benzene_Toluene_tracer'],
            'r_CO_Benzene_Toluene_tracer': r['r_CO_Benzene_Toluene_tracer'],
            'r2_CO_Benzene_Toluene_tracer': r['r2_CO_Benzene_Toluene_tracer'],
            'rmse_CO_Benzene_Toluene_tracer': r['rmse_CO_Benzene_Toluene_tracer'],
            # 'y_data_init_full': r['y_data_init_full'],
            # 'y_data_CO_overlap': r['y_data_CO_overlap'],
            # 'y_data_CO_Benzene_overlap': r['y_data_CO_Benzene_overlap'],
            # 'y_data_CO_Toluene_overlap': r['y_data_CO_Toluene_overlap'],
            # 'y_data_CO_Benzene_Toluene_overlap': r['y_data_CO_Benzene_Toluene_overlap']
        })
    stage3_5_summary_savepath = stage_data_dirs['stage3_5'] + 'csv_storage/stage3_5_summary.csv'
    pd.DataFrame(rows).to_csv(stage3_5_summary_savepath, index=True)
    print('Saved summary to: ' + str(stage3_5_summary_savepath))
################################################################
#Plotting functions
def plot_ml_gap_fill_stage2(x_species_name):
    cache_file_stage2 = f'{CACHE_DIR}/{x_species_name}_stage2.pkl'
    # ---- load if already processed ----
    if os.path.exists(cache_file_stage2):
        with open(cache_file_stage2, 'rb') as f:
            r = pickle.load(f)

    #From Stage2 cache, get ML initial data        
    y_initial = r['y_data_init_full']
    #From Stage2 cache, get newly filled ML data
    y_filled = r['y_data_init_filled']

    #Formats the labels for plots without super long varname
    if x_species_name.startswith('UDAQ_NoVar_'):
        x_species_name = x_species_name[len('UDAQ_NoVar_'):]
    
    fig, (ax1, ax2) = plt.subplots(2,1, figsize = (16,8), tight_layout=True)
    xlim_start_jul = pd.to_datetime('2024-07-15 00:00:00').tz_localize('America/Denver')
    xlim_end_jul = pd.to_datetime('2024-07-31 23:45:00').tz_localize('America/Denver')
    xlim_start_aug = pd.to_datetime('2024-08-01 00:00:00').tz_localize('America/Denver')
    xlim_end_aug = pd.to_datetime('2024-08-18 23:45:00').tz_localize('America/Denver')

    #ax1 is the first row of subplot, for July only
    valid_points_initial = ~np.isnan(y_initial)
    ax1.plot(y_initial.index[valid_points_initial], y_initial[valid_points_initial], linestyle='solid', color = 'm', marker = 'x', label = f'Obs (gaps = {r['y_data_init_gap_count']})')
    valid_points_filled = ~np.isnan(y_filled)
    ax1.plot(y_filled.index[valid_points_filled], y_filled[valid_points_filled], linestyle='solid', color = 'y', marker = '.', label = f'Filled (gaps = {r['gap_count_left_after_y_data_init_filled']})', alpha = 0.7)

    #Set x ticks
    tz_mdt = noaa.index.tz #this time zone should be in Mountain Daylight Time
    ax1.xaxis.set_major_locator(mdates.DayLocator(tz=tz_mdt))
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d', tz=tz_mdt))
    # Minor ticks: every 3 hours
    ax1.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21], tz=tz_mdt))
    # Rotate and format tick labels
    ax1.tick_params(axis='x', which='major')
    ax1.tick_params(axis='x', which='minor', length=3, color='gray')
    

    ax1.set_ylabel(x_species_name + ' (ppb)')
    #ax1.set_xlabel('Date')
    ax1.margins(x=0)
    ax1.set_xlim([xlim_start_jul, xlim_end_jul])

    ax1.legend(loc = 'upper right')

    #ax2 is the second row of subplot, for August only
    ax2.plot(y_initial.index[valid_points_initial], y_initial[valid_points_initial], linestyle='solid', color = 'm', marker = 'x', label = f'Obs (gaps = {r['y_data_init_gap_count']})')
    ax2.plot(y_filled.index[valid_points_filled], y_filled[valid_points_filled], linestyle='solid', color = 'y', marker = '.', label = f'Filled (gaps = {r['gap_count_left_after_y_data_init_filled']})', alpha = 0.7)


    #Set x ticks
    ax2.xaxis.set_major_locator(mdates.DayLocator(tz=tz_mdt))
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d', tz=tz_mdt))
    # Minor ticks: every 3 hours
    ax2.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21], tz=tz_mdt))
    # Rotate and format tick labels
    ax2.tick_params(axis='x', which='major')
    ax2.tick_params(axis='x', which='minor', length=3, color='gray')

    ax2.set_ylabel(x_species_name + ' (ppb)')
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

    plt.savefig(stage_data_dirs['stage2'] + 'plots/timeseries_full/filled_' + str(x_species_name) + '_comparison_july_aug_timeseries.png', dpi =300)
    plt.show()
def scatterplot_obs_pred(odr_fit_results):
    voc_obs = odr_fit_results['y_data_init_overlap']
    voc_pred = odr_fit_results['x_data_corr_overlap']

    voc_name = odr_fit_results['y_species_name']

    xmin = min(voc_obs.min(), voc_pred.min())
    xmax = max(voc_obs.max(), voc_pred.max())


    plt.figure(figsize=(5,5))
    plt.scatter(voc_pred, voc_obs, s=10, alpha=0.3)
    plt.plot([xmin, xmax], [xmin, xmax], '--', color='black', label='1:1 line')
    plt.xlabel(voc_name + ' Predicted (ppb)')
    plt.ylabel(voc_name + ' Obs (ppb)')
    plt.tight_layout()
    #plt.savefig(dirpath + '/Merge_scripts/calibration_adjustments/co_tracer_plots/scatterplot_co_vs_voc/ml_co_with_'+ str(vocname_udaq) + '_scatterplot.png', dpi =300)
    plt.show()

def scatterplot_vocs_to_aromatics(odr_fit_results):
    cache_file_stage2_benzene = f'{CACHE_DIR}/Benzene_stage2.pkl'
    cache_file_stage2_toluene = f'{CACHE_DIR}/Toluene_stage2.pkl'

    with open(cache_file_stage2_benzene, 'rb') as f:
            b = pickle.load(f)
    with open(cache_file_stage2_toluene, 'rb') as f:
            t = pickle.load(f)
    benzene_filled = b['y_data_init_filled']
    toluene_filled = t['y_data_init_filled']

    voc_filled = odr_fit_results['y_data_init_full']
    voc_name = odr_fit_results['y_species_name']

    df_vars, df_vars_overlap = get_overlap(benzene_filled, voc_filled)
    benzene_overlap = df_vars_overlap['xdata']
    voc_overlap = df_vars_overlap['ydata']
    r_val = np.corrcoef(benzene_overlap, voc_overlap)[0, 1]
    print('Benzene Correlation with ', voc_name, ' :', r_val)


    plt.figure(figsize=(5,5))
    plt.scatter(benzene_overlap, voc_overlap, s=10, alpha=0.3)
    plt.xlabel('Benzene (ppb)')
    plt.ylabel(voc_name + ' (ppb)')
    plt.tight_layout()
    #plt.savefig(dirpath + '/Merge_scripts/calibration_adjustments/co_tracer_plots/scatterplot_co_vs_voc/ml_co_with_'+ str(vocname_udaq) + '_scatterplot.png', dpi =300)
    plt.show()

    df_vars, df_vars_overlap = get_overlap(toluene_filled, voc_filled)
    toluene_overlap = df_vars_overlap['xdata']
    voc_overlap = df_vars_overlap['ydata']
    r_val = np.corrcoef(toluene_overlap, voc_overlap)[0, 1]
    print('Toluene Correlation with ', voc_name, ' :', r_val)

    plt.figure(figsize=(5,5))
    plt.scatter(toluene_overlap, voc_overlap, s=10, alpha=0.3)
    plt.xlabel('Toluene (ppb)')
    plt.ylabel(voc_name + ' (ppb)')
    plt.tight_layout()
    #plt.savefig(dirpath + '/Merge_scripts/calibration_adjustments/co_tracer_plots/scatterplot_co_vs_voc/ml_co_with_'+ str(vocname_udaq) + '_scatterplot.png', dpi =300)
    plt.show()

################################################################
#region:scrap
#CO Scatterplot Comparison
# def qc_plots_for_co_correction(odr_fit_results):
#     udaq_voc_name = odr_fit_results['voc_udaq_name']
#     ml_voc_name = odr_fit_results['voc_noaa_name']
    
#     udaq_init_overlap = odr_fit_results['udaq_co_no_corr_overlap']
#     udaq_corr_overlap = odr_fit_results['udaq_co_with_corr_overlap']
#     noaa_overlap = odr_fit_results['noaa_voc_obs_overlap']

#     slope_ml_udaq_init = odr_fit_results['metrics_ml_udaq_init'][0]
#     slope_ml_udaq_corr = odr_fit_results['metrics_ml_udaq_corr'][0]

#     intercept_ml_udaq_init = odr_fit_results['metrics_ml_udaq_init'][1]
#     intercept_ml_udaq_corr = odr_fit_results['metrics_ml_udaq_corr'][1]

#     rmse_ml_udaq_init = odr_fit_results['metrics_ml_udaq_init'][2]
#     rmse_ml_udaq_corr = odr_fit_results['metrics_ml_udaq_corr'][2]

#     r2_ml_udaq_init = odr_fit_results['metrics_ml_udaq_init'][3]
#     r2_ml_udaq_corr = odr_fit_results['metrics_ml_udaq_corr'][3]

#     if udaq_init_overlap.isna().all():
#         print('All vals are NaNs for ' + str(udaq_voc_name))

#     else:
#         #scatterplot of ML vs UDAQ initial
#         fig, ax = plt.subplots(1, 2, figsize=(10,10), tight_layout=True)
#         ax[0].scatter(udaq_init_overlap, noaa_overlap, s=10, alpha=0.5)

#         #To draw regression line, we need a continuous line. Since some species are in ppt, we need
#         #to select an appropriate scale to the step
#         if np.nanmax(udaq_init_overlap) < 0.1:
#             step = 0.001
#         else:
#             step = 0.1

#         xrange_init = np.arange(0,np.nanmax(udaq_init_overlap), step)

#         ax[0].plot(xrange_init, (slope_ml_udaq_init * xrange_init + intercept_ml_udaq_init))
#         ax[0].set_title('Initial')
#         ax[0].set_xlabel('CO (ppb)')
#         ax[0].set_ylabel(str(udaq_voc_name) + ' (ppb)')

#         ax[0].text(0.05, 0.96, "Slope = " + str(round(slope_ml_udaq_init, 3)), transform=ax[0].transAxes)
#         ax[0].text(0.05, 0.94, "Intercept = " + str(round(intercept_ml_udaq_init, 3)), transform=ax[0].transAxes)
#         ax[0].text(0.05, 0.92, "R$^2$= " + str(round(r2_ml_udaq_init, 3)), transform=ax[0].transAxes) 
#         ax[0].text(0.05, 0.90, "RMSE:  " + str(round(rmse_ml_udaq_init, 3)), transform=ax[0].transAxes)


#         ax[1].scatter(udaq_corr_overlap, noaa_overlap, s=10, alpha=0.5)

#         if np.nanmax(udaq_corr_overlap) < 0.1:
#             step = 0.001
#         else:
#             step = 0.1

#         xrange_corr = np.arange(0,np.nanmax(udaq_corr_overlap),step)
        
#         ax[1].plot(xrange_corr, (slope_ml_udaq_corr * xrange_corr + intercept_ml_udaq_corr))
#         ax[1].set_title('Corrected')
#         ax[1].set_xlabel('CO (ppb)')
#         ax[1].set_ylabel(str(udaq_voc_name) + ' (ppb)')

#         ax[1].text(0.05, 0.96, "Slope = " + str(round(slope_ml_udaq_corr, 3)), transform=ax[1].transAxes)
#         ax[1].text(0.05, 0.94, "Intercept = " + str(round(intercept_ml_udaq_corr, 3)), transform=ax[1].transAxes)
#         ax[1].text(0.05, 0.92, "R$^2$= " + str(round(r2_ml_udaq_corr, 3)), transform=ax[1].transAxes) 
#         ax[1].text(0.05, 0.90, "RMSE:  " + str(round(rmse_ml_udaq_corr, 3)), transform=ax[1].transAxes)

#         plt.savefig(stage3_data_dir + '/plots/ml_udaq_scatterplot_comparison_' + str(udaq_voc_name) + '.png', dpi = 300)
#         plt.show()
# def plot_co_scatterplot(udaq_species_name, co_data):
#     cache_file_stage2 = f'{CACHE_DIR}/{udaq_species_name}_stage2.pkl'
#     # ---- load if already processed ----
#     if os.path.exists(cache_file_stage2):
#         with open(cache_file_stage2, 'rb') as f:
#             r = pickle.load(f)

#     #From Stage2 cache, get newly filled ML data
#     ml_filled = r['ml_filled_with_udaq']

#     df_vars, df_vars_overlap = get_overlap(udaq_voc_init_data = co_data, noaa_voc_init_data = ml_filled)
#     co_overlap_vals = df_vars_overlap['UDAQ Data']
#     ml_filled_voc_overlap = df_vars_overlap['ML Data']

#     #ODR fitting
#     (udaq_correction_eq, udaq_vals_corrected, udaq_corrected_overlap, points_considered_in_odr) = run_odr(udaq_voc_init_data = ml_co_raw, noaa_voc_init_data = ml_filled)
#     metrics_co_voc = fitting_metrics(udaq_corrected_overlap, ml_filled_voc_overlap)

#     plt.figure(figsize=(5,5))
#     plt.scatter(udaq_corrected_overlap, ml_filled_voc_overlap, s=10, alpha=0.3)
    
#     #To draw regression line, we need a continuous line. Since some species are in ppt, we need
#     #to select an appropriate scale to the step
#     if np.nanmax(udaq_corrected_overlap) < 0.1:
#         step = 0.001
#     else:
#         step = 0.1

#     xrange_init = np.arange(0,np.nanmax(udaq_corrected_overlap), step)
#     plt.plot(xrange_init, (slope_ml_udaq_init * xrange_init + intercept_ml_udaq_init))
#     plt.xlabel('ML CO (ppb)')
#     plt.ylabel(udaq_species_name + ' (ppb)')
#     plt.tight_layout()
#     plt.savefig(stage2_data_dir + '/plots/co_filled_voc_scatterplot_comparison_' + str(udaq_species_name) + '.png', dpi = 300)
#     plt.show()

# def plot_scatter(voc_name):
#     with open(f"cache/{voc_name}.pkl", "rb") as f:
#         r = pickle.load(f)
    
#     plt.scatter(r["inst1"], r["inst2"], label="raw", alpha=0.5)
#     plt.scatter(r["inst1"], r["inst2_adj"], label="adjusted", alpha=0.5)
#     plt.legend()

#Instead of applying a CO tracer for Alpha pinene, Beta pinene, and Limonene, use Monoterpenes_PTR data to adjust and fill values
#endregion
###############################################################
#Main run block
if __name__ == "__main__":
    #Load UDAQ and NOAA ML VOC files respectively
    udaq_f = dirpath + '/Hawthorne_data/data/script_output/hawthorne_udaq_all_vocs_15min_timezone_carbon_number_updated.csv'
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
    
    formaldehyde_f = hawthorne_data_dir + 'hawthorne_udaq_Formaldehyde_15min_reindexed_timezone_updated.csv'
    df_formaldehyde = pd.read_csv(formaldehyde_f, index_col='time_local', parse_dates=True)
    df_formaldehyde.index = df_formaldehyde.index.tz_localize(None)
    df_formaldehyde.index = df_formaldehyde.index.tz_localize('America/Denver')
    df_formaldehyde['Formaldehyde'] = df_formaldehyde['Formaldehyde'].replace([np.inf, -np.inf], np.nan).mask(df_formaldehyde['Formaldehyde'] < 0, np.nan)
    df_formaldehyde = df_formaldehyde.loc['2024-07-15 00:00:00':'2024-08-18 17:45:00']

    ozone_f = hawthorne_data_dir +'hawthorne_udaq_o3_2024_15min_reindexed_timezone_updated.csv'
    df_o3 = pd.read_csv(ozone_f, index_col='time_local', parse_dates=True)
    df_o3.index = df_o3.index.tz_localize(None)
    df_o3.index = df_o3.index.tz_localize('America/Denver')
    df_o3 = df_o3.loc['2024-07-15 00:00:00':'2024-08-18 17:45:00']
    
    terpenes_f = hawthorne_data_dir + 'hawthorne_udaq_isoprene_alpha_beta_pinene_07152024_08012024_15min_reindexed_timezone_updated.csv'
    df_terpenes =pd.read_csv(terpenes_f, index_col='time_local', parse_dates=True)
    df_terpenes.index = df_terpenes.index.tz_localize(None)
    df_terpenes.index = df_terpenes.index.tz_localize('America/Denver')
    #Set index to only include 2024-07-15 00:00:00 to 2024-08-18 17:45:00 MDT
    df_terpenes = df_terpenes.loc['2024-07-15 00:00:00':'2024-08-18 17:45:00']

    no_f = hawthorne_data_dir + 'hawthorne_udaq_no_07152024_08012024_15min_reindexed_timezone_updated.csv'
    df_no =pd.read_csv(no_f, index_col='time_local', parse_dates=True)
    df_no.index = df_no.index.tz_localize(None)
    df_no.index = df_no.index.tz_localize('America/Denver')
    #Set index to only include 2024-07-15 00:00:00 to 2024-08-18 17:45:00 MDT
    df_no = df_no.loc['2024-07-15 00:00:00':'2024-08-18 17:45:00']

    noy_f = hawthorne_data_dir + 'hawthorne_udaq_noy_07152024_08012024_15min_reindexed_timezone_updated.csv'
    df_noy =pd.read_csv(noy_f, index_col='time_local', parse_dates=True)
    df_noy.index = df_noy.index.tz_localize(None)
    df_noy.index = df_noy.index.tz_localize('America/Denver')
    #Set index to only include 2024-07-15 00:00:00 to 2024-08-18 17:45:00 MDT
    df_noy = df_noy.loc['2024-07-15 00:00:00':'2024-08-18 17:45:00']

    no2_f = hawthorne_data_dir + 'hawthorne_udaq_no2_07152024_08012024_15min_reindexed_timezone_updated.csv'
    df_no2 =pd.read_csv(no2_f, index_col='time_local', parse_dates=True)
    df_no2.index = df_no2.index.tz_localize(None)
    df_no2.index = df_no2.index.tz_localize('America/Denver')
    #Set index to only include 2024-07-15 00:00:00 to 2024-08-18 17:45:00 MDT
    df_no2 = df_no2.loc['2024-07-15 00:00:00':'2024-08-18 17:45:00']

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

    #Give explicit variable to CO that will be used in functions. Turn any negatives into NaNs.
    ml_co_raw = noaa['CO_Piccaro']
    #There seems to be a data issue with the CO Picarro during 07/22 15:00 MDT to 07/24 00:00 MDT so mask that range into NaNs
    ml_co_raw.loc['2024-07-22 11:00:00':'2024-07-22 18:00:00'] = np.nan
    ml_co_raw.loc['2024-07-23 10:30:00':'2024-07-23 12:00:00'] = np.nan

    print(ml_co_raw)
    print(noaa['NO_LIF'])
    

    #Call function to set up VOC mappings for UDAQ:ML data
    (voc_mapping, special_case_vocs) = voc_mapping_setup(udaq)
    print(voc_mapping)

#######################################################
    # STAGE 1

    print('Starting Stage 1')

    stage1_results_all = process_all_vocs_stage1(x_data = udaq, 
                                                 y_data = noaa,
                                                 mapping_species = voc_mapping)
    

    # # Make plots comparing the correlation between ML and UDAQ data, only plot if we have valid ML Data and UDAQ Data
    # # Comment out if unnecessary to save or see plots
    # for result_row in stage1_results_all:
    #     if result_row['y_data_init_full'].isna().all():
    #         print('No ML Data for ' + str(result_row['x_species_name']) + ' or ' + str(result_row['y_species_name']) + ' in Stage 1.')

    #     elif result_row['x_data_init_full'].isna().all():
    #         print('No UDAQ Data for ' + str(result_row['x_species_name']) + ' or ' + str(result_row['y_species_name']) + ' in Stage 1.')
        
    #     else:
    #         scatterplots_for_comparing_init_and_corr_x_to_y(odr_fit_results = result_row, stage_type = 'stage1')

    build_summary_stage1(stage1_results_all)

#######################################################
    # STAGE 2
    print('Starting Stage 2')
    stage2_results_all = process_all_vocs_stage2_udaq_ml(x_data = udaq, mapping_species = voc_mapping)

    build_summary_stage2(stage2_results_all)

    # for result_row in stage2_results_all:
    #     vocname = result_row['x_species_name']
    #     plot_ml_gap_fill_stage2(x_species_name = vocname)

#######################################################
    # STAGE 3
    print('Starting Stage 3')
    stage3_results_all = process_all_vocs_stage3(voc_mapping)

    # for result_row in stage3_results_all:
    #     if result_row['y_data_init_full'].isna().all():
    #         print('No measurements for ' + str(result_row['y_species_name']) + ' or ' + str(result_row['voc_alt_name']) + ' in Stage 3.')
    #     else:    
    #         vocname = result_row['y_species_name']
    #         print(vocname)
            #scatterplots_for_comparing_init_and_corr_x_to_y(odr_fit_results = result_row, stage_type = 'stage3')
            
            #scatterplot_obs_pred(odr_fit_results = result_row)
            #scatterplot_vocs_to_aromatics(odr_fit_results = result_row)

    # build_summary_stage3(stage3_results_all)

    print('Starting Stage 3.5')
    stage3_5_results_all = process_all_vocs_stage3_5(voc_mapping)
    

    # for result_row in stage3_5_results_all:
    #     if result_row['y_data_init_full'].isna().all():
    #         print('No measurements for ' + str(result_row['y_species_name']) + ' or ' + str(result_row['voc_alt_name']) + ' in Stage 3.5.')
    #     else:    
    #         vocname = result_row['y_species_name']
    #         print(vocname)
    #         #scatterplots_for_comparing_init_and_corr_x_to_y(odr_fit_results = result_row, stage_type = 'stage3_5')


    build_summary_stage3_5(stage3_5_results_all)

