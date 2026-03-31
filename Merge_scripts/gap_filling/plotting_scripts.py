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

#region: Load UDAQ and NOAA ML VOC files respectively
#Time index is every 15 min from 07/14/2024 18:00:00 to 8/18/2024 17:45:00 Local Time
noaa_f= dirpath + '/CampaignData_and_Merges/R0/CSL_MobileLab_Parked/merged/rev_15min/all_CSL_MobileLab_Parked_rev15min_iWASupdated.nc'

ds = xr.open_dataset(noaa_f)
noaa = ds.to_dataframe()
noaa = noaa.set_index(['time_local'])
#localize time zone to MDT
noaa.index = noaa.index.tz_localize('America/Denver')
#Set index to only span from 2024-07-15 00:00:00 to 2024-08-18 17:45:00
noaa = noaa.loc['2024-07-15 00:00:00':'2024-08-18 17:45:00']

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
#endregion

ml_udaq_initial_and_filled_ml_load = dirpath + '/Merge_scripts/ml_filling_with_udaq_data/calibration_data/odr_ml_udaq_initial_and_filled_ml_with_udaq_adjusted.csv'
ml_udaq_initial_and_filled_ml_overlap_load = dirpath + '/Merge_scripts/ml_filling_with_udaq_data/calibration_data/odr_ml_udaq_initial_and_filled_ml_with_udaq_adjusted_overlap.csv'
df_ml_udaq_initial_and_filled_ml = pd.read_csv(ml_udaq_initial_and_filled_ml_load, index_col='time_local', parse_dates=True)
df_ml_udaq_initial_and_filled_ml_overlap = pd.read_csv(ml_udaq_initial_and_filled_ml_overlap_load, index_col='time_local', parse_dates=True)

odr_additional_info_path = dirpath + '/Merge_scripts/ml_filling_with_udaq_data/calibration_data/odr_cal_adj_voc_info.csv'
odr_additional_info = pd.read_csv(odr_additional_info_path, index_col=[0])

#For some reason, pandas is reading the time_local as UTC so we convert it back to the UTC-6 time zone
df_ml_udaq_initial_and_filled_ml.index = df_ml_udaq_initial_and_filled_ml.index.tz_localize(None)
df_ml_udaq_initial_and_filled_ml.index = df_ml_udaq_initial_and_filled_ml.index.tz_localize('America/Denver')
df_ml_udaq_initial_and_filled_ml_overlap.index = df_ml_udaq_initial_and_filled_ml_overlap.index.tz_localize(None)
df_ml_udaq_initial_and_filled_ml_overlap.index = df_ml_udaq_initial_and_filled_ml_overlap.index.tz_localize('America/Denver')

def scatterplot_odr(udaq_initial_data, udaq_adjusted_data, ml_data, udaq_voc_name, ml_voc_name, slope_initial, slope_calibrated, intercept_initial, intercept_calibrated, rsquared_initial, rsquared_calibrated, rmse_initial, rmse_calibrated):
    #scatterplot of ML vs UDAQ initial
    fig, ax = plt.subplots(1, 2, figsize=(10,10), tight_layout=True)
    ax[0].scatter(udaq_initial_data,  ml_data, s=10, alpha=0.5)
    if np.nanmax(udaq_initial_data) < 0.1:
        step = 0.001
    else:
        step = 0.1

    xrange_initial = np.arange(0,np.nanmax(udaq_initial_data),step)

    ax[0].plot(xrange_initial, (slope_initial * xrange_initial + intercept_initial))
    ax[0].set_title('Initial UDAQ')
    ax[0].set_xlabel('UDAQ' + str(udaq_voc_name) + ' (ppb)')
    ax[0].set_ylabel('ML ' + str(ml_voc_name) + ' (ppb)')
    # #ADD ANNOTATION OF EQUATION
    # lin_equation_init = 'ML =' + slope_ml_udaq_init + ' UDAQ + ' + intercept_ml_udaq_init
    # #ADD ANNOTATION OF R SQUARED VALUE
    # #rsquared_ml_udaq_init

    

    ax[0].text(0.05, 0.96, "Slope = " + str(slope_initial), transform=ax[0].transAxes)
    ax[0].text(0.05, 0.94, "Intercept = " + str(intercept_initial), transform=ax[0].transAxes)
    ax[0].text(0.05, 0.92, "R^2= " + str(rsquared_initial), transform=ax[0].transAxes) 
    ax[0].text(0.05, 0.90, "RMSE:  " + str(rmse_initial), transform=ax[0].transAxes) 


    ax[1].scatter(udaq_adjusted_data, ml_data, s=10, alpha=0.5)

    if np.nanmax(udaq_initial_data) < 0.1:
        step = 0.001
    else:
        step = 0.1

    xrange_calibrated = np.arange(0,np.nanmax(udaq_adjusted_data),step)

    ax[1].plot(xrange_calibrated, (slope_calibrated * xrange_calibrated + intercept_calibrated))
    ax[1].set_title('Corrected UDAQ')
    ax[1].set_xlabel('UDAQ ' + str(udaq_voc_name) + ' (ppb)')
    ax[1].set_ylabel('ML ' + str(ml_voc_name) + ' (ppb)')

    ax[1].text(0.05, 0.96, "Slope = " + str(slope_calibrated), transform=ax[1].transAxes)
    ax[1].text(0.05, 0.94, "Intercept = " + str(intercept_calibrated), transform=ax[1].transAxes)
    ax[1].text(0.05, 0.92, "R^2= " + str(rsquared_calibrated), transform=ax[1].transAxes) 
    ax[1].text(0.05, 0.90, "RMSE:  " + str(rmse_calibrated), transform=ax[1].transAxes) 


    # #ADD ANNOTATION OF EQUATION
    # lin_equation_corr = 'ML =' + slope_ml_udaq_corr + ' UDAQ + ' + intercept_ml_udaq_corr
    # #ADD ANNOTATION OF R SQUARED VALUE
    # #rsquared_ml_udaq_corr
    plt.savefig(dirpath + '/Merge_scripts/ml_filling_with_udaq_data/adjustment_plots/odr_fitting/scatterplots/ml_udaq_initial_and_calibrated_scatterplot_comparison_'+ str(udaq_voc_name) + '2.png', dpi =300)
    plt.show()

    #Plot 2: ML to UDAQ relationship after applying ODR fit

# for i in range(0, len(df_ml_udaq_initial_and_filled_ml.columns), 4):
#     chunk = df_ml_udaq_initial_and_filled_ml.iloc[:, i:i+4]
#     print(chunk)
print(df_ml_udaq_initial_and_filled_ml.columns)
for i in range(0, len(df_ml_udaq_initial_and_filled_ml.columns), 4):
    cols = df_ml_udaq_initial_and_filled_ml.columns[i:i+4]
    
    ml_col = cols[0]
    udaq_init_col = cols[1]
    udaq_adj_col  = cols[2]
    
    ml_species_match = re.search(r'NOAA_(.*?)_Initial', ml_col)
    ml_species_name = ml_species_match.group(1)
    udaq_species_match = re.search(r"UDAQ_(.*?)_Initial", udaq_init_col)
    udaq_species_name = udaq_species_match.group(1)

    if ml_species_name == 'Formaldehyde':
        ml_species_name = 'HCHO_CRDS'

    idx_of_species_in_additional_info = odr_additional_info.index[odr_additional_info['NOAA_Species_Name'] == ml_species_name]


    ml_udaq_initial_slope = odr_additional_info.loc[idx_of_species_in_additional_info, 'Slope_ML_to_Initial_UDAQ'].iloc[0]
    ml_udaq_calibrated_slope = odr_additional_info.loc[idx_of_species_in_additional_info, 'Slope_ML_to_Calibrated_UDAQ'].iloc[0]

    ml_udaq_initial_intercept = odr_additional_info.loc[idx_of_species_in_additional_info, 'Intercept_ML_to_Initial_UDAQ'].iloc[0]
    ml_udaq_calibrated_intercept = odr_additional_info.loc[idx_of_species_in_additional_info, 'Intercept_ML_to_Calibrated_UDAQ'].iloc[0]

    ml_udaq_initial_rsquared = odr_additional_info.loc[idx_of_species_in_additional_info,'R_Squared_ML_to_Initial_UDAQ'].iloc[0]
    ml_udaq_calibrated_rsquared = odr_additional_info.loc[idx_of_species_in_additional_info,'R_Squared_ML_to_Calibrated_UDAQ'].iloc[0]

    ml_udaq_initial_rmse = odr_additional_info.loc[idx_of_species_in_additional_info,'RMSE_Initial'].iloc[0]
    ml_udaq_calibrated_rmse = odr_additional_info.loc[idx_of_species_in_additional_info,'RMSE_Calibrated'].iloc[0]
    print(ml_udaq_initial_rmse)
    
    # extract data
    ml_initial   = df_ml_udaq_initial_and_filled_ml[ml_col]
    udaq_initial   = df_ml_udaq_initial_and_filled_ml[udaq_init_col]
    udaq_adjusted  = df_ml_udaq_initial_and_filled_ml[udaq_adj_col]
    
    print(ml_species_name, udaq_species_name)

    # print('length of udaq initial:', len(udaq_init_col))
    # print('length of udaq initial * slope + intercept:', len(ml_udaq_initial_slope * udaq_initial + ml_udaq_initial_intercept))
    if ml_initial.isna().all():
        print('No ML Data')
    elif udaq_initial.isna().all():
        print('No UDAQ Data')
    else:
        scatterplot_odr(udaq_initial_data = udaq_initial, 
                        udaq_adjusted_data = udaq_adjusted,
                        ml_data = ml_initial, 
                        udaq_voc_name = udaq_species_name, 
                        ml_voc_name = ml_species_name,
                        slope_initial = ml_udaq_initial_slope,
                        slope_calibrated = ml_udaq_calibrated_slope,
                        intercept_initial = ml_udaq_initial_intercept,
                        intercept_calibrated = ml_udaq_calibrated_intercept,
                        rsquared_initial = ml_udaq_initial_rsquared,
                        rsquared_calibrated = ml_udaq_calibrated_rsquared,
                        rmse_initial = ml_udaq_initial_rmse,
                        rmse_calibrated = ml_udaq_calibrated_rmse
                        )
    


    # if 'Formaldehyde' in col:
    #     df_ml_udaq_initial_and_filled_ml[col]

    # elif col.startswith('Filled_') and col.endswith('_ML_with_UDAQ_Adjusted'):
    #     ml_vocname = col[len('Filled_'):-len('_ML_with_UDAQ_Adjusted')]
    #     print('filled_voc_name: ', ml_vocname)

    #     filled_colname = f'Filled_{ml_vocname}_ML_with_UDAQ_Adjusted'

    #     ml_data = df_ml_udaq_initial_and_filled_ml[col]

    # elif (col.startswith('UDAQ_') and col.endswith('_Initial')) :
    #     udaq_vocname = col[len('UDAQ_'):-len('_Initial')]
    #     udaq_init_data = df_ml_udaq_initial_and_filled_ml[col]

    #     if df_ml_udaq_initial_and_filled_ml[udaq_vocname].isna().all():
    #         print(udaq_vocname, ' has no measurements (UDAQ Data).')
    #     else:
    #         scatterplot_odr(udaq_data, ml_data, udaq_voc_name, noaa_voc_name)

    # elif (col.startswith('UDAQ_') and col.endswith('_Adjusted')):
    #     udaq_vocname = col[len('UDAQ_'):-len('_Adjusted')]
    #     udaq_adj_data = df_ml_udaq_initial_and_filled_ml[col]

    #     if df_ml_udaq_initial_and_filled_ml[udaq_vocname].isna().all():
    #         print(udaq_vocname, ' has no measurements (UDAQ Data).')
    #     else:
    #         scatterplot_odr(udaq_data, ml_data, udaq_voc_name, noaa_voc_name)