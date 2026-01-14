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
import matplotlib.colors as mcolors

#Main datapaths
dirpath = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/'
#define path for Hawthorne data directory
savepath= dirpath + 'Merge_scripts/'
plots_savepath= dirpath + 'Merge_scripts/plots/'

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

VOC_secondary_selected = ['Acetone_WAS', 'Acetaldehyde_PTR','iPropylONO2_WAS', 'nPropylONO2_WAS', 'time_local']
final_VOC_list = VOC_secondary_selected[0:-1]

ml_data = xr.open_dataset('/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/CampaignData_and_Merges/R0/CSL_MobileLab_Parked/merged/rev_15min/all_CSL_MobileLab_Parked_rev15min_iWASupdated_formaldehydeupdated.nc')
ml_data_vars_subset = ml_data[VOC_secondary_selected]
df_ml_data = ml_data_vars_subset.to_dataframe()
df_ml_data.set_index('time_local', inplace=True)

# Select rows between 00:00 and 03:00 (inclusive)
voc_df = df_ml_data.between_time('00:00', '03:00')
# Column-wise mean
voc_means = voc_df.mean()
print(voc_means)

# #Put data into hour-of-day bin
# # df_ml_data['hour'] = df_ml_data.index.hour
# #Group by hour and take mean per species
# # hourly_means = df_ml_data.groupby('hour')[final_VOC_list].mean()
# # for col in final_VOC_list:
# #     plt.figure(figsize=(8,5))
# #     plt.plot(hourly_means.index, hourly_means[col], marker='o')
# #     plt.xlabel('Hour of Day')
# #     plt.ylabel(col + '(ppb)')
# #     # plt.title('Diurnal Profile of Total VOCs')
# #     plt.grid(True)
# #     plt.xticks(range(0,24))
# #     plt.show()

# df_ml_data['hour_2h'] = (df_ml_data.index.hour // 2) * 2
# #Group by hour and take mean per species
# bihourly_means = df_ml_data.groupby('hour_2h')[final_VOC_list].mean()
# for col in final_VOC_list:
#     plt.figure(figsize=(8,5))
#     plt.plot(bihourly_means.index, bihourly_means[col], marker='o')
#     plt.xlabel('Hour')
#     plt.ylabel(col + '(ppb)')
#     # plt.title('Diurnal Profile of Total VOCs')
#     plt.grid(True)
#     plt.xticks(range(0,24))
#     plt.show()

#     print(col + ' mean, 12-4 AM: \n', bihourly_means[col].loc[0:4].mean())

# print(bihourly_means)


def ozone_background_estimate():
    ml_ozone_vars_subset = ml_data[['O3_ppbv', 'time_local']]
    df_ml_ozone = ml_ozone_vars_subset.to_dataframe()
    df_ml_ozone.set_index('time_local', inplace=True)
    #smokefree dates
    smokefree_dates = (
        pd.date_range("2024-07-19", "2024-07-21").tolist() +
        pd.to_datetime([
            "2024-07-24",
            "2024-07-25",
            "2024-07-27",
            "2024-07-29"
        ]).tolist() +
        pd.date_range("2024-08-04", "2024-08-07").tolist() + 
        pd.date_range("2024-08-10", "2024-08-14").tolist() +
        pd.to_datetime(["2024-08-18"]).tolist()
    )

    smokedays = pd.to_datetime([
            "2024-07-31",
            "2024-08-02",
            "2024-08-08",
            "2024-08-09"
        ]).tolist()

    df_smokefree = df_ml_ozone[df_ml_ozone.index.normalize().isin(smokefree_dates)]
    df_smokedays = df_ml_ozone[df_ml_ozone.index.normalize().isin(smokedays)]
    #get median ozone between midnight and 5 AM
    df_smokefree_early_ozone = df_smokefree.between_time('00:00', '5:00')
    df_smokedays_early_ozone = df_smokedays.between_time('00:00', '5:00')
    df_alldays_early_ozone = df_ml_ozone.between_time('00:00', '5:00')
    # Column-wise mean
    smokefree_ozone_median = df_smokefree_early_ozone.median()
    smokedays_ozone_median = df_smokedays_early_ozone.median()
    alldays_ozone_median = df_alldays_early_ozone.median()
    print('smokefree O3 median (12-5 AM): ', smokefree_ozone_median)
    print('smokedays O3 median (12-5 AM): ', smokedays_ozone_median)
    print('All days O3 median (12-5 AM): ', alldays_ozone_median)


ozone_background_estimate()


# species_mean_voc = df_ml_data.mean(skipna=True)
# print('Mean VOC per species:\n', species_mean_voc, '\n')

# mean_total_voc = species_mean_voc.sum()
# print('Mean Total VOC:\n', mean_total_voc, '\n')

# species_sd_voc = df_ml_data.std(skipna=True)
# tvoc_variance = (species_sd_voc**2).sum()
# tvoc_std = np.sqrt(tvoc_variance)

# print('SD Total VOCs:\n', tvoc_std, '\n')

# species_median_voc = df_ml_data.median()
# print('Median VOC per species:\n', species_median_voc, '\n')
# median_total_voc = species_median_voc.sum()
# print('Median Total VOC:\n', median_total_voc, '\n')

# species_q25 = df_ml_data.quantile(0.25)
# species_q75 = df_ml_data.quantile(0.75)
# # Interquartile Range (IQR) Total VOC
# iqr_TVOC = species_q75.sum() - species_q25.sum()
# print('IQR Total VOC:\n', iqr_TVOC, '\n')


# #get percent contributions:
# pct_contributions = (species_mean_voc/mean_total_voc)*100
# pct_contributions_sorted = pct_contributions.sort_values(ascending=False)
