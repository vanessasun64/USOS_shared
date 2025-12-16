import numpy as np 
import os
import sys
import xarray as xr
import pandas as pd
from collections import OrderedDict

global_scripts_path ='/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/global_scripts'
sys.path.insert(0,global_scripts_path)
from dirpath import filepath_source
dirpath = filepath_source('CHPC')
resaved_dir = dirpath + 'USOS_Campaign_analysis/resaved_data/'

#NOT VALID METHOD OF GETTING 1 HOUR MERGES -> can't avg from 30 min merges, must do 1 hr merges again
def resample_hourly_species(rev_30min_merge_filename,species_name,filename):
    #from 30 min merges
    all_days_filepath = dirpath + 'CampaignData_and_Merges/R0/CSL_MobileLab_Parked/merged/rev_30min/' + rev_30min_merge_filename + '.nc'
    all_days_filepath_load = xr.open_dataset(all_days_filepath)
    df_alldays = all_days_filepath_load.to_dataframe()
    df_alldays.reset_index(inplace=True)
    df_alldays.set_index('time_local', inplace=True, drop=False)

    df_species = df_alldays[[species_name]].copy()
    df_species = df_species.rename(columns={'NO2_LIF':'Mobile Lab NO2'})
    df_species_resample = df_species.resample('1h').mean()

    new_start_time = pd.Timestamp('2024-07-15 00:00:00')
    new_end_time = pd.Timestamp('2024-08-18 23:00:00')
    new_index = pd.date_range(start=new_start_time, end=new_end_time, freq='1h')
    df_species = df_species_resample.reindex(new_index)
    df_species.index.name = 'time_local'
    df_species['time_UTC'] = df_species.index + pd.to_timedelta(6, unit='h')

    savepath = resaved_dir + filename
    df_species.to_csv(savepath)
    print('Saved to:' + savepath)
    
resample_hourly_species(
    rev_30min_merge_filename = 'all_CSL_MobileLab_Parked_rev30minv4',
    species_name = 'NO2_LIF',
    filename = 'no2_resampled_hourly.csv')