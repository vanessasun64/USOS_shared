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


def resave_hourly_ozone_as_csv(rev_1hr_merge_filename,filename):
    #from 1 hr merges
    all_days_filepath = dirpath + 'CampaignData_and_Merges/R0/CSL_MobileLab_Parked/merged/rev_1hr/' + rev_1hr_merge_filename + '.nc'
    all_days_filepath_load = xr.open_dataset(all_days_filepath)
    df_alldays = all_days_filepath_load.to_dataframe()
    df_alldays.reset_index(inplace=True)
    df_alldays.set_index('time_local', inplace=True, drop=False)

    df_ozone = df_alldays[['O3_ppbv']].copy()
    df_ozone = df_ozone.rename(columns={'O3_ppbv':'Mobile Lab O3'})

    new_start_time = pd.Timestamp('2024-07-15 00:00:00')
    new_end_time = pd.Timestamp('2024-08-18 23:00:00')
    new_index = pd.date_range(start=new_start_time, end=new_end_time, freq='1h')
    df_ozone = df_ozone.reindex(new_index)
 
    df_ozone.index.name = 'time_local'
    df_ozone['time_UTC'] = df_ozone.index + pd.to_timedelta(6, unit='h')

    savepath = resaved_dir + filename
    df_ozone.to_csv(savepath)
    print('Saved to:' + savepath)


resave_hourly_ozone_as_csv(
    rev_1hr_merge_filename = 'all_CSL_MobileLab_Parked_rev1hr_iWASupdated',
    filename = 'ozone_hourly.csv'
    )