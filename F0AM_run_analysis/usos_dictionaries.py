# For saving dictionaries of different dates for USOS Campaign

import pandas as pd
import numpy as np

dirpath = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/'
interp_vals_filepath = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/CampaignData_and_Merges/R0/CSL_MobileLab_Parked/F0AM_filled/interpolated_values_15min.csv'
df_interp_vals = pd.read_csv(interp_vals_filepath, index_col = 'time_local', parse_dates = True)

whole_campaign_dict = {}
voc_available_days_dict = {}

def mda8_calculate(date_limit, save_filename):
    df_o3_obs = df_interp_vals['O3_ppbv']
    #resample from 15 minute intervals to hourly
    mda8_data = df_o3_obs.resample('1h').mean()

    
    # Calculate the 8-hour rolling average
    rolling_avg_name = '8hr_rolling_avg_O3_obs'
    mda8_data[rolling_avg_name] = mda8_data.rolling(window=8, min_periods=6).mean()
    #print("mda8_data['8hr_rolling_avg_O3_obs']: ", mda8_data['8hr_rolling_avg_O3_obs'])
    df_rolling_8hr_avg_ozone = mda8_data
    df_rolling_8hr_avg_ozone.index.name = 'time_local'
    savepath = dirpath + 'F0AM_run_analysis/' + save_filename + '.csv'
    df_rolling_8hr_avg_ozone.to_csv(savepath)
    print('Saved to:' + savepath)

    # Calculate the maximum 8-hour average for each day
    # Resample to daily frequency, compute the daily max of the rolling averages, drop NA values
    daily_max_8hr_avg_ozone = mda8_data[rolling_avg_name].resample('D').max().dropna()
    daily_max_8hr_avg_ozone.index.name = "date"
    daily8hrmax_ozone_savepath = dirpath + 'F0AM_run_analysis/full_campaign_daily_8hr_max_ozone.csv'
    daily_max_8hr_avg_ozone.to_csv(daily8hrmax_ozone_savepath)

    # Select daytime values only where MD8A > 70
    exceedance_name = 'Exceedance_day_' + 'O3'
    exceedance_day_ozone = daily_max_8hr_avg_ozone[daily_max_8hr_avg_ozone >= 70].index.floor('D')

    whole_campaign_dict['Exceedance_days'] = exceedance_day_ozone
    voc_available_days_dict['Exceedance_days']= exceedance_day_ozone[(exceedance_day_ozone >= '2024-07-19') & (exceedance_day_ozone <= '2024-08-12')]

    return whole_campaign_dict['Exceedance_days'], voc_available_days_dict['Exceedance_days']

# def weekdays_and_weekends():

# def smokedays_and_smokefreedays():

whole_campaign_dict['Exceedance_days'], voc_available_days_dict['Exceedance_days'] = mda8_calculate(date_limit = None, save_filename = 'full_campaign_rolling_8hr_avg_ozone')
print(voc_available_days_dict['Exceedance_days'])