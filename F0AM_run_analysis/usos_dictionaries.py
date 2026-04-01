# For saving dictionaries of different dates for USOS Campaign

import pandas as pd
import numpy as np

dirpath = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/'
interp_vals_filepath = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/CampaignData_and_Merges/R0/CSL_MobileLab_Parked/F0AM_filled/interpolated_values_15min.csv'
df_interp_vals = pd.read_csv(interp_vals_filepath, index_col = 'time_local', parse_dates = True)

# whole_campaign_dict = {}
# voc_available_days_dict = {}

# def voc_available_days_slicing(whole_campaign_input, varname_for_dict_saving):
#     voc_available_days_dict[varname_for_dict_saving]= whole_campaign_input[(whole_campaign_input >= '2024-07-19') & (whole_campaign_input <= '2024-08-12')]
#     return voc_available_days_dict

def mda8_calculate(o3_8hr_rolling_avg_save_filename):
    df_o3_obs = df_interp_vals['O3_ppbv']
    #resample from 15 minute intervals to hourly
    mda8_data = df_o3_obs.resample('1h').mean()

    
    # Calculate the 8-hour rolling average
    rolling_avg_name = '8hr_rolling_avg_O3_obs'
    mda8_data[rolling_avg_name] = mda8_data.rolling(window=8, min_periods=6).mean()
    #print("mda8_data['8hr_rolling_avg_O3_obs']: ", mda8_data['8hr_rolling_avg_O3_obs'])
    df_rolling_8hr_avg_ozone = mda8_data
    df_rolling_8hr_avg_ozone.index.name = 'time_local'
    savepath = dirpath + 'F0AM_run_analysis/' + o3_8hr_rolling_avg_save_filename + '.csv'
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
    
    return exceedance_day_ozone
def weekdays_weekends_dayofweek():
    df_weekdays = df_interp_vals.index[df_interp_vals.index.dayofweek < 5]
    df_weekends = df_interp_vals.index[df_interp_vals.index.dayofweek >= 5]
    df_mondays = df_interp_vals.index[df_interp_vals.index.dayofweek == 0]
    df_tuesdays = df_interp_vals.index[df_interp_vals.index.dayofweek == 1]
    df_wednesdays = df_interp_vals.index[df_interp_vals.index.dayofweek == 2]
    df_thursdays = df_interp_vals.index[df_interp_vals.index.dayofweek == 3]
    df_fridays = df_interp_vals.index[df_interp_vals.index.dayofweek == 4]
    df_saturdays = df_interp_vals.index[df_interp_vals.index.dayofweek == 5]
    df_sundays = df_interp_vals.index[df_interp_vals.index.dayofweek == 6]
    return df_weekdays, df_weekends, df_mondays, df_tuesdays, df_wednesdays, df_thursdays, df_fridays, df_saturdays, df_sundays
def smokedays_and_smokefreedays():
    #only modeled days that were:
    # warm(not rain/cold days: 2024-08-09, 2024-08-12 to 16)
    # days when drives weren't happening during peak photochemical hours (drives during: 2024-07-18, 07-22, 07-23, 07-26, 07-28, 07-30, 08-01, 08-03)

    #compared NOx on weekend vs weekday

    #smoke influenced days: elevated acetonitrile, CO, aerosols
    #looks like should include 2024-07-31, 2024-08-02, and 2024-08-08? Nell modeled 3 smoke days (2024-07-31, 2024-08-02, 2024-08-08)
    #2024-07-19?, 2024-07-24, 2024-07-25,  2024-07-29, 2024-08-05, 2024-08-06, 2024-08-07? Nell modeled 7 non-smoke days
    #saw elevated CO, total VOCs, and NOx on smoke days
    
    #Smoke days -> NOx not elevated -> NOx lost to reservoir species [not sure if this means that when NOx isn't being elevated on smoke days, then NOx is being lost to res species?
    # NOx is elevated on smoke days so my notes aren't very clear]
    #PAN also elevated (but maybe look at all NOx reservoir species? They should be higher when NOx is lower?)


#Get function outputs as variables to put into dictionary
df_weekdays, df_weekends, df_mondays, df_tuesdays, df_wednesdays, df_thursdays, df_fridays, df_saturdays, df_sundays = weekdays_weekends_dayofweek()

#Make dictionary for whole campaign
whole_campaign_dict = {
    'Exceedance_days': mda8_calculate(save_filename = 'full_campaign_rolling_8hr_avg_ozone'),
    'Weekdays': df_weekdays,
    'Weekends': df_weekends,
    'Mondays': df_mondays, 
    'Tuesdays': df_tuesdays, 
    'Wednesdays': df_wednesdays, 
    'Thursdays': df_thursdays, 
    'Fridays': df_fridays, 
    'Saturdays': df_saturdays, 
    'Sundays': df_sundays
    }


# voc_available_days_dict['Exceedance_days'] = voc_available_days_slicing(whole_campaign_input = whole_campaign_dict['Exceedance_days'], 
#                                                     varname_for_dict_saving = 'Exceedance_days')

