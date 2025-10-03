import os 
import sys 
import re 
import yaml
import inspect 
import numpy as np 
import pandas as pd
import xarray as xr
from collections import defaultdict
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt

from sklearn.metrics import r2_score

#get path for USOS_shared directory for correct laptop
#input parameter should be 'CHPC', 'Mac', or 'Windows'
current_dir = os.path.dirname(os.path.abspath(__file__))
global_scripts_path = os.path.abspath(os.path.join(current_dir, "..", "global_scripts"))
# sys.path.insert(0,global_scripts_path)
# from dirpath import filepath_source
# dirpath = filepath_source('CHPC')
dirpath = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/'

#define path for Hawthorne data directory
savepath= dirpath + 'Hawthorne_data/'

def hawthorne_usos_comparison_vocs(udaq_voc_filepath):
    #Read UDAQ VOC measurements file
    df_data = pd.read_csv(udaq_voc_filepath) 
    #Get only the data at Hawthorne
    df_udaq_hawthornedates_usos_only = df_data.loc[(df_data['StationSym'] == 'HW')]
    #set index to datetimeindex
    df_udaq_hawthornedates_usos_only.set_index(['dt'], inplace = True) 
    df_udaq_hawthornedates_usos_only.index = pd.to_datetime(df_udaq_hawthornedates_usos_only.index)

    #Get only the data during the USOS campaign
    df_udaq_hawthornedates_usos_only = df_udaq_hawthornedates_usos_only.sort_index().loc['2024-07-14 00:00:00':'2024-08-18 23:00:00']
    #View Parameter code values as only an int (they print with a .0 at the end initially)
    #Sort and only include the unique parameter codes
    df_udaq_hawthornedates_usos_only.Parameter = df_udaq_hawthornedates_usos_only.Parameter.astype(int)
    # sorted_parameters = sorted(df_hawthorne_usos_only.Parameter.unique())
    # print(sorted_parameters)
    #Correspond to:
    #43243: isoprene
    #45201: Benzene
    #45202: Toluene
    #45220: Styrene
    parameter_codes_needed = [43243, 45201, 45202]
    df_hawthorne_species_needed = df_udaq_hawthornedates_usos_only.loc[df_udaq_hawthornedates_usos_only['Parameter'].isin(parameter_codes_needed)]
    pd.set_option('display.max_rows', 200, 'display.min_rows', 100, 'display.max_columns', None)
    #display(df_hawthorne_species_needed)
    df_reshaped_hawthorne_species = df_hawthorne_species_needed.groupby([df_hawthorne_species_needed.index, 'Parameter'])['Sample Value'].first().unstack()
    print(df_reshaped_hawthorne_species)
    #df_reshaped_hawthorne_species.columns.name = None

    # all_days_filepath = dirpath + 'CampaignData_and_Merges/R0/CSL_MobileLab_Parked/merged/rev_30min/all_CSL_MobileLab_Parked_rev30minv4.nc'
    # all_days_filepath_load = xr.open_dataset(all_days_filepath)
    # df_alldays = all_days_filepath_load.to_dataframe()
    # df_alldays.reset_index(inplace=True)
    # df_alldays.set_index('time_local', inplace=True, drop=False)

    # #adds padding of NaNs for first day and last day of campaign, since they need to have the same length of time as the other days in order to plot
    # new_start_time = pd.Timestamp('2024-07-14 00:00:00')
    # new_end_time = pd.Timestamp('2024-08-18 23:30:00')
    # # Create a new datetime index from new_start to the end of existing index with same frequency
    # new_index = pd.date_range(start=new_start_time, end=new_end_time, freq='1h')

    # # Reindex the dataframe to include new rows
    # df_alldays = df_alldays.reindex(new_index)
    # df_usos_species_needed = ['Isoprene_PTR', 'Benzene_PTR', 
    #                           'Toluene_PTR', 'Styrene_PTR']
    # df_alldays_reshaped = df_alldays[df_usos_species_needed]

    # # df_reshaped_hawthorne_species = df_reshaped_hawthorne_species.dropna(subset = [43243, 45201, 45202, 45220])
    # # df_alldays_reshaped = df_alldays_reshaped.dropna(subset = ['Isoprene_PTR', 'Benzene_PTR', 'Toluene_PTR', 'Styrene_PTR'])
    
    # df_isoprene = pd.concat([df_reshaped_hawthorne_species[43243], df_alldays_reshaped['Isoprene_PTR']], axis=1)
    # hr_avg_isoprene_before_drop = df_isoprene.groupby(df_isoprene.index.hour).mean()
    # display(hr_avg_isoprene_before_drop)

    # df_isoprene = df_isoprene.dropna(subset = [43243, 'Isoprene_PTR'])
    # r2_isoprene = r2_score(df_isoprene[43243], df_isoprene['Isoprene_PTR'])
    # print(r2_isoprene)

    # #obtain m (slope) and b(intercept) of linear regression line
    # m_isoprene, b_isoprene = np.polyfit(df_isoprene[43243], df_isoprene['Isoprene_PTR'], 1)

    # plt.scatter(df_isoprene[43243], df_isoprene['Isoprene_PTR'])
    # plt.plot(df_isoprene[43243], m_isoprene*df_isoprene[43243] + b_isoprene, color = 'k')
    # plt.title('Isoprene')
    # plt.xlabel('UDAQ')
    # plt.ylabel('USOS PTR')
    # plt.show()

    # df_benzene = pd.concat([df_reshaped_hawthorne_species[45201], df_alldays_reshaped['Benzene_PTR']], axis=1)
    # hr_avg_benzene_before_drop = df_benzene.groupby(df_benzene.index.hour).mean()
    # display(hr_avg_benzene_before_drop)
    # df_benzene = df_benzene.dropna(subset = [45201, 'Benzene_PTR'])
    # r2_benzene = r2_score(df_benzene[45201], df_benzene['Benzene_PTR'])
    # print(r2_benzene)

    # m_benzene, b_benzene = np.polyfit(df_benzene[45201], df_benzene['Benzene_PTR'], 1)

    # plt.scatter(df_benzene[45201], df_benzene['Benzene_PTR'])
    # plt.plot(df_benzene[45201], m_benzene*df_benzene[45201] + b_benzene, color = 'k')
    # plt.title('Benzene')
    # plt.xlabel('UDAQ')
    # plt.ylabel('USOS PTR')
    # plt.show()

    # df_toluene = pd.concat([df_reshaped_hawthorne_species[45202], df_alldays_reshaped['Toluene_PTR']], axis=1)
    # hr_avg_toluene_before_drop = df_toluene.groupby(df_toluene.index.hour).mean()
    # display(hr_avg_toluene_before_drop)
    # df_toluene = df_toluene.dropna(subset = [45202, 'Toluene_PTR'])
    # r2_toluene = r2_score(df_toluene[45202], df_toluene['Toluene_PTR'])

    # print(r2_toluene)
    # m_toluene, b_toluene = np.polyfit(df_toluene[45202], df_toluene['Toluene_PTR'], 1)

    # plt.scatter(df_toluene[45202], df_toluene['Toluene_PTR'])
    # plt.plot(df_toluene[45202], m_toluene * df_toluene[45202] + b_toluene, color = 'k')
    # plt.title('Toluene')
    # plt.xlabel('UDAQ')
    # plt.ylabel('USOS PTR')
    # plt.show()


    # ############
    # #display(df_isoprene)
    # hr_avg_isoprene = df_isoprene.groupby(df_isoprene.index.hour).mean()
    # display(hr_avg_isoprene)
    # r2_isoprene_hourly = r2_score(hr_avg_isoprene[43243], hr_avg_isoprene['Isoprene_PTR'])
    # print(r2_isoprene_hourly)
    
    # #obtain m (slope) and b(intercept) of linear regression line
    # m_isoprene, b_isoprene = np.polyfit(hr_avg_isoprene[43243], hr_avg_isoprene['Isoprene_PTR'], 1)

    # plt.scatter(hr_avg_isoprene[43243], hr_avg_isoprene['Isoprene_PTR'])
    # plt.plot(hr_avg_isoprene[43243], m_isoprene*hr_avg_isoprene[43243] + b_isoprene, color = 'k')
    # plt.title('Isoprene')
    # plt.xlabel('UDAQ')
    # plt.ylabel('USOS PTR')
    # plt.show()

def hawthorne_usos_ozone_comparison(udaq_o3_filepath):
    all_days_filepath = dirpath + 'CampaignData_and_Merges/R0/CSL_MobileLab_Parked/merged/rev_30min/all_CSL_MobileLab_Parked_rev30minv4.nc'
    all_days_filepath_load = xr.open_dataset(all_days_filepath)
    df_alldays = all_days_filepath_load.to_dataframe()
    df_alldays.reset_index(inplace=True)
    df_alldays.set_index('time_local', inplace=True, drop=False)

    #adds padding of NaNs for first day and last day of campaign, since they need to have the same length of time as the other days in order to plot
    new_start_time = pd.Timestamp('2024-07-16 00:00:00')
    new_end_time = pd.Timestamp('2024-08-18 23:00:00')
    # Create a new datetime index from new_start to the end of existing index with same frequency
    new_index = pd.date_range(start=new_start_time, end=new_end_time, freq='1h')

    df_udaq_o3_load = pd.read_csv(udaq_o3_filepath)
    # Create a datetime index: 
    date_index = pd.date_range(start='2024-07-16', periods=len(df_udaq_o3_load), freq='1h')
    df_udaq_o3_load.index = date_index

    # Reindex the dataframe to include new rows
    df_alldays = df_alldays.reindex(new_index)
    #print(df_alldays.index)
    print(df_udaq_o3_load.index)

    df_udaq_o3_load['O3'] = df_udaq_o3_load['O3'].where(df_udaq_o3_load['O3']*1000 <= 120)
    print(df_udaq_o3_load['O3'])

    df_ozone = pd.concat([df_udaq_o3_load['O3']*1000, df_alldays['O3_ppbv']], axis=1)
    df_ozone = df_ozone.dropna(subset = ['O3', 'O3_ppbv'])
    m_ozone, b_ozone = np.polyfit(df_ozone['O3'], df_ozone['O3_ppbv'], 1)
    print(m_ozone)
    print(b_ozone)
    display(df_ozone)
    hr_avg_o3_before_drop = df_ozone.groupby(df_ozone.index.hour).mean()
    display(hr_avg_o3_before_drop)

    plt.scatter(df_ozone['O3'], df_ozone['O3_ppbv'])
    plt.plot(df_ozone['O3'], m_ozone * df_ozone['O3'] + b_ozone, color = 'k')
    plt.title('Ozone')
    plt.xlabel('UDAQ Obs')
    plt.ylabel('USOS Obs')
    plt.show()

    ratio_ozone = df_ozone['O3']/df_ozone['O3_ppbv']
    print(ratio_ozone)
    print(ratio_ozone.mean())
    print(ratio_ozone.groupby(ratio_ozone.index.hour).mean())

def hawthorne_ozone_fill_usos(udaq_o3_filepath):
    all_days_filepath = dirpath + 'CampaignData_and_Merges/R0/CSL_MobileLab_Parked/merged/rev_30min/all_CSL_MobileLab_Parked_rev30minv4.nc'
    all_days_filepath_load = xr.open_dataset(all_days_filepath)
    df_alldays = all_days_filepath_load.to_dataframe()
    df_alldays.reset_index(inplace=True)
    df_alldays.set_index('time_local', inplace=True, drop=False)
    df_alldays['USOS O3'] = df_alldays['O3_ppbv']

    #adds padding of NaNs for first day and last day of campaign, since they need to have the same length of time as the other days in order to plot
    new_start_time = pd.Timestamp('2024-07-16 00:00:00')
    new_end_time = pd.Timestamp('2024-08-18 23:00:00')
    # Create a new datetime index from new_start to the end of existing index with same frequency
    new_index = pd.date_range(start=new_start_time, end=new_end_time, freq='1h')

    df_udaq_o3_load = pd.read_csv(udaq_o3_filepath)
    # Create a datetime index: 
    date_index = pd.date_range(start='2024-07-16', periods=len(df_udaq_o3_load), freq='1h')
    df_udaq_o3_load.index = date_index
    print(date_index)

    # Reindex the dataframe to include new rows
    df_alldays = df_alldays.reindex(new_index)

    df_udaq_o3_load['UDAQ O3'] = df_udaq_o3_load['O3']*1000
    df_udaq_o3_load['UDAQ O3'] = df_udaq_o3_load['UDAQ O3'].where(df_udaq_o3_load['UDAQ O3'] <= 120)
    
    df_ozone = pd.concat([df_udaq_o3_load['UDAQ O3'], df_alldays['USOS O3']], axis=1)
    df_ozone.index = date_index
    print(df_ozone.index)

    df_ozone['merged_ozone'] = df_ozone['USOS O3'].combine_first(df_ozone['UDAQ O3'])

    pd.set_option('display.max_rows', 500, 'display.min_rows', 100)
    display(df_ozone)
    
    # Step 1: Calculate the 8-hour rolling average
    df_ozone['8hr_rolling_avg'] = df_ozone['merged_ozone'].rolling(window=8, min_periods=6).mean()
    #rolling average only works for 8 hours with no gaps; at least 6 hours

    # Step 2: Calculate the maximum 8-hour average for each day
    # Resample to daily frequency, compute the daily max of the rolling averages, drop NA values
    daily_max_8hr_avg_ozone = df_ozone['8hr_rolling_avg'].resample('D').max().dropna()
    #print('Max 8 hr avg each day:\n, adjuststart_daily_max_8hr_avg_ozone)

    # Step 3: Map the daily maximum back to the original dataframe
    # Create a new temporary column with the daily max 8-hour average for each timestamp
    df_ozone['MDA8_O3'] = df_ozone.index.floor('D').map(daily_max_8hr_avg_ozone)
    display(df_ozone['MDA8_O3'])
    #print('Daily max 8 hour avg for each timestamp: \n', df_adjuststart_hourly_ozone['MDA8_O3'])

    # Select daytime values only where MD8A > 70
    df_ozone_day_exceedance = df_ozone[(df_ozone.index.hour >=7) & (df_ozone.index.hour<=20) & (df_ozone['MDA8_O3']>=70)]
    #print('When is MDA8 > 70?', df_ozone_day_exceedance)
    df_ozone['Exceedance_day'] = df_ozone['MDA8_O3'] >= 70

    ozone_savepath = dirpath + '/Plotting/USOS_Campaign_analysis/'
    df_ozone.to_csv(ozone_savepath + 'all_ozone_during_usos_with_exceedances.csv')

    threshold = 70
    mda8_exceedance = np.where(df_ozone['MDA8_O3'].values > threshold, df_ozone['MDA8_O3'].values , np.nan)

    # Create boolean masks
    # below_threshold = df_ozone['MDA8_O3'] <= threshold
    # above_threshold = df_ozone['MDA8_O3'] > threshold
    fig, ax = plt.subplots(figsize = (30,10))
    plt.plot(df_ozone.index, df_ozone['MDA8_O3'], color='b', linewidth = 3)
    plt.plot(df_ozone.index, mda8_exceedance, color='r', linewidth = 3)
    plt.margins(x=0)

    #Set x-axis intervals
    ozone_ticks = []
    for ozone_tick_idx in range(6,len(df_ozone.index),6):
        ozone_ticks.append(df_ozone.index[ozone_tick_idx])
    print(ozone_ticks)
    plt.xticks(ticks=ozone_ticks, labels=df_ozone.index.strftime('%m/%d %H%M')[6::6], rotation = 90)

    midnight_vals = []
    for midnight_idx in range(24,len(df_ozone.index),24):
        midnight_vals.append(df_ozone.index[midnight_idx])
    for day_pos in midnight_vals:
        plt.axvline(day_pos, color = 'black')

    plt.plot(df_ozone.index, df_ozone['merged_ozone'], color='k')

    plt.hlines(y=70, xmin = df_ozone['MDA8_O3'].index[0], xmax = df_ozone['MDA8_O3'].index[len(df_ozone['MDA8_O3'])-1], color = 'r')
    plt.xlabel('Time (Date and Hour)')
    plt.ylabel('Ozone (ppbv)')

    plt.title('MDA8 Ozone during USOS')
    
hawthorne_usos_comparison_vocs(
    udaq_voc_filepath = dirpath + 'Hawthorne_data/Verbose.csv'
)

# hawthorne_usos_ozone_comparison(
#     udaq_o3_filepath = dirpath + 'Hawthorne_data/hawthorne_udaq_o3_2024.csv'
# )

# hawthorne_ozone_fill_usos(
#     udaq_o3_filepath = dirpath + 'Hawthorne_data/hawthorne_udaq_o3_2024.csv'
# )