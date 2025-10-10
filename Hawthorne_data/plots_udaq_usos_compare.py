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
import matplotlib.dates as mdates

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

pd.set_option('display.max_rows', 200, 'display.min_rows', 200, 'display.max_columns', None)

# NOTE:
# CHECK WHICH POINTS WERE INVALID FROM QUALITY CONTROL (requires some assistance with color scheme in spreadsheet)

######## FUNCTIONS ########

def hawthorne_ozone_fill_usos(udaq_o3_filepath):
    '''
    This function fills the ozone measurements taken during the USOS campaign at Hawthorne from the parked CSL Mobile Lab's 2B Tech UV Analyzer with the ozone measurements taken by UDAQ also at Hawthorne.
    We alter the dates to be from 2024-07-15 00:00:00 to 2024-08-18 23:00:00 [inclusive] even though there are no data for 2024-07-15, so that the data spans the same time frame
    as other species measured. The merged ozone data is saved as a CSV file.

    '''
    all_days_filepath = dirpath + 'CampaignData_and_Merges/R0/CSL_MobileLab_Parked/merged/rev_30min/all_CSL_MobileLab_Parked_rev30minv4.nc'
    all_days_filepath_load = xr.open_dataset(all_days_filepath)
    df_alldays = all_days_filepath_load.to_dataframe()
    df_alldays.reset_index(inplace=True)
    df_alldays.set_index('time_local', inplace=True, drop=False)
    df_alldays['USOS O3'] = df_alldays['O3_ppbv']

    #adds padding of NaNs for first day and last day of campaign, since they need to have the same length of time as the other days in order to plot
    new_start_time = pd.Timestamp('2024-07-15 00:00:00')
    new_end_time = pd.Timestamp('2024-08-18 23:00:00')
    # Create a new datetime index from new_start to the end of existing index with same frequency
    new_index = pd.date_range(start=new_start_time, end=new_end_time, freq='1h')
    # Reindex the dataframe to include new rows
    df_alldays = df_alldays.reindex(new_index)

    udaq_o3_filepath = dirpath + 'Hawthorne_data/hawthorne_udaq_o3_2024.csv'
    df_udaq_o3_load = pd.read_csv(udaq_o3_filepath, index_col = 'Date')
    index_series = pd.Series(df_udaq_o3_load.index)
    index_fixed = index_series.apply(lambda x: x if ':' in x else f'{x} 00:00')
    df_udaq_o3_load.index = pd.to_datetime(index_fixed, format='%m/%d/%Y %H:%M')
    df_udaq_o3_load = df_udaq_o3_load.reindex(new_index)

    df_udaq_o3_load = df_udaq_o3_load.rename(columns = {'O3':'UDAQ O3'})
    df_udaq_o3_load['UDAQ O3'] = df_udaq_o3_load['UDAQ O3']*1000 #ppm to ppb conversion
    df_udaq_o3_load['UDAQ O3'] = df_udaq_o3_load['UDAQ O3'].where(df_udaq_o3_load['UDAQ O3'] <= 120) #remove values over 120 ppb, likely unusable
    print(df_udaq_o3_load['UDAQ O3'])
    
    df_ozone = pd.concat([df_alldays['USOS O3'], df_udaq_o3_load['UDAQ O3']], axis=1)
 
    #Combine the USOS and UDAQ measurements into merged_ozone column, preferentially filling with USOS vals
    df_ozone['merged_ozone'] = df_ozone['USOS O3'].combine_first(df_ozone['UDAQ O3'])

    for instance in range(0,3):
        # Calculate the 8-hour rolling average
        rolling_avg_name = '8hr_rolling_avg_' + df_ozone.columns[instance]
        df_ozone[rolling_avg_name] = df_ozone[df_ozone.columns[instance]].rolling(window=8, min_periods=6).mean()
        #rolling average only works for 8 hours with no gaps; at least 6 hours

        # Calculate the maximum 8-hour average for each day
        # Resample to daily frequency, compute the daily max of the rolling averages, drop NA values
        daily_max_8hr_avg_ozone = df_ozone[rolling_avg_name].resample('D').max().dropna()
        #print('Max 8 hr avg each day:\n, adjuststart_daily_max_8hr_avg_ozone)

        # Map the daily maximum back to the original dataframe
        # Create a new temporary column with the daily max 8-hour average for each timestamp
        mda8_ozone_name = 'MDA8_O3_' + df_ozone.columns[instance]
        df_ozone[mda8_ozone_name] = df_ozone.index.floor('D').map(daily_max_8hr_avg_ozone)

        #print('Daily max 8 hour avg for each timestamp: \n', df_adjuststart_hourly_ozone['MDA8_O3'])

        # # Select daytime values only where MD8A > 70
        # df_ozone_day_exceedance_usos = df_ozone[(df_ozone.index.hour >=7) & (df_ozone.index.hour<=20) & (df_ozone['MDA8_O3_USOS O3']>=70)]
        # df_ozone_day_exceedance_udaq = df_ozone[(df_ozone.index.hour >=7) & (df_ozone.index.hour<=20) & (df_ozone['MDA8_O3_UDAQ O3']>=70)]
        # df_ozone_day_exceedance_filled = df_ozone[(df_ozone.index.hour >=7) & (df_ozone.index.hour<=20) & (df_ozone['MDA8_O3_merged_ozone']>=70)]
        #print('When is MDA8 > 70?', df_ozone_day_exceedance)
        exceedance_name = 'Exceedance_day_' + df_ozone.columns[instance]
        df_ozone[exceedance_name] = df_ozone[mda8_ozone_name] >= 70
    ozone_savepath = dirpath + '/Plotting/USOS_Campaign_analysis/'
    df_ozone.to_csv(ozone_savepath + 'all_ozone_during_usos_with_exceedances.csv', index = True, index_label = 'dt')

def ozone_plot_comparison_udaq_usos_merged(ozone_saved_filepath):
    '''
    This function makes 2 plots and prints 1 DataFrame:
     Plot 1. A time series of the ozone measurements taken during the USOS campaign at Hawthorne from the parked CSL Mobile Lab's 2B Tech UV Analyzer, the ozone measurements taken by UDAQ also at Hawthorne, 
     and the merge of the two sources which is created by the USOS missing measurements filled with UDAQ measurements when necessary.

     Plot 2. A comparison of the MDA8 exceedance days from using the 3 measurement sources in Plot 1 (USOS, UDAQ, the merge of the two measurement sources)

     DataFrame: counts the NaNs in the USOS, UDAQ, and merged data. This helps to determine how many USOS measurements are being filled by UDAQ measurements.
    '''

    df_all_o3_load = pd.read_csv(ozone_saved_filepath, index_col='dt', parse_dates=True)
    #############
    #Time Series Ozone concentration comparison
    fig, ax = plt.subplots(figsize = (30,10))
    ax.plot(df_all_o3_load.index, df_all_o3_load['UDAQ O3'], linestyle = 'solid', color = 'g', marker = '+', label = 'UDAQ')
    ax.plot(df_all_o3_load.index, df_all_o3_load['merged_ozone'], linestyle = 'solid', color='k', marker='o', label='Merged')
    ax.plot(df_all_o3_load.index, df_all_o3_load['USOS O3'], linestyle = 'solid', color='m', marker='x',label='USOS', alpha = 0.7)
    ax.hlines(y=70, xmin = df_all_o3_load.index[0], xmax = df_all_o3_load.index[len(df_all_o3_load.index)-1], linestyle = 'dashed', color = 'r')
    #Mark midnight for every day
    midnight_vals = []
    for midnight_idx in range(0,len(df_all_o3_load.index),24):
        midnight_vals.append(df_all_o3_load.index[midnight_idx])
    for day_pos in midnight_vals:
        ax.axvline(day_pos, color = 'black', linestyle = 'dotted', alpha = 0.7)

    #Set x ticks
    ax.xaxis.set_major_locator(mdates.DayLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
    # Minor ticks: every 6 hours
    ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 6, 12, 18]))
    # Rotate and format tick labels
    ax.tick_params(axis='x', which='major')
    ax.tick_params(axis='x', which='minor', length=4, color='gray')

    plt.yticks(np.arange(10,100,10))

    plt.ylabel('Ozone (ppbv)')
    plt.xlabel('Date')
    plt.margins(x=0)

    #Put label order to USOS first (primary measurement), then UDAQ (filling holes from USOS data), then merged (final combination)
    handles, labels = plt.gca().get_legend_handles_labels()
    order = [2,0,1]
    plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])

    plt.title('Ozone Comparison Between UDAQ and USOS Measurements')

    #############
    # MDA8 Exceedance Comparison
    fig, ax = plt.subplots(3,1, figsize = (30,15))
    threshold = 70

    instance_name = ['USOS', 'UDAQ', 'Merged']
    instance_varname = ['USOS O3', 'UDAQ O3', 'merged_ozone']
    for instance in range(0,3):
        mda8_ozone_name = 'MDA8_O3_' + df_all_o3_load.columns[instance]
        exceedance_name = 'Exceedance_day_' + df_all_o3_load.columns[instance]
        mda8_exceedance = np.where(df_all_o3_load[mda8_ozone_name].values > threshold, df_all_o3_load[mda8_ozone_name].values , np.nan)
        ax[instance].plot(df_all_o3_load.index, df_all_o3_load[mda8_ozone_name], color='b', linewidth = 3)
        ax[instance].plot(df_all_o3_load.index, mda8_exceedance, color='r', linewidth = 3)
        ax[instance].margins(x=0)
        ax[instance].set_ylabel('Ozone (ppbv)')
        ax[instance].set_xlabel('Date')
        ax[instance].hlines(y=70, xmin = df_all_o3_load.index[0], xmax = df_all_o3_load.index[len(df_all_o3_load.index)-1], linestyle = 'dashed', color = 'r')
        ax[instance].set_title('MDA8 Ozone for ' + instance_name[instance])

        ax[instance].xaxis.set_major_locator(mdates.DayLocator())
        ax[instance].xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))

        # Minor ticks: every 6 hours
        ax[instance].xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 6, 12, 18]))

        # Rotate and format tick labels
        ax[instance].tick_params(axis='x', which='major')
        ax[instance].tick_params(axis='x', which='minor', length=4, color='gray')

        ax[instance].set_yticks(np.arange(30,90,10))

        midnight_vals = []
        for midnight_idx in range(0,len(df_all_o3_load.index),24):
            midnight_vals.append(df_all_o3_load.index[midnight_idx])
        for day_pos in midnight_vals:
            ax[instance].axvline(day_pos, color = 'black', linestyle = 'dotted', alpha = 0.7)

    df_normalized = df_all_o3_load.copy()
    df_normalized['date'] = df_all_o3_load.index.normalize()

    nan_counts = df_normalized.isna().groupby(df_normalized['date']).sum()
    display(nan_counts)

### CALL USABLE FUNCTIONS ########

hawthorne_ozone_fill_usos(
    udaq_o3_filepath = dirpath + 'Hawthorne_data/hawthorne_udaq_o3_2024.csv'
)

ozone_plot_comparison_udaq_usos_merged(
    ozone_saved_filepath = dirpath + 'Plotting/USOS_Campaign_analysis/all_ozone_during_usos_with_exceedances.csv'
)


######## SCRAP FUNCTIONS ########
# def hawthorne_usos_comparison_vocs(udaq_voc_filepath):
#     #Read UDAQ VOC measurements file
#     df_data = pd.read_csv(udaq_voc_filepath) 
#     #Get only the data at Hawthorne
#     df_udaq_hawthornedates_usos_only = df_data.loc[(df_data['StationSym'] == 'HW')]
#     #set index to datetimeindex
#     df_udaq_hawthornedates_usos_only.set_index(['dt'], inplace = True) 
#     df_udaq_hawthornedates_usos_only.index = pd.to_datetime(df_udaq_hawthornedates_usos_only.index)

#     #Get only the data during the USOS campaign
#     df_udaq_hawthornedates_usos_only = df_udaq_hawthornedates_usos_only.sort_index().loc['2024-07-14 00:00:00':'2024-08-18 23:00:00']
#     #View Parameter code values as only an int (they print with a .0 at the end initially)
#     #Sort and only include the unique parameter codes
#     df_udaq_hawthornedates_usos_only.Parameter = df_udaq_hawthornedates_usos_only.Parameter.astype(int)
#     # sorted_parameters = sorted(df_hawthorne_usos_only.Parameter.unique())
#     # print(sorted_parameters)
#     #Correspond to:
#     #43243: isoprene
#     #45201: Benzene
#     #45202: Toluene
#     #45220: Styrene
#     parameter_codes_needed = [43243, 45201, 45202]
#     df_hawthorne_species_needed = df_udaq_hawthornedates_usos_only.loc[df_udaq_hawthornedates_usos_only['Parameter'].isin(parameter_codes_needed)]
#     pd.set_option('display.max_rows', 200, 'display.min_rows', 100, 'display.max_columns', None)
#     #display(df_hawthorne_species_needed)
#     df_reshaped_hawthorne_species = df_hawthorne_species_needed.groupby([df_hawthorne_species_needed.index, 'Parameter'])['Sample Value'].first().unstack()
#     print(df_reshaped_hawthorne_species)
#     #df_reshaped_hawthorne_species.columns.name = None

#     # all_days_filepath = dirpath + 'CampaignData_and_Merges/R0/CSL_MobileLab_Parked/merged/rev_30min/all_CSL_MobileLab_Parked_rev30minv4.nc'
#     # all_days_filepath_load = xr.open_dataset(all_days_filepath)
#     # df_alldays = all_days_filepath_load.to_dataframe()
#     # df_alldays.reset_index(inplace=True)
#     # df_alldays.set_index('time_local', inplace=True, drop=False)

#     # #adds padding of NaNs for first day and last day of campaign, since they need to have the same length of time as the other days in order to plot
#     # new_start_time = pd.Timestamp('2024-07-14 00:00:00')
#     # new_end_time = pd.Timestamp('2024-08-18 23:30:00')
#     # # Create a new datetime index from new_start to the end of existing index with same frequency
#     # new_index = pd.date_range(start=new_start_time, end=new_end_time, freq='1h')

#     # # Reindex the dataframe to include new rows
#     # df_alldays = df_alldays.reindex(new_index)
#     # df_usos_species_needed = ['Isoprene_PTR', 'Benzene_PTR', 
#     #                           'Toluene_PTR', 'Styrene_PTR']
#     # df_alldays_reshaped = df_alldays[df_usos_species_needed]

#     # # df_reshaped_hawthorne_species = df_reshaped_hawthorne_species.dropna(subset = [43243, 45201, 45202, 45220])
#     # # df_alldays_reshaped = df_alldays_reshaped.dropna(subset = ['Isoprene_PTR', 'Benzene_PTR', 'Toluene_PTR', 'Styrene_PTR'])
    
#     # df_isoprene = pd.concat([df_reshaped_hawthorne_species[43243], df_alldays_reshaped['Isoprene_PTR']], axis=1)
#     # hr_avg_isoprene_before_drop = df_isoprene.groupby(df_isoprene.index.hour).mean()
#     # display(hr_avg_isoprene_before_drop)

#     # df_isoprene = df_isoprene.dropna(subset = [43243, 'Isoprene_PTR'])
#     # r2_isoprene = r2_score(df_isoprene[43243], df_isoprene['Isoprene_PTR'])
#     # print(r2_isoprene)

#     # #obtain m (slope) and b(intercept) of linear regression line
#     # m_isoprene, b_isoprene = np.polyfit(df_isoprene[43243], df_isoprene['Isoprene_PTR'], 1)

#     # plt.scatter(df_isoprene[43243], df_isoprene['Isoprene_PTR'])
#     # plt.plot(df_isoprene[43243], m_isoprene*df_isoprene[43243] + b_isoprene, color = 'k')
#     # plt.title('Isoprene')
#     # plt.xlabel('UDAQ')
#     # plt.ylabel('USOS PTR')
#     # plt.show()

#     # df_benzene = pd.concat([df_reshaped_hawthorne_species[45201], df_alldays_reshaped['Benzene_PTR']], axis=1)
#     # hr_avg_benzene_before_drop = df_benzene.groupby(df_benzene.index.hour).mean()
#     # display(hr_avg_benzene_before_drop)
#     # df_benzene = df_benzene.dropna(subset = [45201, 'Benzene_PTR'])
#     # r2_benzene = r2_score(df_benzene[45201], df_benzene['Benzene_PTR'])
#     # print(r2_benzene)

#     # m_benzene, b_benzene = np.polyfit(df_benzene[45201], df_benzene['Benzene_PTR'], 1)

#     # plt.scatter(df_benzene[45201], df_benzene['Benzene_PTR'])
#     # plt.plot(df_benzene[45201], m_benzene*df_benzene[45201] + b_benzene, color = 'k')
#     # plt.title('Benzene')
#     # plt.xlabel('UDAQ')
#     # plt.ylabel('USOS PTR')
#     # plt.show()

#     # df_toluene = pd.concat([df_reshaped_hawthorne_species[45202], df_alldays_reshaped['Toluene_PTR']], axis=1)
#     # hr_avg_toluene_before_drop = df_toluene.groupby(df_toluene.index.hour).mean()
#     # display(hr_avg_toluene_before_drop)
#     # df_toluene = df_toluene.dropna(subset = [45202, 'Toluene_PTR'])
#     # r2_toluene = r2_score(df_toluene[45202], df_toluene['Toluene_PTR'])

#     # print(r2_toluene)
#     # m_toluene, b_toluene = np.polyfit(df_toluene[45202], df_toluene['Toluene_PTR'], 1)

#     # plt.scatter(df_toluene[45202], df_toluene['Toluene_PTR'])
#     # plt.plot(df_toluene[45202], m_toluene * df_toluene[45202] + b_toluene, color = 'k')
#     # plt.title('Toluene')
#     # plt.xlabel('UDAQ')
#     # plt.ylabel('USOS PTR')
#     # plt.show()


#     # ############
#     # #display(df_isoprene)
#     # hr_avg_isoprene = df_isoprene.groupby(df_isoprene.index.hour).mean()
#     # display(hr_avg_isoprene)
#     # r2_isoprene_hourly = r2_score(hr_avg_isoprene[43243], hr_avg_isoprene['Isoprene_PTR'])
#     # print(r2_isoprene_hourly)
    
#     # #obtain m (slope) and b(intercept) of linear regression line
#     # m_isoprene, b_isoprene = np.polyfit(hr_avg_isoprene[43243], hr_avg_isoprene['Isoprene_PTR'], 1)

#     # plt.scatter(hr_avg_isoprene[43243], hr_avg_isoprene['Isoprene_PTR'])
#     # plt.plot(hr_avg_isoprene[43243], m_isoprene*hr_avg_isoprene[43243] + b_isoprene, color = 'k')
#     # plt.title('Isoprene')
#     # plt.xlabel('UDAQ')
#     # plt.ylabel('USOS PTR')
#     # plt.show()

# def hawthorne_usos_ozone_comparison(udaq_o3_filepath):
    # all_days_filepath = dirpath + 'CampaignData_and_Merges/R0/CSL_MobileLab_Parked/merged/rev_30min/all_CSL_MobileLab_Parked_rev30minv4.nc'
    # all_days_filepath_load = xr.open_dataset(all_days_filepath)
    # df_alldays = all_days_filepath_load.to_dataframe()
    # df_alldays.reset_index(inplace=True)
    # df_alldays.set_index('time_local', inplace=True, drop=False)

    # #adds padding of NaNs for first day and last day of campaign, since they need to have the same length of time as the other days in order to plot
    # new_start_time = pd.Timestamp('2024-07-16 00:00:00')
    # new_end_time = pd.Timestamp('2024-08-18 23:00:00')
    # # Create a new datetime index from new_start to the end of existing index with same frequency
    # new_index = pd.date_range(start=new_start_time, end=new_end_time, freq='1h')

    # df_udaq_o3_load = pd.read_csv(udaq_o3_filepath)
    # # Create a datetime index: 
    # date_index = pd.date_range(start='2024-07-16', periods=len(df_udaq_o3_load), freq='1h')
    # df_udaq_o3_load.index = date_index

    # # Reindex the dataframe to include new rows
    # df_alldays = df_alldays.reindex(new_index)
    # #print(df_alldays.index)
    # print(df_udaq_o3_load.index)

    # df_udaq_o3_load['O3'] = df_udaq_o3_load['O3'].where(df_udaq_o3_load['O3']*1000 <= 120)
    # print(df_udaq_o3_load['O3'])

    # df_ozone = pd.concat([df_udaq_o3_load['O3']*1000, df_alldays['O3_ppbv']], axis=1)
    # df_ozone = df_ozone.dropna(subset = ['O3', 'O3_ppbv'])
    # m_ozone, b_ozone = np.polyfit(df_ozone['O3'], df_ozone['O3_ppbv'], 1)
    # print(m_ozone)
    # print(b_ozone)
    # display(df_ozone)
    # hr_avg_o3_before_drop = df_ozone.groupby(df_ozone.index.hour).mean()
    # display(hr_avg_o3_before_drop)

    # plt.scatter(df_ozone['O3'], df_ozone['O3_ppbv'])
    # plt.plot(df_ozone['O3'], m_ozone * df_ozone['O3'] + b_ozone, color = 'k')
    # plt.title('Ozone')
    # plt.xlabel('UDAQ Obs')
    # plt.ylabel('USOS Obs')
    # plt.show()

    # ratio_ozone = df_ozone['O3']/df_ozone['O3_ppbv']
    # print(ratio_ozone)
    # print(ratio_ozone.mean())
    # print(ratio_ozone.groupby(ratio_ozone.index.hour).mean())

# hawthorne_usos_comparison_vocs(
#     udaq_voc_filepath = dirpath + 'Hawthorne_data/Verbose.csv'
# )

# hawthorne_usos_ozone_comparison(
#     udaq_o3_filepath = dirpath + 'Hawthorne_data/hawthorne_udaq_o3_2024.csv'
# )