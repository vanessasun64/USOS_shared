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
from matplotlib.colors import ListedColormap

from scipy.io import savemat
from collections import OrderedDict
dirpath = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/'
mappings_filepath = dirpath + 'Hawthorne_data/mappings/manually_edited/UDAQ_Hawthorne_CRACMM_GEOSCHEM_CB6r5h_mapped_updated_11172025.csv'
merged_data_dir = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/CampaignData_and_Merges/R0/CSL_MobileLab_Parked/merged/'
merged_data_dir_15min = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/CampaignData_and_Merges/R0/CSL_MobileLab_Parked/merged/rev_15min/'
hawthorne_data_dir = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/Hawthorne_data/data/'

#region: NOT WORKING RIGHT NOW, wrong times
def resave_hourly_ozone(filename):
    udaq_o3_filepath = dirpath + 'Hawthorne_data/data/hawthorne_udaq_o3_2024.csv'
    df_udaq_o3_load = pd.read_csv(udaq_o3_filepath, index_col = 'Date')
    index_series = pd.Series(df_udaq_o3_load.index)
    index_fixed = index_series.apply(lambda x: x if ':' in x else f'{x} 00:00')
    df_udaq_o3_load.index = pd.to_datetime(index_fixed, format='%m/%d/%Y %H:%M')
    new_start_time = pd.Timestamp('2024-07-15 00:00:00')
    new_end_time = pd.Timestamp('2024-08-18 23:00:00')
    new_index = pd.date_range(start=new_start_time, end=new_end_time, freq='1h')

    #Shift all values by 1 hour to account for the fact that all values are in MST rather than MDT
    df_udaq_o3_load = df_udaq_o3_load.shift(periods=1)
    df_udaq_o3_load = df_udaq_o3_load.reindex(new_index)
    df_udaq_o3_load = df_udaq_o3_load.rename(columns = {'O3':'UDAQ O3'})
    df_udaq_o3_load['UDAQ O3'] = df_udaq_o3_load['UDAQ O3']*1000 #ppm to ppb conversion
    df_udaq_o3_load['UDAQ O3'] = df_udaq_o3_load['UDAQ O3'].where(df_udaq_o3_load['UDAQ O3'] <= 120) #remove values over 120 ppb, likely unusable
    df_udaq_o3_load.index.name = 'time_local'
    savepath = dirpath + 'Hawthorne_data/data/' + filename + '.csv'
    df_udaq_o3_load.to_csv(savepath)
    print('Saved to:' + savepath)
def resave_one_other_species(parameter_codes, species_names, filename):
    #Quality control file; only has Jul 01 to Aug 01
    udaq_qc_voc_filepath = dirpath + 'Hawthorne_data/data/from_Bart/udaq_QC_07012024_08012024.csv'

    #Read UDAQ QC VOC measurements file
    df_data_qc = pd.read_csv(udaq_qc_voc_filepath) 
    #Get only the data at Hawthorne
    df_udaq_hawthornedates_usos_only_qc = df_data_qc.loc[(df_data_qc['StationSym'] == 'HW')]
    #set index to datetimeindex
    df_udaq_hawthornedates_usos_only_qc.set_index(['dt'], inplace = True) 
    df_udaq_hawthornedates_usos_only_qc.index = pd.to_datetime(df_udaq_hawthornedates_usos_only_qc.index)

    # # #Get only the data during the USOS campaign
    df_udaq_hawthornedates_usos_only_qc = df_udaq_hawthornedates_usos_only_qc.sort_index().loc['2024-07-15 00:00:00':'2024-08-18 23:00:00']
    # #Change dataframe to the parameters for the species we need
    df_udaq_species_needed_qc = df_udaq_hawthornedates_usos_only_qc.loc[df_udaq_hawthornedates_usos_only_qc['Parameter'].isin(parameter_codes)]
    #print(df_udaq_species_needed_qc)

    #Change the sorting so that we only get the sample value and dateindex, take out unnecessary columns
    df_reshaped_udaq_species_qc = df_udaq_species_needed_qc.groupby([df_udaq_species_needed_qc.index, 'Parameter'])['Sample Value'].first().unstack()
    df_reshaped_udaq_species_qc = df_reshaped_udaq_species_qc.rename(columns=dict(zip(parameter_codes,species_names)))
    
    #Shift all values by 1 hour to account for the fact that all values are in MST rather than MDT
    df_reshaped_udaq_species_qc = df_reshaped_udaq_species_qc.shift(periods=1)
    df_reshaped_udaq_species_qc.index = df_reshaped_udaq_species_qc.index.tz_localize('America/Denver')

    new_start_time = pd.Timestamp('2024-07-15 00:00:00')
    new_end_time = pd.Timestamp('2024-08-18 23:00:00')
    new_index = pd.date_range(start=new_start_time, end=new_end_time, freq='1h', tz = 'America/Denver')
    df_reshaped_udaq_species_qc = df_reshaped_udaq_species_qc.reindex(new_index)
    df_reshaped_udaq_species_qc.index.name = 'time_local'

    savepath = dirpath + 'Hawthorne_data/data/' + filename + '.csv'
    df_reshaped_udaq_species_qc.to_csv(savepath)
    print('Saved to:' + savepath)
def resave_all_voc_species_hourly(filename):
    udaq_voc_filepath = dirpath + 'Hawthorne_data/data/Verbose.csv'

    #Read UDAQ VOC measurements file
    df_data = pd.read_csv(udaq_voc_filepath) 
    #Get only the data at Hawthorne
    df_udaq_hawthornedates_usos_only = df_data.loc[(df_data['StationSym'] == 'HW')]
    #set index to datetimeindex
    df_udaq_hawthornedates_usos_only.set_index(['dt'], inplace = True) 
    df_udaq_hawthornedates_usos_only.index = pd.to_datetime(df_udaq_hawthornedates_usos_only.index)

    #Get only the data during the USOS campaign
    df_udaq_hawthornedates_usos_only = df_udaq_hawthornedates_usos_only.sort_index().loc['2024-07-15 00:00:00':'2024-08-18 23:00:00']

    #View Parameter code values as only an int (they print with a .0 at the end initially)
    df_udaq_hawthornedates_usos_only.Parameter = df_udaq_hawthornedates_usos_only.Parameter.astype(int)
    sorted_parameters = sorted(df_udaq_hawthornedates_usos_only.Parameter.unique())
                                
    #Pivot the data
    df_pivot = df_udaq_hawthornedates_usos_only.pivot_table(index=df_udaq_hawthornedates_usos_only.index, columns='Parameter', values='Sample Value', dropna=False)
    #Shift all values by 1 hour to account for the fact that all values are in MST rather than MDT
    df_pivot = df_pivot.shift(periods=1)

    #Resample to hourly to ensure a continuous time index
    new_start_time = pd.Timestamp('2024-07-15 00:00:00')
    new_end_time = pd.Timestamp('2024-08-18 23:00:00')
    #Create a new datetime index from new_start to the end of existing index with same frequency
    new_index = pd.date_range(start=new_start_time, end=new_end_time, freq='1h')
    #Reindex the dataframe to include new rows
    df_pivot = df_pivot.reindex(new_index)
    df_pivot_parameter_sort = df_pivot[sorted(df_pivot.columns)]
    df_pivot_parameter_sort = df_pivot_parameter_sort.drop(43102, axis=1)

    #Get carbon numbers by using the mappings spreadsheet to get the formula and dividing by number of carbons
    df_mapping_parameters = pd.read_csv(mappings_filepath)
    df_mapping_parameters = df_mapping_parameters.drop([0])
    carbon_number = df_mapping_parameters['formula'].str.extract(r'C(\d{1,2})')[0]
    carbon_number = carbon_number.fillna(1)
    carbon_number = carbon_number.astype(int)
    df_udaq_ppbv = df_pivot_parameter_sort.divide(carbon_number.values, axis=1)

    # Create mapping of Parameter Code to UDAQ variable name
    param_to_var = df_mapping_parameters.set_index('Parameter Code')['UDAQ_Variable'].to_dict()

    # # Rename columns in df_data
    df_udaq_vals = df_udaq_ppbv.rename(columns=param_to_var)
    df_udaq_vals.columns.name = None

    df_udaq_vals.index.name = 'time_local'
    savepath = dirpath + 'Hawthorne_data/data/' + filename + '.csv'
    df_udaq_vals.to_csv(savepath)
    print('Saved to:' + savepath)
def compare_qc_and_verbose_vocs():
    #To confirm that the QC values and Verbose.csv values are the same
    udaq_voc_filepath = dirpath + 'Hawthorne_data/data/Verbose.csv'
    #Read UDAQ VOC measurements file
    df_data = pd.read_csv(udaq_voc_filepath) 
    #Get only the data at Hawthorne
    df_udaq_hawthornedates_usos_only = df_data.loc[(df_data['StationSym'] == 'HW')]
    #set index to datetimeindex
    df_udaq_hawthornedates_usos_only.set_index(['dt'], inplace = True) 
    df_udaq_hawthornedates_usos_only.index = pd.to_datetime(df_udaq_hawthornedates_usos_only.index)

    #Get only the data during the USOS campaign
    df_udaq_hawthornedates_usos_only = df_udaq_hawthornedates_usos_only.sort_index().loc['2024-07-15 00:00:00':'2024-08-01 23:00:00']

    #View Parameter code values as only an int (they print with a .0 at the end initially)
    df_udaq_hawthornedates_usos_only.Parameter = df_udaq_hawthornedates_usos_only.Parameter.astype(int)
    sorted_parameters = sorted(df_udaq_hawthornedates_usos_only.Parameter.unique())
                                
    #Pivot the data
    df_pivot = df_udaq_hawthornedates_usos_only.pivot_table(index=df_udaq_hawthornedates_usos_only.index, columns='Parameter', values='Sample Value', dropna=False)
    #Shift all values by 1 hour to account for the fact that all values are in MST rather than MDT
    df_pivot = df_pivot.shift(periods=1)

    #Resample to hourly to ensure a continuous time index
    new_start_time = pd.Timestamp('2024-07-15 00:00:00')
    new_end_time = pd.Timestamp('2024-08-01 23:00:00')
    #Create a new datetime index from new_start to the end of existing index with same frequency
    new_index = pd.date_range(start=new_start_time, end=new_end_time, freq='1h')
    #Reindex the dataframe to include new rows
    df_pivot = df_pivot.reindex(new_index)
    df_pivot_parameter_sort = df_pivot[sorted(df_pivot.columns)]
    df_pivot_parameter_sort = df_pivot_parameter_sort.drop(43102, axis=1)

    #Get carbon numbers by using the mappings spreadsheet to get the formula and dividing by number of carbons
    df_mapping_parameters = pd.read_csv(mappings_filepath)
    df_mapping_parameters = df_mapping_parameters.drop([0])
    carbon_number = df_mapping_parameters['formula'].str.extract(r'C(\d{1,2})')[0]
    carbon_number = carbon_number.fillna(1)
    carbon_number = carbon_number.astype(int)
    df_udaq_ppbv = df_pivot_parameter_sort.divide(carbon_number.values, axis=1)

    # Create mapping of Parameter Code to UDAQ variable name
    param_to_var = df_mapping_parameters.set_index('Parameter Code')['UDAQ_Variable'].to_dict()

    # # Rename columns in df_data
    df_udaq_vals = df_udaq_ppbv.rename(columns=param_to_var)
    df_udaq_vals.columns.name = None

    #Now import the UDAQ QC data
    qc_udaq_voc_filepath = dirpath + 'Hawthorne_data/data/udaq_QC_07012024_08012024.csv'
    #Read UDAQ VOC measurements file
    qc_df_data = pd.read_csv(qc_udaq_voc_filepath) 
    #Get only the data at Hawthorne
    qc_df_udaq_hawthornedates_usos_only = qc_df_data.loc[(qc_df_data['StationSym'] == 'HW')]
    #set index to datetimeindex
    qc_df_udaq_hawthornedates_usos_only.set_index(['dt'], inplace = True) 
    qc_df_udaq_hawthornedates_usos_only.index = pd.to_datetime(qc_df_udaq_hawthornedates_usos_only.index)

    #Get only the data during the USOS campaign
    qc_df_udaq_hawthornedates_usos_only = qc_df_udaq_hawthornedates_usos_only.sort_index().loc['2024-07-15 00:00:00':'2024-08-01 23:00:00']

    #View Parameter code values as only an int (they print with a .0 at the end initially)
    qc_df_udaq_hawthornedates_usos_only.Parameter = qc_df_udaq_hawthornedates_usos_only.Parameter.astype(int)
    qc_sorted_parameters = sorted(qc_df_udaq_hawthornedates_usos_only.Parameter.unique())
    #The QC file has some extra parameters such as carbon monoxide, wind speed, etc. We remove all these non-VOC species that we don't need
    exclude_params = [42101, 42401, 42600, 42601, 42602, 42603, 42612, 61101, 61102, 61301, 62101, 62201, 63301, 64101, 81102, 88101, 88313]
    qc_df_udaq_hawthornedates_usos_only = qc_df_udaq_hawthornedates_usos_only[~qc_df_udaq_hawthornedates_usos_only['Parameter'].isin(exclude_params)]

    # #Pivot the data
    qc_df_pivot = qc_df_udaq_hawthornedates_usos_only.pivot_table(index=qc_df_udaq_hawthornedates_usos_only.index, columns='Parameter', values='Sample Value', dropna=False)
    #Shift all values by 1 hour to account for the fact that all values are in MST rather than MDT
    qc_df_pivot = qc_df_pivot.shift(periods=1)

    # #Resample to hourly to ensure a continuous time index
    qc_new_start_time = pd.Timestamp('2024-07-15 00:00:00')
    qc_new_end_time = pd.Timestamp('2024-08-01 23:00:00')
    #Create a new datetime index from new_start to the end of existing index with same frequency
    qc_new_index = pd.date_range(start=qc_new_start_time, end=qc_new_end_time, freq='1h')
    #Reindex the dataframe to include new rows
    qc_df_pivot = qc_df_pivot.reindex(qc_new_index)
    qc_df_pivot_parameter_sort = qc_df_pivot.drop(43102, axis=1)

    #Get carbon numbers by using the mappings spreadsheet to get the formula and dividing by number of carbons
    qc_df_mapping_parameters = pd.read_csv(mappings_filepath)
    qc_df_mapping_parameters = qc_df_mapping_parameters.drop([0])
    qc_carbon_number = qc_df_mapping_parameters['formula'].str.extract(r'C(\d{1,2})')[0]
    qc_carbon_number = qc_carbon_number.fillna(1).astype(int)
    qc_df_udaq_ppbv = qc_df_pivot_parameter_sort.divide(qc_carbon_number.values, axis=1)

    # Create mapping of Parameter Code to UDAQ variable name
    qc_param_to_var = qc_df_mapping_parameters.set_index('Parameter Code')['UDAQ_Variable'].to_dict()

    # # Rename columns in df_data
    qc_df_udaq_vals = qc_df_udaq_ppbv.rename(columns=qc_param_to_var)
    qc_df_udaq_vals.columns.name = None

    #Plot the first 3 species between the UDAQ Verbose and QC to see how they match up
    for spec in qc_df_udaq_vals.columns[0:3]:
        fig, ax = plt.subplots(figsize = (16,4), constrained_layout=True)
        ax.plot(qc_df_udaq_vals.index, df_udaq_vals[spec], linestyle = 'solid', color = 'g', marker = '+', label = 'Verbose')
        ax.plot(qc_df_udaq_vals.index, qc_df_udaq_vals[spec], linestyle = 'solid', color='m', marker='x',label='UDAQ QC', alpha = 0.7)

        pts1 = ax.scatter(qc_df_udaq_vals.index, df_udaq_vals[spec], marker = '+', color='g', label = 'Verbose', alpha = 0.2)
        pts2 = ax.scatter(qc_df_udaq_vals.index, qc_df_udaq_vals[spec], marker='x', color='m', label = 'UDAQ QC', alpha = 0.2)

        #Mark midnight for every day
        midnight_vals = []
        for midnight_idx in range(0,len(qc_df_udaq_vals.index),24):
            midnight_vals.append(qc_df_udaq_vals.index[midnight_idx])
        for day_pos in midnight_vals:
            ax.axvline(day_pos, color = 'black', linestyle = 'dotted', alpha = 0.7)
        #Set x ticks
        ax.xaxis.set_major_locator(mdates.DayLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
        # Minor ticks: every 3 hours
        ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21]))
        # Rotate and format tick labels
        ax.tick_params(axis='x', which='major')
        ax.tick_params(axis='x', which='minor', length=3, color='gray')
        #ax.grid(True, which='both')
        #plt.yticks(np.arange(10,100,10))
        plt.ylabel(spec + ' Concentration (ppb)')
        plt.xlabel('Date')
        plt.margins(x=0)

        handles, labels = plt.gca().get_legend_handles_labels()
        order = [0,1]
        plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])

        plt.title(spec + ' Comparison Between UDAQ and Mobile Lab Measurements')
        #plt.savefig(dirpath + '/Compare_UDAQ_MobileLab/plots/hawthorne_udaq_mobilelab_o3_comparison_full_campaign.png', dpi =300)
        plt.show()
def formaldehyde_extract(time_interval):
    """
    The USOS Campaign did not have many formaldehyde measurements taken. Here, we consider replacing the formaldhyde measurements with those taken by UDAQ.
    """
    load_15min_iwas_merge = merged_data_dir_15min + 'all_CSL_MobileLab_Parked_rev15min_iWASupdated.nc'
    ds_newmerge = xr.open_dataset(load_15min_iwas_merge)
    #plot existing NOAA HCHO data
    #plot existing UDAQ HCHO data
    #read UDAQ Formaldehyde data, provided by Nell Schafer (CU Boulder/NOAA) and Bart (UDAQ)
    udaq_formaldehyde_load = hawthorne_data_dir + 'hw_zero_corrected_data_formaldehyde.csv'
    df_udaq_formaldehyde = pd.read_csv(udaq_formaldehyde_load, index_col='dt', parse_dates=True)
    
    #shift data by one hour
    df_udaq_formaldehyde = df_udaq_formaldehyde.shift(freq= '1h')
    df_udaq_formaldehyde.index = df_udaq_formaldehyde.index.rename('time_local')
    df_udaq_formaldehyde['time_local']=df_udaq_formaldehyde.index
    #Change index to UTC time
    df_udaq_formaldehyde.index = df_udaq_formaldehyde.index + pd.Timedelta(hours=6)
    df_udaq_formaldehyde['time_UTC']=df_udaq_formaldehyde.index
    df_udaq_formaldehyde.index = df_udaq_formaldehyde.index.rename('time_UTC')
    df_udaq_formaldehyde_usos_only = df_udaq_formaldehyde.sort_index().loc['2024-07-15 00:00:00':'2024-08-18 23:59:00']

    df_udaq_formaldehyde_revised = df_udaq_formaldehyde_usos_only.copy()
    keep_colnames = ['H2CO_Values', 'H2CO_Corrected', 'time_local', 'time_UTC']
    df_udaq_formaldehyde_revised = df_udaq_formaldehyde_revised[keep_colnames]

    # Get the average native sampling frequency in total seconds:
    tseries = df_udaq_formaldehyde_revised.index.to_series()
    min_sep = int(np.round(tseries.diff().median().total_seconds()))
    step_S = time_interval
    
    new_start_time = pd.Timestamp('2024-07-15 00:00:00')
    new_end_time = pd.Timestamp('2024-08-18 23:45:00')
    dts = pd.date_range(new_start_time, new_end_time, freq=str(min_sep) + 's')

    dfn = df_udaq_formaldehyde_revised.reindex(dts, method='nearest', fill_value=np.nan)

    # Take a centered boxcar average around the 900s avg. (for numerical columns only)!    
    #NOTE: .mean() handles Nans like np.nanmean() in this context!!! 
    df_nums=dfn.select_dtypes(exclude=['datetime64'])
    df_nums_new = df_nums.rolling(str(int(step_S)) + 's').mean().resample(str(step_S) + 's').mean()
    dtss=df_nums_new.index.rename('time_UTC')

    df_nonums=dfn.select_dtypes(include=['datetime64'])
    df_nonums_new = df_nonums.reindex(dtss, method='nearest', fill_value=np.nan)
    df_nonums_new.index = df_nonums_new.index.rename('time_UTC')
    
    #our combined df_new dataframe now has our averaged 15 min intervals (in UTC)
    #The column 'time_UTC' is used mostly to check that the index is correct (they should match)
    #while time_local should be UTC -6 hrs to represent Mountain Daylight Time.
    df_new_formaldehyde=pd.concat([df_nums_new, df_nonums_new], axis=1, join="inner")
    df_new_formaldehyde.index = df_new_formaldehyde.index.index.rename('time_UTC')
    print(df_new_formaldehyde)
#endregion

def resave_ozone_15min_reindexed(filename):
    udaq_o3_filepath = dirpath + 'Hawthorne_data/data/from_Bart/hawthorne_udaq_o3_2024.csv'
    df_udaq_o3_load = pd.read_csv(udaq_o3_filepath, index_col = 'Date')
    index_series = pd.Series(df_udaq_o3_load.index)
    index_fixed = index_series.apply(lambda x: x if ':' in x else f'{x} 00:00')
    df_udaq_o3_load.index = pd.to_datetime(index_fixed, format='%m/%d/%Y %H:%M')

    #Shift all values by 1 hour to account for the fact that all values are in MST rather than MDT
    df_udaq_o3_load = df_udaq_o3_load.shift(periods=1)
    df_udaq_o3_load.index = df_udaq_o3_load.index.tz_localize('America/Denver')
    print('df_udaq_o3_load: ', df_udaq_o3_load)
    new_start_time = pd.Timestamp('2024-07-14 18:00:00', tz = 'America/Denver')
    new_end_time = pd.Timestamp('2024-08-18 17:45:00', tz = 'America/Denver')
    new_index = pd.date_range(start=new_start_time, end=new_end_time, freq='15min', tz = 'America/Denver')

    df_udaq_o3 = df_udaq_o3_load.reindex(new_index)
    print('df_udaq_o3: ', df_udaq_o3)
    df_udaq_o3['O3'] = df_udaq_o3['O3']*1000 #ppm to ppb conversion
    df_udaq_o3['O3'] = df_udaq_o3['O3'].mask(df_udaq_o3['O3'] >= 120) #remove values over 120 ppb, likely unusable
    df_udaq_o3.index.name = 'time_local'

    savepath = dirpath + 'Hawthorne_data/data/script_output/' + filename + '.csv'
    df_udaq_o3.to_csv(savepath)
    print('Saved to:' + savepath)
def resave_other_species_15min_reindexed(parameter_codes, species_names, species_names_for_files, filename):
    #Quality control file; only has Jul 01 to Aug 01
    udaq_qc_voc_filepath = dirpath + 'Hawthorne_data/data/from_Bart/udaq_QC_07012024_08012024.csv'

    #Read UDAQ QC VOC measurements file
    df_data_qc = pd.read_csv(udaq_qc_voc_filepath) 
    #Get only the data at Hawthorne
    df_udaq_hawthornedates_usos_only_qc = df_data_qc.loc[(df_data_qc['StationSym'] == 'HW')]
    #set index to datetimeindex
    df_udaq_hawthornedates_usos_only_qc.set_index(['dt'], inplace = True) 
    df_udaq_hawthornedates_usos_only_qc.index = pd.to_datetime(df_udaq_hawthornedates_usos_only_qc.index)

    #Shift all values by 1 hour to account for the fact that all values are in MST rather than MDT
    df_udaq_hawthornedates_usos_only_qc = df_udaq_hawthornedates_usos_only_qc.shift(periods=1)

    # # #Get only the data during the USOS campaign
    df_udaq_hawthornedates_usos_only_qc = df_udaq_hawthornedates_usos_only_qc.sort_index().loc['2024-07-14 18:00:00':'2024-08-18 23:00:00']

    # #Change dataframe to the parameters for the species we need
    df_udaq_species_needed_qc = df_udaq_hawthornedates_usos_only_qc.loc[df_udaq_hawthornedates_usos_only_qc['Parameter'].isin(parameter_codes)]
    #print(df_udaq_species_needed_qc)

    #Change the sorting so that we only get the sample value and dateindex, take out unnecessary columns
    df_reshaped_udaq_species_qc = df_udaq_species_needed_qc.groupby([df_udaq_species_needed_qc.index, 'Parameter'])['Sample Value'].first().unstack()
    df_reshaped_udaq_species_qc = df_reshaped_udaq_species_qc.rename(columns=dict(zip(parameter_codes,species_names)))
    
    df_reshaped_udaq_species_qc.index = df_reshaped_udaq_species_qc.index.tz_localize('America/Denver')

    new_start_time = pd.Timestamp('2024-07-14 18:00:00')
    new_end_time = pd.Timestamp('2024-08-18 17:45:00')
    new_index = pd.date_range(start=new_start_time, end=new_end_time, freq='15min', tz = 'America/Denver')
    df_reshaped_udaq_species_qc = df_reshaped_udaq_species_qc.reindex(new_index)
    df_reshaped_udaq_species_qc.index.name = 'time_local'

    df_reshaped_udaq_species_qc['NO'] = df_reshaped_udaq_species_qc['NO'].mask(df_reshaped_udaq_species_qc['NO'] < 0) #remove negative values

    #df_reshaped_udaq_species_qc['UDAQ CO'] = df_reshaped_udaq_species_qc['UDAQ CO']*1000 #ppm to ppb conversion
    for spec in range(0, len(species_names)):
        savepath = dirpath + 'Hawthorne_data/data/script_output/hawthorne_udaq_' + species_names_for_files[spec] + filename + '.csv'
        df_reshaped_udaq_species_qc[species_names[spec]].to_csv(savepath)
        print('Saved to:' + savepath)

def resave_all_voc_species_15min_reindexed(filename):
    udaq_voc_filepath = dirpath + 'Hawthorne_data/data/from_Bart/Verbose.csv'

    #Read UDAQ VOC measurements file
    df_data = pd.read_csv(udaq_voc_filepath) 
    #Get only the data at Hawthorne
    df_udaq_hawthornedates_usos_only = df_data.loc[(df_data['StationSym'] == 'HW')]
    #set index to datetimeindex
    df_udaq_hawthornedates_usos_only.set_index(['dt'], inplace = True) 
    df_udaq_hawthornedates_usos_only.index = pd.to_datetime(df_udaq_hawthornedates_usos_only.index)

    #Shift all values by 1 hour to account for the fact that all values are in MST rather than MDT
    df_udaq_hawthornedates_usos_only = df_udaq_hawthornedates_usos_only.shift(periods=1)

    #Get only the data during the USOS campaign
    df_udaq_hawthornedates_usos_only = df_udaq_hawthornedates_usos_only.sort_index().loc['2024-07-14 18:00:00':'2024-08-18 17:45:00']

    #View Parameter code values as only an int (they print with a .0 at the end initially)
    df_udaq_hawthornedates_usos_only.Parameter = df_udaq_hawthornedates_usos_only.Parameter.astype(int)
    sorted_parameters = sorted(df_udaq_hawthornedates_usos_only.Parameter.unique())
                                
    #Pivot the data
    df_pivot = df_udaq_hawthornedates_usos_only.pivot_table(index=df_udaq_hawthornedates_usos_only.index, columns='Parameter', values='Sample Value', dropna=False)

    new_start_time = pd.Timestamp('2024-07-14 18:00:00')
    new_end_time = pd.Timestamp('2024-08-18 17:45:00')
    #Create a new datetime index from new_start to the end of existing index with same frequency
    new_index = pd.date_range(start=new_start_time, end=new_end_time, freq='15min')
    #Reindex the dataframe to include new rows
    df_pivot = df_pivot.reindex(new_index, method='nearest', fill_value=np.nan, tolerance = '8min')
    df_pivot_parameter_sort = df_pivot[sorted(df_pivot.columns)]
    #Drop total NMVOCs
    df_pivot_parameter_sort = df_pivot_parameter_sort.drop(43102, axis=1)

    #Get carbon numbers by using the mappings spreadsheet to get the formula and dividing by number of carbons
    df_mapping_parameters = pd.read_csv(mappings_filepath)
    df_mapping_parameters = df_mapping_parameters.drop([0])
    carbon_number = df_mapping_parameters['formula'].str.extract(r'C(\d{1,2})')[0]
    carbon_number = carbon_number.fillna(1)
    carbon_number = carbon_number.astype(int)
    df_udaq_ppbv = df_pivot_parameter_sort.divide(carbon_number.values, axis=1)

    # Create mapping of Parameter Code to UDAQ variable name
    param_to_var = df_mapping_parameters.set_index('Parameter Code')['UDAQ_Variable'].to_dict()

    #Add a column for UTC Time, make it the first column (after the index)
    df_udaq_ppbv.index = df_udaq_ppbv.index.tz_localize('America/Denver')
    df_udaq_ppbv['time_UTC'] = df_udaq_ppbv.index.tz_convert('UTC')
    df_udaq_ppbv.insert(0, 'time_UTC', df_udaq_ppbv.pop('time_UTC'))

    # # Rename columns in df_data
    df_udaq_vals = df_udaq_ppbv.rename(columns=param_to_var)
    df_udaq_vals.columns.name = None
    df_udaq_vals.index.name = 'time_local'
    savepath = dirpath + 'Hawthorne_data/data/' + filename + '.csv'
    df_udaq_vals.to_csv(savepath)
    print('Saved to:' + savepath)
def resave_terpenes_15min_reindexed(filename):
    terpenes_udaq_filepath = hawthorne_data_dir + 'modified_from_Bart_raw/terpenes_HW_20240101_20241231.csv'
    terpenes_udaq_load = pd.read_csv(terpenes_udaq_filepath, index_col = 'Date', parse_dates=True)
    #Shift all values by 1 hour to account for the fact that all values are in MST rather than MDT
    terpenes_udaq_load = terpenes_udaq_load.shift(periods=1)
    terpenes_udaq_load.index = terpenes_udaq_load.index.tz_localize('America/Denver')
    new_start_time = pd.Timestamp('2024-07-14 18:00:00')
    new_end_time = pd.Timestamp('2024-08-18 17:45:00')
    new_index = pd.date_range(start=new_start_time, end=new_end_time, freq='15min', tz = 'America/Denver')
    terpenes_udaq = terpenes_udaq_load.reindex(new_index)
    terpenes_udaq = terpenes_udaq.rename(columns = {'ISOPRENE':'Isoprene', 'ALPHA-PINENE':'Alphapinene', 'BETA-PINENE':'Betapinene'})

    terpenes_udaq.replace('AS', np.nan, inplace=True)
    terpenes_udaq = terpenes_udaq.apply(pd.to_numeric, errors="coerce")
    #Carbon number calculation
    carbon_numbers_terpenes = [5, 10, 10]
    terpenes_udaq = terpenes_udaq.div(carbon_numbers_terpenes, axis=1)

    terpenes_udaq.index.name = 'time_local'
    savepath = dirpath + 'Hawthorne_data/data/script_output/' + filename + '.csv'
    terpenes_udaq.to_csv(savepath)
    print('Saved to:' + savepath)


#Handle formaldehyde data
#If F0AM runs are changed to any interval other than 15 minutes, need to change the step_S value
def formaldehyde_data_15min(filename):
    # read UDAQ Formaldehyde data, provided by and Bart (UDAQ), which is in UTC
    udaq_formaldehyde_load = hawthorne_data_dir + 'from_Bart/hw_zero_corrected_data_formaldehyde.csv'
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
    step_S = 15*60

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
    df_avg_formaldehyde = df_avg_formaldehyde.loc['2024-07-14 18:00:00':'2024-08-18 17:45:00']

    savepath = dirpath + 'Hawthorne_data/data/script_output/' + filename + '.csv'
    df_avg_formaldehyde.to_csv(savepath)
    print('Saved to:' + savepath)

    # # Set any inf, neg. inf, and negative values to NaN
    # df_avg_formaldehyde['Formaldehyde'] = df_avg_formaldehyde['Formaldehyde'].replace([np.inf, -np.inf], np.nan).mask(df_avg_formaldehyde['Formaldehyde'] < 0, np.nan)
    
    
    
# resave_hourly_ozone(
#     filename = 'hawthorne_udaq_o3_2024_timezone_updated'
# )

# resave_ozone_15min_reindexed(
#     filename = 'hawthorne_udaq_o3_2024_15min_reindexed_timezone_updated'
# )

# resave_one_other_species_15min_reindexed(
#     parameter_codes = [42101],
#     species_names = ['UDAQ CO'],
#     filename = 'hawthorne_udaq_co_07152024_08012024_15min_reindexed_timezone_updated'
# )

# resave_other_species_15min_reindexed(
#     parameter_codes = [42601, 42602, 42612],
#     species_names = ['NO', 'NO2', 'NOy'],
#     species_names_for_files = ['no', 'no2', 'noy'],
#     filename = '_07152024_08012024_15min_reindexed_timezone_updated'
# )

# resave_one_other_species_15min_reindexed(
#     parameter_codes = [42602],
#     species_names = ['UDAQ NO2'],
#     filename = 'hawthorne_udaq_no2_07152024_08012024_15min_reindexed_timezone_updated'
# )
# resave_one_other_species_15min_reindexed(
#     parameter_codes = [42612],
#     species_names = ['UDAQ NOy'],
#     filename = 'hawthorne_udaq_noy_07152024_08012024_15min_reindexed_timezone_updated'
# )

resave_terpenes_15min_reindexed(
    filename = 'hawthorne_udaq_isoprene_alpha_beta_pinene_07152024_08012024_15min_reindexed_timezone_updated'
)

# formaldehyde_data_15min(filename = 'hawthorne_udaq_Formaldehyde_15min_reindexed_timezone_updated')

# resave_all_voc_species_hourly(
#     filename = 'hawthorne_udaq_all_vocs_hourly_timezone_carbon_number_updated'
# )

# resave_all_voc_species_15min_reindexed(
#     filename = 'hawthorne_udaq_all_vocs_15min_timezone_carbon_number_updated'
# )

# formaldehyde_extract(
#     time_interval = 15*60
# )