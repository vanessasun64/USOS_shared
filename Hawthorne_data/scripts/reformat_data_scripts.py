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
    udaq_qc_voc_filepath = dirpath + 'Hawthorne_data/data/udaq_QC_07012024_08012024.csv'

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
    new_start_time = pd.Timestamp('2024-07-15 00:00:00')
    new_end_time = pd.Timestamp('2024-08-18 23:00:00')
    new_index = pd.date_range(start=new_start_time, end=new_end_time, freq='1h')
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

def resave_all_voc_species_15min_reindexed(filename):
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
    new_end_time = pd.Timestamp('2024-08-18 23:45:00')
    #Create a new datetime index from new_start to the end of existing index with same frequency
    new_index = pd.date_range(start=new_start_time, end=new_end_time, freq='15min')
    #Reindex the dataframe to include new rows
    df_pivot = df_pivot.reindex(new_index, method='nearest', fill_value=np.nan, tolerance = '8min')
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


# resave_hourly_ozone(
#     filename = 'hawthorne_udaq_o3_2024_timezone_updated'
# )

# resave_one_other_species(
#     parameter_codes = [42602],
#     species_names = ['UDAQ NO2'],
#     filename = 'hawthorne_udaq_no2_07152024_08012024_timezone_updated'
# )

# resave_all_voc_species_hourly(
#     filename = 'hawthorne_udaq_all_vocs_hourly_timezone_carbon_number_updated'
# )

resave_all_voc_species_15min_reindexed(
    filename = 'hawthorne_udaq_all_vocs_15min_timezone_carbon_number_updated'
)