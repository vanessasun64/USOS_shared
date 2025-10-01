import os 
import sys 
import re
import inspect 
import numpy as np 
import pandas as pd
import xarray as xr

current_dir = os.path.dirname(os.path.abspath(__file__))
print(current_dir)
global_scripts_path = os.path.abspath(os.path.join(current_dir, "..", "global_scripts"))
print(global_scripts_path)
sys.path.insert(0,global_scripts_path)
from dirpath import filepath_source
dirpath = filepath_source('CHPC')
savepath= dirpath + '/Hawthorne_data/'
filename='Hawthorne_EPA_match'

def hawthorne_usos_extract_species(udaq_voc_filepath, epa_parameter_info):
    df_data = pd.read_csv(udaq_voc_filepath)
    df_hawthorne_usos_only = df_data.loc[(df_data['StationSym'] == 'HW')] #Get only the data at Hawthorne
    #set index to datetimeindex
    df_hawthorne_usos_only.set_index(['dt'], inplace = True) 
    df_hawthorne_usos_only.index = pd.to_datetime(df_hawthorne_usos_only.index)

    #Get only the data during the USOS campaign
    df_hawthorne_usos_only = df_hawthorne_usos_only.sort_index().loc['2024-07-14 00:00:00':'2024-08-18 23:00:00']

    df_hawthorne_usos_only.Parameter = df_hawthorne_usos_only.Parameter.astype(int)
    sorted_parameters = sorted(df_hawthorne_usos_only.Parameter.unique())

    df_species = pd.DataFrame(sorted_parameters, columns = ['Parameters'])

    df_epa_parameters = pd.read_csv(epa_parameter_info)
    matching_parameters = df_epa_parameters[df_epa_parameters['Parameter Code'].isin(df_species['Parameters'])]

    df_matching_species = matching_parameters.sort_values(by='Parameter Code').reset_index(drop=True)
    df_matching_species = df_matching_species.drop(columns = ['Parameter Abbreviation', 'Still Valid', 'Round or Truncate'])
    df_matching_species = df_matching_species.rename(columns={'Parameter':'EPA Species'})

    excel_file = savepath + filename + '.xlsx'
    # Write the dataframe to an excel file:
    df_matching_species.to_excel(excel_file, index = False)
    print('Output df of parameters saved at: ' + excel_file)


hawthorne_usos_extract_species(
    udaq_voc_filepath = dirpath + 'Hawthorne_data/Verbose.csv',
    epa_parameter_info = dirpath + 'Hawthorne_data/epa_parameters_official.csv'
)