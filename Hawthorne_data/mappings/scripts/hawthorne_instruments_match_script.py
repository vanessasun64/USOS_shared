#This program takes spreadsheets given by UDAQ on all the instruments set up at Hawthorne and matches the EPA Parameter Code to
# give us a spreadsheet containing the Parameter Code, Species, Method Code, info about instrument used to take measurement, and
# units for all species measured at Hawthorne by UDAQ.

import pandas as pd

# All instruments UDAQ uses for measurements at Hawthorne
hawthorne_instruments_data_load = pd.read_csv('/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/Hawthorne_data/manually_edited/hawthorne_instruments_with_parameter_code_and_method_code.csv')
all_instruments_info = pd.read_csv('/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/Hawthorne_data/given_data/all_instrument_methods.csv')
# Merge df1 with df2 on Parameters and Methods
merged_instruments_data = pd.merge(hawthorne_instruments_data_load, all_instruments_info, on=["ParameterCode", "MethodCode"], how="left")

# Select only the desired columns
merged_instruments_data = merged_instruments_data[['ParameterCode', 'Parameter', 'MethodCode', 'RecordingMode', 'CollectionDescription', 'AnalysisDescription', 'CollectionType', 'Units']]

merged_instruments_data.to_csv('/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/Hawthorne_data/script_output/hawthorne_all_species_measured_by_udaq_instruments_info.csv',index=False)