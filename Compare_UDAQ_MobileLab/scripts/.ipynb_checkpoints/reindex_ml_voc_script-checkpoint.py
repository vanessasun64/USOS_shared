#Resave VOC data for comparing UDAQ and Mobile Lab
#ML data from '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/CampaignData_and_Merges/R0/CSL_MobileLab_Parked/merged/rev_1hr/all_CSL_MobileLab_Parked_rev1hr_iWASupdated.nc'
#has index of 2024-07-14 18:00:00 to 2024-08-18 17:00:00
#UDAQ data from '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/Hawthorne_data/data/hawthorne_udaq_all_vocs_hourly_timezone_carbon_number_updated.csv'
#has index of 2024-07-15 00:00:00 to 2024-08-18 23:00:00
#Reindexed to 2024-07-15 00:00:00 to 2024-08-18 23:00:00 for comparing the overlap

import xarray as xr
import numpy as np
import pandas as pd

#Filepaths for loading
dirpath = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/'
resaved_dir = dirpath + 'Compare_UDAQ_MobileLab/resaved_data/'

#Get names of variables needed from Mobile Lab and UDAQ dataset
voc_mappings_filepath = dirpath + 'Hawthorne_data/mappings/manually_edited/UDAQ_Hawthorne_CRACMM_GEOSCHEM_CB6r5h_mapped_updated_11172025.csv'
df_voc_mapping_parameters = pd.read_csv(voc_mappings_filepath)
df_voc_mapping_parameters = df_voc_mapping_parameters.drop([0]) #drop Total NMVOCs
usos_voc_vars = [str(spec) for spec in df_voc_mapping_parameters['USOS Mapping'].dropna()]
usos_voc_vars.append('time_local')
mapping = dict(zip(df_voc_mapping_parameters['USOS Mapping'], df_voc_mapping_parameters['UDAQ_Variable']))


ml_voc_data = xr.open_dataset('/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/CampaignData_and_Merges/R0/CSL_MobileLab_Parked/merged/rev_1hr/all_CSL_MobileLab_Parked_rev1hr_iWASupdated.nc')
ml_voc_data_vars_subset = ml_voc_data[usos_voc_vars]
df_ml_voc_data = ml_voc_data_vars_subset.to_dataframe()
df_ml_voc_data.set_index('time_local', inplace=True)

#Turn any negatives into NaNs
for spec in df_ml_voc_data.columns:
    df_ml_voc_data[spec] = df_ml_voc_data[spec].mask(df_ml_voc_data[spec] < 0, np.nan)

df_ml_voc_data.columns = [
    mapping.get(col) if not pd.isna(col) else mapping.get(np.nan, col)
    for col in df_ml_voc_data.columns
]

new_start_time = pd.Timestamp('2024-07-15 00:00:00')
new_end_time = pd.Timestamp('2024-08-18 23:00:00')
#Create a new datetime index from new_start to the end of existing index with same frequency
new_index = pd.date_range(start=new_start_time, end=new_end_time, freq='1h')
#Reindex the dataframe to include new rows
df_ml_voc_data = df_ml_voc_data.reindex(new_index)
df_ml_voc_data.index.name = 'time_local'

filename = 'ml_hourly_voc_overlap_reindexed'

savepath = resaved_dir + filename + '.csv'
df_ml_voc_data.to_csv(savepath)
print('Saved to:' + savepath)