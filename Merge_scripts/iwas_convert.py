# When the iWAS took measurements for the USOS campaign, it accidentally compiled all days of data into one file (for 07152024).
# Measurements were VOCs in units of ppbv, taken at approximately once every two hours but sometimes more frequently
# This program:
# - reads the raw icartt files for iWAS measurements
# - takes the hourly average, fills any missing times with NaN values
# - renames the variable name to what's used in the revised USOS parked spreadsheet for the revised merges of all other measurements
# - reads the revised 1 hr merges from running revise_USOS_merges.py and replaces those incorrect values with the correct hourly average iWAS VOCs

# Must have updated 1 hr revised merges from running revise_USOS_merges.py first!

import os
import sys
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import matplotlib as mpl
import xarray as xr
import matplotlib.dates as mdates

sys.path.insert(0,'/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/Merge_scripts/icartt_read_and_merge/')
# from icartt_read_and_merge import icartt_merger
import ict_utils 

def read_icartt(icartt_file: str, flt_num: int = -99, meta: dict = {},
                instr_name_prefix: bool = False, add_file_no: bool = False,
                delimiter: str = ',', line1_delim:str =',',
                meta_dict_replace_with: str = '_'):
    """Parse a single ICARTT file to a pandas dataframe."""
   
    # Get the header row number from the ICARTT.
    with open(icartt_file, "r", errors="ignore") as f : # ,encoding='ANSI') as f:
        header_row = int(f.readlines()[0].split(line1_delim)[0]) - 1

    # Parse the table starting where data begins (e.g. after the header).engine='python',
    print(icartt_file,header_row)
    df = pd.read_csv(icartt_file, header=header_row, delimiter=delimiter)#, encoding='utf-8')
    
    # Set possible error values to NaNs.
    df.replace(-9, np.nan, inplace=True)
    df.replace(-99, np.nan, inplace=True)
    df.replace(-999, np.nan, inplace=True)
    df.replace(-9999, np.nan, inplace=True)
    df.replace(-99999, np.nan, inplace=True)
    df.replace(-999999, np.nan, inplace=True)
    
    # Strip leading/tailing white space around variable names
    df.columns = [c.strip() for c in list(df.columns)]

    # Build/ append metadata from ICARTT to a dictionary file, add prefix to columns 
    df,meta = ict_utils.build_meta_dict(df, icartt_file, meta=meta,  flt_num=flt_num,
                            instr_name_prefix=instr_name_prefix,
                            line1_delim=line1_delim, add_file_no=add_file_no,
                            replace_with=meta_dict_replace_with)

    if add_file_no is True:
        # Create a column same length as data that contains the file #
        sz = len(df[df.columns[0]])  # get appropriate length
        fnum_arr = np.full(shape=sz, fill_value=flt_num, dtype=int)
        df['Flight_N'] = fnum_arr

    return df, meta  # dataframe with data, and df with metadata

#Original starting time of midnight UTC on 07-15-2024 
df_iwas, meta_iwas = read_icartt(icartt_file = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/CampaignData_and_Merges/R0/CSL_MobileLab_Parked/raw/20240715/USOS-iWAS_MobileLabGround_20240715_R1.ict')
base_dt = pd.Timestamp('2024-07-15 00:00:00')
df_iwas['Start_dt_UTC'] = base_dt + pd.to_timedelta(df_iwas['iWAS_Start_UTC'], unit = 's')
df_iwas['Stop_dt_UTC'] = base_dt + pd.to_timedelta(df_iwas['iWAS_Stop_UTC'], unit = 's')
df_iwas['Mid_dt_UTC'] = base_dt + pd.to_timedelta(df_iwas['iWAS_Mid_UTC'], unit = 's')

start_rounded = (df_iwas['iWAS_Start_UTC']).round()
df_iwas['Start_dt_UTC_rounded'] = base_dt + pd.to_timedelta(start_rounded, unit = 's')
#print(df_iwas['Start_dt_UTC_rounded'].values)
#df_iwas['Start_dt_UTC_rounded'] = df_iwas['Start_dt_UTC_rounded'].dt.round('min')
#df_iwas.set_index(['Start_dt_UTC_rounded'], inplace = True)

#There are four obs taken quickly at 07/23/2024 23:00:00
# Re-saves only the first value taken
jul23 = pd.Timestamp("2024-07-23")
df_jul23_hr23 = (
    (df_iwas['Start_dt_UTC_rounded'].dt.date == jul23.date()) &
    (df_iwas['Start_dt_UTC_rounded'].dt.hour == 23))
first_val_jul23 = df_iwas.loc[df_jul23_hr23].index.min()
df_iwas_jul23_hr23_oneval = pd.concat([
    df_iwas.loc[~df_jul23_hr23],
    df_iwas.loc[[first_val_jul23]]
]).sort_index()
#resave rounded time as index
df_iwas_jul23_hr23_oneval.set_index('Start_dt_UTC_rounded', inplace = True)
#print(df_iwas_jul23_hr23_oneval['Acetone_ppbv'].loc['2024-08-14 23'])

# #get plots for how the various species look for 07/23, where we have the biggest issue of multiple measurements taken
# #around 07/23 23:00:00
# jul23 = pd.Timestamp("2024-07-23")
# jul22 = jul23 - pd.Timedelta(days=1)
# jul24   = jul23 + pd.Timedelta(days=2)   # +2 because end is exclusive

# df_window = df_iwas[(df_iwas.index >= jul22) & (df_iwas.index < jul24)]

# noplot_colnames = ['iWAS_Start_UTC', 'iWAS_Stop_UTC', 'iWAS_Mid_UTC', 'Start_dt_UTC', 'Stop_dt_UTC', 'Mid_dt_UTC']
# df_window_copy = df_window.copy()
# df_window_copy = df_window_copy.drop(columns=noplot_colnames)
# print(df_window_copy.index)
# for spec in df_window_copy.columns:
#     fig, ax = plt.subplots(figsize = (8,4), constrained_layout=True)
#     ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
#     plt.plot(df_window_copy.index, df_window_copy[spec], color = 'r', marker = 'o')
#     plt.xlabel('Date/Time (UTC)')
#     plt.ylabel('Concentration (ppb)')
#     plt.margins(x=0)
#     plt.title(spec)
#     plt.xlim([pd.Timestamp("2024-07-23 14:00:00"), pd.Timestamp("2024-07-24 01:00:00")])

full_index = pd.date_range(start = pd.Timestamp('2024-07-15 00:00:00'),
                           end=pd.Timestamp('2024-08-18 23:45:00'),
                           freq='15min')

# tseries = df_iwas_jul23_hr23_oneval.index.to_series()
# min_sep = int(np.round(tseries.diff().median().total_seconds()))
# lim = max(1, np.round(4350 / min_sep))
df_iwas_new = df_iwas_jul23_hr23_oneval.reindex(full_index, method='nearest', fill_value=np.nan, tolerance = '8min')
#print(df_iwas_new['Acetone_ppbv'].loc['2024-08-14'])

noplot_colnames = ['iWAS_Start_UTC', 'iWAS_Stop_UTC', 'iWAS_Mid_UTC', 'Start_dt_UTC', 'Stop_dt_UTC', 'Mid_dt_UTC']
df_iwas_original = df_iwas_jul23_hr23_oneval.copy()
df_iwas_original = df_iwas_original.drop(columns=noplot_colnames)
df_iwas_updated= df_iwas_new.copy()
df_iwas_updated = df_iwas_updated.drop(columns=noplot_colnames)

# df_iwas_updated.index.name = 'time_UTC'
# df_iwas_updated['time_UTC'] = df_iwas_updated.index

#Rename the iWAS variables into the correct ones used in the original merges
rev_spreadsheet = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/Merge_spreadsheets/revised_USOS_parked_vars_updated_062025_r1_was.xlsx'
meta_df= pd.read_excel(rev_spreadsheet, index_col=0).fillna('')

#Read the Original_varname in excel file. For columns of our dataframe, rename dataframe column name
map_df = pd.DataFrame({'ORIGINAL_VARNAME':meta_df['ORIGINAL_VARNAME'], 'NEW_VARNAME':meta_df['NEW_VARNAME']})
mapping = dict(zip(map_df['ORIGINAL_VARNAME'], map_df['NEW_VARNAME']))
df_iwas_updated = df_iwas_updated.rename(columns=mapping)
#rename index to 'time_UTC' to match netcdf file index name
df_iwas_updated.index.rename('time_UTC', inplace=True)
#print(df_iwas_updated)

# for spec in df_iwas_updated.columns[7:12]:
#     fig, ax = plt.subplots(figsize = (8,4), constrained_layout=True)
#     #ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
#     # plt.plot(df_iwas_original.index,df_iwas_original[spec], color = 'k', marker = 'o', label = 'Original')
#     plt.plot(df_iwas_updated.index, df_iwas_updated[spec], color = 'r', marker = 'x', label = 'Reindexed')
#     plt.xlabel('Date/Time (UTC)')
#     plt.ylabel('Concentration (ppb)')
#     plt.margins(x=0)
#     plt.title(spec)
#     plt.xlim([pd.Timestamp("2024-07-15"), pd.Timestamp("2024-07-16")])


# #iWAS values are saving with a new dimension called 'index' so make sure it is instead called 'time_UTC'
# ##Replace iWAS values from the merge with our new avgs
ds_nc = xr.open_dataset('/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/CampaignData_and_Merges/R0/CSL_MobileLab_Parked/merged/rev_15min/all_CSL_MobileLab_Parked_rev15min.nc')

for var in df_iwas_updated.columns:
    if var in ds_nc.data_vars:
        # Align the DataFrame to the NetCDF time
        values = df_iwas_updated[var].reindex(ds_nc.time_UTC).to_numpy()

        # Remove old variable completely
        ds_nc = ds_nc.drop_vars(var)

        # Assign as a DataArray with explicit dims and coords
        da_new = xr.DataArray(
            data=values,
            dims=["time_UTC"],           # specify the correct dimension
            coords={"time_UTC": ds_nc.time_UTC},  # assign coordinate explicitly
            name=var
        )
        ds_nc[var] = da_new


# for var in df_iwas_updated.columns:
#     if var in ds_nc.data_vars:
#         if 'index' in ds_nc[var].dims:
#             ds_nc = ds_nc.drop_vars(var)
#         values = df_iwas_updated[var].reindex(ds_nc.time_UTC).to_numpy()
#         ds_nc[var] = ('time_UTC',values)
#print(ds_nc)
# print(ds_nc.dims)
# print(ds_nc.coords)
savedir_path = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/CampaignData_and_Merges/R0/CSL_MobileLab_Parked/merged/rev_15min/'
save_ncfilename = 'all_CSL_MobileLab_Parked_rev15min_iWASupdated.nc'
full_savepath = savedir_path + save_ncfilename
ds_nc.to_netcdf(full_savepath, format = 'NETCDF4',mode='w')
print('Saved netCDF file with updated iWAS measurements to: ', full_savepath)

# ds_new = xr.open_dataset(full_savepath)
# print(ds_new)

