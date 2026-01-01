
import os
import sys
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import matplotlib as mpl
import xarray as xr
import matplotlib.dates as mdates

from sklearn.linear_model import LinearRegression
import statsmodels.api as sm

sys.path.insert(0,'/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/Merge_scripts/icartt_read_and_merge/')
# from icartt_read_and_merge import icartt_merger
import ict_utils 

#region: filepaths
#Filepaths for loading
dirpath = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/'
merged_data_dir = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/CampaignData_and_Merges/R0/CSL_MobileLab_Parked/merged/'
merged_data_dir_15min = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/CampaignData_and_Merges/R0/CSL_MobileLab_Parked/merged/rev_15min/'
hawthorne_data_dir = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/Hawthorne_data/data/'
merge_scripts_dir = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/Merge_scripts/'
merge_scripts_plots_dir = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/Merge_scripts/plots/'

#endregion

#region: Plot formatting
mpl.rcParams['xtick.labelsize'] = 15
mpl.rcParams['ytick.labelsize'] = 15
mpl.rcParams['legend.fontsize'] = 16
mpl.rcParams['axes.labelsize'] = 18
mpl.rcParams['axes.titlesize'] = 28
# Set the font family to 'serif'
mpl.rcParams['font.family'] = 'serif'
# Specify preferred serif font (Computer Modern Roman is 'cmr10')
mpl.rcParams['font.serif'] = 'Lato' 
# Optionally, configure mathtext to use Computer Modern fonts as well
mpl.rcParams['mathtext.fontset'] = 'cm'
# Ensure minus signs are rendered correctly with CM fonts
mpl.rcParams['axes.unicode_minus'] = False
#endregion

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

def iwas_convert(savefilename):
    """
    # When the iWAS took measurements for the USOS campaign, it accidentally compiled all days of data into one file (for 07152024).
    # Measurements were VOCs in units of ppbv, taken at approximately once every two hours but sometimes more frequently
    # This program:
    # - reads the raw icartt files for iWAS measurements
    # - takes the measurement for every 15 minutes, fills any missing times with NaN values
    # - renames the variable name to what's used in the revised USOS parked spreadsheet for the revised merges of all other measurements
    # - reads the revised 15 min merges from running revise_USOS_merges.py and replaces those incorrect values with the correct hourly average iWAS VOCs

    # Must have updated 15 min revised merges from running revise_USOS_merges.py first!
    """

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


    #Remove four obs taken quickly at 07/23/2024 23:00:00, thanks to instrument issues
    jul23 = pd.Timestamp("2024-07-23")
    df_jul23_hr23 = (
        (df_iwas['Start_dt_UTC_rounded'].dt.date == jul23.date()) &
        (df_iwas['Start_dt_UTC_rounded'].dt.hour == 23))
    df_iwas = df_iwas[~(df_jul23_hr23)]
    df_iwas = df_iwas.set_index(['Start_dt_UTC_rounded'])
    #To check that reindexing was correct, use acetone as an example; compare with Line 95
    #print(df_iwas['Acetone_ppbv'].loc['2024-08-14'])

    #Change index to 15 min frequency, with each value going into its closest 15 minute interval
    full_index = pd.date_range(start = pd.Timestamp('2024-07-15 00:00:00'),
                            end=pd.Timestamp('2024-08-18 23:45:00'),
                            freq='15min')

    df_iwas_reindex = df_iwas.reindex(full_index, method='nearest', fill_value=np.nan, tolerance = '8min')
    #To check that reindexing was correct, use acetone as an example; compare with Line 86
    #print(df_iwas_new['Acetone_ppbv'].loc['2024-08-14'])

    #Drop columns not needed
    noplot_colnames = ['iWAS_Start_UTC', 'iWAS_Stop_UTC', 'iWAS_Mid_UTC', 'Start_dt_UTC', 'Stop_dt_UTC', 'Mid_dt_UTC']
    df_iwas_updated = df_iwas_reindex.drop(columns=noplot_colnames)
    #Rename the iWAS variables into the correct ones used in the original merges
    rev_spreadsheet = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/Merge_spreadsheets/revised_USOS_parked_vars_updated_062025_r1_was.xlsx'
    meta_df= pd.read_excel(rev_spreadsheet, index_col=0).fillna('')

    #Read the Original_varname in excel file. For columns of our dataframe, rename dataframe column name
    map_df = pd.DataFrame({'ORIGINAL_VARNAME':meta_df['ORIGINAL_VARNAME'], 'NEW_VARNAME':meta_df['NEW_VARNAME']})
    mapping = dict(zip(map_df['ORIGINAL_VARNAME'], map_df['NEW_VARNAME']))
    df_iwas_updated = df_iwas_updated.rename(columns=mapping)
    #rename index to 'time_UTC' to match netcdf file index name
    df_iwas_updated.index.rename('time_UTC', inplace=True)

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

    save_ncfilename = savefilename + '.nc'
    full_savepath = merged_data_dir_15min + save_ncfilename
    ds_nc.to_netcdf(full_savepath, format = 'NETCDF4',mode='w')
    print('Saved netCDF file with updated iWAS measurements to: ', full_savepath)

    ds_new = xr.open_dataset(full_savepath)
    print(ds_new)

def formaldehyde_averaging(time_interval):
    """
    The USOS Campaign did not have many formaldehyde measurements taken. Here, we consider replacing the formaldehyde measurements with those taken by UDAQ.
    """
    load_15min_iwas_merge = merged_data_dir_15min + 'all_CSL_MobileLab_Parked_rev15min_iWASupdated.nc'
    ds_newmerge = xr.open_dataset(load_15min_iwas_merge)
    df_ml_data = ds_newmerge[['HCHO_CRDS', 'time_local']].to_dataframe()
    df_ml_data['HCHO_CRDS'] = df_ml_data['HCHO_CRDS'].mask(df_ml_data['HCHO_CRDS'] < 0, np.nan)

    # #read UDAQ Formaldehyde data, provided by Nell Schafer (CU Boulder/NOAA) and Bart (UDAQ)
    udaq_formaldehyde_load = hawthorne_data_dir + 'hw_zero_corrected_data_formaldehyde.csv'
    df_udaq_formaldehyde = pd.read_csv(udaq_formaldehyde_load, index_col='dt', parse_dates=True)

    df_udaq_formaldehyde_usos_only = df_udaq_formaldehyde.sort_index().loc['2024-07-15 00:00:00':'2024-08-18 23:59:00']
    df_udaq_formaldehyde_usos_only.index = df_udaq_formaldehyde_usos_only.index.rename('time_UTC')
    keep_colnames = ['H2CO_Values', 'H2CO_Corrected']
    df_udaq_formaldehyde_usos_only = df_udaq_formaldehyde_usos_only[keep_colnames]

    #This averaging in this section is based off the same code used for averaging the icartt data
    #as in the function align2master_timeline in time_utils.py of the icartt_read_and_merge

    # Get the average native sampling frequency in total seconds:
    tseries = df_udaq_formaldehyde_usos_only.index.to_series()
    #Avg native sampling frequency
    min_sep = int(np.round(tseries.diff().median().total_seconds()))
    #Intended interval
    step_S = time_interval

    #Reindex to avg native sampling frequency (in this case, 1 minute)
    new_start_time = pd.Timestamp('2024-07-15 00:00:00')
    new_end_time = pd.Timestamp('2024-08-18 23:59:00')
    dts = pd.date_range(new_start_time, new_end_time, freq=str(min_sep) + 's')
    dfn = df_udaq_formaldehyde_usos_only.reindex(dts, method='nearest', fill_value=np.nan)

    # Take a centered boxcar average around the 900s (15 min) avg. (for numerical columns only)
    #NOTE: .mean() handles Nans like np.nanmean() in this context!!! 
    df_nums=dfn.select_dtypes(exclude=['datetime64'])
    df_nums_new = df_nums.rolling(str(int(step_S)) + 's').mean().resample(str(step_S) + 's').mean()

    #Part of original code but seems unnecessary here?
    # #Make sure our index is the correct interval:
    # dtss=df_nums_new.index
    # df_nonums=dfn.select_dtypes(include=['datetime64'])
    # df_nonums_new = df_nonums.reindex(dtss, method='nearest', fill_value=np.nan)
    # df_nonums_new.index = df_nonums_new.index.rename('time_UTC')
    
    # our new dataframe now has our averaged 15 min intervals (in UTC)
    # df_new_formaldehyde=pd.concat([df_nums_new, df_nonums_new], axis=1, join="inner")
    df_new_formaldehyde = df_nums_new
    df_new_formaldehyde.index = df_new_formaldehyde.index.rename('time_UTC')
    #print(df_new_formaldehyde)
    df_new_formaldehyde = df_new_formaldehyde.rename(columns = {'H2CO_Values':'H2CO_UDAQ_NoCorrection', 'H2CO_Corrected': 'H2CO_UDAQ_Corrected'})
    return df_ml_data, df_new_formaldehyde

def formaldehyde_correction_compare(time_interval_used):
    averaging_func_output = formaldehyde_averaging(time_interval = time_interval_used)
    df_ml_data = averaging_func_output[0]
    df_formaldehyde_averaged = averaging_func_output[1]
    fig, (ax1, ax2) = plt.subplots(2,1, figsize = (16,8), tight_layout=True)
    xlim_start_jul = pd.to_datetime('2024-07-15 00:00:00')
    xlim_end_jul = pd.to_datetime('2024-07-31 23:00:00')

    #ax1 is the first row of subplot, for July only
    ax1.plot(df_formaldehyde_averaged.index, df_formaldehyde_averaged['H2CO_UDAQ_NoCorrection'],linestyle = 'solid', color = 'g', marker = '+', label = 'UDAQ No Correction')
    ax1.plot(df_formaldehyde_averaged.index, df_formaldehyde_averaged['H2CO_UDAQ_Corrected'],linestyle = 'solid', color = 'b', marker = '+', label = 'UDAQ Corrected')
    ax1.plot(df_ml_data.index, df_ml_data['HCHO_CRDS'], linestyle = 'solid', color='m', marker='x',label='Mobile Lab', alpha = 0.7)

    #Set x ticks
    ax1.xaxis.set_major_locator(mdates.DayLocator())
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
    # Minor ticks: every 3 hours
    ax1.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21]))
    # Rotate and format tick labels
    ax1.tick_params(axis='x', which='major')
    ax1.tick_params(axis='x', which='minor', length=3, color='gray')
    #ax.grid(True, which='both')
    
    ax1.set_ylabel('Formaldehyde (ppb)')
    #ax1.set_xlabel('Date')
    ax1.margins(x=0)
    ax1.set_xlim([xlim_start_jul, xlim_end_jul])

    ax1.legend(loc = 'upper right')
    
    #ax2 is the second row of subplot, for August only
    ax2.plot(df_formaldehyde_averaged.index, df_formaldehyde_averaged['H2CO_UDAQ_NoCorrection'],linestyle = 'solid', color = 'g', marker = '+', label = 'UDAQ No Correction')
    ax2.plot(df_formaldehyde_averaged.index, df_formaldehyde_averaged['H2CO_UDAQ_Corrected'],linestyle = 'solid', color = 'b', marker = '+', label = 'UDAQ Corrected')
    ax2.plot(df_ml_data.index, df_ml_data['HCHO_CRDS'], linestyle = 'solid', color='m', marker='x',label='Mobile Lab', alpha = 0.7)

    #Set x ticks
    ax2.xaxis.set_major_locator(mdates.DayLocator())
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
    # Minor ticks: every 3 hours
    ax2.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21]))
    # Rotate and format tick labels
    ax2.tick_params(axis='x', which='major')
    ax2.tick_params(axis='x', which='minor', length=3, color='gray')

    ax2.set_ylabel('Formaldehyde (ppb)')
    ax2.set_xlabel('Date (UTC)')
    ax2.margins(x=0)

    ax2.set_xlim([pd.to_datetime('2024-08-01 00:00:00'), pd.to_datetime('2024-08-18 23:00:00')])
    ax2.legend(loc = 'upper right')

    #Mark midnight for every day
    midnight_vals = []
    for midnight_idx in range(0,len(df_ml_data.index),96):
        midnight_vals.append(df_ml_data.index[midnight_idx])
    for day_pos in midnight_vals:
        ax1.axvline(day_pos, color = 'black', linestyle = 'dotted')
        ax2.axvline(day_pos, color = 'black', linestyle = 'dotted')

    plt.savefig(merge_scripts_plots_dir + 'hawthorne_udaq_ml_' + 'comparison_july_aug_formaldehyde.png', dpi =300)
    plt.show()

    #estimate mean bias for No Correction vs Corrected
    corrected_varname = ['H2CO_UDAQ_NoCorrection', 'H2CO_UDAQ_Corrected']
    mb_avg_both = []
    mb_hrly_both = []
    mnb_avg_both = []
    mnb_hrly_both = []
    for corr in corrected_varname:
        df_form = pd.DataFrame({'obs':df_ml_data['HCHO_CRDS'], 'model':df_formaldehyde_averaged[corr]})
        df_form['difference']= df_form['model']-df_form['obs']
        mb_total = df_form['difference'].dropna().mean()
        mb_avg_both.append(mb_total)

        df_form['hour']=df_form.index.hour 
        hourly_MB = df_form.groupby(df_form.index.hour)['difference'].mean()
        mb_hrly_both.append(hourly_MB)

        df_mnb_div = df_form['difference']/df_form['obs']
        mnb_total = df_mnb_div.dropna().mean()
        mnb_avg_both.append(mnb_total)

        df_mnb_div = pd.DataFrame({'mnb_div':df_mnb_div})
        df_mnb_div['hour']=df_mnb_div.index.hour 
        mnb_hrly = df_mnb_div.groupby('hour')['mnb_div'].mean()
        mnb_hrly_both.append(mnb_hrly)


    plt.figure(figsize=(10,6), tight_layout = True)

    plt.plot(mb_hrly_both[0].index, mb_hrly_both[0], color='g', marker='.', label=f"Hrly MB No Correction (Min.={np.min(mb_hrly_both[0]):.2f} ppb)")
    plt.plot(mb_hrly_both[0].index, np.ones(len(mb_hrly_both[0].index))*mb_avg_both[0], linestyle = 'dashed', color='g',label=f"Avg. MB No Correction={mb_avg_both[0]:.2f} ppb")
    plt.plot(mb_hrly_both[1].index, mb_hrly_both[1], color='b', marker='.', label=f"Hrly MB Corrected (Min.={np.min(mb_hrly_both[1]):.2f} ppb)")
    plt.plot(mb_hrly_both[1].index, np.ones(len(mb_hrly_both[1].index))*mb_avg_both[1], linestyle = 'dashed', color='b',label=f"Avg. MB Corrected={mb_avg_both[1]:.2f} ppb")

    plt.hlines(y=0, xmin=mb_hrly_both[0].index[0], xmax= mb_hrly_both[0].index[len(mb_hrly_both[0])-1], linestyle='solid', color = 'k')

    plt.ylabel('Mean Bias (ppb)')
    plt.xlabel('Hour (UTC)')
    plt.xlim([0, 23])
    
    hour_range = np.arange(0,24,1)
    plt.xticks(hour_range)
    plt.yticks(np.arange(0, 0.4, 0.05))
    plt.ylim([0,0.45])
    
    plt.grid()
    plt.legend(loc='upper left')
    plt.savefig(merge_scripts_plots_dir + 'hawthorne_udaq_ml_' + 'meanbias_formaldehyde.png', dpi =300)
    plt.show()

    #Plotting Mean Normalized Bias
    plt.figure(figsize=(10,6), tight_layout = True)
    plt.plot(mnb_hrly_both[0].index, mnb_hrly_both[0]*100, color='g', marker='.', label=f"Hrly MNB No Correction")
    plt.plot(mnb_hrly_both[0].index, np.ones(len(mnb_hrly_both[0].index))*mnb_avg_both[0]*100, linestyle='dashed', color='g',label=f"Avg. MNB No Correction")
    plt.plot(mnb_hrly_both[1].index, mnb_hrly_both[1]*100, color='b', marker='.', label=f"Hrly MNB Corrected")
    plt.plot(mnb_hrly_both[1].index, np.ones(len(mnb_hrly_both[1].index))*mnb_avg_both[1]*100, linestyle='dashed', color='b',label=f"Avg. MNB Corrected")

    plt.hlines(y=0, xmin=mnb_hrly_both[0].index[0], xmax= mnb_hrly_both[0].index[len(mnb_hrly_both[0])-1], linestyle='solid', color = 'k')
            
    plt.ylabel('Mean Normalized Bias (%)')
    plt.xlabel('Hour (UTC)')
    plt.xlim([0, 23])
    plt.xticks(hour_range)
    plt.yticks(np.arange(0, 20, 1))
    plt.ylim([0,20])
    plt.grid()
    plt.legend(loc='upper right')
    plt.savefig(merge_scripts_plots_dir + 'hawthorne_udaq_ml_' + 'meannormalizedbias_formaldehyde.png', dpi =300)
    plt.show()

def bias_and_error_analysis(time_interval_used):
    averaging_func_output = formaldehyde_averaging(time_interval = time_interval_used)
    df_ml_data = averaging_func_output[0]
    df_formaldehyde_averaged = averaging_func_output[1]
    udaq_form_nocorr = df_formaldehyde_averaged['H2CO_UDAQ_NoCorrection']
    udaq_form_corr = df_formaldehyde_averaged['H2CO_UDAQ_Corrected']
    ml_form = df_ml_data['HCHO_CRDS']

    diff_nocorr = udaq_form_nocorr-ml_form
    diff_corr = udaq_form_corr-ml_form
    
    abs_diff_nocorr = abs(udaq_form_nocorr-ml_form)
    abs_diff_corr = abs(udaq_form_corr-ml_form)

    df_obs_modeled_nocorr= pd.DataFrame({'Observed': ml_form, 'Modeled': udaq_form_nocorr, 'Residual':udaq_form_nocorr-ml_form, 'Abs. Residual':abs(udaq_form_nocorr-ml_form)})
    df_clean_nocorr = df_obs_modeled_nocorr.dropna()
    df_obs_modeled_corr= pd.DataFrame({'Observed': ml_form, 'Modeled': udaq_form_corr, 'Residual':udaq_form_corr-ml_form, 'Abs. Residual':abs(udaq_form_corr-ml_form)})
    df_clean_corr = df_obs_modeled_corr.dropna()
    
    #Scatter Plot of ML vs UDAQ
    fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize = (12,6), tight_layout=True)
    ax1.scatter(x = ml_form, y = udaq_form_nocorr, color = 'g')
    ax1.set_xlabel('ML HCHO (ppb)')
    ax1.set_ylabel('UDAQ HCHO No Correction (ppb)')
    ax1.set_title('Model vs Obs')

    #Plot regression line
    O_reshaped_nocorr = df_clean_nocorr['Observed'].values.reshape(-1,1)
    model_nocorr = LinearRegression().fit(O_reshaped_nocorr, df_clean_nocorr['Modeled'].values)
    slope_nocorr = model_nocorr.coef_[0]
    intercept_nocorr = model_nocorr.intercept_
    # Predicted M based on regression
    M_pred_nocorr = model_nocorr.predict(O_reshaped_nocorr)
    # Regression line
    ax1.plot(df_clean_nocorr['Observed'].values, M_pred_nocorr, color='k', linewidth=2, label=f'Regression line: y={intercept_nocorr:.2f}+{slope_nocorr:.2f}x')
    # 1:1 line
    ax1.plot([df_clean_nocorr['Observed'].values.min(), df_clean_nocorr['Observed'].values.max()], [df_clean_nocorr['Observed'].values.min(), df_clean_nocorr['Observed'].values.max()], color='y', linestyle='--', label='1:1 line')
    ax1.legend(fontsize = 10)

    #Scatter Plot of Residuals vs Obs
    ax2.scatter(x = ml_form, y = diff_nocorr, color = 'g')
    ax2.set_xlabel('ML HCHO (ppb)')
    ax2.set_ylabel('Diff. M - O (ppb)')
    ax2.set_title('Residuals vs Obs')

    #Give Ordinary Least Squares regression summary
    #Looking at Coefficient of Observed which is the slope. If the p-value (P>|t|) for the Observed > 0.05 then the slope is not significantly different from zero and suggests that the residuals are constant 
    col_ones = sm.add_constant(df_clean_nocorr['Observed']) #adds column of ones to ind var
    ols_nocorr = sm.OLS(df_clean_nocorr['Residual'], col_ones).fit()
    print(ols_nocorr.summary())

    #Scatter Plot of Residuals vs Modeled
    ax3.scatter(x = udaq_form_nocorr, y = diff_nocorr, color = 'g')
    ax3.set_xlabel('UDAQ HCHO No Correction (ppb)')
    ax3.set_ylabel('Diff. M - O (ppb)')
    ax3.set_title('Residuals vs Modeled')
    plt.savefig(merge_scripts_plots_dir + 'hawthorne_udaq_ml_' + 'scatter_formaldehyde_nocorrection_ModelvsObs_and_Residuals.png', dpi =300)
    plt.show()

    #Scatter Plot of Abs. Residuals vs Obs
    fig, (ax1, ax2) = plt.subplots(1,2, figsize = (12,6), tight_layout=True)
    ax1.scatter(x = ml_form, y = abs_diff_nocorr, color = 'g')
    ax1.set_xlabel('ML HCHO (ppb)')
    ax1.set_ylabel('Abs. Diff. M - O (ppb)')
    ax1.set_title('Abs. Residuals vs Obs')

    #Give Absolute Residual binning summary
    df_bin_copy_nocorr = df_clean_nocorr.copy()
    df_bin_copy_nocorr['O_bin'] = pd.cut(df_bin_copy_nocorr['Observed'], bins=10)
    abs_summary_nocorr = df_bin_copy_nocorr.groupby('O_bin')['Abs. Residual'].agg(['mean','median','std'])
    print(abs_summary_nocorr)

    #Scatter Plot of Absolute Residuals vs Modeled
    ax2.scatter(x = udaq_form_nocorr, y = abs_diff_nocorr, color = 'g')
    ax2.set_xlabel('UDAQ HCHO No Correction (ppb)')
    ax2.set_ylabel('Abs. Diff. M - O (ppb)')
    ax2.set_title('Abs. Residuals vs Modeled')
    plt.savefig(merge_scripts_plots_dir + 'hawthorne_udaq_ml_' + 'scatter_formaldehyde_nocorrection_AbsResiduals.png', dpi =300)
    plt.show()

    #Scatter Plot of Relative Error vs Obs
    df_clean_nocorr['Rel. Error'] = df_clean_nocorr['Abs. Residual']/df_clean_nocorr['Observed']
    plt.figure(figsize=(12,6), tight_layout = True)
    plt.scatter(x=df_clean_nocorr['Observed'], y=df_clean_nocorr['Rel. Error'], color = 'c')
    plt.xlabel('ML HCHO (ppb)')
    plt.ylabel('Relative Error')
    plt.show()

    #Relative Error binned
    df_bin_copy_nocorr = df_clean_nocorr.copy()
    df_bin_copy_nocorr['O_bin'] = pd.cut(df_bin_copy_nocorr['Observed'], bins=10)
    print(df_bin_copy_nocorr.groupby('O_bin')['Rel. Error'].median())

    #Corrected Formaldehyde Plots
    #Scatter Plot of ML vs UDAQ
    fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize = (12,6), tight_layout=True)
    ax1.scatter(x = ml_form, y = udaq_form_corr, color = 'b')
    ax1.set_xlabel('ML HCHO (ppb)')
    ax1.set_ylabel('UDAQ HCHO Corrected (ppb)')
    ax1.set_title('Model vs Obs')

    #Plot regression line
    O_reshaped_corr = df_clean_corr['Observed'].values.reshape(-1,1)
    model_corr = LinearRegression().fit(O_reshaped_corr, df_clean_corr['Modeled'].values)
    slope_corr = model_corr.coef_[0]
    intercept_corr = model_corr.intercept_
    # Predicted M based on regression
    M_pred_corr = model_corr.predict(O_reshaped_corr)
    # Regression line
    ax1.plot(df_clean_corr['Observed'].values, M_pred_corr, color='k', linewidth=2, label=f'Regression line: y={intercept_corr:.2f}+{slope_corr:.2f}x')
    # 1:1 line
    ax1.plot([df_clean_corr['Observed'].values.min(), df_clean_corr['Observed'].values.max()], [df_clean_corr['Observed'].values.min(), df_clean_corr['Observed'].values.max()], color='c', linestyle='--', label='1:1 line')
    ax1.legend(fontsize = 10)

    #Scatter Plot of Residuals vs Obs
    ax2.scatter(x = ml_form, y = diff_corr, color = 'b')
    ax2.set_xlabel('ML HCHO (ppb)')
    ax2.set_ylabel('Diff. M - O (ppb)')
    ax2.set_title('Residuals vs Obs')

    #Give Ordinary Least Squares regression summary
    #Looking at Coefficient of Observed which is the slope. If the p-value (P>|t|) for the Observed > 0.05 then the slope is not significantly different from zero and suggests that the residuals are constant 
    ols_corr = sm.OLS(df_clean_corr['Residual'], col_ones).fit()
    print(ols_corr.summary())

    #Scatter Plot of Residuals vs Modeled
    ax3.scatter(x = udaq_form_corr, y = diff_corr, color = 'b')
    ax3.set_xlabel('UDAQ HCHO Corrected (ppb)')
    ax3.set_ylabel('Diff. M - O (ppb)')
    ax3.set_title('Residuals vs Modeled')
    plt.savefig(merge_scripts_plots_dir + 'hawthorne_udaq_ml_' + 'scatter_formaldehyde_corrected_ModelvsObs_and_Residuals.png', dpi =300)
    plt.show()

    #Scatter Plot of Absolute Residuals vs Obs
    fig, (ax1, ax2) = plt.subplots(1,2, figsize = (12,6), tight_layout=True)
    ax1.scatter(x = ml_form, y = abs_diff_corr, color = 'b')
    ax1.set_xlabel('ML HCHO (ppb)')
    ax1.set_ylabel('Abs. Diff. M - O (ppb)')
    ax1.set_title('Abs. Residuals vs Obs')

    #Give Absolute Residual binning summary
    df_bin_copy_corr = df_clean_corr.copy()
    #Give Absolute Residual binning summary
    df_bin_copy_corr['O_bin'] = pd.cut(df_bin_copy_corr['Observed'], bins=10)
    abs_summary_corr = df_bin_copy_corr.groupby('O_bin')['Abs. Residual'].agg(['mean','median','std'])
    print(abs_summary_corr)

    #Scatter Plot of Absolute Residuals vs Modeled
    ax2.scatter(x = udaq_form_corr, y = abs_diff_corr, color = 'b')
    ax2.set_xlabel('UDAQ HCHO Corrected (ppb)')
    ax2.set_ylabel('Abs. Diff. M - O (ppb)')
    ax2.set_title('Abs. Residuals vs Modeled')
    plt.savefig(merge_scripts_plots_dir + 'hawthorne_udaq_ml_' + 'scatter_formaldehyde_corrected_AbsResiduals.png', dpi =300)
    plt.show()

    #Give Absolute Residual binning summary
    df_bin_copy_corr = df_clean_corr.copy()
    df_bin_copy_corr['M_bin'] = pd.cut(df_bin_copy_corr['Modeled'], bins=10)
    abs_summary_corr = df_bin_copy_corr.groupby('M_bin')['Abs. Residual'].agg(['mean','median','std'])
    print(abs_summary_corr)

    #Scatter Plot of Relative Error vs Obs
    df_clean_corr['Rel. Error'] = df_clean_corr['Abs. Residual']/df_clean_corr['Observed']
    plt.figure(figsize=(12,6), tight_layout = True)
    plt.scatter(x=df_clean_corr['Observed'], y=df_clean_corr['Rel. Error'], color = 'c')
    plt.xlabel('ML HCHO (ppb)')
    plt.ylabel('Relative Error')
    plt.show()

    #Relative Error binned
    df_bin_copy_corr = df_clean_corr.copy()
    df_bin_copy_corr['O_bin'] = pd.cut(df_bin_copy_corr['Observed'], bins=10)
    print(df_bin_copy_corr.groupby('O_bin')['Rel. Error'].median())

def calibration_adjustment(time_interval_used):
    averaging_func_output = formaldehyde_averaging(time_interval = time_interval_used)
    df_ml_data = averaging_func_output[0]
    df_formaldehyde_averaged = averaging_func_output[1]
    udaq_form_nocorr = df_formaldehyde_averaged['H2CO_UDAQ_NoCorrection']
    udaq_form_corr = df_formaldehyde_averaged['H2CO_UDAQ_Corrected']
    ml_form = df_ml_data['HCHO_CRDS']

    diff_nocorr = udaq_form_nocorr-ml_form
    diff_corr = udaq_form_corr-ml_form
    
    abs_diff_nocorr = abs(udaq_form_nocorr-ml_form)
    abs_diff_corr = abs(udaq_form_corr-ml_form)

    df_obs_modeled_nocorr= pd.DataFrame({'Observed': ml_form, 'Modeled': udaq_form_nocorr, 'Residual':udaq_form_nocorr-ml_form, 'Abs. Residual':abs(udaq_form_nocorr-ml_form)})
    df_clean_nocorr = df_obs_modeled_nocorr.dropna()
    df_obs_modeled_corr= pd.DataFrame({'Observed': ml_form, 'Modeled': udaq_form_corr, 'Residual':udaq_form_corr-ml_form, 'Abs. Residual':abs(udaq_form_corr-ml_form)})
    df_clean_corr = df_obs_modeled_corr.dropna()
    df_bin_copy_corrected = df_clean_corr.copy()

    # Add constant term for intercept
    M_with_const = sm.add_constant(df_bin_copy_corrected['Modeled'])

    # Fit OLS regression: O = a + b*M
    ols_calibration = sm.OLS(df_bin_copy_corrected['Observed'], M_with_const).fit()

    # Get calibration coefficients
    a = ols_calibration.params['const']  # additive offset
    b = ols_calibration.params[df_bin_copy_corrected['Modeled'].name]   # multiplicative factor

    print(f"Calibration equation: O_hat = {a:.4f} + {b:.4f} * M")
    print(ols_calibration.summary())

    modeled_vals_adjusted = (-0.1498+1.0271*udaq_form_corr)

    fig, (ax1, ax2) = plt.subplots(2,1, figsize = (16,8), tight_layout=True)
    xlim_start_jul = pd.to_datetime('2024-07-15 00:00:00')
    xlim_end_jul = pd.to_datetime('2024-07-31 23:00:00')

    #ax1 is the first row of subplot, for July only
    ax1.plot(df_formaldehyde_averaged.index, df_formaldehyde_averaged['H2CO_UDAQ_Corrected'],linestyle = 'solid', color = 'b', marker = '+', label = 'UDAQ Corrected')
    ax1.plot(df_ml_data.index, df_ml_data['HCHO_CRDS'], linestyle = 'solid', color='m', marker='x',label='Mobile Lab', alpha = 0.7)
    ax1.plot(modeled_vals_adjusted.index, modeled_vals_adjusted, linestyle='solid', color = 'y', marker = '.', label = 'Calibration-adjusted UDAQ')

    #Set x ticks
    ax1.xaxis.set_major_locator(mdates.DayLocator())
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
    # Minor ticks: every 3 hours
    ax1.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21]))
    # Rotate and format tick labels
    ax1.tick_params(axis='x', which='major')
    ax1.tick_params(axis='x', which='minor', length=3, color='gray')
    #ax.grid(True, which='both')
    
    ax1.set_ylabel('Formaldehyde (ppb)')
    #ax1.set_xlabel('Date')
    ax1.margins(x=0)
    ax1.set_xlim([xlim_start_jul, xlim_end_jul])

    ax1.legend(loc = 'upper right')
    
    #ax2 is the second row of subplot, for August only
    ax2.plot(df_formaldehyde_averaged.index, df_formaldehyde_averaged['H2CO_UDAQ_Corrected'],linestyle = 'solid', color = 'b', marker = '+', label = 'UDAQ Corrected')
    ax2.plot(df_ml_data.index, df_ml_data['HCHO_CRDS'], linestyle = 'solid', color='m', marker='x',label='Mobile Lab', alpha = 0.7)
    ax2.plot(modeled_vals_adjusted.index, modeled_vals_adjusted, linestyle='solid', color = 'y', marker = '.', label = 'Calibration-adjusted UDAQ')

    #Set x ticks
    ax2.xaxis.set_major_locator(mdates.DayLocator())
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
    # Minor ticks: every 3 hours
    ax2.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21]))
    # Rotate and format tick labels
    ax2.tick_params(axis='x', which='major')
    ax2.tick_params(axis='x', which='minor', length=3, color='gray')

    ax2.set_ylabel('Formaldehyde (ppb)')
    ax2.set_xlabel('Date (UTC)')
    ax2.margins(x=0)

    ax2.set_xlim([pd.to_datetime('2024-08-01 00:00:00'), pd.to_datetime('2024-08-18 23:00:00')])
    ax2.legend(loc = 'upper right')

    #Mark midnight for every day
    midnight_vals = []
    for midnight_idx in range(0,len(df_ml_data.index),96):
        midnight_vals.append(df_ml_data.index[midnight_idx])
    for day_pos in midnight_vals:
        ax1.axvline(day_pos, color = 'black', linestyle = 'dotted')
        ax2.axvline(day_pos, color = 'black', linestyle = 'dotted')

    plt.savefig(merge_scripts_plots_dir + 'hawthorne_udaq_ml_' + 'comparison_july_aug_formaldehyde_calibration_adjusted.png', dpi =300)
    plt.show()

    #To determine if calibration adjustment aligns closer with ML data, use Root Mean Squared Error
    # Original RMSE
    rmse_original = np.sqrt(np.mean((udaq_form_corr - ml_form)**2))

    # Calibrated RMSE
    rmse_calibrated = np.sqrt(np.mean((modeled_vals_adjusted - ml_form)**2))

    print(f"Original RMSE: {rmse_original:.4f}")
    print(f"Calibrated RMSE: {rmse_calibrated:.4f}")

    #Calculate percent improvement
    rmse_percent_improvement = 100*((rmse_original - rmse_calibrated) / (rmse_original))
    print(f'Percent improvement of calibration adjusted RMSE: {rmse_percent_improvement:.4f}')
    

def fill_formaldehyde_holes_with_calibrated_adj_vals(time_interval_used, savefilename):
    averaging_func_output = formaldehyde_averaging(time_interval = time_interval_used)
    df_ml_data = averaging_func_output[0]
    df_formaldehyde_averaged = averaging_func_output[1]
    udaq_form_corr = df_formaldehyde_averaged['H2CO_UDAQ_Corrected']
    ml_form = df_ml_data['HCHO_CRDS']

    modeled_vals_adjusted = (-0.1498+1.0271*udaq_form_corr)
    print(ml_form.index)
    print(modeled_vals_adjusted.index)

    df_merged = ml_form.fillna(modeled_vals_adjusted)
    print(df_merged)

    ds_nc = xr.open_dataset('/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/CampaignData_and_Merges/R0/CSL_MobileLab_Parked/merged/rev_15min/all_CSL_MobileLab_Parked_rev15min_iWASupdated.nc')
    print(ds_nc['HCHO_CRDS'].values)
    # Remove old variable completely
    ds_nc = ds_nc.drop_vars('HCHO_CRDS')

    # Assign as a DataArray with explicit dims and coords
    da_new = xr.DataArray(
        data=df_merged,
        dims=["time_UTC"],           # specify the correct dimension
        coords={"time_UTC": ds_nc.time_UTC},  # assign coordinate explicitly
        name='HCHO_CRDS'
    )
    ds_nc['HCHO_CRDS'] = da_new

    save_ncfilename = savefilename + '.nc'
    full_savepath = merged_data_dir_15min + save_ncfilename
    ds_nc.to_netcdf(full_savepath, format = 'NETCDF4',mode='w')
    print('Saved netCDF file with updated Formaldehyde measurements to: ', full_savepath)

    ds_new = xr.open_dataset(full_savepath)
    print(ds_new['HCHO_CRDS'].values)




## CALL FUNCTIONS HERE
# iwas_convert(
#     savefilename = 'all_CSL_MobileLab_Parked_rev15min_iWASupdated'
# )
# formaldehyde_averaging(
#     time_interval = 15*60
# )

formaldehyde_correction_compare(
    time_interval_used = 15*60

)

bias_and_error_analysis(
    time_interval_used = 15*60
)

calibration_adjustment(
    time_interval_used =  15*60
)
# fill_formaldehyde_holes_with_calibrated_adj_vals(
#     time_interval_used = 15*60,
#     savefilename = 'all_CSL_MobileLab_Parked_rev15min_iWASupdated_formaldehydeupdated'
# )