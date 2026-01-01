#Compare iWAS and PTR measurements for several species
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.lines as mlines
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerTuple
from matplotlib.patches import Patch

#region: Settings
#region: filepaths
#Filepaths for loading
dirpath = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/'
#Filepaths for saving
compare_udaq_ml_savepath = dirpath + 'USOS_Campaign_analysis/plots/'
#endregion

nc_filepath = dirpath + 'CampaignData_and_Merges/R0/CSL_MobileLab_Parked/merged/rev_15min/all_CSL_MobileLab_Parked_rev15min_iWASupdated_formaldehydeupdated.nc'    
nc_load = xr.open_dataset(nc_filepath)
df_ml_data = nc_load[['Alpha_Pinene_WAS','Beta_Pinene_WAS', 'Monoterpenes_PTR', "Isoprene_WAS", 'Isoprene_PTR', 'time_local']].to_dataframe()
df_ml_data = df_ml_data.set_index(['time_local'])
display(df_ml_data)

df_ml_data['pinene_sum'] = df_ml_data['Alpha_Pinene_WAS'] + df_ml_data['Beta_Pinene_WAS']
display(df_ml_data)

f0am_filled_path = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/CampaignData_and_Merges/R0/CSL_MobileLab_Parked/' + 'F0AM_filled/' + '20240805_20240807_15min_CSL_mobile_lab_parked_with_interp_nell_match_with_formaldehyde'+ '.csv'
df_f0am_data = pd.read_csv(f0am_filled_path, index_col='time_local', parse_dates=True)

new_start_time_vocs = pd.Timestamp('2024-07-15 00:00:00')
new_end_time_vocs = pd.Timestamp('2024-08-18 23:00:00')
#Create a new datetime index from new_start to the end of existing index with same frequency
new_index_vocs = pd.date_range(start=new_start_time_vocs, end=new_end_time_vocs, freq='1h')
udaq_voc_file = dirpath + 'Hawthorne_data/data/hawthorne_udaq_all_vocs_hourly_timezone_carbon_number_updated.csv'
df_udaq_voc_data = pd.read_csv(udaq_voc_file, index_col='time_local', parse_dates=True)
df_udaq_voc_data = df_udaq_voc_data.reindex(new_index_vocs)


def plot_species_compare():
    fig, (ax1, ax2) = plt.subplots(2,1, figsize = (16,8), tight_layout=True)
    xlim_start_jul = pd.to_datetime('2024-07-15 00:00:00')
    xlim_end_jul = pd.to_datetime('2024-07-31 23:00:00')

    #ax1 is the first row of subplot, for July only
    ax1.plot(df_ml_data.index, df_ml_data['Isoprene_WAS'],linestyle = 'solid', color = 'g', marker = '+', label = 'WAS')
    ax1.plot(df_ml_data.index, df_ml_data['Isoprene_PTR'],linestyle = 'solid', color = 'b', marker = '+', label = 'PTR')
    ax1.plot(df_f0am_data.index, df_f0am_data['Isoprene_WAS'], linestyle = 'solid', color = 'y', marker = '*', label = 'F0AM filled')
    ax1.plot(df_udaq_voc_data.index, df_udaq_voc_data['Isoprene'],linestyle = 'solid', color = 'm', marker = '.', label = 'UDAQ')
    
    #Set x ticks
    ax1.xaxis.set_major_locator(mdates.DayLocator())
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
    # Minor ticks: every 3 hours
    ax1.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21]))
    # Rotate and format tick labels
    ax1.tick_params(axis='x', which='major')
    ax1.tick_params(axis='x', which='minor', length=3, color='gray')
    #ax.grid(True, which='both')

    ax1.set_ylabel('Isoprene (ppb)')
    #ax1.set_xlabel('Date')
    ax1.margins(x=0)
    ax1.set_xlim([xlim_start_jul, xlim_end_jul])

    ax1.legend(loc = 'upper right')

    #ax2 is the second row of subplot, for August only
    ax2.plot(df_ml_data.index, df_ml_data['Isoprene_WAS'],linestyle = 'solid', color = 'g', marker = '*', label = 'WAS')
    ax2.plot(df_ml_data.index, df_ml_data['Isoprene_PTR'],linestyle = 'solid', color = 'b', marker = '+', label = 'PTR')
    ax2.plot(df_f0am_data.index, df_f0am_data['Isoprene_WAS'], linestyle = 'solid', color = 'y', marker = '.', label = 'F0AM filled')
    ax2.plot(df_udaq_voc_data.index, df_udaq_voc_data['Isoprene'],linestyle = 'solid', color = 'm', marker = '.', label = 'UDAQ', alpha = 0.6)

    #Set x ticks
    ax2.xaxis.set_major_locator(mdates.DayLocator())
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
    # Minor ticks: every 3 hours
    ax2.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21]))
    # Rotate and format tick labels
    ax2.tick_params(axis='x', which='major')
    ax2.tick_params(axis='x', which='minor', length=3, color='gray')

    ax2.set_ylabel('Isoprene (ppb)')
    ax2.set_xlabel('Date (MDT)')
    ax2.margins(x=0)

    ax2.set_xlim([pd.to_datetime('2024-08-01 00:00:00'), pd.to_datetime('2024-08-18 23:00:00')])
    ax2.legend(loc = 'upper right')

    #Mark midnight for every day
    midnight_vals = []
    for midnight_idx in range(24,len(df_ml_data.index),96):
        midnight_vals.append(df_ml_data.index[midnight_idx])
    for day_pos in midnight_vals:
        ax1.axvline(day_pos, color = 'black', linestyle = 'dotted')
        ax2.axvline(day_pos, color = 'black', linestyle = 'dotted')

    #plt.savefig(merge_scripts_plots_dir + 'hawthorne_udaq_ml_' + 'comparison_july_aug_formaldehyde.png', dpi =300)
    plt.show()

def plot_species_compare_short():
    fig, (ax1) = plt.subplots(1,1, figsize = (16,8), tight_layout=True)
    xlim_start_jul = pd.to_datetime('2024-08-05 00:00:00')
    xlim_end_jul = pd.to_datetime('2024-08-07 23:00:00')

    #ax1 is the first row of subplot, for July only
    ax1.plot(df_ml_data.index, df_ml_data['Isoprene_WAS'],linestyle = 'solid', color = 'g', marker = '*', label = 'WAS')
    ax1.plot(df_ml_data.index, df_ml_data['Isoprene_PTR'],linestyle = 'solid', color = 'b', marker = '+', label = 'PTR')
    ax1.plot(df_f0am_data.index, df_f0am_data['Isoprene_WAS'], linestyle = 'solid', color = 'y', marker = '.', label = 'F0AM filled', alpha=0.3)
    ax1.plot(df_udaq_voc_data.index, df_udaq_voc_data['Isoprene'],linestyle = 'solid', color = 'm', marker = '.', label = 'UDAQ')
    #Set x ticks
    ax1.xaxis.set_major_locator(mdates.DayLocator())
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
    # Minor ticks: every 3 hours
    ax1.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21]))
    # Rotate and format tick labels
    ax1.tick_params(axis='x', which='major')
    ax1.tick_params(axis='x', which='minor', length=3, color='gray')
    #ax.grid(True, which='both')

    ax1.set_ylabel('Isoprene (ppb)')
    #ax1.set_xlabel('Date')
    ax1.margins(x=0)
    ax1.set_xlim([xlim_start_jul, xlim_end_jul])

    ax1.legend(loc = 'upper right')

    # #Mark midnight for every day
    # midnight_vals = []
    # for midnight_idx in range(24,len(df_ml_data.index),96):
    #     midnight_vals.append(df_ml_data.index[midnight_idx])
    # for day_pos in midnight_vals:
    #     ax1.axvline(day_pos, color = 'black', linestyle = 'dotted')

    #plt.savefig(merge_scripts_plots_dir + 'hawthorne_udaq_ml_' + 'comparison_july_aug_formaldehyde.png', dpi =300)
    plt.show()

def jno2_ratio():
    f0am_filled_path = dirpath+'CampaignData_and_Merges/R0/CSL_MobileLab_Parked/F0AM_filled/20240805_20240807_15min_CSL_mobile_lab_parked_with_interp_nell_match_with_formaldehyde.csv'
    df_f0am_filled = pd.read_csv(f0am_filled_path,index_col = 'time_local', parse_dates=True)
    
    start_time = pd.to_datetime('2024-08-05 00:00:00')
    end_time = pd.to_datetime('2024-08-05 23:45:00')
    df_f0am_filled = df_f0am_filled.loc[start_time:end_time]
    obs_div_calc = df_f0am_filled['jNO2_meas']/df_f0am_filled['jNO2']
    print(obs_div_calc.values)
    print(df_f0am_filled['jNO2_ratio'].values)
    print(df_f0am_filled['jNO2_meas'].values)
    print(df_f0am_filled['jNO2'].values)

#CALL FUNCTIONS
plot_species_compare()
plot_species_compare_short()
#jno2_ratio()