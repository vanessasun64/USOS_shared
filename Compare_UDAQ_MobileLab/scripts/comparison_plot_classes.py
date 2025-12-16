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
compare_udaq_ml_savepath = dirpath + 'Compare_UDAQ_MobileLab/plots/'
#endregion


#region: resave ozone concentration to get an index that has the ozone overlap for both UDAQ and ML
#Also calculate MDA8 ozone and associated variables

mobilelab_hourly_o3_filepath = dirpath + 'USOS_Campaign_analysis/resaved_data/ozone_hourly.csv'
df_ml_o3_data = pd.read_csv(mobilelab_hourly_o3_filepath,index_col = 'time_local', parse_dates=True)
udaq_o3_filepath = dirpath + 'Hawthorne_data/data/hawthorne_udaq_o3_2024_timezone_updated.csv'
df_udaq_o3_data = pd.read_csv(udaq_o3_filepath, index_col = 'time_local', parse_dates=True)

def mda8_calculate(save_filename):
    mda8_index_start = pd.to_datetime('2024-07-16 12:00:00')
    mda8_index_end = pd.to_datetime('2024-08-18 18:00:00')
    mda8_new_index = pd.date_range(start=mda8_index_start, end=mda8_index_end, freq='1h')
    mda8_udaq_data = df_udaq_o3_data.reindex(mda8_new_index)
    mda8_ml_data = df_ml_o3_data.reindex(mda8_new_index)
    #MDA8 Calculation

    # Calculate the 8-hour rolling average
    rolling_avg_name_UDAQ = '8hr_rolling_avg_' + 'UDAQ O3'
    rolling_avg_name_ML = '8hr_rolling_avg_' + 'ML O3'
    mda8_udaq_data[rolling_avg_name_UDAQ] = mda8_udaq_data['UDAQ O3'].rolling(window=8, min_periods=6).mean()
    mda8_ml_data[rolling_avg_name_ML] = mda8_ml_data['Mobile Lab O3'].rolling(window=8, min_periods=6).mean()
    #rolling average only works for 8 hours with no gaps; at least 6 hours

    # Calculate the maximum 8-hour average for each day
    # Resample to daily frequency, compute the daily max of the rolling averages, drop NA values
    daily_max_8hr_avg_ozone_UDAQ = mda8_udaq_data[rolling_avg_name_UDAQ].resample('D').max().dropna()
    daily_max_8hr_avg_ozone_ML = mda8_ml_data[rolling_avg_name_ML].resample('D').max().dropna()

    df_combined_daily8hrmax_ozone = pd.concat([daily_max_8hr_avg_ozone_UDAQ, daily_max_8hr_avg_ozone_ML], axis = 1)
    df_combined_daily8hrmax_ozone.index.name = "date"
    daily8hrmax_ozone_savepath = dirpath + 'Compare_UDAQ_MobileLab/resaved_data/daily8hrmax_ozone_udaq_ml.csv'
    df_combined_daily8hrmax_ozone.to_csv(daily8hrmax_ozone_savepath)

    # Map the daily maximum back to the original dataframe
    # Create a new temporary column with the daily max 8-hour average for each timestamp
    mda8_ozone_name_UDAQ = 'MDA8_O3_' + 'UDAQ O3'
    mda8_ozone_name_ML = 'MDA8_O3_' + 'ML O3'
    mda8_udaq_data[mda8_ozone_name_UDAQ] = mda8_udaq_data.index.floor('D').map(daily_max_8hr_avg_ozone_UDAQ)
    mda8_ml_data[mda8_ozone_name_ML] = mda8_ml_data.index.floor('D').map(daily_max_8hr_avg_ozone_ML)

    # # Select daytime values only where MD8A > 70
    exceedance_name_UDAQ = 'Exceedance_day_' + 'UDAQ O3'
    exceedance_name_ML = 'Exceedance_day_' + 'ML O3'
    mda8_udaq_data[exceedance_name_UDAQ] = mda8_udaq_data[mda8_ozone_name_UDAQ] >= 70
    mda8_ml_data[exceedance_name_ML] = mda8_ml_data[mda8_ozone_name_ML] >= 70

    df_combined_mda8_ozone = mda8_udaq_data.join(mda8_ml_data)
    df_combined_mda8_ozone.index.name = 'time_local'
    savepath = dirpath + 'Compare_UDAQ_MobileLab/resaved_data/' + save_filename + '.csv'
    df_combined_mda8_ozone.to_csv(savepath)
    print('Saved to:' + savepath)
#mda8_calculate('mda8_udaq_ml_combined_df')

#endregion

# region: file loading
#MDA8 Ozone - has index from 2024-07-16 12:00:00 to 2024-08-18 18:00:00
mda8_ozone_filepath = dirpath + 'Compare_UDAQ_MobileLab/resaved_data/mda8_udaq_ml_combined_df.csv'
df_mda8_o3_data = pd.read_csv(mda8_ozone_filepath, index_col = 'time_local', parse_dates=True)

daily8hrmax_ozone_filepath = dirpath + 'Compare_UDAQ_MobileLab/resaved_data/daily8hrmax_ozone_udaq_ml.csv'
df_daily8hrmax_ozone = pd.read_csv(daily8hrmax_ozone_filepath, index_col = 'date', parse_dates=True)

#VOCs 
ml_voc_file = dirpath + 'Compare_UDAQ_MobileLab/resaved_data/ml_hourly_voc_overlap_reindexed.csv'
df_ml_voc_data = pd.read_csv(ml_voc_file, index_col='time_local', parse_dates=True)
#Reindex to hourly to ensure a continuous time index: 2024-07-15 03:00:00 and 2024-08-12 23:00:00
new_start_time_vocs = pd.Timestamp('2024-07-15 03:00:00')
new_end_time_vocs = pd.Timestamp('2024-08-12 23:00:00')
#Create a new datetime index from new_start to the end of existing index with same frequency
new_index_vocs = pd.date_range(start=new_start_time_vocs, end=new_end_time_vocs, freq='1h')
#Reindex the dataframe to include new rows
df_ml_voc_data = df_ml_voc_data.reindex(new_index_vocs)

udaq_voc_file = dirpath + 'Hawthorne_data/data/hawthorne_udaq_all_vocs_hourly_timezone_carbon_number_updated.csv'
df_udaq_voc_data = pd.read_csv(udaq_voc_file, index_col='time_local', parse_dates=True)
df_udaq_voc_data = df_udaq_voc_data.reindex(new_index_vocs)
# endregion

#region: Global variables, used by multiple plots / repeatedly
udaq_o3_original = df_udaq_o3_data['UDAQ O3']
udaq_o3_index_original = df_udaq_o3_data.index
ml_o3_original = df_ml_o3_data['Mobile Lab O3']
ml_o3_index_original = df_ml_o3_data.index

daily_max_8hr_avg_ozone_udaq = df_daily8hrmax_ozone['8hr_rolling_avg_UDAQ O3']
daily_max_8hr_avg_ozone_ml = df_daily8hrmax_ozone['8hr_rolling_avg_ML O3']

df_mda8_o3_data['hour']=df_mda8_o3_data.index.hour
df_udaq_voc_data['hour']=df_udaq_voc_data.index.hour
df_ml_voc_data['hour']=df_ml_voc_data.index.hour

#hrly mean ozone
mean_hrly_udaq_ozone=df_mda8_o3_data.groupby('hour')['UDAQ O3'].apply(np.nanmean)
mean_hrly_ml_ozone=df_mda8_o3_data.groupby('hour')['Mobile Lab O3'].apply(np.nanmean)
#hrly standard deviation ozone
std_hrly_udaq_ozone=df_mda8_o3_data.groupby('hour')['UDAQ O3'].apply(np.nanstd)
std_hrly_ml_ozone=df_mda8_o3_data.groupby('hour')['Mobile Lab O3'].apply(np.nanstd)

#hrly median ozone
median_hrly_udaq_ozone=df_mda8_o3_data.groupby('hour')['UDAQ O3'].apply(np.nanmedian)
median_hrly_ml_ozone=df_mda8_o3_data.groupby('hour')['Mobile Lab O3'].apply(np.nanmedian)
#hrly 25th and 75th quartile
udaq_percentile_25_ozone=df_mda8_o3_data.groupby('hour')['UDAQ O3'].apply(lambda x: np.nanpercentile(x, 25))
udaq_percentile_75_ozone=df_mda8_o3_data.groupby('hour')['UDAQ O3'].apply(lambda x: np.nanpercentile(x, 75))
ml_percentile_25_ozone=df_mda8_o3_data.groupby('hour')['Mobile Lab O3'].apply(lambda x: np.nanpercentile(x, 25))
ml_percentile_75_ozone=df_mda8_o3_data.groupby('hour')['Mobile Lab O3'].apply(lambda x: np.nanpercentile(x, 75))

hour_range = np.arange(0,24,1)

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

#endregion

#region: Notes
#Any stats use only the overlap in index between ml and udaq
#endregion


class time_series_udaq_ml_comparison:
    def __init__(self, species_name, var_name_modification=None):
        self.species_name = species_name
        self.var_name_modification = var_name_modification
    def species_full_time_series_comparison_no_uncertainties(self, legend_loc, SavePlotSpeciesName):
        fig, (ax1, ax2) = plt.subplots(2,1, figsize = (16,8), tight_layout=True)

        if self.species_name == 'Ozone':
            species_udaq_index = udaq_o3_index_original
            species_udaq_var = udaq_o3_original
            species_ml_index = ml_o3_index_original
            species_ml_var = ml_o3_original
            species_yticks = np.arange(10,120,10)

            species_name_var = 'Ozone'
            species_unit = '(ppb)'
            xlim_start_jul = pd.to_datetime('2024-07-16 00:00:00')
            xlim_end_jul = pd.to_datetime('2024-07-31 23:00:00')


            #Add 70 ppb horizontal line
            ax1.hlines(y=70, xmin = species_udaq_index[0], xmax = species_udaq_index[len(species_udaq_index)-1], linestyle = 'dashed', color = 'r')
            ax2.hlines(y=70, xmin = species_udaq_index[0], xmax = species_udaq_index[len(species_udaq_index)-1], linestyle = 'dashed', color = 'r')

        else:
            species_udaq_index = df_udaq_voc_data.index
            species_udaq_var = df_udaq_voc_data[self.species_name]
            species_ml_index = df_ml_voc_data.index
            species_ml_var = df_ml_voc_data[self.species_name]
            if self.var_name_modification is not None:
                species_name_var = self.var_name_modification
            else:
                species_name_var = self.species_name
            species_unit = '(ppb)' #will be modified in future if we end up using any non-ppb units
            xlim_start_jul = pd.to_datetime('2024-07-15 00:00:00')
            xlim_end_jul = pd.to_datetime('2024-07-31 23:00:00')

        #ax1 is the first row of subplot, for July only
        ax1.plot(species_udaq_index, species_udaq_var, linestyle = 'solid', color = 'g', marker = '+', label = 'UDAQ')
        ax1.plot(species_ml_index, species_ml_var, linestyle = 'solid', color='m', marker='x',label='Mobile Lab', alpha = 0.7)
            
        #Mark midnight for every day
        midnight_vals = []
        for midnight_idx in range(0,len(species_udaq_index),24):
            midnight_vals.append(species_udaq_index[midnight_idx])
        for day_pos in midnight_vals:
            ax1.axvline(day_pos, color = 'black', linestyle = 'dotted', alpha = 0.7)

        #Set x ticks
        ax1.xaxis.set_major_locator(mdates.DayLocator())
        ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
        # Minor ticks: every 3 hours
        ax1.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21]))
        # Rotate and format tick labels
        ax1.tick_params(axis='x', which='major')
        ax1.tick_params(axis='x', which='minor', length=3, color='gray')
        #ax.grid(True, which='both')
        if self.species_name == 'Ozone':
            ax1.set_yticks(species_yticks)
        else:
            pass
        ax1.set_ylabel(species_name_var + ' ' + species_unit)
        #ax1.set_xlabel('Date')
        ax1.margins(x=0)
        ax1.set_xlim([xlim_start_jul, xlim_end_jul])

        ax1.legend(loc = legend_loc)
        
        #ax2 is the second row of subplot, for August only
        ax2.plot(species_udaq_index, species_udaq_var, linestyle = 'solid', color = 'g', marker = '+', label = 'UDAQ')
        ax2.plot(species_ml_index, species_ml_var, linestyle = 'solid', color='m', marker='x', label='Mobile Lab', alpha = 0.7)
        
        #Mark midnight for every day
        midnight_vals = []
        for midnight_idx in range(0,len(species_udaq_index),24):
            midnight_vals.append(species_udaq_index[midnight_idx])
        for day_pos in midnight_vals:
            ax2.axvline(day_pos, color = 'black', linestyle = 'dotted', alpha = 0.7)
        
        #Set x ticks
        ax2.xaxis.set_major_locator(mdates.DayLocator())
        ax2.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
        # Minor ticks: every 3 hours
        ax2.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21]))
        # Rotate and format tick labels
        ax2.tick_params(axis='x', which='major')
        ax2.tick_params(axis='x', which='minor', length=3, color='gray')

        if self.species_name == 'Ozone':
            ax2.set_yticks(species_yticks)
        else:
            pass

        ax2.set_ylabel(species_name_var + ' ' + species_unit)
        ax2.set_xlabel('Date')
        ax2.margins(x=0)

        ax2.set_xlim([pd.to_datetime('2024-08-01 00:00:00'), pd.to_datetime('2024-08-18 23:00:00')])
        ax2.legend(loc = legend_loc)
        plt.savefig(compare_udaq_ml_savepath + 'hawthorne_udaq_ml_' + 'comparison_july_aug_no_uncertainties.png', dpi =300)
        plt.show()
    def species_full_time_series_comparison_with_uncertainties(self, legend_loc, SavePlotSpeciesName,
        InstrumentUncertaintyUDAQ=None,
        InstrumentUncertaintyML=None
        ):

        #Make sure to input a float as function arg for instrument uncertainty if it exists
        #Ozone set to 1.5 ppb uncertainty
        #If uncertainty is unknown, can be left off

        fig, (ax1, ax2) = plt.subplots(2,1, figsize = (16,8), tight_layout=True)

        if self.species_name == 'Ozone':
            species_udaq_index = udaq_o3_index_original
            species_udaq_var = udaq_o3_original
            species_ml_index = ml_o3_index_original
            species_ml_var = ml_o3_original

            #Set uncertainty for 2B Tech
            udaq_instr_uncertainty = 1.5
            ml_instr_uncertainty = 1.5

            species_yticks = np.arange(10,120,10)
            species_name_var = 'Ozone'
            species_unit = '(ppb)'
            xlim_start_jul = pd.to_datetime('2024-07-16 00:00:00')
            xlim_end_jul = pd.to_datetime('2024-07-31 23:00:00')

            #Add 70 ppb horizontal line
            ax1.hlines(y=70, xmin = species_udaq_index[0], xmax = species_udaq_index[len(species_udaq_index)-1], linestyle = 'dashed', color = 'r')
            ax2.hlines(y=70, xmin = species_udaq_index[0], xmax = species_udaq_index[len(species_udaq_index)-1], linestyle = 'dashed', color = 'r')
            

        else:
            species_udaq_index = df_udaq_voc_data.index
            species_udaq_var = df_udaq_voc_data[self.species_name]
            species_ml_index = df_ml_voc_data.index
            species_ml_var = df_ml_voc_data[self.species_name]

            udaq_instr_uncertainty = InstrumentUncertaintyUDAQ
            ml_instr_uncertainty = InstrumentUncertaintyML

            species_yticks = None
            if self.var_name_modification is not None:
                species_name_var = self.var_name_modification
            else:
                species_name_var = self.species_name
            species_unit = '(ppb)' #will be modified in future if we end up using any non-ppb units
            xlim_start_jul = pd.to_datetime('2024-07-15 00:00:00')
            xlim_end_jul = pd.to_datetime('2024-07-31 23:00:00')

        #ax1 is the first row of subplot, for July only
        ax1.plot(species_udaq_index, species_udaq_var, linestyle = 'solid', color = 'g', marker = '+', label = 'UDAQ')
        ax1.plot(species_ml_index, species_ml_var, linestyle = 'solid', color='m', marker='x',label='Mobile Lab', alpha = 0.7)
        
        if udaq_instr_uncertainty is not None and ml_instr_uncertainty is not None:
            #plot uncertainty in shaded
            ax1.fill_between(species_udaq_index, y1 = species_udaq_var - udaq_instr_uncertainty, y2 = species_udaq_var + udaq_instr_uncertainty,  color = 'g', alpha = 0.4)
            ax1.fill_between(species_ml_index, y1 = species_ml_var - ml_instr_uncertainty, y2 = species_ml_var + ml_instr_uncertainty, color = 'm', alpha = 0.4)


        #Mark midnight for every day
        midnight_vals = []
        for midnight_idx in range(0,len(species_udaq_index),24):
            midnight_vals.append(species_udaq_index[midnight_idx])
        for day_pos in midnight_vals:
            ax1.axvline(day_pos, color = 'black', linestyle = 'dotted', alpha = 0.7)

        #Set x ticks
        ax1.xaxis.set_major_locator(mdates.DayLocator())
        ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
        # Minor ticks: every 3 hours
        ax1.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21]))
        # Rotate and format tick labels
        ax1.tick_params(axis='x', which='major')
        ax1.tick_params(axis='x', which='minor', length=3, color='gray')
        #ax.grid(True, which='both')

        if self.species_name == 'Ozone':
            ax1.set_yticks(species_yticks)
        else:
            pass

        ax1.set_xlim([xlim_start_jul, xlim_end_jul])

        ax1.set_ylabel(species_name_var + ' ' + species_unit)
        #ax1.set_xlabel('Date')
        ax1.margins(x=0)
        
        ax1.legend(loc = legend_loc)
        
        #ax2 is the second row of subplot, for August only

        ax2.plot(species_udaq_index, species_udaq_var, linestyle = 'solid', color = 'g', marker = '+', label = 'UDAQ')
        ax2.plot(species_ml_index, species_ml_var, linestyle = 'solid', color='m', marker='x', label='Mobile Lab', alpha = 0.7)
        
        if udaq_instr_uncertainty is not None and ml_instr_uncertainty is not None:
            #plot uncertainty in shaded
            ax2.fill_between(species_udaq_index, y1 = species_udaq_var - udaq_instr_uncertainty, y2 = species_udaq_var + udaq_instr_uncertainty,  color = 'g', alpha = 0.4)
            ax2.fill_between(species_ml_index, y1 = species_ml_var - ml_instr_uncertainty, y2 = species_ml_var + ml_instr_uncertainty, color = 'm', alpha = 0.4)

        #Mark midnight for every day
        midnight_vals = []
        for midnight_idx in range(0,len(species_udaq_index),24):
            midnight_vals.append(species_udaq_index[midnight_idx])
        for day_pos in midnight_vals:
            ax2.axvline(day_pos, color = 'black', linestyle = 'dotted', alpha = 0.7)
        
        #Set x ticks
        ax2.xaxis.set_major_locator(mdates.DayLocator())
        ax2.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
        # Minor ticks: every 3 hours
        ax2.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21]))
        # Rotate and format tick labels
        ax2.tick_params(axis='x', which='major')
        ax2.tick_params(axis='x', which='minor', length=3, color='gray')
      
        if self.species_name == 'Ozone':
            ax2.set_yticks(species_yticks)
        else:
            pass

        ax2.set_xlim([pd.to_datetime('2024-08-01 00:00:00'), pd.to_datetime('2024-08-18 23:00:00')])

        ax2.set_ylabel(species_name_var + ' ' + species_unit)
        ax2.set_xlabel('Date')
        ax2.margins(x=0)

        ax2.legend(loc = legend_loc)
        plt.savefig(compare_udaq_ml_savepath + 'hawthorne_udaq_ml_' + 'comparison_july_aug_with_uncertainties.png', dpi =300)
        plt.show()
    def species_full_time_series_comparison_with_inset_and_uncertainties(self, legend_loc, SavePlotSpeciesName,
        InstrumentUncertaintyUDAQ=None,
        InstrumentUncertaintyML=None
        ):

        #Make sure to input a float as function arg for instrument uncertainty if it exists
        #Ozone set to 1.5 ppb uncertainty
        #If uncertainty is unknown, can be left off

        fig, (ax1, ax2) = plt.subplots(2,1, figsize = (16,8), tight_layout=True)

        if self.species_name == 'Ozone':
            species_udaq_index = udaq_o3_index_original
            species_udaq_var = udaq_o3_original
            species_ml_index = ml_o3_index_original
            species_ml_var = ml_o3_original

            #Set uncertainty for 2B Tech
            udaq_instr_uncertainty = 1.5
            ml_instr_uncertainty = 1.5

            species_yticks = np.arange(10,130,10)
            species_name_var = 'Ozone'
            species_unit = '(ppb)'
            xlim_start_jul = pd.to_datetime('2024-07-16 00:00:00')
            xlim_end_jul = pd.to_datetime('2024-07-31 23:00:00')

            #Add 70 ppb horizontal line
            ax1.hlines(y=70, xmin = species_udaq_index[0], xmax = species_udaq_index[len(species_udaq_index)-1], linestyle = 'dashed', color = 'r')
            ax2.hlines(y=70, xmin = species_udaq_index[0], xmax = species_udaq_index[len(species_udaq_index)-1], linestyle = 'dashed', color = 'r')
            

        else:
            species_udaq_index = df_udaq_voc_data.index
            species_udaq_var = df_udaq_voc_data[self.species_name]
            species_ml_index = df_ml_voc_data.index
            species_ml_var = df_ml_voc_data[self.species_name]

            udaq_instr_uncertainty = InstrumentUncertaintyUDAQ
            ml_instr_uncertainty = InstrumentUncertaintyML

            species_yticks = None
            if self.var_name_modification is not None:
                species_name_var = self.var_name_modification
            else:
                species_name_var = self.species_name
            species_unit = '(ppb)' #will be modified in future if we end up using any non-ppb units
            xlim_start_jul = pd.to_datetime('2024-07-15 00:00:00')
            xlim_end_jul = pd.to_datetime('2024-07-31 23:00:00')

        #ax1 is the first row of subplot, for July only
        ax1.plot(species_udaq_index, species_udaq_var, linestyle = 'solid', color = 'g', marker = '+', label = 'UDAQ')
        ax1.plot(species_ml_index, species_ml_var, linestyle = 'solid', color='m', marker='x',label='Mobile Lab', alpha = 0.7)
        
        if udaq_instr_uncertainty is not None and ml_instr_uncertainty is not None:
            #plot uncertainty in shaded
            ax1.fill_between(species_udaq_index, y1 = species_udaq_var - udaq_instr_uncertainty, y2 = species_udaq_var + udaq_instr_uncertainty,  color = 'g', alpha = 0.4)
            ax1.fill_between(species_ml_index, y1 = species_ml_var - ml_instr_uncertainty, y2 = species_ml_var + ml_instr_uncertainty, color = 'm', alpha = 0.4)

            axin = ax1.inset_axes([0.6, 0.75, 0.2, 0.2])
            axin.set_xlim([pd.to_datetime('2024-07-24 09:00:00'), pd.to_datetime('2024-07-24 20:00:00')])
            axin.set_ylim(70, 90)
            axin.plot(species_udaq_index, species_udaq_var, linestyle = 'solid', color = 'g', marker = '.', label = 'UDAQ')
            axin.plot(species_ml_index, species_ml_var, linestyle = 'solid', color='m', marker='x',label='Mobile Lab', alpha = 0.7)
            axin.fill_between(species_udaq_index, y1 = species_udaq_var - udaq_instr_uncertainty, y2 = species_udaq_var + udaq_instr_uncertainty,  color = 'g', alpha = 0.3)
            axin.fill_between(species_ml_index, y1 = species_ml_var - ml_instr_uncertainty, y2 = species_ml_var + ml_instr_uncertainty, color = 'm', alpha = 0.3)
            axin.tick_params(axis='x', labelsize=10)
            axin.tick_params(axis='y', labelsize=10)
            axin.xaxis.set_major_locator(mdates.HourLocator())
            axin.xaxis.set_major_formatter(mdates.DateFormatter('%H'))

            ax1.indicate_inset_zoom(axin)

        #Mark midnight for every day
        midnight_vals = []
        for midnight_idx in range(0,len(species_udaq_index),24):
            midnight_vals.append(species_udaq_index[midnight_idx])
        for day_pos in midnight_vals:
            ax1.axvline(day_pos, color = 'black', linestyle = 'dotted', alpha = 0.7)

        #Set x ticks
        ax1.xaxis.set_major_locator(mdates.DayLocator())
        ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
        # Minor ticks: every 3 hours
        ax1.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21]))
        # Rotate and format tick labels
        ax1.tick_params(axis='x', which='major')
        ax1.tick_params(axis='x', which='minor', length=3, color='gray')
        #ax.grid(True, which='both')

        if self.species_name == 'Ozone':
            ax1.set_yticks(species_yticks)
        else:
            pass

        ax1.set_xlim([xlim_start_jul, xlim_end_jul])

        ax1.set_ylabel(species_name_var + ' ' + species_unit)
        #ax1.set_xlabel('Date')
        ax1.margins(x=0)
        
        ax1.legend(loc = legend_loc)
        
        #ax2 is the second row of subplot, for August only

        ax2.plot(species_udaq_index, species_udaq_var, linestyle = 'solid', color = 'g', marker = '+', label = 'UDAQ')
        ax2.plot(species_ml_index, species_ml_var, linestyle = 'solid', color='m', marker='x', label='Mobile Lab', alpha = 0.7)
        
        if udaq_instr_uncertainty is not None and ml_instr_uncertainty is not None:
            #plot uncertainty in shaded
            ax2.fill_between(species_udaq_index, y1 = species_udaq_var - udaq_instr_uncertainty, y2 = species_udaq_var + udaq_instr_uncertainty,  color = 'g', alpha = 0.4)
            ax2.fill_between(species_ml_index, y1 = species_ml_var - ml_instr_uncertainty, y2 = species_ml_var + ml_instr_uncertainty, color = 'm', alpha = 0.4)

        #Mark midnight for every day
        midnight_vals = []
        for midnight_idx in range(0,len(species_udaq_index),24):
            midnight_vals.append(species_udaq_index[midnight_idx])
        for day_pos in midnight_vals:
            ax2.axvline(day_pos, color = 'black', linestyle = 'dotted', alpha = 0.7)
        
        #Set x ticks
        ax2.xaxis.set_major_locator(mdates.DayLocator())
        ax2.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
        # Minor ticks: every 3 hours
        ax2.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21]))
        # Rotate and format tick labels
        ax2.tick_params(axis='x', which='major')
        ax2.tick_params(axis='x', which='minor', length=3, color='gray')
      
        if self.species_name == 'Ozone':
            ax2.set_yticks(species_yticks)
        else:
            pass

        ax2.set_xlim([pd.to_datetime('2024-08-01 00:00:00'), pd.to_datetime('2024-08-18 23:00:00')])

        ax2.set_ylabel(species_name_var + ' ' + species_unit)
        ax2.set_xlabel('Date')
        ax2.margins(x=0)

        ax2.legend(loc = legend_loc)
        plt.savefig(compare_udaq_ml_savepath + 'hawthorne_udaq_ml_' + 'comparison_july_aug_with_inset_and_uncertainties.png', dpi =300)
        plt.show()

class mda8_ozone:
    def __init__(self):
        self.data = []
    def mda8_ozone_one_per_day(self, legend_loc, SavePlotSpeciesName):
        fig, ax = plt.subplots(2,1, figsize = (24,8), tight_layout = True)
        threshold = 70
        ax[0].margins(x=0)
        ax[0].plot(daily_max_8hr_avg_ozone_udaq.index, daily_max_8hr_avg_ozone_udaq, color='k',marker='.', label = 'Non-Exceedance Day')
        # ax[0].plot(df_ozone.index, daily_max_8hr_avg_ozone_UDAQ, color='r', linewidth = 3, label = 'Exceedance')
        mda8_exceedance_udaq = np.where(daily_max_8hr_avg_ozone_udaq > threshold, daily_max_8hr_avg_ozone_udaq, np.nan)
        ax[0].plot(daily_max_8hr_avg_ozone_udaq.index, mda8_exceedance_udaq, linestyle = '', marker = 'o', color='r', label = 'Exceedance Day')
        print(mda8_exceedance_udaq)

        ax[0].set_ylabel('MDA8 Ozone (ppb)')

        ax[0].hlines(y=70, xmin = daily_max_8hr_avg_ozone_udaq.index[0], xmax = daily_max_8hr_avg_ozone_udaq.index[len(daily_max_8hr_avg_ozone_udaq)-1], color = 'r', linestyle = 'dashed')

        ax[0].set_xticks(daily_max_8hr_avg_ozone_udaq.index)
        ax[0].set_xticklabels(daily_max_8hr_avg_ozone_udaq.index.strftime('%m/%d'))
        ax[0].set_ylim(35,90)

        # midnight_vals_30min = []
        # for midnight_idx_30min in range(16,len(df_adjuststart_30min_ozone.index),48):
        #     midnight_vals_30min.append(df_adjuststart_30min_ozone.index[midnight_idx_30min])
        for day_pos_30min in daily_max_8hr_avg_ozone_udaq.index:
            ax[0].axvline(day_pos_30min, color = 'black', linestyle = 'dotted')


        lines, labels = ax[0].get_legend_handles_labels()
        #lines2, labels2 = ax2.get_legend_handles_labels()
        ax[0].legend(lines, labels, loc=legend_loc)
        ax[0].set_title('UDAQ')


        ax[1].margins(x=0)

        ax[1].plot(daily_max_8hr_avg_ozone_ml.index, daily_max_8hr_avg_ozone_ml, color='k',marker='.', label = 'Non-Exceedance Day')
        mda8_exceedance_ml = np.where(daily_max_8hr_avg_ozone_ml > threshold, daily_max_8hr_avg_ozone_ml, np.nan)
        # ax[1].plot(df_ozone.index, daily_max_8hr_avg_ozone_UDAQ, color='r', linewidth = 3, label = 'Exceedance')
        ax[1].plot(daily_max_8hr_avg_ozone_ml.index, mda8_exceedance_ml, linestyle = '', marker='o', color='r', label = 'Exceedance Day')

        ax[1].set_ylabel('MDA8 Ozone (ppb)')
        ax[1].set_xlabel('Date')

        ax[1].hlines(y=70, xmin = daily_max_8hr_avg_ozone_ml.index[0], xmax = daily_max_8hr_avg_ozone_ml.index[len(daily_max_8hr_avg_ozone_ml)-1], color = 'r', linestyle = 'dashed')

        ax[1].set_xticks(daily_max_8hr_avg_ozone_ml.index)
        ax[1].set_xticklabels(daily_max_8hr_avg_ozone_ml.index.strftime('%m/%d'))
        ax[1].set_ylim(35,90)

        # midnight_vals_30min = []
        # for midnight_idx_30min in range(16,len(df_adjuststart_30min_ozone.index),48):
        #     midnight_vals_30min.append(df_adjuststart_30min_ozone.index[midnight_idx_30min])
        for day_pos_30min in daily_max_8hr_avg_ozone_ml.index:
            ax[1].axvline(day_pos_30min, color = 'black', linestyle = 'dotted')

        lines, labels = ax[1].get_legend_handles_labels()
        #lines2, labels2 = ax2.get_legend_handles_labels()
        ax[1].legend(lines, labels, loc=legend_loc)
        ax[1].set_title('Mobile Lab')
        plt.savefig(compare_udaq_ml_savepath + 'hawthorne_udaq_ml_' + SavePlotSpeciesName + '.png', dpi =300)
        plt.show()
    def mda8_ozone_one_per_day_same_plot(self, legend_loc, SavePlotSpeciesName):
        fig, ax = plt.subplots(figsize = (24,4), tight_layout = True)
        threshold = 70
        ax.margins(x=0)
        plt.plot(daily_max_8hr_avg_ozone_udaq.index, daily_max_8hr_avg_ozone_udaq, color='g',marker='.', label = 'UDAQ Non-Exceedance Day')
        # ax[0].plot(df_ozone.index, daily_max_8hr_avg_ozone_UDAQ, color='r', linewidth = 3, label = 'Exceedance')
        mda8_exceedance_udaq = np.where(daily_max_8hr_avg_ozone_udaq > threshold, daily_max_8hr_avg_ozone_udaq, np.nan)
        mda8_exceedance_ml = np.where(daily_max_8hr_avg_ozone_ml > threshold, daily_max_8hr_avg_ozone_ml, np.nan)
        plt.plot(daily_max_8hr_avg_ozone_udaq.index, mda8_exceedance_udaq, linestyle = '', marker = 'o', color='r')
        plt.plot(daily_max_8hr_avg_ozone_ml.index, daily_max_8hr_avg_ozone_ml, color='m',marker='.', label = 'Mobile Lab Non-Exceedance Day')
        plt.plot(daily_max_8hr_avg_ozone_ml.index, mda8_exceedance_ml, linestyle = '', marker = 'o', color='r', label = 'Exceedance Day')

        print(mda8_exceedance_udaq)

        ax.set_ylabel('MDA8 Ozone (ppb)')

        ax.hlines(y=70, xmin = daily_max_8hr_avg_ozone_udaq.index[0], xmax = daily_max_8hr_avg_ozone_udaq.index[len(daily_max_8hr_avg_ozone_udaq)-1], color = 'r', linestyle = 'dashed')

        ax.set_xticks(daily_max_8hr_avg_ozone_udaq.index)
        ax.set_xticklabels(daily_max_8hr_avg_ozone_udaq.index.strftime('%m/%d'))
        ax.set_ylim(35,90)

        # midnight_vals_30min = []
        # for midnight_idx_30min in range(16,len(df_adjuststart_30min_ozone.index),48):
        #     midnight_vals_30min.append(df_adjuststart_30min_ozone.index[midnight_idx_30min])
        for day_pos_30min in daily_max_8hr_avg_ozone_udaq.index:
            ax.axvline(day_pos_30min, color = 'black', linestyle = 'dotted')


        lines, labels = ax.get_legend_handles_labels()
        #lines2, labels2 = ax2.get_legend_handles_labels()
        ax.legend(lines, labels, loc=legend_loc)
        #ax.set_title('Exceedance Days')


        # ax[1].margins(x=0)

        # ax[1].plot(daily_max_8hr_avg_ozone_ML.index, daily_max_8hr_avg_ozone_ML, color='k',marker='.', label = 'Non-Exceedance Day')
        # mda8_exceedance_ml = np.where(daily_max_8hr_avg_ozone_ML > threshold, daily_max_8hr_avg_ozone_ML, np.nan)
        # # ax[1].plot(df_ozone.index, daily_max_8hr_avg_ozone_UDAQ, color='r', linewidth = 3, label = 'Exceedance')
        # ax[1].plot(daily_max_8hr_avg_ozone_ML.index, mda8_exceedance_ml, linestyle = '', marker='o', color='r', label = 'Exceedance Day')

        # ax[1].set_ylabel('MDA8 Ozone (ppb)')
        # ax[1].set_xlabel('Date')

        # ax[1].hlines(y=70, xmin = daily_max_8hr_avg_ozone_ML.index[0], xmax = daily_max_8hr_avg_ozone_ML.index[len(daily_max_8hr_avg_ozone_ML)-1], color = 'r', linestyle = 'dashed')

        # ax[1].set_xticks(daily_max_8hr_avg_ozone_ML.index)
        # ax[1].set_xticklabels(daily_max_8hr_avg_ozone_ML.index.strftime('%m/%d'))
        # ax[1].set_ylim(35,90)

        # # midnight_vals_30min = []
        # # for midnight_idx_30min in range(16,len(df_adjuststart_30min_ozone.index),48):
        # #     midnight_vals_30min.append(df_adjuststart_30min_ozone.index[midnight_idx_30min])
        # for day_pos_30min in daily_max_8hr_avg_ozone_ML.index:
        #     ax[1].axvline(day_pos_30min, color = 'black', linestyle = 'dotted')

        # lines, labels = ax[1].get_legend_handles_labels()
        # #lines2, labels2 = ax2.get_legend_handles_labels()
        # ax[1].legend(lines, labels, loc=legend_loc)
        # ax[1].set_title('Mobile Lab')
        plt.savefig(compare_udaq_ml_savepath + 'hawthorne_udaq_ml_' + SavePlotSpeciesName + '.png', dpi =300)
        plt.show()
    def mda8_ozone_frequency_histogram_udaq(self, SavePlotSpeciesName):
        idx_max_8hr_rolling_avg = df_mda8_o3_data['8hr_rolling_avg_UDAQ O3'].groupby(df_mda8_o3_data.index.date).idxmax()
        hrs_for_max8hr_rollingavg = []
        for val in idx_max_8hr_rolling_avg:
            hrs_for_max8hr_rollingavg.append(val.hour)

        # Count the occurrences of each hour
        hour_counts_max8hr_rollingavg = pd.Series(hrs_for_max8hr_rollingavg).value_counts().sort_index()
        # print(hour_counts_max8hr_rollingavg)

        all_hours = pd.Index(range(24), name='hour')

        hour_counts_max8hr_rollingavg = hour_counts_max8hr_rollingavg.reindex(all_hours, fill_value=0)

        hour_range = np.arange(0,24,1)
        plt.figure(figsize=(10,6), tight_layout=True)
        plt.hist(hour_counts_max8hr_rollingavg.index, 
                 bins=hour_range, 
                 weights = hour_counts_max8hr_rollingavg.values,
                 edgecolor='black')
        plt.xlabel('Hour')
        plt.xticks(hour_range)
        plt.ylabel('Contribution to MDA8')
        plt.show()

        ranges = idx_max_8hr_rolling_avg.to_frame(name='end')

        hours = pd.Series(8, index=ranges.index)
        hours.iloc[0] = 6  # first row exception
        ranges['start'] = ranges['end'] - pd.to_timedelta(hours, unit='h')

        # reorder columns
        ranges = ranges[['start', 'end']]
        # print(ranges)

        hours = []

        for start, end in zip(ranges['start'], ranges['end']):
            hours += pd.date_range(start=start, end=end, freq='h').hour.tolist()

        # Count the occurrences of each hour
        hour_counts = pd.Series(hours).value_counts().sort_index()
        # print(hour_counts)

        all_hours = pd.Index(range(24), name='hour')

        hour_counts = hour_counts.reindex(all_hours, fill_value=0)

        plt.figure(figsize=(10,6), tight_layout=True)
        plt.hist(
            hour_counts.index, 
            bins=hour_range, 
            weights = hour_counts.values, 
            edgecolor='black')
        plt.xlabel('Hour')
        plt.xticks(hour_range)
        plt.ylabel('Contribution to MDA8')
        #plt.title('Frequency of Hour Counts (Including Zeros)')
        plt.show()
    def mda8_ozone_frequency_histogram_ml(self, SavePlotSpeciesName):
        idx_max_8hr_rolling_avg = df_mda8_o3_data['8hr_rolling_avg_ML O3'].groupby(df_mda8_o3_data.index.date).idxmax()
        hrs_for_max8hr_rollingavg = []
        for val in idx_max_8hr_rolling_avg:
            hrs_for_max8hr_rollingavg.append(val.hour)

        # Count the occurrences of each hour
        hour_counts_max8hr_rollingavg = pd.Series(hrs_for_max8hr_rollingavg).value_counts().sort_index()
        # print(hour_counts_max8hr_rollingavg)

        all_hours = pd.Index(range(24), name='hour')

        hour_counts_max8hr_rollingavg = hour_counts_max8hr_rollingavg.reindex(all_hours, fill_value=0)

        hour_range = np.arange(0,24,1)
        plt.figure(figsize=(10,6), tight_layout=True)
        plt.hist(hour_counts_max8hr_rollingavg.index, 
                 bins=hour_range, 
                 weights = hour_counts_max8hr_rollingavg.values,
                 edgecolor='black')
        plt.xlabel('Hour')
        plt.xticks(hour_range)
        plt.ylabel('Contribution to MDA8')
        plt.show()

        ranges = idx_max_8hr_rolling_avg.to_frame(name='end')

        hours = pd.Series(8, index=ranges.index)
        hours.iloc[0] = 6  # first row exception
        ranges['start'] = ranges['end'] - pd.to_timedelta(hours, unit='h')

        # reorder columns
        ranges = ranges[['start', 'end']]
        # print(ranges)

        hours = []

        for start, end in zip(ranges['start'], ranges['end']):
            hours += pd.date_range(start=start, end=end, freq='h').hour.tolist()

        # Count the occurrences of each hour
        hour_counts = pd.Series(hours).value_counts().sort_index()
        # print(hour_counts)

        all_hours = pd.Index(range(24), name='hour')

        hour_counts = hour_counts.reindex(all_hours, fill_value=0)

        plt.figure(figsize=(10,6), tight_layout=True)
        plt.hist(
            hour_counts.index, 
            bins=hour_range, 
            weights = hour_counts.values, 
            edgecolor='black')
        plt.xlabel('Hour')
        plt.xticks(hour_range)
        plt.ylabel('Contribution to MDA8')
        #plt.title('Frequency of Hour Counts (Including Zeros)')
        plt.show()

class diurnal_comparison:
#    def avg_diurnal_with_std():
    def __init__(self, species_name, var_name_modification=None):
        self.species_name = species_name
        self.var_name_modification = var_name_modification
    def diurnal_mean_median_with_quartiles(self, legend_loc, SavePlotSpeciesName):
        if self.var_name_modification is not None:
            species_name_var = self.var_name_modification
        else:
            species_name_var = self.species_name
            
        if self.species_name == 'Ozone':
            mean_hrly_udaq = mean_hrly_udaq_ozone
            median_hrly_udaq =  median_hrly_udaq_ozone
            udaq_percentile_25 = udaq_percentile_25_ozone
            udaq_percentile_75 = udaq_percentile_75_ozone

            mean_hrly_ml = mean_hrly_ml_ozone
            median_hrly_ml = median_hrly_ml_ozone
            ml_percentile_25 = ml_percentile_25_ozone
            ml_percentile_75 = ml_percentile_75_ozone

        else:
            mean_hrly_udaq=df_udaq_voc_data.groupby('hour')[self.species_name].apply(np.nanmean)
            median_hrly_udaq = df_udaq_voc_data.groupby('hour')[self.species_name].apply(np.nanmedian)
            udaq_percentile_25=df_udaq_voc_data.groupby('hour')[self.species_name].apply(lambda x: np.nanpercentile(x, 25))
            udaq_percentile_75=df_udaq_voc_data.groupby('hour')[self.species_name].apply(lambda x: np.nanpercentile(x, 75))

            mean_hrly_ml=df_ml_voc_data.groupby('hour')[self.species_name].apply(np.nanmean)
            median_hrly_ml = df_ml_voc_data.groupby('hour')[self.species_name].apply(np.nanmedian)
            ml_percentile_25=df_ml_voc_data.groupby('hour')[self.species_name].apply(lambda x: np.nanpercentile(x, 25))
            ml_percentile_75=df_ml_voc_data.groupby('hour')[self.species_name].apply(lambda x: np.nanpercentile(x, 75))

        plt.figure(figsize=(10,6), tight_layout = True)
        udaq_plot, = plt.plot(hour_range, median_hrly_udaq, '-', label='UDAQ Median', color='g')
        ml_plot, = plt.plot(hour_range, median_hrly_ml, '-', label='Mobile Lab Median', color='m')

        quartiles_udaq = plt.fill_between(x=hour_range, y1=udaq_percentile_25, y2=udaq_percentile_75, label= 'UDAQ 25th & 75 Percentile',  color = 'g', alpha = 0.2)
        quartiles_ml = plt.fill_between(x=hour_range, y1=ml_percentile_25, y2=ml_percentile_75, label= 'Mobile Lab 25th & 75 Percentile',  color = 'm', alpha = 0.2)

        mean_dashed_udaq, = plt.plot(hour_range, mean_hrly_udaq, linestyle = 'dashed', color = 'g')
        mean_dashed_ml, = plt.plot(hour_range, mean_hrly_ml, linestyle = 'dashed', color = 'm')
        # Custom legend handle with two lines
        combined_dashed = (mean_dashed_udaq, mean_dashed_ml)
        shade_handle1 = Patch(facecolor='g', alpha=0.2)
        shade_handle2 = Patch(facecolor='m', alpha=0.2)
        combined_shades = (shade_handle1, shade_handle2)

        plt.margins(x=0)
        plt.xlabel('Hour (MDT)')
        plt.ylabel(species_name_var + ' (ppb)')
        plt.xticks(hour_range)
        plt.legend(
            handles = [udaq_plot, ml_plot,
            combined_shades, combined_dashed],
            labels = ['UDAQ Median','ML Median', '25th & 75 Percentile', 'Mean'],
            handler_map = {tuple: HandlerTuple(ndivide=None)},
            loc = legend_loc
        )
        plt.savefig(compare_udaq_ml_savepath + 'hawthorne_udaq_ml_' + 'diurnal_' + SavePlotSpeciesName + '.png', dpi =300)
        plt.show()

class mean_bias_plots:
    def __init__(self, species_name):
        self.species_name = species_name
    def scatter_calculation(self):
        if self.species_name == "Ozone":
            df_species = pd.DataFrame({'obs':df_mda8_o3_data['Mobile Lab O3'], 'model':df_mda8_o3_data['UDAQ O3']})
        else:
            df_species = pd.DataFrame({'obs':df_ml_voc_data[self.species_name], 'model':df_udaq_voc_data[self.species_name]})
        
        print('mean observed value:', df_species['obs'].mean())
        # 1. Compute pointwise error
        df_species['error'] = df_species['model'] - df_species['obs']

        # 2. Overall scatter (standard deviation of error)
        overall_std = df_species['error'].std()  # pandas std ignores NaNs by default
        print("Overall standard deviation of error:", overall_std)

        # 3. Optional: Coefficient of Variation (relative scatter)
        mean_observed = df_species['obs'].mean()
        cv_error = overall_std / mean_observed
        print("Coefficient of variation of error:", cv_error)

        # 4. Hourly scatter (standard deviation per hour of day)
        hourly_std = df_species.groupby(df_species.index.hour)["error"].std()
        print("Hourly standard deviation of error:\n", hourly_std)

    def drift_plot(self):
        if self.species_name == "Ozone":
            df_o3 = pd.DataFrame({'obs O3':df_mda8_o3_data['Mobile Lab O3'], 'model O3':df_mda8_o3_data['UDAQ O3']})
            fig, (ax1, ax2) = plt.subplots(2,1, figsize = (16,8), tight_layout=True)
            df_o3['O3_diff']= df_o3['model O3']-df_o3['obs O3']
            xlim_start_jul = pd.Timestamp('2024-07-16 00:00:00')
            xlim_end_jul = pd.to_datetime('2024-07-31 23:00:00')

            ax1.set_xlim([xlim_start_jul, xlim_end_jul])
            ax1.plot(df_o3.index, df_o3['O3_diff'], color = 'c')

            #Set x ticks
            ax1.xaxis.set_major_locator(mdates.DayLocator())
            ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
            # Minor ticks: every 3 hours
            ax1.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21]))
            # Rotate and format tick labels
            ax1.tick_params(axis='x', which='major')
            ax1.tick_params(axis='x', which='minor', length=3, color='gray')
            ax1.grid(which='both', axis='both')
    
            ax1.set_ylabel('Difference (ppb)')

            ax2.plot(df_o3.index, df_o3['O3_diff'], color = 'c')
            ax2.set_xlim([pd.to_datetime('2024-08-01 00:00:00'), pd.to_datetime('2024-08-18 23:00:00')])
            
            #Set x ticks
            ax2.xaxis.set_major_locator(mdates.DayLocator())
            ax2.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
            # Minor ticks: every 3 hours
            ax2.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21]))
            # Rotate and format tick labels
            ax2.tick_params(axis='x', which='major')
            ax2.tick_params(axis='x', which='minor', length=3, color='gray')
            ax2.grid(which='both', axis='both')
            
            ax2.set_ylabel('Difference (ppb)')

            plt.show()

        else:
            df_voc = pd.DataFrame({'obs VOC':df_ml_voc_data[self.species_name], 'model VOC':df_udaq_voc_data[self.species_name]})
            fig, (ax1, ax2) = plt.subplots(2,1, figsize = (16,8), tight_layout=True)
            df_voc['VOC_diff']= df_voc['model VOC']-df_voc['obs VOC']
            xlim_start_jul = pd.Timestamp('2024-07-15 03:00:00')
            xlim_end_jul = pd.to_datetime('2024-07-31 23:00:00')

            ax1.set_xlim([xlim_start_jul, xlim_end_jul])
            ax1.plot(df_voc.index, df_voc['VOC_diff'], color = 'c')

            #Set x ticks
            ax1.xaxis.set_major_locator(mdates.DayLocator())
            ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
            # Minor ticks: every 3 hours
            ax1.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21]))
            # Rotate and format tick labels
            ax1.tick_params(axis='x', which='major')
            ax1.tick_params(axis='x', which='minor', length=3, color='gray')
            ax1.grid(which='both', axis='both')
            
            ax1.set_ylabel('Difference (ppb)')

            ax2.plot(df_voc.index, df_voc['VOC_diff'], color = 'c')
            ax2.set_xlim([pd.to_datetime('2024-08-01 00:00:00'), pd.to_datetime('2024-08-12 23:00:00')])
            
            #Set x ticks
            ax2.xaxis.set_major_locator(mdates.DayLocator())
            ax2.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
            # Minor ticks: every 3 hours
            ax2.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21]))
            # Rotate and format tick labels
            ax2.tick_params(axis='x', which='major')
            ax2.tick_params(axis='x', which='minor', length=3, color='gray')
            ax2.grid(which='both', axis='both')
            
            ax2.set_ylabel('Difference (ppb)')
            
            plt.show()

    def mean_bias(self, legend_loc, SavePlotSpeciesName, InstrumentUncertaintyUDAQ=None):
        if self.species_name == "Ozone":
            df_o3 = pd.DataFrame({'obs O3':df_mda8_o3_data['Mobile Lab O3'], 'model O3':df_mda8_o3_data['UDAQ O3']})
            df_o3['O3_diff']= df_o3['model O3']-df_o3['obs O3']
            mb_total = df_o3['O3_diff'].dropna().mean()

            df_o3['hour']=df_o3.index.hour 
            hourly_MB = df_o3.groupby(df_o3.index.hour)['O3_diff'].mean()

            # Plot the hourly Mean Bias and the average over time: 
            plt.figure(figsize=(10,6), tight_layout = True)
            plt.plot(hourly_MB.index, hourly_MB, color='g', marker='.', label=f"Hrly MB (Min.={np.min(hourly_MB):.2f} ppb)")
            plt.hlines(y=0, xmin=hourly_MB.index[0], xmax= hourly_MB.index[len(hourly_MB)-1], linestyle='solid', color = 'k')
            plt.plot(hourly_MB.index, np.ones(len(hourly_MB.index))*mb_total, linestyle = 'dashed', color='g',label=f"Avg. MB={mb_total:.2f} ppb")
            plt.fill_between(hourly_MB.index, -1.5, 1.5,  color = 'm', alpha = 0.2, label = 'Instrum. Uncertainty = $\pm$ 1.5 ppb')

            #plt.ylim([-0.14, 0])
            #plt.title('Mean Bias')
            plt.ylabel('Mean Bias (ppb)')
            plt.xlabel('Hour (MDT)')
            plt.xlim([0, 23])
            plt.xticks(hour_range)
            #plt.yticks(np.arange(-7, 3, 1))
            plt.grid()
            plt.legend(loc=legend_loc)
            plt.savefig(compare_udaq_ml_savepath + 'hawthorne_udaq_ml_' + 'mean_bias_' + SavePlotSpeciesName + '.png', dpi =300)
            plt.show()
        else:
            df_voc = pd.DataFrame({'obs VOC':df_ml_voc_data[self.species_name], 'model VOC':df_udaq_voc_data[self.species_name]})
            df_voc['VOC_diff']= df_voc['model VOC']-df_voc['obs VOC']
            mb_total = df_voc['VOC_diff'].dropna().mean()

            df_voc['hour']=df_voc.index.hour 
            hourly_MB = df_voc.groupby(df_voc.index.hour)['VOC_diff'].mean()
            
            # Plot the hourly Mean Bias and the average over time: 
            plt.figure(figsize=(10,6), tight_layout = True)
            plt.plot(hourly_MB.index, hourly_MB, color='g', marker='.', label=f"Hrly MB (Min.={np.min(hourly_MB):.2f} ppb)")
            plt.plot(hourly_MB.index, np.ones(len(hourly_MB.index))*mb_total, linestyle = 'dashed', color='g',label=f"Avg. MB={mb_total:.2f} ppb")

            udaq_instr_uncertainty = InstrumentUncertaintyUDAQ
            if udaq_instr_uncertainty is not None:
            #plot uncertainty in shaded
                plt.fill_between(hourly_MB.index, -1*udaq_instr_uncertainty, udaq_instr_uncertainty,  color = 'm', alpha = 0.2, label = 'Instrum. Uncertainty = $\pm$ '+ str(udaq_instr_uncertainty)+' ppb')
                if any(val > 0 for val in hourly_MB.values) or udaq_instr_uncertainty >= 0: 
                    plt.hlines(y=0, xmin=hourly_MB.index[0], xmax= hourly_MB.index[len(hourly_MB)-1], linestyle='solid', color = 'k')
                else:
                    pass
            else:
                if any(val > 0 for val in hourly_MB.values):
                    plt.hlines(y=0, xmin=hourly_MB.index[0], xmax= hourly_MB.index[len(hourly_MB)-1], linestyle='solid', color = 'k')
                else:
                    pass

            #plt.ylim([-0.14, 0])
            #plt.title('Mean Bias')
            plt.ylabel('Mean Bias (ppb)')
            plt.xlabel('Hour (MDT)')
            plt.xlim([0, 23])
            plt.xticks(hour_range)
            #plt.yticks(np.arange(-7, 3, 1))
            plt.grid()
            plt.legend(loc=legend_loc)
            plt.savefig(compare_udaq_ml_savepath + 'hawthorne_udaq_ml_' + 'mean_bias_' + SavePlotSpeciesName + '.png', dpi =300)
            plt.show()

    def normalized_mean_bias(self, legend_loc, SavePlotSpeciesName, InstrumentUncertaintyUDAQ=None):
        if self.species_name == "Ozone":
            df_o3 = pd.DataFrame({'obs O3':df_mda8_o3_data['Mobile Lab O3'], 'model O3':df_mda8_o3_data['UDAQ O3']})
            mask = df_o3['obs O3'].notna() & df_o3['model O3'].notna() & (df_o3['obs O3'] > 0)
            overall_nmb = (
                (df_o3.loc[mask, 'model O3'] - df_o3.loc[mask, 'obs O3']).sum() /
                df_o3.loc[mask, 'obs O3'].sum()
            )

            # Group by hour of day and compute NMB
            df_valid=df_o3.loc[mask]
            hourly_nmb = df_valid.groupby(df_valid.index.hour).apply(
                lambda x: (x['model O3'] - x['obs O3']).sum() / x['obs O3'].sum()
            )

            plt.figure(figsize=(10,6), tight_layout = True)
            plt.plot(hourly_nmb.index, hourly_nmb*100, color='b', marker='.', label=f"Hrly MNB (Min. ={np.min(hourly_nmb)*100:.2f}%)")
            plt.plot(hourly_nmb.index, np.ones(len(hourly_nmb.index))*overall_nmb*100, linestyle='dashed', color='b',label=f"Avg. MNB ={overall_nmb*100:.2f}%")
        
            plt.fill_between(hourly_nmb.index, -2, 2,  color = 'c', alpha = 0.2, label = 'Instrum. Uncertainty = $\pm$2%')
            plt.hlines(y=0, xmin=hourly_nmb.index[0], xmax= hourly_nmb.index[len(hourly_nmb)-1], linestyle='solid', color = 'k')

            #plt.ylim([-13, 6])
            #plt.yticks(np.arange(-13,6,1))
        else:
            df_voc = pd.DataFrame({'obs VOC':df_ml_voc_data[self.species_name], 'model VOC':df_udaq_voc_data[self.species_name]})
            mask = df_voc['obs VOC'].notna() & df_voc['model VOC'].notna() & (df_voc['obs VOC'] > 0)
            overall_nmb = (
                (df_voc.loc[mask, 'model VOC'] - df_voc.loc[mask, 'obs VOC']).sum() /
                df_voc.loc[mask, 'obs VOC'].sum()
            )

            # Group by hour of day and compute NMB
            df_valid=df_voc.loc[mask]
            hourly_nmb = df_valid.groupby(df_valid.index.hour).apply(
                lambda x: (x['model VOC'] - x['obs VOC']).sum() / x['obs VOC'].sum()
            )

            plt.figure(figsize=(10,6), tight_layout = True)
            plt.plot(hourly_nmb.index, hourly_nmb*100, color='b', marker='.', label=f"Hrly MNB (Min. ={np.min(hourly_nmb)*100:.2f}%)")
            plt.plot(hourly_nmb.index, np.ones(len(hourly_nmb.index))*overall_nmb*100, linestyle='dashed', color='b',label=f"Avg. MNB ={overall_nmb*100:.2f}%")

            udaq_instr_uncertainty = InstrumentUncertaintyUDAQ
            if udaq_instr_uncertainty is not None:
            #plot uncertainty in shaded
                plt.fill_between(hourly_nmb.index, -1*udaq_instr_uncertainty, udaq_instr_uncertainty,  color = 'b', alpha = 0.2, label = 'Instrum. Uncertainty = $\pm$ '+ str(udaq_instr_uncertainty)+'%')
                if any(val > 0 for val in hourly_nmb.values) or udaq_instr_uncertainty >= 0: 
                    plt.hlines(y=0, xmin=hourly_nmb.index[0], xmax= hourly_nmb.index[len(hourly_nmb)-1], linestyle='solid', color = 'k')
                else:
                    pass
            else:
                if any(val > 0 for val in hourly_nmb.values):
                    plt.hlines(y=0, xmin=hourly_nmb.index[0], xmax= hourly_nmb.index[len(hourly_nmb)-1], linestyle='solid', color = 'k')
                else:
                    pass

        plt.ylabel('Normalized Mean Bias (%)')
        plt.xlabel('Hour (MDT)')
        plt.xlim([0, 23])
        plt.xticks(hour_range)
        plt.grid()
        plt.legend(loc=legend_loc)
        plt.savefig(compare_udaq_ml_savepath + 'hawthorne_udaq_ml_' + 'normalized_mean_bias_' + SavePlotSpeciesName + '.png', dpi =300)
        plt.show()

    def mean_normalized_bias(self, legend_loc, SavePlotSpeciesName, InstrumentUncertaintyUDAQ=None):
        if self.species_name == "Ozone":
            df_mnb_species = pd.DataFrame({'obs':df_mda8_o3_data['Mobile Lab O3'], 'model':df_mda8_o3_data['UDAQ O3']})
        else:
            df_mnb_species = pd.DataFrame({'obs':df_ml_voc_data[self.species_name], 'model':df_udaq_voc_data[self.species_name]})

        df_mnb_species['species_diff']= df_mnb_species['model']-df_mnb_species['obs']
        df_mnb_div = df_mnb_species['species_diff']/df_mnb_species['obs']
        mnb_total = df_mnb_div.dropna().mean()
        print(mnb_total)
        df_mnb_div = pd.DataFrame({'mnb_div':df_mnb_div})
        df_mnb_div['hour']=df_mnb_div.index.hour 
        mnb_hrly = df_mnb_div.groupby('hour')['mnb_div'].mean()
        print(mnb_hrly)

        plt.figure(figsize=(10,6), tight_layout = True)
        plt.plot(mnb_hrly.index, mnb_hrly*100, color='b', marker='.', label=f"Hrly MNB (Min. ={np.min(mnb_hrly)*100:.2f}%)")
        plt.plot(mnb_hrly.index, np.ones(len(mnb_hrly.index))*mnb_total*100, linestyle='dashed', color='b',label=f"Avg. MNB ={mnb_total*100:.2f}%")
    
        if self.species_name == "Ozone":
            plt.fill_between(mnb_hrly.index, -2, 2,  color = 'c', alpha = 0.2, label = 'Instrum. Uncertainty = $\pm$2%')
        else:
            pass
        
        udaq_instr_uncertainty = InstrumentUncertaintyUDAQ
        if udaq_instr_uncertainty is not None:
        #plot uncertainty in shaded
            plt.fill_between(mnb_hrly.index, -1*udaq_instr_uncertainty, udaq_instr_uncertainty,  color = 'b', alpha = 0.2, label = 'Instrum. Uncertainty = $\pm$ '+ str(udaq_instr_uncertainty)+'%')
            if any(val > 0 for val in mnb_hrly.values) or udaq_instr_uncertainty >= 0: 
                plt.hlines(y=0, xmin=mnb_hrly.index[0], xmax= mnb_hrly.index[len(mnb_hrly)-1], linestyle='solid', color = 'k')
            else:
                pass
        else:
            if any(val > 0 for val in mnb_hrly.values):
                plt.hlines(y=0, xmin=mnb_hrly.index[0], xmax= mnb_hrly.index[len(mnb_hrly)-1], linestyle='solid', color = 'k')
            else:
                pass
        
        plt.ylabel('Normalized Mean Bias (%)')
        plt.xlabel('Hour (MDT)')
        plt.xlim([0, 23])
        plt.xticks(hour_range)
        plt.grid()
        plt.legend(loc=legend_loc)
        plt.savefig(compare_udaq_ml_savepath + 'hawthorne_udaq_ml_' + 'mean_normalized_bias_' + SavePlotSpeciesName + '.png', dpi =300)
        plt.show()

        