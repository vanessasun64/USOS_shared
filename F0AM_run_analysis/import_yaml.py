import yaml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator
#import cmasher as cmr
import matplotlib.colors as mcolors
import matplotlib.dates as mdates

import IPython
notebook_name = "/".join(
        IPython.extract_module_locals()[1]["__vsc_ipynb_file__"].split("/")[-5:]
    )

#region: Plot formatting
mpl.rcParams['xtick.labelsize'] = 15
mpl.rcParams['ytick.labelsize'] = 15
mpl.rcParams['legend.fontsize'] = 16
mpl.rcParams['axes.labelsize'] = 18
mpl.rcParams['axes.titlesize'] = 28
mpl.rcParams['axes.xmargin'] = 0

# Set the font family to 'serif'
mpl.rcParams['font.family'] = 'serif'
# Specify preferred serif font (Computer Modern Roman is 'cmr10')
mpl.rcParams['font.serif'] = 'Lato' 
# Optionally, configure mathtext to use Computer Modern fonts as well
mpl.rcParams['mathtext.fontset'] = 'cm'
# Ensure minus signs are rendered correctly with CM fonts
mpl.rcParams['axes.unicode_minus'] = False

#Colorblind friendly colors in RGB:
forest_green_color = np.array([51,117,56])
bright_green_color = np.array([51, 117, 56])
indigo_color = np.array([46, 37, 133])
bright_blue_color = np.array([0, 114, 178])
light_blue_color = np.array([86, 180, 233])
pale_blue_color = np.array([148, 203, 236])
yellow_orange_color = np.array([230, 159, 0])
sunny_yellow_color = np.array([230, 159, 20])
lemon_yellow_color = np.array([240, 228, 66])
bright_orange_color = np.array([93, 168, 153])
desaturated_pink_color = np.array([194, 106, 119])
saturated_orchid_color = np.array([194, 106, 119])
purple_orchid_color = np.array([159, 74, 150])


def rgb_range_correction(color_names):
    mpl_formatted_colors = []
    for c in range(0,len(color_names)):
        rgb_division = color_names[c]/255
        mpl_formatted_colors.append(rgb_division)
    return mpl_formatted_colors

four_colorset1 = rgb_range_correction([forest_green_color, bright_blue_color, sunny_yellow_color, desaturated_pink_color])

#endregion

#read YAML files
dirpath = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/'

run_month = '8'
run_number = '48'
mech_name = ['CB6r5h', 'CRACMM', 'GEOS-Chem', 'MCM']

df_CB6r5h_species_conc = pd.DataFrame()
df_CRACMM_species_conc = pd.DataFrame()
df_GEOSChem_species_conc = pd.DataFrame()
df_MCM_species_conc = pd.DataFrame()

df_CB6r5h_analysis = pd.DataFrame()
df_CRACMM_analysis = pd.DataFrame()
df_GEOSChem_analysis = pd.DataFrame()
df_MCM_analysis = pd.DataFrame()

for mech in mech_name:
    for run_day in range(4,8):
        run_date = run_month + '_' + str(run_day) + '_' + '2024'
        yaml_usos_dirpath = dirpath + 'F0AM-4.4.2/Runs/USOS_' + run_date
        yaml_usos_full_filepath_species_conc = yaml_usos_dirpath + '/' + 'Run' + run_number + '/' + mech +'/USOS_' + run_date + '_' + mech + '_' + 'Run' + run_number + '_species_conc_output' + '.yaml'
        yaml_usos_full_filepath_analysis = yaml_usos_dirpath + '/' + 'Run' + run_number + '/' + mech +'/USOS_' + run_date + '_' + mech + '_' + 'Run' + run_number + '_analysis_output' + '.yaml'
        with open(yaml_usos_full_filepath_species_conc, 'r') as f:
            yaml_usos_species_conc = yaml.full_load(f)
        with open(yaml_usos_full_filepath_analysis, 'r') as f:
            yaml_usos_analysis = yaml.full_load(f)

        processed_data_species_conc = {k: [item[0] if isinstance(item, list) else item for item in v] for k, v in yaml_usos_species_conc.items()}
        df_yaml_species_conc = pd.DataFrame(processed_data_species_conc)

        
        if mech == 'CB6r5h':
            df_CB6r5h_species_conc = pd.concat([df_CB6r5h_species_conc, df_yaml_species_conc], ignore_index=True)
        #     print(df_CB6r5h_analysis)
        #     print(df_yaml_analysis)
        #     df_CB6r5h_analysis = pd.concat([df_CB6r5h_analysis, df_yaml_analysis], ignore_index=True)
        elif mech == 'CRACMM':
            df_CRACMM_species_conc = pd.concat([df_CRACMM_species_conc, df_yaml_species_conc], ignore_index=True)
        #     print(df_CRACMM_analysis)
        #     #print(df_yaml_analysis)
        #     #df_CRACMM_analysis = pd.concat([df_CRACMM_analysis, df_yaml_analysis], ignore_index=True)

        elif mech == 'GEOS-Chem':
            df_GEOSChem_species_conc = pd.concat([df_GEOSChem_species_conc, df_yaml_species_conc], ignore_index=True)
        #     print(df_GEOSChem_analysis)
        #     #df_GEOSChem_analysis = pd.concat([df_GEOSChem_analysis, df_yaml_analysis], ignore_index=True)
        elif mech == 'MCM':
            df_MCM_species_conc = pd.concat([df_MCM_species_conc, df_yaml_species_conc], ignore_index=True)
        #     print(df_MCM_analysis)
        #     #df_MCM_analysis = pd.concat([df_MCM_analysis, df_yaml_analysis], ignore_index=True)
        else:
            print('Mechanism name invalid')

new_start_time = pd.Timestamp('2024-08-04 00:00:00')
new_end_time = pd.Timestamp('2024-08-07 23:45:00')
#Create a new datetime index from new_start to the end of existing index with same frequency
new_index = pd.date_range(start=new_start_time, end=new_end_time, freq='15min')

df_CB6r5h_hrly_species_conc = df_CB6r5h_species_conc.copy()
df_CRACMM_hrly_species_conc = df_CRACMM_species_conc.copy()
df_GEOSChem_hrly_species_conc = df_GEOSChem_species_conc.copy()
df_MCM_hrly_species_conc = df_MCM_species_conc.copy()

df_CB6r5h_hrly_species_conc.set_index(new_index, inplace=True)
df_CRACMM_hrly_species_conc.set_index(new_index, inplace=True)
df_GEOSChem_hrly_species_conc.set_index(new_index, inplace=True)
df_MCM_hrly_species_conc.set_index(new_index, inplace=True)

hour_range = np.arange(0,24,1)

def plot_ozone_conc_time_series():
    fig, ax = plt.subplots(figsize = (10,7), layout = 'tight')
    plt.rcParams['figure.figsize'] = [1.618 * i for i in plt.rcParams['figure.figsize']]

    plt.plot(new_index, df_CB6r5h_species_conc['O3Initconc'], label = 'Obs.', linestyle = 'solid', linewidth = 2, color = 'black')
    plt.plot(new_index, df_CB6r5h_species_conc['O3conc'], label = 'CB6r5h', linestyle = 'dashed', linewidth = 2, color = four_colorset1[0].tolist())
    plt.plot(new_index,df_CRACMM_species_conc['O3conc'], label = 'CRACMM', linestyle = 'dashed', linewidth = 2, color = four_colorset1[1].tolist())
    plt.plot(new_index,df_GEOSChem_species_conc['O3conc'], label = 'GEOS-Chem', linestyle = 'dashed', linewidth = 2, color = four_colorset1[2].tolist())
    plt.plot(new_index,df_MCM_species_conc['O3conc'], label = 'MCM', linestyle = 'dashed', linewidth = 2, color = four_colorset1[3].tolist())

    ax.xaxis.set_major_locator(mdates.DayLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
    # Minor ticks: every 3 hours
    ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21]))
    # Rotate and format tick labels
    ax.tick_params(axis='x', which='major')
    ax.tick_params(axis='x', which='minor', length=3, color='gray')

    plt.xlabel('Hour (Mountain Time)')
    # plt.xticks(np.arange(0,24))
    plt.margins(x=0)

    plt.ylabel('Ozone Concentration (ppb)')
    plt.legend()
    #plt.title('Ozone Between Chemical Mechanisms \nfor July 19, 2024\n')
    plt.savefig(dirpath + '/F0AM_run_analysis/plots/ozone_conc_mechanism_comparison_timeseries_08042024_08072024.png', dpi =150)
    plt.show()

    fig, ax = plt.subplots(figsize = (10,7), layout = 'tight')
    plt.rcParams['figure.figsize'] = [1.618 * i for i in plt.rcParams['figure.figsize']]

    plt.plot(new_index, df_CB6r5h_species_conc['O3Initconc']+df_CB6r5h_species_conc['NO2Initconc'], label = 'Obs.', linestyle = 'solid', linewidth = 2, color = 'black')
    plt.plot(new_index, df_CB6r5h_species_conc['O3conc']+df_CB6r5h_species_conc['NO2conc'], label = 'CB6r5h', linestyle = 'dashed', linewidth = 2, color = four_colorset1[0].tolist())
    plt.plot(new_index, df_CRACMM_species_conc['O3conc']+df_CRACMM_species_conc['NO2conc'], label = 'CRACMM', linestyle = 'dashed', linewidth = 2, color = four_colorset1[1].tolist())
    plt.plot(new_index, df_GEOSChem_species_conc['O3conc']+df_GEOSChem_species_conc['NO2conc'], label = 'GEOS-Chem', linestyle = 'dashed', linewidth = 2, color = four_colorset1[2].tolist())
    plt.plot(new_index, df_MCM_species_conc['O3conc']+df_MCM_species_conc['NO2conc'], label = 'MCM', linestyle = 'dashed', linewidth = 2, color = four_colorset1[3].tolist())

    ax.xaxis.set_major_locator(mdates.DayLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
    # Minor ticks: every 3 hours
    ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21]))
    # Rotate and format tick labels
    ax.tick_params(axis='x', which='major')
    ax.tick_params(axis='x', which='minor', length=3, color='gray')

    plt.xlabel('Hour (Mountain Time)')
    # plt.xticks(np.arange(0,24))
    plt.margins(x=0)

    plt.ylabel('O$_x$ Concentration (ppb)')
    plt.legend()
    #plt.title('Ozone Between Chemical Mechanisms \nfor July 19, 2024\n')
    plt.savefig(dirpath + '/F0AM_run_analysis/plots/ox_conc_mechanism_comparison_timeseries_08042024_08072024.png', dpi =150)
    plt.show()

def plot_ozone_hourly_mean():
    #Put data into hour-of-day bin
    df_CB6r5h_hrly_species_conc['hour'] = df_CB6r5h_hrly_species_conc.index.hour
    df_CRACMM_hrly_species_conc['hour'] = df_CRACMM_hrly_species_conc.index.hour
    df_GEOSChem_hrly_species_conc['hour'] = df_GEOSChem_hrly_species_conc.index.hour
    df_MCM_hrly_species_conc['hour'] = df_MCM_hrly_species_conc.index.hour

    #Group by hour and take mean per species
    df_CB6r5h_hourly_means_species_conc = df_CB6r5h_hrly_species_conc.groupby('hour').mean()
    df_CRACMM_hourly_means_species_conc = df_CRACMM_hrly_species_conc.groupby('hour').mean()
    df_GEOSChem_hourly_means_species_conc = df_GEOSChem_hrly_species_conc.groupby('hour').mean()
    df_MCM_hourly_means_species_conc = df_MCM_hrly_species_conc.groupby('hour').mean()

    fig, ax = plt.subplots(figsize = (10,7), layout = 'tight')

    plt.plot(hour_range, df_CB6r5h_hourly_means_species_conc['O3Initconc'], label = 'Obs.', linestyle = 'solid', linewidth = 2, color = 'black')
    plt.plot(hour_range, df_CB6r5h_hourly_means_species_conc['O3conc'], label = 'CB6r5h', linestyle = 'dashed', linewidth = 2, color = four_colorset1[0].tolist())
    plt.plot(hour_range, df_CRACMM_hourly_means_species_conc['O3conc'], label = 'CRACMM', linestyle = 'dashed', linewidth = 2, color = four_colorset1[1].tolist())
    plt.plot(hour_range, df_GEOSChem_hourly_means_species_conc['O3conc'], label = 'GEOS-Chem', linestyle = 'dashed', linewidth = 2, color = four_colorset1[2].tolist())
    plt.plot(hour_range, df_MCM_hourly_means_species_conc['O3conc'], label = 'MCM', linestyle = 'dashed', linewidth = 2, color = four_colorset1[3].tolist())
    
    #ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21]))
    ax.xaxis.set_minor_locator(MultipleLocator(3))
    ax.xaxis.set_major_locator(MultipleLocator(6))
    plt.xlabel('Hour (Mountain Time)')
    plt.xticks(np.arange(0,24,6))
    ax.grid(which='both', linestyle='--', linewidth=0.5)
    plt.margins(x=0)

    plt.ylabel('Avg Ozone Concentration (ppb)')
    plt.legend()
    #plt.title('Ozone Between Chemical Mechanisms \nfor July 19, 2024\n')
    plt.savefig(dirpath + '/F0AM_run_analysis/plots/ozone_conc_mechanism_comparison_08042024_08072024.png', dpi =150)
    plt.show()

def plot_ozone_hourly_median():

    #Put data into hour-of-day bin
    df_CB6r5h_hrly_species_conc['hour'] = df_CB6r5h_hrly_species_conc.index.hour
    df_CRACMM_hrly_species_conc['hour'] = df_CRACMM_hrly_species_conc.index.hour
    df_GEOSChem_hrly_species_conc['hour'] = df_GEOSChem_hrly_species_conc.index.hour
    df_MCM_hrly_species_conc['hour'] = df_MCM_hrly_species_conc.index.hour

    #Group by hour and take mean per species
    df_CB6r5h_hourly_means_species_conc = df_CB6r5h_hrly_species_conc.groupby('hour').median()
    df_CRACMM_hourly_means_species_conc = df_CRACMM_hrly_species_conc.groupby('hour').median()
    df_GEOSChem_hourly_means_species_conc = df_GEOSChem_hrly_species_conc.groupby('hour').median()
    df_MCM_hourly_means_species_conc = df_MCM_hrly_species_conc.groupby('hour').median()

    hour_range = np.arange(0,24,1)

    fig, ax = plt.subplots(figsize = (10,7), layout = 'tight')

    plt.plot(hour_range, df_CB6r5h_hourly_means_species_conc['O3Initconc'], label = 'Obs.', linestyle = 'solid', linewidth = 2, color = 'black')
    plt.plot(hour_range, df_CB6r5h_hourly_means_species_conc['O3conc'], label = 'CB6r5h', linestyle = 'dashed', linewidth = 2, color = four_colorset1[0].tolist())
    plt.plot(hour_range,df_CRACMM_hourly_means_species_conc['O3conc'], label = 'CRACMM', linestyle = 'dashed', linewidth = 2, color = four_colorset1[1].tolist())
    plt.plot(hour_range,df_GEOSChem_hourly_means_species_conc['O3conc'], label = 'GEOS-Chem', linestyle = 'dashed', linewidth = 2, color = four_colorset1[2].tolist())
    plt.plot(hour_range,df_MCM_hourly_means_species_conc['O3conc'], label = 'MCM', linestyle = 'dashed', linewidth = 2, color = four_colorset1[3].tolist())
    
    #ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21]))
    ax.xaxis.set_minor_locator(MultipleLocator(3))
    ax.xaxis.set_major_locator(MultipleLocator(6))
    plt.xlabel('Hour (Mountain Time)')
    plt.xticks(np.arange(0,24,6))
    ax.grid(which='both', linestyle='--', linewidth=0.5)
    plt.margins(x=0)

    plt.ylabel('Median Ozone Concentration (ppb)')
    plt.legend()
    #plt.title('Ozone Between Chemical Mechanisms \nfor July 19, 2024\n')
    plt.savefig(dirpath + '/F0AM_run_analysis/plots/ozone_conc_mechanism_comparison.png', dpi =150)
    plt.show()

def mean_bias():
    df_o3 = pd.DataFrame({'obs O3':df_CB6r5h_hrly_species_conc['O3Initconc'], 'CB6r5h O3':df_CB6r5h_hrly_species_conc['O3conc'], 'CRACMM O3':df_CRACMM_hrly_species_conc['O3conc'], 'GEOS-Chem O3':df_GEOSChem_hrly_species_conc['O3conc'], 'MCM O3':df_MCM_hrly_species_conc['O3conc']})
    mech_varlist = ['CB6r5h O3', 'CRACMM O3', 'GEOS-Chem O3', 'MCM O3']

    plt.figure(figsize=(10,6), tight_layout = True)
    for mech_var in mech_varlist:
        df_o3['O3_diff']= df_o3[mech_var]-df_o3['obs O3']
        mb_total = df_o3['O3_diff'].dropna().mean()
        print(mb_total)

        df_o3['hour']=df_o3.index.hour
        hourly_MB = df_o3.groupby(df_o3.index.hour)['O3_diff'].mean()

        if mech_var == 'CB6r5h O3':
            color_val = four_colorset1[0].tolist()
        elif mech_var == 'CRACMM O3':
            color_val = four_colorset1[1].tolist()
        elif mech_var == 'GEOS-Chem O3':
            color_val = four_colorset1[2].tolist()
        elif mech_var == 'MCM O3':
            color_val = four_colorset1[3].tolist()
        else:
            print('Mechanism name invalid')

        print(color_val)
        # Plot the hourly Mean Bias and the average over time: 
        plt.plot(hourly_MB.index, hourly_MB, color=color_val, marker='.', label=f"{mech_var} Hrly MB (Min.={np.min(hourly_MB):.2f} ppb)")
        plt.hlines(y=0, xmin=hourly_MB.index[0], xmax= hourly_MB.index[len(hourly_MB)-1], linestyle='solid', color = 'k')
        # plt.plot(hourly_MB.index, np.ones(len(hourly_MB.index))*mb_total, linestyle = 'dashed', color=color_val,label=f"Avg. MB={mb_total:.2f} ppb")
        #plt.fill_between(hourly_MB.index, -1.5, 1.5,  color = 'm', alpha = 0.2, label = 'Instrum. Uncertainty = $\pm$ 1.5 ppb')

        #plt.ylim([-0.14, 0])
        #plt.title('Mean Bias')
        plt.ylabel('Mean Bias (ppb)')
        plt.xlabel('Hour (MDT)')
        plt.xlim([0, 23])
        plt.xticks(hour_range)
        #plt.yticks(np.arange(-7, 3, 1))
        plt.grid()
        #plt.legend(loc=legend_loc)
        #plt.savefig(compare_udaq_ml_savepath + 'hawthorne_udaq_ml_' + 'mean_bias_' + SavePlotSpeciesName + '.png', dpi =300)
    plt.show()

def plot_OPE_total():
    print(len(df_CB6r5h_analysis['ope_total']))
    print(len(df_CRACMM_analysis['ope_total']))
    print(len(df_GEOSChem_analysis['ope_total']))
    print(len(df_MCM_analysis['ope_total']))
    
    # fig, ax = plt.subplots(figsize = (12,10), layout = 'tight')
    # plt.plot(new_index, df_CB6r5h_analysis['ope_total'], label = 'CB6r5h', linewidth = 2, color = four_colorset1[0].tolist())
    # plt.plot(new_index, df_CRACMM_analysis['ope_total'], label = 'CRACMM', linewidth = 2, color = four_colorset1[1].tolist())
    # plt.plot(new_index, df_GEOSChem_analysis['ope_total'], label = 'GEOS-Chem', linewidth = 2, color = four_colorset1[2].tolist())
    # plt.plot(new_index, df_MCM_analysis['ope_total'], label = 'MCM',linewidth = 2, color = four_colorset1[3].tolist())
    # plt.legend()
    # plt.margins(x=0,y=0)
    # plt.xlabel('Hour')
    # plt.ylabel('OPE')
    # plt.ylim([0,2])
    # plt.savefig(dirpath + '/F0AM_run_analysis/plots/ope_comparison_timeseries_08042024_08072024.png', dpi =150)
    # plt.show()

#CALL FUNCTIONS
# plot_ozone_conc_time_series()
# plot_ozone_hourly_mean()
# # plot_ozone_hourly_median()
# mean_bias()
plot_OPE_total()