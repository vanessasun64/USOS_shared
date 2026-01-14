import numpy as np 
import os 
import xarray as xr
import pandas as pd

from scipy.io import savemat
from collections import OrderedDict
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
import matplotlib as mpl

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

#region: filepaths
#Filepaths for loading
dirpath = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/'
merged_data_dir = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/CampaignData_and_Merges/R0/CSL_MobileLab_Parked/merged/'
merged_data_dir_15min = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/CampaignData_and_Merges/R0/CSL_MobileLab_Parked/merged/rev_15min/'
#endregion

#Read in data
#ml_voc_data = xr.open_dataset(merged_data_dir_15min +'all_CSL_MobileLab_Parked_rev15min_iWASupdated_formaldehydeupdated.nc')
ml_voc_data = xr.open_dataset(merged_data_dir_15min +'all_CSL_MobileLab_Parked_rev15min_iWASupdated_formaldehydeupdated_raw_noadjustment.nc')
df_ml_data = ml_voc_data.to_dataframe()
df_ml_data.reset_index(inplace=True)
df_ml_data.set_index('time_local', inplace=True)


blh_filepath = dirpath+'Boundary_layer_height/BLH_Nell.csv'
df_blh = pd.read_csv(blh_filepath)
#Convert from Igor Pro Time to local time
df_blh['time_local'] = (pd.to_datetime("1904-01-01") + pd.to_timedelta(df_blh['Time_start_15min_local'], unit="s"))
df_blh = df_blh.rename(columns={'UDAQ_mixing_height_m_15min_avg':'BLH_Nell_m', 'UDAQ_mixing_height_m_15min_avg_smoothed':'BLH_Nell_m_smoothed'})
df_blh = df_blh.set_index(['time_local'])
df_blh = df_blh.drop(columns=['Time_start_15min_local'])
df_ml_data = df_ml_data.join(df_blh, how="inner")   # only matching index values
print(df_ml_data.columns)



# Define which variables we need to make sure don't have Nans/ negs since we'll be using then as constraints in F0AM: 
# True means we can't have holes in them!
need2fill= {
    #Nell constrained
    'O3_ppbv':True,
    'CO_Piccaro':True,
    'CH4_Piccaro':True,
    'NO_LIF':True,
    'NO2_LIF':True,
    'HONO_CIMS':True,
    'Alpha_Pinene_WAS':True,
    'Benzene_PTR':True,
    'Toluene_PTR':True,
    'Styrene_PTR':True,
    'Methanol_PTR':True,
    'Acetaldehyde_PTR':True,
    'Ethanol_PTR':True,
    'HCOOH_CIMS':True,
    'Acrolein_PTR':True,
    'Monoterpenes_PTR':True,
    'Br2_CIMS':True,
    'Cl2_CIMS':True,
    'Isoprene_WAS':True,
    'nButane_WAS':True,
    'iButane_WAS':True,
    'nPentane_WAS':True,
    'iPentane_WAS':True,
    'nHexane_WAS':True,
    'nOctane_WAS':True,
    'nNonane_WAS':True,
    'nDecane_WAS':True,
    'Ethane_WAS':True,
    'Ethene_WAS':True,
    'Propane_WAS':True,
    'Propene_WAS':True,
    't2Butene_WAS':True,
    'c2Butene_WAS':True,
    'x1_Pentene_WAS':True,
    't2Pentene_WAS':True,
    'c2Pentene_WAS':True,
    'x1_Butene_WAS':True,
    'x2_MethylPentane_WAS':True,
    'Ethyl_WAS':True,
    'm_p_Xylene_WAS':True,
    'o_Xylene_WAS':True,
    'C2Cl4_WAS':True,
    'Furan_WAS':True,
    'Benzaldehyde_PTR':True,
    'MVK_MACR_PTR':True,
    'iPropylONO2_WAS':True,
    'ISOPN_CIMS':True,
    'HCHO_CRDS':True,
    'BLH_Nell_m_smoothed':True,
    'BLH_Nell_m':True,
    ###### End of Nell Matching
    ###### Need to be true for F0AM to run:
    'Temp_K':True,
    'Pressure_mb':True,
    'Altitude_m':True,
    'WindSpd_ms':True,
    'jNO2_meas':True,
    'RH_percent':True,
    'H2O_Piccaro':True,
    'jNO2':True,
    'jBrCl':True,
    'jBr2':True,
    'jCCl4':True,
    'jCH2Oa':True,
    'jCH2Ob':True,
    'jClNO2':True,
    'jClOa':True,
    'jClOb':True,
    'jCl2':True,
    'jHNO2':True,
    'jHNO3':True,
    'jI2':True,
    'jNO3a':True,
    'jNO3b':True,
    'jN2O5':True,
    'jO3':True,
    'Surface_area_conc_POPS':True,
    'Lat':True,
    'Lon':True,
    #####
    'BrO_CIMS':False,
    'ClNO2_CIMS':False,
    'BrCl_CIMS':False,
    'NCl3_CIMS':False,
    'File_Index':False,
    'NOy_LIF':False,
    'ACCURACY_1sigma_NO':False,
    'ACCURACY_1sigma_NO2':False,
    'ACCURACY_1sigma_NOy':False,
    'N2O5_CIMS':False,
    'C4H7NO5_CIMS':False,
    'C5H10O3_CIMS':False,
    'C5H9NO5_CIMS':False,
    'C10H17NO4_CIMS':False,
    'HNO3_CIMS':False,
    'Acetonitrile_PTR':False,
    'Methanethiol_PTR':False,
    'Acetone_Propanal_PTR':False,
    'DMS_PTR':False,
    'C8Aromatics_PTR':False,
    'C9Aromatics_PTR':False,
    'Isoprene_PTR':False,
    'Naphthalene_PTR':False,
    'Octanal_PTR':False,
    'Nonanal_PTR':False,
    'C7H4ClF3_PTR':False,
    'D5_siloxane_PTR':False,
    'PAN_CIMS':False,
    'APAN_CIMS':False,
    'PPN_CIMS':False,
    'Course_deg':False,
    'GndSpd_ms':False,
    'Heading_deg':False,
    'WindDir_deg':False,
    'CO2_Piccaro':False,
    'H2O_CRDS':False,
    'CH4_CRDS':False,
    'Time_Start_WAS':False,
    'Time_Stop_WAS':False,
    'Time_Mid_WAS':False,
    'Acetone_WAS':False,
    'Acrolein_WAS':False,
    'Benzene_WAS':False,
    'C2HCl3_WAS':False,
    'CCl4_WAS':False,
    'CF2Cl2_WAS':False,
    'CFCl3_WAS':False,
    'CH2Cl2_WAS':False,
    'CH3Br_WAS':False,
    'CycloPentane_WAS':False,
    'Ethyne_WAS':False,
    'Limonene_WAS':False,
    'MACR_WAS':False,
    'MethylCycloHexane_WAS':False,
    'MethylCycloPentane_WAS':False,
    'Toluene_WAS':False,
    'Beta_Pinene_WAS':False,
    'iPropylBenzene_WAS':False,
    'nHeptane_WAS':False,
    'nPropylBenzene_WAS':False,
    'nPropylONO2_WAS':False,
    'x123_TriMethylBenzene_WAS':False,
    'x124_TriMethylBenzene_WAS':False,
    'x135_TriMethylBenzene_WAS':False,
    'x13_Butadiene_WAS':False,
    'x224_TriMethylPentane_WAS':False,
    'x22_DiMethylButane_WAS':False,
    'x2_EthylToluene_WAS':False,
    'x3_MethylPentane_WAS':False,
    'x3_x4_EthylToluene_WAS':False,
    'AOD':False}
    # 'Time_Start_POPS':False,
    # 'Time_Stop_POPS':False,
    # 'Press_mb_POPS':False,
    # 'Temp_C_POPS':False,
    # 'Aerosol_conc_POPS':False,
    # 'Surface_area_conc_POPS':False,
    # 'Volume_density_POPS':False,
    # 'Effective_radius_POPS':False,
    # 'Extinction_POPS':False,
    # 'Particle_conc_bin01_POPS':False,
    # 'Particle_conc_bin02_POPS':False,
    # 'Particle_conc_bin03_POPS':False,
    # 'Particle_conc_bin04_POPS':False,
    # 'Particle_conc_bin05_POPS':False,
    # 'Particle_conc_bin06_POPS':False,
    # 'Particle_conc_bin07_POPS':False,
    # 'Particle_conc_bin08_POPS':False,
    # 'Particle_conc_bin09_POPS':False,
    # 'Particle_conc_bin10_POPS':False,
    # 'Particle_conc_bin11_POPS':False,
    # 'Particle_conc_bin12_POPS':False,
    # 'Particle_conc_bin13_POPS':False,
    # 'Particle_conc_bin14_POPS':False,
    # 'Particle_conc_bin15_POPS':False,
    # 'Particle_conc_bin16_POPS':False,
    # 'Particle_conc_bin17_POPS':False,
    # 'Particle_conc_bin18_POPS':False,
    # 'Particle_conc_bin19_POPS':False,
    # 'Particle_conc_bin20_POPS':False,
    # 'Particle_conc_bin21_POPS':False,
    # 'Particle_conc_bin22_POPS':False,
    # 'Particle_conc_bin23_POPS':False,
    # 'Particle_conc_bin24_POPS':False,
    # 'Particle_conc_bin25_POPS':False,
    # 'Particle_conc_bin26_POPS':False,
    # 'Particle_conc_bin27_POPS':False,
    # 'Particle_conc_bin28_POPS':False,
    # 'Particle_conc_bin29_POPS':False,
    # 'Particle_conc_bin30_POPS':False,
    # 'Particle_conc_bin31_POPS':False,
    # 'Particle_conc_bin32_POPS':False,
    # 'Particle_conc_bin33_POPS':False,
    # 'Particle_conc_bin34_POPS':False,
    # 'Particle_conc_bin35_POPS':False,
    # 'Particle_conc_bin36_POPS':False}

# Get names of vars we need to fill nans in: 
vars2fill=[key for key,value in need2fill.items() if value ==True]


def plot_interps():
    """
    This function shows the interpolation plots for our 15 minute revised merges. Currently, any time there is data missing from the USOS Campaign, 
    there is a NaN as the value. These plots will show what it looks like to fill a NaN with an interpolated value.
    """
    df_interp=df_ml_data.copy()

    #Nell turned all benzaldehyde, styrene, HONO negative values into 0
    neg_species = ['Benzaldehyde_PTR','Styrene_PTR', 'HONO_CIMS']

    for i,col in enumerate(vars2fill):
        if col in neg_species:
            df_interp[col] = df_interp[col].mask(df_interp[col] < 0, 0)
        else:
            # Set any negative values to NaN so we can interp them... 
            df_interp[col] = df_interp[col].mask(df_interp[col] < 0, np.nan)
                                
        # Calc number of points that are negative or Nans: 
        n_baddies= len([item for item in df_ml_data[col] if item <0 or np.isnan(item)]) 
        
        if n_baddies > 0: 
            #apply the linear interpolation
            df_interp[col] = df_interp[col].interpolate(method='linear')
            
            #Plot it so we can take a look at it... 
            fig, (ax1, ax2, ax3, ax4) = plt.subplots(4,1, figsize = (30,20), tight_layout=True)

            #Each subplot shows approximately 1 week
            xlim_start_w1 = pd.to_datetime('2024-07-15 00:00:00')
            xlim_end_w1 = pd.to_datetime('2024-07-23 23:00:00')
            xlim_start_w2 = pd.to_datetime('2024-07-24 00:00:00')
            xlim_end_w2 = pd.to_datetime('2024-07-31 23:00:00')
            xlim_start_w3 = pd.to_datetime('2024-08-01 00:00:00')
            xlim_end_w3 = pd.to_datetime('2024-08-08 23:00:00')
            xlim_start_w4 = pd.to_datetime('2024-08-09 00:00:00')
            xlim_end_w4 = pd.to_datetime('2024-08-18 23:00:00')

            ax1.plot(df_ml_data.index, df_ml_data[col], color='k', marker='o',label=f'Original (baddies={n_baddies})')
            ax1.plot(df_interp.index, df_interp[col], color='r', marker='x', label='Interpolated')

            #Set x ticks
            ax1.xaxis.set_major_locator(mdates.DayLocator())
            ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
            # Minor ticks: every 3 hours
            ax1.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21]))
            # Rotate and format tick labels
            ax1.tick_params(axis='x', which='major')
            ax1.tick_params(axis='x', which='minor', length=3, color='gray')
            ax1.set_xlim([xlim_start_w1, xlim_end_w1])
            ax1.legend(loc = 'upper right')

            ax2.plot(df_ml_data.index, df_ml_data[col], color='k', marker='o',label=f'Original (baddies={n_baddies})')
            ax2.plot(df_interp.index, df_interp[col], color='r', marker='x', label='Interpolated')

            # #Set x ticks
            ax2.xaxis.set_major_locator(mdates.DayLocator())
            ax2.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
            # Minor ticks: every 3 hours
            ax2.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21]))
            # Rotate and format tick labels
            ax2.tick_params(axis='x', which='major')
            ax2.tick_params(axis='x', which='minor', length=3, color='gray')
            ax2.set_xlim([xlim_start_w2, xlim_end_w2])
            ax2.legend(loc = 'upper right')

            ax3.plot(df_ml_data.index, df_ml_data[col], color='k', marker='o',label=f'Original (baddies={n_baddies})')
            ax3.plot(df_interp.index, df_interp[col], color='r', marker='x', label='Interpolated')
  
            # #Set x ticks
            ax3.xaxis.set_major_locator(mdates.DayLocator())
            ax3.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
            # Minor ticks: every 3 hours
            ax3.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21]))
            # Rotate and format tick labels
            ax3.tick_params(axis='x', which='major')
            ax3.tick_params(axis='x', which='minor', length=3, color='gray')
            ax3.set_xlim([xlim_start_w3, xlim_end_w3])
            ax3.legend(loc = 'upper right')

            ax4.plot(df_ml_data.index, df_ml_data[col], color='k', marker='o',label=f'Original (baddies={n_baddies})')
            ax4.plot(df_interp.index, df_interp[col], color='r', marker='x', label='Interpolated')

            # #Set x ticks
            ax4.xaxis.set_major_locator(mdates.DayLocator())
            ax4.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
            # Minor ticks: every 3 hours
            ax4.xaxis.set_minor_locator(mdates.HourLocator(byhour=[0, 3, 6, 9, 12, 15, 18, 21]))
            # Rotate and format tick labels
            ax4.tick_params(axis='x', which='major')
            ax4.tick_params(axis='x', which='minor', length=3, color='gray')
            ax4.set_xlim([xlim_start_w4, xlim_end_w4])
            ax4.legend(loc = 'upper right')

            plt.suptitle(col)
            plt.show()
#region: extra functions needed to export to matlab
def find_files(directory, search_string, ext):
    matches = []
    for root, dirs, files in os.walk(directory):
        for filename in files:
            if search_string in filename and filename.endswith(ext):
                matches.append(os.path.join(root, filename))
    return matches
def dataframe_to_nested_dict(df):
    nested_dict = {}
    for column in df.columns:
        nested_dict[column] = df[column][0:-1].to_numpy()
        
    # Add the index as a key-value pair
    nested_dict['time_local'] = df.index[0:-1].to_numpy()
    return nested_dict
#endregion
def subset_day(date_time_start, date_time_stop,file_subset_name, var_name):
    """
    This function is for if you want to interpolate a subset of the days for the campaign.

    INPUTS:
        date_time_start: A string that sets the beginning of your date and time subset, in format YYYY-MM-DD HH:MM:SS
        date_time_stop: A string that sets the end of your date and time subset, in format YYYY-MM-DD HH:MM:SS. 
                        NOTE: THIS IS AN INCLUSIVE VALUE so date_time_stop= "2024-08-08 23:30:00" includes the 23:30:00 value.
        file_subset_name: A string with format separated by underscore:
            date start
            date end
            time for averaging (such as parked data with 30 min averages)
            CSL_mobile_lab
            parked / driving
            with_interp (to indicate that it includes the interpolation)
        var_name: A string that sets the name of the MATLAB variable when you import the file into MATLAB.
    """
    df_subsetdays = df_ml_data.sort_index().loc[date_time_start:date_time_stop]
    df_interp_subset=df_subsetdays.copy()

    #Nell turned all benzaldehyde, styrene, HONO negative values into NaNs = 0
    neg_species = ['Benzaldehyde_PTR','Styrene_PTR', 'HONO_CIMS']

    for i,col in enumerate(vars2fill):
        if col in neg_species:
            df_interp_subset[col] = df_interp_subset[col].mask(df_interp_subset[col] < 0, 0)
        elif col == 'Lon': #avoids making longitude points into NaNs
            df_interp_subset['Lon'] = df_interp_subset['Lon'].interpolate(method='linear')
        else:
            # Set any negative values to NaN so we can interp them... 
            df_interp_subset[col] = df_interp_subset[col].mask(df_interp_subset[col] < 0, np.nan)
        #Benzaldehyde has all NaNs so substitute with zeros instead
        if col == 'Benzaldehyde_PTR':
            df_interp_subset[col] =  df_interp_subset[col].mask(np.isnan(df_interp_subset[col]), 0)
        else:
            pass
        # Calc number of points that are negative or Nans: 
        n_baddies= len([item for item in df_ml_data[col] if item <0 or np.isnan(item)]) 
        
        if n_baddies > 0: 
            #apply the linear interpolation
            df_interp_subset[col] = df_interp_subset[col].interpolate(method='linear')

    #get a ratio for jNO2 measured to TUV
    df_interp_subset['jNO2_ratio'] = df_interp_subset['jNO2_meas']/df_interp_subset['jNO2']
    #level out the inf values and values that are too high for the jNO2 ratio
    msk = ((df_interp_subset['jNO2_ratio'] ==np.inf)  | (df_interp_subset['jNO2_ratio'] >10) )
    df_interp_subset.loc[msk,'jNO2_ratio'] = 1.0

    savepath = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/CampaignData_and_Merges/R0/CSL_MobileLab_Parked/' + 'F0AM_filled/' + file_subset_name + '.csv'
    df_interp_subset.to_csv(savepath)
    print('Saved CSV to:' + savepath)

    # Convert the dataframe to a nested dictionary (so scipy can output to a matlab structure!) 
    ddict=dataframe_to_nested_dict(df_interp_subset)

    # Sort alphabetically so not annoying in MATLAB...  
    ddict= OrderedDict(sorted(ddict.items())) 

    # Save the USOS data in an output .mat file: 
    outpath = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/F0AM-4.4.2/Campaign_Data/matlab_merge/parked/original/'
    matfilename = file_subset_name + '.mat'
    savemat(outpath+matfilename,{var_name: ddict})
    print('Saved MATLAB file to:' + outpath + matfilename)

def all_days(file_name, var_name):
    df_interp_alldays=df_ml_data.copy()

    #Nell turned all benzaldehyde, styrene, HONO negative values into NaNs = 0
    neg_species = ['Benzaldehyde_PTR','Styrene_PTR', 'HONO_CIMS']

    for i,col in enumerate(vars2fill):
        if col in neg_species:
            df_interp_alldays[col] = df_interp_alldays[col].mask(df_interp_alldays[col] < 0, 0)
        elif col == 'Lon':
            df_interp_alldays['Lon'] = df_interp_alldays['Lon'].interpolate(method='linear')
        else:
            # Set any negative values to NaN so we can interp them... 
            df_interp_alldays[col] = df_interp_alldays[col].mask(df_interp_alldays[col] < 0, np.nan)
        #Benzaldehyde has all NaNs so substitute with zeros instead
        if col == 'Benzaldehyde_PTR':
            df_interp_alldays[col] =  df_interp_alldays[col].mask(np.isnan(df_interp_alldays[col]), 0)
        else:
            pass
        # Calc number of points that are negative or Nans: 
        n_baddies= len([item for item in df_ml_data[col] if item <0 or np.isnan(item)]) 
        
        if n_baddies > 0: 
            #apply the linear interpolation
            df_interp_alldays[col] = df_interp_alldays[col].interpolate(method='linear')
        
        else:
            pass

    #get a ratio for jNO2 measured to TUV
    df_interp_alldays['jNO2_ratio'] = df_interp_alldays['jNO2_meas']/df_interp_alldays['jNO2']
    #level out the inf values and values that are too high for the jNO2 ratio
    msk = ((df_interp_alldays['jNO2_ratio'] ==np.inf)  | (df_interp_alldays['jNO2_ratio'] >10) )
    df_interp_alldays.loc[msk,'jNO2_ratio'] = 1.0

    print(df_interp_alldays['Lon'].values)
    savepath = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/CampaignData_and_Merges/R0/CSL_MobileLab_Parked/' + 'F0AM_filled/' + file_name + '.csv'
    df_interp_alldays.to_csv(savepath)
    print('Saved CSV to:' + savepath)

    # Convert the dataframe to a nested dictionary (so scipy can output to a matlab structure!) 
    ddict=dataframe_to_nested_dict(df_interp_alldays)

    # Sort alphabetically so not annoying in MATLAB...  
    ddict= OrderedDict(sorted(ddict.items())) 

    # Save the USOS data in an output .mat file: 
    outpath = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/F0AM-4.4.2/Campaign_Data/matlab_merge/parked/original/'
    matfilename = file_name + '.mat'
    savemat(outpath+matfilename,{var_name: ddict})
    print('Saved MATLAB file to:' + outpath + matfilename)

###### CALL FUNCTIONS #########
# plot_interps()

# start_month = '08'
# start_day = '05'
# end_month = '08'
# end_day = '07'
# # subset_day(
# #      date_time_start = '2024-' + start_month + '-' + start_day + ' 00:00:00', 
# #      date_time_stop = '2024-' + end_month + '-' + end_day + ' 23:45:00',
# #      file_subset_name = '2024' + start_month + start_day +'_' + '2024' + end_month + end_day + '_15min_CSL_mobile_lab_parked_with_interp_nell_match_with_formaldehyde_raw_noadjustment', 
# #      var_name = 'USOS'
# # )

all_days(
    file_name = 'alldays_15min_CSL_mobile_lab_parked_with_interp_nell_match_with_formaldehyde_raw_noadjustment',
    var_name = 'USOS'
)