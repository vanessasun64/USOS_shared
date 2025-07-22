import numpy as np 
import os 
import xarray as xr
import pandas as pd
from scipy.io import savemat
from collections import OrderedDict

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

def subset_days(date_time_start, date_time_stop,file_subset_name, var_name):
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

    df_subsetdays = df_alldays.sort_index().loc[date_time_start:date_time_stop]

    df_interp_subset=df_subsetdays.copy()
    for i,col in enumerate(vars2fill):
        # Set any negative values to NaN so we can iterp them... 
        df_interp_subset[col] = df_interp_subset[col].mask(df_interp_subset[col] < 0, np.nan)
                                
        # Calc number of points that are negative or Nans: 
        n_baddies= len([item for item in df_subsetdays[col] if item <0 or np.isnan(item)])
                                
        if n_baddies > 0: 
                #apply the linear interpolation
            df_interp_subset[col] = df_interp_subset[col].interpolate(method='linear')

    # Convert the dataframe to a nested dictionary (so scipy can output to a matlab structure!) 
    ddict=dataframe_to_nested_dict(df_interp_subset)

    # Sort alphabetically so not annoying in MATLAB...  
    ddict= OrderedDict(sorted(ddict.items())) 

    # Save the USOS data in an output .mat file: 
    outpath = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/F0AM-4.3.0.1/Campaign_Data/matlab_merge/parked/original/'
    matfilename = file_subset_name + '.mat'
    savemat(outpath+matfilename,{var_name: ddict})
    print('Saved MATLAB file to:' + outpath + matfilename)

def all_days(file_alldays_name):
    fixed_ct=0; 

    df_interp_alldays=df_alldays.copy()
    for i,col in enumerate(vars2fill):
        # Set any negative values to NaN so we can iterp them... 
        df_interp_alldays[col] = df_interp_alldays[col].mask(df_interp_alldays[col] < 0, np.nan)
                                
        # Calc number of points that are negative or Nans: 
        n_baddies= len([item for item in df_alldays[col] if item <0 or np.isnan(item)])
                                
        if n_baddies > 0: 
            if fixed_ct<10: 
                #apply the linear interpolation
                df_interp_alldays[col] = df_interp_alldays[col].interpolate(method='linear')
            fixed_ct=fixed_ct+1

    # Convert the dataframe to a nested dictionary (so scipy can output to a matlab structure!) 
    ddict=dataframe_to_nested_dict(df_interp_alldays)

    # Sort alphabetically so not annoying in MATLAB...  
    ddict= OrderedDict(sorted(ddict.items())) 

    # Save the USOS data in an output .mat file: 
    outpath = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/F0AM-4.3.0.1/Setups/'
    matfilename = file_alldays_name + '.mat'
    savemat(outpath+matfilename,{'dat': ddict})
    print('Saved MATLAB file to:' + outpath + matfilename)

# %% BEGIN MAIN 

#import NetCDF file for all days
all_days_filepath = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/CampaignData_and_Merges/R0/CSL_MobileLab_Parked/merged/rev_30min/all_CSL_MobileLab_Parked_rev30minv4.nc'
all_days_filepath_load = xr.open_dataset(all_days_filepath)
df_alldays = all_days_filepath_load.to_dataframe()
df_alldays.reset_index(inplace=True)
df_alldays.set_index('time_local', inplace=True)

# Define which variables we need to make sure don't have Nans/ negs since we'll be using then as constraints in F0AM: 
# True means constrain
need2fill= {'Br2_CIMS':True,
        'BrO_CIMS':True,
        'Cl2_CIMS':True,
        'ClNO2_CIMS':True,
        'BrCl_CIMS':True,
        'NCl3_CIMS':False,
        'File_Index':False,
        'NO_LIF':True,
        'NO2_LIF':True,
        'NOy_LIF':False,
        'ACCURACY_1sigma_NO':False,
        'ACCURACY_1sigma_NO2':False,
        'ACCURACY_1sigma_NOy':False,
        'N2O5_CIMS':True,
        'C4H7NO5_CIMS':False,
        'ISOPN_CIMS':False,
        'C5H10O3_CIMS':False,
        'C5H9NO5_CIMS':False,
        'C10H17NO4_CIMS':False,
        'HNO3_CIMS':True,
        'Methanol_PTR':True,
        'Acetonitrile_PTR':True,
        'Acetaldehyde_PTR':True,
        'Ethanol_PTR':True,
        'Methanethiol_PTR':False,
        'Acrolein_PTR':True,
        'Acetone_Propanal_PTR':False,
        'DMS_PTR':True,
        'Isoprene_PTR':True,
        'MVK_MACR_PTR':True,
        'Benzene_PTR':True,
        'Toluene_PTR':True,
        'Styrene_PTR':True,
        'Benzaldehyde_PTR':False,
        'C8Aromatics_PTR':False,
        'C9Aromatics_PTR':False,
        'Naphthalene_PTR':False,
        'Octanal_PTR':False,
        'Monoterpenes_PTR':True,
        'Nonanal_PTR':False,
        'C7H4ClF3_PTR':False,
        'D5_siloxane_PTR':False,
        'PAN_CIMS':True,
        'APAN_CIMS':False,
        'PPN_CIMS':True,
        'HCOOH_CIMS':True,
        'HONO_CIMS':True,
        'O3_ppbv':True,
        'Altitude_m':False,
        'Course_deg':False,
        'GndSpd_ms':False,
        'Heading_deg':False,
        'Lat':False,
        'Lon':False,
        'Temp_K':True,
        'Pressure_mb':True,
        'WindDir_deg':False,
        'WindSpd_ms':True,
        'jNO2_meas':True,
        'CO2_Piccaro':False,
        'CH4_Piccaro':True,
        'CO_Piccaro':True,
        'H2O_Piccaro':False,
        'HCHO_CRDS':True,
        'H2O_CRDS':True,
        'CH4_CRDS':True,
        'Time_Start_WAS':False,
        'Time_Stop_WAS':False,
        'Time_Mid_WAS':False,
        'Acetone_WAS':False,
        'Acrolein_WAS':False,
        'Benzene_WAS':False,
        'C2Cl4_WAS':False,
        'C2HCl3_WAS':False,
        'CCl4_WAS':False,
        'CF2Cl2_WAS':False,
        'CFCl3_WAS':False,
        'CH2Cl2_WAS':False,
        'CH3Br_WAS':False,
        'CycloPentane_WAS':False,
        'Ethyl_WAS':False,
        'Ethane_WAS':False,
        'Ethene_WAS':False,
        'Ethyne_WAS':False,
        'Furan_WAS':False,
        'Isoprene_WAS':False,
        'Limonene_WAS':False,
        'MACR_WAS':False,
        'MethylCycloHexane_WAS':False,
        'MethylCycloPentane_WAS':False,
        'Propane_WAS':False,
        'Propene_WAS':False,
        'Toluene_WAS':False,
        'Alpha_Pinene_WAS':False,
        'Beta_Pinene_WAS':False,
        'c2Butene_WAS':False,
        'c2Pentene_WAS':False,
        'iButane_WAS':False,
        'iPentane_WAS':False,
        'iPropylBenzene_WAS':False,
        'iPropylONO2_WAS':False,
        'm_p_Xylene_WAS':False,
        'nButane_WAS':False,
        'nDecane_WAS':False,
        'nHeptane_WAS':False,
        'nHexane_WAS':False,
        'nNonane_WAS':False,
        'nOctane_WAS':False,
        'nPentane_WAS':False,
        'nPropylBenzene_WAS':False,
        'nPropylONO2_WAS':False,
        'o_Xylene_WAS':False,
        't2Butene_WAS':False,
        't2Pentene_WAS':False,
        'x123_TriMethylBenzene_WAS':False,
        'x124_TriMethylBenzene_WAS':False,
        'x135_TriMethylBenzene_WAS':False,
        'x13_Butadiene_WAS':False,
        'x1_Butene_WAS':False,
        'x1_Pentene_WAS':False,
        'x224_TriMethylPentane_WAS':False,
        'x22_DiMethylButane_WAS':False,
        'x2_EthylToluene_WAS':False,
        'x2_MethylPentane_WAS':False,
        'x3_MethylPentane_WAS':False,
        'x3_x4_EthylToluene_WAS':False,
        'RH_percent':True,
        'AOD':False}
        # 'jNO2':True,
        # 'jBrCl':True,
        # 'jBr2':True,
        # 'jCCl4':True,
        # 'jCH2Oa':True,
        # 'jCH2Ob':True,
        # 'jClNO2':True,
        # 'jClOa':True,
        # 'jClOb':True,
        # 'jCl2':True,
        # 'jHNO2':True,
        # 'jHNO3':True,
        # 'jI2':True,
        # 'jNO3a':True,
        # 'jNO3b':True,
        # 'jN2O5':True,
        # 'jO3':True,
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

#07/07/2025 Subset for smokefree days for first F0AM run, with varname set to Campaign Name.
subset_days(date_time_start = "2024-08-04 00:00:00", 
            date_time_stop= "2024-08-08 23:00:00", 
            file_subset_name ='20240804_20240808_30min_CSL_mobile_lab_parked_with_interp_with_pan_interp',
            var_name = 'USOS')

# subset_days(date_time_start = "2024-08-04 00:00:00", 
#             date_time_stop= "2024-08-08 23:30:00", 
#             file_subset_name ='20240804_20240808_30min_CSL_mobile_lab_parked_with_interp',
#             var_name = 'USOS')
# all_days(file_alldays_name = "30min_USOS_CSL_mobile_lab_parked_with_interp")