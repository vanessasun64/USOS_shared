import numpy as np 
import os 
import xarray as xr
import pandas as pd
#from scipy.io import savemat
from collections import OrderedDict
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
import matplotlib as mpl

def single_day_full_vars_availability_plot(date_to_plot):
    file_path='/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_merges/R0/CSL_MobileLab_Parked/merged/rev_30min/USOS_rev30min_Merged_Parked_2024'+ str(date_to_plot) +'.nc'

    # Load in the USOS 30min parked data on 1 day: 
    ds=xr.open_dataset(file_path)

    # Convert the Dataset to a DataFrame
    df_all = ds.to_dataframe()

    # Setting local 'time' as the index 
    df_all.set_index('time_local', inplace=True)

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
            'Benzaldehyde_PTR':True,
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
            'AOD':False,
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
            'Time_Start_POPS':False,
            'Time_Stop_POPS':False,
            'Press_mb_POPS':False,
            'Temp_C_POPS':False,
            'Aerosol_conc_POPS':False,
            'Surface_area_conc_POPS':False,
            'Volume_density_POPS':False,
            'Effective_radius_POPS':False,
            'Extinction_POPS':False,
            'Particle_conc_bin01_POPS':False,
            'Particle_conc_bin02_POPS':False,
            'Particle_conc_bin03_POPS':False,
            'Particle_conc_bin04_POPS':False,
            'Particle_conc_bin05_POPS':False,
            'Particle_conc_bin06_POPS':False,
            'Particle_conc_bin07_POPS':False,
            'Particle_conc_bin08_POPS':False,
            'Particle_conc_bin09_POPS':False,
            'Particle_conc_bin10_POPS':False,
            'Particle_conc_bin11_POPS':False,
            'Particle_conc_bin12_POPS':False,
            'Particle_conc_bin13_POPS':False,
            'Particle_conc_bin14_POPS':False,
            'Particle_conc_bin15_POPS':False,
            'Particle_conc_bin16_POPS':False,
            'Particle_conc_bin17_POPS':False,
            'Particle_conc_bin18_POPS':False,
            'Particle_conc_bin19_POPS':False,
            'Particle_conc_bin20_POPS':False,
            'Particle_conc_bin21_POPS':False,
            'Particle_conc_bin22_POPS':False,
            'Particle_conc_bin23_POPS':False,
            'Particle_conc_bin24_POPS':False,
            'Particle_conc_bin25_POPS':False,
            'Particle_conc_bin26_POPS':False,
            'Particle_conc_bin27_POPS':False,
            'Particle_conc_bin28_POPS':False,
            'Particle_conc_bin29_POPS':False,
            'Particle_conc_bin30_POPS':False,
            'Particle_conc_bin31_POPS':False,
            'Particle_conc_bin32_POPS':False,
            'Particle_conc_bin33_POPS':False,
            'Particle_conc_bin34_POPS':False,
            'Particle_conc_bin35_POPS':False,
            'Particle_conc_bin36_POPS':False}

    # Get names of vars we need to fill nans in: 
    vars2fill=[key for key,value in need2fill.items() if value ==True]

    df_foam=pd.DataFrame(index=df_all.index)

    for col in vars2fill:
        if col in df_all.columns:
        # Everything that was a NaN or negative is assigned a 0, positive values are assigned a 1
            df_foam[col] = np.where((df_all[col] > 0) & (~df_all[col].isna()), 1, 0 ) 
        else:
            df_foam[col]= 0

    df_transposed=df_foam.T

    # Create the heatmap using imshow
    plt.figure(figsize=(10, 10))
    plt.imshow(df_transposed, aspect='auto', cmap='viridis')

    # Set x-ticks and y-ticks
    plt.xticks(ticks=np.arange(len(df_transposed.columns)), labels=df_transposed.columns.strftime('%H:%M:%S'))
    plt.yticks(ticks=np.arange(len(df_transposed.index)), labels=df_transposed.index)

    # Add color bar
    plt.colorbar(label='Value')

    # Set labels and title
    plt.xlabel('Time')
    plt.ylabel('Series')
    plt.title('Data availability for: ' + str(date_to_plot))

    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45)

    # Show plot
    plt.tight_layout()
    #plt.show()

    save_filepath_merger = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_merges/data_availability_plots/full_var_availability_'+ str(date_to_plot) + '.png'
    plt.savefig(fname = save_filepath_merger,format='png')
    print('Saved figure as:' + str(save_filepath_merger))

def three_day_important_availability_plot(dates_to_plot = []):
    all_ds=[]
    for date_parked in dates_to_plot:
        file_path='/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_merges/R0/CSL_MobileLab_Parked/merged/rev_30min/USOS_rev30min_Merged_Parked_2024'+ str(date_parked) +'.nc'

    # Load in the USOS 30min parked data on 1 day: 
        ds=xr.open_dataset(file_path)
        all_ds.append(ds)
    combo_ds = xr.concat(all_ds, dim='time_UTC')

    # Convert the Dataset to a DataFrame
    df_all = combo_ds.to_dataframe()

    # Setting local 'time' as the index 
    df_all.set_index('time_local', inplace=True)

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
            'Benzaldehyde_PTR':True,
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
            'AOD':False,
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
            'Time_Start_POPS':False,
            'Time_Stop_POPS':False,
            'Press_mb_POPS':False,
            'Temp_C_POPS':False,
            'Aerosol_conc_POPS':False,
            'Surface_area_conc_POPS':False,
            'Volume_density_POPS':False,
            'Effective_radius_POPS':False,
            'Extinction_POPS':False,
            'Particle_conc_bin01_POPS':False,
            'Particle_conc_bin02_POPS':False,
            'Particle_conc_bin03_POPS':False,
            'Particle_conc_bin04_POPS':False,
            'Particle_conc_bin05_POPS':False,
            'Particle_conc_bin06_POPS':False,
            'Particle_conc_bin07_POPS':False,
            'Particle_conc_bin08_POPS':False,
            'Particle_conc_bin09_POPS':False,
            'Particle_conc_bin10_POPS':False,
            'Particle_conc_bin11_POPS':False,
            'Particle_conc_bin12_POPS':False,
            'Particle_conc_bin13_POPS':False,
            'Particle_conc_bin14_POPS':False,
            'Particle_conc_bin15_POPS':False,
            'Particle_conc_bin16_POPS':False,
            'Particle_conc_bin17_POPS':False,
            'Particle_conc_bin18_POPS':False,
            'Particle_conc_bin19_POPS':False,
            'Particle_conc_bin20_POPS':False,
            'Particle_conc_bin21_POPS':False,
            'Particle_conc_bin22_POPS':False,
            'Particle_conc_bin23_POPS':False,
            'Particle_conc_bin24_POPS':False,
            'Particle_conc_bin25_POPS':False,
            'Particle_conc_bin26_POPS':False,
            'Particle_conc_bin27_POPS':False,
            'Particle_conc_bin28_POPS':False,
            'Particle_conc_bin29_POPS':False,
            'Particle_conc_bin30_POPS':False,
            'Particle_conc_bin31_POPS':False,
            'Particle_conc_bin32_POPS':False,
            'Particle_conc_bin33_POPS':False,
            'Particle_conc_bin34_POPS':False,
            'Particle_conc_bin35_POPS':False,
            'Particle_conc_bin36_POPS':False}

    # Get names of vars we need to fill nans in: 
    vars2fill=[key for key,value in need2fill.items() if value ==True]

    df_foam=pd.DataFrame(index=df_all.index)

    for col in vars2fill:
        if col in df_all.columns:
        # Everything that was a NaN or negative is assigned a 0, positive values are assigned a 1
            df_foam[col] = np.where((df_all[col] > 0) & (~df_all[col].isna()), 1, 0 ) 
        else:
            df_foam[col]= 0

    df_transposed=df_foam.T

    # Create the heatmap using imshow
    plt.figure(figsize=(30, 10))
    plt.imshow(df_transposed, aspect='auto', cmap='viridis')

    # Set x-ticks and y-ticks
    plt.xticks(ticks=np.arange(len(df_transposed.columns)), labels=df_transposed.columns.strftime('%H:%M:%S'))
    plt.yticks(ticks=np.arange(len(df_transposed.index)), labels=df_transposed.index)

    # Add color bar
    plt.colorbar(label='Value')

    # Set labels and title
    plt.xlabel('Time')
    plt.ylabel('Series')
    plt.title('Data availability for: ' + str(dates_to_plot))

    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45)

    # Show plot
    plt.tight_layout()
    #plt.show()

    #save_filepath_merger = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_merges/data_availability_plots/three_day_availability_'+ str(date_to_plot) + '.png'
    #plt.savefig(fname = save_filepath_merger,format='png')
    #print('Saved figure as:' + str(save_filepath_merger))
# parked_dates_list = ['0715','0716','0717','0718','0719','0720','0721','0722',
# '0723','0724','0725','0726','0727','0728','0729','0730','0731',
# '0801','0802','0803','0804','0805','0806','0807','0808',
# '0809','0810','0811','0812','0813','0814','0815','0816','0817','0818']

# for date_parked in parked_dates_list:
#     availability_plot(date_parked)

smoke_free_days = ['0804','0805','0806']
three_day_important_availability_plot(smoke_free_days)