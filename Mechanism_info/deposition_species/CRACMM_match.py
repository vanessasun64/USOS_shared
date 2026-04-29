import pandas as pd # ok with rdkit 
import pickle
import numpy as np
# import ast
import sys
from rdkit import Chem

from rdkit.Chem import  Descriptors, rdMolDescriptors, Fragments
from scipy.io import savemat
from collections import OrderedDict, defaultdict
import csv
import re

dirpath = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/'


smart_groups= dict({
    'Is_Radical':	'[#8H0;v1+0]',
    'RO2s':	'[#6,S;!$([#6](=[#8]))]~[#8][#8;X1v1+0]',
    'ROs':	'[#6,S,c;!$([#6](=[#8]))]~[$([#8X2+0]),$([#8;X1v1+0]);!$([#8X2H1]);!$([#8X2][#8]);!$([#8X2][#7]);!$([#8X2]([#6])[#6]);!$([#8X2][#16]);!$([#8v+](=[#6])[#8v-])]',
    'Acyl_RO2s':	'[#6](=[#8])~[#8][#8;X1v1+0]',
    'Acyl_ROs':	'[#6](=[#8])[$([#8X2]),$([#8;X1v1+0]);!$([#8X2H1]);!$([#8X2][#8]);!$([#8X2][#7]);!$([#8X2]([#6])[#6]);!$([#8X2][#16])]',
    'All_NO3s':	'[$([#7X3](=[#8X1])(=[#8X1])[#8]),$([#7X3+]([#8X1-])(=[#8X1])[#8]),$([#7X3+](=[#8X1][#8X1-])[#8])]',
    'Non_PAN_NO3s':	'[$([#7X3](=[#8X1])(=[#8X1])[#8;!$([#8]~[#8])]),$([#7X3+]([#8X1-])(=[#8X1])[#8;!$([#8]~[#8])])]',
    'Tertiary_NO3s':	'[$([#7X3](=[#8X1])(=[#8X1])([#8]-[#6X4H0])),$([#7X3+]([#8X1-])(=[#8X1])([#8]-[#6X4H0]))]',
    'RONO2s':	'[$([#7X3](=[#8X1])(=[#8X1])[#8]~[*;!$([#8])]),$([#7X3+]([#8X1-])(=[#8X1])[#8]~[*;!$([#8])])]',
    'RO2NO2s':	'[$([#7X3](=[#8X1])(=[#8X1])[#8]-[#8][#6;!$([#6]=[#8])]),$([#7X3+]([#8X1-])(=[#8X1])[#8]-[#8][#6;!$([#6]=[#8])])]',
    'PANs':	'[$([#7X3](=[#8X1])(=[#8X1])[#8]-[#8]-[#6](=[#8])),$([#7X3+]([#8X1-])(=[#8X1])[#8]-[#8]-[#6](=[#8]))]',
    'Nitros':	'[#6][$([#7X3](=[#8])=[#8]),$([#7X3+](=[#8])[#8-])](~[#8])(~[#8])',
    'Peroxides':	'[#8X2;!$([#8]-[#6](=[#8]));!$([#8]-[#7])][$([#8X2H1v2+0]),$([#8X2;!$([#8]-[#6](=[#8]));!$([#8]-[#7])])]',
    'OrganicPeroxides':	'[#6X4;!$([#6](=[#8]))][#8X2H0][#8X2H0][#6X4]',
    'HydroPeroxides':	'[#8X2;!$([#8]-[#6](=[#8]));!$([#8]-[#7])][#8X2H1v2+0]',
    'Peracids':	'[*;!$([#8]);!$([#7]);!$([#16])](=[#8])-[#8]-[#8X2H1v2+0]',
    'Carbonyls':	'[$([#6X3]=[#8X1]),$([#6X3+]-[#8X1-])]',
    'Carboxylic_Acids':	'[#6X3;!$([#6X3H])](=[#8])[#8X2H]',
    'Ketones':	'[#6X3;!$([#6X3H1]);!$([#6X3H2]);!$([#6]-[#8])](=[#8])',
    'Aldehydes':	'[#6X3H1;!$([#6X3H1]-[#8X2])](=[#8;!$([#8][#8])])',
    'Esters':	'[$([#6H1](=[#8])(-[#8][#6])),$([#6X3](=[#8])(-[#8][#6])[*;!$([#8])]),$([#6X3](=[#8])(-[#8][#6])[#8])]',
    'Ethers':	'[#8D2]([#6;!$([#6]=[#8])])[#6;!$([#6]=[#8])]',
    'Carbonates':	'[#6X3](=[#8X1])(-[$([#8X2H0]),$([#8X1-1])])(-[$([#8X2H0]),$([#8X1-1])])',
    'All_OHs':	'[#6]-[#8X2H1]',
    'Hydroxyls':	'[$([#8X2H]);!$([#8X2H][#6X3](=[#8]));!$([#8X2H]-[$([#8]),$([#16](=[#8])(=[#8]))])]',
    'Enols':	'[$([#6](=[#6])(-[#8X2H1]));!$([$([c](-[#8X2H1]));!$([c](=[#8])(-[#8X2H1]))])]',
    'Phenols': '[$([c](-[#8X2H1]));!$([c](=[#8])(-[#8X2H1]))]', 
    'Thiols':	'[#6]-[#16X2H1]',
    'Dihydroxys':	'[#6X4;!$([#6]=[#8])](-[#8X2H1])(-[#8X2H1])',
    'Aliphatic_Alcohols':	'[$([CX4R0](-[#8X2H1]));!$([#6X3;!$([#6X3H])](=[#8])[#8X2H])!$([#6X4;!$([#6]=[#8])](-[#8X2H1])(-[#8X2H1]))]',
    'Other_Alcohols':'[$([#6;!$([$([#6](=[#6])(-[#8X2H1]));!$([$([c](-[#8X2H1]));!$([c](=[#8])(-[#8X2H1]))])]);!$([#6]-[#16X2H1]);!$([#6X4;!$([#6]=[#8])](-[#8X2H1])(-[#8X2H1]))][#8X2H1]);!$([$([c](-[#8X2H1]));!$([c](=[#8])(-[#8X2H1]))]);!$([#6X3;!$([#6X3H])](=[#8])[#8X2H]);!$([$([CX4R0](-[#8X2H1]));!$([#6X3;!$([#6X3H])](=[#8])[#8X2H])!$([#6X4;!$([#6]=[#8])](-[#8X2H1])(-[#8X2H1]))])]',
    #Not enol, Not carbon thiol, not a dihydroxy OH carbon, But that carbon does gotta be attached to an OH group somewhere, no phenol, no carb acid, and no aliphatic alcohols
    'C':	'[#6]',
    'H':	'[#1]',
    'O':	'[#8]',
    'N':	'[#7]',
    'S':	'[#16]',
    'Cl':	'[#17]',
    'Br':	'[#35]',
    'Toluene': '[#6;$([c]1ccccc1[CH3])]',
    'Benzene': '[c;r6]1[c;r6][c;r6][c;r6][c;r6][c;r6]1',
    'Methyl': '[CH3]',
    'Cresol': '[c]-[CH3]'
})
def group_all_species(df_in, smart_groups):
    df_out=pd.DataFrame()
    # df_out[use]= df_in[use].copy()
    df_out['Species'] = df_in['Species']
    df_out['Description'] = df_in['Description']
    df_out['Phase'] = df_in['Phase']
    df_out['Stable'] = df_in['Stable']
    df_out['Molecular Weight (g/mol)'] = df_in['Molecular Weight (g/mol)']
    df_out['Explicit/Lumped'] = df_in['Explicit/Lumped']
    df_out['Representative'] = df_in['Representative']
    df_out['SMILES'] = df_in['SMILES']
    df_out['DTXSID'] = df_in['DTXSID']
    df_out['H Law (M/atm)'] = df_in['H Law (M/atm)']
    df_out['Enthalpy of solution (K)'] = df_in['Enthalpy of solution (K)']
    df_out['Aerosol density (kg/m3)'] = df_in['Aerosol density (kg/m3)']
    df_out['Kappa_org'] = df_in['Kappa_org']
    df_out['C* (microg/m3)'] = df_in['C* (microg/m3)']
    df_out['Enthalpy of vaporization (J/mol)'] = df_in['Enthalpy of vaporization (J/mol)']
    df_out['OM to OC (g/g)'] = df_in['OM to OC (g/g)']

    df_out['Epoxides']=list(np.full([len(df_out),1], np.nan)) 
    for i, species in enumerate(df_allspecies['Species']):
        if df_in.loc[i, 'SMILES'] is not np.nan:
            print(i, ' Species: ', df_allspecies.loc[i,'Species']) # Print off name of what you're parsing. 
            print('SMILES: ', df_in.loc[i, 'SMILES'])
            df_out.at[i, 'SMILES'] = df_in.loc[i, 'SMILES']
            molec = Chem.MolFromSmiles(df_in.loc[i, 'SMILES'])  #Turn this molec it into an RDKit molecule object. 
            #molec = Chem.AddHs(molec)
            #molec = Chem.SetAromaticity(molec)
            if molec is not None:
                molec.UpdatePropertyCache(strict=True)  # for radicals!
                for key in smart_groups:  # Loop over every functional group you want to search for. 
                    print('Obtaining info on Functional Group ', key, '\n')
                    # Turn the SMARTs string for this functional group into a RDKit molec fragment. 
                    frag = Chem.MolFromSmarts(smart_groups[key])
                    
                    # Get a list of the indices of atom #s in molecule that match this fragment 
                    inds=list(molec.GetSubstructMatches(frag))
                    
                    # Save the len of this list as the # of functional group matches you found!)
                    df_out.at[i,key]=np.int64(len(inds))
                
                    # # Use RDKit to get fragments(not always as specific as our group matches...) 
                    # rd_frags=get_rdkit_frags(df_in, use, molec)

                df_out.at[i, 'Epoxides'] = Fragments.fr_epoxide(molec) # Number of epoxide rings 
        else:
            print('No SMILES for ', species)
            df_out.at[i, 'SMILES'] = np.nan
            df_out.at[i, 'Is_Radical'] = np.nan
            df_out.at[i, 'RO2s'] = np.nan
            df_out.at[i, 'ROs'] = np.nan
            df_out.at[i, 'Acyl_RO2s'] = np.nan
            df_out.at[i, 'Acyl_ROs'] = np.nan
            df_out.at[i, 'All_NO3s'] = np.nan
            df_out.at[i, 'Non_PAN_NO3s'] = np.nan
            df_out.at[i, 'Tertiary_NO3s'] = np.nan
            df_out.at[i, 'RONO2s'] = np.nan
            df_out.at[i, 'RO2NO2s'] = np.nan
            df_out.at[i, 'PANs'] = np.nan
            df_out.at[i, 'Nitros'] = np.nan
            df_out.at[i, 'Peroxides'] = np.nan
            df_out.at[i, 'OrganicPeroxides'] = np.nan
            df_out.at[i, 'HydroPeroxides'] = np.nan
            df_out.at[i, 'Peracids'] = np.nan
            df_out.at[i, 'Carbonyls'] = np.nan
            df_out.at[i, 'Carboxylic_Acids'] = np.nan
            df_out.at[i, 'Ketones'] = np.nan
            df_out.at[i, 'Aldehydes'] = np.nan
            df_out.at[i, 'Esters'] = np.nan
            df_out.at[i, 'Ethers'] = np.nan
            df_out.at[i, 'Carbonates'] = np.nan
            df_out.at[i, 'All_OHs'] = np.nan
            df_out.at[i, 'Hydroxyls'] = np.nan
            df_out.at[i, 'Enols'] = np.nan
            df_out.at[i, 'Phenols'] = np.nan
            df_out.at[i, 'Thiols'] = np.nan
            df_out.at[i, 'Dihydroxys'] = np.nan
            df_out.at[i, 'Aliphatic_Alcohols'] = np.nan
            df_out.at[i, 'Other_Alcohols'] = np.nan
            df_out.at[i, 'C'] = np.nan
            df_out.at[i, 'H'] = np.nan
            df_out.at[i, 'O'] = np.nan
            df_out.at[i, 'N'] = np.nan
            df_out.at[i, 'S'] = np.nan
            df_out.at[i, 'Cl'] = np.nan
            df_out.at[i, 'Br'] = np.nan
            df_out.at[i, 'Epoxides'] = np.nan
            df_out.at[i, 'Toluene'] = np.nan
            df_out.at[i, 'Benzene'] = np.nan
            df_out.at[i, 'Methyl'] = np.nan
            df_out.at[i, 'Cresol'] = np.nan
    # Add a column that has the number of OH groups on a compounds that might cause it to be an organic acid... 
    Organic_Acid_OHs=['Enols', 'Phenols', 'Thiols', 'Carboxylic_Acids']
    df_out['Organic_Acid_OHs']=df_out[Organic_Acid_OHs].sum(axis=1)

    #Get all the information in Excel file
    print('df_merged_total: ', df_out)
    df_out.to_excel('cracmm_allspecies_info_rdkit.xlsx')
    print('Saved spreadsheet for all CRACMM RDKit info: ', 'cracmm_allspecies_info_rdkit.xlsx')

def choose_smarts_classifications():
    # df_updated_info = pd.read_excel(dirpath + 'Mechanism_info/geoschem_allspecies_info_rdkit.xlsx', index_col=0)
    df_updated_info = pd.read_excel('cracmm_allspecies_info_rdkit.xlsx', index_col=0)

    #Get species for groupings needed for deposition files in MCM.
    #These are RONO2s, carboxylic acids, hydroperoxides, ovocs

    alkyl_nitrates_excluding_peroxy_nitrates_and_pans = df_updated_info.loc[(df_updated_info['RONO2s'] > 0) & (df_updated_info['PANs'] == 0) & (df_updated_info['Is_Radical'] == 0)]
    carboxylic_acids = df_updated_info.loc[(df_updated_info['Carboxylic_Acids'] > 0)]
    hydroperoxides = df_updated_info.loc[(df_updated_info['HydroPeroxides'] > 0)]

    #Alfie counts any carbonyls and any alcohols as OVOCs. 
    #Our carbonyls functional SMARTs includes all carbonyls including carboxylic acids. Since we have a separate depositional velocity, we want all the carbonyls excluding carboxylic acids.
    ovocs = df_updated_info.loc[((df_updated_info['Carbonyls'] > 0) & (df_updated_info['Carboxylic_Acids'] == 0)) | ((df_updated_info['Dihydroxys'] > 0)) | ((df_updated_info['Aliphatic_Alcohols'] > 0)) | ((df_updated_info['Other_Alcohols'] > 0))]

    # Check for duplicates in the OVOC list from carbonyls, dihidroxys, aliphatic alcohols, and other alcohols setup
    ovocs_arr = ovocs['Species'].values
    print('OVOCs: ', ovocs_arr)

    u, c = np.unique(ovocs_arr, return_counts=True)
    duplicates_carbonyls_alcohols = u[c > 1]
    print('duplicates_carbonyls_alcohols: ', duplicates_carbonyls_alcohols)
    #We have no duplicates from the carbonyls and alcohols list. 

    # Now we want to check if there are any duplicates in multifunctional species. We get arrays of all the species in each grouping.
    rono2_arr = alkyl_nitrates_excluding_peroxy_nitrates_and_pans['Species'].values
    carboxylic_acids_arr = carboxylic_acids['Species'].values
    hydroperoxides_arr = hydroperoxides['Species'].values

    #First we tackle duplicates related to the deposition, keeping the duplicate in the functional group with the largest depositional velocity
    arrays_to_check_duplicates_deposition = [rono2_arr, carboxylic_acids_arr, hydroperoxides_arr, ovocs_arr]
    arrays_to_check_duplicates_deposition_names = ['rono2_arr', 'carboxylic_acids_arr', 'hydroperoxides_arr', 'ovocs_arr']
    #Order of depositional velocities: RONO2, carboxylic acids, hydroperoxides, ovocs
    deposition_velocities = [2.0, 1.0, 1.8, 1.2]

    # Step 1: map value -> list of array names
    locations = defaultdict(list)

    for name, arr in zip(arrays_to_check_duplicates_deposition_names, arrays_to_check_duplicates_deposition):
        for val in arr:
            locations[val].append(name)

    # Step 2: find duplicates
    duplicates_deposition = {val: names_list for val, names_list in locations.items() if len(names_list) > 1}

    print("Duplicates found:")
    for val, names_list in duplicates_deposition.items():
        print(f"'{val}' appears in arrays: {names_list}")

    # Step 3: decide best owner (highest rate)
    best_owner = {}
    move_log = []

    name_to_rate = dict(zip(arrays_to_check_duplicates_deposition_names, deposition_velocities))  # helper mapping

    for val, names_list in locations.items():
        best = max(names_list, key=lambda n: name_to_rate[n])
        best_owner[val] = best

        # log removals
        for n in names_list:
            if n != best:
                move_log.append({
                    "value": val,
                    "removed_from": n,
                    "kept_in": best
                })

    # Step 4: build new arrays with "_dep" naming
    arrays_dep = {}

    for name, arr in zip(arrays_to_check_duplicates_deposition_names, arrays_to_check_duplicates_deposition):
        filtered = [val for val in arr if best_owner[val] == name]
        arrays_dep[f"{name}_dep"] = np.array(filtered, dtype=arr.dtype)

    # Step 5: print updated arrays
    print("\nUpdated arrays:")
    for name, arr in arrays_dep.items():
        print(f"{name}: {arr}")

    # Optional: print move log
    print("\nMove log:")
    for entry in move_log:
        print(f"{entry['value']} removed from {entry['removed_from']} -> kept in {entry['kept_in']}")

    #Log what was changed from the duplicates in CSV file
    #with open(dirpath + 'Mechanism_info/GEOSChem_duplicate_log.csv', 'w', newline='') as f:
    with open('CRACMM_duplicate_log.csv', 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['value', 'removed_from', 'kept_in'])
        writer.writeheader()
        writer.writerows(move_log)

    #################################
    # Now let's organize other groupings not related to deposition
    # We need to classify nitrocresols, nitrophenols, nitrotoluenes, and nitrobenzenes as NOx reservoir species

    #Nitrocresol: Has nitro, has only one phenol, has methyl group bonded to aromatic ring (cresol)
    nitro_cresols =  df_updated_info.loc[(df_updated_info['Nitros'] > 0) & (df_updated_info['Phenols'] == 1) & (df_updated_info['Is_Radical'] == 0) & (df_updated_info['Cresol'] > 0)]
    #Nitrophenol: Has nitro, has a phenol, does not have a methyl group
    nitro_phenols =  df_updated_info.loc[(df_updated_info['Nitros'] > 0) & (df_updated_info['Phenols'] == 1) & (df_updated_info['Is_Radical'] == 0) & (df_updated_info['Cresol'] == 0)]
    #Nitrotoluene: Has nitro, has toluene
    nitro_toluene =  df_updated_info.loc[(df_updated_info['Nitros'] > 0) & (df_updated_info['Toluene'] > 0) & (df_updated_info['Is_Radical'] == 0) & (df_updated_info['Phenols'] == 0)]
    #Nitrobenzene: Has nitro, has benzene ring, is not a phenol, does not have toluene
    nitro_benzene =  df_updated_info.loc[(df_updated_info['Nitros'] > 0) & (df_updated_info['Benzene'] > 0) &  (df_updated_info['Is_Radical'] == 0) & (df_updated_info['Phenols'] == 0) & (df_updated_info['Toluene'] == 0)]
    #nitrocatechol
    nitro_catechol =  df_updated_info.loc[(df_updated_info['Nitros'] > 0) & (df_updated_info['Phenols'] == 2) & (df_updated_info['Is_Radical'] == 0)]

    nitro_cresols_arr = nitro_cresols['Species'].values
    nitro_phenols_arr = nitro_phenols['Species'].values
    nitro_toluene_arr = nitro_toluene['Species'].values
    nitro_benzene_arr = nitro_benzene['Species'].values
    nitro_catechol_arr = nitro_catechol['Species'].values

    arrays_to_check_duplicates_nitros = [nitro_cresols_arr, nitro_phenols_arr, nitro_toluene_arr, nitro_benzene_arr, nitro_catechol_arr]
    arrays_to_check_duplicates_nitros_names = ['nitro_cresols_arr', 'nitro_phenols_arr', 'nitro_toluene_arr', 'nitro_benzene_arr', 'nitro_catechol_arr']
    # Step 1: map value -> list of array names
    locations = defaultdict(list)

    for name, arr in zip(arrays_to_check_duplicates_nitros_names, arrays_to_check_duplicates_nitros):
        for val in arr:
            locations[val].append(name)

    # Step 2: find duplicates
    duplicates_nitros = {val: names_list for val, names_list in locations.items() if len(names_list) > 1}

    print("Duplicates found for nitro duplicates:")
    for val, names_list in duplicates_nitros.items():
        print(f"'{val}' appears in arrays: {names_list}")

    #No duplicates across the different nitro groupings

    #############
    # Now we check for duplicates in other groupings.
    pans_species_names = df_updated_info.loc[(df_updated_info['PANs'] > 0) & (df_updated_info['Is_Radical'] == 0)]
    non_pan_peroxy_nitrates = df_updated_info.loc[(df_updated_info['RO2NO2s'] > 0) & (df_updated_info['Is_Radical'] == 0)]
    total_sum_nitrates = pd.concat([alkyl_nitrates_excluding_peroxy_nitrates_and_pans, pans_species_names, non_pan_peroxy_nitrates], ignore_index = True)
    pans_arr = pans_species_names['Species'].values
    ro2no2_arr = non_pan_peroxy_nitrates['Species'].values
    total_sum_nitrates_arr = total_sum_nitrates['Species'].values

    # #Check if there is any overlap between 
    # # RONO2s: Alkyl Nitrates - Organonitrates, catches R-NO3, including ions, excluding peroxy nitrates and excluding PANs
    # # RO2NO2s: Non-PAN Peroxy Nitrates 
    # # PANs: Peroxy-acyl Nitrates 
    overlap_colcheck = ['RONO2s', 'RO2NO2s', 'PANs']
    mask = (df_updated_info[overlap_colcheck] >= 1).sum(axis=1) >= 2
    overlap_results =  df_updated_info.loc[mask]
    overlap_results_print = df_updated_info.loc[mask, ['Species'] + overlap_colcheck]
    print('Overlap RONO2s, RO2NO2s, PANs: \n', overlap_results_print)
    #We find several overlaps between RONO2s and PANs. These species are RONO2s but are also PANs. 

    # Creates dictionary by passing Series objects as values
    den_groupings = {
        # Save the RONO2s that exclude PANs
        'RONO2s': rono2_arr,
        'PANs': pans_arr,
        'RO2NO2s': ro2no2_arr,
        'Total_sum_nitrates': total_sum_nitrates_arr,
        'Carboxylic_acids': carboxylic_acids_arr,
        'HydroPeroxides': hydroperoxides_arr,
        'Nitrocresols': nitro_cresols_arr,
        'Nitrophenols': nitro_phenols_arr,
        'Nitrotoluene': nitro_toluene_arr,
        'Nitrobenzene': nitro_benzene_arr,
        'Nitrocatechol': nitro_catechol_arr,
        'RONO2s_dep': arrays_dep['rono2_arr_dep'], 
        'CarboxylicAcids_dep': arrays_dep['carboxylic_acids_arr_dep'], 
        'Hydroperoxides_dep': arrays_dep['hydroperoxides_arr_dep'],
        'OVOCs': arrays_dep['ovocs_arr_dep']
    }

    # Save the dictionary in an output .mat file: 
    savemat("CRACMM_species_classifications.mat", {"CRACMM_species_classifications": den_groupings})
    
    return den_groupings

def deposition_file_CRACMM(den_groupings):
    #Copy and paste the text file made with this function to use in the deposition_MCM.m file
    #make a text file

    #Organic Hydroperoxides: ROOH
    # with open(dirpath + 'Mechanism_info/deposition_GEOSChem_text_for_matlab_hydroperoxides.txt', 'w', encoding="utf-8") as f:
    with open('deposition_CRACMM_text_for_matlab_hydroperoxides.txt', 'w', encoding="utf-8") as f:
        f.write('%----------------------------------------------\n')
        f.write('% Organic Hydroperoxides\n')
        f.write('%----------------------------------------------\n')
        for rooh_spec in den_groupings['Hydroperoxides_dep']:
            f.write('i=i+1;\n')
            rooh_rnames_string = "Rnames{i} = '" + str(rooh_spec) + " = DEP';\n"
            f.write(rooh_rnames_string)

            f.write('k(:,i) = ROOH_dep./(BLH_cm); %  s-1\n')

            rooh_gstr_string = "Gstr{i,1} = '" + str(rooh_spec) + "';\n"
            f.write(rooh_gstr_string)

            rooh_fspecies_string = "f" + str(rooh_spec) + "(i) = -1;\n\n"
            f.write(rooh_fspecies_string)
        print('Saved text file to filepath: ', 'deposition_CRACMM_text_for_matlab_hydroperoxides.txt')

    #RONO2s
    # with open(dirpath + 'Mechanism_info/deposition_GEOSChem_text_for_matlab_rono2s.txt', 'w', encoding="utf-8") as f:
    with open('deposition_CRACMM_text_for_matlab_rono2s.txt', 'w', encoding="utf-8") as f:
        f.write('%----------------------------------------------\n')
        f.write('% RONO2s (Organic nitrates)\n')
        f.write('%----------------------------------------------\n')

        for rono2_spec in den_groupings['RONO2s_dep']:
            f.write('i=i+1;\n')
            rono2_rnames_string = "Rnames{i} = '" + str(rono2_spec) + " = DEP';\n"
            f.write(rono2_rnames_string)

            f.write('k(:,i) = RONO2_dep./(BLH_cm); %  s-1\n')

            rono2_gstr_string = "Gstr{i,1} = '" + str(rono2_spec) + "';\n"
            f.write(rono2_gstr_string)

            rono2_fspecies_string = "f" + str(rono2_spec) + "(i) = -1;\n\n"
            f.write(rono2_fspecies_string)
        print('Saved text file to filepath: ', 'deposition_CRACMM_text_for_matlab_rono2s.txt')

    #CarboxylicAcids
    # with open(dirpath + 'Mechanism_info/deposition_GEOSChem_text_for_matlab_CarboxylicAcids.txt', 'w', encoding="utf-8") as f:
    with open('deposition_CRACMM_text_for_matlab_CarboxylicAcids.txt', 'w', encoding="utf-8") as f:
        f.write('%----------------------------------------------\n')
        f.write('% Carboxylic Acids\n')
        f.write('%----------------------------------------------\n')

        for carboxylicacids_spec in den_groupings['CarboxylicAcids_dep']:
            f.write('i=i+1;\n')
            carboxylicacids_rnames_string = "Rnames{i} = '" + str(carboxylicacids_spec) + " = DEP';\n"
            f.write(carboxylicacids_rnames_string)

            f.write('k(:,i) = RCOOH_dep./(BLH_cm); %  s-1\n')

            carboxylicacids_gstr_string = "Gstr{i,1} = '" + str(carboxylicacids_spec) + "';\n"
            f.write(carboxylicacids_gstr_string)

            carboxylicacids_fspecies_string = "f" + str(carboxylicacids_spec) + "(i) = -1;\n\n"
            f.write(carboxylicacids_fspecies_string)
        print('Saved text file to filepath: ', 'deposition_CRACMM_text_for_matlab_CarboxylicAcids.txt')

    #OVOCs
    # with open(dirpath + 'Mechanism_info/deposition_GEOSChem_text_for_matlab_ovocs.txt', 'w', encoding="utf-8") as f:
    with open('deposition_CRACMM_text_for_matlab_ovocs.txt', 'w', encoding="utf-8") as f:
        f.write('%----------------------------------------------\n')
        f.write('% OVOCs (Oxidized VOCs)\n')
        f.write('%----------------------------------------------\n')
        for ovocs_spec in den_groupings['OVOCs']:
            f.write('i=i+1;\n')
            ovocs_rnames_string = "Rnames{i} = '" + str(ovocs_spec) + " = DEP';\n"
            f.write(ovocs_rnames_string)

            f.write('k(:,i) = OVOC_dep./(BLH_cm); %  s-1\n')

            ovocs_gstr_string = "Gstr{i,1} = '" + str(ovocs_spec) + "';\n"
            f.write(ovocs_gstr_string)

            ovocs_fspecies_string = "f" + str(ovocs_spec) + "(i) = -1;\n\n"
            f.write(ovocs_fspecies_string)
        print('Saved text file to filepath: ', 'deposition_CRACMM_text_for_matlab_ovocs.txt')

### CALL FUNCTIONS ###

# path_all_species = '../../CRACMM/CRACMM_20251014update/metadata/cracmm3m/cracmm3m_metadata.csv'
# df_allspecies = pd.read_csv(path_all_species)

# group_all_species(
#     df_in = df_allspecies, 
#     smart_groups = smart_groups)

den_groupings = choose_smarts_classifications()
deposition_file_CRACMM(den_groupings)
