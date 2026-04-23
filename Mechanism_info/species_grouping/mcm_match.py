import pandas as pd # ok with rdkit 
import pickle
import numpy as np
import ast
import sys
from rdkit import Chem

from rdkit.Chem import  Descriptors, rdMolDescriptors, Fragments
from scipy.io import savemat
from collections import OrderedDict

dirpath = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/'

#Read in excel file of MCM species
path_to_mcm_file = dirpath + 'Mechanism_info/mcm_species.xlsx'
df_mcm_species = pd.read_excel(path_to_mcm_file)
#only include Mech == MCM
df_mcm_species  = df_mcm_species.loc[(df_mcm_species['Mech'] == 'MCM')]
#print(df_mcm_species)

path_to_mcm_using_file = dirpath + 'Mechanism_info/mcm_bb_sherwen_species.xlsx'
df_mcm_species_using = pd.read_excel(path_to_mcm_using_file, index_col=0)
df_mcm_species_using = df_mcm_species_using.drop_duplicates(subset="Name")
df_mcm_species_using.to_excel(dirpath + 'Mechanism_info/mcm_bb_sherwen_species.xlsx')

#Species in our mechanism not in Jessica's MCM Species spreadsheet
not_in_df_mcm_species = df_mcm_species_using[~df_mcm_species_using["Name"].isin(df_mcm_species["Name"])]
#print(not_in_df_mcm_species)

#Species in Jessica's MCM Species spreadsheet, not in our mechanism
not_in_df_mcm_species_using = df_mcm_species[~df_mcm_species["Name"].isin(df_mcm_species_using["Name"])]
#print(not_in_df_mcm_species_using)

#Took species missing and put them into spreadsheet:
path_to_mcm_missing_species = dirpath + 'Mechanism_info/mcm_species_missing.xlsx'
df_mcm_missing_species = pd.read_excel(path_to_mcm_missing_species, index_col=0)
# print(df_mcm_missing_species)

df_mcm_missing_species_copy = df_mcm_missing_species.copy()

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
    'Enols':	'[$([#6](=[#6])(-[#8X2H1]));!$([$([$(c1ccccc1)](-[#8X2H1]));!$([c](=[#8])(-[#8X2H1]))])]',
    #'Phenols': '[c][OX2H]',
    'Phenols': '[$([$(c1ccccc1)](-[#8X2H1]));!$([c](=[#8])(-[#8X2H1]))]', # -OH on benzene ring (but no carboxylic acids)
    'Thiols':	'[#6]-[#16X2H1]',
    'Dihydroxys':	'[#6X4;!$([#6]=[#8])](-[#8X2H1])(-[#8X2H1])',
    'Aliphatic_Alcohols':	'[$([CX4R0](-[#8X2H1]));!$([#6X3;!$([#6X3H])](=[#8])[#8X2H])!$([#6X4;!$([#6]=[#8])](-[#8X2H1])(-[#8X2H1]))]',
    'Other_Alcohols':	'[$([#6;!$([$([#6](=[#6])(-[#8X2H1]));!$([$([$(c1ccccc1)](-[#8X2H1]));!$([c](=[#8])(-[#8X2H1]))])]);!$([#6]-[#16X2H1]);!$([#6X4;!$([#6]=[#8])](-[#8X2H1])(-[#8X2H1]))][#8X2H1]);!$([$([$(c1ccccc1)](-[#8X2H1]));!$([c](=[#8])(-[#8X2H1]))]);!$([#6X3;!$([#6X3H])](=[#8])[#8X2H]);!$([$([CX4R0](-[#8X2H1]));!$([#6X3;!$([#6X3H])](=[#8])[#8X2H])!$([#6X4;!$([#6]=[#8])](-[#8X2H1])(-[#8X2H1]))])]',
    'C':	'[#6]',
    'H':	'[#1]',
    'O':	'[#8]',
    'N':	'[#7]',
    'S':	'[#16]',
    'Cl':	'[#17]',
    'Br':	'[#35]',
    'Toluene': '[$(Cc1ccccc1)]',
    'BenzeneRing': 'c1ccccc1',
    'Methyl': '[#6](-[#1])(-[#1])(-[#1])'
})


def get_rdkit_frags(df_in, use, molec): 
    df_out=pd.DataFrame()
    df_out[use]= df_in[use].copy()
    
    fr_cols=['fr_Hydroxyls','fr_Carbonyls','fr_Carboxylic_Acids','fr_Aldehydes',
    'fr_Esters','fr_Ketones','fr_Epoxides','fr_Nitros','fr_Ethers','fr_Phenols',
    'fr_Thiols']
    
    for col in fr_cols: df_out[col]=list(np.full([len(df_out),1], np.nan)) 
    
    for  i, sp in enumerate(df_in[use]): 
        molec = Chem.inchi.MolFromInchi(df_in.loc[i,'InChI'],sanitize=False, removeHs=False,logLevel=None)  #Turn this molec it into an RDKit molecule object. 
        
        if molec is not None:
            molec.UpdatePropertyCache(strict=True)  # for radicals!
            
        df_out.at[i,'fr_Hydroxyls'] = Fragments.fr_Al_OH(molec) # Number of aliphatic hydroxyl groups
        df_out.at[i,'fr_Carbonyls'] = Fragments.fr_C_O(molec) # Number of Carbonyl O 
        df_out.at[i,'fr_Carboxylic_Acids'] = Fragments.fr_COO2(molec) #N umber of carboxylic acids 
        df_out.at[i,'fr_Aldehydes']=Fragments.fr_aldehyde(molec) # Number of aldehydes
        df_out.at[i,'fr_Esters']=Fragments.fr_ester(molec) # Number of esters
        df_out.at[i,'fr_Ketones']=Fragments.fr_ketone(molec) # Number of ketones
        df_out.at[i,'fr_Nitros']= Fragments.fr_nitro(molec) # Number of nitro groups
        df_out.at[i,'fr_Ethers']=Fragments.fr_ether(molec) # Number of ether oxygens (including phenoxy)
        df_out.at[i, 'fr_Phenols']=Fragments.fr_phenol_noOrthoHbond(molec) # Number of phenolic OH excluding ortho intramolecular Hbond substituents
        df_out.at[i, 'fr_Thiols']=Fragments.fr_SH(molec) # Number of thiol groups
    
    return df_out

#Get the information for all the species in the NOAA BB and Halogens Sherwen mechanism that aren't in the MCM Species List from Jessica
def get_info_on_missing_species():
    for i, species in enumerate(df_mcm_missing_species_copy['Name']):
        #print('Checking for: ', species)
        if species == 'GUAIACOLO':
            print('Addressing GUAIACOLO species')
            # Your radical SMILES
            rad_smiles_guaiacolo = 'COC1=C(O)C=C[C.]C1O'

            # Replace [C.] with [C] temporarily so RDKit can parse
            neutral_smiles_guaiacolo = rad_smiles_guaiacolo.replace('[C.]', '[C]')

            # Create a Mol object from the neutral SMILES
            guaiacolo_mol = Chem.MolFromSmiles(neutral_smiles_guaiacolo)
            if guaiacolo_mol is None:
                raise ValueError("Failed to parse SMILES. Check your input.")

            # Identify the radical atom(s)
            # In this simple example, we know it's the carbon previously written as [C.]
            # Let's find its index by matching atoms: C with degree 2 and not part of OH/OMe
            for atom in guaiacolo_mol.GetAtoms():
                if atom.GetSymbol() == 'C':
                    # crude check: match the atom that corresponds to radical position
                    # here we just take the 5th carbon in the ring (0-indexed)
                    if atom.GetIdx() == 5:  # adjust index if needed
                        atom.SetNumRadicalElectrons(1)

            inchi = Chem.MolToInchi(guaiacolo_mol)
            form=rdMolDescriptors.CalcMolFormula(guaiacolo_mol)
            mw=Descriptors.MolWt(guaiacolo_mol)
            canon_smiles =Chem.MolToSmiles(guaiacolo_mol, canonical = True)
            df_mcm_missing_species_copy.loc[i, 'InChI'] = inchi
            df_mcm_missing_species_copy.loc[i, 'Formula'] = form
            df_mcm_missing_species_copy.loc[i, 'Molecular_Weight'] = mw
            df_mcm_missing_species_copy.loc[i, 'Canonical_SMILES'] = canon_smiles

        #Enter loop if 'InChI' is not empty and not NaN
        elif (pd.notna(df_mcm_missing_species_copy['InChI'][i]) and df_mcm_missing_species_copy['InChI'][i] != ''):
            #print('first loop entered for species ', species)
            #strip any spaces
            df_mcm_missing_species_copy.loc[i,'InChI'] = df_mcm_missing_species_copy['InChI'][i].strip(' ')
            # Use each string to make an RDKit molecule object. 
            molec = Chem.MolFromInchi(df_mcm_missing_species_copy['InChI'][i])

            if molec is not None:
                molec.UpdatePropertyCache(strict=True)  # for radicals!
                # Pull formula, molecular weight, its canonical smiles, and InChI strings. 
                form=rdMolDescriptors.CalcMolFormula(molec)
                mw=Descriptors.MolWt(molec)
                smiles =Chem.MolToSmiles(molec)
                df_mcm_missing_species_copy.loc[i, 'Formula'] = form
                df_mcm_missing_species_copy.loc[i, 'Molecular_Weight'] = mw
                df_mcm_missing_species_copy.loc[i, 'Canonical_SMILES'] = smiles

        #Enter loop if 'SMILES' is not empty and not NaN, and InChi is empty or NaN
        elif (pd.notna(df_mcm_missing_species_copy['SMILES'][i]) and df_mcm_missing_species_copy['SMILES'][i] != '' and (df_mcm_missing_species_copy['InChI'][i] == '' or pd.isna(df_mcm_missing_species_copy['InChI'][i]))):
            molec = Chem.MolFromSmiles(df_mcm_missing_species_copy['SMILES'][i])
            #print('second loop entered for species ', species)

            if molec is not None:
                molec.UpdatePropertyCache(strict=True)  # for radicals!
                # Pull formula, molecular weight, its canonical smiles, and InChI strings. 
                inchi = Chem.MolToInchi(molec)
                form=rdMolDescriptors.CalcMolFormula(molec)
                mw=Descriptors.MolWt(molec)
                smiles =Chem.MolToSmiles(molec)
                df_mcm_missing_species_copy.loc[i, 'InChI'] = inchi
                df_mcm_missing_species_copy.loc[i, 'Formula'] = form
                df_mcm_missing_species_copy.loc[i, 'Molecular_Weight'] = mw
                df_mcm_missing_species_copy.loc[i, 'Canonical_SMILES'] = smiles

    df_mcm_missing_species_copy.to_excel(dirpath + 'Mechanism_info/mcm_species_missing_updated.xlsx')
def subtract_mcm_species_minus_extra():
    path_to_mcm_missing_species_updated = dirpath + 'Mechanism_info/mcm_species_missing_updated.xlsx'
    df_mcm_missing_species_updated = pd.read_excel(path_to_mcm_missing_species_updated, index_col=0)

    df_diff = df_mcm_species.merge(not_in_df_mcm_species_using, how="left", indicator=True)
    df_diff = df_diff[df_diff["_merge"] == "left_only"].drop(columns="_merge")
    df_diff.drop(columns=df_diff.columns[0], axis=1, inplace=True)
    
    df_mcm_bb_halogens_sherwen_total = pd.concat([df_diff, df_mcm_missing_species_updated], ignore_index=True, sort=False)
    #dupes = df_mcm_bb_halogens_sherwen_total[df_mcm_bb_halogens_sherwen_total["Name"].duplicated(keep=False)]
    df_mcm_bb_halogens_sherwen_total = df_mcm_bb_halogens_sherwen_total.drop_duplicates(subset="Name", ignore_index=True)
    df_mcm_bb_halogens_sherwen_total.to_excel(dirpath + 'Mechanism_info/mcm_species_bb_sherwen_total_info.xlsx')
def group_all_species(df_in, use, smart_groups):
    df_out=pd.DataFrame()
    df_out[use]= df_in[use].copy()
    df_out['Epoxides']=list(np.full([len(df_out),1], np.nan)) 
    print('df_in: ', df_in)
    # df_out['Toluene']=list(np.full([len(df_out),1], np.nan)) 
    # df_out['BenzeneRing']=list(np.full([len(df_out),1], np.nan)) 
    # df_out['Methyl']=list(np.full([len(df_out),1], np.nan)) 
    
    for i, species in enumerate(df_allspecies['Name']):
        #print(df_allspecies.loc[i,'Name']) # Print off name of what you're parsing. 
        df_allspecies.loc[i,'InChI'] = str(df_allspecies.loc[i,'InChI']).strip(' ')
        molec = Chem.MolFromInchi(df_allspecies.loc[i,'InChI'],sanitize=False, removeHs=False,logLevel=None)  #Turn this molec it into an RDKit molecule object. 
        
        if molec is not None:
            molec.UpdatePropertyCache(strict=True)  # for radicals!
            
            for key in smart_groups:  # Loop over every functional group you want to search for. 
                print('Testing ', key, '\n')
                # Turn the SMARTs string for this functional group into a RDKit molec fragment. 
                frag = Chem.MolFromSmarts(smart_groups[key])
                
                # Get a list of the indices of atom #s in molecule that match this fragment 
                inds=list(molec.GetSubstructMatches(frag))
                
                # Save the len of this list as the # of functional group matches you found!)
                df_out.at[i,key]=np.int64(len(inds))
                
                # # Use RDKit to get fragments(not always as specific as our group matches...) 
                #rd_frags=get_rdkit_frags(df_in, use, molec)
            
            df_out.at[i,'Epoxides'] = Fragments.fr_epoxide(molec) # Number of epoxide rings 
            # df_out.at[i, 'Toluene'] = Fragments.fr_toluene(molec)
            # df_out.at[i, 'BenzeneRing'] = Fragments.fr_benzenering(molec)
            # df_out.at[i, 'Methyl'] = Fragments.fr_methyl(molec)
            
    # Add a column that has the number of OH groups on a compounds that might cause it to be an organic acid... 
    Organic_Acid_OHs=['Enols','Phenols','Thiols', 'Carboxylic_Acids']
    df_out['Organic_Acid_OHs']=df_out[Organic_Acid_OHs].sum(axis=1)

    print('df_out: ', df_out)
    print('df_allspecies: ', df_allspecies)

    #If we merge, there will be a double counting of some of the parameters like "Is_Radical". Instead, we replace the older
    #duplicate columns with the newer calculations, just in case the older code is wrong. We'll have all the older and newer
    #columns in the merged dataframe, but only one column of the duplicates, taking in the new values.

    duplicate_columns = ['Is_Radical', 'RO2s', 'ROs', 'Acyl_RO2s', 'Acyl_ROs', 
                       'All_NO3s', 'Non_PAN_NO3s', 'Tertiary_NO3s', 'RONO2s', 
                       'RO2NO2s', 'PANs', 'Nitros', 'Peroxides', 'OrganicPeroxides',
                       'HydroPeroxides', 'Peracids', 'Carbonyls', 'Carboxylic_Acids', 
                       'Ketones', 'Aldehydes', 'Esters', 'Ethers', 'Carbonates', 'All_OHs', 
                       'Hydroxyls', 'Enols', 'Phenols', 'Thiols', 'Dihydroxys', 'Aliphatic_Alcohols',
                         'Other_Alcohols', 'C', 'H', 'O', 'N', 'S', 'Cl', 'Br', 'Epoxides', 'Organic_Acid_OHs']

    df_merged_total = (
        df_allspecies
        .merge(
            df_out[['Name'] + duplicate_columns  + ['Toluene'] + ['BenzeneRing'] + ['Methyl']],
            on='Name',
            how='left',
            suffixes=('', '_right')
        )
    )

    for col in duplicate_columns:
        df_merged_total[col] = df_merged_total[f'{col}_right']
        df_merged_total.drop(columns=f'{col}_right', inplace=True)

    print('df_merged_total: ', df_merged_total)
    print(list(df_merged_total.columns.values))
    df_merged_total.to_excel(dirpath + 'Mechanism_info/mcm_allspecies_bb_sherwen_info_updated.xlsx')

def choose_smarts_classifications():
    df_updated_info = pd.read_excel(dirpath + 'Mechanism_info/mcm_allspecies_bb_sherwen_info_updated.xlsx', index_col=0)

    # #Check if there is any overlap between 
    # # RONO2s: Alkyl Nitrates - Organonitrates, catches R-NO3, including ions, excluding peroxy nitrates and excluding PANs
    # # RO2NO2s: Non-PAN Peroxy Nitrates 
    # # PANs: Peroxy-acyl Nitrates 
    # overlap_colcheck = ['RONO2s', 'RO2NO2s', 'PANs']
    # mask = (df_updated_info[overlap_colcheck] >= 1).sum(axis=1) >= 2
    # overlap_results =  df_updated_info.loc[mask]
    # overlap_results_print = df_updated_info.loc[mask, ['Name'] + overlap_colcheck]
    # print(overlap_results_print)

    #We find several overlaps between RONO2s and PANs. These species are RONO2s but are also PANs. 
    #We want to categorize them as PANs when we separate into subcategories. 
    alkyl_nitrates_excluding_peroxy_nitrates_and_pans = df_updated_info.loc[(df_updated_info['RONO2s'] > 0) & (df_updated_info['PANs'] == 0) & (df_updated_info['Is_Radical'] == 0)]
    pans_species_names = df_updated_info.loc[(df_updated_info['PANs'] > 0) & (df_updated_info['Is_Radical'] == 0)]
    non_pan_peroxy_nitrates = df_updated_info.loc[(df_updated_info['RO2NO2s'] > 0) & (df_updated_info['Is_Radical'] == 0)]
    total_sum_nitrates = pd.concat([alkyl_nitrates_excluding_peroxy_nitrates_and_pans, pans_species_names, non_pan_peroxy_nitrates], ignore_index = True)
    carboxylic_acids = df_updated_info.loc[(df_updated_info['Carboxylic_Acids'] > 0)]
    hydroperoxides = df_updated_info.loc[(df_updated_info['HydroPeroxides'] > 0)]

    nitro_cresols =  df_updated_info.loc[(df_updated_info['Nitros'] > 0) & (df_updated_info['Phenols'] > 0) & (df_updated_info['Is_Radical'] == 0) & (df_updated_info['Methyl'] > 0)]
    nitro_phenols =  df_updated_info.loc[(df_updated_info['Nitros'] > 0) & (df_updated_info['Phenols'] > 0) & (df_updated_info['Is_Radical'] == 0) & (df_updated_info['Methyl'] == 0)]
    nitro_toluene =  df_updated_info.loc[(df_updated_info['Nitros'] > 0) & (df_updated_info['Toluene'] > 0) &  (df_updated_info['Is_Radical'] == 0)]
    nitro_benzene =  df_updated_info.loc[(df_updated_info['Nitros'] > 0) & (df_updated_info['BenzeneRing'] > 0) &  (df_updated_info['Is_Radical'] == 0) & (df_updated_info['Phenols'] == 0)] #May also have to set condition of no toluene
    print('Nitro Cresols:', nitro_cresols)
    print('Nitro Phenols:', nitro_phenols)
    print('Nitro Toluene:', nitro_toluene)
    print('Nitro Benzene:', nitro_benzene)

    #Alfie counts any carbonyls and any alcohols as OVOCs. 
    #Our carbonyls functional SMARTs includes all carbonyls including carboxylic acids. Since we have a separate depositional velocity, we want all the carbonyls excluding carboxylic acids.
    carbonyls_group = (df_updated_info['Carbonyls'] > 0) & (df_updated_info['Carboxylic_Acids'] == 0)
    dihydroxys_group = (df_updated_info['Dihydroxys'] > 0)
    aliphatic_alcohols_group = (df_updated_info['Aliphatic_Alcohols'] > 0)
    other_alcohols_group = (df_updated_info['Other_Alcohols'] > 0)
    # ovocs = df_updated_info.loc[(carbonyls_group) | (dihydroxys_group) | (aliphatic_alcohols_group) | (other_alcohols_group)]
    ovocs = df_updated_info.loc[((df_updated_info['Carbonyls'] > 0) & (df_updated_info['Carboxylic_Acids'] == 0)) | ((df_updated_info['Dihydroxys'] > 0)) | ((df_updated_info['Aliphatic_Alcohols'] > 0)) | ((df_updated_info['Other_Alcohols'] > 0))]

    # Check for duplicates from carbonyls, dihidroxys, aliphatic alcohols, and other alcohols setup
    print('OVOCs: ', ovocs)
    # duplicates_carbonyls_alcohols = [i for i in set(ovocs) if ovocs.count(i) > 1]
    # print('duplicates_carbonyls_alcohols: ', duplicates_carbonyls_alcohols)

    # # Oxidized VOCs are all the species that are not primary, not inorganic, not radicals, not RO2s, and not duplicates
    # # STEP 1: Determine Primary VOCs (so that we can subtract them from the list of total species)
    # # Used the primary VOCs exporting from the MCM website to primary_vocs_and_precursors_mcm.xlsx file, along with the precursors in the function precursor_ro2_classification() below (with alphabetized sorting)
    # # Added the furan relevant primary VOCs to the MCM website's VOC precursors to make primary_vocs_mcm_list.csv and primary_vocs_mcm_list.txt
    # # where the CSV file is a column with heading MCM_primary_with_furans and the text file is a list called primary_vocs_mcm_list

    # with open(dirpath + 'Mechanism_info/species_grouping/ovocs/primary_vocs_mcm_list.txt') as f: #read primary vocs
    #     lines = f.readlines()

    # # This primary VOCs list unfortunately includes a few species that we want in our OVOCs list. The OVOCs from Rickly et al., 2023 include:
    # # Acetaldehyde, Acrolein, Formaldehyde, Ethanol, Formic Acid, Butanedione [not identified as primary VOC by MCM website], Glycolaldehyde [not identified as primary VOC by MCM website], 
    # # Isopropanol, MACR, Methyl acetate, Methyl ethyl ketone, Methanol, MVK

    # rickly_ovocs_to_include = ['CH3CHO', 'ACR', 'HCHO', 'C2H5OH', 'HCOOH', 'IPROPOL', 'MACR', 'METHACET', 'MEK', 'CH3OH', 'MVK']

    # # Exclude all PANs and peroxy nitrates
    # pans_species_names
    # non_pan_peroxy_nitrates
    # alkyl_nitrates_excluding_peroxy_nitrates_and_pans

    # # Exclude duplicates, as in anything in the other categories of 
    # hydroperoxides
    # carboxylic_acids    


    #Creates dictionary by passing Series objects as values
    den_groupings = {
        'RONO2s': alkyl_nitrates_excluding_peroxy_nitrates_and_pans['Name'].values,
        'PANs': pans_species_names['Name'].values,
        'RO2NO2s': non_pan_peroxy_nitrates['Name'].values,
        'Total_sum_nitrates': total_sum_nitrates['Name'].values,
        'Carboxylic_acids': carboxylic_acids['Name'].values,
        'HydroPeroxides': hydroperoxides['Name'].values
    }
    print(den_groupings)

    # Save the dictionary in an output .mat file: 
    #savemat(dirpath + "Mechanism_info/species_grouping/MCM_species_classifications.mat", {"MCM_species_classifications": den_groupings})
def precursor_ro2_classification():
    df_updated_info = pd.read_excel(dirpath + 'Mechanism_info/mcm_allspecies_bb_sherwen_info_updated.xlsx', index_col=0)
    #Now we want to create groupings for the precursor species for RO2s. 
    ro2_species = df_updated_info.loc[(df_updated_info['PANs'] > 0)]
    df_updated_info['Precursors'] = df_updated_info['Precursors'].apply(lambda x: ast.literal_eval(x) if pd.notna(x) else x)
    print(type(df_updated_info['Precursors'][0]))

    #Go through all the precursors and print out all the possible species
    all_precursors = [spec for sublist in df_updated_info['Precursors'] if isinstance(sublist, list) for spec in sublist]
    # Get unique values
    unique_precursors = list(set(all_precursors))
    # print(unique_precursors)

    #region: output precursors
    precursors_listed = ['M2HEX', 'CHEX2ENE', 'IC4H10', 'OXYL', 'DIETETHER', 'PETHTOL', 
    'NPROACET', 'DMM', 'CH3BR', 'PBENZ', 'NC8H18', 'TOLUENE', 'NC5H12', 
    'C5H8_RCIM', 'ETHOX', 'MEFURFURAL', 'CCL2CH2', 'M3HEX', 'CBUT2ENE', 
    'PXYL', 'ME3BUOL', 'DIBRET', 'MTBE', 'NC9H20', 'IC5H12', 'CHCL3', 
    'CH2CL2', 'HEX1ENE', 'FURFURAL', 'BOX2PROL', 'MEPROPENE', 'TRICLETH', 
    'CHCL2CH3', 'MBO', 'DM23BU2ENE', 'CL12PROP', 'NC6H14', 'TM135B', 'DIMEFURAN', 
    'CH2CLCH2CL', 'TDICLETH', 'BCARY', 'STYRENE', 'DMS', 'M3PE', 'TM123B', 'MXYL', 
    'CHCL2CHCL2', 'C4ALDB', 'TPENT2ENE', 'MO2EOL', 'GUAIACOL', 'APINENE', 'CPENT2ENE',
    'DIME35EB', 'TM124B', 'C2H2', 'ME2BUT2ENE', 'FURAN', 'C5H8', 'TCE', 'VINCL',
    'DIET35TOL', 'DIIPRETHER', 'M22C4', 'CH3CL', 'NC7H16', 'SBUTACET', 'SYRINGOL', 
    'PR2OHMOX', 'CH2CLCHCL2', 'OETHTOL', 'NC10H22', 'NC11H24', 'NBUTACET', 'CH3CH2CL',
    'EBENZ', 'M23C4', 'CH3CCL3', 'NEOP', 'CDICLETH', 'CHEX', 'IPBENZ', 'TBUT2ENE', 
    'M2PE', 'THEX2ENE', 'BPINENE', 'LIMONENE', 'NC12H26', 'ME3BUT1ENE', 'M2F', 
    'ETBE', 'METHTOL', 'ME2BUT1ENE', 'BUOX2ETOH']
    print('sorted: ', sorted(precursors_listed))
    #endregion


    #All groupings:
    # biomass_burning_grouping, xylene_grouping, eth_benz_grouping, methbenz_grouping, toluene_grouping, other_alkyl_benz_grouping, styrene_grouping, summed_aromatic_grouping, traffic_aromatic_grouping, 
    # major_traffic_aromatic_grouping, all_fuel_aromatic_grouping, major_fuel_aromatic_grouping, all_industrial_aromatic_grouping, major_industrial_aromatic_grouping, acetate_esters,
    # fuel_evaporative_ethers_grouping, industrial_solvent_voc_ether_grouping, monoterpenes_grouping, isoprene_grouping, greater_than_c4_alkenes_mostly_from_isoprene_grouping, isoprene_products_grouping,
    # biogenics_total_grouping, biogenic_OVOCs_grouping, greater_than_c4_alkenes_from_biogenics_and_anthro_grouping, primarily_fossil_fuel_combustion_products_grouping, primarily_fuel_evaporative_sources_grouping,
    # other_petroleum_related_vocs_grouping, halogenated_vocs_grouping, dms_other_grouping

    #BB includes Furans, guaiacol, syringol
    biomass_burning_grouping = ['FURAN', 'M2F', 'DIMEFURAN', 'FURFURAL', 'MEFURFURAL', 'GUAIACOL', 'SYRINGOL']
    
    xylene_grouping = ['OXYL','MXYL','PXYL']
    eth_benz_grouping = ['EBENZ']
    methbenz_grouping = ['TM135B','TM123B','TM124B', 'DIME35EB', 'DIET35TOL']
    toluene_grouping = ['TOLUENE', 'METHTOL','PETHTOL', 'OETHTOL']
    other_alkyl_benz_grouping = ['IPBENZ','PBENZ'] #cumene, propylbenzene
    styrene_grouping = ['STYRENE']
    summed_aromatic_grouping = ['OXYL','MXYL','PXYL', 'EBENZ', 'TM135B','TM123B','TM124B','DIME35EB','DIET35TOL',
                                'TOLUENE', 'METHTOL','PETHTOL','OETHTOL', 'IPBENZ', 'PBENZ', 'STYRENE']
    traffic_aromatic_grouping = ['OXYL','MXYL','PXYL', 'EBENZ', 'TOLUENE', 'METHTOL','PETHTOL', 'OETHTOL', 'IPBENZ','PBENZ', 'TM135B','TM123B','TM124B']
    major_traffic_aromatic_grouping = ['OXYL','MXYL','PXYL', 'EBENZ', 'TOLUENE', 'METHTOL','PETHTOL', 'OETHTOL']
    all_fuel_aromatic_grouping = ['OXYL','MXYL','PXYL', 'EBENZ', 'TM135B','TM123B','TM124B','DIME35EB','DIET35TOL', 'TOLUENE', 'METHTOL','PETHTOL','OETHTOL', 'IPBENZ', 'PBENZ', 'STYRENE']
    major_fuel_aromatic_grouping = ['OXYL','MXYL','PXYL','EBENZ','TOLUENE', 'METHTOL','PETHTOL', 'OETHTOL', 'PBENZ', 'TM135B','TM123B','TM124B']
    all_industrial_aromatic_grouping = ['OXYL','MXYL','PXYL', 'EBENZ', 'TM135B','TM123B','TM124B','DIME35EB','DIET35TOL','TOLUENE', 'METHTOL','PETHTOL','OETHTOL', 'IPBENZ', 'PBENZ', 'STYRENE']
    major_industrial_aromatic_grouping = ['IPBENZ','TM135B','TM123B','TM124B', 'DIME35EB', 'DIET35TOL', 'STYRENE']
    #Minor would be all minus major
    
    #solvent-related OVOCs; typically from solvents or industrial processes
    acetate_esters = ['NPROACET', 'NBUTACET', 'SBUTACET']

    fuel_evaporative_ethers_grouping = [ 'DIIPRETHER', 'ETBE','MTBE','DIETETHER']
    #DIIPRETHER: Diisopropyl ether
    #ETBE: tert-Butyl ethyl ether
    #MTBE: Tert butyl methyl ether
    #DIETETHER: diethyl ether

    industrial_solvent_voc_ether_grouping = ['BOX2PROL','MO2EOL', 'PR2OHMOX', 'BUOX2ETOH', 'ETHOX', 'DMM']
    #ether from industrial/solvent
    #'BOX2PROL' 1-Butoxy-2-propanol
    #'MO2EOL' 2 -methoxyethanol
    #'PR2OHMOX' 1-Methoxy-2-propanol
    #'BUOX2ETOH' 2 -Butoxyethanol
    # 'ETHOX',  #ETHYLENE Oxide
    # 'DMM' #Dimethoxymethane; 1, 2-dibromoethane

    #Biogenics
    monoterpenes_grouping = ['APINENE', 'BPINENE', 'BCARY', 'LIMONENE'] 
    isoprene_grouping = ['C5H8']
    greater_than_c4_alkenes_mostly_from_isoprene_grouping = ['ME2BUT2ENE', 'ME3BUT1ENE', 'ME2BUT1ENE', 'MEPROPENE']
    isoprene_products_grouping = ['C5H8', 'ME2BUT2ENE', 'ME3BUT1ENE', 'ME2BUT1ENE', 'MEPROPENE']
    biogenics_total_grouping =  ['APINENE', 'BPINENE', 'BCARY', 'LIMONENE', 'C5H8', 'ME2BUT2ENE', 'ME3BUT1ENE', 'ME2BUT1ENE', 'MEPROPENE', 'MBO', 'ME3BUOL']
    biogenic_OVOCs_grouping = ['ME2BUT2ENE', 'ME3BUT1ENE', 'ME2BUT1ENE', 'MEPROPENE', 'MBO', 'ME3BUOL']
    #region
    #BCARY is (4Z)-4,11,11-trimethyl-8-methylidenebicyclo(7.2.0)undec-4-ene
    # 'MBO', #2-Methyl-3-buten-2-ol, both branched alcohol & alkene
    # #Primarily biogenic
    # 'ME3BUOL', #3-Methyl-1-butanol, branched alcohol
    # #sometimes group ME3BUOL and MBO as biogenic OVOC or OVOC, or OVOC-derived alkyl RO₂ precursor
    #formed more from isoprene oxidation, all >C4 alkenes
    # 'ME2BUT2ENE',#2-Methyl-2-butene
    # 'ME3BUT1ENE', # 3-METHYL-1-BUTENE; Isopentene
    # 'ME2BUT1ENE', #2-METHYL-1-BUTENE
    # 'MEPROPENE',  #Isobutene
    #endregion

    #both biogenic & anthro >C4 alkenes
    greater_than_c4_alkenes_from_biogenics_and_anthro_grouping = ['TBUT2ENE', 'CBUT2ENE', 'CPENT2ENE', 'TPENT2ENE', 'HEX1ENE', 'CHEX2ENE', 'THEX2ENE', 'TBUT2ENE']
    #region
    # 'CBUT2ENE',#cis-2-butene
    # 'CPENT2ENE',#cis-2-PENTENE
    # 'TPENT2ENE',#trans-2-pentene 
    # 'HEX1ENE', #1-HEXENE;
    # 'CHEX2ENE', #cis-2-Hexene
    # 'THEX2ENE', # trans-2-Hexene
    # 'TBUT2ENE',#trans-2-butene
    #endregion

    #primarily fossil fuel/combustion
    primarily_fossil_fuel_combustion_products_grouping = ['M2PE', 'M3PE', 'M2HEX', 'M3HEX', 'C2H2', 'NC9H20', 'NC10H22', 'NC11H24', 'NC12H26', 'C4ALDB']
    #region
    # 'M2PE', #methylpentane
    # 'M3PE', #3-METHYLPENTANE
    # 'M2HEX', #2-METHYLHEXANE
    # 'M3HEX',  #3-METHYLHEXANE
    # 'C2H2'#acetylene/ethyne; alyne; combustion
    # 'NC9H20', #nonane
    # 'NC10H22', #decane
    # 'NC11H24', #undecane
    # 'NC12H26', #dodecane
    #'C4ALDB'#Butenal; crotonaldehyde ; combustion
    #endregion

    primarily_fuel_evaporative_sources_grouping = ['NC5H12', 'NC6H14', 'NC7H16', 'NC8H18', 'M23C4', 'M22C4', 'IC4H10', 'IC5H12', 'NEOP']
    #region
    # 'NC5H12', #pentane
    # 'NC6H14', #hexane
    # 'NC7H16',  #heptane
    # 'NC8H18', #octane
    # 'M23C4'   #2,3-DIMETHYLBUTANE; fuel-related/evaporative
    # 'M22C4' #2,2-DIMETHYLBUTANE; fuel-related/evaporative
    # 'IC4H10',  #isobutane
    # 'IC5H12', #isopentane / 2-Methylbutane
    # 'NEOP', neopentane; gasoline
    #endregion
    
    other_petroleum_related_vocs_grouping = ['CHEX', 'DM23BU2ENE']
    # 'CHEX', #cyclohexane; fuel and solvent
    # 'DM23BU2ENE', #2,3-DIMETHYL-2-BUTENE; fuel evaporation and petrochemical

    halogenated_vocs_grouping = ['DIBRET', 'CH3BR', 'CH2CL2', 'CCL2CH2', 'CHCL3', 'CH3CL', 'CHCL2CH3', 'CH2CLCH2CL', 'CH3CH2CL', 'TRICLETH', 'TCE', 'CH2CLCHCL2', 'CL12PROP', 'VINCL', 'TDICLETH', 'CDICLETH']
    #region
    #'DIBRET', 1, 2-dibromoethane
    #'CH3BR',  Bromomethane
    # #dichloromethane
    # 'CH2CL2',
    # #Vinylidene chloride 
    # 'CCL2CH2'
    # #chloroform (trichloromethane)
    # 'CHCL3', 
    # #chloromethane/methyl chloride
    # 'CH3CL',
    # #Methylchloroform
    # 'CH3CCL3'
    # #DICHLOROETHANE
    # 'CHCL2CH3',
    # #1,2-dichloroethane
    # 'CH2CLCH2CL'
    # #1,1,2,2-Tetrachloroethane
    # 'CHCL2CHCL2',
    # #Chloroethane
    # 'CH3CH2CL',
    # #TRICHLOROETHYLENE
    # 'TRICLETH',
    # # Tetrachloroethene / Tetrachloroethylene
    # 'TCE',
    # #1,1,2-TRICHLOROETHANE
    # 'CH2CLCHCL2',
    # #1,2 dichloropropane
    # 'CL12PROP',
    # #VINYL CHLORIDE
    # 'VINCL',
    # #TRANS-1,2-DICHLOROETHYLENE
    #  'TDICLETH',
    # # CIS-1,2-DICHLOROETHYLENE
    # 'CDICLETH',
    #endregion

    dms_other_grouping = ['DMS']

    #Checking to make sure there is no overlap between RO2s and Acyl_RO2s column (there shouldn't be)
    # overlap_colcheck_ro2s = ['RO2s', 'Acyl_RO2s']
    # mask_ro2s = (df_updated_info[overlap_colcheck_ro2s] >= 1).sum(axis=1) >= 2
    # overlap_results_ro2s =  df_updated_info.loc[mask_ro2s]
    # overlap_results_ro2s_print = df_updated_info.loc[mask_ro2s, ['Name'] + overlap_colcheck_ro2s]
    # print(overlap_results_ro2s_print)

    #Get list of RO2s and compare to the one in the official NOAA Biomass Burning Mechanism
    ro2_species = df_updated_info.loc[(df_updated_info['RO2s'] > 0) | (df_updated_info['Acyl_RO2s'] > 0)]
    ro2_species_smarts_list = ro2_species['Name'].tolist()
    
    #region: mcm_bb_ro2_list
    mcm_bb_ro2_list = ['HOCH2CH2O2', 'HO1C3O2', 'HYPROPO2', 'IPROPOLO2', 'NBUTOLAO2', 'NBUTOLBO2', 'BUT2OLO2', 'IBUTOLBO2', 'IBUTOLCO2', 'TBUTOLO2',
    'HO3C5O2', 'PE2ENEBO2', 'HM2C43O2', 'M2BUOL2O2', 'HM33C3O2', 'ME3BUOLO2', 'HO2M2C4O2', 'ME2BU2OLO2', 'PROL11MO2', 'H2M3C4O2',
    'ME2BUOLO2', 'CYHEXOLAO2', 'MIBKAOHAO2', 'MIBKAOHBO2', 'MIBKHO4O2', 'CH3CO3', 'NMBOAO2', 'NMBOBO2', 'MBOAO2', 'MBOBO2',
    'CH3O2', 'HCOCH2O2', 'C2H5CO3', 'C2H5O2', 'C3H7CO3', 'BUTALO2', 'NC3H7O2', 'IPRCO3', 'IBUTALBO2', 'IBUTALCO2',
    'IC3H7O2', 'C4H9CO3', 'C4CHOBO2', 'NC4H9O2', 'ACO3', 'ACRO2', 'OCCOHCO2', 'CH3C2H2O2', 'MACO3', 'MACRO2',
    'MACROHO2', 'C3DBCO3', 'C4CONO3O2', 'C4NO3COO2', 'C4OCCOHCO2', 'COCCOH2CO2', 'SC4H9O2', 'IC4H9O2', 'TC4H9O2', 'PEAO2',
    'PEBO2', 'PECO2', 'IPEAO2', 'IPEBO2', 'IPECO2', 'NEOPO2', 'HEXAO2', 'HEXBO2', 'HEXCO2', 'M2PEAO2',
    'M2PEBO2', 'M2PECO2', 'M2PEDO2', 'M3PEAO2', 'M3PEBO2', 'M3PECO2', 'M22C43O2', 'M22C4O2', 'M33C4O2', 'M23C43O2',
    'M23C4O2', 'HEPTO2', 'M2HEXAO2', 'M2HEXBO2', 'M3HEXAO2', 'M3HEXBO2', 'OCTO2', 'NONO2', 'DECO2', 'UDECO2',
    'DDECO2', 'CHEXO2', 'ETHENO3O2', 'PRONO3AO2', 'PRONO3BO2', 'BU1ENO3O2', 'C43NO34O2', 'HO3C4O2', 'C42NO33O2', 'MPRANO3O2',
    'MPRBNO3O2', 'C51NO32O2', 'C52NO31O2', 'PE1ENEAO2', 'PE1ENEBO2', 'C52NO33O2', 'C53NO32O2', 'PE2ENEAO2', 'C4NO32M1O2', 'C4NO32M2O2',
    'C4M3NO31O2', 'C4M3NO32O2', 'IPRCHOOA', 'ME3BU2OLO2', 'C4M2NO32O2', 'C4M2NO33O2', 'C65NO36O2', 'C66NO35O2', 'C6OH5O2', 'HO5C6O2',
    'C62NO33O2', 'C63NO32O2', 'C64OH5O2', 'C65OH4O2', 'C4ME2NO3O2', 'C4ME2OHO2', 'BZBIPERO2', 'C6H5CH2O2', 'TLBIPERO2', 'OXYLO2',
    'OXYBIPERO2', 'MXYLO2', 'MXYBIPERO2', 'PXYLO2', 'PXYBIPERO2', 'C6H5C2H4O2', 'EBZBIPERO2', 'PHC3O2', 'PBZBIPERO2', 'PHIC3O2',
    'IPBZBIPRO2', 'TM123BO2', 'TM123BPRO2', 'TM124BO2', 'TM124BPRO2', 'TMBO2', 'TM135BPRO2', 'ETOLO2', 'OETLBIPRO2', 'METLBIPRO2',
    'PETLBIPRO2', 'DM35EBO2', 'DMEBIPRO2', 'DE35TO2', 'DETLBIPRO2', 'NSTYRENO2', 'STYRENO2', 'C6H5CO3', 'C6H5O2', 'CH2CLO2',
    'CHCL2O2', 'CCL3O2', 'CCL3CH2O2', 'TCEOHO2', 'C2CL3OHAO2', 'C2CL3OHBO2', 'C2CL2OHO2', 'DICLETO2', 'CH2OHCL2O2', 'CL2OHCH2O2',
    'CL12PRAO2', 'CL12PRBO2', 'CL12PRCO2', 'CH3CCL2O2', 'CHCL2CH2O2', 'CH2CLCH2O2', 'CH3CHCLO2', 'CHCL2CL2O2', 'CH2CL3O2', 'CHCL3O2',
    'CCLNO3O2', 'CNO3CLO2', 'CCLOHCH2O2', 'CH2OHCCLO2', 'NBUTDAO2', 'NBUTDBO2', 'BUTDAO2', 'BUTDBO2', 'BUTDCO2', 'NISOPO2',
    'ISOP34O2', 'CHOOCH2O2', 'METHACETO2', 'MOCOCH2O2', 'ACETC2H4O2', 'EOCOCH2O2', 'ETHACETO2', 'NPROACEAO2', 'NPROACEBO2', 'NPROACECO2',
    'IPRACBO2', 'IPROACETO2', 'NBUACETAO2', 'NBUACETBO2', 'NBUACETCO2', 'SBUACETAO2', 'SBUACETBO2', 'MCOOTBO2', 'TBOCOCH2O2', 'CH3OCH2O2',
    'DIETETO2', 'ETOC2O2', 'MTBEAO2', 'MTBEBO2', 'DIIPRETO2', 'IPROMC2O2', 'ETBEAO2', 'ETBEBO2', 'ETBECO2', 'MO2EOLAO2',
    'MO2EOLBO2', 'EOX2EOLAO2', 'EOX2EOLBO2', 'PR2OHMOXO2', 'H2C3OCO2', 'BOX2EOHAO2', 'BOX2EOHBO2', 'BOXPROLAO2', 'BOXPROLBO2', 'CH2BRO2',
    'DIBRETO2', 'CH3COCH2O2', 'MEKAO2', 'MEKBO2', 'MEKCO2', 'CO2C54O2', 'MPRKAO2', 'DIEKAO2', 'DIEKBO2', 'MIPKAO2',
    'MIPKBO2', 'HEX2ONAO2', 'HEX2ONBO2', 'HEX2ONCO2', 'HEX3ONAO2', 'HEX3ONBO2', 'HEX3ONCO2', 'HEX3ONDO2', 'MIBKAO2', 'MIBKBO2',
    'MTBKO2', 'CYHXONAO2', 'NAPINAO2', 'NAPINBO2', 'APINAO2', 'APINBO2', 'APINCO2', 'NBPINAO2', 'NBPINBO2', 'BPINAO2',
    'BPINBO2', 'BPINCO2', 'NLIMO2', 'LIMAO2', 'LIMBO2', 'LIMCO2', 'NBCO2', 'BCAO2', 'BCBO2', 'BCCO2',
    'DMMAO2', 'DMMBO2', 'DMCO2', 'CH3SCH2O2', 'HODMSO2', 'ETHOXO2', 'BUT2CO3', 'C3ME3CO3', 'C3ME3CHOO2', 'HOCH2CO3',
    'CH3CHOHCO3', 'IPRHOCO3', 'IPRCHOO', 'BZEMUCCO3', 'BZEMUCO2', 'C5DIALO2', 'NPHENO2', 'PHENO2', 'CRESO2', 'NCRESO2',
    'TLEMUCCO3', 'TLEMUCO2', 'C615CO2O2', 'OXYMUCCO3', 'OXYMUCO2', 'MC6CO2O2', 'NOXYOLO2', 'OXYOLO2', 'MXYMUCCO3', 'MXYMUCO2',
    'C726CO5O2', 'MXYOLO2', 'NMXYOLO2', 'PXYMUCCO3', 'PXYMUCO2', 'C6M5CO2O2', 'NPXYOLO2', 'PXYOLO2', 'EBENZOLO2', 'NEBNZOLO2',
    'EBZMUCCO3', 'EBZMUCO2', 'C715CO2O2', 'NPBNZOLO2', 'PBENZOLO2', 'PBZMUCCO3', 'PBZMUCO2', 'C815CO2O2', 'IPBENZOLO2', 'NIPBNZOLO2',
    'IPBZMUCCO3', 'IPGLOOB', 'IPBZMUCO2', 'C7M15CO2O2', 'NTM123OLO2', 'TM123OLO2', 'TM123MUCO2', 'NTM124OLO2', 'TM124OLO2', 'TM124MUCO3',
    'TM124MUCO2', 'C7CO2M5O2', 'NTM135OLO2', 'TM135OLO2', 'TM135MUCO3', 'TM135MUCO2', 'C7M2CO5O2', 'OETLMUCCO3', 'OETLMUCO2', 'MC7CO2O2',
    'NOETOLO2', 'OETOLO2', 'METLMUCCO3', 'METLMUCO2', 'C826CO3O2', 'METOLO2', 'NMETOLO2', 'PETLMUCCO3', 'PETLMUCO2', 'C7M6CO2O2',
    'NPETOLO2', 'PETOLO2', 'DMEBMUCO3', 'DMEBMUCO2', 'C8M2CO6O2', 'NDMEPHOLO2', 'DMEPHOLO2', 'NDEMPHOLO2', 'DEMPHOLO2', 'DETLMUCO3',
    'DETLMUCO2', 'C9M2CO6O2', 'HMVKAO2', 'HMVKBO2', 'MVKO2', 'CISOPAO2', 'ISOPBO2', 'CISOPCO2', 'ISOPDO2', 'NC526O2',
    'C530O2', 'M3BU3ECO3', 'C45O2', 'NC51O2', 'C51O2', 'CH2CHCH2O2', 'ISOPAO2', 'ISOPCO2', 'MEMOXYCO3', 'EOX2MECO3',
    'ETOMEO2', 'PRONEMOXO2', 'BOXMCO3', 'BOX2MO2', 'BOXPRONAO2', 'BOXPRONBO2', 'C107O2', 'C109O2', 'C96O2', 'NOPINAO2',
    'NOPINBO2', 'NOPINCO2', 'NOPINDO2', 'LIMALAO2', 'LIMALBO2', 'C923O2', 'BCALAO2', 'BCALBO2', 'C136O2', 'BCALCO2',
    'C141O2', 'HOC2H4CO3', 'HOIPRCO3', 'HO13C5O2', 'HO3C4CO3', 'C54O2', 'H2M2C3CO3', 'PROL1MCO3', 'C56O2', 'HO2C43CO3',
    'MIBKCOOHO2', 'NC4OHCO3', 'C4OH2CO3', 'CO2C3CO3', 'HO2C3CO3', 'IBUDIALCO3', 'PROPALO2', 'CO3C4CO3', 'HO1C4O2', 'A2PANOO',
    'HCOCOHCO3', 'HCOCO3', 'MACRNCO3', 'MACRNBCO3', 'CHOMOHCO3', 'CO2H3CO3', 'HO1C5O2', 'HO2C5O2', 'C52O2', 'TBUTCO3',
    'HO1C6O2', 'C5H11CO3', 'HO2C6O2', 'HO3C6O2', 'HO1MC5O2', 'C54CO3', 'HO2MC5O2', 'EIPKAO2', 'EIPKBO2', 'HO2M2C5O2',
    'H1MC5O2', 'M3C4CO3', 'H2MC5O2', 'M2BKAO2', 'M2BKBO2', 'HM33C4O2', 'M22C3CO3', 'HM22C4O2', 'M33C3CO3', 'HM23C4O2',
    'M2C43CO3', 'HO3C76O2', 'CO3C75O2', 'H2M5C65O2', 'C75O2', 'H2M2C65O2', 'H2M4C65O2', 'C710O2', 'H3M3C6O2', 'HO3C86O2',
    'CO3C85O2', 'HO3C96O2', 'C91O2', 'HO3C106O2', 'C101O2', 'HO3C116O2', 'C111O2', 'HO3C126O2', 'C121O2', 'CO1C6O2',
    'NO3CH2CO3', 'PRNO3CO3', 'CO3C4NO3O2', 'HO3C3CO3', 'MPRBNO3CO3', 'C5NO3COAO2', 'C4NO3CO3', 'C5OH2CO4O2', 'C4OHCO3', 'C5NO3CO4O2',
    'C5CONO34O2', 'C43NO3CO3', 'C4MCONO3O2', 'C3MNO3CO3', 'C3M3OH2CO3', 'MC4CONO3O2', 'C65NO36CO3', 'MNO3COC4O2', 'C4COMOH3O2', 'HO5C5CO3',
    'C6NO3CO5O2', 'C6CONO34O2', 'MALDIALCO3', 'EPXDLCO3', 'C3DIALO2', 'MALDIALO2', 'OXYL1O2', 'C5CO14O2', 'OXYLCO3', 'EPXM2DLCO3',
    'C4MCO2O2', 'DM123O2', 'MXYLCO3', 'MXYL1O2', 'C3MCODBCO3', 'EPXMDLCO3', 'C3MDIALO2', 'MXY1O2', 'PXYLCO3', 'PXYL1O2',
    'PXY1O2', 'C6H5CH2CO3', 'EBENZO2', 'C6DCARBBO2', 'PHCOETO2', 'PBENZO2', 'C7DCCO3', 'IPBENZO2', 'IC7DCCO3', 'IPGLOO',
    'TM123BCO3', 'TM123O2', 'EPXKTMCO3', 'C4CO2O2', 'TM124BCO3', 'DM124O2', 'TM124O2', 'TMBCO3', 'DMPHO2', 'C4MCODBCO3',
    'EPXMKTCO3', 'CO24C53O2', 'MPHCOMEO2', 'EPXMEDLCO3', 'C4ECO2O2', 'OET1O2', 'MET1O2', 'PET1O2', 'DMPHCOMO2', 'EMPHCOMO2',
    'EMPHCO3', 'C7CODBCO3', 'EPXEKTCO3', 'C3EDIALO2', 'CO24C63O2', 'CCL3CO3', 'CLETO3', 'CL2OHCO3', 'CL12CO3', 'CLCOCLMEO2',
    'CHCL2CO3', 'CLCOCH2O2', 'CLCOCLO2', 'CCLOHCO3', 'HNMVKO2', 'NC3CO3', 'C42O2', 'HC3CO3', 'C41O2', 'MVKOHAO2',
    'MVKOHBO2', 'HC3CCO3', 'INCO2', 'NC4CO3', 'C510O2', 'C536O2', 'C537O2', 'INAO2', 'C58O2', 'HC4CO3',
    'CHOCOMOXO2', 'ACETMECO3', 'HOACETETO2', 'MECOACETO2', 'ACPRONEO2', 'ACCOETO2', 'ACETC2CO3', 'IPRACBCO3', 'ACBUONEAO2', 'ACBUONEBO2',
    'ACCOC3H6O2', 'SBUACONEO2', 'TBUACCO3', 'MTBEACHOO2', 'MTBEBCO3', 'IPROC21O2', 'IPROMCCO3', 'EIPEO2', 'ETBEACO3', 'ETBECCO3',
    'BOXCOEOLO2', 'BRETO3', 'HO1CO3C4O2', 'BIACETO2', 'HO2CO4C5O2', 'CO23C54O2', 'HOCO3C54O2', 'C53O2', 'C41CO3', 'CO2HOC61O2',
    'CO24C6O2', 'CO25C6O2', 'HO2C4O2', 'C61O2', 'CO23C65O2', 'C6CO3OH5O2', 'C6CO34O2', 'C6HO1CO3O2', 'C3COCCO3', 'PEN2ONE1O2',
    'MIBK3COO2', 'C612O2', 'CO2M33CO3', 'C6COCHOO2', 'CY6DIONO2', 'NC101O2', 'C96CO3', 'C720O2', 'NC91CO3', 'C8BCO2',
    'C918CO3', 'C923CO3', 'C141CO3', 'NBCALO2', 'BCALO2', 'BCSOZO2', 'C151O2', 'C152O2', 'MMFO2', 'MMCFO2',
    'DMSO2O2', 'CHOC4CO3', 'C6DIALO2', 'CHOC4O2', 'CYC6DIONO2', 'CONM2CO3', 'NBZFUO2', 'BZFUO2', 'CATEC1O2', 'MCATEC1O2',
    'MC3CODBCO3', 'C4M2ALOHO2', 'C5DICARBO2', 'NTLFUO2', 'TLFUO2', 'MC4CODBCO3', 'MC5CO2OHO2', 'NOXYFUO2', 'C6OTKETO2', 'OXYFUO2',
    'OXCATEC1O2', 'C5MCO2OHO2', 'NMXYFUO2', 'C23O3MO2', 'MXYFUO2', 'NPXYFUO2', 'MCOCOMOXO2', 'PXYFUO2', 'MXCATEC1O2', 'DMKOHO2',
    'PXCATEC1O2', 'ECATEC1O2', 'C6DICARBO2', 'NEBFUO2', 'BUTALAO2', 'EBFUO2', 'C7CO3OHO2', 'PCATEC1O2', 'C7DCO2', 'NPBFUO2',
    'C4CHOAO2', 'PBFUO2', 'C8CO3OHO2', 'PHCOMEO2', 'IPCATEC1O2', 'IC7DCO2', 'NIPBFUO2', 'IC4CHOAO2', 'IPBFUO2', 'C7MCO3OHO2',
    'T123CAT1O2', 'C7ADCCO3', 'C7ADCO2', 'NTMB1FUO2', 'TMB1FUO2', 'NTMB2FUO2', 'MC6OTKETO2', 'TMB2FUO2', 'C7BDCO2', 'T124CAT1O2',
    'OTCATEC1O2', 'MTCATEC1O2', 'C7EDCO2', 'PTCATEC1O2', 'C7DDCCO3', 'C7DDCO2', 'NMEBFUO2', 'C23O3EO2', 'MEBFUO2', 'EMPHO2',
    'CH3COCCLO2', 'CLCOCCL2O2', 'C527O2', 'C526O2', 'HC4ACO3', 'C58AO2', 'INB1O2', 'INB2O2', 'HPC52O2', 'HC4CCO3',
    'C57AO2', 'C57O2', 'INDO2', 'C59O2', 'C524O2', 'ETHFORMO2', 'IPRMEETO2', 'CHOOMCO3', 'PRONFORMO2', 'PRCOOMCO3',
    'PRCOOMO2', 'BOXCOCHOO2', 'BOXFORMO2', 'PRONOCOPO2', 'BOXCOCOMO2', 'PINALO2', 'C108O2', 'C89CO3', 'C920CO3', 'C920O2',
    'C97O2', 'C85CO3', 'C85O2', 'C719O2', 'C918O2', 'C9DCO2', 'C915O2', 'C917O2', 'NLIMALO2', 'LIMALO2',
    'C729CO3', 'C822CO3', 'C924O2', 'C816CO3', 'NORLIMO2', 'C816O2', 'NLMKAO2', 'LMKAO2', 'LMKBO2', 'C146O2',
    'C131CO3', 'BCLKAO2', 'BCLKBO2', 'BCLKCO2', 'C131O2', 'C147O2', 'C126CO3', 'C136CO3', 'C148O2', 'C1311O2',
    'NC1313O2', 'C1313O2', 'C126O2', 'C144O2', 'C142O2', 'NBCKO2', 'BCKAO2', 'BCKBO2', 'CH3SOO', 'H13C43CO3',
    'C42CO3', 'HOC3H6CO3', 'C3DIOLO2', 'HO2C4CO3', 'HOIBUTCO3', 'C63O2', 'HO3C5CO3', 'C64O2', 'HO2C54O2', 'HO2C54CO3',
    'C66O2', 'CO3C54CO3', 'H2M2C4CO3', 'C67O2', 'C610O2', 'H2M3C4CO3', 'C68O2', 'C69O2', 'C611O2', 'HM33C3CO3',
    'HM22C3O2', 'HM22C3CO3', 'HM2C43CO3', 'C71O2', 'C76O2', 'C77O2', 'C78O2', 'C711O2', 'H3M3C5O2', 'H3M3C5CO3',
    'C82O2', 'C81O2', 'C93O2', 'C92O2', 'HO6C7O2', 'C103O2', 'C102O2', 'HO7C8O2', 'C113O2', 'C112O2',
    'HO8C9O2', 'C123O2', 'C122O2', 'CO1H63O2', 'C3NO3COO2', 'NPHEN1O2', 'NNCATECO2', 'NCATECO2', 'NBZQO2', 'PBZQO2',
    'NPTLQO2', 'PTLQO2', 'NCRES1O2', 'MNNCATECO2', 'MNCATECO2', 'NOXYOL1O2', 'NOXYQO2', 'OXYQO2', 'OXNNCATCO2', 'OXNCATECO2',
    'C534O2', 'NMXYOL1O2', 'NMXYQO2', 'MXYQO2', 'MXNNCATCO2', 'MXNCATECO2', 'NPXYOL1O2', 'NPXYQO2', 'PXYQO2', 'PXNNCATCO2',
    'PXNCATECO2', 'NEBNZ1O2', 'NPEBQO2', 'PEBQO2', 'ENNCATECO2', 'ENCATECO2', 'CO3H4CO3', 'PHCOCOCO2', 'NPBNZ1O2', 'NPPRBQO2',
    'PPRBQO2', 'PNNCATECO2', 'PNCATECO2', 'C5O45OHCO3', 'NIPBNZ1O2', 'NIPRBQO2', 'IPRBQO2', 'IPNNCATCO2', 'IPNCATECO2', 'C4MOHOCO3',
    'NT123L1O2', 'T123NNCTO2', 'T123NCATO2', 'NT124L1O2', 'NTM124QO2', 'TM124QO2', 'T124NNCTO2', 'T124NCATO2', 'C5CO234O2', 'NOETOL1O2',
    'NOETLQO2', 'OETLQO2', 'OTNNCATCO2', 'OTNCATECO2', 'NMETOL1O2', 'NMETLQO2', 'METLQO2', 'MTNNCATCO2', 'MTNCATECO2', 'NPETOL1O2',
    'NPETLQO2', 'PETLQO2', 'PTNNCATCO2', 'PTNCATECO2', 'CO234C65O2', 'H13CO2CO3', 'CO2N3CO3', 'C535O2', 'C58NO3CO3', 'ACCOCOMEO2',
    'ACEETOHO2', 'ACCOMCOMO2', 'ACCOCOETO2', 'MTBEAALCO3', 'C62O2', 'HO13C4O2', 'HM22CO3', 'C6COCHOCO3', 'C5COCHOO2', 'CHOC2H4O2',
    'HCOCH2CO3', 'CY6TRIONO2', 'C6CYTONO2', 'NC102O2', 'C512CO3', 'C89O2', 'C926O2', 'C817CO3', 'C817O2', 'NC826O2',
    'C826O2', 'C729O2', 'LMLKAO2', 'LMLKBO2', 'C116CO3', 'C116O2', 'C129O2', 'C1210O2', 'CH3SOO2', 'C1H4C5CO3',
    'CHOC4OHO2', 'HOC4CHOO2', 'C6145COO2', 'COHM2CO3', 'CO2C4CO3', 'HOBUT2CO3', 'CO3C5CO3', 'CO2C54CO3', 'C65O2', 'CO2M3C4CO3',
    'C72O2', 'CO25C73O2', 'CO25C74O2', 'C712O2', 'C713O2', 'C714O2', 'C84O2', 'C94O2', 'C104O2', 'C114O2',
    'C6H13CO3', 'C124O2', 'MALANHYO2', 'NDNPHENO2', 'DNPHENO2', 'NDNCRESO2', 'DNCRESO2', 'C6O4KETO2', 'NDNOXYOLO2', 'DNOXYOLO2',
    'MMALANHYO2', 'CH3COCO3', 'NDNMXYOLO2', 'DNMXYOLO2', 'TL4OHNO2O2', 'NDNPXYOLO2', 'DNPXYOLO2', 'NDNEBNZLO2', 'DNEBNZLO2', 'NDNPBNZLO2',
    'DNPBNZLO2', 'C61CO3', 'NDNIPBZLO2', 'DNIPBNZLO2', 'C62CO3', 'NDNT123LO2', 'DNT123LO2', 'TM124NO2O2', 'NDNT124LO2', 'DNT124LO2',
    'MXYOHNO2O2', 'NDNOETOLO2', 'DNOETOLO2', 'NDNMETOLO2', 'DNMETOLO2', 'NDNPETOLO2', 'DNPETOLO2', 'CO356OCO2', 'C531O2', 'INCNCO3',
    'IEACO3', 'IECCO3', 'HPC52CO3', 'INDHCO3', 'C57NO3CO3', 'INAHPCO3', 'INANCO3', 'INAHCO3', 'NC524O2', 'C525O2',
    'HMACO3', 'HMACRO2', 'ACCOMECO3', 'IPRFORMO2', 'PRCOFORMO2', 'PRONOCOMO2', 'CO23C4CO3', 'C5CO34CO3', 'C106O2', 'C717O2',
    'C811CO3', 'C921O2', 'C98O2', 'C86O2', 'C919O2', 'C914O2', 'C916O2', 'C88CO3', 'C88O2', 'C512O2',
    'C619O2', 'C626CO3', 'C626O2', 'C735O2', 'C822O2', 'C823CO3', 'C925O2', 'C622CO3', 'C1011CO3', 'C1210CO3',
    'C132O2', 'C137CO3', 'C1013CO3', 'C1312O2', 'C127O2', 'C143O2', 'CH3SO2O2', 'HO24C5O2', 'C55O2', 'C67CO3',
    'H3M2C4CO3', 'C79O2', 'H3M3C4CO3', 'H13M3C5O2', 'HO4C5CO3', 'HO5C6CO3', 'HO6C7CO3', 'HO7C8CO3', 'HO8C9CO3', 'C5CO2OHCO3',
    'C6CO2OHCO3', 'C5M2OHOCO3', 'C4COMOHCO3', 'C23O3MCO3', 'C23O3CCO3', 'C7CO2OHCO3', 'C6MOHCOCO3', 'C7OHCO2CO3', 'ECO3CO3', 'C8OHCO2CO3',
    'C8CO2OHCO3', 'NDMMALYO2', 'DMMALYO2', 'C7MOHCOCO3', 'C5MEJCO3', 'C6EO2OHCO3', 'C7MJPCO3', 'C23O3ECO3', 'EMPOHNO2O2', 'C47CO3',
    'INB1HPCO3', 'INB1NACO3', 'INB1NBCO3', 'MMALNACO3', 'MMALNBCO3', 'INDHPCO3', 'INANCOCO3', 'HIEB1O2', 'HIEB2O2', 'HO13C3CO3',
    'C5CO23O2', 'CHOC2CO3', 'CHOC3COCO3', 'C5124COCO3', 'CO235C6CO3', 'NC71O2', 'C811O2', 'CHOC3COO2', 'H3C25C6CO3', 'H3C25C6O2',
    'C810O2', 'C818O2', 'C727CO3', 'NC728O2', 'C728O2', 'C622O2', 'C823O2', 'C819O2', 'C731CO3', 'C1011O2',
    'C137O2', 'C1013O2', 'C1010O2', 'C117O2', 'C830CO3', 'C145O2', 'C927O2', 'C1214O2', 'CHOC4DOLO2', 'C6TRONOHO2',
    'C23C54CO3', 'C73O2', 'C74O2', 'C715O2', 'C83O2', 'C95O2', 'C105O2', 'C115O2', 'C125O2', 'C4CO2DBCO3',
    'C5CO2DBCO3', 'C4DBM2CO3', 'C5DBCO2CO3', 'C7CO2DBCO3', 'C8CO2DBCO3', 'C8DBCO2CO3', 'C4DBMECO3', 'C5DBECO3', 'C5EDBCO3', 'C31CO3',
    'C533O2', 'MECOFORMO2', 'C5124COO2', 'CO235C6O2', 'C716O2', 'C922O2', 'C614O2', 'C511O2', 'C620O2', 'C87CO3',
    'C616O2', 'C718CO3', 'C513O2', 'CO25C6CO3', 'C627O2', 'C727O2', 'C511CO3', 'C517CO3', 'C517O2', 'C628O2',
    'C824O2', 'C1211CO3', 'C133O2', 'C830O2', 'C128O2', 'HO24C4CO3', 'C613O2', 'CO2OH3MCO3', 'C812O2', 'C721CO3',
    'C721O2', 'H3C2C4CO3', 'C87O2', 'C718O2', 'C514O2', 'C820O2', 'C518CO3', 'NC623O2', 'C623O2', 'C825O2',
    'C731O2', 'C732CO3', 'C1012O2', 'C1211O2', 'C139O2', 'C1014O2', 'C736O2', 'C118O2', 'C928CO3', 'C630O2',
    'C1215O2', 'EMALANHYO2', 'PMALANHYO2', 'IPMALNHYO2', 'C312COCO3', 'CHOCOCH2O2', 'NC72O2', 'C621O2', 'C515CO3', 'C515O2',
    'C821O2', 'HMVKBCO3', 'C520O2', 'C624CO3', 'C732O2', 'C829O2', 'C134O2', 'C827CO3', 'C522CO3', 'C831O2',
    'C928O2', 'C46CO3', 'C930O2', 'C813O2', 'C722O2', 'C615CO3', 'C617CO3', 'C618CO3', 'C617O2', 'C618O2',
    'NC730O2', 'C730O2', 'C624O2', 'C733O2', 'C1212O2', 'C827O2', 'C1310O2', 'NC61CO3', 'C615O2', 'C519CO3',
    'C629O2', 'C734O2', 'C521O2', 'C135O2', 'COO2C4CO3', 'COO2C4O2', 'C929O2', 'C516O2', 'C44O2', 'H1C23C4CO3',
    'H1C23C4O2', 'CO1M22CO3', 'C519O2', 'C625O2', 'C1213O2', 'COO2C3CO3', 'C828CO3', 'C828O2', 'HYDFURANO2', 'HYDMEFURANO2',
    'C5HYDCARBO2', 'HYDDIMEFURANO2', 'FURFURALO2', 'HYDFURFURALO2', 'MEFURFURALO2', 'HYDMEFURFURALO2', 'ALDFURFURALO2', 'FURANO2', 'CARBFURANO2', 'GUAIACOLO']
    #endregion
    # #Check the differences in my list vs official RO2ToAdd from mechanism
    # diff_ro2_lists = list(set(mcm_bb_ro2_list) ^ set(ro2_species_smarts_list))
    # print(diff_ro2_lists)
    #Only GUAIACOLO is different between the lists. I think that whoever made the mechanism file put GUAIACOLO as an RO2 incorrectly? So we are good to use our ro2_species_smarts_list


#print(nitrates_and_nox_reservoirs)
#     print(len(nitrates_and_nox_reservoirs))
#     
		
#CALL FUNCTIONS
get_info_on_missing_species()
subtract_mcm_species_minus_extra()

path_all_species = dirpath + 'Mechanism_info/mcm_species_bb_sherwen_total_info.xlsx'
df_allspecies = pd.read_excel(path_all_species, index_col = 0)
group_all_species(
    df_in = df_allspecies, 
    use = 'Name', 
    smart_groups = smart_groups)

choose_smarts_classifications()