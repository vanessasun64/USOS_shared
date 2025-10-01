import os 
import sys 
import re 
import yaml
import inspect 
import numpy as np 
import pandas as pd
import xarray as xr
import pubchempy as pcp
import rdkit
from rdkit import Chem
from rdkit.Chem import Fragments
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem

from collections import defaultdict

#get path for USOS_shared directory for correct laptop
#input parameter should be 'CHPC', 'Mac', or 'Windows'
current_dir = os.path.dirname(os.path.abspath(__file__))
global_scripts_path = os.path.abspath(os.path.join(current_dir, "..", "global_scripts"))
sys.path.insert(0,global_scripts_path)
from dirpath import filepath_source
dirpath = filepath_source('CHPC')

#define path for Hawthorne data directory
savepath= dirpath + 'Hawthorne_data/'

#import CRACMM utilities and mapper
utildir = dirpath +'CRACMM/utilities/'
sys.path.append(utildir)

import cracmm1_mapper as cracmm1   # includes: get_cracmm_roc(smiles,koh,log10cstar) (Version 1)
import cracmm2_mapper as cracmm2   # includes: get_cracmm_roc(smiles,koh,log10cstar) (Version 2)

##########################################################################################################################################################
#FUNCTIONS

def hawthorne_usos_extract_species(udaq_voc_filepath, epa_parameter_info, filename):
    """Function to match UDAQ's 2024 VOC measurements for Hawthorne during the USOS Campaign to the species using EPA AQS Parameter codes.
    
    INPUTS: 
    -------  
    (1) udaq_voc_filepath - STR of full path to the file from UDAQ (called Verbose.csv as of Sept. 2025)

    (2) epa_parameter_info - STR of full path to CSV file of EPA's AQS Parameter codes, downloaded from:
                         https://aqs.epa.gov/aqsweb/documents/codetables/parameters.html - last downloaded on Sept. 27, 2025

    (3) filename - STR of full path of the file to save as output

    OUTPUTS
    ------- 
    Dataframe of EPA 'Parameter Codes', 'EPA Species', 'Parameter Alternate Name', 'CAS Number', and 'Standard Units', saved as file types:
    (1) Excel file
    (2) CSV file

    Keep in mind that these are ORDERED DataFrames

    """
    #Read UDAQ VOC measurements file
    df_data = pd.read_csv(udaq_voc_filepath) 
    #Get only the data at Hawthorne
    df_hawthorne_usos_only = df_data.loc[(df_data['StationSym'] == 'HW')]
    #set index to datetimeindex
    df_hawthorne_usos_only.set_index(['dt'], inplace = True) 
    df_hawthorne_usos_only.index = pd.to_datetime(df_hawthorne_usos_only.index)

    #Get only the data during the USOS campaign
    df_hawthorne_usos_only = df_hawthorne_usos_only.sort_index().loc['2024-07-14 00:00:00':'2024-08-18 23:00:00']
    #View Parameter code values as only an int (they print with a .0 at the end initially)
    #Sort and only include the unique parameter codes
    df_hawthorne_usos_only.Parameter = df_hawthorne_usos_only.Parameter.astype(int)
    sorted_parameters = sorted(df_hawthorne_usos_only.Parameter.unique())

    #Save the Parameter codes as a new DataFrame
    df_species = pd.DataFrame(sorted_parameters, columns = ['Parameters'])

    #Read CSV file of the EPA Parameters
    df_epa_parameters = pd.read_csv(epa_parameter_info)

    #Match the parameter codes from the UDAQ VOC measurements to the EPA Parameter Code, drop unnecessary columns
    matching_parameters = df_epa_parameters[df_epa_parameters['Parameter Code'].isin(df_species['Parameters'])]
    df_matching_species = matching_parameters.sort_values(by='Parameter Code').reset_index(drop=True)
    df_matching_species = df_matching_species.drop(columns = ['Parameter Abbreviation', 'Still Valid', 'Round or Truncate'])
    df_matching_species = df_matching_species.rename(columns={'Parameter':'EPA Species'})

    # Write the dataframe to an excel file and csv file:
    excel_file = filename + '.xlsx'
    csv_file = filename + '.csv'
    df_matching_species.to_excel(excel_file, index = False)
    df_matching_species.to_csv(csv_file, index = False)
    return print('Output df of parameters saved at: ' + excel_file)
    return print('Output df of parameters saved at: ' + csv_file)

def csv_retrieve_compounds(data_filepath, parameter_column, species_column, cas_column):
    """Function to read CSV file of species from UDAQ's 2024 VOC measurements for Hawthorne during the USOS Campaign. 
    Requires the CSV output (called Hawthorne_EPA_match file.csv) from the hawthorne_usos_extract_species function that matches the original UDAQ file without species name data to the EPA parameters.

    INPUTS: 
    -------  
    (1) data_filepath - STR of full path to the file from hawthorne_usos_extract_species output (Hawthorne_EPA_match file.csv)
    (2) parameter_column - STR of the column name for Parameter Codes in Hawthorne_EPA_match.csv file
    (3) species_column - STR of the column name for EPA Species in Hawthorne_EPA_match.csv file
    (4) cas_column - STR of the column name for CAS Numbers in Hawthorne_EPA_match.csv file

    OUTPUTS
    ------- 
    (1) parameter_codes_list - LIST of Parameter Codes that we have for UDAQ VOC measurements taken at Hawthorne during 2024 USOS Campaign
    (2) species_list - LIST of EPA species that we have for UDAQ VOC measurements taken at Hawthorne during 2024 USOS Campaign
    (3) cas_list - LIST of CAS Numbers that we have for UDAQ VOC measurements taken at Hawthorne during 2024 USOS Campaign

    Keep in mind that these are ORDERED lists

    """
    parameter_codes_list = []
    species_list = []
    cas_list = []
    df_data = pd.read_csv(data_filepath)
    #for only testing a few species:
    # for few_spec in range(50,53):
    #     parameter_codes_list.append(df_data[parameter_column][few_spec])
    #     species_list.append(df_data[species_column][few_spec])
    #     cas_list.append(df_data[cas_column][few_spec])
    for spec in range(0,len(df_data)):
        parameter_codes_list.append(df_data[parameter_column][spec])
        species_list.append(df_data[species_column][spec])
        cas_list.append(df_data[cas_column][spec])
    return parameter_codes_list, species_list, cas_list

def get_chem_info_epa(cmpd_list, parameters_list, cas_hawthorne_list, return_df, quiet,
                  out2yaml, out2excel, filename, overwrite):
    
    """Function to retrieve PubChem info on UDAQ's 2024 VOC measurements for Hawthorne during the USOS Campaign and save a new Excel & CSV file of info.
    Includes:
        - EPA Species Name
     	- query_results: confirms that our species was found by CID match
     	- iupac_name
     	- synonyms: synonyms in PubChem for each species
     	- cid: PubChem chemical ID number
     	- molecular weight
     	- chemical formula
     	- inchi
     	- inchi_key
     	- SMILES
     	- cas_number

    
    INPUTS: 
    -------  
      (1) cmpd_list - name of LIST containing chemical compound names from UDAQ's 2024 VOC measurements for Hawthorne during the USOS Campaign (from csv_retrieve_compounds output)

      (2) parameters_list - name of LIST containing parameter codes from UDAQ's 2024 VOC measurements for Hawthorne during the USOS Campaign (from csv_retrieve_compounds output)
                    
      (3) cas_hawthorne_list - name of LIST containing CAS numbers from UDAQ's 2024 VOC measurements for Hawthorne during the USOS Campaign (from csv_retrieve_compounds output)
      
      (4) return_df - (OPTIONAL) BOOL indicating if output dict, should be 
                  returned as a pandas dataframe rather than dictionary (default). 
                  Default is  set to FALSE, as dictionary is easier for parsing later, 
                  but as df, output is much more read-able for humans...  

      (5) quiet    - (OPTIONAL) BOOL of whether or not to print additional output. 
                     Default is True. 
                     
      (6) out2yaml - (OPTIONAL) BOOL indicating if outputs should be saved to a 
                     yaml file. Default is FALSE.  
       
      (7) out2excel - (OPTIONAL) BOOL indicating if outputs should be saved to 
                      an excel .xlsx sheet. Also saves CSV. Default is FALSE. 


       ------ below only used if either out2yaml or out2excel is true ---------
          
      
      (8) filename - (OPTIONAL) STR with name of the output file to write. Do NOT include extension.
      
      (9) overwrite - (OPTIONAL) BOOL indicating whether or not output files 
                        should overwrite any exisiting files at output_dir
                        with output_fname or not. If FALSE, a version # is added. 
       
    OUTPUTS:  
    -------- 
    (1) info - DICT or DF with info on all compounds in list including NA for 
               items w/o matches. 
               If out2excel is True, info will output to Excel and CSV file.

    #AFTER CALLING, MAKE SURE TO EDIT M AND P XYLENE MANUALLY!
    Index   Parameter   EPA Species query_results       iupac_name  synonyms    cid     mw      formula inchi                                           inchi_key                   SMILES          cas_number
    45      45109	    m-Xylene	Entered manually	1,3-xylene	            7929	106.16	C8H10	InChI=1S/C8H10/c1-7-4-3-5-8(2)6-7/h3-6H,1-2H3	IVSZLXZYQVIEFR-UHFFFAOYSA-N	CC1=CC(=CC=C1)C	108-38-3
    46      45109	    p-Xylene	Entered manually	1,4-xylene	            7809	106.16	C8H10	InChI=1S/C8H10/c1-7-3-5-8(2)6-4-7/h3-6H,1-2H3	URLKBWYHVLBVBO-UHFFFAOYSA-N	CC1=CC=C(C=C1)C	106-42-3


    AUTHOR: 
    -------
       Vanessa Sun (vanessa.sun@utah.edu)
       Github:@vanessasun64

       Adapted from get_SMILES_form_name.py file from:
       Prof. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD
    """
    
    # Intialize output vars: 
    info={}           # Output Dict to store info on all compounds 
    no_matches=[]     # List to store names of cmpds with no matches. 
    multi_matches={}; # Dict to store query results about cmpds with multiple matches...

    for i, name in enumerate(cmpd_list): 
        # Initialize output dict w/ blanks: 
        info[name]= {'Parameter Code':'', 'EPA Species':'','query_results':'', 'iupac_name':'',
                'synonyms':'', 'cid':None, 'mw':np.nan,
                'formula':'', 'inchi':'', 'inchi_key':'', 'SMILES':'', 'cas_number':''}
        
        #Fill Parameter Codes and EPA Species names
        info[name]['Parameter Code'] =  parameters_list[i]
        info[name]['EPA Species']=str(name)   

        # Print progress to the screen (it's slow). 
        if quiet is False:  print(f'Parsing {i+1}/{len(cmpd_list)}... {name}')
        
        # First, get the chemical IDs of all potential matches to this name. 
        cids= pcp.get_cids(name, 'name', 'substance', list_return='flat')
        
        # And use that to get a list of all potential compound objects: 
        matches=[pcp.Compound.from_cid(cid) for cid in cids]

        print(matches)
        # PubChem has a filtered "whitelist" with  human-chosen CIDS for popular 
        # names. Go ahead and figure out what they think the info is/should be: 
        pref_match = pcp.get_compounds(name, 'name')
        
        #No Match from EPA Species Name to PubChem entry scenario: 
        if len(matches)==0 or len(pref_match)==0: 
            print('Could not find a match by EPA Species name ' + str(name) + '. \n Might have to enter manually')
        
            # If you didn't find any matches, then print/store that info. 
            no_matches.append(name)
            info[name]['query_results']='No Matches. Try CAS manually or invalid.'

            if quiet is False: 
                print(f'\tNo matches found for: "{name}"')

        #Scenario if there are matches
        else: 
            #Check if the CAS Number of the preferred match is the same as the species
            #Since Pubchempy doesn't allow you to match with CAS numbers, we loop through the synonyms for each match to find the CAS number
            #and add to info. Save info like query_results to confirm that we matched by CAS number.
            for syn in pref_match[0].synonyms:
                if syn == cas_hawthorne_list[i]:
                    print('Found CAS Number in preferred match')
                    info[name]['query_results']='Found by CAS Number in preferred match.'
                    info[name]['iupac_name']=pref_match[0].iupac_name
                    info[name]['cid']=str(pref_match[0].cid)
                    info[name]['synonyms']=str(pref_match[0].synonyms)
                    info[name]['cid']=int(pref_match[0].cid)
                    info[name]['mw']=float(pref_match[0].molecular_weight)
                    info[name]['formula']=str(pref_match[0].molecular_formula)
                    info[name]['inchi'] =str(pref_match[0].inchi)
                    info[name]['inchi_key'] =str(pref_match[0].inchikey)
                    info[name]['SMILES']=str(pref_match[0].smiles)
                    info[name]['cas_number'] = str(syn)

                    preferred_match_bool = True
                    break
            
            if preferred_match_bool == False:
                for cmpd in matches:
                    for syn in cmpd.synonyms:
                        if syn == cas_hawthorne_list[i]:
                            print('Found CAS Number in non-preferred match')
                            info[name]['query_results']='Found by CAS Number in non-preferred match.'
                            info[name]['iupac_name']=cmpd.iupac_name
                            info[name]['cid']=str(cmpd.cid)
                            info[name]['synonyms']=str(cmpd.synonyms)
                            info[name]['cid']=int(cmpd.cid)
                            info[name]['mw']=float(cmpd.molecular_weight)
                            info[name]['formula']=str(cmpd.molecular_formula)
                            info[name]['inchi'] =str(cmpd.inchi)
                            info[name]['inchi_key'] =str(cmpd.inchikey)
                            info[name]['SMILES']=str(cmpd.smiles)
                            info[name]['cas_number'] = str(syn)
                            break

    ###########################################################################
    #          Store & return all that nice output how they want... 
    ###########################################################################
    # if out2yaml==True: 
    #     # Check that the desired filename is available & return the full path to write the file: 
    #     yaml_file=check_filename(filename=filename, default_name= 'species_info', ext='.yml', 
    #                              savepath=savepath, overwrite=overwrite, return_full=True, quiet=quiet)
    #     # Write the dictionary to a YAML file
    #     with open(yaml_file, 'w') as f:   
    #         yaml.dump(info, f)
    #     print(f'Output dictionary of all Pubchem matches saved at: \n\t{yaml_file}')
    print(info)
    if return_df== True or out2excel== True: 
        #Our info is saved as a dictionary. However, if we convert to a DataFrame now, the index will be the species name.
        #We convert the index to be numbered instead
        combined_dict = defaultdict(list)
        print('info.values: ')
        print(info.values())

        for inner_dict in info.values():
            print('inner_dict: ')
            print(inner_dict)
            for key, value in inner_dict.items():
                combined_dict[key].append(value)
                print('combined_dict[key]: ')
                print(combined_dict[key])

        combined_dict = dict(combined_dict)
        print(combined_dict)
        
        # Convert dict to dataframe with all this information... 
        info_df = pd.DataFrame.from_dict(combined_dict)

        #Save info into Excel and CSV file
        if out2excel==True: 
            excel_file = savepath + filename + '.xlsx'
            csv_file =  savepath + filename + '.csv'
            # Write the dataframe to an excel file:
            info_df.to_excel(excel_file, index = False)
            info_df.to_csv(csv_file, index = False)
            print('Output df of PubChem matches saved at:' +  excel_file)
            print('Output df of PubChem matches saved at: ' + csv_file)
    
    if return_df== True:
        return info_df
    else: 
        return info

def comptox_extract(comptox_file, pubchem_match_file, save_filename):
    """
    #BEFORE CALLING, MAKE SURE TO EDIT M AND P XYLENE MANUALLY IN EPA_Hawthorne_pubchem_match.xlsx and CSV
    Index   Parameter   EPA Species query_results       iupac_name  synonyms    cid     mw      formula inchi                                           inchi_key                   SMILES          cas_number
    45      45109	    m-Xylene	Entered manually	1,3-xylene	            7929	106.16	C8H10	InChI=1S/C8H10/c1-7-4-3-5-8(2)6-7/h3-6H,1-2H3	IVSZLXZYQVIEFR-UHFFFAOYSA-N	CC1=CC(=CC=C1)C	108-38-3
    46      45109	    p-Xylene	Entered manually	1,4-xylene	            7809	106.16	C8H10	InChI=1S/C8H10/c1-7-3-5-8(2)6-4-7/h3-6H,1-2H3	URLKBWYHVLBVBO-UHFFFAOYSA-N	CC1=CC=C(C=C1)C	106-42-3
    
    -----

    Function to get kOH and log10 C star from EPA CompTox dashboard for UDAQ's 2024 VOC measurements for Hawthorne during the USOS Campaign
    then use the CRACMM Mapper to figure out which CRACMM species the measurements correspond to. Saves Excel and CSV file.

    INPUTS
    -----
    (1) comptox - STR of the path to the CSV file where we exported data from EPA CompTox by copying InChI keys of our species from UDAQ's 2024 VOC measurements 
    for Hawthorne during the USOS Campaign (See Step 5 of Hawthorne VOC Workflow explanation file)
    (2) pubchem_match_file - STR of the path to the CSV file where we extracted data from PubChem on our species from 
    UDAQ's 2024 VOC measurements for Hawthorne during the USOS Campaign. This is the output from running the get_chem_info_epa function. 
    (3) save_filename - STR of only the file name we want to save CRACMM mapping for UDAQ's 2024 VOC measurements for Hawthorne during the USOS Campaign.

    OUTPUTS
    ------
    Excel and CSV file including the following for each species that we have measurements for, from UDAQ at Hawthorne, during the USOS Campaign:
        - Parameter Code
        - EPA Species
        - query_results
        - iupac_name
        - synonyms
        - cid	
        - mw
        - formula
        - inchi
        - inchi_key
        - SMILES
        - cas_number
        - CompTox Name
        - kOH (cm3/molecule*s)
        - Vapor Pressure 25C (mmHg)
        - Vapor Pressure 25C (Pa)
        - Cstar
        - log10Cstar
        - CRACMM Mapping

    """
    #Read in the EPA CompTox CSV file
    df_comptox_data = pd.read_csv(comptox_file)
    #Add a row to act as Total NMVOCs row
    new_row= {'PREFERRED_NAME': "Total NMVOCs"}
    new_row_df = pd.DataFrame([new_row])
    df_comptox_data = pd.concat([new_row_df, df_comptox_data], ignore_index=True)

    #Select which columns of data we need, keep Preferred name for manually checking that merge between CompTox info and EPA Species worked properly
    columns_to_keep = ['PREFERRED_NAME', 'ATMOSPHERIC_HYDROXYLATION_RATE_(AOH)_CM3/MOLECULE*SEC_OPERA_PRED','VAPOR_PRESSURE_MMHG_OPERA_PRED']
    #renaming columns to some shorter names
    df_comptox_trimmed = df_comptox_data[columns_to_keep]
    #convert vapor pressure from mmHg to Pa
    df_comptox_trimmed = df_comptox_trimmed.rename(columns = {'PREFERRED_NAME':'CompTox Name','ATMOSPHERIC_HYDROXYLATION_RATE_(AOH)_CM3/MOLECULE*SEC_OPERA_PRED':'kOH (cm3/molecule*s)', 'VAPOR_PRESSURE_MMHG_OPERA_PRED':'Vapor Pressure 25C (mmHg)'})
    df_comptox_trimmed['Vapor Pressure 25C (Pa)'] = df_comptox_trimmed['Vapor Pressure 25C (mmHg)'].multiply(133.322)

    #Read in our PubChem matches
    df_pubchem_match_data = pd.read_csv(pubchem_match_file)
    #Merge PubChem matches with EPA CompTox file to new dataframe
    df_merge_comptox_pubchem = pd.merge(df_pubchem_match_data, df_comptox_trimmed, how='outer', left_index=True, right_index=True)

    #Calculate the CStar and log 10 CStar from CRACMM utilities using vapor pressure and molecular weight, then add to our dataframe
    R = 8.314
    T = 298
    df_merge_comptox_pubchem['Cstar'] = (df_merge_comptox_pubchem['Vapor Pressure 25C (Pa)'] * df_merge_comptox_pubchem['mw']*10**6) / (R*T)
    df_merge_comptox_pubchem['log10Cstar'] = np.log10(df_merge_comptox_pubchem['Cstar'])

    #Use CRACMM Mapper to get the suggested name in CRACMM mechanism using SMILES, kOH, and log10 CStar values
    cracmm_mapping_list = [np.nan] #First row of DataFrame is Total NMVOCs so we want the first CRACMM mapping value to be blank or NaN
    for idx in range(1,len(df_merge_comptox_pubchem)):
        smiles = df_merge_comptox_pubchem['SMILES'][idx]
        kOH = df_merge_comptox_pubchem['kOH (cm3/molecule*s)'][idx]
        log10Cstar = df_merge_comptox_pubchem['log10Cstar'][idx]
        cracmm_mapping = cracmm2.get_cracmm_roc(smiles, kOH, log10Cstar)
        print("EPA Species is " + df_merge_comptox_pubchem['EPA Species'][idx] + ', mapped to CRACMM Species: ')
        print(cracmm_mapping)
        cracmm_mapping_list.append(cracmm_mapping)
    #Add CRACMM Mapping to our DataFrame
    df_merge_comptox_pubchem['CRACMM Mapping'] = cracmm_mapping_list
    excel_file = savepath + save_filename + '.xlsx'
    csv_file = savepath +  save_filename + '.csv'

    #Write the dataframe to an excel file and csv file:
    df_merge_comptox_pubchem.to_excel(excel_file, index = False)
    df_merge_comptox_pubchem.to_csv(csv_file, index = False)
    print('Output df of parameters saved at: ' + excel_file)
    print('Output df of parameters saved at: ' + csv_file)


###############################################################################################
#CALL FUNCTIONS

# hawthorne_usos_extract_species(
#     udaq_voc_filepath = dirpath + 'Hawthorne_data/Verbose.csv',
#     epa_parameter_info = dirpath + 'Hawthorne_data/epa_parameters_official.csv',
#     filename='Hawthorne_EPA_match'
# )

# parameter_codes_list, species_list, cas_list = csv_retrieve_compounds(
#     data_filepath = dirpath + 'Hawthorne_data/Hawthorne_EPA_match.csv',
#     parameter_column = 'Parameter Code',
#     species_column = 'EPA Species',
#     cas_column = 'CAS Number'
# )

# get_chem_info_epa(
#     cmpd_list = species_list,
#     parameters_list= parameter_codes_list,
#     cas_hawthorne_list= cas_list,
#     return_df=False,
#     quiet = False,
#     out2yaml = False,
#     out2excel = True,
#     filename = 'EPA_Hawthorne_pubchem_match',
#     overwrite = False)

comptox_extract(
    comptox_file = savepath + 'comptox_batch_search_all.csv',
    pubchem_match_file = savepath + 'EPA_Hawthorne_pubchem_match.csv',
    save_filename= 'EPA_CRACMM_mapped'     
)

