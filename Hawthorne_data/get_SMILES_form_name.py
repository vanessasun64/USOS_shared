import os 
import sys 
import re 
import yaml
import inspect 
import numpy as np 
import pandas as pd 

import pubchempy as pcp
import rdkit
from rdkit import Chem
from rdkit.Chem import Fragments
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem

###############################################################################
# Functions used to actually query info: 
############################################################################### 

def csv_retrieve_compounds(data_filepath):
    df_data = pd.read_csv(data_filepath)

csv_retrieve_compounds(
    data_filepath = '/Hawthorne_data/'
)

cmpd_list=['chlorine'] 
return_df=False
quiet=False
out2yaml=False
out2excel=False
savepath='/uufs/chpc.utah.edu/common/home/u6044586/'
filename='BYU_species_info'
overwrite=False

def get_chem_info(cmpd_list, return_df:bool=False, quiet:bool=False,
                  out2yaml:bool=False, out2excel:bool=False, 
                  savepath:str='',filename:str='species_info', overwrite:bool=False):
    
    """Function to take a list of chemical compound's names  & try to retrieve
    info like formula, InChI, molecular weight, etc. from PubChem's database.
    Can return output of all matches as dict of df and save as either yaml or .xlsx.
    
    INPUTS: 
    -------  
      (1) cmpd_list - Iterable object containing chemical compound names
                    
      (2) return_df - (OPTIONAL) BOOL indicating if output dict, should be 
                  returned as a pandas dataframe rather than dictionary (default). 
                  Default is  set to FALSE, as dictionary is easier for parsing later, 
                  but as df, output is much more read-able for humans...  

      (3) quiet    - (OPTIONAL) BOOL of whether or not to print additional output. 
                     Default is True. 
                     
      (3) out2yaml - (OPTIONAL) BOOL indicating if outputs should be saved to a 
                     yaml file. Default is FALSE.  
       
      (4) out2excel - (OPTIONAL) BOOL indicating if outputs shoudl be saved to 
                      an excel .xlsx sheet. Default is FALSE. 

       ------ below only used if either out2yaml or out2excel is true ---------
          
      (5) savepath - (OPTIONAL) STR containing the path where the output file
                        should be saved.Default is location of this file. 
      
      (6) filename - (OPTIONAL) STR with name of the output file to write. 
                         Default is 'species_fino.{ext}'
      
      (7) overwrite - (OPTIONAL) BOOL indicating whether or not output files 
                        should overwrite any exisiting files at output_dir
                        with output_fname or not. If FALSE, a version # is added. 
       
    OUTPUTS:  
    -------- 
    (1) info - DICT or DF with info on all compounds in list including NA for 
               items w/o matches. 
               
    (2) no_matches -LIST of all compounds in input list with no match found. 
       
    (3) multi_matches - DICT with keys as names of compounds multiple matches
                        were found for & values as the results of the PubChem 
                        queries for all potential matches 

    AUTHOR: 
    -------
       Prof. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD
    """
    
    # Intiatlize output vars: 
    info={}           # Output Dict to store info on all compounds 
    no_matches=[]     # List to store names of cmpds with no matches. 
    multi_matches={}; # Dict to store query results about cmpds with multiple matches...

    for i, name in enumerate(cmpd_list): 
        
        # Initialize output dict w/ blanks: 
        info[name]= {'query_results':'', 'iupac_name':'',
                'synonyms':'', 'cid':None, 'mw':np.nan,'charge':np.nan,
                'formula':'', 'inchi':''} 
            
        # Print progress to the screen (its slow). 
        if quiet is False:  print(f'Parsing {i+1}/{len(cmpd_list)}... {name}')
        
        # First, get the chemical IDs of all potential matches to this name. 
        cids= pcp.get_cids(name, 'name', 'substance', list_return='flat')
        
        # And use that to get a list of all potential compound objects: 
        matches=[pcp.Compound.from_cid(cid) for cid in cids]
        
        # PubChem has a filtered "whitelist" with  human-chosen CIDS for popular 
        # names. Go ahead and figure out what they think the info is/should be: 
        pref_match = pcp.get_compounds(name, 'name')
        
        if len(matches)==0 or len(pref_match)==0: 
            # If you didn't find any matches, then print/ store that info. 
            no_matches.append(name)
            info[name]['query_results']='No Matches'
            
            if quiet is False: 
                print(f'\tNo matches found for: "{name}"')
        else: 
            for cmpd in matches: 
                if cmpd.inchi: 
                    # Create an RDKit molecule from this INCHI. 
                    mol = Chem.MolFromInchi(cmpd.inchi)
                                      
            if cmpd.cid in cids: print('CID of pref match in items')
            
            # disp=[f"Match {i}/{len(matches)} for INPUT_CMPD='{name}'\n",
            # f"\t{'iupac_name :':>10} {cmpd.iupac_name}\n",
            # f"\t{'synonyms :':>10} {';'.join(cmpd.synonyms)}\n",
            # f"\t{'mw :':>10} {float(cmpd.molecular_weight)}\n",
            # f"\t{'charge :':>10} {float(cmpd.charge)}\n",
            # f"\t{'formula :':>10} {cmpd.molecular_formula}\n",
            # f"\t{'inchi :':>10}  {cmpd.inchi}\n"]
            
            # print(s for s in disp)

            elif len(results) > 1: 
                # If you got multiple matches also print/store that info. 
                multi_matches[name] =results
                
                # And summarize/ place info about this result in info: 
                matched_to=';'.join(['{m.iupac_name} ({m.formula})' for m in results])
                info[name]['query_results']='Multiple matches: {matched_to}'
            
                if quiet is False: 
                    print(f'\tMultiple potential matches found for: "{name} \n\t {matched_to}"')
                        
            elif len(results)==1: 
                # If you ONLY got 1 match for this name... then store all its info: 
                cmpd=results[0]
                
                # Pull about it from built-in methods & store in output diiionary 
                info[name]['query_results']='Uniquely Matched'
                info[name]['iupac_name']=cmpd.iupac_name
                info[name]['synonyms']=cmpd.synonyms
                info[name]['cid']=int(cmpd.cid)
                info[name]['mw']=float(cmpd.molecular_weight)
                info[name]['charge']=float(cmpd.charge)
                info[name]['formula']=cmpd.molecular_formula
                info[name]['inchi'] =cmpd.inchi


    ###########################################################################
    #          Store & return all that nice output how they want... 
    ###########################################################################
    if out2yaml==True: 
        # Check that the desired filename is available & return the full path to write the file: 
        yaml_file=check_filename(filename=filename, default_name= 'species_info', ext='.yml', 
                                 savepath=savepath, overwrite=overwrite, return_full=True, quiet=quiet)
        # Write the dictionary to a YAML file
        with open(yaml_file, 'w') as f:   
            yaml.dump(info, f)
        print(f'Output dictioary of all Pubchem matches saved at: \n\t{yaml_file}')
    
    if return_df== True or out2excel== True: 
   
        # First convert from a record based dict to an index based on so df 
        # formatting makes sense... 
        info = _convert_to_category_dict(info) 
        
        # Convert dict to dataframe with all this information... 
        info_df = pd.DataFrame(info)
        
        if out2excel==True: 
            # Check that the filename is available & return the full path to write the file: 
            excel_file=check_filename(filename=filename, default_name= 'species_info', ext='.xlsx', 
                                savepath=savepath, overwrite=overwrite, return_full=True, quiet=quiet)
            
            # Write the dataframe to an excel file: 
            info_df.to_excel(excel_file)
            print(f'Output df of all PubChem matches saved at: \n\t{excel_file}')
    
    if return_df== True:
        return info_df, no_matches, multi_matches
    else: 
        return info, no_matches,multi_matches