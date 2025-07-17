#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 16 19:10:13 2025

@author: u6044586
"""
import os
import sys
import time 
import numpy as np 
import pandas as pd 
import xarray as xr
from rdkit import Chem
from rdkit.Chem import  Descriptors, rdMolDescriptors, Fragments
from datetime import datetime, timedelta 

sys.path.append('/uufs/chpc.utah.edu/common/home/u6044586/python_scripts/modules/pyMCM/')
from pyMCM_utils import dict2df


def make_initial_var_info_spreadsheet(file, savepath): 
    ds=xr.open_dataset(file)
    
    # Get a list of all unique attrs of any variable in dataset: 
    unq_attrs=[] 
    for var in ds.data_vars: 
        for attr in list( ds[var].attrs.keys()): 
            if attr not in unq_attrs: 
                unq_attrs.append(attr)
    
    # Build intial dict holding info aon all vars in parked datasets: 
    dat={'Var':[]}
    for key in unq_attrs: dat[key]=[] 
    
    # Fill dict with info on each variable: 
    for var in ds.data_vars: 
        dat['Var'].append(var)
        for key in unq_attrs: 
            if key in list(ds[var].attrs.keys()):
                dat[key].append(ds[var].attrs[key]) 
            else: 
                dat[key].append('')
    
    # Turn dict into pandas df: 
    df=pd.DataFrame(dat) 
    
    # save dataframe as output .xlsx spreadsheet: 
    df.to_excel(savepath)
    print(f'Output var info spreadsheet saved at: \n\t{savepath}')
    
    return df 

def is_alkane(molec):

    # Iterate over all atoms in the molecule
    for atom in molec.GetAtoms():
        # Check if an atom is not carbon or hydrogen
        if atom.GetSymbol() not in ['C', 'H']:
            return False

    # Check for double or triple bonds
    for bond in molec.GetBonds():
        if bond.GetBondType() not in [Chem.BondType.SINGLE]:
            return False

    return True

def is_alkene(molecule):

    contains_double_bond = False
    
    # Iterate over all atoms in the molecule
    for atom in molecule.GetAtoms():
        # Check if an atom is not carbon or hydrogen
        if atom.GetSymbol() not in ['C', 'H']:
            return False

    # Check for at least one double bond
    for bond in molecule.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            contains_double_bond = True
        elif bond.GetBondType() not in [Chem.BondType.SINGLE, Chem.BondType.DOUBLE]:
            return False

    return contains_double_bond

def is_aromatic(molecule):
    
    # Check for at least one aromatic bond
    for bond in molecule.GetBonds():
        if bond.GetIsAromatic():
            return True

    return False

def is_halogen(molecule):
    
    # Check for the presence of halogen atoms
    halogens = {'F', 'Cl', 'Br', 'I'}
    for atom in molecule.GetAtoms():
        if atom.GetSymbol() in halogens:
            return True

    return False

def is_saturated(molecule):
    
    # Check for any double or triple bonds
    for bond in molecule.GetBonds():
        if bond.GetBondType() in (Chem.BondType.DOUBLE, Chem.BondType.TRIPLE):
            return False

    # Alternatively, check for the presence of rings, which might indicate unsaturation
    if molecule.GetRingInfo().NumRings() > 0:
        return False

    # If no multiple bonds and no rings, the molecule is saturated
    return True

def is_PAN(molecule):
    
    if molecule is None:
        return False

    # Define the SMARTS pattern for a peroxyacetyl nitrate
    # A simplified pattern capturing the essential connections in PAN
    pan_smarts = '[$([#7X3](=[#8X1])(=[#8X1])[#8]-[#8]-[#6](=[#8])),$([#7X3+]([#8X1-])(=[#8X1])[#8]-[#8]-[#6](=[#8]))]'

    # Create a molecule object from the SMARTS pattern
    pan_pattern = Chem.MolFromSmarts(pan_smarts)

    # Check if the pattern matches the molecule
    if molecule.HasSubstructMatch(pan_pattern):
        return True

    return False

def query_rdkit_info(df_in, use ,add_functional_groups:bool=False,
                        save= True, savepath: str = '', filename:str='', verbose:bool=True, 
                        nm_col:str= 'MCM_Name', overwrite:bool=False):
       """Function that takes a pandas dataframe with a column named 'InChI' or 'SMILES' and 
       uses rdkit to extract its canonical SMILES string, Formula, Molecular weight, and (optionally)
       how many functional groups each compound has. It outputs that info as a dataframe 
       and saves it to a .xlsx NOTE: Comma delimited files will NOT work
       because InChI strings contain commas.
       
       Inputs: 
       -------
           df_in - A pandas DataFrame with a column of either 'InChI', 'SMILES' or 'Canonical_SMILES' 
                   that tells us what molecules we want to query info about. Optional (If you 
                   have a column like "MCM_name" then you can set verbose to True for it to tell you which 

             
           use - column name in df_in that has INCHI string... 
           
           add_functional_groups - Option of whether you'd like to add counts of functional groups 
                      to the output dataframe. 
        
                      
           verbose - (optional) Boolean. Set True to see warnings/ errors 
           
       Outputs: 
       --------  
           df - Pandas dataframe with all the original data and new columns iwth 
                data added
       
       Author: 
       -------
           Dr. Jessica D. Haskins (jhaskins@alum.mit.edu) GitHub: @jhaskinsPhD
       
       Change Log: 
       ----------
           10/29/2021    JDH Created for pyMCM
           03/16/2025    JDH modified for USOS
           
       """
          
       if add_functional_groups is True: 
           filename='/uufs/chpc.utah.edu/common/home/u6044586/python_scripts/modules/pyMCM/data/Functional_Group_SMARTs'
           path='/uufs/chpc.utah.edu/common/home/u6044586/'
           groups= dict2df(savepath=path, full_file=filename, parse_chars=False, reverse=True)

       df= df_in.copy() # Just make a copy so you're not changing stuff in the input df.

       # Initialize lists to hold stuff and warnings of inconsistencies! 
       formulas=[]; mws=[]; smiles=[]; is_alka=[]; is_alke =[];  is_arom=[];   is_hal=[]; is_sat=[]; is_pan=[]          
       
       print ('\n', '---------- BEGIN RDKIT WARNINGS ----------')
       for iind,ind in enumerate(df.index): # Loop over all compounds in the input dataframe.    

           if df.loc[ind,use] !='':
               # Use each string to make an RDKit molecule object. 
               molec = Chem.inchi.MolFromInchi(df.loc[ind,use], logLevel=None)
               
               if molec is not None:
                   molec.UpdatePropertyCache(strict=True)  # for radicals!
                   
                   # Pull formula, molecular weight, its canonical smiles,
                   form=rdMolDescriptors.CalcMolFormula(molec)
                   mw=Descriptors.MolWt(molec)
                   alka=is_alkane(molec)
                   alke=is_alkene(molec)
                   arom=is_aromatic(molec)
                   hal=is_halogen(molec)
                   sat=is_saturated(molec)
                   pan=sat=is_PAN(molec)
                   
                   if add_functional_groups is True: 
                        molec_i=Chem.rdmolops.AddHs(molec) # Make sure the molec has H's before we search for matches! 
                        
                        for key in groups:  # Loop over ever functional group you want to search for. 
                            # Turn the SMARTs string for this functional group into a RDKit molec fragment. 
                            frag = Chem.MolFromSmarts(groups[key])
                            
                            # Get a list of the indices of atom #s in molecule that match this fragment 
                            inds=list(molec_i.GetSubstructMatches(frag))
                            
                            # Save the len of this list as the # of functional group matches you found!)
                            df.at[ind,key]=np.float64(len(inds))
                        df.at[ind,'Epoxides'] = Fragments.fr_epoxide(molec_i) # Add Number of epoxide rings 
                        df.at[ind,'Alk4s']=Fragments.fr_unbrch_alkane(molec_i) #Number of unbranched alkanes of at least 4 members (excludes halogenated alkanes)
                                      
                   sm=Chem.MolToSmiles(molec)         
               else: 
                   form= '??'; sm='??';  mw=np.nan; alka='??';  alke='??';  arom='??'; hal='??'; sat='??'; pan='??'
           else: 
               form= '??'; sm='??';  mw=np.nan; alka='??';  alke='??';  arom='??'; hal='??'; sat='??'; pan='??'
               
           #Append all the info into lists to fill output columns with later. 
           formulas.append(form)
           mws.append(mw)
           smiles.append(sm)
           is_alka.append(alka)
           is_alke.append(alke)
           is_arom.append(arom)
           is_hal.append(hal)
           is_sat.append(sat)
           is_pan.append(pan)
              
       # Add each column of info to the df if you need it... 
       df['FORMULA']=formulas 
       df['MOLECULAR_WEIGHT']=mws
       df['SMILES']=smiles
       df['IS_ALKANE']= is_alka
       df['IS_ALKENE']=is_alke
       df['IS_AROMATIC']=is_arom
       df['IS_HALOGEN']=is_hal
       df['IS_SATURATED']= is_sat
       df['IS_PAN']=is_pan
       

       print ('---------- END RDKIT WARNINGS ----------', '\n')
       
       return df     
   
def list_nc_files(directory):
    nc_files = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith('.nc'):
                nc_files.append(os.path.join(root, file))
    return nc_files

def revise_merges_n_add_better_meta(rev_spreadsheet, rev, type_mrg, time_mrg):
    
    # Read in revised spreadsheet with meta data updated / more readable & includes INCHIs of all data: 
    meta_df= pd.read_excel(rev_spreadsheet, index_col=0).fillna('')
    for col in ['IS_LUMPED','IS_VOC','IS_RONO2','IS_ATTR']:
        meta_df[col] = meta_df[col].astype(bool)
        
    # Make sure INCHI's don't have any spaces or new line chars accidentally attached to them: 
    meta_df['InChI']=[var.strip().replace(' ','').replace('/n','') for var in meta_df['InChI']]
    meta_df.to_excel(rev_spreadsheet )
    
    # Use RDKit to pull info based off the InChI codes defined in df
    meta_df=query_rdkit_info(meta_df,'InChI',add_functional_groups=False)
    
    # Get list of all raw merged .nc files: 
    data_dir='/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_merges/'
    dirpath=os.path.join(data_dir, rev, type_mrg, 'merged', time_mrg) 
    # 07072025: Vanessa added sorting because file list was out of order, messing up the dates for the combined merge
    file_list= [os.path.join(dirpath,f) for f in sorted(os.listdir(dirpath)) if f.endswith('.nc')]
    
    # Initialize list to hold all data for "master merge": 
    all_ds=[]
    
    #  Loop over each merged file in directory & add metadata: 
    for file_i in file_list: 
        
        # Load in raw merged data file: 
        ds_in=xr.open_dataset(file_i) 
        
        # Rename time as time_UTC for clarity... 
        ds_in = ds_in.rename({"time": "time_UTC"})
    
        # Create dataframe form old raw merged dataset: 
        df_in= ds_in.to_dataframe() 
        
        # Loop over all vars we have updated metadata info about and make dict holding 
        # metadata for each var to use in assigning attributes to updated dataset data variables: 
        meta={}
        
        # Define the old meta data attrs we want to keep info from raw merge and new attr names that will hold those: 
        keep_old_meta={'Original_Filename':'ORIGINAL_FILENAME','Revision':'REVISION', 
                       'Stipulations':'STIPULATIONS', 'Units':'ORIGINAL_UNITS'}
        
        # Loop over all vars in meta_df we have info about: 
        for ind in meta_df.index:
            # Make sure we're ripping data about something that should be a data var 
            # not an attribute of a data var. 
            if meta_df.loc[ind,'IS_ATTR']==False:
                # Pull out new varname and intialize meta[newvar] as a dictionary: 
                newvar=meta_df.loc[ind, 'NEW_VARNAME']
                meta[newvar]={}
                # Loop over all columns in metadata dataframe: 
                for col in meta_df.columns: 
                    # Make sure we're not pulling data on the attributes... 
                    if col not in ['IS_ATTR', 'ATTR_OF', 'ORIGINAL_FILENAME', 'NEW_VARNAME','NEW_UNITS']: 
                        # Don't assign keys/values to blank or NaN quantities in metadata dataframe: 
                        if ((type(meta_df.loc[ind, col])==str) and (len(meta_df.loc[ind, col]) >0)) or \
                            ((type(meta_df.loc[ind, col])!=str) and (np.isnan(meta_df.loc[ind, col])==0)): 
                            
                            # Assign info about this var to metadata dict: 
                            meta[newvar][col]=meta_df.loc[ind, col]
                            
                # Also add meta data from original ds that may differ with each file... 
                oldvar=meta_df.loc[ind, 'ORIGINAL_VARNAME']
                for attr_i in list(keep_old_meta.keys()):
                    if oldvar in list(ds_in.data_vars) and attr_i in list(ds_in[oldvar].attrs.keys()):
                        value=ds_in[oldvar].attrs[attr_i]
                        if attr_i=='Stipulations': value=value.replace('_', ' ') 
                        meta[newvar][keep_old_meta[attr_i]]=value
            else: 
                # Assign vars that should instead be attrbiutes to dict of what they should be attrs of... 
                var=meta_df.loc[ind, 'NEW_VARNAME']
                attr_of=meta_df.loc[ind, 'ATTR_OF']
                
                meta[attr_of][var]=ds_in[meta_df.loc[ind, 'ORIGINAL_VARNAME']]
                
        # Rename columns in dataset to new varname: 
        names= meta_df.set_index('ORIGINAL_VARNAME')['NEW_VARNAME'].to_dict()
        df_in.rename(columns=names, inplace=True)
        
        # Drop the vars that should be attributes from the df so they don't get assigned. 
        attr_vars={meta_df.loc[ind,'NEW_VARNAME']:meta_df.loc[ind,'ATTR_OF'] for ind in meta_df.index if meta_df.loc[ind, 'IS_ATTR']==True}
        df_in.drop(columns=list(attr_vars.keys()), inplace=True)
        
        # Print warning about cols we didn't re-name/ have info about in spreadsheet: 
        missing_info=[col for col in list(df_in.columns) if col not in list(names.values())]
        if len(missing_info)> 0: 
            print('Missing info about columns:'+','.join(missing_info))
        
        # Create new xarray dataset from dataframe: 
        ds_out = xr.Dataset.from_dataframe(df_in)
        
        # Assign attr info about time coordinate: 
        ds_out['time_UTC'].attrs['LONG_NAME']='Time in UTC. Offset to local timezone (MDT) is -6h in July/Aug 2024.'
        
        # Assign stuff in metadata dict as attributes of new xarray: 
        for array_name, attrs in meta.items():
            if array_name in ds_out:
                # Loop through each attribute in the attrs dictionary
                for key, value in attrs.items():
                    if type(value) == str and value != '??': 
                        ds_out[array_name].attrs[key] = value
                    elif isinstance(value, np.bool_) or isinstance(value, bool):
                        ds_out[array_name].attrs[key] = str(value)  
                    else:
                        ds_out[array_name].attrs[key] = value  
    
        # Define local time var: 
        coords=ds_out['Temp_K'].coords
        dims= ds_out['Temp_K'].dims
        ds_out['time_local']=xr.DataArray(ds_out.time_UTC.values-pd.Timedelta(hours=6), coords=coords, dims=dims)
        ds_out['time_local'].attrs['LONG_NAME']='Local Time in MDT (UTC offset =-6h in July/Aug 2024).'
        
        # Convert units to what "new units" are defined in spreadsheet: 
        for ind in meta_df.index: 
            var=meta_df.loc[ind,'NEW_VARNAME']
            new_unit=meta_df.loc[ind,'NEW_UNITS'].strip().replace(' ','')
            old_unit=meta_df.loc[ind,'ORIGINAL_UNITS'].strip().replace(' ','')
            
            if new_unit=='ppbv' and old_unit=='pptv':fac=1e-3   # 1pptv = 1e-3 ppbv 
            elif new_unit=='ppbv' and old_unit=='ppmv':fac=1e3  # 1 ppmv = 1e3 ppbv
            else: fac=1
        
            if var in list(ds_out.data_vars): 
                # Convert units using factor.
                ds_out[var].values=ds_out[var].values*fac
                # Update UNITS attribute: 
                ds_out[var].attrs['UNITS']=new_unit
            else: 
                # ALso update units of vars in attrbiutes: 
                for key, value in attr_vars.items():     
                    if var in list(ds_out[value].attrs.keys()): 
                        # Define new attr with conversion & new units
                        ds_out[value].attrs[var+'_'+new_unit]=ds_out[value].attrs[var].values*fac
                        # Delete old atribute with old units: 
                        del ds_out[value].attrs[var]

        
        # Convert temp froôm C to K (despite being named correctly already is actually in C).  
        ds_out['Temp_K'] = ds_out['Temp_K'] +273.15 
        ds_out['Temp_K'].attrs['UNITS']='degrees Kelvin (K)'
        
        # Save revised merge as its own dataset in same dir/same file name, but with 'rev_' in front of the dirpath. 
        savepath=os.path.join('/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_merges/', rev, type_mrg, 'merged', 'rev_'+time_mrg)

        # Check if the directory exists
        if not os.path.exists(savepath):
            # Create the directory if it doesn't exist: 
            os.makedirs(savepath)
        outfile= os.path.join(savepath, os.path.basename(file_i).replace(time_mrg, 'rev'+time_mrg))
        ds_out.to_netcdf(outfile)
        print(f'Revised file saved at: {outfile}')
        
        # Append dataset to list of all for "master merge". 
        all_ds.append(ds_out.isel(time_UTC=slice(0, -1)))
        
    # Concatenate all  along the 'time' dimension
    combo_ds = xr.concat(all_ds, dim='time_UTC')
    
    combo_outfile= os.path.join(savepath, 'all_'+type_mrg+'_rev'+time_mrg)+'.nc'
    combo_ds.to_netcdf(combo_outfile)
    print(f'Master merge of all saved at: {combo_outfile}')

    return combo_ds 

###############################################################################
#    Revise 30min parked USOS merges: 
###############################################################################
# Make intial var info spreadsheet for parked USOS drives: 
#March or April 2025: Jessica originally created the spreadsheet with RA data, which was missing variables
#Early June 2025: Vanessa remade this spreadsheet to add relative humidity
#Late June 2025: Vanessa had to remake the spreadsheet and redo merges, as the POPS instrument R0 data had extra variables for each aerosol size bin. 


#parked_merged_file= '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_merges/R0/CSL_MobileLab_Parked/merged/30min/USOS_30min_Merged_Parked_20240805.nc'
#parked_savepath= '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_merges/initial_USOS_parked_vars_updated_062025.xlsx'
#make_initial_var_info_spreadsheet(parked_merged_file, parked_savepath)

## JDH MADE BY-HAND EDITS TO THAT INTIAL SPREADSHEET RESULTING IN:
# "/uufs/chpc.utah.edu/common/home/haskins-group1/data/Campaign_Data/Raw_Data/USOS_2024/RA/parked/merged/revised_USOS_parked_var_info_spreadsheet.xlsx" 

# Vanessa made By-Hand edits to the initial spreadsheet, to merge Jessica's by-hand edits with new variables in (jvalues from TUV; remove iWAS):
# June 23, 2025; Vanessa deleted Particle Count variable, Added all new POPS variables. Added new iWAS variables from the new R1 data for July 15, 2024 USOS date.
# '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_merges/revised_USOS_parked_vars_updated_062025.xlsx'

# Read in revised spreadsheet with meta data updated / more readable & includes INCHIs of all data: 
parked_rev_spreadsheet='/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_merges/revised_USOS_parked_vars_updated_062025_r1_was.xlsx'

# # Call function that does the heavy lifting to revise & combine all merges: 
revise_merges_n_add_better_meta(parked_rev_spreadsheet, rev='R0', type_mrg='CSL_MobileLab_Parked', time_mrg='1hr')

  
# ###############################################################################
# #    Revise 1min Driving USOS merges: 
# # ###############################################################################
# # # Make intial var info spreadsheet for parked USOS drives: 
# # driving_merged_file= '/uufs/chpc.utah.edu/common/home/haskins-group1/data/Campaign_Data/Raw_Data/USOS_2024/RA/driving/merged/10s/USOS_10s_Merged_Drive_20240728.nc'
# # driving_savepath= '/uufs/chpc.utah.edu/common/home/haskins-group1/data/Campaign_Data/Raw_Data/USOS_2024/RA/driving/merged/initial_USOS_driving_var_info_spreadsheet.xlsx'
# # #make_initial_var_info_spreadsheet(driving_merged_file, driving_savepath)

# # ## JDH MADE BY-HAND EDITS TO THAT INTIAL SPREADSHEET RESULTING IN:
# # # "/uufs/chpc.utah.edu/common/home/haskins-group1/data/Campaign_Data/Raw_Data/USOS_2024/RA/driving/merged/revised_USOS_driving_var_info_spreadsheet.xlsx" 

# # # Read in revised spreadsheet with meta data updated / more readable & includes INCHIs of all data: 
# # driving_rev_spreadsheet='/uufs/chpc.utah.edu/common/home/haskins-group1/data/Campaign_Data/Raw_Data/USOS_2024/RA/driving/merged/revised_USOS_driving_var_info_spreadsheet.xlsx'


# # #Call function that does the heavy lifting to revise & combine all merges: 
# # revise_merges_n_add_better_meta(driving_rev_spreadsheet, rev='RA', type_mrg='driving', time_mrg='1min')
