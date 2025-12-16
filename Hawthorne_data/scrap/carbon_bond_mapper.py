import os
import sys
import time 
import numpy as np 
import pandas as pd 
import xarray as xr
from rdkit import Chem
from rdkit.Chem import  Descriptors, rdMolDescriptors, Fragments
from datetime import datetime, timedelta 

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

def is_aldehyde():
#Check for C=O double bond, that C bonded to single hydrogen, C also bonded to C unless it's formaldehyde


def is_internal_olefin

def is_ketone

def is_terminal_olefin

def is_monoterpene

def is_toluene_or_monoalkyl_aromatic

def is_xylene_or_polyalkyl_aromatic

def is_unreactive_carbon
#quaternary alkyl groups
#carboxylic acid
#ester groups
#halogenated carbons
#carbons of nitrile groups


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