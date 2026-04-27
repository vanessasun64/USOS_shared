# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 15:29:36 2022

@author: jhask
"""

import pandas as pd # ok with rdkit 
import pickle
import numpy as np 
import sys
from rdkit import Chem

from rdkit.Chem import  Descriptors, rdMolDescriptors, Fragments

def get_rdkit_frags(molec): 
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

def find_groups(df_in, use, groups, verbose=False):

    # Create an output dataframe and make a column for epoxides... 
    df_out=pd.DataFrame(); df_out[use]= df_in[use].copy()
    df_out['Epoxides']=list(np.full([len(df_out),1], np.nan)) 
    
    for  i, sp in enumerate(df_in[use]): 
        if verbose is True: print(df_out.loc[i,'MCM_Name']) # Print off name of what you're parsing. 
        molec = Chem.inchi.MolFromInchi(df_in.loc[i,'InChI'],sanitize=False, removeHs=False,logLevel=None)  #Turn this molec it into an RDKit molecule object. 
        
        if molec is not None:
            molec.UpdatePropertyCache(strict=True)  # for radicals!
            
            for key in groups:  # Loop over ever functional group you want to search for. 
                
                # Turn the SMARTs string for this functional group into a RDKit molec fragment. 
                frag = Chem.MolFromSmarts(groups[key])
                
                # Get a list of the indices of atom #s in molecule that match this fragment 
                inds=list(molec.GetSubstructMatches(frag))
                
                # Save the len of this list as the # of functional group matches you found!)
                df_out.at[i,key]=np.int64(len(inds))
                
                # Use RDKit to get fragments(not always as specific as our group matches...) 
                rd_frags=get_rdkit_frags(molec)
            
            df_out.at[i,'Epoxides'] = Fragments.fr_epoxide(molec) # Number of epoxide rings 
            
    # Add a column that has the number of OH groups on a compounds that might cause it to be an organic acid... 
    Organic_Acid_OHs=['Enols','Phenols','Thiols', 'Carboxylic_Acids']
    df_out['Organic_Acid_OHs']=df_out[Organic_Acid_OHs].sum(axis=1)
    
    return df_out    

# Dictionary of groups and their smart search strings....      
groups= dict({
    'All_Radicals': '[#8H0;v1+0]', # Includes peroxy radicals, acyl peroxy radicals, alkoxy and acyloxy radicals
    'RO2s':'[#6,S;!$([#6](=[#8]))]~[#8][#8;X1v1+0]', # RO2 peroxy radicals 
    'ROs':'[#6,S,c;!$([#6](=[#8]))]~[$([#8X2+0]),$([#8;X1v1+0]);!$([#8X2H1]);!$([#8X2][#8]);!$([#8X2][#7]);!$([#8X2]([#6])[#6]);!$([#8X2][#16]);!$([#8v+](=[#6])[#8v-])]', # All RO radicals. 
    'Acyl_RO2s':'[#6](=[#8])~[#8][#8;X1v1+0]', # acylperoxy radical
    'Acyl_ROs':'[#6](=[#8])[$([#8X2]),$([#8;X1v1+0]);!$([#8X2H1]);!$([#8X2][#8]);!$([#8X2][#7]);!$([#8X2]([#6])[#6]);!$([#8X2][#16])]',

    'All_NO3s':'[$([#7X3](=[#8X1])(=[#8X1])[#8]),$([#7X3+]([#8X1-])(=[#8X1])[#8])]', # Catches ANY NO3 group including ions, PANS, RO2NO2, RONO2,
    'Non_PAN_NO3s': '[$([#7X3](=[#8X1])(=[#8X1])[#8;!$([#8]~[#8])]),$([#7X3+]([#8X1-])(=[#8X1])[#8;!$([#8]~[#8])])]', # Non-PAN nitrates (Should count NO3 groups)
    'RONO2s': '[$([#7X3](=[#8X1])(=[#8X1])[#8]~[*;!$([#8])]),$([#7X3+]([#8X1-])(=[#8X1])[#8]~[*;!$([#8])])]', # Non-PAN Organonitrates (Same as NO3s, but requires attached R chain)
    'RO2NO2s':'[$([#7X3](=[#8X1])(=[#8X1])[#8]-[#8][#6;!$([#6]=[#8])]),$([#7X3+]([#8X1-])(=[#8X1])[#8]-[#8][#6;!$([#6]=[#8])])]', # Peroxy Nitrates
    'PANs':'[$([#7X3](=[#8X1])(=[#8X1])[#8]-[#8]-[#6](=[#8])),$([#7X3+]([#8X1-])(=[#8X1])[#8]-[#8]-[#6](=[#8]))]',  # Peroxy Acyl Nitrates
    'Nitros':'[#6][$([#7X3](=[#8])=[#8]),$([#7X3+](=[#8])[#8-])](~[#8])(~[#8])', # Nitro_groups
    'Tertiary_Nitrates': '[$([#7X3](=[#8X1])(=[#8X1])([#8]-[#6X4H0])),$([#7X3+]([#8X1-])(=[#8X1])([#8]-[#6X4H0]))]',
    
    'ALK_Chains_GT4:':'[$([CX4][CX4][CX4][CX4])]',
    'Alkenes:':'[$([CX2](=C)=C),$([CX3]=[CX3])]',
    'Alkynes':'[$([CX2]#C)]',
    'Peroxides': '[#8X2;!$([#8]-[#6](=[#8]));!$([#8]-[#7])][$([#8X2H1v2+0]),$([#8X2;!$([#8]-[#6](=[#8]));!$([#8]-[#7])])]', #ROOR or ROOH, not peracids
    'OrganicPeroxides': '[#6X4;!$([#6](=[#8]))][#8X2H0][#8X2H0][#6X4]', # R-O-O-R but not R-O-O-H            
    'HydroPeroxides':'[#8X2;!$([#8]-[#6](=[#8]));!$([#8]-[#7])][#8X2H1v2+0]', #R-O-O-H  but not C(=O)OH
    'Peracids': '[*;!$([#8]);!$([#7]);!$([#16])](=[#8])-[#8]-[#8D1H1v2+0]',    # R(=O)-O-OH  
    
    'Carbonyls': '[$([#6X3]=[#8X1]),$([#6X3+]-[#8X1-])]', # Not selective. Hits C=O
    'Carboxylic_Acids':'[#6X3;!$([#6X3H])](=[#8])[#8X2H]',    # Hits C(=O)-OH
    'Ketones': '[#6X3;!$([#6X3H1]);!$([#6X3H2]);!$([#6]-[#8])](=[#8])',    # Hits C(=O)(-R)(-R)
    'Aldehydes': '[#6X3H1;!$([#6X3H1]-[#8X2])](=[#8;!$([#8][#8])])',    # Hits R-CH(=O), but not # Hits O-CH(=O)
    'Esters':   '[$([#6H1](=[#8])(-[#8][#6])),$([#6X3](=[#8])(-[#8][#6])[*;!$([#8])]),$([#6X3](=[#8])(-[#8][#6])[#8])]', # Hits O-CH(=O), O-C(=O)(-C), or  O-C(=O)(-O) 
    'Ethers':'[#8D2]([#6;!$([#6]=[#8])])[#6;!$([#6]=[#8])]',
    'Carbonates': '[#6X3](=[#8X1])(-[$([#8X2H0]),$([#8X1-1])])(-[$([#8X2H0]),$([#8X1-1])])',
    
    'All_OHs':'[#6]-[#8X2H1]', # All C-OHs
    'Hydroxyls':'[$([#8X2H]);!$([#8X2H][#6X3](=[#8]));!$([#8X2H]-[$([#8]),$([#16](=[#8])(=[#8]))])]', # No carboxylic acids or sulfuric acids
    'Enols': '[$([#6](=[#6])(-[#8X2H1]));!$([$([$(c1ccccc1)](-[#8X2H1]));!$([c](=[#8])(-[#8X2H1]))])]', # Enols that are not phenols. 
    'Phenols':'[$([$(c1ccccc1)](-[#8X2H1]));!$([c](=[#8])(-[#8X2H1]))]', # -OH on benzene ring (but no carboxylic acids)
    'Thiols':'[#6]-[#16X2H1]',
    'Dihydroxys':'[#6X4;!$([#6]=[#8])](-[#8X2H1])(-[#8X2H1])', # C(-OH) (-OH) , not carboxylic acids not on rings. 
    'Aliphatic_Alcohols':'[$([CX4R0](-[#8X2H1]));!$([#6X3;!$([#6X3H])](=[#8])[#8X2H])!$([#6X4;!$([#6]=[#8])](-[#8X2H1])(-[#8X2H1]))]', #C-OH not on rings.
    'Other_Alcohols': ''.join(['[$([#6;'
        '!$([$([#6](=[#6])(-[#8X2H1]));!$([$([$(c1ccccc1)](-[#8X2H1]));!$([c](=[#8])(-[#8X2H1]))])]);' # Not Carbon with a Phenol
        '!$([#6]-[#16X2H1]);' # Not carbon thiol
        '!$([#6X4;!$([#6]=[#8])](-[#8X2H1])(-[#8X2H1]))];' #not a dihydroxy OH carbon. 
        '[#8X2H1]);' # But that carbon does gotta be attached to an OH group somewhere. 
        '!$([$([$(c1ccccc1)](-[#8X2H1]));!$([c](=[#8])(-[#8X2H1]))]);' #no phenol 
        '!$([#6X3;!$([#6X3H])](=[#8])[#8X2H]);' # no carb acid
        '!$([$([CX4R0](-[#8X2H1]));!$([#6X3;!$([#6X3H])](=[#8])[#8X2H])!$([#6X4;!$([#6]=[#8])](-[#8X2H1])(-[#8X2H1]))])]']), # and no aliphatic acohols. 
   
    'C':'[#6]', #  Number of Carbon atoms in molecule
    'H':'[#1]', #  Number of Hydrogen atoms in molecule
    'O':'[#8]', #  Number of Oxygen atoms in molecule
    'N':'[#7]', #  Number of Nitrogen atoms in molecule
    'S':'[#16]', #  Number of Sulfur atoms in molecule
    'Cl':'[#17]', #  Number of Chlorine atoms in molecule
    'Br':'[#35]', #  Number of Bromine atoms in molecule
    })

# Load in all the data about the MCM comopuds... 
mypath_to_pyMCM_IO_data= 'C:/Users/jhask/OneDrive/Documents/Python/'
MCM_df=pd.load_dataframe(mypath_to_pyMCM_IO_data+'v331_AllRxn_scrape_all.xlsx')

#Loop over MCM using the SMART fragment search strings above... 
out=find_groups(MCM_df, 'MCM_Name', groups)
    
