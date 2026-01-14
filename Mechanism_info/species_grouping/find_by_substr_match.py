import numpy as np  
import pandas as pd  
from rdkit import Chem
from rdkit.Chem import  Descriptors, rdMolDescriptors, Fragments

# Use each string to make an RDKit molecule object. 
molec = Chem.inchi.MolFromInchi(df.loc[ind,use], logLevel=None)

if molec is not None:
    molec.UpdatePropertyCache(strict=True)  # for radicals!
	
# Pull formula, molecular weight, its canonical smiles, and InChI strings. 
            form=rdMolDescriptors.CalcMolFormula(molec)
            mw=Descriptors.MolWt(molec)
			
 molec=Chem.rdmolops.AddHs(molec) # Make sure the molec has H's before we search for matches! 
 
 
for group, smarts_string_i in smarts_strings.items():  # Loop over ever functional group you want to search for. 
                
        # Turn the SMARTs string for this functional group into a RDKit molec fragment. 
        frag = Chem.MolFromSmarts(smarts_string_i)
        
        # Get a list of the indices of atom #s in molecule that match this fragment 
        inds=list(molec.GetSubstructMatches(frag))
        
        # Save the len of this list as the # of functional group matches you found!)
        np.float64(len(inds))
		