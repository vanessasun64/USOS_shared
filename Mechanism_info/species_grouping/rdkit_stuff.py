#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 17 13:55:27 2025

@author: u6044586
"""

# -*- coding: utf-8 -*-
"""
Created on Wed May 26 17:07:35 2021

@author: phavala
contributors: Havala Pye, Karl Seltzer USEPA, Nash Skipper

Nash Skipper - updates for CRACMM2
    1. Standardize input SMILES string to rdkit canonical SMILES.
    2. Use canonical SMILES for explicit species mapping. This eliminates the need to account
       for different valid SMILES for the same species.
    3. Map MEK, MVK, and TOL explicitly with SMILES. Previously these were mapped based on atom
       counts and functional groups.
    4. Map STY (styrene) and EBZ (ethylbenzene) explicitly for CRACMM2.
    5. Remove dependence on kOH for mapping xylenes. All xylenes are now XYL instead of being
       split into XYE or XYM depending on kOH.
    6. Update alkene (OLI and OLT) mapping to use rdkit properties instead of being based on
       the position of the double bond in the SMILES string.
    7. Update ROC[N|P]#ALK species mapping to include mapping to ROC[N|P]#OXY# species based on 
       O:C ratio.
    8. Add V or A to the start of ROC* species to indicate gas or particle phase.
    9. Add warning if input was not mapped to a CRACMM species (mapped to UNKCRACMM).

"""
import sys 
import numpy  as np 
import warnings
import rdkit
from rdkit import Chem
from rdkit.Chem import Fragments
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem
from rdkit.Chem import GetPeriodicTable
from rdkit.Chem import BondType

#######################################################################
# Define *MUTUALLY EXCLUSIVE* SMARTS Patterns to ID & classify RO2s/ROs
#######################################################################
# Both strings below match **neutral, organic** peroxy radicals (R-O-O*) or 
# alkoxy radicals  (R-O*)attached to any type of carbon sp3, sp2, and aromatic). 
# These strings do NOT match peroxy/alkoxy radical species that contain more 
# than 1 unbound electron or those without the radical on a terminal O atom. 
# Strings are defined exclusively by including an explicity check of whther the 
# carbon atom directly bonded to the peroxy/alkoxy group is part of a carbonyl 
# group. When it is, that's an acyl RO2/RO, otherwise, it goes in other_RO2/RO 
# representingmatches to any alkyl, aryl, arylalkyl, heteroatom-substituted & 
# any other non-acyl organic RO2s/ROs. 

# Note: Using [O;X1v1+0] to make sure radicals are not charged, and that they
#		only have 1 bond. 
smarts = {
    'acyl RO2'    : Chem.MolFromSmarts('[#6](=[O])[OX2][O;X1v1+0]'),
    'non-acyl RO2': Chem.MolFromSmarts('[#6;!$([#6](=[O]))][OX2][O;X1v1+0]'),
    'acyl RO'     : Chem.MolFromSmarts('[#6](=[O])[O;X1v1+0]'),
    'non-acyl RO' : Chem.MolFromSmarts('[#6;!$([#6](=[O])[O;X1v1+0])][O;X1v1+0]'),
}

def valid_element(element_symbol):
    """
    Checks if the input is a valid element symbol on the periodic table (case-insensitive).

    Args:
        element_symbol (str): The element symbol to check.

    Returns:
        bool: True if the element symbol is valid, False otherwise.
    """
    periodic_table = GetPeriodicTable()

    try:
        if not isinstance(element_symbol, str):
            return False  # Input must be a string

        element_symbol = element_symbol.strip()  # Remove leading/trailing whitespace
        if not element_symbol:
            return False  # Handle empty strings

        # Try to get the atomic number from the symbol (case-insensitive)
        atomic_num = periodic_table.GetAtomicNumber(element_symbol.capitalize())  # Element symbols are capitalized

        return True  # Valid symbol

    except KeyError:
        return False  # Invalid symbol
    except Exception as e:
        # Handle unexpected errors (e.g., RDKit internal errors)
        print(f"An error occurred: {e}")  # Log the error if it occurs.
        return False  # Treat as invalid for safety

def count_elements(mol): 
    """Function to return the count of various defined elements & any other elements 
    within the input rdkit molecule object."""     
    
    # Define a list of the typical atoms we def are going to need count: 
    get_count_of=['C','H','O','N','S','Cl','Br','I']
    
    # And add any other unique atoms that make up this molecule to the list too. 
    get_count_of=np.unique(get_count_of+[atom.GetSymbol() for atom in mol.GetAtoms()]) 
    
    # Make sure the molec has H's before we search for the # of Hs... 
    mol=Chem.AddHs(mol)
    
    ct={} # Initalize the output dictionary 

    # Now loop over all our elements of interest & count how many there are in the molec. 
    for elem in get_count_of: 
        
        ct[f"{elem}"] = int(sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == elem))
    
    return ct

def classify_RO(mol):
    """
    Analyzes a molecule to determine if it contains an *organic* alkoxy radical
    (acyl or non-acyl) and excludes inorganic alkoxy radicals.
    """
    # Check whether the input molecule contains a match to either of these 
    # substructures: 
    has_acyl_RO= mol.HasSubstructMatch(smarts['acyl RO'])
    has_othr_RO= mol.HasSubstructMatch(smarts['non-acyl RO'])
    
    # Use results to decide if its an RO and if so, what kind of RO is is. 
    if has_acyl_RO and not has_othr_RO: 
        tf_is_RO=True; RO_type='acyl RO'
    elif has_othr_RO and not has_acyl_RO: 
        tf_is_RO=True; RO_type='non-acyl RO'
    elif has_acyl_RO and has_othr_RO: 
        smiles = Chem.MolToSmiles(mol, canonical=True)
        print(f'WARNING: Molecule for {smiles} matched as both an acyl RO & as a non-acyl RO/',
              'Issues with the SMARTS strings defining their patterns are likely to blame.',
              'The RO identification & classification might be messed up.\n')
        
        tf_is_RO=True;  RO_type=None
    else: 
        tf_is_RO=False; RO_type=None
    
    return tf_is_RO, RO_type
    
def classify_RO2(mol):
    """
    Analyzes a molecule to determine if it contains an *organic* peroxy radical
    (acyl or non-acyl) and excludes inorganic peroxy radicals.
    """
    # Update the property cache  to ensure that RDKit's internal representation of
    # the molecule is accurate / reccomended before SMARTS matching. 
    mol.UpdatePropertyCache(strict=True)  
    
    # Check whether the input molecule contains a match to either of these 
    # substructures: 
    has_acyl_RO2= mol.HasSubstructMatch(smarts['acyl RO2'])
    has_othr_RO2= mol.HasSubstructMatch(smarts['non-acyl RO2'])
    
    # Use results to decide if its an RO2 and if so, what kind of RO2 is is. 
    if has_acyl_RO2 and not has_othr_RO2: 
        tf_is_RO2=True; RO2_type='acyl RO2'
    elif has_othr_RO2 and not has_acyl_RO2: 
        tf_is_RO2=True; RO2_type='non-acyl RO2'
    elif has_acyl_RO2 and has_othr_RO2: 
        smiles = Chem.MolToSmiles(mol, canonical=True)
        print(f'WARNING: Molecule for {smiles} matched as both an acyl RO2 & as a non-acyl RO2.',
              'Issues with the SMARTS strings defining their patterns are likely to blame.',
              'The RO2 identification & classification might be messed up.\n') 
        tf_is_RO2=True;  RO2_type=None
    else: 
        tf_is_RO2=False; RO2_type=None
    
    return tf_is_RO2, RO2_type



class atmos_mol:
    def __init__(self, mol):
        """
        Initializes an atmos_mol object from an existing RDKit Mol object.

        Args:
            mol: An RDKit Mol object.  If None, an error is raised.
        """
        if not isinstance(mol, Chem.Mol):
            raise TypeError("Input must be an RDKit Mol object.")

        self.mol = mol
        
        # Add a list of all the elements inside this molecule: 
        ls=[]
        [ls.append(str(atom.GetSymbol())) for atom in self.GetAtoms() if str(atom.GetSymbol()) not in ls]
        self.made_of=ls
        
        # Add counts of all the major atoms & the ones in this molecule: 
        self.atom_counts=count_elements(mol)
                                    
        
    def __getattr__(self, name):
        """
        Forwards attribute/method lookups from the underlying RDKit Mol object.
        """
        # Whenever you try to access an attribute or method of atmos_mol that 
        # doesn't already exist directly on the atmos_mol class, it checks if the 
        # underlying RDKit Mol object DOES have that method & returns it instead.
        if hasattr(self.mol, name):
            return getattr(self.mol, name)
        else:
            raise AttributeError(f"'atmos_mol' object has no attribute '{name}'")

    
    def get_count_of(self,elem):
        """ Return the count of any element in the molecule (e.g. carbon number ,etc).""" 
        if self is None: return None  # Handle None input gracefully
        
         # Make sure they have passed a valid element though... 
        if valid_element(elem) is True: 
            
            # Must actually add Hydrogens to the molecule befor we can count them. 
            mol2use=Chem.AddHs(self) if elem =='H' else self 

            # Update the property cache to ensure that RDKit's internal 
            # representation of the molecule is accurate 
            mol2use.UpdatePropertyCache(strict=True)  
            
            ct= int(sum(1 for atom in mol2use.GetAtoms() if atom.GetSymbol() == elem))
            return ct 

        else: 
            return None 
        
    def get_OtoC(self): 
        """Get O to C ratio of molecule. Returns None if # of C's =0""" 
        
        if self.atom_counts['C'] >0: 
            o2c=self.atom_counts['O']/self.atom_counts['C']
        else:
            o2c=None
        return o2c 
    
    def is_organic(self):
        return self.atom_counts['C'] >0 
   
    def is_RO2(self):
        tf_is_RO2, _ = classify_RO2(self)
        return tf_is_RO2
    
    def get_RO2_type(self):
        _, RO2_type = classify_RO2(self)
        return RO2_type
    
    def is_RO(self):
        tf_is_RO, _ = classify_RO(self)
        return tf_is_RO
    
    def get_RO_type(self):
        _, RO_type = classify_RO(self)
        return RO_type
        
     
    def count_CdblC(self, verbose:bool=False):
        """
        Counts the number of C=C double bonds in an RDKit molecule,
        including those RDKit stores as "double" or "aromatic... but ensuring each 
        bond is counted only once.
        
        Args:
            mol (rdkit.Chem.Mol): The RDKit molecule to analyze.
        Returns:
            int: The number of C=C double bonds in the molecule.
        """
        if self is None:
            return False  # Handle None input gracefully
        
        counted_bonds = set()
        n_CdblC = 0; n_CdblC_arom=0
        for bond in self.GetBonds():
            
            # Query RDKit to get the bond type & the two atoms its binding: 
            btype=bond.GetBondType()
            atom1=bond.GetBeginAtom().GetSymbol()
            atom2=bond.GetEndAtom().GetSymbol()
            atom1_idx = bond.GetBeginAtomIdx()
            atom2_idx = bond.GetEndAtomIdx()
            
            # Use sorted tuple to ensure (1, 2) and (2, 1) are treated the same
            bond_tuple = tuple(sorted((atom1_idx, atom2_idx)))
            
            # Consider double AND aromatic bonds between two C atoms, but make sure 
            # not to double count. 
            if ((btype == BondType.DOUBLE)  and atom1 == 'C' and atom2 == 'C'):
               
               # Figure out if you already counted the bond between these atoms in 
               # your count or not: 
               if bond_tuple not in counted_bonds:
                   n_CdblC += 1
                   counted_bonds.add(bond_tuple)
            elif (btype == BondType.AROMATIC and atom1 == 'C' and atom2 == 'C'):
                # Figure out if you already counted the bond between these atoms
                if bond_tuple not in counted_bonds:
                    n_CdblC_arom += 1
                    counted_bonds.add(bond_tuple)
        
        # Total # of double bondsin the molecule is the sum of what RDKit calls doubles 
        # and 1/2 the # it specifies as aromatic bonds... 
        total_dbls= n_CdblC+n_CdblC_arom/2
        if verbose is True: 
            print(f'Dbl:{n_CdblC}  +  Arom:{n_CdblC_arom/2} = Total Dbls:{total_dbls}')
        
        return total_dbls

    def has_radical(self, verbose:bool=False):
        """Highly generic check for if any atom in molecule has an unpaired electon"""
        if self is None:
            return False  # Handle None input gracefully
        
        # Update the property cache  to ensure that RDKit's internal representation of
        # the molecule is accurate / reccomended before SMARTS matching. 
        self.UpdatePropertyCache(strict=True)  
        
        rad=False
        
        # Check if the molecule has any radicals
        for atom in self.GetAtoms():
            num_radical_e = atom.GetNumRadicalElectrons()
            rad=True if num_radical_e > 0 else rad 
            
            if verbose is True: 
                print(f"Atom {atom.GetSymbol()} is a radical with {num_radical_e} unpaired electrons.")
        return rad       
    
    def is_alkane(self):
        """
        Checks if an RDKit molecule is an alkane with a number of carbon atoms
        within the specified range (inclusive).

        Args:
            mol: An RDKit molecule object.
            min_carbons (int): The minimum number of carbon atoms required (default is 4).
            max_carbons (int): The maximum number of carbon atoms allowed (default is infinity).

        Returns:
            bool: True if the molecule is an alkane with a number of carbons in the range,
                  False otherwise.
        """

        if self is None:
            return False  # Handle None input gracefully

        # Check for only carbon and hydrogen atoms
        if not all([item  in['C','H'] for item in self.made_of]): 
            return False

        # Check for single bonds only
        for bond in self.GetBonds():
            if bond.GetBondType() != Chem.BondType.SINGLE:
                return False

        # Check for fully saturated carbons
        for atom in self.GetAtoms():
            if atom.GetAtomicNum() == 6:
                num_Hs = atom.GetTotalNumHs()
                if num_Hs + atom.GetDegree() != 4:  # A carbon's degree + number of Hs must equal 4.
                    return False

        return True  # All conditions met, it's an alkane!
    
    def count_alcohols(self):
        """Count # of alphatic & aromatic alcohols (RDKit's method)"""
        nalcohol  = Chem.Fragments.fr_Al_OH(self,countUnique=True) + \
                    Chem.Fragments.fr_Ar_OH(self,countUnique=True)    
        return nalcohol
                      
                      
    def is_alkene(self):
        """
        Checks if an RDKit molecule is an alkene (contains at least one C=C double bond)
        but is not aromatic (does not contain any aromatic rings).
    
        Args:
            mol: An RDKit molecule object.
    
        Returns:
            bool: True if the molecule is an alkene and not aromatic, False otherwise.
        """
        if self is None:
            return False  # Handle None input gracefully
    
        # Check if it has ANY C=C double bonds (including aromatic ones). 
        if self.count_CdblC()==0 :
          return False # No double bonds, so can't be an alkene.
    
        # Check for aromaticity (if there are any aromatic bonds, it is aromatic)
        for bond in self.GetBonds():
            if bond.GetIsAromatic():
                return False  # It is aromatic, so it can't be the result.
            
        # Ok, now also check that it doesn't have a smaller ring (cyclo, etc). 
        if self.has_ring():
            return False 

        return True  # It's an alkene and not aromatic


    def has_ring(self): 
        got_ring=None
        for atom in range(len(self.GetAtoms())):
            if self.GetAtomWithIdx(atom).IsInRing():
                got_ring = True
                break
            else: 
                got_ring = False 
        return got_ring
    
    def is_TMB(self):
        """
        Checks if the input molecule is a trimethylbenzene (a benzene ring with 
        exactly three methyl groups) or not. Returns BOOL
        """
    
        if self is None:
            return False  # Handle None input gracefully
    
        # Check for a benzene ring
        benzene = Chem.MolFromSmiles('c1ccccc1')
        if not self.HasSubstructMatch(benzene):
            return False
    
        # Find the attachment points on the benzene ring
        matches = self.GetSubstructMatches(benzene)
    
        if len(matches) != 1:
            return False # Not one benzene ring.
    
        # Count the number of methyl groups attached to the benzene ring
        methyl_count = 0
        benzene_atoms = matches[0]
        for atom_index in benzene_atoms:
            atom = self.GetAtomWithIdx(atom_index)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor.GetTotalNumHs() == 3:  # Methyl group
                    methyl_count += 1
    
        # It's a trimethylbenzene if there are exactly three methyl groups attached
        return methyl_count == 3
     


                 

    
    
                                                                    
def count_functional_groups(mol): 
    # Count functional groups (http://rdkit.org/docs/source/rdkit.Chem.Fragments.html)
    n_acid     = rdkit.Chem.Fragments.fr_COO(mol,countUnique=True)     # carboxylic acid
    n_ketone   = rdkit.Chem.Fragments.fr_ketone(mol,countUnique=True)
    n_aldehyde = rdkit.Chem.Fragments.fr_aldehyde(mol,countUnique=True)
    n_carbonyl = n_ketone + n_aldehyde
    n_benzene  = rdkit.Chem.Fragments.fr_benzene(mol,countUnique=True)
    n_alcohol  = rdkit.Chem.Fragments.fr_Al_OH(mol,countUnique=True) + \
                  rdkit.Chem.Fragments.fr_Ar_OH(mol,countUnique=True)# aliphatic and aromatic
    n_furan    = rdkit.Chem.Fragments.fr_furan(mol,countUnique=True) # number of furan rings

    n_nitrate = smiles_upper.count('O[N+](=O)[O-]') + smiles_upper.count('O[N+]([O-])=O')
    n_peroxide = smiles_upper.count('COO') +  smiles_upper.count('OOC') - smiles_upper.count('COOC')
    #namine = smiles_upper.count('CN') + smiles_upper.count('NC') + smiles_upper.count('C(N') + smiles_upper.count('N(C')
    
#######################################################################
# CODE FOR DEBUGGING THE ABOVE STUFF TO MAKE SURE IT WORKS RIGHT: 
#######################################################################
def compare_dicts(smiles, dict1, dict2, print_me):
    """
    Compares two dictionaries and prints the differences.
    
    Args:
        dict1 (tuple): A tuple containing the name of the first dictionary and the actual dictionary.
        dict2 (tuple): A tuple containing the name of the second dictionary and the actual dictionary.
    """
    
    dict1_name, dict1_data = dict1
    dict2_name, dict2_data = dict2
    
    # Get the keys that are in both dictionaries
    common_keys = set(dict1_data.keys()) & set(dict2_data.keys())
    
    # Get the keys that are only in dict1
    dict1_only_keys = set(dict1_data.keys()) - common_keys
    
    # Get the keys that are only in dict2
    dict2_only_keys = set(dict2_data.keys()) - common_keys
    
    # Print the differences
    pass_test=True
    if dict1_only_keys:
        if pass_test is True: 
            print(print_me)
            pass_test=False
        these1='  ,  '.join([f"'{key}': {dict1_data[key]}" for key in dict1_only_keys])
        print(f"\tFAIL:The following keys are found only in {dict1_name}... {these1}")


    if dict2_only_keys:
        if pass_test is True: 
            print(print_me)
            pass_test=False
        these2='  ,  '.join([f"{key}: {dict2_data[key]}" for key in dict2_only_keys])
        print(f"\tFAIL:The following keys are found only in {dict2_name}... {these2}")
        pass_test=False

    for key in common_keys:
        if dict1_data[key] != dict2_data[key]:
            if pass_test is True: 
                print(print_me)
                pass_test=False
                
            print(f"\tFAIL: Different values for key '{key}' found in dicts:")
            print(f"\t\t['{key}']={dict1_data[key]} in {dict1_name}")
            print(f"\t\t['{key}']={dict2_data[key]} in {dict2_name}")
    
    return pass_test
      
def check_radical_matches(smiles, rtype):
    mol = Chem.MolFromSmiles(smiles)
    matches1 = mol.GetSubstructMatches(smarts['all RO2'])
   
    print(f"\t# of any RO2 SMARTS pattern matches: {len(matches1)}")
    for match in matches1:
        print(f"\t\t on Atoms: {match}")
    
    matches2 = mol.GetSubstructMatches(smarts['acyl '+rtype])
    print(f"\t# of acyl {rtype} SMARTS pattern matches: {len(matches2)}")
    for match in matches2:
        print(f"\t\t on Atoms: {match}")
    
    matches3 = mol.GetSubstructMatches(smarts['non-acyl '+rtype])
    print(f"\t# of non-acyl {rtype} SMARTS pattern matches: {len(matches3)}")
    for match in matches3:
        print(f"\t\t on Atoms: {match}")
    print('\n')
    return 

def test_radical_id_n_classification(): 
    # Define some test  cases w/ real answers to compare to function results: 
    right4rads = {'Is_RO2': False, 'RO2_type': None, 
                  'Is_RO': False, 'RO_type': None, 'Is_Rad':True}
    right4nrads={**right4rads, 'Is_Rad':False} 

    # Define dict with input SMILES strings as keys.. Expected answers as values. 
    test_cases = {
        # Some dif RO2s where we should indeed get a match: 
        "[O]OC(=O)C=C(C)CON(=O)=O":   {**right4rads,'Is_RO2': True,'RO2_type':'acyl RO2',    'reason': "It's an Acyl RO2."},
        "CC(O[N+](=O)[O-])(C(=O)CO[N+](=O)[O-])C(=O)O[O]": 
                                      {**right4rads,'Is_RO2': True,'RO2_type':'acyl RO2','reason': "It's an Acyl RO2(MCM name=MMALNACO3)."},
        "CC(O[N+](=O)[O-])(C(=O)CC(=O)C=O)C(=O)O[O]": 
                                      {**right4rads,'Is_RO2': True,'RO2_type':'acyl RO2','reason': "It's an Acyl RO2(MCM name=NC61CO3)."},
        "CC1(OC1C(=O)C)C(=O)O[O]":    {**right4rads,'Is_RO2': True,'RO2_type':'acyl RO2','reason': "It's an aryl-Acyl RO2(MCM name=EPXEKTCO3)."},
        "c1ccccc1O[O]":               {**right4rads,'Is_RO2': True,'RO2_type':'non-acyl RO2','reason': "It's an Aryl RO2."},
        "CO[O]":                      {**right4rads,'Is_RO2': True,'RO2_type':'non-acyl RO2','reason': "It's an Alkyl RO2."},
        "C1CCCCC1O[O]":               {**right4rads,'Is_RO2': True,'RO2_type':'non-acyl RO2','reason': "It's a Cycloalkyl RO2."},
        "[O]OCC=C(C)CON(=O)=O":       {**right4rads,'Is_RO2': True,'RO2_type':'non-acyl RO2','reason': "It's an Alkyl RO2 from isoprene+NO3."},
        # Some dif ROs where we should get a match: 
        "CC(=O)C=C(C)C([O])=O":        {**right4rads,'Is_RO': True,'RO_type':'acyl RO',     'reason': "It's an Acyl RO (MCM name=C4MCODBCO2)."},
        "[O]C(=O)C=C(C=O)[N+](=O)[O-]":{**right4rads,'Is_RO': True,'RO_type':'acyl RO',     'reason': "It's an Acyl RO (MCM name=NC4DCO2)."},
        "CC1(C)CC(C(=O)CCC=O)C1CCC([O])=O": {**right4rads,'Is_RO': True,'RO_type':'acyl RO', 'reason': "It's a aryl-acyl RO (MCM name=C1210CO2)."},
        "C(C)(C)(C)[O]":              {**right4rads,'Is_RO': True,'RO_type':'non-acyl RO', 'reason': "It's a non-acyl, sterically hindered RO (tert-Butoxy)"},
        "[O]CC=C(C)CON(=O)=O":         {**right4rads,'Is_RO': True,'RO_type':'non-acyl RO', 'reason': "It's a non-acyl RO from isoprene+NO3."},
        "c1ccccc1[O]":                 {**right4rads,'Is_RO': True,'RO_type':'non-acyl RO', 'reason': "It's an Aryl RO."},
        "C[O]":                        {**right4rads,'Is_RO': True,'RO_type':'non-acyl RO', 'reason': "It's an Alkyl RO."},
        # Test it against some radicals it SHOULD fail. 
        "[N+](=O)O[O-]": {**right4rads, 'reason': "It's an *INORGANIC* peroxy nitrate (has no carbon)"},
        "[O]Cl": {**right4rads, 'reason': "It's an *INORGANIC* radical (has no carbon)"},
        "SO[O]": {**right4rads, 'reason': "It's an *INORGANIC* RO (has no carbon)"},
        "CO[O-]": {**right4nrads, 'reason': "It's a *CHARGED* RO2 (anion)"},
        "C[O+]": {**right4rads, 'reason': "It's a *CHARGED* RO (cation)"},
        "CC(=O)[O+]": {**right4rads, 'reason': "It's a *CHARGED* acyl RO"},
        "[O][O]": {**right4rads, 'reason': "It's a diradical with an RO2 group, but has no carbon."},
        "[Cl]": {**right4rads, 'reason': "It's a radical, but NOT a RO/RO2 since it's *INORGANIC* (has no carbon)"},
        # Test it against some non-radicals w/ similar underlying functional groups: 
        "OC(=O)CC1CC(C(=O)C)C1(C)C":    {**right4nrads, 'reason': "its NOT a radical. It's a pinonic acid."},
        "[Cl-]":{**right4nrads, 'reason': "It's NOT a radical & is *INORGANIC* (has no carbon)"},
        'O=S=O':  {**right4nrads,'reason': "its not a radical & is *INORGANIC* (has no carbon)"},
        "CC[O][O]CC": {**right4nrads, 'reason': "It's NOT a radical. It's a (stable) aliphatic organic peroxide"},
        "c1ccccc1[O][O]c1ccccc1": {**right4nrads, 'reason': "It's NOT a radical. It's a (stable) aromatic organic peroxide"},
        "CC(=O)OOC(=O)CC": {**right4nrads, 'reason': "It's NOT a radical. It's a (stable) acyl peroxide"},
        "CCO": {**right4nrads, 'reason': "It's NOT a radical. It's a (stable) aliphatic alcohol (aliphatic)"},
        "c1ccccc1O": {**right4nrads, 'reason': "It's NOT a radical. It's a (stable) Phenol (aromatic)"},
        "CC(=O)OC": {**right4nrads, 'reason': "It's NOT a radical. It's a (stable) Ester."},
        "CC(=O)O": {**right4nrads, 'reason': "It's NOT a radical. It's a (stable) Carboxylic acid."},
    }

    # Test the is_RO2, is_RO and is_radical functions: 
    failed=False
    for smiles, expected in test_cases.items():
        
        # Convert smiles to a rdkit molecule. 
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Warning: Invalid SMILES: {smiles}")
            sys.exit()
            continue
        
        # Run functions: 
        #tf_is_rad=is_radical(mol)
        tf_is_RO2, RO2_type = classify_RO2(mol)
        tf_is_RO, RO_type = classify_RO(mol)
        
        # Print off messages noting if stuff failed for Is_RO2(): 
        if tf_is_RO2!=expected['Is_RO2']: 
            failed=True
            print(f"is_RO2() FAILED for SMILES='{smiles}' on BOOL identification:\n",
                  f"\tInfo on compound: {expected['reason']}\n",
                  f"\tGot : {tf_is_RO2}\n",
                  f"\tExpected :{expected['Is_RO2']}\n\n")
        if ((expected['Is_RO2'] and RO2_type != expected['RO2_type']) or 
            (not expected['Is_RO2'] and RO2_type is not None)):
            failed=True
            print(f"is_RO2() FAILED for SMILES='{smiles}' on TYPE classification:\n",
                  f"\tInfo on compound: {expected['reason']}\n",
                  f"\tGot : {RO2_type}\n",
                  f"\tExpected :{expected['RO2_type']}")
            check_radical_matches(smiles, 'RO2')

        # Print off messages noting if stuff failed for Is_RO(): 
        if tf_is_RO2!=expected['Is_RO2']: 
            failed=True
            print(f"is_RO() FAILED for SMILES='{smiles}' on BOOL identification:\n",
                  f"\tInfo on compound: {expected['reason']}\n",
                  f"\tGot : {tf_is_RO}\n",
                  f"\tExpected :{expected['Is_RO']}\n\n")
        if ((expected['Is_RO'] and RO_type != expected['RO_type']) or 
            (not expected['Is_RO'] and RO_type is not None)):
            failed=True
            print(f"is_RO() FAILED for SMILES='{smiles}' on TYPE classification:\n",
                  f"\tInfo on compound: {expected['reason']}\n",
                  f"\tGot : {RO_type}\n",
                  f"\tExpected :{expected['RO_type']}")
            check_radical_matches(smiles, 'RO')
                
        # # Print off messages noting if stuff failed for Is_Radical(): 
        # if tf_is_rad!=expected['Is_Rad']: 
        #     failed=True
        #     print(f"is_Radical() FAILED for SMILES='{smiles}' on BOOL identification:\n",
        #           f"\tInfo on compound: {expected['reason']}\n",
        #           f"\tGot : {tf_is_rad}\n",
        #           f"\tExpected :{expected['Is_Rad']}\n\n")
    
    if failed ==False: 
        print('ALL TESTS PASSED for is_RO2(), is_RO(), is_radical()') 
    else: 
        print('Tests completed for is_RO2(), is_RO(), is_radical().')

    return 

def test_count_elements(): 
    
    # Define some smiles to test & the actual answers: 
    right_ans={elem:0 for elem in ['C','H','O','N','S','Cl','Br','I']}
    test_cases ={
           'CC(=C)C=C'      : {**right_ans, 'C':5, 'H':8},              # Isoprene (C5H8)
           'CO[O]'          : {**right_ans, 'C':1, 'H':3, 'O':2},       # MO2 (CH3O2) 
           'C/C(=C/C=O)/COO': {**right_ans, 'C':5, 'H':8, 'O':3,'R':2}, # HPALD 1 (C5H8O3)
           'Br[Hg]Br'       : {**right_ans, 'Br':2 ,'Hg':2, },          # HgBr2  
           'CC(O)(CO[N+](=O)[O-])C(O[N+](=O)[O-])C(=O)O[O]':            # INCNCO3  (C5H7N2O10)
                              {**right_ans, 'C':5,'H':7,'N':2,'O':10}, 
           'CC(=O)CCC1C(C(CO)(CCC=O)O[N+](=O)[O-])CC1(C)C':             # BCALNO3 (C15H25NO6)
                              {**right_ans, 'C':15,'H':25,'N':1,'O':6}, 
           'ClO[O]'         : {**right_ans, 'Cl':1,'O':2}               # CLOO*  (ClO2)
           }
        
    # Initialize dir to hold test results: 
    results=[] 
    
    for smiles,expected in test_cases.items(): 
        # Convert smiles to a rdkit molecule. 
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Warning: Invalid SMILES: {smiles}")
            sys.exit()    
            
        # Test function to count number of atoms: 
        atom_nums=count_elements(mol);  
        
        # Compare to results we expected to get... print info on diffs. 
        print_me=f'Testing count_elements() using SMILES="{smiles}"'
        pass_test=compare_dicts( ('count_elements()',atom_nums),('Expected Answer', expected), print_me)
        
        # Store info if we passed in dict: 
        results.append(pass_test) 
        
        if all(results): 
            print('ALL TESTS PASSED! for count_elements()')
        else: 
            print('Tests completed for count_elements()')
        
    return
    
def test_count_CdblC():
    """
    Tests the count_CdblC() function with several test cases.
    """
    # Specify test cases & the right answer to # of double bonds they contain: 
    test_cases = {
        # Test some common atmospheric compounds first: 
        'C=CC(=C)C': 2,                  # isoprene 
        'Cc1ccccc1': 3,                  # Toluene
        'C=Cc1ccccc1': 4,                # styrene
        'O1CC1': 0,                      # ethylene oxide (an epoxide)
        'C=Cc1ccccc1':4 ,                # styrene (added in CRACMM2)
        'C1=CC=CC=C1': 3,                # benzene
        'CCc1ccccc1': 3 ,                # ethylbenzene (added in CRACMM2)
        'C=CC(C)=O' :  1 ,               # methly vinyl ketone
        'C=CC=C' :  2,                   # 1,3 butadiene   
        'C=CC=O' :  1 ,                  # acrolein
        'C=C': 1,                        # ethene aka ethylene
        'C#C':0 ,                        # acetylene
        # And now just go ham. 
        'CC1(C)C2CC(OO)C(C)(O)C1C2': 0, # APINBOOH 
        'O=CCCC(OO)(CO[N+](=O)[O-])C1CC(C)(C)C1CCC(=O)C': 0,  # NBCALOOH
        'OC(=O)c1ccccc1': 3,            # PHCOOH
        'C1CC=CC1': 1,                  # cyclopentene
        'C=CC=C': 2,
        'C=C(C)C=C': 2,
        'C=CC#N': 1,
        'C1CC1': 0 ,
        'C=CC(=O)C': 1,
        'C=CC=CC=C': 3,
        'C=C[N+](=O)[O-]': 1,
        'C=CC(Cl)=C': 2,
        'C=CC=CC=CC=C[O]': 4,
        'C=C1C=CC=C1': 3,
    }

    for smiles, expected_count in test_cases.items():
        mol = atmos_mol(Chem.MolFromSmiles(smiles))
        actual_count = mol.count_CdblC()
        assert actual_count == expected_count, f"Test case failed for SMILES '{smiles}': expected {expected_count}, got {actual_count}"

    print("All test cases passed for count_CdblC()!")
    return


def test_is_trimethylbenzene():

    # Define some SMILES strings to test and the expected answers (True/False)
    test_cases = {
        'Cc1cc(C)cc(C)c1': True,  # 1,3,5-trimethylbenzene
        'Cc1ccccc1': False,       # Toluene (one methyl)
        'c1ccccc1': False,         # Benzene (no methyls)
        'Cc1cc(C)c(C)cc1': True,  # 1,2,4-trimethylbenzene
        'Cc1c(C)c(C)ccc1': True,  # 1,2,3-trimethylbenzene
        'CC1=CC=CC=C1': False,    # Toluene (Kekulé form)
        'CC1=C(C)C=C(C)C=C1': True, # 1,3,5-trimethylbenzene (Kekulé form)
        'C=C': False,           # Ethylene (alkene)
        'C1CCCCC1': False,      # Cyclohexane
        'CC(C)C': False        # Isobutane
    }

    # Initialize list to hold test results:
    results = []

    for smiles, expected in test_cases.items():
        # Convert SMILES to an RDKit molecule
        mol = atmos_mol(Chem.MolFromSmiles(smiles))
        if mol is None:
            print(f"Warning: Invalid SMILES: {smiles}")
            sys.exit()

        # Test the is_trimethylbenzene function
        actual = mol.is_trimethylbenzene()

        # Compare the actual result to the expected result
        print_me = f'Testing is_trimethylbenzene() using SMILES="{smiles}"'
        pass_test = (actual == expected)
        print(f"{print_me}: {'PASSED' if pass_test else 'FAILED'} (Expected: {expected}, Actual: {actual})") 


        # Store the test result
        results.append(pass_test)

    if all(results):
        print('ALL TESTS PASSED! for is_TMB()')
    else:
        print('Tests completed for is_TMB() - Some FAILED :(')
        
        
def debug_all(): 
    # Test function count_elements()
    test_count_elements()
    
    # Test functions is_RO2(), is_RO(), and is_radical for ID'ing & classifying type. 
    test_radical_id_n_classification() 
    
    # Test our abaility to count the  # of double bonds in a rdkit molecule: 
    test_count_CdblC() 
        
    return 
    
#######################################################################
#                                     FROM CRACMM
#######################################################################

def count_substruct_matches(mol, pattern):
    """Utility Function to count occurrences of a SMARTS pattern in a molecule"""
    if pattern is None:
        print(f"Invalid SMARTS pattern: {pattern}")
        return 0
    return len(mol.GetSubstructMatches(pattern))



def standardized_SMILES(smiles_input): 
    if smiles_input == '-':
        unksmiles_msg = (
            f'SMILES {smiles_input} not recognized. Mapping to UNKSMILES.'
        )
        warnings.warn(unksmiles_msg)
        return 'UNKSMILES'
    smiles       = Chem.CanonSmiles(smiles_input) # standardize input SMILES string
    m            = Chem.MolFromSmiles(smiles)
    smiles_upper = smiles_input.upper()
    return m

def get_cracmm_roc(smiles_input,koh,log10cstar,phase=None):
    '''
    This function is part of CRACMM
    
    :param smiles_input: smiles string
    :type smiles_input: string
    :param koh: kOH
    :type koh: float
    :param log10cstar:
    :type log10cstar: float
    :param phase: This is used to label species that can exist in both gas and particle phases.
    :type phase: Default is None, other options are "gas" and "particle"
    :returns: mechspecies (string) 
    '''
    
    # CRACMM2

    # Prep inputs
    if smiles_input == '-':
        unksmiles_msg = (
            f'SMILES {smiles_input} not recognized. Mapping to UNKSMILES.'
        )
        warnings.warn(unksmiles_msg)
        return 'UNKSMILES'
    smiles       = Chem.CanonSmiles(smiles_input) # standardize input SMILES string
    m            = Chem.MolFromSmiles(smiles)
    smiles_upper = smiles_input.upper()
    
    # Count C=C and atoms
    nCdblC  = smiles_upper.count('=C')-smiles_upper.count('O=C')
    if nCdblC < 0:
        nCdblC = 0
    nC      = smiles_upper.count('C')-smiles_upper.count('CL')
    nO      = smiles_upper.count('O')
    nN      = smiles_upper.count('N')
    nSi     = smiles_upper.count('SI')
    nH      = 0
    for atom in m.GetAtoms():
        nH += atom.GetTotalNumHs()
    # O:C ratio
    if nC > 0:
        OtoC = nO/nC
    
    # Count functional groups (http://rdkit.org/docs/source/rdkit.Chem.Fragments.html)
    nacid     = rdkit.Chem.Fragments.fr_COO(m,countUnique=True)     # carboxylic acid
    nketone   = rdkit.Chem.Fragments.fr_ketone(m,countUnique=True)
    naldehyde = rdkit.Chem.Fragments.fr_aldehyde(m,countUnique=True)
    ncarbonyl = nketone + naldehyde
    nbenzene  = rdkit.Chem.Fragments.fr_benzene(m,countUnique=True)
    nalcohol  = rdkit.Chem.Fragments.fr_Al_OH(m,countUnique=True) + \
                  rdkit.Chem.Fragments.fr_Ar_OH(m,countUnique=True)      # aliphatic and aromatic
    nfuran    = rdkit.Chem.Fragments.fr_furan(m,countUnique=True) # number of furan rings
    # gotring variable is never used so this could be removed
    #     it may be useful to keep it here commented out in case it is needed in the future
    #for atom in range(len(m.GetAtoms())):
        #if m.GetAtomWithIdx(atom).IsInRing():
        #    gotring = 1 # 0 = no ring, 1 = ring
        #    break
        #else: gotring = 0
    nnitrate = smiles_upper.count('O[N+](=O)[O-]') + smiles_upper.count('O[N+]([O-])=O')
    nperoxide = smiles_upper.count('COO') +  smiles_upper.count('OOC') - smiles_upper.count('COOC')
    #namine = smiles_upper.count('CN') + smiles_upper.count('NC') + smiles_upper.count('C(N') + smiles_upper.count('N(C')
    tfmonoterpene = (nC == 10 and nH == 18 and nO == 1) or (nC == 10 and nH == 16) 

    # Mapper is for ROC only and not elemental carbon
    if   ( nC <= 0 ):                 mechspecies = 'UNKCRACMM'
    elif ( smiles == '[C]' ):         mechspecies = 'UNKCRACMM'
    
    # Map CO to UNKCRACMM; CO will be mapped to SLOWROC if not handled explicitly
    elif ( smiles == '[C-]#[O+]' ):   mechspecies = 'UNKCRACMM'
    # The same applies to CO2
    elif ( smiles == 'O=C=O' ):       mechspecies = 'UNKCRACMM'

     # Explicit species
    elif ( smiles == 'CC=O' ):        mechspecies = 'ACD'   # acetaldehyde
    elif ( smiles == 'C#C' ):         mechspecies = 'ACE'   # acetylene
    elif ( smiles == 'CC(C)=O' ):     mechspecies = 'ACT'   # acetone
    elif ( nC==6 and nH==6 and nO==0 and nbenzene==1 ):
                                      mechspecies = 'BEN'   # benzene
    elif ( smiles == 'C'  ):          mechspecies = 'ECH4'  # methane
    elif ( smiles == 'CCO'):          mechspecies = 'EOH'   # ethanol
    elif ( smiles == 'C=C'):          mechspecies = 'ETE'   # ethene aka ethylene
    elif ( smiles == 'OCCO'):         mechspecies = 'ETEG'  # ethylene glycol
    elif ( smiles == 'CC' ):          mechspecies = 'ETH'   # ethane
    elif ( smiles == 'C=O'):          mechspecies = 'HCHO'  # formaldehyde
    elif ( smiles == 'C=CC(=C)C' ):   mechspecies = 'ISO'   # isoprene (output of Chem.CanonSmiles)
    elif ( smiles == 'CO'):           mechspecies = 'MOH'   # methanol
    elif ( smiles == 'O=CO'):         mechspecies = 'ORA1'  # formic acid
    elif ( smiles == 'COO'):          mechspecies = 'OP1'   # methyl hydrogen peroxide
    elif ( smiles == 'C=Cc1ccccc1'):  mechspecies = 'STY'   # styrene (added in CRACMM2)
    elif ( smiles == 'CCc1ccccc1'):   mechspecies = 'EBZ'   # ethylbenzene (added in CRACMM2)
    elif ( smiles == 'CCC(C)=O' ):    mechspecies = 'MEK'   # methyl ethyl ketone
    elif ( smiles == 'C=CC(C)=O' ):   mechspecies = 'MVK'   # methly vinyl ketone
    elif ( smiles == 'Cc1ccccc1' ):   mechspecies = 'TOL'   # toluene
    elif ( smiles == 'C=CC=C' ):      mechspecies = 'BDE13' # 1,3 butadiene   
    elif ( smiles == 'C=CC=O' ):      mechspecies = 'ACRO'  # acrolein

    # Glyoxal and glycoaldehyde (here due to solubility alt mappings: ACD, ETEG)
    elif ( nC==2 and nO==2 and naldehyde>=1 ):     
                                      mechspecies = 'GLY'   # glyoxal
    # Propylene glycol
    elif ( nC==3 and nO==2 and nalcohol==2): mechspecies = 'PROG' # propylene glycol and other 3 carbon dialcohols                                

    # Low reactivity species (new to v0.1)
    elif ( koh < 3.5e-13 ):           mechspecies = 'SLOWROC' # low reactivity gas

    # Lumped terpene species 
    elif ( nC == 15 and nH == 24 ) and ( nCdblC>=1 ): mechspecies = 'SESQ' # sesquiterpenes
    elif ( tfmonoterpene and nCdblC==1 ): mechspecies = 'API' # a-pinene monoterpenes
    elif ( tfmonoterpene and nCdblC>=2 ): mechspecies = 'LIM' # limonene monoterpenes

    # Furans and dienes other than 1,3 BDE
    elif ( nfuran > 0 ):               mechspecies = 'FURAN' # furans and other dienes

    # Multi-ring aromatics (PAH and NAPH can be collapsed together if necessary)  
    # elif ( nbenzene >= 1 and log10cstar < 3.5 and (nO/nC) == 0 ): mechspecies = 'PAH'  # PAH and other lower-volatility aromatics (v0.1)
    elif ( nbenzene > 1 and nO/nC == 0 ): mechspecies = 'NAPH' # Naphthalene-like, PAH with 2 rings

    # SVOC species binned
    elif ( log10cstar < -1.5 ): # C* bin centered on 0.01 ug/m3 
        if OtoC > 0.6:
            mechspecies = 'ROCN2OXY8'
        elif OtoC > 0.3:
            mechspecies = 'ROCN2OXY4'
        elif OtoC > 0.1:
            mechspecies = 'ROCN2OXY2'
        else:
            mechspecies = 'ROCN2ALK'
    elif ( log10cstar < -0.5 ): # C* bin centered on 0.1 ug/m3
        if OtoC > 0.45:
            mechspecies = 'ROCN1OXY6'
        elif OtoC > 0.2:
            mechspecies = 'ROCN1OXY3'
        elif OtoC > 0.05:
            mechspecies = 'ROCN1OXY1'
        else:
            mechspecies = 'ROCN1ALK'
    elif ( log10cstar < 0.5 ): # C* bin centered on 1
        if OtoC > 0.3:
            mechspecies = 'ROCP0OXY4'
        elif OtoC > 0.1:
            mechspecies = 'ROCP0OXY2'
        else:
            mechspecies = 'ROCP0ALK'
    elif ( log10cstar < 1.5 ): # C* bin centered on 10^1
        if OtoC > 0.2:
            mechspecies = 'ROCP1OXY3'
        elif OtoC > 0.05:
            mechspecies = 'ROCP1OXY1'
        else:
            mechspecies = 'ROCP1ALK'
    elif ( log10cstar < 2.5 ): # C* bin centered on 10^2
        if OtoC > 0.1:
            mechspecies = 'ROCP2OXY2'
        else:
            mechspecies = 'ROCP2ALK'

    # Single-ring aromatics (excluding explicit species)
    elif ( nbenzene > 0 ): # Single-ring aromatics
        if ( naldehyde > 0 ):                mechspecies = 'BALD'     # Benzaldehyde and arom. aldehydes
        elif ( nC>=7 and nalcohol>=2 ):      mechspecies = 'MCT'      # methylcatechol
        elif ( nC>=7 and nalcohol>=1 ):      mechspecies = 'CSL'      # cresol
        elif ( nC==6 and nalcohol>=1 ):      mechspecies = 'PHEN'     # phenol
        elif ( log10cstar < 5.5 ):           mechspecies = 'VROCP5ARO' # C* bin centered on 10^5 (v0.1)
        elif ( log10cstar < 6.5 ):           mechspecies = 'VROCP6ARO' # C* bin centered on 10^6 (v0.1)
        # any single-ring aromatics that have not been mapped by rules above
        else:                                mechspecies = 'XYL'      # xylenes and other aromatics (CRACMM2)
    
    # Species with double bonds, not aromatic
    elif ( nCdblC>=1 and log10cstar < 5.5 ): mechspecies = 'VROCP5ARO' # C* bin centered on 10^5 (v0.1)
    elif ( nCdblC>=1 and log10cstar < 6.5 ): mechspecies = 'VROCP6ARO' # C* bin centered on 10^6 (v0.1) 
    elif ( nCdblC>=2 ):                      mechspecies = 'FURAN'    # some additional dienes (nC>=4 gauranteed)
    elif ( nCdblC==1 and ncarbonyl>=2 ):     mechspecies = 'DCB1'     # unsaturated dicarbonyls
    elif ( nCdblC==1 and nC==4 and naldehyde==1 ): 
                                             mechspecies = 'MACR'     # methacrolein (2.9e-11) (and crotonaldehyde)
    elif ( nCdblC==1 and naldehyde>=1 ):     mechspecies = 'UALD'     # unsaturated aldehydes (3.4e-11)
    # map OLI and OLT
    elif ( nCdblC==1 ):
        atoms = [a for a in m.GetAtoms()]
        bonds = [b for b in m.GetBonds()]
        bondtype = [b.GetBondType() for b in bonds]
        #        double bond at end                         bond isn't part of a ring    the double bond is to a carbon
        if bondtype[0]==rdkit.Chem.rdchem.BondType.DOUBLE and not bonds[0].IsInRing() and atoms[0].GetSymbol()=='C':
                                               mechspecies = 'OLT'
        elif bondtype[-1]==rdkit.Chem.rdchem.BondType.DOUBLE and not bonds[-1].IsInRing() and atoms[-1].GetSymbol()=='C':
                                               mechspecies = 'OLT'
        else:                                  mechspecies = 'OLI'
    
    # IVOC species binned
    elif ( ( log10cstar < 6.5 and nO/nC >= 0.1 ) or nSi > 0 ):  
                                       mechspecies = 'VROCIOXY' # Oxygenated IVOCs and any silanes/siloxanes
    elif ( log10cstar < 3.5 ): # C* bin centered on 1000 ug/m3
        if OtoC > 0.1:
            mechspecies = 'ROCP3OXY2'
        else:
            mechspecies = 'ROCP3ALK'
    # gas phase only for C* bins above 1000 ug/m3 so it can only be VROC (not AROC)
    elif ( log10cstar < 4.5 ): # C* bin centered on 10^4
        if OtoC > 0.1:
            mechspecies = 'VROCP4OXY2'
        else:
            mechspecies = 'VROCP4ALK'
    elif ( log10cstar < 5.5 ): # C* bin centered on 10^5
        if OtoC > 0.05:
            mechspecies = 'VROCP5OXY1'
        else:
            mechspecies = 'VROCP5ALK'
    elif ( log10cstar < 6.5 ): # C* bin centered on 10^6
        if OtoC > 0.05:
            mechspecies = 'VROCP6OXY1'
        else:
            mechspecies = 'VROCP6ALK'
    
    # Oxygenated species without double bonds (mapped in order of decreasing koh)
    elif ( naldehyde>=1 and nketone>=1 ):mechspecies = 'MGLY' # methylglyoxal and similar like C4H6O2 (1.5e-11)
    elif ( naldehyde>=1 ):               mechspecies = 'ALD'  # higher aldehydes (C>3) (1.98e-11)  
    elif ( nperoxide>=1 ):               mechspecies = 'OP2'  # higher organic peroxides (6.4e-12)  
    elif ( nalcohol>=1 and nketone>=1 ): mechspecies = 'HKET' # hydroxy ketone (3.0e-12) 
    elif ( nketone>=1 ):                 mechspecies = 'KET'  # all other ketones (2.9e-12) 
    elif ( nnitrate>=1 ):                mechspecies = 'ONIT' # organic nitrates (2.2e-12)
    elif ( nalcohol>=1 ):                mechspecies = 'ROH'  # C3 and higher alcohols (1.3e-12)  
    elif ( nacid>=1 ):                   mechspecies = 'ORA2' # acetic acid and higher acids (C>=2) (6.5e-13) 

    # HC Series, koh in cm3/s, 298 K, 1 atm
    elif ( koh < 3.4E-12 ):                     mechspecies = 'HC3'  # slow "alkanes"
    elif ( koh >= 3.4E-12 and koh <= 6.8E-12 ): mechspecies = 'HC5'  # medium "alkanes"   
    elif ( koh > 6.8E-12 ):                     mechspecies = 'HC10' # fast "alkanes"

    else: mechspecies = 'UNKCRACMM' # Species is unknown to CRACMM
    
    # append V or A to ROC* species to indicate gas (V) or particle (A)
    # except log10(C*)>3 are only in gas phase and already have the V added above
    if mechspecies[:3]=='ROC':
        # only 'gas', 'particle', or None allowed; otherwise exit with error
        try:
            assert phase in ['gas', 'particle', None]
        except AssertionError as e:
            err_msg = f'Invalid phase option {phase} supplied. Valid options are "gas", "particle", or None.'
            raise Exception(err_msg).with_traceback(e.__traceback__)
        if phase=='particle':
            append = 'A'
        else:
            append = 'V' # assumes gas if phase=None was supplied
        mechspecies = append + mechspecies
    
    # warning if species was not mapped
    if mechspecies == 'UNKCRACMM':
        unkcracmm_msg = (
            f'Species with SMILES {smiles_input} is unknown in CRACMM'
            + ' and has been mapped to UNKCRACMM.'
        )
        warnings.warn(unkcracmm_msg)
    
    return mechspecies

