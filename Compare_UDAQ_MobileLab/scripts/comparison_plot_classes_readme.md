Inputs:

### Species names 

* species_name: A string with species to plot
       Formatting notes:

       * Ozone should be in format 'Ozone'
       * VOC Options:
              'Ethane', 'Ethylene', 'Propane', 'Propylene', 'Acetylene', 
              'nButane', 'Isobutane', 'trans2Butene', 'cis2Butene', 
              'buta13diene', 'nPentane', 'Isopentane', 'Pent1ene', 
              'trans2Pentene', 'cis2Pentene', 'x3_Methylpentane', 'nHexane', 
              'nHeptane', 'nOctane', 'nNonane', 'nDecane', 
              'Cyclopentane', 'Isoprene', 'x22_Dimethylbutane', 'x224_Trimethylpentane', 
              'Methylcyclohexane', 'Methylcyclopentane', 'But1ene', 'x2_Methylpentane', 
              'mpXylene', 'Benzene', 'Toluene', 'Ethylbenzene', 
              'oXylene', 'x135_Trimethylbenzene', 'x124_Trimethylbenzene', 'nPropylbenzene',
              'Isopropylbenzene', 'oEthyltoluene', 'pEthyltoluene', 'Styrene', 
              'x123_Trimethylbenzene'

### Labels for plots 
* var_name_modifications: A string with any modification to the string for species_udaq. Will be used for plot axes or labels. For example, you could change 'oXylene' input for species_udaq to the more formal 'o-Xylene' for var_name_modifications.

Variable List:
mobilelab_o3
df_udaq_o3_load
usos_vars: a list of VOC species measured during USOS with Formaldehyde from Piccaro and Benzene, Toluene, Styrene, Isoprene from PTR. 
df_ml_data: a pandas dataframe


