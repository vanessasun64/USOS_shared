#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 24 15:04:34 2023

@author: u6044586
"""
import pandas as pd 
import numpy as np 

def unit_cleaner(meta):
    """Function to take strangely formatted units in icartts and make nice."""
    bads=['ppb', 'ppm', 'ppt', 
            'ugm3','ug/m3', ' ug m3', 'ug_m3', 'ug_m3_', 'ug_m_3', 'ug_sm-3', 'ugm-3', 'ug  std  m-3',
            'ug std m-3', 'ng std m-3',
            'micrograms  per  cubicmeter', 
            'micrograms per ambient cubicmeter', 
            'micrograms  per  ambient  cubicmeter',
            'parts per billion by volume',
            'parts  per  billion  by  volume',
            'ng_per_m3', 
            'inverse  megameter', 
            'microns', 'nanometers',
            'um2/cm3',
            'um3/cm3','um3 cm-3',
            'nanometers'
            '/s']
    
    goods= ['ppbv', 'ppmv', 'pptv',
        'μg/m3', 'μg/m3','μg/m3','μg/m3', 'μg/m3', 'μg/m3', 'μg/m3','μg/m3','μg/m3',
        'μg/m3',' ng/m3',
        'μg/m3', 
        'μg/m3_amb','μg/m3_amb',
        'ppbv', 'ppbv',
        'ng/m3',
        'Mm^{-1}', 
        'μm', 'nm',
        'μm2/cm3',
        'μm3/cm3','μm3/cm3',
        's^{-1}']
  
    for key in meta['Units']: 
        value= meta['Units'][key]
        for i in range(0, len(bads)): 
                
            if bads[i] in value: 
                meta['Units'][key]=goods[i]
                
    return meta    


def units_from_datainfo(df, meta, instr_name_prefix: bool = False  ):
    """Fill in units from Data Info if no units were found in the ICARTT file."""
    keys=meta['Units'].keys()
    for cols in df.columns: # Loop over columns in df
        if cols in keys: 
            current_units= meta['Units'][cols] # get existing unit for this 
            this_instr= meta['Vars_Instrument'][cols] # and its instrument.
                
            if current_units  =='': # doesn't have assigned units. 
                meta['Units'][cols]=meta['Data_Info'][this_instr]
                    
    return df, meta 

