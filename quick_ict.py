#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 17 14:11:47 2025

@author: u6044586
"""
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
from datetime import date 

p='/uufs/chpc.utah.edu/common/home/haskins-group1/data/Campaign_Data/Raw_Data/USOS_2024/R0/CSL_MobileLab_Parked/20240717/USOS-PTR_MobileLabGround_20240717_RA.ict'

with open(p) as file:
    lines=file.readlines() 
    
hd_idx= int(lines[0].split(',')[0].strip())-1

headers=[h.strip() for h in lines[hd_idx].split(',')]

data={var: [] for var in headers}

for ln in range(hd_idx+1, len(lines)): 
    ln=lines[ln].strip() 
    items=ln.split(',') 
    
    for i,var in enumerate(headers): 
        try:
            val=float(items[i])
            if round(val) == -9999:
                val=np.nan 
        except ValueError: 
            val= items[i]
            
        data[var]= data[var]+[val]

    
specific_date = date(2024, 7, 17)
deltas = pd.to_timedelta(data['Time_Start'], unit='s')
base_timestamp = pd.to_datetime(specific_date)
datetimes = base_timestamp + deltas
data['Datetimes']= datetimes


plt.plot(data['Datetimes'], data['C5H8_NOAAPTR_ppbv'])

