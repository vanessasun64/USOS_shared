#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 14:35:49 2025

@author: u6044586
"""

import sys
import numpy as np
import pandas as pd 
import xarray as xr
import os
import re
import glob
from datetime import datetime, timedelta

sys.path.insert(0,'/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_merges/icartt_read_and_merge/')
from icartt_read_and_merge import icartt_merger
#Note, line above has changed from importing ict to importing icartt_merger so the otter/flight functions need to be changed to reflect this and have  as of June 2025, 
#since Vanessa is only using parked and drives data

def list_directories(path):
    directories = []
    with os.scandir(path) as entries:
        for entry in entries:
            if entry.is_dir() and entry.name.isdigit():
                directories.append(entry.name)
    return directories

def get_master_timeline(date_str, time_int):
    # Parse the date string
    date_obj = datetime.strptime(date_str, '%Y%m%d')
    
    # Format the date into the desired string format
    formatted_date = date_obj.strftime('%Y-%m-%d 00:00:00')
    
    # Calculate the next day
    next_day = date_obj + timedelta(days=1)
    
    # Format the next day into the desired string format
    formatted_next_day = next_day.strftime('%Y-%m-%d 00:00:00')
    
    master_timeline=[formatted_date, formatted_next_day,time_int]
    
    return master_timeline

def remove_meta_pkl_files(dirpath): 
    
    # Use glob to find all .pkl files in the directory and subdirectories
    pkl_files = glob.glob(os.path.join(dirpath, '**', '*.pickle'), recursive=True)

    # Delete each .pkl file
    for file in pkl_files:
        try:
            os.remove(file)
            parts = file.split("USOS_2024" + os.sep, 1)
            print(f"Deleted: {parts[1]}")
        except Exception as e:
            print(f"Error deleting {parts[1]}: {e}")
    return
        
def merge_drives(res='10s',merge_path=''):
    
    # Path to raw USOS Mobile lab data for each drive date: 
    usos_path='/uufs/chpc.utah.edu/common/home/haskins-group1/data/Campaign_Data/Raw_Data/USOS_2024/R0/CSL_MobileLab_Driving/'
    
    # Set Merged data path: 
    if len(merge_path)==0: 
        merge_path=os.path.join(usos_path,'/merged/',res,'/')
        print('Merged ', res, 'data will be saved at: \n\t', merge_path) 
        
    if res=='1s': 
        interval_s=1 # 1 s
    elif res=='10s': 
        interval_s=10 # 10s
    elif res=='30s': 
        interval_s=30 # 30s
    elif res=='1min': 
        interval_s=60 # 60s  
        
    # Get a list of all dates of drives: 
    drive_dates=list_directories(usos_path)
    
    # Loop over each drive
    for drive_i in drive_dates: 
        
        # Create the "master timeline" for that day. in XXs intervals... 
        master_timeline_i= get_master_timeline(drive_i, interval_s)
        
        # Merge individual data beside one another for this date: 
        ds, meta, filename= icartt_merger(data_directory=usos_path+drive_i+'/',
                                              master_timeline=master_timeline_i,
                                              mode_input="Merge_Beside",
                                              output_directory=merge_path,
                                              output_filename= 'USOS_'+res+'_Merged_Drive_'+drive_i, # What output file will be named... 
                                              prefix_instr_name=False, # instr name in prefix already! 
                                              units_in_data_info= False,
                                              tz_default='UTC')
        
    # Remove all the meta data pickle files: 
    remove_meta_pkl_files(merge_path)
        
    return 

def merge_parked(res='30min',merge_path=''):
    # Path to raw USOS Mobile lab data for each: 
    usos_path=  '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_merges/R0/CSL_MobileLab_Parked/raw/'
    
    # Set Merged data path: 
    if len(merge_path)==0: 
        merge_path=os.path.join(usos_path,'/merged/',res,'/')
        print('Merged ', res, 'data will be saved at: \n\t', merge_path) 
        
    if res=='30min': 
        interval_s=30*60 # 30mins * 60s = 1800s. 
    elif res=='1hr': 
        interval_s=60*60 # 60mins * 60s = 3600s. 
   
    # Get a list of all dates of drives: 
    drive_dates=list_directories(usos_path)
    
    # Loop over each drive
    for drive_i in drive_dates: 
        
        # Create the "master timeline" for that day. in XXs intervals... 
        master_timeline_i= get_master_timeline(drive_i, interval_s)
    
        # Merge individual data beside one another for this date: 
        ds, meta, filename= icartt_merger(data_directory=usos_path+drive_i+'/',
                                              mode_input="Merge_Beside",
                                              master_timeline=master_timeline_i,
                                              output_directory='/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_merges/R0/CSL_MobileLab_Parked/merged/',
                                              output_filename= res+'/USOS_'+res+'_Merged_Parked_'+drive_i, # What output file will be named... 
                                              prefix_instr_name=False, # instr name in prefix already! 
                                              units_in_data_info= False,
                                              tz_default='UTC')
    
    # Remove all the meta data pickle files: 
    remove_meta_pkl_files(merge_path)
        
    return merge_path

def organize_twin_otter_by_date(base_dir, destination_dir):
    """Function to organize files by date and leg instead of by instrument regardless of revision #"""    
    
    # Define the pattern to extract the flight date and leg/revision #
    pattern = re.compile(r'_(\d{8})_(R\w+)(_L\d+)?\.ict$')
    
    # List all subdirectories within the base directory that don't start with a # (already organized by date if so!)
    subdirs = [d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d)) and not d[0].isdigit() and d!='merged']

    for subdir in subdirs:
        subdir_path = os.path.join(base_dir, subdir)
        for filename in os.listdir(subdir_path):
            if filename.endswith('.ict'):
               
                # Extract the date and leg using the regular expression
                match = pattern.search(filename)
                if match:
                    date_str = match.group(1) # get date 
                    rev_str =match.group(2) # revision & leg combined
                    if '_' in rev_str: # seperate into leg if it has it in that group! 
                        pts=rev_str.split('_')
                        rev_str=pts[0]
                        leg_str='_'+pts[1]
                    else: 
                        leg_str='' # otherwise leave leg as a blank. 
            
                    # Construct the target directory name based on date and leg
                    target_dir_name = f'{date_str}{leg_str}'
                    target_dir_path = os.path.join(destination_dir, target_dir_name)
                    
                    # Create target directory if it doesn't exist
                    os.makedirs(target_dir_path, exist_ok=True)
                    
                    # Source and target file paths
                    source_file_path = os.path.join(subdir_path, filename)
                    target_file_path = os.path.join(target_dir_path, filename)
                    
                    # Copy the file to the target directory
                    shutil.copy2(source_file_path, target_file_path)
                    
    return 

def get_flt_start_stop_times(base_dir, avg_step=10):
    """Function to pull out min start/max stop time & time resolution of all indv files on each flight
    contained in a dir organized by date underneath the input base_dir. Used to build master timeline 
    for merge beside on flights when you want to keep the start/stop of the flights but don't know it... """ 
    
    # Initialize output dict to hold info: 
    flt_times={}
    
    # Get a list of all date/leg subdirs in folder: 
    flts=[os.path.join(base_dir,d) for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d)) and  d[0].isdigit() and d!='merged']
    
    for subdir_path in flts: 
        # Get the the date from the folder name: 
        fldr = os.path.basename(subdir_path)
        date_str=fldr.split('_')[0]if '_' in fldr else fldr 
        date = pd.to_datetime(date_str, format='%Y%m%d')
            
        # Get a list of all files in that folder (from each instr): 
        files= [os.path.join(subdir_path,filename) for filename in os.listdir(subdir_path) if filename.endswith('.ict')]
        
        # Loop over all files with this date & get min start time/max stop time:
        start_tm=[];stop_tm=[]; res=[]; time_info={} 
        for file_i in files: 
            
            df, meta=ict.read_icartt(file_i, instr_name_prefix=True, add_file_no= False)    
            time_var=df.columns[0] # Get independent time variable 
            start_tm.append(df[time_var][0]) # pull out start time 
            stop_tm.append(df[time_var][len(df)-1])# pull out stop time 
            res.append(df[time_var].diff().mean())# pull out average time resolution (s) 
            
        # Pull min/max from all files & convert from s since midnight to datetime: 
        min_start=date+pd.to_timedelta(np.nanmin(start_tm), unit='s')
        max_stop=date+pd.to_timedelta(np.nanmax(stop_tm), unit='s') 
        min_res=np.nanmin(res); max_res=np.nanmax(res)
        
        # Store time info about this date/leg in output dict: 
        flt_times[fldr]={}
        flt_times[fldr]=[min_start.strftime('%Y-%m-%d %H:%M:%S'),
                        max_stop.strftime('%Y-%m-%d %H:%M:%S'),
                        avg_step]
    
    return flt_times

def merge_twin_otter(avg_step='10s'): 
    """Function to organize, merge beside, and merge together all Twin Otter data""" 
    
    base_dir='/uufs/chpc.utah.edu/common/home/haskins-group1/data/Campaign_Data/Raw_Data/USOS_2024/R0/Twin_Otter/'

    # Create dir to store data organized by date if it doesn't yet exist: 
    #os.makedirs(base_dir+'/by_date/', exist_ok=True)

    # Twin Otter data was organized by instrument... re-organize it by date/leg instead
    # so we can merge beside for each before we merge on top for all! 
    #organize_twin_otter_by_date(base_dir+'/by_instrument/',base_dir+'/by_date/')
    
    # Get the start/stop times for each flight from the data to use in merge beside 
    # to create input master_timeline list (at input for avg_Step XXs) invervals
    #avg_step_num=np.int64(avg_step.replace('s',''))
    #flt_times=get_flt_start_stop_times(base_dir+'by_date/', avg_step=avg_step_num)
    
    # Create output merge beside dir if it doesn't already exist: 
    mrg_beside_dir=base_dir+'merged/merged_beside_'+avg_step+'/'
    #os.makedirs(mrg_beside_dir, exist_ok=True)
    
    # Merge data from each flight BESIDE one another first on each date/leg:  
    # for flt_i in list(flt_times.keys()): 
    #       ds, meta, filename= ict.icartt_merger(data_directory=base_dir+'by_date/'+flt_i,
    #                                       mode_input="Merge_Beside",
    #                                       master_timeline=flt_times[flt_i], 
    #                                       output_directory=mrg_beside_dir,
    #                                       output_filename= 'USOS_'+avg_step+'MrgBeside_TwinOtter_'+flt_i, # What output file will be named... 
    #                                       prefix_instr_name=False, # instr name in prefix already! 
    #                                       units_in_data_info= False,
    #                                       save_meta=False,
    #                                       tz_default='UTC')
    
    # Get a list of all the individual merged beside files: 
    nc_files=[os.path.join(mrg_beside_dir,filename) for filename in os.listdir(mrg_beside_dir) if filename.endswith('.nc')]
            
    # Define pattern to extract the flight date and leg using a regular expression
    pattern = re.compile(r'_(\d{8})(_L\d+)?\.nc$')
    
    # Open each merge beside dataset and add info about flight leg/ flight num. 
    # Then concatenate them all together to make the final merged file: 
    ds_list=[]
    for i,file in enumerate(nc_files):
        # Pull out date & flight leg #: 
        match = pattern.search(file)
        if match:
            date_str = match.group(1) # get date
            if match.group(2): 
                leg= np.int64(match.group(2).replace('L','').replace('_',''))
            else: 
                leg=np.nan
            
        # Open Dataset: 
        ds_i=xr.open_dataset(file) 
        
        #Grab coords/dims of one var: 
        coords=ds_i['Flight_N'].coords
        dims= ds_i['Flight_N'].dims
        
        # Add info on flt_leg and flt num this data came from: 
        ds_i['Flt_Leg']=xr.DataArray(np.ones(len(ds_i['Flight_N']))*leg, coords=coords, dims=dims)
        ds_i['Flt_Num']=xr.DataArray(np.ones(len(ds_i['Flight_N']))*i, coords=coords, dims=dims)
        ds_i['Flt_Date']=xr.DataArray(date_str*len(ds_i['Flight_N']), coords=coords, dims=dims)

        # Remove the old dummy flight_N var created when we merged beside... meaningless: 
        ds_i=ds_i.drop_vars('Flight_N')
        
        # Append to list of datasets to concatenate: 
        ds_list.append(ds_i)
    
    # Concatenate all data from every flight into a single merged file: 
    ds= xr.concat(ds_list, dim='time')
       
    # Save output file: 
    outfile= base_dir+'merged/USOS_'+avg_step+'Merged_TwinOtter.nc'
    ds.to_netcdf(outfile)
    
    # Tell them where it saved it: 
    print(f'Merged Twin Otter data ({avg_step}) saved at: \n\t {outfile}')
    
    return ds

# # Merge all parked indv icartts into single day files w/ 30min average: 
merge_parked(res='1hr')

# # Merge all parked indv icartts into single day files w/ 10s average: 
#merge_drives(res='10s')

# # Merge all twin otter indv icartts into merged beside & then merge all at 10s average: 
#ds = merge_twin_otter(avg_step='10s')



