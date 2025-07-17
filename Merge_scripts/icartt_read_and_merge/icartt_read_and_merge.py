"""Set of Functions that allow you to read in Icartt Files and Merge them.

#  DESCRIPTION
#
#     This script provides all the functionality you need to read in
#     ICARTT data as a xarray dataset & save it as a .nc file for
#     quick usage in python scripts. Options are provided to merge
#     data from multiple icartt files into a single file, and to
#     remap the time averaging. Meta data is attached from the icartt file.
#
#  NOTES
#
#  * parser is based on icartt v2.0 spec. most convenient source, imo:
#    https://www-air.larc.nasa.gov/missions/etc/IcarttDataFormat.htm
#
#  * and spec reference on earthdata;s page:
#    https://cdn.earthdata.nasa.gov/conduit/upload/6158/ESDS-RFC-029v2.pdf
#
#  ACKNOWLEDGEMENTS
#    This module was modified from a really nice & much fancier module
#    (ornldaac_icartt_to_netcdf) written by mcnelisjj@ornl.gov
#    that was intended to convert ICARTT v2 files into netCDF.
#    Thanks to them for doing the bulk of the work.
#
#  WRITTEN BY:
#    Dr. Jessica D.Haskins 
#    Email: jessica.haskins@utah.edu
#    GitHub: @jdhask
"""

import os
import numpy as np
import pandas as pd
import pickle


import ict_utils 
import unit_utils
import time_utils


def read_icartt(icartt_file: str, flt_num: int = -99, meta: dict = {},
                instr_name_prefix: bool = False, add_file_no: bool = False,
                delimiter: str = ',', line1_delim:str =',',
                meta_dict_replace_with: str = '_'):
    """Parse a single ICARTT file to a pandas dataframe."""
   
    # Get the header row number from the ICARTT.
    with open(icartt_file, "r", errors="ignore") as f : # ,encoding='ANSI') as f:
        header_row = int(f.readlines()[0].split(line1_delim)[0]) - 1

    # Parse the table starting where data begins (e.g. after the header).engine='python',
    print(icartt_file,header_row)
    df = pd.read_csv(icartt_file, header=header_row, delimiter=delimiter)#, encoding='utf-8')
    
    # Set possible error values to NaNs.
    df.replace(-9, np.nan, inplace=True)
    df.replace(-99, np.nan, inplace=True)
    df.replace(-999, np.nan, inplace=True)
    df.replace(-9999, np.nan, inplace=True)
    df.replace(-99999, np.nan, inplace=True)
    df.replace(-999999, np.nan, inplace=True)
    
    # Strip leading/tailing white space around variable names
    df.columns = [c.strip() for c in list(df.columns)]

    # Build/ append metadata from ICARTT to a dictionary file, add prefix to columns 
    df,meta = ict_utils.build_meta_dict(df, icartt_file, meta=meta,  flt_num=flt_num,
                            instr_name_prefix=instr_name_prefix,
                            line1_delim=line1_delim, add_file_no=add_file_no,
                            replace_with=meta_dict_replace_with)

    if add_file_no is True:
        # Create a column same length as data that contains the file #
        sz = len(df[df.columns[0]])  # get appropriate length
        fnum_arr = np.full(shape=sz, fill_value=flt_num, dtype=int)
        df['Flight_N'] = fnum_arr

    return df, meta  # dataframe with data, and df with metadata


def read_icartt_multileg(icartt_files: list = [], flt_num: int = -99, 
                          meta:dict = {}, instr_name_prefix: bool = False,
                          add_file_no: bool = False):
    """Parse multi-leg icartt files, combine into single df."""
    icartt_files=icartt_files[0] # pop list from 2d list
    
    # Sort the list of input ICARTTs.
    icartts = sorted(icartt_files)

    df = None  # Set an empty merged flight data data frame for this leg

    for ict in icartts:  # Loop over the dif ICARTT Legs
        # Parse individual file
        df_i, meta_i = read_icartt(ict, flt_num=flt_num, meta=meta,
                                   instr_name_prefix=instr_name_prefix,
                                   add_file_no=add_file_no)
        
        meta = meta_i  # update metadata file... gets appended upstream.
        
        # If df for merged is still None, set to first iter. Otherwise append.
        if df is None:
            df = df_i
        else:
            #df = df.append(df_i, ignore_index=True)
            df = pd.concat([df, df_i], ignore_index=True)

    return df, meta  # Return the merged df


def main_loop_parse_flights(DATA: dict):
    """Looper for parsing indv flights in a directory."""
    
    # Make groupings of standard and multileg flights.
    DATA['FLIGHTS'] = ict_utils.organize_standard_and_multileg_flights(DATA)

    ct = 1  # Initialize number of flights you're looping over.

    # Loop over the ICARTT files in the data directory.
    for flight, icartt in DATA['FLIGHTS'].items():

        # Tell the people which file you're processing:
        print("\n - [ {} / {} ] {} ".format(
            ct, len(DATA['FLIGHTS']), os.path.basename(flight)))

        if int(ct) == 1:  # Initialize empty dictionary to store metadata
            meta ={'Original_Filename':{}, 'Instruments': {}, 'Data_Info': {}, 'Instrument_Info': {},
                'PI_Info': {}, 'Uncertainty': {}, 'Revision': {}, 'Long_Name':{},
                'Stipulations': {}, 'Institution_Info': {}, 'Units': {},'Vars_Instrument': {},
                'Platform' : {} ,'PI_Contact_Info': {} , 'Associated Data': {} , 
                'Location': {} , 'ULOD_Flag': {} , 'ULOD_Value': {} , 
                'LLOD_Flag': {} , 'LLOD_Value': {} ,'Project_Info': {} , 
                'Other_Comments': {}}
        
        # If the type(icartt) is a list, must be MULTI-LEG flight.
        if type(icartt) is list:
            add_file_no = True
            # Handle multileg flight so it is merged to single df.
            df_data, new_meta = read_icartt_multileg(icartt_files=[icartt], flt_num=int(ct),
                                                      meta=meta, add_file_no=add_file_no, 
                                                      instr_name_prefix=DATA['PREFIX_OPT'])
            meta = new_meta  # Update the metadata dict to be the new one.
            icartt=icartt[0] # only need one file from here on out, not list of files. 
            
        # Else if the type(icartt) is a string, parse SINGLE flight (not leg)
        elif type(icartt) is str:
            add_file_no= True
            # Call the parse_icartt_table directly, pass existing meta dict.
            df_data, new_meta = read_icartt(icartt, flt_num=int(ct),
                                            meta=meta, instr_name_prefix=DATA['PREFIX_OPT'],
                                            add_file_no=add_file_no)
            meta = new_meta  # Update the metadata dict to be the new one.
            
        # Parse the string date columns in the indv icartt, and pick which
        # one to use, then convert it to a datetime object.
        df_data, times = time_utils.master_icartt_time_parser(df_data, icartt,
                                                   quiet=True,
                                                   remove_old_time=True, 
                                                   tz_default=DATA['TZ_DEFAULT'])

        if DATA['MODE'] == 'Merge_Beside':
            # To merge beside, must get all dfs must be on the SAME
            # time axis... so align them to a master timeline.
            # Aligns each file after its opened.
            tmln = DATA['MSTR_TMLN']  # wouldn't want to index twice!
            print(tmln)
            df_data = time_utils.align2master_timeline(df_data, tmln[0], tmln[1],
                                            tmln[2], quiet=True,
                                            datetime_index=False)
            df_data.index.name='datetime' # Assign name to index (need for xarray)
            # Should come out with a datetime index we can merge along!

        if ct == 1:  # If first icartt file, make the larger dataframe!
            df_all = df_data
        else:
            if DATA['MODE'] == 'Stack_On_Top':
                # For all subsequent loops append the new one UNDER the old df
                df_all = pd.concat([df_all, df_data], ignore_index=True)
            else:
                # For all subsequent loops append new columns to existing df
                # NEXTT To each other but same index. e.g. MEGE_BESIDE
                df_all = pd.concat([df_all, df_data], axis=1)

        ct += 1  # Update the counting variable.

    # Check if the User wants us to align the Stacked data to a master timeline
    # Aligns AFTER all icartts have been loaded in.
    if (DATA['MODE'] == 'Stack_On_Top') and (bool(DATA.get('MSTR_TMLN'))
                                             is True):
        tmln = DATA['MSTR_TMLN']  # wouldn't want to index twice!
        df_all = df_all.set_index(['datetime'])
        
        df_all = time_utils.align2master_timeline(df_all, tmln[0], tmln[1], tmln[2],
                                       quiet=True, datetime_index=True)
        
        
    elif (DATA['MODE'] == 'Stack_On_Top'):  # Make file # and datetime indexes.
        df_all = df_all.set_index(['datetime'])

    # Clean column names from weird chars. 
    cur_cols=df_all.columns
    for i in range(0, len(cur_cols)):
        good_name= ict_utils.char_cleaner(cur_cols[i],  ignore= '-')
        df_all.rename(columns={cur_cols[i]: good_name}, inplace=True)
    
    return df_all, meta


def icartt_merger(data_directory: str,
                  mode_input: str,
                  master_timeline: list = [],
                  output_directory: str = '',
                  output_filename: str = 'icartt_merge_output',
                  prefix_instr_name: bool = True, 
                  units_in_data_info: bool = False,
                  tz_default:str ='',# preffered TZ when none found. 
                  save_meta:bool=True): 
    """Merge a directory of icarrts into an Xarray .nc file.
    
    # ========================================================================
    # ========================    INPUTS   ===================================
    # ========================================================================
    #
    #   (1) data_directory - A string containing the absolute path to a folder
    #                      which contains all the individual icartt files that
    #                      you wish to merge together.
    #
    #   (2) mode_input - A string describing HOW you would like to merge these
    #                    icartt files in the data_directory together. Only 2
    #                    valid options are supported right now.  Either
    #                    "Stack_On_Top" or "Merge_Beside".
    #
    #     "Stack_On_Top": Each icartt file is for a different date,
    #      but contains data from multiple instruments or mutltiple
    #      measurements and  you want that data in a single
    #      file (e.g. indexed by time, and File/Flight #). Contents of
    #      individual icarrt files will be "stacked on top" of one another.
    #
    #     "Merge_Beside": Each icartt file is for the entire
    #      sampling period, but contains different measurements.
    #      You want to have all of these differnt measurments
    #      on the same time base, throughout the whole period. The contents
    #      of each icartt file will be "merged beside" one another.
    #
    #   (3) master_timeline - OPTIONAL if "Stack_On_Top", required if
    #                        "Merge_Beside". It is list with 3 items:
    #
    #       -  Startdate_str:  A string containing the start date of the
    #                        "mastertimeline" that all data  will be
    #                        merged to.Format is 'YYYY-MM-DD HH:MM:SS'
    #       -  Enddate_str:  A string containing the end date of the
    #                        "mastertimeline" that all data  will be
    #                        merged to. Format is 'YYYY-MM-DD HH:MM:SS'
    #       - Averaging_Step: An integer that is the number of seconds
    #                        for each timestep in between startdate and
    #                        end date. So 120 for a 2 minute average.
    #
    #    (4) output_drectory - OPTIONAL string containing the  abs path where
    #                         the output file will be written. If not set, the
    #                        output will be to stored in the  input data_dir.
    #
    #    (5) output_filename - OPTIONAL string containing what you'd like the
    #                         output file to be called (not including its
    #                         extension). Default is 'icartt_merge_output'
    #
    #    (6) prefix_instr_name - OPTIONAL boolean value indicating whether you
    #                         would like to append the instrument name
    #                         contained in the icartt file to all the var
    #                         names. Default is True since when merging
    #                         icartt files it is common to have some PI's
    #                         measuring the same items & naming them the same.
    #
    #    (7) units_in_data_info - OPTIONAL boolean value indicating whether the 
    #                         units of values are contained in the "DATA INFO" 
    #                         section of the ICARTT File (as is sometimes case 
    #                         for merged files... ). Default is false. 
    
    #    (8) tz_default - OPTIONAL str indicating the default timezone to assume 
    #                         measurements are in if none is found. 
    #
    #    (9) save_meta - OPTIONAL bool indicating whetherto save the output metadata 
    #                         in a pickle file or not. Default is TRUE (to save it). 
    #
    # ========================================================================
    """
    # Format the input for easier referencing.
    inputs = {'DIR_ICARTT': data_directory,
              'DIR_OUTPUT': output_directory,
              'O_FILENAME': output_filename,
              'MODE': mode_input,
              'PREFIX_OPT': prefix_instr_name,
              'MSTR_TMLN': master_timeline,
              'TZ_DEFAULT': tz_default} # Tz to assume if none provided by default. 

    # Make sure you got appropriate inputs from the user, retrieve icartt files
    DATA = ict_utils.handle_input_configuration(inputs)
    
    # Loop through parsing the flights & collecting them in a single dataframe.
    df, meta = main_loop_parse_flights(DATA)
    # # remove duplicate datetime cols. (messes up xarray conversion.)
    df = df.loc[:, ~df.columns.duplicated()]
    
    if units_in_data_info is True:  # attach units from data info if they ask. 
        df, meta = unit_utils.units_from_datainfo(df, meta, instr_name_prefix= prefix_instr_name)
    
    if df.index.name is None:
        df.index.name = 'datetime'
    
    meta= unit_utils.unit_cleaner(meta) # Clean up the metadata unit names to be correct formats. 
    print('Converting output to X-Array, please wait!') # This can take a while.
    ds_i0 = df.to_xarray() # Convert the dataframe to an xarray
    
    ds_i=ict_utils.attach_meta(ds_i0, meta) # Attach info from metadata file to xarray
    
    # Datetime as an index is formatted weird after xarray convert. so fix it
    # Pull time from the DF index, convert to something we can save as a series/var
    t = pd.to_datetime(df.index.get_level_values('datetime')).to_series().reset_index(drop=True)
    
    # Assign a new variable to xarray that is time.
    ds_i['time'] = ("datetime", t)
    
    # Tell xarray that we want time to be a coordinate.
    ds_i = ds_i.set_coords('time')
    
    # And tell it to replace index # with time as the preferred dimension.
    ds_i = ds_i.swap_dims({"datetime": "time"})
    
    # And now drop the old coordinate index.
    ds_i= ds_i.reset_coords('datetime', drop=True)
    
    # Save the output xarray with re-mapped names to a netcdf
    filename=DATA['DIR_OUTPUT'] + DATA['O_FILENAME']+'.nc'
    ds_i.to_netcdf(filename)
    
    # Save the rest of metadata to a pickle if asked. 
    if save_meta==True: 
        filename_meta = DATA['DIR_OUTPUT'] + DATA['O_FILENAME'] + '_meta.pickle'    
        with open(filename_meta, 'wb') as handle:
            pickle.dump(meta, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    # Tell the people where you saved it.
    print('Output nc file and additional metadata saved at:' + filename)
    
    return ds_i, meta, filename

