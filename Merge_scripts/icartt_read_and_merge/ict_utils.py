#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 24 15:00:26 2023

@author: u6044586
"""

import os
import re
import sys
import csv
import numpy as np


def attach_meta(ds, meta, names_flexible: bool = False):
    """Function to take meta data and append units to attribute of ds"""
    print('Attaching MetaData to X-Array.')
    count=0  # Attached Units to meta data 
    printed_already=[]  # list of printed % s to terminal to prevent dupes
    

    tlen=len(ds.data_vars)
    for varname in ds.data_vars:
        
        # Get case insensitive list of units  dict keys. 
        units_lower = {k.lower(): v for k,v in meta['Units'].items()}
        if varname.lower()  in units_lower:
               value = units_lower[varname.lower()]
               ds[varname]= ds[varname].assign_attrs({'Units': value})
        else: 
            ds[varname]= ds[varname].assign_attrs({'Units': 'No units found for this variable...'})
        
        # Get case insensitive list of units  dict keys. 
        long_names_lower = {k.lower(): v for k,v in meta['Long_Name'].items()}
        if varname.lower()  in long_names_lower:
               value = long_names_lower[varname.lower()]
               ds[varname]= ds[varname].assign_attrs({'Long_Name': value})
        else: 
            ds[varname]= ds[varname].assign_attrs({'Long_Name': 'No long name found for this variable...'})

        # Get case insensitive list of units  dict keys. 
        instr_lower = {k.lower(): v for k,v in meta['Vars_Instrument'].items()}
        if varname.lower()  in instr_lower:
            value = instr_lower[varname.lower()]
            ds[varname]= ds[varname].assign_attrs({'Instr_from_file': value})
            
            # Other meta data is stored by instrument name so... attach here. 
            other_metas=list(meta.keys())
            other_metas.remove('Vars_Instrument')
            other_metas.remove('Units')
            other_metas.remove('Long_Name')

            for title in other_metas:
                title=str(title)
                if title != 'Instruments':
                    ldict = list(meta[title].keys())
                    if value in ldict:
                        ds[varname]= ds[varname].assign_attrs({title: meta[title][value]})
                    else: 
                        ds[varname]= ds[varname].assign_attrs({title: 'Not found'})
                          
        else: 
            ds[varname]= ds[varname].assign_attrs({'Instr_from_file': 'No Instrument found for this variable...'})

        # Print off someoutput for the people to see some progress!
        percent=  np.round((count/ tlen) * 100)
        if (percent % 5 == 0) and (percent != 0): # % processed is multitple of 5
            if count==0:
                print('Progress attaching units: ', np.round((count/ tlen) * 100),' %')
            else: 
                if percent not in printed_already: #Check if you printed this % to temrianl
                    print(np.round((count/ tlen) * 100),' %')
                    printed_already.append(percent) # append this percent to list of printed vals. 
       
        count=count+1
    
    return ds


def warn(message: str = "( miscellaneous warning )"):
    """Print the input error message and exits the script with a failure."""
    # Args: message (str): warning message gets printed during eval and ignored
    print("   WARN: {}. Skipping.".format(message))


def exit_with_error(message=str):
    """Print the input error message and exits the script with a failure."""
    sys.exit(print("ERROR: {} -- Abort.".format(message)))


def crawl_directory(path: str, extension: str = None):
    """Crawl an input directory for a list of ICARTT files.\
     Parameters:------------------\
        path (str): full path to an input directory.\
        ext (str): An optional extension to limit search.\
    Returns:  A list of paths to data files (strings)."""
    selected_files = []  # Create an empty list

    for root, dirs, files in os.walk(path):  # Walk directory.
        for f in files:  # Loop over files,
            fext = os.path.splitext(f)[1]  # Get the extension.

            # If file matches input extension or if no extension given,
            if extension is None or extension == fext:

                # Join to root for the full path.
                fpath = os.path.join(root, f)
                # Add to list.
                selected_files.append(fpath)

    # Return the complete list.
    return selected_files


def organize_standard_and_multileg_flights(DATA: dict):
    """Organize the Multileg flights & parse them."""
    # A regular expression catches the multi leg flight suffix.
    multileg_regex = re.compile('_L[0-9].ict')

    # A dictionary stores the output filename and legs as child list.
    flights = {}

    for ict in DATA['ICARTT_FILES']:  # Loop over all files in the directory.
        # If regular expression is not matched anywhere in string,
        if re.search(multileg_regex, ict) is None:
            # Add to list of "standard" flights (e.g. not a leg)
            flights[ict] = ict

        else:  # Else if regular expression is matched in string.
            # The output file won't have the suffix.
            output_filename = ict[:-7] + ".ict"

            # Add this file to the dict of multi-leg flights.
            if output_filename not in flights:
                flights[output_filename] = [ict]
            else:
                flights[output_filename].append(ict)

    return flights  # Return the organized flights as a dictionary.


def char_cleaner(mystring, ignore: list = [], replace_with: str ='_'):
    """Clean up gross strings from weird characters."""
    after = mystring.strip()  # strip all leading/trailing whitespace

    # Then, replace common representations with a word.
    after = after.replace('%', 'percent')
    after = after.replace('_+_', '+')
    after = after.replace('-->', '_to_')
    after = after.replace('->', '_to_')

    # A list of bad chars we don't want in our string.
    bad_chars = [' ', ',', '.', '"', '*', '!', '@', '#', '$', '^', '&',
                  '=', '?', '/', '\\', ':', ';', '~', '`', '<',
                 '>', '{', '}']

    for i in range(0, len(bad_chars)):
        if bad_chars[i] not in ignore:  # don't replace chars they want
            after = after.replace(bad_chars[i], replace_with)
    
    after = after.replace('_to_', '->')
    
    
    return after


def build_meta_dict(df, icartt_file: str, meta: dict = {}, flt_num: int = -99, 
                     instr_name_prefix: bool = False, line1_delim: str = ',', 
                     replace_with: str = '_', add_file_no: bool = False):
    """Take and combines metadata from different\
    icartt files into a dictionary. So you can access metadata from an\
    individual icartt file by typing in its instrument name & flt #."""
    
    # If meta empty, initialize the dict.
    if bool(meta) is False:
        meta = {'Original_Filename': {}, 'Long_Name':{},
                'Instruments': {}, 'Data_Info': {}, 'Instrument_Info': {},
                'PI_Info': {}, 'Uncertainty': {}, 'Revision': {},'Vars_Instrument': {},
                'Stipulations': {}, 'Institution_Info': {}, 'Units': {},
                'Platform' : {} ,'PI_Contact_Info': {} , 'Associated Data': {} , 
                'Location': {} , 'ULOD_Flag': {} , 'ULOD_Value': {} , 
                'LLOD_Flag': {} , 'LLOD_Value': {} ,'Project_Info': {} , 
                'Other_Comments': {}}
        
    # Open the file.
    with open(icartt_file, "r",errors="ignore") as f:  # Get number of header rows
        header_row = int(f.readlines()[0].split(line1_delim)[0]) - 1
    df_vars=[] 
    for i in range(0, len(df.columns)): 
        df_vars.append(char_cleaner(df.columns[i], ignore= '-', replace_with=replace_with))
    still_units=True # boolean to see if we're still in the part of the header that has units
    
    with open(icartt_file, "r",errors="ignore") as f: # read line by line.
        reader = csv.reader(f)
        ln_num = 0  # intitalize line counting var.
        for row in reader:
            line = " ".join(row)  # read line by line.

            # Icartt splits headers with a ":", use that to split them.
            before, sep, after = line.rpartition(":")
            after = char_cleaner(after, ignore=':')  # Pass to string cleaner.
            
            # Line 9 usually holds the "time Variable. ICARTT standard files should 
            # have time in seconds after midnight, but some dummies use fractional day of year 
            # instead. Here while reading the header, we'll check if they've done that 
            # and if they have, convert this column to seconds after midnight. 
            if (('fractional ' in line.lower()) or ('julian ' in line.lower())):
                col_name, sep, after = line.rpartition("  ")
                # from julian day to fractional day... 
                df[col_name]= df[col_name]-np.floor(df[col_name][0])
    
                df[col_name]=df[col_name]*24*60*60 # from day to seconds. 
            
            # First 3 lines have set parameters in ICARTT Files.
            if ln_num == 1:
                PI = after
            if ln_num == 2:
                Institution = after
            if ln_num == 3:
                Instrument = after
            
            # Once you know the instrument, you can start to build the dict
            # (becase we are using the instrument part of the dict index we
            # can't do it until we get to this line number. )
            if ln_num == 3:
                meta['PI_Info'][Instrument] = PI
                meta['Institution_Info'][Instrument] = Institution
                meta['Original_Filename'][Instrument]=icartt_file
                
                vals= meta['Instruments'].values()
                if Instrument not in vals: 
                    meta['Instruments'][len(vals)+1]= Instrument
                
            # The rest of the meta data is on arbitrary line #s based on how
            # much info the author of the ICARTT included, so just parse the
            # string to ID which row that is contained on. Then,  append info
            # from  this file into the meta dictionary indexed on the
            # instrument name           
            meta_icartt= ['Data_Info', 'Uncertainty', 'Revision', 
                                'Instrument_Info', 'Stipulations', 'Platform',
                                'PI_Contact_Info', 'Associated Data', 'Location', 
                                'ULOD_Flag', 'ULOD_Value', 'LLOD_Flag', 'LLOD_Value',
                                'Project_Info', 'Other_Comments']
            
            if 'PI_CONTACT_INFO' in before: 
                still_units= False 
            
            for t in range (0, len(meta_icartt)): 
                 if meta_icartt[t].upper() in before:
                     meta[meta_icartt[t]][Instrument] = after
            
            # --------------  Keep UNITS of EACH df varname:  ----------------
            # Icartt splits variable names from units with a space, so split name of var from unit.
            line_list= str(line).split(' ',1)

           
            # Check to see that it is in a spot where we expect varname sapce unit. 
            if  (ln_num < header_row-1) and (ln_num > 11) and still_units is True: 
                
                # Rip the varname 
                intro_line_ripped_varname=char_cleaner(line_list[0], ignore= ['-'] , replace_with=replace_with)
            
                # Do lower case comparison of df_vars, line with varname ripped. 
                df_vars_lower= [v.lower() for v in df_vars]
                intro_line_ripped_varname_lower= intro_line_ripped_varname.lower()

                for v in range(0, len(df_vars_lower)):
                    keys=meta['Units'].keys() # Get deys of "units" at each step... 

                    vkey = df_vars[v] # decide what the unit key will be named.
                    if instr_name_prefix is True:
                        vkey= str(Instrument)+'_'+vkey
                    
                    # Debugging line to see if what they have in header is same as var names. 
                    #print(intro_line_ripped_varname_lower, df_vars_lower[v])
                    
                    # Check to see if you have a partial string match to this var 
                    if intro_line_ripped_varname_lower in df_vars_lower[v] and intro_line_ripped_varname_lower.isnumeric() is False: 
                        if (len(line_list)> 1): #Check to see that it has a unit
                            sep = ', ' # rip the unit. 
                            unit= sep.join(line_list[1:]).strip().replace(',', ' ')
                            longname=' '.join((' '.join(line_list[1:]).strip()).split(' ')[1:])
                            
                            # Check to see if you've already saved info about this var & its units: 
                            if str(df_vars[v]) not in keys:
                                meta['Units'][vkey]= unit.strip()
                            if str(df_vars[v]) not in list(meta['Long_Name'].keys()):
                                meta['Long_Name'][vkey]= longname.strip()
                        else: # have partial string match but no units were found on line. 
                            meta['Units'][vkey]= '' # make unit empty
                            meta['Long_Name'][vkey]= '' # make unit empty
                            

                        
            if ln_num > header_row  -1:
                break  # top once you reach data.
            ln_num = ln_num + 1

    # Keep info on each dfvar of which isntrument each variable came from
    # rename columsn if appropriate. 
    for k in df_vars:
        if instr_name_prefix is True:
            meta['Vars_Instrument'][str(Instrument)+'_'+k]= Instrument
            df = df.rename(columns= {k : Instrument + '_'+ k })
        else:         
            meta['Vars_Instrument'][k]= Instrument
    
    # Fill meta of things we didn't find with blanks. 
    # list of unique instruments we got data variables from: 
    instrs = list( set( val for dic in [meta['Vars_Instrument']] for val in dic.values()))
    for h in range(0, len(meta_icartt)):
        for i in range(0, len(instrs)): 
            if instrs[i] not in meta[meta_icartt[h]].keys(): 
                meta[meta_icartt[h]][instrs[i]]= '' # make blanks. 

    return df, meta


def handle_input_configuration(DATA: dict):
    """Make sure user passed appropriate input."""
    print('1. Input ICARTT directory:' + DATA['DIR_ICARTT'])

    # 1. Ensure ICARTT directory is valid.
    if not os.path.isdir(DATA['DIR_ICARTT']):
        exit_with_error("Input ICARTT directory is invalid.")

    # 2. Ensure that directory  has ICARTT files in it. return  list.
    DATA['ICARTT_FILES'] = crawl_directory(DATA['DIR_ICARTT'],extension='.ict')

    if len(DATA['ICARTT_FILES']) == 0:
        # If no icartts were found, exit and notify user.
        exit_with_error("No ICARTT files found in the input directory.")
    else:  # Else, inform on the number of ICARTT files
        print(" - Found [ {} ] ICARTTs.\n".format(len(DATA['ICARTT_FILES'])))

    # 3. Check that a valid mode has been passed
    valid_modes = ['Stack_On_Top', 'Merge_Beside']
    if DATA['MODE'] not in valid_modes:
        exit_with_error(("Input mode entered is invalid."),
                         ("Valid Options are:" + valid_modes))

    # 4. Check that master_timeline info has been provided if necessary.
    if DATA['MODE'] == 'Merge_Beside':
        if bool(DATA.get('MSTR_TMLN')) is False:
            exit_with_error(("For input mode 'Merge_Beside', input for "),
                             ("MSTR_TMLN is also needed."))

    return DATA  # give input back with the list of icartts now included.


