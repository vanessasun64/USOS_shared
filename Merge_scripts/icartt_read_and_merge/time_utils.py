#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 24 15:04:06 2023

@author: u6044586
"""

import re
import datetime
import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
import warnings

import ict_utils # import custom submodule ict_utils 


def align2master_timeline(df: pd.DataFrame, startdt: str, enddt: str,
                          step_S: int, quiet: bool = False,
                          lim: int = None, datetime_index: bool = False,
                          tzf: str = 'UTC', use_TZ: bool = True, 
                          ):
    """Resample dataframes to appropriate timelines."""
    # Function to take a dataframe and appropriately remap it to a new time
    # index, considering the native sampling frequency as it is relative
    # to the desired new time. Writte 2/6/21, jessica d. haskins

    # --------------------------Inputs:-------------------------------------
    # df - dataframe, mustcontain column 'datetime'.
    # startdt, enddt = '2006-03-01 00:00:00', '2006-04-01 00:00:00'
    # step_S= 120 averaging step in seconds (120 for a 2 minute average).
    # quiet - Set to False to show sanity check on averaging.
    # lim -manually set the limit of # of points to include in an avg.
    # datetime index - Set to True if datetime is already an index of the df.
    
    if (datetime_index is False) and 'datetime' not in df:
        ict_utils.exit_with_error(("Dataframe passed to align2master_timeline()",
                         "does not contain a column 'datetime'. "))

    # Make datetime an index and remove duplicates.
    if datetime_index is False:
        df = df.dropna(subset=['datetime'])  # if datetime is nan drop whole row.
        df.set_index('datetime', inplace=True)  # Make the datetime an index.

    df.sort_index(inplace=True) # Sort all columns by datetime index. 
    df = df[~df.index.duplicated()]  # remove any duplicates rows

    # Get the average native sampling frequency in total seconds:
    tseries = df.index.to_series()
    min_sep = int(np.round(tseries.diff().median().total_seconds()))

    if quiet is False:
        print('Native Mean Time Sep. (s): ', str(min_sep) + 's')

    # If the native time seperation is less than X seconds, you'll
    # reindex it to our full date range, as close to native  freq as you can,
    # then take a roliing avg to get the X second avg on the time base we want
    if min_sep < step_S:
        if use_TZ is True: 
            dts = pd.date_range(startdt, enddt, freq=str(min_sep) + 's', tz=tzf)
        else: 
            dts = pd.date_range(startdt, enddt, freq=str(min_sep) + 's') # just do native time!

        if not lim:
            lim = np.round(4350 / min_sep)  # don't fill if collect>than 1H out
        

        dfn = df.reindex(dts, method='nearest', fill_value=np.nan,limit=int(lim))
       
        # Take a centered boxcar average around the X s avg. (for numerical columns only)!    
        #NOTE: .mean() handles Nans like np.nanmean() in this context!!! 
        df_nums=dfn.select_dtypes(exclude=['object'])
        df_nums_new = df_nums.rolling(str(int(step_S)) + 's').mean().resample(str(step_S) + 's').mean()
        dtss=df_nums_new.index
        df_nonums=dfn.select_dtypes(include=['object'])
        df_nonums_new = df_nonums.reindex(dtss, method='nearest', fill_value=np.nan)
        
        df_new=pd.concat([df_nums_new, df_nonums_new], axis=1, join="inner")
                               
        
    else:
        # The native sampling frequency is longer than X seconds, so just pull
        # the closest values along to fill our array.
        print(('WARNING: You have input an averaging frequency that is LESS'),
              ('than this instruments average native sampling frequency.'),
              ('This is dangerous and can lead to errors because it is NOT'),
              ('interpolating data, but rather "filling from nearest".'),
              ('Consider raising your input averaging time step.'))
        if use_TZ is True: 
            dts = pd.date_range(startdt, enddt, freq=str(step_S) + 's', tz=tzf)
        else: 
            dts = pd.date_range(startdt, enddt, freq=str(step_S) + 's')
        
        if not lim:
            # Don't fill if collected > than 1H 15 mins out from here.
            # If lim not set and too small, just use 1 point.
            lim = max(1, np.round(4350 / min_sep))
       
        df_new = df.reindex(dts, method='nearest', fill_value=np.nan ,limit=int(lim))


    # Plot the Original Data & the Re- Mapped stuff so you can see if its good:
    if quiet is False:
        one = df.columns[0]
        fig, ax = plt.subplots(1, figsize=(8, 8))
        ax.plot(df[one], label="original",
                color='orange', linewidth=1)
        ax.plot(df_new[one],
                label="re-mapped", color='blue', linewidth=1)
        # ax.set_xlim(df.index[200], df.index[400])
        # ax.set_ylim(34, 39)
        ax.legend()
        plt.show()

    return df_new


def find_datelike_cols(df: pd.DataFrame, icartt_file: str,
                        quiet: bool = True, more_names: list = list()):
    """Identify the date like columns in a dataframe from an icartt file."""
    # List of partial Strings to search columnNames for that will ID them as
    # time cols. ***NOTE:  If you are getting a persistent error about this
    # function, try adding to the list of tm_names. E.g. add other timezone
    # denotations. Setting quiet to False will print your col names and
    # what it is finding so you can decide what to add to this list.
    tm_names = ['utc', 'cst', 'cdt', 'local', 'lst', 'est', 'pst','mst',
                'gmt', 'time_mid', 'central_standard',
                'eastern_standard', 'pacific_standard'] +more_names

    print('All Column Names:', df.columns) if quiet is False else None

    # Identify all the names of time related columns in the dataframe.
    times = list()  # Empty list to contain columns with time-like names.

    #Edited June 2025 by Jessica Haskins, for excluding humidity as a time: was reading any variable with "mid" as a time.
    for col in df.columns:
        if any(nm in col.lower() for nm in tm_names if (('humidity' not in col.lower()) and (nm !='mid'))):
            
            times.append(col.lower())  # fill the list with those names.
            # Rename all time cols in lowercase to make string matches easier
            df.rename(columns={col: col.lower()}, inplace=True)
    
    # Make sure you haven't grabbed a day column for time by accident. Drop it.
    for h in times:
        times.remove(h) if 'day' in h else None

    print('Original Time Columns Found:', times) if quiet is False else None
    # Return the dataframe with lowercase time names, a list of the time
    # columns found, and a list of time strings you searched for in the col
    # names (this is just for consistency when looping...)
    return df, times, tm_names


def make_time_midpoint_cols(df: pd.DataFrame, tm_names: list, times: list,
                             quiet: bool = True):
    """Take start/stop times from datelike cols and turn them into midpts."""
    # Get start/stop pairs for time to make a midpoint if you can find both.

    # Create a dict for time names & their ID'd "type" for scanning later.
    nn_times = {}
    for j in range(0, len(tm_names)):  # Check each time zone name
        start_j = None  # reset on each loop through a dif time_nm
        stop_j = None
        has_tm = None
        for i in range(0, len(times)):  # Check each time col name.
            # Assign start/stop pair variables.
            if (tm_names[j] in times[i]) and ('start' in times[i]):
                start_j = times[i]

            if (tm_names[j] in times[i]) and \
               ('stop' in times[i] or 'end' in times[i]):
                stop_j = times[i]

            # Or just let us know if this "time" col exists.
            if (tm_names[j] in times[i]):
                has_tm = times[i]

        # You found a start/stop pair for this time name.
        if (start_j is not None) and (stop_j is not None):
            # Make a new column in the data frame that is the midpoint.
            df[tm_names[j] + '_mid'] = (df[start_j] + df[stop_j]) / 2

            # Remove the start & stop pair of time from the larger df
            df = df.drop(columns=[start_j, stop_j])

            # Update list of column names with times to reflect above.
            times.append(tm_names[j] + '_mid')  # add midpoint name
            times.remove(start_j)  # remove oldies
            times.remove(stop_j)

        elif has_tm is not None:
            # Populate dict associating its actual timename and its "type" so
            # we can choose between times to use as master later if we want.
            nn_times[has_tm] = 'time_' + tm_names[j]

            # Didn't find a start/stop pair but do have a col with this tm nm.
            # df.rename(columns={has_tm: 'time_' + tm_names[j]},
            #          inplace=True)  # rename it so we know what it is called

            # Update list of column names with times to reflect above.
            # times.append('time_' + tm_names[j])  # add new name
            # times.remove(has_tm)  # remove old name

    print('Time cols After Mid_Point Assign:',
          times) if quiet is False else None

    return df, times, tm_names, nn_times


def pick_a_single_time_col(times: list, nn_times: dict, quiet: bool = True, tz_default:str=''):
    """Decide which time column you prefer to use a indx if you got lots."""
   
    # Create dictionary for "picking" which time variable we prefer to use.
    # Change the rank for preferences here. Lower = more desirable.
    zones=['utc','est','edt','cst','cdt','mst','mdt','pst','pdt']
    pref_time=dict({}); ct=1 
    for z in zones: 
        pref_time['time_'+z]=ct
        pref_time[z+'_mid']=ct+1 
        ct=ct+2 
    others=['time_mid','time_time','time_time_mid','time_mid_mid','time_start','time_mid',
            'time_igortime','igortime','igor_time','time_local','local_mid','time_lt',
            'lt_mid','time_lst','lst_mid']
    for o in others: pref_time[o]= ct; ct=ct+1 

    # Create dictionary to associate a "time" column name with its timezone.
    tz_info = {'time_utc': 'UTC', 'utc_mid': 'UTC',
               'time_est': 'EST5EDT', 'est_mid': 'EST5EDT',
               'time_edt': 'EST5EDT', 'edt_mid': 'EST5EDT',
               'time_cst': 'CST6CDT', 'cst_mid': 'CST6CDT',
               'time_cdt': 'CST6CDT', 'cdt_mid': 'CST6CDT',
               'time_mst': 'MST7MDT', 'mst_mid': 'MST7MDT',
               'time_mdt': 'MST7MDT', 'mdt_mid': 'MST7MDT',
               'time_pdt': 'PST8PDT', 'pdt_mid': 'PST8PDT',
               'time_pst': 'PST8PDT', 'pst_mid': 'PST8PDT'}

    pref_arr = list()  # empty list to contain pref rank of "time cols"
    tz_arr = list()  # Empty list to contain the timezone of the "time cols"
    
    for n in range(0, len(times)):  # Get an array of # preferences for "times"
        # Get the nickname of the time col if its not a mid point.
        check_if_NOT_mid = nn_times.get(times[n], None)
        if check_if_NOT_mid is None:  # has no nickname,it is a midpoint.
            print('Found time cols:', times[n], 'Rank of time pref is:', pref_time.get(times[n]), \
                  'Timezone is : ', tz_info.get(times[n])) if quiet is False else None
            pref_arr.append(pref_time.get(times[n], 100))
            tz_arr.append(tz_info.get(times[n], 100))
        else:  # has a nickname, pass that to get preference.
            print('Found time cols:', times[n], 'has nickname: ',nn_times.get(times[n]),\
                  'Rank of time pref is:', pref_time.get(check_if_NOT_mid), \
                  'Timezone is: ', tz_info.get(times[n])) if quiet is False else None
            pref_arr.append(pref_time.get(check_if_NOT_mid, 100))
            tz_arr.append(tz_info.get(check_if_NOT_mid, 100))

    if min(pref_arr) == 100:
        ict_utils.exit_with_error('Time columnName in the "pick_a_single_time_col()" '+\
                         'function could not be properly identified. Try '+\
                         'running the call to function with quiet=False '+\
                         'to see debug output. ')

    # Pick one of the time columns based on your ranked preferences.
    time_pref = times[(pref_arr == min(pref_arr))]  # Name of the pref time col
    tz_pref = tz_arr[(pref_arr == min(pref_arr))]  # Time zone of pre time col

    if tz_pref == 100: # esp if using prefixes, sometimes won't catch timezone using dict. 
        # So if tz not found, then search time names for common timezone indicators. 
        check_for_tz_str=list() # empty list to check time pref names for timezone info.
        check_for_tz_str= [t for t in ['est', 'cst', 'cdt', 'utc', 'pst', 'pdt', 'mst', 'mdt'] if t in time_pref]
        if len(check_for_tz_str)> 0:
            if check_for_tz_str[0] =='est': tz_pref= 'EST5EDT'
            if (check_for_tz_str[0] =='cst') or (check_for_tz_str[0] =='cdt'): tz_pref= 'CST6CDT'
            if (check_for_tz_str[0] =='mst') or (check_for_tz_str[0] =='mdt'): tz_pref= 'MST7MDT'
            if (check_for_tz_str[0] =='pst') or (check_for_tz_str[0] =='pdt'): tz_pref= 'PST8PDT'
            if check_for_tz_str[0] =='utc': tz_pref= 'UTC'  

    if tz_pref == 100:
        if tz_default=='':
            print('The timezone of the preferred time array is ambiguious & could not be identified. '+ 
                      'Please input the timezone abbreviation we should use for this time array.'+ 
                      'Valid options include: UTC, EST5EDT, CST6CDT, MST7MDT, PST8PDT, etc.')
            tz_pref = input ('... ')
            print('You entered: ', tz_pref)
    
            if tz_pref == 100: # if tz pref still undefined... 
                ict_utils.exit_with_error('Still has no time zone. Exiting.')      
        else: 
            print('The timezone of the preferred time array is ambiguious & could not be identified. '+\
                  'Using inputted preffered timezone: '+tz_default)
            tz_pref=tz_default
            
    print('Pref time col:', time_pref, ' in', tz_pref,
          'Timezone') if quiet is False else None

    bad_times = times  # duplicate list of all times
    bad_times.remove(time_pref)  # and drop all non pref times.

    # Return name of preffered column its timezone and the names of all the
    # columns that you don't want to use.
    return time_pref, tz_pref, bad_times


def icartt_time_to_datetime(df: pd.DataFrame, yr, mon, day, time_col: str,
                            tz_pref: str, remove_old_time: bool = True):
    """Convert seconds since midnight (icartt time col) 2 datetime obj col."""
    # Takes yr mon day in string or int form. Get to ints if strings passed.
    yr = int(yr) if type(yr) == str else yr
    mon = int(mon) if type(mon) == str else mon
    day = int(day) if type(day) == str else day

    # Add "timedelta" of seconds since midnight to the date the icarrt
    # file started on (typically in the file name).
    datetime_col_i = datetime.datetime(year=yr, month=mon, day=day) + \
        pd.to_timedelta(df[time_col], unit='s')
            
    # Teach this tz unaware pandas series its native timezone.
    datetime_col = datetime_col_i.dt.tz_localize(tz_pref,nonexistent='shift_forward')  # now it has a tz

    # Convert from native tz to UTC time.
    datetimecol_in_UTC = datetime_col.dt.tz_convert('UTC')

    df['datetime'] = datetimecol_in_UTC  # Create new column with date in UTC.

    if remove_old_time is True:  # Drop the old col from the df if asked
        df = df.drop(columns=time_col)

    return df


def master_icartt_time_parser(df: pd.DataFrame, icartt_file: str,
                               quiet: bool = True, remove_old_time:
                                   bool = True, tz_default:str=''):
    """Identify the date like columns, convert to TZ aware datetimes)."""

    # Identify the time-like columns in this dataframe:
    df, times, tm_names = find_datelike_cols(df, icartt_file, quiet)
    
    if len(times)==0: # try again with more names for "time" (will pick up too many if try first time.)
        df, times, tm_names = find_datelike_cols(df, icartt_file, quiet, more_names =['Time','time', 'StartTime', 'mid', 'start', 'stop'])
         
    # Take Start/Stop pairs of time cols and convert them to midpoint times.
    df, times, tm_names, nn_times = make_time_midpoint_cols(df, tm_names,
                                                             times, quiet)
    
    # Fix "IgorTime" if you got it in this Icartt, then just overwrite var in df 
    igor_times= [x for x in times if 'igor' in x.lower()]
    for bt in range(0, len(igor_times)):  # needs to be in seconds since midnight... 
        df[igor_times[bt]]= pd.to_datetime(df[igor_times[bt]], origin="1904-01-01", unit='s') - datetime.datetime(1904,1,1)

    # Pick which time var to use and which time zone its in.
    time_pref, tz_pref, bad_times = pick_a_single_time_col(times,nn_times, quiet=quiet, tz_default=tz_default)
                        
    if remove_old_time is True:  # Remove the non-preferred times
        df = df.drop(columns=bad_times)
        # you've dropped all other names from the df.
        times = list([time_pref])
    else:
        # you haven't dropped anything.
        times = list([time_pref, bad_times])

    # Get the date this data was collected on from the icarttfile name passed.
    date_full = re.search(r'\d{4}\d{2}\d{2}', icartt_file).group(0)
    yr = date_full[0:4]
    mm = date_full[4:6]
    dd = date_full[6:8]
    
    # Tell the people which variable was chosen as "datetime".
    print('The time variable chosen to be converted to "datetime" is:',
          time_pref+' assumed to be in '+tz_pref)
    
    # Convert the preffered time column to a column named 'datetime' and drop
    # all the other time columns from the larger dataframe.
    df = icartt_time_to_datetime(
        df, yr, mm, dd, time_pref, tz_pref, remove_old_time)

    # And update the times list to remove the old time and add "datetime"
    times.append('datetime')
    if remove_old_time is True:
        times = times.remove(time_pref)

    return df, times
