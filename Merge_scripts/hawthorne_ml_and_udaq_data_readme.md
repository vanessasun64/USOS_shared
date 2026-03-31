# All measurements together
* Has datetime index 2024-07-14 18:00:00 to 2024-08-18 17:45:00 local time (MDT) with index called 'time_local'

`/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/Merge_scripts/gap_filling/all_measured_species.csv`
    * Notes: 
        * UDAQ mEthyltoluene and pEthyltoluene replaced for combined 'mpEthyltoluene' variable


* Filled with duplicate measurements and tracer species, interpolation:
`/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/CampaignData_and_Merges/R0/CSL_MobileLab_Parked/merged/rev_15min/all_CSL_MobileLab_Parked_rev15min_iWASupdated_filled.nc`
    * ** NOTE: In UTC **

# NOAA ML measurements
* merged & has metadata including iWAS data indexed to 15 minutes, has holes: 
`/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/CampaignData_and_Merges/R0/CSL_MobileLab_Parked/merged/rev_15min/all_CSL_MobileLab_Parked_rev15min_iWASupdated.nc`

# UDAQ measurements:

## Most VOCs
* All have datetime index 2024-07-14 18:00:00 to 2024-08-18 17:45:00 local time (MDT) with index called 'time_local'

* Original UDAQ GC data is in ppb of carbon & 1 hour off (in MST instead of MDT), converted to ppb & shifted index to MDT using `/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/Hawthorne_data/scripts/reformat_data_scripts.py`, indexed to 15 minutes: `/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/Hawthorne_data/data/script_output/hawthorne_udaq_all_vocs_15min_timezone_carbon_number_updated.csv`
    * Index is in column time_local; pandas might read time_local as UTC.  I read in the file as:
        ```
        df = pd.read_csv(filename, index_col='time_local', parse_dates=True)
        df.index = df.index.tz_localize(None)
        df.index = df.index.tz_localize('America/Denver')
        ```

## Terpenes (Isoprene, Alpha & Beta pinene)
* Renamed variables/sliced from raw files in
`/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/Hawthorne_data/data/modified_from_Bart_raw/terpenes_HW_20240101_20241231.csv`, then shifted time from MST to MDT (1 hour) in `/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/Hawthorne_data/data/script_output/hawthorne_udaq_isoprene_alpha_beta_pinene_07152024_08012024_15min_reindexed_timezone_updated.csv`
    * Index is in column time_local; pandas might read time_local as UTC.  I read in the file as:
        ```
        df = pd.read_csv(filename, index_col='time_local', parse_dates=True)
        df.index = df.index.tz_localize(None)
        df.index = df.index.tz_localize('America/Denver')
        ```

## Formaldehyde
* raw data in `/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/Hawthorne_data/from_Bart/hw_zero_corrected_data_formaldehyde.csv`
    * **NOTE: Data in UTC (unlike most other UDAQ measurements)**
* averaged for 15 minutes index in MDT (local time) by script in `/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/Hawthorne_data/scripts/reformat_data_scripts.py`, data saved in:
`/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/Hawthorne_data/data/script_output/hawthorne_udaq_Formaldehyde_15min_reindexed_timezone_updated.csv`
    * Index is in column time_local; pandas might read time_local as UTC.  I read in the file as:
        ```
        df = pd.read_csv(filename, index_col='time_local', parse_dates=True)
        df.index = df.index.tz_localize(None)
        df.index = df.index.tz_localize('America/Denver')  


## Ozone
* raw data: `/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/Hawthorne_data/data/from_Bart/hawthorne_udaq_o3_2024.csv`
* reformatted date format, shifted time zone by 1 hour from MST original to MDT and localized time zone, reindexed to 15 minutes, removed values over 120 ppb, and converted from ppm to ppb in script `/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/Hawthorne_data/scripts/reformat_data_scripts.py`, data saved in:
`/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/Hawthorne_data/data/script_output/hawthorne_udaq_o3_2024_15min_reindexed_timezone_updated.csv`

## NO, NOx, and NO2
* raw data: `/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/Hawthorne_data/data/from_Bart/udaq_QC_07012024_08012024.csv`
* shifted time zone by 1 hour from MST original to MDT and localized time zone, reindexed to 15 minutes by script: `/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/Hawthorne_data/scripts/reformat_data_scripts.py`
    Data saved in:
    * `/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/Hawthorne_data/data/script_output/hawthorne_udaq_no_07152024_08012024_15min_reindexed_timezone_updated.csv`
    * `/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/Hawthorne_data/data/script_output/hawthorne_udaq_no2_07152024_08012024_15min_reindexed_timezone_updated.csv`
    * `/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/Hawthorne_data/data/script_output/hawthorne_udaq_noy_07152024_08012024_15min_reindexed_timezone_updated.csv`
