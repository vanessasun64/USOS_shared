import numpy as np 
import os 
import xarray as xr
import pandas as pd

#from scipy.io import savemat
from collections import OrderedDict
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
import matplotlib as mpl
from matplotlib import colors
from matplotlib import colormaps
import matplotlib.colors as mcolors
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import cm
import matplotlib

import mplcursors

import folium

#load file for all dates
all_days_filepath = '../../CampaignData_and_Merges/R0/CSL_MobileLab_Parked/merged/rev_30min/all_CSL_MobileLab_Parked_rev30minv4.nc'
all_days_filepath_load = xr.open_dataset(all_days_filepath)
df_alldays = all_days_filepath_load.to_dataframe()
df_alldays.reset_index(inplace=True)
df_alldays.set_index('time_local', inplace=True, drop=False)

first_val_for_date = df_alldays.groupby(df_alldays.index.date).apply(lambda x: x.iloc[[0]])

map = folium.Map(location = [40.73, -111.8], tiles="OpenStreetMap", zoom_start=9)
for time_len in range(2, len(first_val_for_date)-1):
    folium.Marker(
        location = [first_val_for_date['Lat'].iloc[time_len], first_val_for_date['Lon'].iloc[time_len]],
        popup=first_val_for_date.index,
    ).add_to(map)

map