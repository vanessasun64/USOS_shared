import os 
import sys 
import re 
import yaml
import inspect 
import numpy as np 
import pandas as pd
import xarray as xr
from collections import defaultdict
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.colors import ListedColormap

from scipy.io import savemat
from collections import OrderedDict
import matplotlib.colors as mcolors

#Plot formatting
mpl.rcParams['xtick.labelsize'] = 15
mpl.rcParams['ytick.labelsize'] = 15
mpl.rcParams['legend.fontsize'] = 16
mpl.rcParams['axes.labelsize'] = 18
mpl.rcParams['axes.titlesize'] = 28
mpl.rcParams['axes.xmargin'] = 0

# Set the font family to 'serif'
mpl.rcParams['font.family'] = 'serif'
# Specify preferred serif font (Computer Modern Roman is 'cmr10')
mpl.rcParams['font.serif'] = 'Lato' 
# Optionally, configure mathtext to use Computer Modern fonts as well
mpl.rcParams['mathtext.fontset'] = 'cm'
# Ensure minus signs are rendered correctly with CM fonts
mpl.rcParams['axes.unicode_minus'] = False

# Main datapaths
dirpath = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/'
#define path for Hawthorne data directory
savepath= dirpath + 'Merge_scripts/'
plots_savepath= dirpath + 'Merge_scripts/plots/'

#Get names of variables needed from Mobile Lab and UDAQ dataset
mappings_filepath = dirpath + 'Hawthorne_data/mappings/manually_edited/UDAQ_Hawthorne_CRACMM_GEOSCHEM_CB6r5h_mapped_updated_11172025.csv'
df_mapping_parameters = pd.read_csv(mappings_filepath)
df_mapping_parameters = df_mapping_parameters.drop([0]) #drop Total NMVOCs

usos_vars = [str(spec) for spec in df_mapping_parameters['USOS Mapping'].dropna()]
usos_vars.append('time_local')
print(usos_vars)

ml_data = xr.open_dataset('/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/CampaignData_and_Merges/R0/CSL_MobileLab_Parked/merged/rev_1hr/all_CSL_MobileLab_Parked_rev1hr_iWASupdated.nc')
ml_data_vars_subset = ml_data[usos_vars]
df_ml_data = ml_data_vars_subset.to_dataframe()
df_ml_data.set_index('time_local', inplace=True)

udaq_voc_file = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/Hawthorne_data/data/hawthorne_udaq_all_vocs_hourly_timezone_carbon_number_updated.csv'
df_udaq_data = pd.read_csv(udaq_voc_file, index_col='time_local', parse_dates=True)

#Get only the data during the overlapping times
df_ml_data = df_ml_data.sort_index().loc['2024-07-15 12:00:00':'2024-08-18 18:00:00']
df_udaq_data = df_udaq_data.sort_index().loc['2024-07-15 12:00:00':'2024-08-18 18:00:00']

#Map the same species in ML and UDAQ data
mapping = dict(zip(df_mapping_parameters['USOS Mapping'], df_mapping_parameters['UDAQ_Variable']))
#print(mapping)

#remove Formaldehyde from VOC species
#print(df_ml_data)
df_ml_data = df_ml_data.drop(columns = 'HCHO_CRDS')
df_udaq_data = df_udaq_data.drop(columns ='Formaldehyde')
del mapping['HCHO_CRDS']
print(df_ml_data.columns)
print(df_udaq_data.columns)
print(mapping)

#Turn any negatives into NaNs
for spec in df_ml_data.columns:
    df_ml_data[spec] = df_ml_data[spec].mask(df_ml_data[spec] < 0, np.nan)

#Rename columns to a more common name by using the mappings
df_ml_data.columns = [
    mapping.get(col) if not pd.isna(col) else mapping.get(np.nan, col)
    for col in df_ml_data.columns
]

# Fraction of each species per hour per day
hourly_fraction = df_ml_data.div(df_ml_data.sum(axis=1), axis=0)
print(hourly_fraction)
# Average fractions across all days grouped by hour
mean_per_hour = hourly_fraction.groupby(hourly_fraction.index.hour).mean()

print(mean_per_hour)
# # Select the typical time fractions
fractions_9 = mean_per_hour.loc[9]
fractions_12 = mean_per_hour.loc[12]
fractions_15 = mean_per_hour.loc[15]
fractions_18 = mean_per_hour.loc[18]
# print(fractions_9)

ml_data = xr.open_dataset('/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/CampaignData_and_Merges/R0/CSL_MobileLab_Parked/merged/rev_1hr/all_CSL_MobileLab_Parked_rev1hr_iWASupdated.nc')
ml_data_vars_subset = ml_data[usos_vars]
df_ml_data = ml_data_vars_subset.to_dataframe()
df_ml_data.set_index('time_local', inplace=True)

udaq_voc_file = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/Hawthorne_data/data/hawthorne_udaq_all_vocs_hourly_timezone_carbon_number_updated.csv'
df_udaq_data = pd.read_csv(udaq_voc_file, index_col='time_local', parse_dates=True)

#Get only the data during the overlapping times
df_ml_data = df_ml_data.sort_index().loc['2024-07-15 12:00:00':'2024-08-18 18:00:00']
df_udaq_data = df_udaq_data.sort_index().loc['2024-07-15 12:00:00':'2024-08-18 18:00:00']

#Map the same species in ML and UDAQ data
mapping = dict(zip(df_mapping_parameters['USOS Mapping'], df_mapping_parameters['UDAQ_Variable']))
#print(mapping)

#Turn any negatives into NaNs
for spec in df_ml_data.columns:
    df_ml_data[spec] = df_ml_data[spec].mask(df_ml_data[spec] < 0, np.nan)

#Rename columns to a more common name by using the mappings
df_ml_data.columns = [
    mapping.get(col) if not pd.isna(col) else mapping.get(np.nan, col)
    for col in df_ml_data.columns
]
# Step 1: Resample to hourly sums
hourly_sum = df_ml_data.resample('h').sum()

# Step 2: Compute fractional contribution per hour
hourly_fraction = hourly_sum.div(hourly_sum.sum(axis=1), axis=0) * 100  # Percent

# # Step 3: Group by hour of the day and take the mean
hourly_fraction['hour'] = hourly_fraction.index.hour
mean_fraction_by_hour = hourly_fraction.groupby('hour').mean()
# mean_fraction_by_hour.drop(columns='hour', inplace=True)

# Get the contributions for 12 PM (or any hour)
contributions = mean_fraction_by_hour.loc[12]
# Sort from largest to smallest
contributions_sorted = contributions.sort_values(ascending=False)

# Now mean_fraction_by_hour contains the average fractional contribution for each hour (0-23)
# Set threshold
threshold = 3

# Create labels: only show car names if contribution >= threshold, else empty string
labels = [car if val >= threshold else '' for car, val in zip(contributions_sorted.index, contributions_sorted)]

# Custom autopct function
def autopct_threshold(pct, threshold=3):
    return f'{pct:.1f}%' if pct >= threshold else ''  # Show label only if >= threshold

plt.figure(figsize=(6,6))
plt.pie(
    contributions_sorted,
    labels=labels,
    autopct=lambda pct: autopct_threshold(pct, threshold=3),
    startangle=90
)

# Add a single 'Other' annotation
plt.annotate(
    'Other', 
    xy=(0, 0),  # Center of the pie chart
    xytext=(0.9, 0.8),  # Position text below the pie
    ha='center', 
    fontsize=15
)
plt.title('12 PM')
plt.show()

# Step 1: Resample to hourly sums
hourly_sum = df_ml_data.resample('h').sum()

# Step 2: Compute fractional contribution per hour
hourly_fraction = hourly_sum.div(hourly_sum.sum(axis=1), axis=0) * 100  # Percent

# # Step 3: Group by hour of the day and take the mean
hourly_fraction['hour'] = hourly_fraction.index.hour
mean_fraction_by_hour = hourly_fraction.groupby('hour').mean()
# mean_fraction_by_hour.drop(columns='hour', inplace=True)

hours_to_plot = [6, 9, 12, 15, 18, 21]
times_for_plot = ['6 AM', '9 AM', '12 PM', '3 PM', '6 PM', '9 PM']
df_plot = mean_fraction_by_hour.loc[hours_to_plot]

print(df_plot)
# # Sort from largest to smallest
# df_plot_sorted = df_plot.sort_values(by=df_plot.columns, ascending=False)
# print(df_plot_sorted)

#print(df_plot)
series_list = []
for idx, row in df_plot.iterrows():
    row_sorted = row.sort_values(ascending=False)
    display_labels = []
    other_series = pd.Series([0], index = ['Other'])
    row_sorted_with_other = pd.concat([row_sorted, other_series])
    #print(row_sorted_with_other) #This gives us a pandas series with index as the species and value as its fractional contribution, one for each hour at 6 AM, 9 AM, etc in 3 hour intervals
    
    for val, col in zip(row_sorted_with_other.values, row_sorted_with_other.index):
        #Go through all the values in these dataseries
        #If it is less than the threshold and not in the other category, add the label to display_labels
        if val >= threshold and col != 'Other' :
            display_labels.append(col)
            #total_percentages_excluding_other.append(val)
        #If the value is less than 3% contribution to the total VOCs, then add sum it into the Other category, replace the Other value with it, and remove the species from the Series index
        elif val < threshold and col != 'Other':
            spec_plus_other = row_sorted_with_other.loc[col] + row_sorted_with_other.loc['Other']
            #print(spec_plus_other)
            row_sorted_with_other = row_sorted_with_other.replace(row_sorted_with_other.loc['Other'], spec_plus_other)
            row_sorted_with_other = row_sorted_with_other.drop(col)
        else:
            display_labels.append('')  # hide small slice
        #This loop will complete for Series 1, which represents the full series for 6 AM then go on to altering 9 AM, then 12 PM.
    
    #print(row_sorted_with_other)
    series_list.append(row_sorted_with_other)
#print(series_list)
all_species = pd.Index([])

for s in series_list:
    all_species = all_species.union(s.index)

#print(all_species)

# hues = np.linspace(0, 1, 11, endpoint=False)
# colors_list = [mcolors.hsv_to_rgb([h, 0.3, 0.95]) for h in hues]
cmap = plt.get_cmap("tab10")
colors_list = [cmap(i) for i in range(11)]

# Map each column to a color
colors_dict = {col: color for col, color in zip(all_species, colors_list)}
print(colors_dict)

fig, axes = plt.subplots(2,3, figsize=(16, 8), constrained_layout=False)

# Adjust spacing between subplots
plt.subplots_adjust(
    left=0.05,    # space on left of figure
    right=0.97,   # space on right of figure
    top=0.95,     # space at top
    bottom=0.001,  # space at bottom
    wspace=0.1,   # horizontal space between subplots
    hspace=0.04    # vertical space between subplots
)

axes = axes.flatten()
threshold = 3  # percentage threshold

# hues = np.linspace(0, 1, 11, endpoint=False)
# colors_list = [mcolors.hsv_to_rgb([h, 0.4, 0.9]) for h in hues]

# # Map each column to a color
# colors_dict = {col: color for col, color in zip(mean_fraction_by_hour.columns, colors_list)}


# for idx, row in df_plot.iterrows():
#     row_sorted = row.sort_values(ascending=False)
#     display_labels = []
#     other_series = pd.Series([0], index = ['Other'])
#     row_sorted_with_other = pd.concat([row_sorted, other_series])

#     for val, col in zip(row_sorted_with_other.values, row_sorted_with_other.index):
#         if val >= threshold and col != 'Other' :
#             display_labels.append(col)
#             #total_percentages_excluding_other.append(val)
#         elif val < threshold and col != 'Other':
#             spec_plus_other = row_sorted_with_other.loc[col] + row_sorted_with_other.loc['Other']
#             #print(spec_plus_other)
#             row_sorted_with_other = row_sorted_with_other.replace(row_sorted_with_other.loc['Other'], spec_plus_other)
#             row_sorted_with_other = row_sorted_with_other.drop(col)
#         else:
#             display_labels.append('')  # hide small slice


for ax, (idx, row), label in zip(axes, df_plot.iterrows(),  times_for_plot):
    # # Predefine a color for each column
    # cmap = plt.get_cmap('tab20')  # or any other colormap
    # colors_dict = {col: cmap(i / len(df_plot.columns)) for i, col in enumerate(df_plot.columns)}
    row_sorted = row.sort_values(ascending=False)
    #print(row_sorted)

 # Prepare labels: hide if percentage < threshold
    display_labels = []
    # total_percentages_excluding_other = []
    #print(row_sorted.index)
    other_series = pd.Series([0], index = ['Other'])
    row_sorted_with_other = pd.concat([row_sorted, other_series])
    #print(row_sorted_with_other.index)
    for val, col in zip(row_sorted_with_other.values, row_sorted_with_other.index):
        if val >= threshold and col != 'Other' :
            display_labels.append(col)
            #total_percentages_excluding_other.append(val)
        elif val < threshold and col != 'Other':
            spec_plus_other = row_sorted_with_other.loc[col] + row_sorted_with_other.loc['Other']
            #print(spec_plus_other)
            row_sorted_with_other = row_sorted_with_other.replace(row_sorted_with_other.loc['Other'], spec_plus_other)
            row_sorted_with_other = row_sorted_with_other.drop(col)
        else:
            display_labels.append('Other')  # hide small slice

    print(display_labels)
            
    # total_percentages_excluding_other = np.array(total_percentages_excluding_other)
    # sum_percentages_excluding_other = np.sum(total_percentages_excluding_other)
    # other_percentage = round(100-sum_percentages_excluding_other, ndigits=1)
    # print(other_percentage)
    
    # Prepare colors: match the column name
    row_colors = [colors_dict[col] for col in row_sorted_with_other.index]
    
    # autopct function: hide percentages < threshold
    autopct = lambda pct: f'{pct:.1f}%' if pct >= threshold else ''

    wedges, texts, autotexts = ax.pie(
        row_sorted_with_other,
        labels=display_labels,
        colors=row_colors,
        autopct=autopct,
        startangle=90
    )

    # Outer labels (display_labels)
    for text in texts:
        text.set_fontsize(14)  # adjust outer labels size

    # Percentages inside slices
    for autotext in autotexts:
        autotext.set_fontsize(10)  # adjust percentages size
    
    ax.set_title(label, y=0.9, fontsize = 24)

    # ax.annotate(
    #     'Other', 
    #     xy=(0, 0),  # Center of the pie chart
    #     xytext=(0.48, 0.97),  # Position text around the pie
    #     ha='center', 
    #     fontsize=14
    # )
    # ax.annotate(
    #     str(other_percentage) + '%', 
    #     xy=(0, 0),  # Center of the pie chart
    #     xytext=(0.4, 0.6),  # Position text around the pie
    #     ha='center', 
    #     fontsize=10
    # )



# # Hide unused axes
# for ax in axes[n_charts:]:
#     ax.axis('off')

plt.savefig(dirpath + '/UDAQ_reports/q1_report/plots/mobilelab_voc_speciation.png', dpi =300)
plt.show()

