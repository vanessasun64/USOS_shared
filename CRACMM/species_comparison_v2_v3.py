#This program finds all the new species in CRACMM3M (expanded halogen chemistry) that were not in CRACMM2

import pandas as pd

dirpath = '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/'
metadata_dirpath = dirpath+'CRACMM/CRACMM_20251014update/metadata/'

pd.set_option('display.max_rows', None)

cracmm2_load = pd.read_csv(metadata_dirpath+'cracmm2/cracmm2_metadata.csv')
cracmm3_load = pd.read_csv(metadata_dirpath+'cracmm3m/cracmm3m_metadata.csv')

#select columns from each dataframe
v2_cols = cracmm2_load[["Species", "Description"]]
v3_cols = cracmm3_load[["Species", "Description"]]

# # Combine them into a new dataframe
# combined_df = pd.concat([v2_cols, v3_cols], axis=1)
# display(combined_df)

# Add a source label and preserve original index
v2_labeled = v2_cols.assign(source='V2', original_index=v2_cols.index)
v3_labeled = v3_cols.assign(source='V3', original_index=v3_cols.index)

# Combine only the identifying columns
combined = pd.concat([v2_labeled, v3_labeled], ignore_index=True)

# Count occurrences of each species across both dataframes
type_counts = combined['Species'].value_counts()

# Keep only those species that appear in only one dataframe
unique_types = type_counts[type_counts == 1].index

# Filter combined dataframe to keep only unique species
result = combined[combined['Species'].isin(unique_types)].reset_index(drop=True)
display(result)