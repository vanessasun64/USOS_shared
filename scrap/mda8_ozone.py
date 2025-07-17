# Filter the DataFrame to only include summer dates between June 1st and August 1st with O3 < 120ppb

#data had hourly data
df = df[(df.index.month >= 6) & (df.index.month <= 8) & (df['Ox']< 100)]
 
# Step 1: Calculate the 8-hour rolling average
df['O3_8hr_avg'] = df['O3'].rolling(window=8, min_periods=8).mean()
#rolling average only works for 8 hours with no gaps; at least 6 hours
 
# Step 2: Calculate the maximum 8-hour average for each day
# Resample to daily frequency, compute the daily max of the rolling averages, drop NA values
daily_max_8hr_avg = df['O3_8hr_avg'].resample('D').max().dropna()
 
# Step 3: Map the daily maximum back to the original dataframe
# Create a new temporary column with the daily max 8-hour average for each timestamp
df['MD8A_O3'] = df.index.floor('D').map(daily_max_8hr_avg)
 
# Select daytime values only where MD8A > 65
df_day = df[(df.index.hour >=7) & (df.index.hour<=20) & (df['MD8A_O3']>65)]