#%% 0. Start
import pandas as pd


#%% 1. Open data
#%%% 1.1. Blocking catalogue
blocks_daily = pd.read_csv('../Data/Output_data/03-Blocking_daily_catalogue_2019_2020_NH.csv')
blocks_event = pd.read_csv('../Data/Output_data/03-Blocking_event_catalogue_2019_2020_NH.csv')


#%% 3. Filter by season
winter_All = blocks_event[blocks_event.MONTH_START.isin([12,1,2])]
spring_All = blocks_event[blocks_event.MONTH_START.isin([3,4,5])]
summer_All = blocks_event[blocks_event.MONTH_START.isin([6,7,8])]
fall_All = blocks_event[blocks_event.MONTH_START.isin([9,10,11])]


#%% 4. Numbers for table
print(f'All events (winter)             |     nmr: {len(winter_All)}  |  duration: {round(winter_All.DURATION.mean(),2)}')
print(f'All events (spring)             |     nmr: {len(spring_All)}  |  duration: {round(spring_All.DURATION.mean(),2)}')
print(f'All events (summer)             |     nmr: {len(summer_All)}  |  duration: {round(summer_All.DURATION.mean(),2)}')
print(f'All events (fall)               |     nmr: {len(fall_All)}  |  duration: {round(fall_All.DURATION.mean(),2)}')
print(f'All events (total)              |     nmr: {len(blocks_event)} |  duration: {round(blocks_event.DURATION.mean(),2)}')