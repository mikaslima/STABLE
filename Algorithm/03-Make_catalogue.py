##################################################################
#
#
# Script that computes several statistics of each daily event and overall events during their lifetimes
# Takes Z500 and tracked continuous structure information as input
# Outputs a netcdf with struct masks with their respective id and two csv files with their statistics per day and per event
#
#
##################################################################
#%% 0. Start and initialize variables
import numpy as np
import xarray as xr
import pandas as pd
import pickle
from tqdm import tqdm

#### Import namelist to be used with the initiaization variables
namelist_input = pd.read_csv('../Data/Input_data/namelist_input.txt', sep=' ', header=0)
def get_namelist_var(name):
    return namelist_input[namelist_input.variable == name].value.values[0]

use_subset = int(get_namelist_var('use_subset'))    # 1 - Use the full input data, 2 - Use a time-cut of the input data
if use_subset == 1:
    year_i = int(get_namelist_var('year_i')); year_f = int(get_namelist_var('year_f'))
    year_file_i = year_i; year_file_f = year_f
elif use_subset == 2:
    year_i = int(get_namelist_var('year_i')); year_f = int(get_namelist_var('year_f'))
    year_file_i = int(get_namelist_var('year_file_i')); year_file_f = int(get_namelist_var('year_file_f'))
    
res = float(get_namelist_var('res'))                                   # Data resolution
region = get_namelist_var('region')                                    # Hemisphere to be analysed
data_type = get_namelist_var('data_type')                              # Data origin (ERA5 or NCAR, atm)
n_days_before = int(get_namelist_var('n_days_before'))                 # Number of days to be captured to compute LATmin


#%% 1. Functions
#%% 1.1. Check if year is leap or not
def is_leap(year):
    return year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)

#%%% 1.2. Geographical distance function
def dist(lat1, lat2, lon1, lon2):
    R = 6371.0
    lat1 = np.radians(lat1)
    lat2 = np.radians(lat2)
    lon1 = np.radians(lon1)
    lon2 = np.radians(lon2)
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    distance = R * c
    return distance

#%%% 1.3. Compute area matrix
def area_matrix(lon, lat, res):
    lat_2d = np.array([lat for i in lon]).T
    lon_2d = np.array([lon for i in lat])

    a = dist(lat_2d+res/2, lat_2d+res/2, lon_2d-res/2, lon_2d+res/2)
    b = dist(lat_2d-res/2, lat_2d-res/2, lon_2d-res/2, lon_2d+res/2)
    h = dist(lat_2d+res/2, lat_2d-res/2, lon_2d, lon_2d)

    return np.around((a+b)*h/2,1)

#%%% 1.4. Compute blocking index according to Wiedenmann et al. (2002) but extending 2D as in Davini et al. (2012)
def Blocking_intensity_index(smask, data_in_day):
    bmask = np.zeros(np.shape(smask)); bmask[:] = np.nan
    for loni in range(len(lon)):

        #### Get local geopotential
        Z_local = data_in_day[:,loni]

        #### Compute the minimum geopotential upstream
        if lon[loni]+60 >= 180:
            Z_upstream1 = data_in_day[:, loni:]
            Z_upstream2 = data_in_day[:, :np.where(lon == lon[loni]+60-360)[0][0]+1]
            Z_upstream = np.min([np.min(Z_upstream1, axis=1), np.min(Z_upstream2, axis=1)], axis=0)
        else:
            Z_upstream = np.min(data_in_day[:, loni:np.where(lon == lon[loni]+60)[0][0]+1], axis=1)

        #### Compute the minimum geopotential downstream
        if lon[loni]-60 < -180:
            Z_downstream1 = data_in_day[:, :loni+1]
            Z_downstream2 = data_in_day[:, np.where(lon == lon[loni]-60+360)[0][0]:]
            Z_downstream = np.min([np.min(Z_downstream1, axis=1), np.min(Z_downstream2, axis=1)], axis=0)
        else:
            Z_downstream = np.min(data_in_day[:, np.where(lon == lon[loni]-60)[0][0]:loni+1], axis=1)

        #### Compute the index
        RC = (((Z_upstream+Z_local)/2)+((Z_downstream+Z_local)/2))/2
        bmask[:,loni] = 100 * (Z_local/RC-1)
        
    BI_mask = np.copy(bmask)
    BI_mask[smask == 0] = np.nan
    
    return BI_mask


#%% 2. Open data
#%%% 2.1. Open tracked structures data
with open(f'../Data/Output_data/02-TrackedStruct_{year_i}_{year_f}_{region}.pkl', 'rb') as handle:
    data_dict = pickle.load(handle)

#%%% 2.2. Open original z500 data
original_data = xr.open_dataset(f'../Data/Input_data/Z500_{year_file_i}_{year_file_f}_{region}_{data_type}.nc')

#### Cut the subset if needed
if use_subset == 2:
    original_data = original_data.sel(time=slice(str(year_i), str(year_f)))
  
#### Invert the data array if needed and open the array
if region == 'NH':
    original_array = original_data.z.values
elif region == 'SH':
    original_array = original_data.z.values[:,::-1,:]

#%% 2.3. Lat and Lon and years
if region == 'NH':
    lat = original_data.latitude.values
elif region == 'SH':
    lat = -original_data.latitude.values[::-1]

lon = original_data.longitude.values

lons, lats = np.meshgrid(lon,lat)

time = original_data.time.values


#%% 3. Compute area matrix for weighted average
area = area_matrix(lon, lat, res)


#%% 4. Catalogue all structures
#%%% 4.1. Initialize arrays
Year = []                    # Year
Month = []                   # Month
Day = []                     # Day
Jul = []                     # Julian Day
Types = []                   # Type
Step = []                    # structure current duration
Duration = []                # Duration of event
Struct_ID = []               # Structure ID (entire catalogue)
Areas = []                   # area (km^2)
Area_per_prev = []           # % overlap area previous day
Area_per_dur = []            # % overlap area current day
Cent_desl = []               # mass center deslocation (km)
Cent_lats = []               # mass center latitude
Cent_lons = []               # mass center longitude
Max_block_index = []         # max daily intensity index
Mean_block_index = []        # mean daily intensity index
min_lats = []                # latitude minimum
max_lats = []                # latitude maximum
min_lons = []                # longitude minimum
max_lons = []                # longitude maximum
z500_max = []                # std Z500 maximum
z500_mean = []               # std Z500 mean

#### Masks for structs and 2D intensity
struct_array_tosave = np.zeros((len(time),len(lat),len(lon))).astype(np.float32)
BI_array_tosave = np.zeros((len(time),len(lat),len(lon))).astype(np.float32)

#%%% 4.2. Get data
for day_n in tqdm(range(len(time))):
    if day_n < n_days_before:
        continue
    else:
        data_string = time[day_n].astype(str)[:10]
        dict_day = data_dict[data_string]

        struct_array_tosave[day_n] = dict_day['Struct_array']
        BI_array_tosave[day_n] = Blocking_intensity_index(dict_day['Struct_array'], original_array[day_n])

        for strct_id in dict_day.keys():
            if strct_id != 'Struct_array' and strct_id not in Struct_ID:

                #### Invert the struct mask if it is in the SH (the resulting lat values will be reverted ahead)
                if region == 'NH':
                    struct_mask = np.where(dict_day['Struct_array'] == int(strct_id),1,0)
                elif region == 'SH':
                    struct_mask = np.where(dict_day['Struct_array'] == int(strct_id),1,0)[::-1,:]
                area_start = np.sum(struct_mask*area)

                data_in_day = original_array[day_n]

                Struct_ID.append(strct_id)
                Types.append(dict_day[strct_id]['type'])
                min_lats.append(lats[np.where(struct_mask==1)].min())
                max_lats.append(lats[np.where(struct_mask==1)].max())
                min_lons.append(lons[np.where(struct_mask==1)].min())
                max_lons.append(lons[np.where(struct_mask==1)].max())
                Year.append(int(data_string[:4]))
                Month.append(int(data_string[5:7]))
                Day.append(int(data_string[8:10]))
                Jul.append(int(pd.to_datetime(data_string).strftime('%j')))
                Areas.append(int(area_start))
                Area_per_prev.append(np.nan)
                Area_per_dur.append(np.nan)
                Step.append(dict_day[strct_id]['step'])
                Duration.append(dict_day[strct_id]['duration'])

                #### Maximum z500 value within the structure
                z500 = data_in_day*struct_mask
                z500[z500 == 0] = np.nan
                zmax = np.nanmax(z500)
                z500_max.append(round(np.nanmax(z500),1))
                z500_mean.append(round(np.nanmean(z500),1))

                #### Compute blocking index
                BI_mask = Blocking_intensity_index(struct_mask, data_in_day)
                Max_block_index.append(np.nanmax(BI_mask))
                Mean_block_index.append(np.nanmean(BI_mask))

                #### Compute the weighted center value of lon
                if np.sum(~np.isnan(BI_mask[:,0])) != 0 and np.sum(~np.isnan(BI_mask[:,-1])) != 0:
                    lons_anti = np.append(lons[:,np.shape(lons)[1]//2:]-360,
                                          lons[:,:np.shape(lons)[1]//2],
                                          axis = 1)
                    BI_anom_anti = np.append(BI_mask[:,np.shape(BI_mask)[1]//2:],
                                             BI_mask[:,:np.shape(BI_mask)[1]//2],
                                             axis = 1)
                    lon_temp = np.average(np.ma.masked_array(lons_anti, np.isnan(BI_anom_anti)),
                                          weights=np.ma.masked_array(area*(BI_anom_anti+abs(np.nanmin(BI_anom_anti))), np.isnan(BI_anom_anti)))
                    if lon_temp < -180:
                        lon_temp = lon_temp + 360
                else:
                    lon_temp = np.average(np.ma.masked_array(lons, np.isnan(BI_mask)),
                                          weights=np.ma.masked_array(area*(BI_mask+abs(np.nanmin(BI_mask))), np.isnan(BI_mask)))
                Cent_lons.append(lon_temp)
                #############################################

                #### Compute the weighted value of lat
                lat_temp = np.average(np.ma.masked_array(lats, np.isnan(BI_mask)),
                                            weights=np.ma.masked_array(area*(BI_mask+abs(np.nanmin(BI_mask))), np.isnan(BI_mask)))
                Cent_lats.append(lat_temp)

                #### Translation from one step to the next (the first is NaN)
                Cent_desl.append(np.nan)

                #### Repeat the analysis for the following days of each event
                for t in np.arange(1, dict_day[strct_id]['duration']):

                    day_after = str(time[day_n+t])[:10]

                    dict_day_after = data_dict[day_after]

                    struct_mask_after = np.where(dict_day_after['Struct_array'] == int(strct_id),1,0)

                    area_after = np.sum(struct_mask_after*area)

                    overlap_mask = struct_mask * struct_mask_after
                    overlap_area = np.sum(overlap_mask*area)

                    data_in_day = original_array[day_n+t]

                    Struct_ID.append(strct_id)
                    Types.append(dict_day_after[strct_id]['type'])
                    min_lats.append(lats[np.where(struct_mask_after==1)].min())
                    max_lats.append(lats[np.where(struct_mask_after==1)].max())
                    min_lons.append(lons[np.where(struct_mask_after==1)].min())
                    max_lons.append(lons[np.where(struct_mask_after==1)].max())
                    Year.append(int(day_after[:4]))
                    Month.append(int(day_after[5:7]))
                    Day.append(int(day_after[8:10]))
                    Jul.append(int(pd.to_datetime(day_after).strftime('%j')))
                    Areas.append(int(area_after))
                    Area_per_prev.append(np.round(overlap_area/area_start,3))
                    Area_per_dur.append(np.round(overlap_area/area_after,3))
                    Step.append(dict_day_after[strct_id]['step'])
                    Duration.append(dict_day_after[strct_id]['duration'])

                    z500 = data_in_day*struct_mask_after
                    z500[z500 == 0] = np.nan
                    zmax = np.nanmax(z500)
                    z500_max.append(round(np.nanmax(z500),1))
                    z500_mean.append(round(np.nanmean(z500),1))

                    ### Compute blocking index
                    BI_mask = Blocking_intensity_index(struct_mask_after, data_in_day)
                    Max_block_index.append(np.nanmax(BI_mask))
                    Mean_block_index.append(np.nanmean(BI_mask))

                    #### Compute the weighted center value of lon
                    if np.sum(~np.isnan(BI_mask[:,0])) != 0 and np.sum(~np.isnan(BI_mask[:,-1])) != 0:
                        lons_anti = np.append(lons[:,np.shape(lons)[1]//2:]-360,
                                              lons[:,:np.shape(lons)[1]//2],
                                              axis = 1)
                        BI_anom_anti = np.append(BI_mask[:,np.shape(BI_mask)[1]//2:],
                                                 BI_mask[:,:np.shape(BI_mask)[1]//2],
                                                 axis = 1)
                        lon_temp_after = np.average(np.ma.masked_array(lons_anti, np.isnan(BI_anom_anti)),
                                                    weights=np.ma.masked_array(area*(BI_anom_anti+abs(np.nanmin(BI_anom_anti))), np.isnan(BI_anom_anti)))
                        if lon_temp < -180:
                            lon_temp_after = lon_temp + 360
                    else:
                        lon_temp_after = np.average(np.ma.masked_array(lons, np.isnan(BI_mask)),
                                              weights=np.ma.masked_array(area*(BI_mask+abs(np.nanmin(BI_mask))), np.isnan(BI_mask)))
                    Cent_lons.append(lon_temp_after)
                    #############################################

                    #### Compute the weighted value of lat
                    lat_temp_after = np.average(np.ma.masked_array(lats, np.isnan(BI_mask)),
                                          weights=np.ma.masked_array(area*(BI_mask+abs(np.nanmin(BI_mask))), np.isnan(BI_mask)))
                    Cent_lats.append(lat_temp_after)

                    Cent_desl.append(dist(lat_temp, lat_temp_after, lon_temp, lon_temp_after))

                    lon_temp = lon_temp_after
                    lat_temp = lat_temp_after
                    struct_mask = struct_mask_after
                    area_start = np.sum(struct_mask*area)

#%%% 4.3. Create data in dataframe and save
if region == 'NH':
    dict_to_save = {'SID': np.array(Struct_ID),
                    'YEAR': np.array(Year),
                    'MONTH': np.array(Month),
                    'DAY': np.array(Day),
                    'TYPE': np.array(Types),
                    'CLON': np.array(Cent_lons),
                    'CLAT': np.array(Cent_lats),
                    'DESL': np.array(Cent_desl),
                    'LATMIN': np.array(min_lats),
                    'LATMAX': np.array(max_lats),
                    'LONMIN': np.array(min_lons),
                    'LONMAX': np.array(max_lons),
                    'JUL': np.array(Jul),
                    'AREA': np.array(Areas),
                    'AREA_PER_PREV': np.array(Area_per_prev),
                    'AREA_PER_DUR': np.array(Area_per_dur),
                    'STEP': np.array(Step),
                    'DURATION': np.array(Duration),
                    'BI_MAX': np.array(Max_block_index),
                    'BI_MEAN': np.array(Mean_block_index),
                    'Z500_MAX': np.array(z500_max),
                    'Z500_MEAN': np.array(z500_mean)}
elif region == 'SH':
    dict_to_save = {'SID': np.array(Struct_ID),
                    'YEAR': np.array(Year),
                    'MONTH': np.array(Month),
                    'DAY': np.array(Day),
                    'TYPE': np.array(Types),
                    'CLON': np.array(Cent_lons),
                    'CLAT': -np.array(Cent_lats),
                    'DESL': np.array(Cent_desl),
                    'LATMIN': -np.array(min_lats),
                    'LATMAX': -np.array(max_lats),
                    'LONMIN': np.array(min_lons),
                    'LONMAX': np.array(max_lons),
                    'JUL': np.array(Jul),
                    'AREA': np.array(Areas),
                    'AREA_PER_PREV': np.array(Area_per_prev),
                    'AREA_PER_DUR': np.array(Area_per_dur),
                    'STEP': np.array(Step),
                    'DURATION': np.array(Duration),
                    'BI_MAX': np.array(Max_block_index),
                    'BI_MEAN': np.array(Mean_block_index),
                    'Z500_MAX': np.array(z500_max),
                    'Z500_MEAN': np.array(z500_mean)}

blocks = pd.DataFrame(data = dict_to_save)
blocks.to_csv(f'../Data/Output_data/03-Blocking_daily_catalogue_{year_i}_{year_f}_{region}.csv',index=False)

#%% 4.4. Create netcdfs for masks
#### Create dataset considering the region chosen
if region == 'NH':
    data = xr.Dataset(
        data_vars = dict(Structs=(['time', 'lat', 'lon'], struct_array_tosave),
                         Intensity=(['time', 'lat', 'lon'], BI_array_tosave)),
        coords = dict(time=time, lat=lat, lon=lon),
        attrs = dict(description='Structures and Blocking intensity based on Sousa et al 2021')
        )
elif region == 'SH':
    data = xr.Dataset(
        data_vars = dict(Structs=(['time', 'lat', 'lon'], struct_array_tosave[:,::-1,:]),
                         Intensity=(['time', 'lat', 'lon'], BI_array_tosave[:,::-1,:])),
        coords = dict(time=time, lat=-lat[::-1], lon=lon),
        attrs = dict(description='Structures and Blocking intensity based on Sousa et al 2021')
        )

data.to_netcdf(f'../Data/Output_data/03-CatalogueMasks_{year_i}_{year_f}_{region}.nc')


#%% 5. Create dataframe for event catalogue (this is simpler, just means, maxs, mins, and counts)
Event_ID = np.unique(blocks.SID.values)
Event_dom = []
Event_year = []
Event_month = []
Event_day = []
Event_duration = []
Event_intmean = []
Event_intmax = []
Event_latmean = []
Event_lonmean = []
Event_areamean = []
Event_z500absmax = []
Event_z500absmean = []
Event_areamax = []
Event_overlapmean = []
Event_deslocmean = []
Event_percridge = []
Event_percomega = []
Event_perchybrid = []
Event_percrex = []
Event_percpolar = []

for strct_id in tqdm(Event_ID):
    strct_to_eval = blocks[blocks.SID == strct_id]

    ########## Other attributes
    Event_year.append(strct_to_eval.YEAR.values[0])
    Event_month.append(strct_to_eval.MONTH.values[0])
    Event_day.append(strct_to_eval.DAY.values[0])
    Event_latmean.append(strct_to_eval.CLAT.values.mean())
    Event_lonmean.append(strct_to_eval.CLON.values.mean())
    Event_deslocmean.append(strct_to_eval.DESL.values.mean())
    Event_duration.append(strct_to_eval.STEP.values[-1])
    Event_z500absmax.append(strct_to_eval.Z500_MAX.max())
    Event_z500absmean.append(strct_to_eval.Z500_MAX.mean())
    Event_overlapmean.append(strct_to_eval.AREA_PER_DUR.mean())
    Event_areamean.append(strct_to_eval.AREA.mean())
    Event_areamax.append(strct_to_eval.AREA.max())
    Event_intmean.append(np.nanmean(strct_to_eval.BI_MEAN.values))
    Event_intmax.append(np.nanmax(strct_to_eval.BI_MAX.values))

    ########### Check dominant type and percentages
    count_strcts = strct_to_eval.groupby('TYPE').size()
    prevalent = count_strcts.keys()[count_strcts == count_strcts.max()].values

    # Ridge
    try:
        Event_percridge.append(count_strcts['Ridge']/len(strct_to_eval))
    except:
        Event_percridge.append(0)

    # Omega
    try:
        Event_percomega.append(count_strcts['Omega block']/len(strct_to_eval))
    except:
        Event_percomega.append(0)

    # Hybrid
    try:
        Event_perchybrid.append(count_strcts['Rex block (block)']/len(strct_to_eval))
    except:
        Event_perchybrid.append(0)

    # Rex
    try:
        Event_percrex.append(count_strcts['Rex block']/len(strct_to_eval))
    except:
        Event_percrex.append(0)

    # Polar
    try:
        Event_percpolar.append(count_strcts['Rex block (polar)']/len(strct_to_eval))
    except:
        Event_percpolar.append(0)

    # Dominant type
    if len(prevalent) == 1:
            prevalent = prevalent[0]
            Event_dom.append(prevalent)
    else:
        if 'Rex block (polar)' in prevalent:
            prevalent = 'Rex block (polar)'
            Event_dom.append(prevalent)
            continue
        elif 'Rex block' in prevalent:
            prevalent = 'Rex block'
            Event_dom.append(prevalent)
            continue
        elif 'Rex block (hybrid)' in prevalent:
            prevalent = 'Rex block (hybrid)'
            Event_dom.append(prevalent)
            continue
        elif 'Omega block' in prevalent:
            prevalent = 'Omega block'
            Event_dom.append(prevalent)
            continue
        else:
            break

#%%% 5.1. Create data in dataframe and save
dict_to_save = {'SID': np.array(Event_ID),
                'DOM_TYPE': np.array(Event_dom),
                'DESL_MEAN': np.array(Event_deslocmean),
                'CLAT_MEAN': np.array(Event_latmean),
                'CLON_MEAN': np.array(Event_lonmean),
                'YEAR_START': np.array(Event_year),
                'MONTH_START': np.array(Event_month),
                'DAY_START': np.array(Event_day),
                'AREA_MAX': np.array(Event_areamax),
                'AREA_MEAN': np.array(Event_areamean),
                'OVERLAP_MEAN': np.array(Event_overlapmean),
                'DURATION': np.array(Event_duration),
                'BI_MAX': np.array(Event_intmax),
                'BI_MEAN': np.array(Event_intmean),
                'Z500_MAX': np.array(Event_z500absmax),
                'Z500_MEAN': np.array(Event_z500absmean),
                'PERC_RIDGE': np.array(Event_percridge),
                'PERC_OMEGA': np.array(Event_percomega),
                'PERC_HYBRID': np.array(Event_perchybrid),
                'PERC_REX': np.array(Event_percrex),
                'PERC_POLAR': np.array(Event_percpolar),}
blocks = pd.DataFrame(data = dict_to_save)
blocks.to_csv(f'../Data/Output_data/03-Blocking_event_catalogue_{year_i}_{year_f}_{region}.csv',index=False)