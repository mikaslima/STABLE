##################################################################
#
#
# Script that tracks the structures in time according to the conditions in Sousa et al. (2021)
# Takes LATmin and continuous structure information as input
# Outputs a pkl with continuous structures in space and time respecting structure thresholds, as well as their information
#
#
##################################################################
#%% 0. Start
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
area_threshold = float(get_namelist_var('area_threshold'))             # Overlap from structures in following days
persistence = float(get_namelist_var('persistence'))                   # Minimum number of days a structure needs to exist

# Tracking method parameters (refer to REAME or ppt for further explanation):
# 1 considers excedance of the area_threshold in the day of the analysis or next (as in Sousa et al., 2021);
# 2 considers excedance in both days;
# 3 considers only on the next day and 4 on the day of the analysis.
tracking_method = int(get_namelist_var('tracking_method'))


#%% 1. Functions
#%%% 1.1. Geographical distance function
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

#%%% 1.2. Compute area matrix
def area_matrix(lon, lat, res):
    lat_2d = np.array([lat for i in lon]).T
    lon_2d = np.array([lon for i in lat])

    a = dist(lat_2d+res/2, lat_2d+res/2, lon_2d-res/2, lon_2d+res/2)
    b = dist(lat_2d-res/2, lat_2d-res/2, lon_2d-res/2, lon_2d+res/2)
    h = dist(lat_2d+res/2, lat_2d-res/2, lon_2d, lon_2d)

    return np.around((a+b)*h/2,1)


#%% 2. Open data
with open(f'../Data/Output_data/01-StructTypes_{year_i}_{year_f}_{region}.pkl', 'rb') as handle:
    struct_type = pickle.load(handle)

struct_mask = xr.open_dataset(f'../Data/Output_data/01-StructMasks_{year_i}_{year_f}_{region}.nc')
struct_mask_array = struct_mask.Structs.values

#%%% 2.1. Lat, Lon and time
lat = struct_mask.lat.values
lon = struct_mask.lon.values
time = struct_mask.time.values


#%% 3. Compute area matrix for weighted average
area = area_matrix(lon, lat, res)


#%% 4. Track structures per day of start
data_dict = {}
struct_array_tosave = np.zeros((len(time),len(lat),len(lon)))

for str_dates in pd.Series(time.astype(str)).str[:10].values:
    data_dict[str_dates] = {}
    data_dict[str_dates]['Struct_array'] = np.zeros(np.shape(area), dtype=np.float32)

for day_n in tqdm(range(len(time))):
    data_string = time[day_n].astype(str)[:10]

    #### Initialize counter of structures for each year
    if time[day_n].astype(str)[5:10] == '01-01':
        structs_in_year = 0

    #### Don't count the first 15 days
    if day_n < 15:
        continue

    #### Start analysis
    else:
        struct_mask_day = struct_mask_array[day_n]
        struct_info_day = struct_type[data_string]

        #### Check for this day if there are any already analysed blockings
        list_keys = []
        for existing_ids in data_dict[data_string].keys():
            if existing_ids != 'Struct_array':
                list_keys.append(existing_ids)
        list_keys = np.unique(list_keys)

        #### If any blockings have already been analysed save the "parent" on a list so they are not taken into account again (this approach saves a lot of analysis time)
        #### FLAG: THIS MAY DIFFER SIGNIFICANTLY FROM PEDRO's APPROACH
        #### cont. Not sure how Pedro dealt with previously analysed and discarded structures
        #### Structures that separate in two and then could be considered two blockings in the right conditions are not taken in account using this approach
        analysed_parents = []
        if len(list_keys) != 0:
            for id_key in list_keys:
                analysed_parents.append(data_dict[data_string][id_key]['parent'])

        #### Analyse structures in day
        for struct_num in list(struct_info_day.keys()):
            ######### Determine if the structure has been analysed before
            if struct_num in analysed_parents:
                continue                                   #### If it has been analysed before, ignore and continue to the next one
            elif struct_num not in analysed_parents:       #### Otherwise start analysis on tentative new blocking
                struct_mask = np.zeros(np.shape(area))
                struct_mask[struct_mask_day == int(struct_num)] = 1    #### Mask of starting day for the structure alone

                struct_3d_arr = [struct_mask]    #### Initialize "cube" of masks for the blocking

                continue_y = True    #### Flag to continue the analysis
                count_n = 1          #### Count the number of days of the tentative blocking
                saved_struct_nums = [(data_string,int(struct_num))]   #### Save the day and the number of the structure in that day

                while continue_y:    #### While the flag is active, keep going with the analysis, this advances the analysis in time

                    if day_n+count_n < len(time):    #### If the day being analised is not greater than the available dataset, keep the analysis going

                        day_after = str(pd.to_datetime(time)[day_n+count_n])[:10]    #### String of the day being analysed
                        struct_info_after = struct_type[day_after]                   #### Information about the structures in the day
                        struct_mask_after = struct_mask_array[day_n+count_n]         #### Mask for the structures in the day being analysed

                        contained_structs = struct_3d_arr[-1]*struct_mask_after      #### Multiply to see how many and which structures overlap with the one in the previous day
                        area_before = int(np.sum(struct_3d_arr[-1]*area))            #### Area of structure in the day being analysed

                        #### Compute the possible several overlapping areas (most of the time there will be either 0 or 1, and we always keep the highest)
                        contained_areas_after = []; contained_areas_before = []; contained_areas_n = []            #### Store the contained areas inside a list to check if the analysis goes on or not
                        for cont_n in np.unique(contained_structs).astype(int):
                            if cont_n != 0:
                                contained_array = np.zeros(np.shape(area))
                                contained_array[contained_structs == cont_n] = 1     #### Mask of overlapping block
                                contained_area = int(np.sum(contained_array*area))   #### Overlap area

                                temp_mask_after = np.zeros(np.shape(area))
                                temp_mask_after[struct_mask_after == cont_n] = 1     #### Mask of next day block
                                temp_area = int(np.sum(temp_mask_after*area))        #### Area of next day block

                                perc_contained_after = contained_area/temp_area      #### Percentage of area of overlap with the next day block
                                perc_contained_before = contained_area/area_before   #### Percentage of area of overlap with the present day block

                                contained_areas_after.append(perc_contained_after)
                                contained_areas_before.append(perc_contained_before)
                                contained_areas_n.append(cont_n)

                        contained_areas_after = np.array(contained_areas_after)
                        contained_areas_before = np.array(contained_areas_before)
                        contained_areas_n = np.array(contained_areas_n)

                        if len(contained_areas_after) != 0:          #### If there are at least one structure overlapping then proceed with the analysis
                            biggest = max(contained_areas_after)

                            #### Check if there are valid structures depending on the tracking method
                            valid_structs = []; valid_area_after = []; valid_area_before = []      #### Duplicates of "contained_areas_..." but for valid conditions
                            for i, j, k in zip(contained_areas_after, contained_areas_before, contained_areas_n):
                                if tracking_method == 1:           #### Checks if overlapping structure exceeds the threshold (typically 50%) compared to the next day OR the day of the analysis (Same as in Sousa et al., 2021)
                                    if i >= area_threshold or j >= area_threshold:
                                        valid_structs.append(k)
                                        valid_area_after.append(i)
                                        valid_area_before.append(j)

                                elif tracking_method == 2:           #### Checks if overlapping structure exceeds the threshold (typically 50%) compared to the next day AND the day of the analysis (Technically more correct approach)
                                    if i >= area_threshold and j >= area_threshold:
                                        valid_structs.append(k)
                                        valid_area_after.append(i)
                                        valid_area_before.append(j)

                                elif tracking_method == 3:           #### Checks if overlapping structure exceeds the threshold (typically 50%) compared to the next day
                                    if i >= area_threshold:
                                        valid_structs.append(k)
                                        valid_area_after.append(i)
                                        valid_area_before.append(j)

                                elif tracking_method == 4:           #### Checks if overlapping structure exceeds the threshold (typically 50%) compared to the day of the analysis
                                    if j >= area_threshold:
                                        valid_structs.append(k)
                                        valid_area_after.append(i)
                                        valid_area_before.append(j)

                            #### Check if there were any found and proceed with the analysis
                            if len(valid_structs) >= 1:
                                if len(valid_structs) == 1:      #### In the case of just a single valid structure there are no problems
                                    which_next = valid_structs[0]    #### Structure that "heirs" the blocking to the next step

                                elif len(valid_structs) > 1:     #### For more than one, we need to consider the tracking method to choose the structure to follow

                                    if tracking_method == 1:             #### Checks if overlapping structure exceeds the threshold (typically 50%) compared to the next day OR the day of the analysis (Same as in Sousa et al., 2021)
                                        temp_list_and = []; temp_structs_and = [] #### Check if there are any that overlap behyond the threshold both on the day and after
                                        temp_list_or = []; temp_structs_or = []   #### Otherwise, we only have the "or" condition

                                        for i, j, k in zip(valid_area_after, valid_area_before, valid_structs):
                                            if i >= area_threshold and j >= area_threshold:
                                                temp_list_and.append(i*j); temp_structs_and.append(k)
                                            else:
                                                temp_list_or.append(i*j); temp_structs_or.append(k)

                                        if len(temp_list_and) != 0:   #### If both the overlap in the day of the analysis and the day after exceed the threshold then they are prioritized
                                            biggest = max(temp_list_and)
                                            which_next = temp_structs_and[np.where(temp_list_and == biggest)[0][0]]   #### Which structure "heirs" the blocking to the next step
                                        else:                         #### Otherwise consider the missing ones
                                            biggest = max(temp_list_or)
                                            which_next = temp_structs_or[np.where(temp_list_or == biggest)[0][0]]   #### Which structure "heirs" the blocking to the next step

                                    elif tracking_method == 2:           #### Checks if overlapping structure exceeds the threshold (typically 50%) compared to the next day AND the day of the analysis (Technically more correct approach)
                                        biggest = max(valid_area_after*valid_area_before)
                                        which_next = valid_structs[np.where(valid_area_after*valid_area_before == biggest)[0][0]]   #### Which structure "heirs" the blocking to the next step

                                    elif tracking_method == 3:           #### Checks if overlapping structure exceeds the threshold (typically 50%) compared to the next day
                                        biggest = max(valid_area_after)
                                        which_next = valid_structs[np.where(valid_area_after == biggest)[0][0]]   #### Which structure "heirs" the blocking to the next step

                                    elif tracking_method == 4:           #### Checks if overlapping structure exceeds the threshold (typically 50%) compared to the day of the analysis
                                        biggest = max(valid_area_before)
                                        which_next = valid_structs[np.where(valid_area_before == biggest)[0][0]]   #### Which structure "heirs" the blocking to the next step


                                #### Check for this day if there are any already analysed blockings
                                list_keys = []
                                for existing_ids in data_dict[day_after].keys():
                                    if existing_ids != 'Struct_array':
                                        list_keys.append(existing_ids)
                                list_keys = np.unique(list_keys)

                                #### If any blockings have already been analysed save the "parent" on a list so they are not taken into account again
                                analysed_parents = []
                                if len(list_keys) != 0:
                                    for id_key in list_keys:
                                        analysed_parents.append(data_dict[day_after][id_key]['parent'])

                                if str(which_next) not in analysed_parents:
                                    struct_to_save = np.zeros(np.shape(area))
                                    struct_to_save[struct_mask_after == which_next] = 1    #### Mask of this "heir" structure
                                    saved_struct_nums.append((day_after,which_next))       #### Save the day and number of the structure in that day

                                    struct_3d_arr.append(struct_to_save)    #### Append to the blocking cube

                                    continue_y = True      #### Flag is updated to keep the analysis going
                                    count_n += 1           #### Increment the day by one

                                elif str(which_next) in analysed_parents:
                                    if count_n < persistence:
                                        continue_y = False

                                    elif count_n >= persistence:

                                        struct_3d_arr = np.array(struct_3d_arr)

                                        ############ Block to save data
                                        structs_in_year += 1
                                        give_id = f'{data_string[:4]}{str(structs_in_year).zfill(3)}'

                                        for i, (date_to_save, num_of_struct) in enumerate(saved_struct_nums):
                                            n_to_save = np.where(pd.Series(time.astype(str)).str[:10].values == date_to_save)[0][0]
                                            temp_info = struct_type[date_to_save]
                                            temp_mask = struct_mask_array[n_to_save]

                                            data_dict[date_to_save][give_id] = {'type': temp_info[str(num_of_struct)],
                                                                                'parent': str(num_of_struct),
                                                                                'step': i+1,
                                                                                'duration': len(saved_struct_nums)}

                                            data_dict[date_to_save]['Struct_array'][temp_mask == num_of_struct] = np.float32(give_id)
                                        ############
                                        continue_y = False

                            #### If there are no valid overlapping structures then end that structure's analysis
                            elif len(valid_structs) == 0:
                                if count_n < persistence:          #### If the biggest overlapping area does not exceed the threshold and there are not enough days in the analysis drop the tentative blocking
                                    continue_y = False

                                elif count_n >= persistence:       #### If the biggest overlapping area does not exceed the threshold but the structure is semi-stationary save the data, new blocking found!
                                    struct_3d_arr = np.array(struct_3d_arr)

                                    ############ Block to save data
                                    structs_in_year += 1
                                    give_id = f'{data_string[:4]}{str(structs_in_year).zfill(3)}'

                                    for i, (date_to_save, num_of_struct) in enumerate(saved_struct_nums):
                                        n_to_save = np.where(pd.Series(time.astype(str)).str[:10].values == date_to_save)[0][0]
                                        temp_info = struct_type[date_to_save]
                                        temp_mask = struct_mask_array[n_to_save]

                                        data_dict[date_to_save][give_id] = {'type': temp_info[str(num_of_struct)],
                                                                            'parent': str(num_of_struct),
                                                                            'step': i+1,
                                                                            'duration': len(saved_struct_nums)}

                                        data_dict[date_to_save]['Struct_array'][temp_mask == num_of_struct] = np.float32(give_id)

                                    ############
                                    continue_y = False

                        else:     #### if no structure shares area then end and save data if conditions are met
                            if count_n < persistence:
                                continue_y = False

                            elif count_n >= persistence:

                                struct_3d_arr = np.array(struct_3d_arr)

                                ############ Block to save data
                                structs_in_year += 1
                                give_id = f'{data_string[:4]}{str(structs_in_year).zfill(3)}'

                                for i, (date_to_save, num_of_struct) in enumerate(saved_struct_nums):
                                    n_to_save = np.where(pd.Series(time.astype(str)).str[:10].values == date_to_save)[0][0]
                                    temp_info = struct_type[date_to_save]
                                    temp_mask = struct_mask_array[n_to_save]

                                    data_dict[date_to_save][give_id] = {'type': temp_info[str(num_of_struct)],
                                                                        'parent': str(num_of_struct),
                                                                        'step': i+1,
                                                                        'duration': len(saved_struct_nums)}

                                    data_dict[date_to_save]['Struct_array'][temp_mask == num_of_struct] = np.float32(give_id)

                                ############
                                continue_y = False

                    elif day_n+count_n >= len(time):    #### If the day being analysed exceeds the dataset end the analysis and check if data is to save or not

                        if count_n < persistence:
                            continue_y = False

                        elif count_n >= persistence:

                            struct_3d_arr = np.array(struct_3d_arr)

                            ############ Block to save data
                            structs_in_year += 1
                            give_id = f'{data_string[:4]}{str(structs_in_year).zfill(3)}'

                            for i, (date_to_save, num_of_struct) in enumerate(saved_struct_nums):
                                n_to_save = np.where(pd.Series(time.astype(str)).str[:10].values == date_to_save)[0][0]
                                temp_info = struct_type[date_to_save]
                                temp_mask = struct_mask_array[n_to_save]

                                data_dict[date_to_save][give_id] = {'type': temp_info[str(num_of_struct)],
                                                                    'parent': str(num_of_struct),
                                                                    'step': i+1,
                                                                    'duration': len(saved_struct_nums)}

                                data_dict[date_to_save]['Struct_array'][temp_mask == num_of_struct] = np.float32(give_id)

                            ############
                            continue_y = False


#%% 5. Save data
with open(f'../Data/Output_data/02-TrackedStruct_{year_i}_{year_f}_{region}.pkl', 'wb') as handle:
        pickle.dump(data_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)