##################################################################
#
#
# Script that identifies individual structures based on the charachteristics imposed in Sousa et al. (2021)
# Takes Z500 as input
# Outputs a netcdf with type of structure in each grid cell and a pkl file with the type of structure associated to each days structures
#
#
##################################################################
#%% 0. Start and initialize variables
import numpy as np
import xarray as xr
from tqdm import tqdm
from scipy import ndimage
from scipy.ndimage import label, binary_dilation
import pickle
from scipy.signal import savgol_filter
import pandas as pd

#### Import namelist to be used with the initiaization variables
namelist_input = pd.read_csv('../Data/Input_data/namelist_input.txt', sep=' ', header=0)
def get_namelist_var(name):
    return namelist_input[namelist_input.variable == name].value.values[0]

use_subset = int(get_namelist_var('use_subset'))                       # 1 - Use the full input data, 2 - Use a time-cut of the input data
if use_subset == 1:
    year_i = int(get_namelist_var('year_i')); year_f = int(get_namelist_var('year_f'))
    year_file_i = year_i; year_file_f = year_f
elif use_subset == 2:
    year_i = int(get_namelist_var('year_i')); year_f = int(get_namelist_var('year_f'))
    year_file_i = int(get_namelist_var('year_file_i')); year_file_f = int(get_namelist_var('year_file_f'))

res = float(get_namelist_var('res'))                                   # Data resolution
region = get_namelist_var('region')                                    # Hemisphere to be analysed
data_type = get_namelist_var('data_type')                              # Data origin (ERA5 or NCAR, atm)
min_struct_area = float(get_namelist_var('min_struct_area'))           # Minimum structure area to be captured
n_days_before = int(get_namelist_var('n_days_before'))                 # Number of days to be captured to compute LATmin
horizontal_LATmin = int(get_namelist_var('horizontal_LATmin'))         # 1 - makes the computed LATmin be the same as Sousa et al., 2021; 2 - makes the new moving LATmin based on data from the last n_days_before
omega_hybrid_method = int(get_namelist_var('omega_hybrid_method'))     # 1 - Simple distinction between omega and hybrid as described in Sousa et al. (2021), 2 - New distinction method with hybrids in mixed systems
consider_polar = int(get_namelist_var('consider_polar'))               # 1 - Consider regions poleward of 90-delta as Sousa et al. (2021), 2 - Consider regions poleward of 90-delta, 3 - Do not consider regions poleward of 90-delta
lat_polar_circle = float(get_namelist_var('lat_polar_circle'))         # Latitude (decimal) of the polar circle  to be considered

#### Establish latitude delta of analysis for gradient computation
delta = int(get_namelist_var('delta')) # degrees south or north


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

#%%% 1.3. Find nearest value in array
def find_nearest(array, value, which):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    if which == 'value':
        return array[idx]
    elif which == 'index':
        return idx

#%%% 1.4. Function to clean holes in the data
def clean_holes(array, which):
    out = np.copy(array)
    if which == 'hole':
        out[~np.isnan(out)] = 1
        out[np.isnan(out)] = 0
    elif which == 'excess':
        out[~np.isnan(out)] = 0
        out[np.isnan(out)] = 1
    out = ndimage.binary_fill_holes(out).astype(float)
    out[out == 0] = np.nan
    return out

#%%% 1.5. Identify isolated structures
def identify_structures(array):
    grow = binary_dilation(array, structure=np.ones((1, 1), dtype=int)) # If you want to change the degrees of freedom change the numbers inside np.ones(), 1 degree was found as ideal in this approach
    lbl, npatches = label(grow)
    lbl[array==0] = 0
    return lbl

#%%% 1.6. Compute gradient arrays
def gradients_function(day_n, LATmin_day, data_in_day):
    #### Location of LATmin as index
    LATmin_loc = [find_nearest(lat, lati, 'index') for lati in LATmin_day]

    #### Initialize arrays in each day
    GHG_lat = np.zeros((len(lat),len(lon)),dtype=np.float32); GHG_lat[:] = np.nan
    GHG_lon = np.zeros((len(lat),len(lon)),dtype=np.float32); GHG_lon[:] = np.nan
    GHGS = np.zeros((len(lat),len(lon)),dtype=np.float32); GHGS[:] = np.nan
    GHGS2 = np.zeros((len(lat),len(lon)),dtype=np.float32); GHGS2[:] = np.nan
    GHGN = np.zeros((len(lat),len(lon)),dtype=np.float32); GHGN[:] = np.nan

    #### Find the boudary of the analysis
    limit_pos_north_others = find_nearest(lat, 90-delta, 'index')
    limit_pos_north_GHG = find_nearest(lat, 90-delta/2, 'index')
    limit_pos_south_GHG = find_nearest(lat, delta/2, 'index')
    limit_pos_south_others = find_nearest(lat, delta, 'index')

    #### Produce the analysis at latitude of the desired hemisphere
    for lat_pos in range(len(lat)):

        #### Compute the lat component of GHG for each lat between the boundary (bounded both in the north and south)
        #### FLAG: THIS MAY DIFFER SIGNIFICANTLY FROM PEDRO's APPROACH (Not clear how GHG was computed in the paper)
        if lat_pos > limit_pos_north_GHG and lat_pos < limit_pos_south_GHG:
            lat_pos_below_h = find_nearest(lat, lat[lat_pos]-delta/2, 'index')
            lat_pos_above_h = find_nearest(lat, lat[lat_pos]+delta/2, 'index')
            GHG_lat[lat_pos] = (data_in_day[lat_pos_above_h]-data_in_day[lat_pos_below_h])/delta

        #### Compute the lat component of GHG for the lats above the north boundary towards the pole
        elif lat_pos <= limit_pos_north_GHG:
            lat_pos_below_h = find_nearest(lat, lat[lat_pos]-delta/2, 'index')
            lat_pos_above_h = find_nearest(lat, 180-lat[lat_pos]-delta/2, 'index')
            GHG_lat[lat_pos] = (data_in_day[lat_pos_above_h,lon_opposite_arr]-data_in_day[lat_pos_below_h])/delta

        #### Compute GHGS, GHGS2, and GHGN for each lat
        if lat_pos <= limit_pos_south_others:
            lat_pos_below = find_nearest(lat, lat[lat_pos]-delta, 'index')
            GHGS[lat_pos] = (data_in_day[lat_pos]-data_in_day[lat_pos_below])/delta

            lat_pos_below2 = find_nearest(lat, lat[lat_pos]-delta*2, 'index')
            GHGS2[lat_pos] = (data_in_day[lat_pos_below]-data_in_day[lat_pos_below2])/delta

            #### GHGN can just be computed up to a certain lat, like GHG before
            if lat_pos >= limit_pos_north_others:
                lat_pos_above = find_nearest(lat, lat[lat_pos]+delta, 'index')
                GHGN[lat_pos] = (data_in_day[lat_pos_above]-data_in_day[lat_pos])/delta

            #### To compute the GHGN over 75º up to the pole we make the inversion to the other side of the hemisphere
            elif lat_pos < limit_pos_north_others:
                #### Retrieve the opposite latitude
                opposite_lat = find_nearest(lat, 180-lat[lat_pos]-delta, 'index')

                #### Compute the gradient
                GHGN[lat_pos] = (data_in_day[opposite_lat, lon_opposite_arr]-data_in_day[lat_pos])/delta

    #### Compute the lon component of GHG for each lon, we use the same loop to clean values equatorward of LATmin
    #### FLAG: THIS MAY DIFFER SIGNIFICANTLY FROM PEDRO's APPROACH (Not clear how GHG was computed in the paper)
    for lon_pos in np.arange(len(lon)):
        #### Save the before and after positions to compute the gradient and to check for the boundary
        lon_pos_before = find_nearest(lon, lon[lon_pos]-delta/2, 'index')
        lon_pos_after = find_nearest(lon, lon[lon_pos]+delta/2, 'index')

        #### West boundary correction (cyclic)
        if lon[lon_pos]-delta/2 < -180:
            lon_pos_before = find_nearest(lon, lon[lon_pos]-delta/2+360, 'index')

        #### East boundary correction (cyclic)
        if lon[lon_pos]+delta/2 > 180:
            lon_pos_after = find_nearest(lon, lon[lon_pos]+delta/2-360, 'index')

        #### Compute the GHG lon component
        GHG_lon[:,lon_pos] = (data_in_day[:,lon_pos_before]-data_in_day[:,lon_pos_after])/(delta*np.cos(np.deg2rad(lat)))

        #### Clean values below LATmin
        GHG_lat[LATmin_loc[lon_pos]:, lon_pos] = np.nan
        GHG_lon[LATmin_loc[lon_pos]:, lon_pos] = np.nan
        GHGN[LATmin_loc[lon_pos]:, lon_pos] = np.nan
        GHGS[LATmin_loc[lon_pos]:, lon_pos] = np.nan
        GHGS2[LATmin_loc[lon_pos]:, lon_pos] = np.nan

    return GHG_lat, GHG_lon, GHGS, GHGS2, GHGN

#%%% 1.7. Static LATmin
def compute_static_LATmin(day_n):
    ############ Sousa et al., 2021 methodology: Fixed LATmin with longitude
    #### Agregate data from last "n_days_before" days
    data_for_mean = original_array[day_n-n_days_before:day_n]

    #### Compute the zonal mean of the aggregated data
    zonal = savgol_filter(np.mean(data_for_mean, axis=(0,2)), 6, 3)

    auxMED = []; auxSUM = []

    for i in range(zonal.shape[0]):
        MED = zonal[i] * np.cos(lat_rad[i])
        SUM = np.cos(lat_rad[i])

        auxMED.append(MED)
        auxSUM.append(SUM)

    z500_thres = np.sum(auxMED)/np.sum(auxSUM)
    diff = np.abs(zonal - z500_thres)

    return [lat[np.argmin(diff)] for i in lon], z500_thres

#%%% 1.8. Moving LATmin method
def compute_moving_LATmin(day_n):
    ############ Compute LATmin based on data from last 15 days
    #### Agregate data from last "n_days_before" days
    data_for_mean = original_array[day_n-n_days_before:day_n]

    #### Compute the zonal mean of the aggregated data
    zonal = savgol_filter(np.mean(data_for_mean, axis=(0,2)), 6, 3)

    auxMED = []; auxSUM = []

    for i in range(zonal.shape[0]):
        MED = zonal[i] * np.cos(lat_rad[i])
        SUM = np.cos(lat_rad[i])

        auxMED.append(MED)
        auxSUM.append(SUM)

    mean_val = np.sum(auxMED)/np.sum(auxSUM)

    #### Mean array for the last "n_days_before" days
    data_in_day = np.average(data_for_mean, axis=0)
    
    #### Keep data below the hemispheric mean val for the last "n_days_before" days
    data_in_day[data_in_day < mean_val] = np.nan

    #### Clean "holes" in data with a binary dilation algorithm, this dilation has one degree of freedom (i.e., only immediately touching pixels are accounted, no diagonals)
    new_data_in_day = clean_holes(data_in_day, 'hole')

    #### Compute the min lat for each lon based on daily values from last 15 days
    LATarray_moving = []
    for i in range(np.shape(new_data_in_day)[1]):
        where_data = np.where(~np.isnan(new_data_in_day[:,i]))[0]
        if np.max(np.diff(where_data)) != 1:
            where = np.where(np.diff(where_data) != 1)[0][0]
            where_data = where_data[where+1:]
        LATarray_moving.append(lat[where_data[0]])

    return LATarray_moving, mean_val


#%% 2. Open data
#%%% 2.1. Open z array
original_data = xr.open_dataset(f'../Data/Input_data/Z500_{year_file_i}_{year_file_f}_{region}_{data_type}.nc')

#### Cut the subset if needed
if use_subset == 2:
    original_data = original_data.sel(time=slice(str(year_i), str(year_f)))

#### Invert the data array if needed and open the array
if region == 'NH':
    original_array = original_data.z.values
elif region == 'SH':
    original_array = original_data.z.values[:,::-1,:]

#%%% 2.3. Lat and Lon and years
if region == 'NH':
    lat = original_data.latitude.values
elif region == 'SH':
    lat = -original_data.latitude.values[::-1]
lat_rad = np.deg2rad(lat) # Factor to use later

closest_lat_polar = np.argmin(abs(lat-lat_polar_circle))

lon = original_data.longitude.values
#### Define opposite longitude for calculation in the pole region
lon_opposite_arr = [find_nearest(lon, lon[loni]-180, 'index') if lon[loni] >= 0 else find_nearest(lon, lon[loni]+180, 'index') for loni in range(len(lon))]

time = original_data.time.values


#%% 3. Compute area matrix for weighted average
area = area_matrix(lon, lat, res)


#%% 4. Daily data retrieval and save
data_dict = {}
struct_array_tosave = np.zeros((len(time),len(lat),len(lon)))
for day_n in tqdm(range(len(time))):

    #### First n_days_before days do not count since there is no LATmin data
    if day_n < n_days_before:
        struct_info = {}

    #### Full analysis
    else:

        #### Retrieve the LATmin value for each day
        if horizontal_LATmin == 1:
            LATmin_day, mean_val = compute_static_LATmin(day_n)
        elif horizontal_LATmin == 2:
            # LATmin_day, mean_val = compute_moving_LATmin(day_n)
            LATmin_day, mean_val = compute_moving_LATmin(day_n)

        #### Retrieve daily data
        data_in_day = original_array[day_n]

        #### Retrieve the gradient data for each day
        GHG_lat_day, GHG_lon_day, GHGS_day, GHGS2_day, GHGN_day = gradients_function(day_n, LATmin_day, data_in_day)

        #### Compute and clean Jet regions
        Jet = np.sqrt(GHG_lat_day**2+GHG_lon_day**2)
        Jet[Jet < 20] = np.nan
        Jet[~np.isnan(Jet)] = 1
        if consider_polar == 1:
            limit_pos_pole = find_nearest(lat, 90-delta/2, 'index')
            Jet[:limit_pos_pole,:][np.sqrt(GHG_lon_day[:limit_pos_pole,:]**2) >= 20] = 1
        if consider_polar == 3:
            limit_pos_pole = find_nearest(lat, 90-delta/2, 'index')
            Jet[:limit_pos_pole,:] = np.nan

        #### Blocked locations
        Block = np.ones((len(lat),len(lon)))
        Block = np.where((GHGS_day > 0) & (GHGN_day < 0), 1, np.nan)

        #### Deal with polar structures (1 considers only the southern gradient as in Sousa et al., 2021, 2 considers the area between 90-delta and 85º but computes the gradient in the between region, 3 does not consider pole regions)
        if consider_polar == 1:
            limit_pos_pole = find_nearest(lat, 90-delta, 'index')
            Block[:limit_pos_pole,:][GHGS_day[:limit_pos_pole,:] > 0] = 1
            Block[0] = np.nan
        elif consider_polar == 2:
            limit_pos_pole = find_nearest(lat, 85, 'index')
            Block[:limit_pos_pole,:] = np.nan
        elif consider_polar == 3:
            limit_pos_pole = find_nearest(lat, 90-delta, 'index')
            Block[:limit_pos_pole,:] = np.nan

        #### Clean where jet affects
        Block[Jet == 1] = np.nan

        #### Clean ridge pixel locations
        data_anom = (data_in_day - mean_val)
        data_anom[data_anom < 0] = np.nan

        #### Clean "holes" in data with a binary dilation algorithm, this dilation has one degree of freedom (i.e., only immediately touching pixels are accounted, no diagonals)
        #### FLAG: THIS MAY DIFFER SIGNIFICANTLY FROM PEDRO's APPROACH
        new_data_anom = clean_holes(data_anom, 'excess')
        new_data_anom = np.where(np.isnan(new_data_anom),1,np.nan)

        Ridge = np.copy(new_data_anom)

        LATmin_day_locs = [find_nearest(lat, k, 'index') for k in LATmin_day]
        for lon_loc in range(len(lon)):
            Ridge[np.array(LATmin_day_locs[lon_loc]):,lon_loc] = np.nan
        Ridge[~np.isnan(Ridge)] = 1
        Ridge[Jet == 1] = np.nan
        Ridge[Block == 1] = np.nan

        #### Rex pixel identifications
        Rex = np.where((Block == 1) & (GHGS_day - GHGS2_day > 20), 1, np.nan)

        #### Omega pixel identification
        Omega = np.where((Block == 1) & (GHGS_day - GHGS2_day <= 20), 1, np.nan)

        #### Change struct array to save pixel types
        struct_array = np.zeros((len(lat),len(lon)), dtype = int)
        struct_array[Ridge == 1] = 1
        struct_array[Rex == 1] = 2
        struct_array[Omega == 1] = 3

        #### Check and save structures that obey the daily rules
        #### Identify continuous structures
        bin_struct = np.copy(struct_array)
        bin_struct[bin_struct != 0] = 1   #### Create binary array for the binary dilation algorithm

        #### FLAG: THIS MAY DIFFER SIGNIFICANTLY FROM PEDRO's APPROACH (Not sure how they grouped the structures)
        iso_struct = identify_structures(bin_struct) #### Binary dilation algorithm to identify structures, with two degree of freedoms (pixels touching in the diagonal count as the same structure)

        #### Initialize the masks array to be used in the analysis loop
        antimeridian_structs = np.zeros(np.shape(iso_struct))
        antimeridian_type = np.zeros(np.shape(iso_struct))

        #### Number of structures in the day
        total_structures = np.max(iso_struct)

        #### Initialize the dictionary that will store the data in that day
        info = {}  # Keep center type, and mask
        info['Struct_array'] = np.zeros(np.shape(iso_struct))

        #### Count the valid structures to be saved
        struct_num = 1

        #### Run all structures and filter between normal ones, antimeridians (those that touch the boundary), and polar ones
        for num in range(1,total_structures+1):

            #### Individual structure mask (either one or NaN)
            indv = np.zeros(np.shape(iso_struct))
            indv[iso_struct == num] = 1
            indv[iso_struct != num] = np.nan

            #### Verify if structure is touching the antimeridian
            low_lon = min(np.where(indv == 1)[1]); high_lon = max(np.where(indv == 1)[1])
            if low_lon == 0 or high_lon == len(lon)-1:
                antimeridian_structs[indv == 1] = num
                #### Carry the structure type to the new antimeridian array
                antimeridian_type[indv*struct_array == 1] = 1
                antimeridian_type[indv*struct_array == 2] = 2
                antimeridian_type[indv*struct_array == 3] = 3
                continue

            #### If structure is neither on pole or antimeridian, store its info
            struct_area = int(np.nansum(indv*area))
            if struct_area >= min_struct_area:
                indv_iso_struct = indv*struct_array

                #### Latitudinally extend rex region
                rex_longs = np.unique(np.where(indv_iso_struct == 2)[1])
                indv_coords = np.where(indv == 1)
                indv_lat_tofill = []; indv_lon_tofill = []
                for i,j in zip(indv_coords[0], indv_coords[1]):
                    if j in rex_longs:
                        indv_lat_tofill.append(i)
                        indv_lon_tofill.append(j)
                indv_iso_struct[indv_lat_tofill, indv_lon_tofill] = 2

                #### Check if Ridge or block and which type
                indv_block = np.zeros(np.shape(indv_iso_struct))
                indv_block[(indv_iso_struct == 2) | (indv_iso_struct == 3)] = 1
                blocked_area = int(np.sum(indv_block*area))

                indv_omega = np.zeros(np.shape(indv_iso_struct))
                indv_omega[(indv_iso_struct == 3)] = 1
                omega_area = int(np.sum(indv_omega*area))

                #### Check distances to LATmin (as LATmin varies with lon) and the max latitudinal extension of structure
                high_diff_latmin = []; low_diff_latmin = []
                highest_lat = []; lowest_lat = []
                for i in range(np.shape(indv)[1]):
                    where_data = np.where(~np.isnan(indv[:,i]))[0]
                    if len(where_data) != 0:
                        high_diff_latmin.append(lat[np.min(where_data)]-LATmin_day[i])
                        low_diff_latmin.append(lat[np.max(where_data)]-LATmin_day[i])

                        highest_lat.append(lat[np.min(where_data)])
                        lowest_lat.append(lat[np.max(where_data)])

                high_diff_latmin = np.max(high_diff_latmin)
                low_diff_latmin = np.min(low_diff_latmin)

                max_ext = np.max(highest_lat)-np.min(lowest_lat)

                ### Check for the individual charachteristcs as set by Sousa et al. (2021) with new modifications
                if blocked_area < min_struct_area or high_diff_latmin < 15:        #### Ridge if blocked area bellow the min struct area and the maximum struct lat does not exceed LATmin+15º
                    struct_type = 'Ridge'

                elif blocked_area >= min_struct_area and omega_area/blocked_area >= 0.5 and low_diff_latmin < 15:

                    if omega_hybrid_method == 1:        #### If omega is the largest then all structure is omega (Same as Sousa et al., 2021)
                        struct_type = 'Omega block'

                    elif omega_hybrid_method == 2:      #### Altered method, may consider hybrid in this step

                        if blocked_area-omega_area > min_struct_area:   #### If Rex is substantially present, then it is considered an hybrid
                            struct_type = 'Rex block (hybrid)'

                        else:
                            struct_type = 'Omega block'     #### If not, then it is omega

                else:
                    if low_diff_latmin == res:          #### Hybrid if Rex and connected to the subtropical belt
                        struct_type = 'Rex block (hybrid)'

                    else:                               #### If none of the other conditions were valid, then all that remains is that the structure is Rex (either Pure, Polar, or hybrid between the two)

                        high_lat = min(np.where(indv == 1)[0])   #### Compute the structures latitude closest to the pole

                        #### FLAG: THIS MAY DIFFER SIGNIFICANTLY FROM PEDRO's APPROACH (Not sure how pole structures were considered in the paper)
                        if lat[high_lat] > lat_polar_circle:            #### Check if any part of the structure is above the polar circle

                            indv_above_polar = np.copy(indv)
                            indv_above_polar[closest_lat_polar+1:,:] = np.nan    #### Structure area poleward of the polar circle

                            #### If the area poleward of the polar circle is the majority and substantial (above the minimum threshold), then it is considered a polar block
                            if int(np.nansum(indv_above_polar*area)) > min_struct_area and int(np.nansum(indv_above_polar*area))/struct_area > 0.5:
                                struct_type = 'Rex block (polar)'

                            else:                            #### If the area poleward of the polar circle is not substantial then it is a regular Rex block
                                struct_type = 'Rex block'

                        else:                                #### If there is no part poleward of the polar circle then it is a regular Rex block
                            struct_type = 'Rex block'

                #### Save struct type and mask, update the number of valid structures
                info[f'{struct_num}'] = struct_type
                info['Struct_array'][indv == 1] = struct_num
                struct_num += 1

        #### Analyse antimeridian structures (Same process as before but use arrays cut in half and joined at the antimeridian for a cyclic boundary)
        #### FLAG: THIS MAY DIFFER SIGNIFICANTLY FROM PEDRO's APPROACH (Unclear if the structures considered a cyclic boundary, but makes sense to consider this)
        antim_structs_cut = np.append(antimeridian_structs[:,np.shape(antimeridian_structs)[1]//2:],
                                      antimeridian_structs[:,:np.shape(antimeridian_structs)[1]//2],
                                      axis = 1)
        antim_iso_struct = identify_structures(antim_structs_cut)

        antim_Structure = np.append(struct_array[:,np.shape(struct_array)[1]//2:],
                                    struct_array[:,:np.shape(struct_array)[1]//2],
                                    axis = 1)

        antim_total_structures = np.max(antim_iso_struct)

        for num in range(1,antim_total_structures+1):
            indv = np.zeros(np.shape(antim_iso_struct))
            indv[antim_iso_struct == num] = 1
            indv[antim_iso_struct != num] = np.nan

            #### Check distances to LATmin (as LATmin varies with lon) and the max latitudinal extension of structure
            high_diff_latmin = []; low_diff_latmin = []
            highest_lat = []; lowest_lat = []
            for i in range(np.shape(indv)[1]):
                where_data = np.where(~np.isnan(indv[:,i]))[0]
                if len(where_data) != 0:
                    high_diff_latmin.append(lat[np.min(where_data)]-LATmin_day[i])
                    low_diff_latmin.append(lat[np.max(where_data)]-LATmin_day[i])

                    highest_lat.append(lat[np.min(where_data)])
                    lowest_lat.append(lat[np.max(where_data)])

            high_diff_latmin = np.max(high_diff_latmin)
            low_diff_latmin = np.min(low_diff_latmin)

            max_ext = np.max(highest_lat)-np.min(lowest_lat)

            #### If structure has bigger area than threshold, store its info
            struct_area = int(np.nansum(indv*area))
            if struct_area >= min_struct_area:
                indv_iso_struct = indv*antim_Structure

                #### Latitudinally extend rex region
                rex_longs = np.unique(np.where(indv_iso_struct == 2)[1])
                indv_coords = np.where(indv == 1)
                indv_lat_tofill = []; indv_lon_tofill = []
                for i,j in zip(indv_coords[0], indv_coords[1]):
                    if j in rex_longs:
                        indv_lat_tofill.append(i)
                        indv_lon_tofill.append(j)
                indv_iso_struct[indv_lat_tofill, indv_lon_tofill] = 2

                #### Check if Ridge or block and which type
                indv_block = np.zeros(np.shape(indv_iso_struct))
                indv_block[(indv_iso_struct == 2) | (indv_iso_struct == 3)] = 1
                blocked_area = int(np.sum(indv_block*area))

                indv_omega = np.zeros(np.shape(indv_iso_struct))
                indv_omega[(indv_iso_struct == 3)] = 1
                omega_area = int(np.sum(indv_omega*area))

                ### Check for the individual charachteristcs as set by Sousa et al. (2021) with new modifications
                if blocked_area < min_struct_area or high_diff_latmin < 15:        #### Ridge if blocked area bellow the min struct area and the maximum struct lat does not exceed LATmin+15º
                    struct_type = 'Ridge'

                elif blocked_area >= min_struct_area and omega_area/blocked_area >= 0.5 and low_diff_latmin < 15:
                    if omega_hybrid_method == 1:        #### If omega is the largest then all structure is omega (Same as Sousa et al., 2021)
                        struct_type = 'Omega block'

                    elif omega_hybrid_method == 2:      #### Altered method, may consider hybrid in this step

                        if blocked_area-omega_area > min_struct_area:   #### If Rex is substantially present, then it is considered an hybrid
                            struct_type = 'Rex block (hybrid)'

                        else:
                            struct_type = 'Omega block'     #### If not, then it is omega

                else:
                    if low_diff_latmin == res:          #### Hybrid if Rex and connected to the subtropical belt
                        struct_type = 'Rex block (hybrid)'

                    else:                               #### If none of the other conditions were valid, then all that remains is that the structure is Rex (either Pure, Polar, or hybrid between the two)
                        high_lat = min(np.where(indv == 1)[0])   #### Compute the structures latitude closest to the pole

                        #### FLAG: THIS MAY DIFFER SIGNIFICANTLY FROM PEDRO's APPROACH (Not sure how pole structures were considered in the paper)
                        if lat[high_lat] > lat_polar_circle:            #### Check if any part of the structure is above the polar circle

                            indv_above_polar = np.copy(indv)
                            indv_above_polar[closest_lat_polar+1:,:] = np.nan     #### Structure area poleward of the polar circle

                            #### If the area poleward of the polar circle is the majority and substantial (above the minimum threshold), then it is considered a polar block
                            if int(np.nansum(indv_above_polar*area)) > min_struct_area and int(np.nansum(indv_above_polar*area))/struct_area > 0.5:
                                struct_type = 'Rex block (polar)'

                            else:                            #### If the area poleward of the polar circle is not substantial then it is a regular Rex block
                                struct_type = 'Rex block'

                        else:      #### If there is no part poleward of the polar circle then it is a regular Rex block
                            struct_type = 'Rex block'

                indv_to_save = np.append(indv[:,np.shape(indv)[1]//2:],
                                          indv[:,:np.shape(indv)[1]//2],
                                          axis = 1)

                info[f'{struct_num}'] = struct_type
                info['Struct_array'][indv_to_save == 1] = struct_num
                struct_num += 1

        #########################
        ##### Copy output
        struct_info = info.copy()

        #### Save structure data in each day
        struct_array_tosave[day_n] = struct_info['Struct_array']

        #### Delete structure array to save space
        del struct_info['Struct_array']

    #### Save information about each structure in a dictionary
    data_dict[time[day_n].astype(str)[:10]] = struct_info


#%% 5. Save final netcdf
data = xr.Dataset(
    data_vars = dict(Structs=(['time', 'lat', 'lon'], struct_array_tosave.astype(np.float32))),
    coords = dict(time=time, lat=lat, lon=lon),
    attrs = dict(description='Structures inspired by Sousa et al 2021, altered by Miguel Lima',
                  units='Structure type in associated pkl file')
    )

data.to_netcdf(f'../Data/Output_data/01-StructMasks_{year_i}_{year_f}_{region}.nc')

#### Save pkl file with the dictionary to see the type of each structure
with open(f'../Data/Output_data/01-StructTypes_{year_i}_{year_f}_{region}.pkl', 'wb') as handle:
        pickle.dump(data_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
