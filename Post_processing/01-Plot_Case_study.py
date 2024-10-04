#%% 0. Start
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch


#%% 1. Functions
#%% 1.1. Find nearest value in array

#%% 2. Open data
#%%% 2.1. Open event catalogue
blocks_daily = pd.read_csv('../Data/Output_data/03-Blocking_daily_catalogue_2019_2020_NH.csv')

#%%% 2.2. Open blocking masks
masks_data = xr.open_dataset('../Data/Output_data/03-CatalogueMasks_2019_2020_NH.nc')
masks_array = masks_data.Structs.values

#%%% 2.3. Open Z500 data
original_data = xr.open_dataset('../Data/Input_data/Z500_2019_2020_NH_NCAR.nc')
original_array = original_data.z.values

#%%% 2.4. Lat, Lon and time
lat = masks_data.lat.values
lon = masks_data.lon.values
lons, lats = np.meshgrid(lon,lat)

time = pd.Series(masks_data.time.astype(str)).str[:10].values


#%% 3. Chose case studies
#%%%% 3.1. 2010 summer case
blocks_2019 = blocks_daily[blocks_daily.YEAR == 2019]
start_date_2019 = '2019-06-19'
end_date_2019 = '2019-07-28'


#%% 4. Plot case studies
#%%% 4.0. Colors by type
colors = ['red', 'limegreen', 'purple', 'aqua', 'blue']

#%%% 4.1. Eulerian case studies
#%%% 4.1.1. Function
def plt_hms(fig, lines, columns, which, date_index):
    date = time[date_index]
    data_in_day = masks_array[date_index]
    events = np.unique(data_in_day)
    arr_types = np.copy(data_in_day)
    
    if len(events) > 1:
        for ev in events:
            if ev != 0:
                obs = blocks_daily[(blocks_daily.SID == ev) &
                                    (blocks_daily.YEAR == int(date[:4])) &
                                    (blocks_daily.MONTH == int(date[5:7])) &
                                    (blocks_daily.DAY == int(date[8:10]))]
                stype = obs.TYPE.values[0]
                
                if stype == 'Ridge':
                    arr_types[arr_types == ev] = 1
                elif stype == 'Omega block':
                    arr_types[arr_types == ev] = 2
                elif stype == 'Rex block (hybrid)':
                    arr_types[arr_types == ev] = 3
                elif stype == 'Rex block':
                    arr_types[arr_types == ev] = 4
                elif stype == 'Rex block (polar)':
                    arr_types[arr_types == ev] = 5
    
    arr_types[arr_types == 0] = np.nan
    
    #### Start figure
    lambert = ccrs.Orthographic(central_longitude=0,central_latitude=45)
    lambert._threshold /= 100.0
    ax = fig.add_subplot(lines, columns, which, projection=lambert)
    ax.set_title(f'{date}', x = 1, y = -0.135, ha = 'right', fontsize = 20)
    ax.set_extent([-35, 40, 30, 80], ccrs.PlateCarree())
    
    # Add gridlines with labels
    gl = ax.gridlines(draw_labels=True, linestyle='--', lw=0.5, color='k', alpha=0.75)
    gl.bottom_labels = False
    gl.left_labels = False
    
    coast = NaturalEarthFeature(category='physical', scale='10m',
                                facecolor='none', name='coastline')
    ax.add_feature(coast, edgecolor='k', linewidth=0.9, alpha=0.9, zorder=2)
    
    land = NaturalEarthFeature(category='physical', scale='10m',
                                facecolor='lightgray', name='land')
    ax.add_feature(land, edgecolor='k', linewidth=0.9, alpha=0.9, zorder=0)
    borders = NaturalEarthFeature(category='cultural',
                                  name='admin_0_countries',
                                  scale='10m',
                                  facecolor='none')
    ax.add_feature(borders, edgecolor = 'k', linewidth=0.6, alpha=0.9, zorder=2)
    
    ####
    im = ax.pcolormesh(lon, lat, arr_types,
                        cmap = LinearSegmentedColormap.from_list('Blockings', colors, N=6),
                        vmin = 1, vmax = 5,
                        transform = ccrs.PlateCarree())
    
    ################ Colorbar
    legend_elements = [Patch(facecolor='red', edgecolor='k', label='Ridge'),
                        Patch(facecolor='limegreen', edgecolor='k', label='Omega'),
                        Patch(facecolor='purple', edgecolor='k', label='Rex (hybrid)'),
                        Patch(facecolor='aqua', edgecolor='k', label='Rex (pure)'),
                        Patch(facecolor='blue', edgecolor='k', label='Rex (polar)')]
    plt.legend(handles=legend_elements, loc=(0,-0.13), fontsize=12, ncol=3)
    
    #### Z500
    Z500_lvls = np.arange(4500,6000,20)
    CS = ax.contour(lon, lat, original_array[date_index],
                      levels=Z500_lvls, linewidths=np.where(Z500_lvls%100 == 0, 1.5, 0.2),
                      colors = 'black', transform = ccrs.PlateCarree())
    ax.clabel(CS, CS.levels[::5], inline = True, fontsize = 15)
    
    return ax, im

#%%% 4.1.2. 2010 warm case study
# date_index = np.where(time == start_date_2019)[0][0]
# date_in_study = time[date_index]

# while date_in_study != end_date_2019:
#     plt.close('all')
#     fig = plt.figure(figsize=(14, 7))
#     plt_hms(fig, 1, 1, 1, date_index)
#     plt.savefig(f'../Figures/01-Case_study/{time[date_index]}.jpg',dpi=300, bbox_inches = 'tight')
    
#     date_index += 1
#     date_in_study = time[date_index]

date_index = np.where(time == '2019-07-23')[0][0]
plt.close('all')
fig = plt.figure(figsize=(14, 7))
plt_hms(fig, 1, 1, 1, date_index)
plt.savefig('../Figures/01-Case_study.jpg',dpi=300, bbox_inches = 'tight')

date_index += 1
date_in_study = time[date_index]
