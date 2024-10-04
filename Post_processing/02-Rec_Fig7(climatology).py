#%% 0. Start
import numpy as np
import xarray as xr
import pandas as pd
import pickle
from tqdm import tqdm
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature
import matplotlib.pyplot as plt
import matplotlib.path as mpath


#%% 1. Functions
#%% 1.1. Find nearest value in array
def find_nearest(array, value, which):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    if which == 'value':
        return array[idx]
    elif which == 'index':
        return idx


#%% 2. Open data
year_i = 1950; year_f = 2020
region = 'NH'

#%%% 2.1. Open tracked structures data
obs_NH = pd.read_csv('../Data/Output_data/03-Blocking_daily_catalogue_2019_2020_NH.csv')
masks_NH = xr.open_dataset('../Data/Output_data/03-CatalogueMasks_2019_2020_NH.nc')


#%%% 2.3. Lat, Lon and time
lat = masks_NH.lat.values
lon = masks_NH.lon.values
lons, lats = np.meshgrid(lon,lat)

time = pd.DatetimeIndex(masks_NH.time.values)

seasons = {'winter': [1,2,12],
           'spring': [3,4,5],
           'summer': [6,7,8],
           'fall': [9,10,11]}


#%% 3. Function to compute climatology per year or season
def give_climatology(season_to_eval, masks, masks_data):

    Ridge_days = np.zeros(np.shape(lons))
    Omega_days = np.zeros(np.shape(lons))
    Rex_days = np.zeros(np.shape(lons))
    Rex_pure = np.zeros(np.shape(lons))
    Rex_hybrid = np.zeros(np.shape(lons))
    Rex_polar = np.zeros(np.shape(lons))
    
    indexes = np.where(time.month.isin(seasons[season_to_eval]))
    
    dates_to_an = time.astype(str).str[:10].values[indexes]
    
    n_days = 0
    for data_string0 in tqdm(dates_to_an):
        
        n_days += 1
        day_n = np.where(time == data_string0)[0][0]
        
        if day_n < 15:
            continue
        else:
            mask_day = masks[day_n]
            
            for strct_id in np.unique(mask_day):
                
                if strct_id != 0:
                    
                    strct_type = masks_data[(masks_data.SID == strct_id) &
                                            (masks_data.YEAR == int(data_string0[:4])) &
                                            (masks_data.MONTH == int(data_string0[5:7])) &
                                            (masks_data.DAY == int(data_string0[8:]))].TYPE.values[0]
                    
                    daily_struct_array = np.zeros(np.shape(lons))
                    daily_struct_array[mask_day == strct_id] = 1
                    
                    if strct_type == 'Ridge':
                        Ridge_days[daily_struct_array == 1] += 1
                    elif strct_type == 'Omega block':
                        Omega_days[daily_struct_array == 1] += 1
                    elif strct_type == 'Rex block (hybrid)':
                        Rex_hybrid[daily_struct_array == 1] += 1
                        Omega_days[daily_struct_array == 1] += 1 
                    elif strct_type == 'Rex block':
                        Rex_pure[daily_struct_array == 1] += 1
                        Rex_days[daily_struct_array == 1] += 1
                    elif strct_type == 'Rex block (polar)':
                        Rex_polar[daily_struct_array == 1] += 1
                        Rex_days[daily_struct_array == 1] += 1
        
        season_info = {'Ridge': Ridge_days,
                        'Omega': Omega_days,
                        'Rex': Rex_days,
                        'Rex (pure)': Rex_pure,
                        'Rex (hybrid)': Rex_hybrid,
                        'Rex (polar)': Rex_polar}
        
    return season_info, n_days


#%% 4. Retrieve climatologies per season
winter_info, n_days_winter = give_climatology('winter', masks_NH.Structs.values, obs_NH)
spring_info, n_days_spring = give_climatology('spring', masks_NH.Structs.values, obs_NH)
summer_info, n_days_summer = give_climatology('summer', masks_NH.Structs.values, obs_NH)
fall_info, n_days_fall = give_climatology('fall', masks_NH.Structs.values, obs_NH)


#%% 5. Figure
#%%% 5.1. Function to draw hemispheric maps
def plt_hms(fig, lines, columns, which, to_draw, n_days, alinea, title, title2, title2str):
    # ax for Northern Hemisphere
    ax = fig.add_subplot(lines, columns, which, projection=ccrs.NorthPolarStereo())
    ax.set_extent([-180, 180, 90, 35], ccrs.PlateCarree())

    theta = np.linspace(0, 2*np.pi, 100)
    map_circle = mpath.Path(np.vstack([np.sin(theta), np.cos(theta)]).T * 0.5 + [0.5, 0.5])

    ax.set_boundary(map_circle, transform=ax.transAxes)

    coast = NaturalEarthFeature(category='physical', scale='110m',
                                facecolor='none', name='coastline')
    ax.add_feature(coast, edgecolor='k', linewidth=0.9, alpha=0.9, zorder=2)

    land = NaturalEarthFeature(category='physical', scale='110m',
                                facecolor='lightgray', name='land')
    ax.add_feature(land, edgecolor='k', linewidth=0.9, alpha=0.9, zorder=0)

    ax.gridlines()

    # Draw Ridge location
    vmin=0; vmax=35; step=1
    to_plot = to_draw/n_days*100

    to_plot[to_plot < 2] = np.nan

    im = ax.contourf(lon, lat, to_plot, levels = np.arange(vmin, vmax+step, step),
                      cmap = 'hot_r', alpha = 0.8, transform=ccrs.PlateCarree(), zorder=1, extend = 'max')
    ax.contour(lon, lat, to_plot, levels = [2,5,10,20,30],
                alpha = 0.8, transform=ccrs.PlateCarree(), zorder=2, colors='k', linewidths = 0.7)

    ##
    ax.set_title(title, fontsize=18, pad = 10)
    if title2 == 'yes':
        ax.annotate(title2str, (-0.35, 0.5), xycoords='axes fraction', fontsize = 24,
                    transform=ccrs.PlateCarree(), rotation=90, va='center', ma='center')

    ax.annotate(alinea, (0.005, 0.95), xycoords='axes fraction', fontsize = 22,
                transform=ccrs.PlateCarree())

    return ax, im

#%%% 5.2. Make all figure
plt.close('all')
fig = plt.figure(figsize=(16, 12))

################ Ridges
############ Ax1
ax1, im = plt_hms(fig, 3, 4, 1, winter_info['Ridge'], n_days_winter, 'a)', 'DJF', 'yes', 'Ridges')

############ Ax2
ax2, im = plt_hms(fig, 3, 4, 2, spring_info['Ridge'], n_days_spring, 'b)', 'MAM', 'no', 'Ridges')

############ Ax3
ax3, im = plt_hms(fig, 3, 4, 3, summer_info['Ridge'], n_days_summer, 'c)', 'JJA', 'no', 'Ridges')

############ Ax4
ax4, im = plt_hms(fig, 3, 4, 4, fall_info['Ridge'], n_days_fall, 'd)', 'SON', 'no', 'Ridges')

################ Omega
############ Ax5
ax5, im = plt_hms(fig, 3, 4, 5, winter_info['Omega'], n_days_winter, 'e)', 'DJF', 'yes', 'Omega\nblocks')

############ Ax6
ax6, im = plt_hms(fig, 3, 4, 6, spring_info['Omega'], n_days_spring, 'f)', 'MAM', 'no', 'Omega')

############ Ax7
ax7, im = plt_hms(fig, 3, 4, 7, summer_info['Omega'], n_days_summer, 'g)', 'JJA', 'no', 'Omega')

############ Ax8
ax8, im = plt_hms(fig, 3, 4, 8, fall_info['Omega'], n_days_fall, 'h)', 'SON', 'no', 'Omega')

################ Rex
############ Ax5
ax9, im = plt_hms(fig, 3, 4, 9, winter_info['Rex'], n_days_winter, 'i)', 'DJF', 'yes', 'Rex\nblocks')

############ Ax6
ax10, im = plt_hms(fig, 3, 4, 10, spring_info['Rex'], n_days_spring, 'j)', 'MAM', 'no', 'Rex')

############ Ax7
ax11, im = plt_hms(fig, 3, 4, 11, summer_info['Rex'], n_days_summer, 'k)', 'JJA', 'no', 'Rex')

############ Ax8
ax12, im = plt_hms(fig, 3, 4, 12, fall_info['Rex'], n_days_fall, 'l)', 'SON', 'no', 'Rex')

################ Colorbar
cax = fig.add_axes([0.2, 0.075, 0.62, 0.02])
cbar = fig.colorbar(im, cax = cax, ticks = np.arange(0, 40, 5),
                    orientation='horizontal', extend = 'max', pad = 0.45, aspect = 25)
cbar.set_label('% of days in season', fontsize=22, labelpad=8)
cbar.ax.tick_params(labelsize=20)

################ Complete
plt.subplots_adjust(wspace=0.08)
plt.savefig('../Figures/02-Rec_Fig7(climatology).jpg', dpi=300, bbox_inches = 'tight')