#%% 0. Start
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
from scipy.signal import butter,filtfilt

#%% 1. Open data
#%%% 1.1. Blocking catalogue
blocks_daily = pd.read_csv('../Data/Output_data/03-Blocking_daily_catalogue_2019_2020_NH.csv')


#%% 2. Retrieve data by structure type
#%%% 2.1. Create month day pair on original catalogue
blocks_daily['MONTH_DAY'] = [f'{str(mm).zfill(2)}-{str(dd).zfill(2)}' for mm, dd in zip(blocks_daily.MONTH.values, blocks_daily.DAY.values)]

#%%% 2.2. Function to order the dates and respective values
days_in_month = {1: 31,   # January
                 2: 28,   # February
                 3: 31,   # March
                 4: 30,   # April
                 5: 31,   # May
                 6: 30,   # June
                 7: 31,   # July
                 8: 31,   # August
                 9: 30,   # September
                 10: 31,  # October
                 11: 30,  # November
                 12: 31}  # December
dates = np.array([f'{str(month).zfill(2)}-{str(day).zfill(2)}' for month, max_day in days_in_month.items() for day in range(1, max_day + 1)])
def reorder(series):
    ordered_vals = []
    for date_temp in dates:
        try:
            ordered_vals.append(series.values[np.where(series.index == date_temp)[0][0]])
        except:
            ordered_vals.append(0)
    return np.array(ordered_vals)
    
#%%% 2.3. Get distributions
all_dist = reorder(blocks_daily.MONTH_DAY.value_counts())
ridge_dist = reorder(blocks_daily[blocks_daily.TYPE == 'Ridge'].MONTH_DAY.value_counts())
omega_dist = reorder(blocks_daily[blocks_daily.TYPE == 'Omega block'].MONTH_DAY.value_counts())
hybrid_dist = reorder(blocks_daily[blocks_daily.TYPE == 'Rex block (hybrid)'].MONTH_DAY.value_counts())
rex_dist = reorder(blocks_daily[blocks_daily.TYPE == 'Rex block'].MONTH_DAY.value_counts())

hybpolar_dist = reorder(blocks_daily[blocks_daily.TYPE == 'Rex block (hybrid polar)'].MONTH_DAY.value_counts())
polar_dist = reorder(blocks_daily[blocks_daily.TYPE == 'Rex block (polar)'].MONTH_DAY.value_counts())

polar_dist = polar_dist + hybpolar_dist


#%% 3. Figure
#%%% 3.1. Function to draw axes
def round_up_to_nearest(number, near):
    return int(math.ceil(number / near)) * near

def round_down_to_nearest(number, near):
    return int(math.floor(number / near)) * near

def butter_lowpass_filter(data, cutoff, fs, order):
    nyq=0.5*fs
    normal_cutoff = cutoff / nyq
    # Get the filter coefficients 
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    y = filtfilt(b, a, data)
    return y

def draw_axes(lines, columns, which, title, inlabel, data, color):
    # Filter requirements.
    T = 24*60*60        # Sample Period
    fs = 1/T            # sample rate, Hz
    cutoff = 1/(T*30)    # desired cutoff frequency of the filter, Hz, slightly higher than actual fs (Filter all variation below a month)
    order = 2       # sin wave can be approx represented as quadratic
    dataf = butter_lowpass_filter(data, cutoff, fs, order)
    
    ## Fig
    ax = fig.add_subplot(lines,columns,which)
    ax.grid(zorder=0, linestyle='--', alpha = 0.6, lw=0.5)

    ax.plot(np.arange(len(dates)), data, c=color, lw = 1)
    ax.plot(np.arange(len(dates)), dataf, color='black', lw=2, ls='--')

    ax.set_ylabel('# of events', fontsize = 15)
    ax.set_title(title, fontsize = 16, pad = 10)
    ax.set_xticks([0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334][1::2],
                   ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'][1::2], fontsize=13)

    top = round_up_to_nearest(np.max(data)+np.max(data)*0.05,10)
    bot = round_down_to_nearest(np.min(data)-np.min(data)*0.05,10)
    if bot < 0:
        bot = 0
    ax.set_ylim(bot, top); ax.set_xlim(0,365)
    
    ax.annotate(inlabel, (0.87, 0.88), xycoords='axes fraction', fontsize = 18)
    
    return ax

#%%% 3.2. Produce image
plt.close('all')
fig = plt.figure(figsize=(16, 9))

### All types
ax1 = draw_axes(2, 3, 1, 'All', 'a)', all_dist,'darkgray')

### Ridges
ax2 = draw_axes(2, 3, 2, 'Ridge', 'b)', ridge_dist, 'tab:red')

### Omega types
ax3 = draw_axes(2, 3, 3, 'Omega blocks', 'c)', omega_dist, 'limegreen')

### Hybrid types
ax4 = draw_axes(2, 3, 4, 'Hybrid blocks', 'd)', hybrid_dist, 'purple')

### Rex types
ax5 = draw_axes(2, 3, 5, 'Rex blocks', 'e)', rex_dist, 'tab:cyan')

### Polar types
ax6 = draw_axes(2, 3, 6, 'Polar blocks', 'f)', polar_dist, 'aqua')

############ Adjust
plt.subplots_adjust(wspace=0.23, hspace=0.27)

############ Save
plt.savefig('../Figures/04-Year_climatology_NH.jpg', dpi=300, bbox_inches = 'tight')