#%% 0. Start
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math


#%% 1. Open data
#%%% 1.1. Blocking catalogue
blocks_event = pd.read_csv('../Data/Output_data/03-Blocking_event_catalogue_2019_2020_NH.csv')


#%% 2. Prepare arrays to plot
#%%% 2.1. Histogram function
def histo(array):
    clean_arr = array[~np.isnan(array)]
    start = 0; end = 10; step = 0.5
    y, x = np.histogram(clean_arr, np.arange(start,end+step,step))
    return np.diff(x)+x[:-1]-step/2, y, np.nanmean(array)

#%%% 2.2. Create series
all_vals, all_dist, all_mean = histo(blocks_event.BI_MAX.values)
ridge_vals, ridge_dist, ridge_mean = histo(blocks_event[blocks_event.DOM_TYPE == 'Ridge'].BI_MAX.values)
omega_vals, omega_dist, omega_mean = histo(blocks_event[blocks_event.DOM_TYPE == 'Omega block'].BI_MAX.values)
hybrid_vals, hybrid_dist, hybrid_mean = histo(blocks_event[blocks_event.DOM_TYPE == 'Rex block (hybrid)'].BI_MAX.values)
rex_vals, rex_dist, rex_mean = histo(blocks_event[blocks_event.DOM_TYPE == 'Rex block'].BI_MAX.values)
polar_vals, polar_dist, polar_mean = histo(blocks_event[(blocks_event.DOM_TYPE == 'Rex block (polar)') |
                                                        (blocks_event.DOM_TYPE == 'Rex block (hybrid polar)')].BI_MAX.values)


#%% 3. Figure
#%%% 3.1. Function to draw axes
def round_up_to_nearest(number, near):
    return int(math.ceil(number / near)) * near

def round_down_to_nearest(number, near):
    return int(math.floor(number / near)) * near

def draw_axes(lines, columns, which, title, inlabel, mean, abcs, data, color):
    ax = fig.add_subplot(lines,columns,which)
    ax.grid(zorder=0, linestyle='--', alpha = 0.6, lw=0.5)

    ax.bar(abcs, data, width = 0.4, color=color)
    
    ax.set_ylabel('# of events', fontsize = 15)
    ax.set_title(title, fontsize = 16, pad = 10)
    ax.set_xticks(np.arange(0,11,1), fontsize=13)
    ax.set_xlim(min(abcs)-0.25,max(abcs))

    top = round_up_to_nearest(np.max(data)+np.max(data)*0.05,10)
    bot = 0
    ax.set_ylim(bot, top)
    
    ax.annotate(f'{inlabel}\n$\mu$={round(mean,1)}', (0.94, 0.76), xycoords='axes fraction', fontsize = 18,horizontalalignment='right')
    
    return ax

#%%% 3.2. Produce image
plt.close('all')
fig = plt.figure(figsize=(16, 9))

### All types
ax1 = draw_axes(2, 3, 1, 'All', 'a)', all_mean, all_vals, all_dist,'darkgray')

### Ridges
ax2 = draw_axes(2, 3, 2, 'Ridge', 'b)', ridge_mean, ridge_vals, ridge_dist, 'tab:red')

### Omega types
ax3 = draw_axes(2, 3, 3, 'Omega blocks', 'c)', omega_mean, omega_vals, omega_dist, 'limegreen')

### Hybrid types
ax4 = draw_axes(2, 3, 4, 'Hybrid blocks', 'd)', hybrid_mean, hybrid_vals, hybrid_dist, 'purple')

### Rex types
ax5 = draw_axes(2, 3, 5, 'Rex blocks', 'e)', rex_mean, rex_vals, rex_dist, 'tab:cyan')

### Polar types
ax6 = draw_axes(2, 3, 6, 'Polar blocks', 'f)', polar_mean, polar_vals, polar_dist, 'aqua')

############ Adjust
plt.subplots_adjust(wspace=0.27, hspace=0.27)

############ Save
plt.savefig('../Figures/06-Intensity_climatology_NH.jpg', dpi=300, bbox_inches = 'tight')