#%% 0. Start
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import scipy


#%% 1. Open data
#%%% 1.1. Blocking catalogue
blocks_event = pd.read_csv('../Data/Output_data/03-Blocking_event_catalogue_2019_2020_NH.csv')


#%% 2. Prepare arrays to plot
#%%% 2.1. Reorder function
def reorder(series):
    ordered_vals = []
    # ordered_fill_index = np.arange(4,series.index.max()+1)
    ordered_fill_index = np.arange(5,30)
    for indx in ordered_fill_index:
        if indx in series.index.values:
            ordered_vals.append(series.values[np.where(series.index == indx)[0][0]])
        else:
            ordered_vals.append(0)
    return ordered_fill_index, np.array(ordered_vals)

#%%% 2.2. Create series
all_vals, all_dist = reorder(blocks_event.DURATION.value_counts())
ridge_vals, ridge_dist = reorder(blocks_event[blocks_event.DOM_TYPE == 'Ridge'].DURATION.value_counts())
omega_vals, omega_dist = reorder(blocks_event[blocks_event.DOM_TYPE == 'Omega block'].DURATION.value_counts())
hybrid_vals, hybrid_dist = reorder(blocks_event[blocks_event.DOM_TYPE == 'Rex block (hybrid)'].DURATION.value_counts())
rex_vals, rex_dist = reorder(blocks_event[blocks_event.DOM_TYPE == 'Rex block'].DURATION.value_counts())

polar_vals, polar_dist = reorder(blocks_event[(blocks_event.DOM_TYPE == 'Rex block (polar)') |
                                              (blocks_event.DOM_TYPE == 'Rex block (hybrid polar)')].DURATION.value_counts())


#%% 3. Figure
#%%% 3.1. Function to draw axes
def round_up_to_nearest(number, near):
    return int(math.ceil(number / near)) * near

def round_down_to_nearest(number, near):
    return int(math.floor(number / near)) * near

def monoExp(x, m, t, b):
    return m * np.exp(-t * x) + b

def draw_axes(lines, columns, which, title, inlabel, abcs, data, color):
    ## Poly fit
    # params, cv = scipy.optimize.curve_fit(monoExp,abcs,data,(10000,0.3,30))
    # m, t, b = params
    
    ## Fig
    ax = fig.add_subplot(lines,columns,which)
    ax.grid(zorder=0, linestyle='--', alpha = 0.6, lw=0.5)

    ax.bar(abcs, data, width = 0.9, color=color)
    # ax.plot(abcs, monoExp(abcs, m, t, b), color='black', lw=1.5, ls='--')

    ax.set_ylabel('# of events', fontsize = 15)
    ax.set_title(title, fontsize = 16, pad = 10)
    ax.set_xticks(np.arange(0,30,5), fontsize=13)
    ax.set_xlim(min(abcs)-0.75,max(abcs))

    top = round_up_to_nearest(np.max(data)+np.max(data)*0.05,10)
    bot = 0
    ax.set_ylim(bot, top)
    
    ax.annotate(inlabel, (0.87, 0.88), xycoords='axes fraction', fontsize = 18)
    
    return ax, #params

#%%% 3.2. Produce image
plt.close('all')
fig = plt.figure(figsize=(16, 9))

### All types
ax1 = draw_axes(2, 3, 1, 'All', 'a)', all_vals, all_dist,'darkgray')

### Ridges
ax2 = draw_axes(2, 3, 2, 'Ridge', 'b)', ridge_vals, ridge_dist, 'tab:red')

### Omega types
ax3 = draw_axes(2, 3, 3, 'Omega blocks', 'c)', omega_vals, omega_dist, 'limegreen')

### Hybrid types
ax4 = draw_axes(2, 3, 4, 'Hybrid blocks', 'd)', hybrid_vals, hybrid_dist, 'purple')

### Rex types
ax5 = draw_axes(2, 3, 5, 'Rex blocks', 'e)', rex_vals, rex_dist, 'tab:cyan')

### Polar types
ax6 = draw_axes(2, 3, 6, 'Polar blocks', 'f)', polar_vals, polar_dist, 'aqua')

############ Adjust
plt.subplots_adjust(wspace=0.27, hspace=0.27)

############ Save
plt.savefig('../Figures/05-Duration_climatology_NH.jpg', dpi=300, bbox_inches = 'tight')