#%% 0. Start
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


#%% 1. Functions


#%% 2. Open data
#%%% 2.1. Blocking catalogue
blocks = pd.read_csv('../Data/Output_data/03-Blocking_event_catalogue_2019_2020_NH.csv')
nmr_years = 2

#%% 3. Retrieve info from each season individually
#%%% 3.1. Initialize
seasons = {'winter': [1,2,12],
           'spring': [3,4,5],
           'summer': [6,7,8],
           'fall': [9,10,11]}

#%%% 3.2. Function to iterate over each SID and save info
def stats_season(season_to_eval):
    
    #### Dataframe to evaluate by the prevalent type of structure
    blocks_inseason = blocks[blocks.MONTH_START.isin(seasons[season_to_eval])]
    SIDs_inseason = np.unique(blocks_inseason.SID.values)
    
    nmr_all = 0; dur_all = []
    nmr_ridge = 0; dur_ridge = []
    nmr_omega = 0; dur_omega = []
    nmr_rex = 0; dur_rex = []
    nmr_pole = 0; dur_pole = []
    
    for strct_id in SIDs_inseason:
        strct_to_eval = blocks_inseason[blocks_inseason.SID == strct_id]
        strct_type = np.unique(strct_to_eval.DOM_TYPE.values)
        strct_dur = np.unique(strct_to_eval.DURATION.values)
        
        nmr_all += 1
        dur_all.append(strct_dur)
        if strct_type == 'Ridge':
            nmr_ridge += 1; dur_ridge.append(strct_dur)
        elif strct_type == 'Omega block':
            nmr_omega += 1; dur_omega.append(strct_dur)
        elif strct_type == 'Rex block':
            nmr_rex += 1; dur_rex.append(strct_dur)
        elif strct_type == 'Rex block (hybrid)':
            nmr_rex += 1; dur_rex.append(strct_dur)
        elif strct_type == 'Rex block (polar hybrid)':
            nmr_pole += 1; dur_pole.append(strct_dur)
        elif strct_type == 'Rex block (polar)':
            nmr_pole += 1; dur_pole.append(strct_dur)
    
    info_dict = {'nmr_all': round(nmr_all/nmr_years,1), 'dur_all': round(np.mean(dur_all),1),
                 'nmr_ridge': round(nmr_ridge/nmr_years,1), 'dur_ridge': round(np.mean(dur_ridge),1),
                 'nmr_omega': round(nmr_omega/nmr_years,1), 'dur_omega': round(np.mean(dur_omega),1),
                 'nmr_rex': round(nmr_rex/nmr_years,1), 'dur_rex': round(np.mean(dur_rex),1),
                 'nmr_pole': round(nmr_pole/nmr_years,1), 'dur_pole': round(np.mean(dur_pole),1)}
    
    return info_dict

#%%% 3.3. Retrieve info for all seasons
info_winter = stats_season('winter')
info_spring = stats_season('spring')
info_summer = stats_season('summer')
info_fall = stats_season('fall')


#%% 4. Figure
seasons = ('DJF', 'MAM', 'JJA', 'SON')
info_to_plot = {'All': (info_winter['dur_all'], info_spring['dur_all'], info_summer['dur_all'], info_fall['dur_all']),
                'Ridge': (info_winter['dur_ridge'], info_spring['dur_ridge'], info_summer['dur_ridge'], info_fall['dur_ridge']),
                'Omega': (info_winter['dur_omega'], info_spring['dur_omega'], info_summer['dur_omega'], info_fall['dur_omega']),
                'Rex': (info_winter['dur_rex'], info_spring['dur_rex'], info_summer['dur_rex'], info_fall['dur_rex']),
                'Rex (polar)': (info_winter['dur_pole'], info_spring['dur_pole'], info_summer['dur_pole'], info_fall['dur_pole'])
                }
info_to_label = {'All': (info_winter['nmr_all'], info_spring['nmr_all'], info_summer['nmr_all'], info_fall['nmr_all']),
                'Ridge': (info_winter['nmr_ridge'], info_spring['nmr_ridge'], info_summer['nmr_ridge'], info_fall['nmr_ridge']),
                'Omega': (info_winter['nmr_omega'], info_spring['nmr_omega'], info_summer['nmr_omega'], info_fall['nmr_omega']),
                'Rex': (info_winter['nmr_rex'], info_spring['nmr_rex'], info_summer['nmr_rex'], info_fall['nmr_rex']),
                'Rex (polar)': (info_winter['nmr_pole'], info_spring['nmr_pole'], info_summer['nmr_pole'], info_fall['nmr_pole'])
                }

plt.close('all')
fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(1,1,1)

ax.grid(zorder=0, linestyle='--')

width = 0.15
multiplier = -1
x = np.arange(len(seasons))
colors = ['grey', 'tab:red', 'tab:green', 'tab:cyan', 'tab:pink']

#### All
for i, (attribute, measurement) in enumerate(info_to_plot.items()):
    offset = width * multiplier
    rects = ax.bar(x + offset, measurement, width, label=attribute, zorder = 3, color = colors[i])
    ax.bar_label(rects, labels=info_to_label[attribute], padding=4, fontsize=15)
    multiplier += 1

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Duration (days)', fontsize = 18)
ax.set_title('Events and Duration per dominant type (N.H.)', fontsize = 18, pad = 10)
ax.set_xticks(x + width, seasons, fontsize=16)
ax.set_yticks(np.arange(0,38,4), np.arange(0,38,4), fontsize=16)
ax.legend(loc='upper left', ncol=3, fontsize = 16)
ax.set_ylim(0, 38)
ax.annotate('Events\nper season', (2.9, 9.4),
        xytext=(0.7, 0.88), textcoords='axes fraction',
        arrowprops=dict(facecolor='black', shrink=0.1, width = 1.5, headwidth = 10),
        fontsize=20,
        horizontalalignment='center', verticalalignment='top')

### Savefig
plt.savefig('../Figures/03-Rec_Fig9(bars).jpg', dpi=300, bbox_inches = 'tight')