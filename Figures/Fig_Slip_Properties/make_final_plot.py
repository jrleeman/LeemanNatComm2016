import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Rectangle

def load_event_properties(experiment):
    """
    Load event property file picks for a given experiment number and return
    that data as an array
    """
    return np.loadtxt('%s_event_properties.txt'%experiment,delimiter=',',skiprows=1)

def load_blacklist(experiment):
    """
    Load event numbers from the blacklist file for each experiment and
    return them as an array
    """
    blacklist = np.loadtxt('%s_blacklist.txt'%experiment)
    return blacklist

def load_events(experiment):
    """
    Loads all events from a given experiment that are not on the blacklist
    file for that experiment. Returns array of event properties.
    """
    event_properties = load_event_properties(experiment)
    blacklist = load_blacklist(experiment)
    return np.delete(event_properties,blacklist,axis=0)

def filter(data,col,low,high):
    """
    Take array, filter out rows in which the element in the given column
    is not in the range low-high (inclusive)
    """
    inds = np.where(data[:,col]>=low)
    data_trim = data[inds]
    inds = np.where(data_trim[:,col]<=high)
    data_trim = data_trim[inds]
    return data_trim

def get_kc(disp):
    slope = (7.e-4-2.6e-6)/(16.-6.)
    if disp >= 16:
        return 7e-4
    else:
        return slope*disp - 0.0004

# Tuple of experiments we'll consider for plotting even data from
experiments_with_event_data = ('p4342','p4343','p4344','p4345','p4346',
                               'p4347','p4348','p4350','p4351')

# Tuple of experiments we'll plot unload/reload stiffness from
experiments_with_unload_reload = ('p4267','p4268','p4269','p4270','p4271',
                                  'p4272','p4273','p4309','p4310','p4311',
                                  'p4312','p4313','p4314','p4316','p4317',
                                  'p4327','p4328','p4329','p4330')

# Read those experiments into a dictionary of event data
experiment_event_data = dict()

for experiment in experiments_with_event_data:
    experiment_event_data[experiment] = load_events(experiment)

# Make the plot
# Setup figure and axes
# Generally plots is ~1.33x width to height (10,7.5 or 12,9)
fig = plt.figure(figsize=(12,12))
ax1 = plt.subplot(211)
ax2 = plt.subplot(212)


#
# Top Plot
#

exps = ['p4267','p4268','p4269','p4270','p4271','p4272','p4273',
        'p4309','p4310','p4311','p4312','p4313','p4314','p4316','p4317',
        'p4327','p4328','p4329','p4330']

ax1.text(-0.1,0.9,'A',transform = ax1.transAxes,fontsize=24)

# Set labels and tick sizes
ax1.set_xlabel(r'Average LP Displacement [mm]',fontsize=18)
ax1.set_ylabel(r'Stiffness [1/um]x1000',fontsize=18)
ax1.tick_params(axis='both', which='major', labelsize=16)

# Turns off chart clutter

# Turn off top and right tick marks
ax1.get_xaxis().tick_bottom()
ax1.get_yaxis().tick_left()

# Turn off top and right splines
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)

# Plotting

for exp in exps:

    df = pd.read_csv('/Users/jleeman/Dropbox/PennState/BiaxExperiments/%s/%s_stiffness_cycles.txt'%(exp,exp))

    temp = df[df['Behavior']=='stable']
    ax1.scatter(temp['AvgDisp']/1000.,temp['Slope']*1000,color='k',s=50,alpha=0.6)

    #temp = df[df['Behavior']=='slow']
    #ax1.scatter(temp['AvgDisp']/1000.,temp['Slope'],color='r',s=50,alpha=0.6)

    #temp = df[df['Behavior']=='fast']
    #ax1.scatter(temp['AvgDisp']/1000.,temp['Slope'],color='r',s=50,alpha=0.6)

# Add rectangle for where figure B comes from
rect_x1 = 10.
rect_x2 = 50.
rect_y1 = 0.
rect_y2 = 0.0009*1000
rect_width = rect_x2-rect_x1
rect_height = rect_y2-rect_y1
ax1.add_patch(Rectangle((rect_x1,rect_y1),rect_width,rect_height,alpha=0.2, zorder=0,facecolor="k"))

# Set limits
ax1.set_xlim(0,50)
ax1.set_ylim(0,0.004*1000)

# Plot Kc
df = pd.read_excel('/Users/jleeman/Dropbox/PennState/BiaxExperiments/p4309/p4309_rsf_fits.xlsx')




for i,fit in df.iterrows():

    if fit['Grade'] == 'A':
        #color='#000066'
        #color='#FFFFFF'
        color='#0000FF'
    elif fit['Grade'] == 'B':
        color='#0066CC'
        color='#0000FF'
        #color='#FFFFFF'
    elif fit['Grade'] == 'C':
        #color='#00CCFF'
        color='#FFFFFF'
        continue
    elif fit['Grade'] == 'D':
        #color='#00FFFF'
        color='#FFFFFF'
        continue

    if fit['Type']=='Down' and fit['Law']=='r' and fit['k']==0.0055:
        ax1.scatter(fit['LP_Disp']/1000.,fit['Kc']*1000,c=color,s=60,marker='v',zorder=50)
        print fit['LP_Disp']/1000.,fit['Kc']

    elif fit['Type']=='Up' and fit['Law']=='r' and fit['k']==0.0055:
        ax1.scatter(fit['LP_Disp']/1000.,fit['Kc']*1000,c=color,s=60,marker='^',zorder=50)
        print fit['LP_Disp']/1000.,fit['Kc']
    else:
        pass


low_color = 10./1000.
high_color = 4000./1000.
color_map = plt.get_cmap('rainbow_r')
marker_size = 40
marker_alpha=0.5
color_col=11

for key in experiment_event_data:
    event_data = experiment_event_data[key]
    sc = ax1.scatter(event_data[:,9]/1000.,event_data[:,5]*1000,c=event_data[:,color_col]/1000.,s=marker_size,alpha=marker_alpha,vmin=low_color,vmax=high_color,cmap=color_map)
    print key,np.min(event_data[:,color_col]), np.max(event_data[:,color_col])

# Plot line for kc definition
ax1.plot([6,16,50],[2.6e-6*1000,7e-4*1000,7e-4*1000],color='k',linewidth=2)

# Set labels and tick sizes
ax2.set_xlabel(r'Load Point Displacement [$\mu m$]',fontsize=18,labelpad=15)
ax2.set_ylabel(r'Stiffness [1/$\mu m$]',fontsize=18)
ax2.tick_params(axis='both', which='major', labelsize=16)

# Turns off chart clutter

# Turn off top and right tick marks
ax2.get_xaxis().tick_bottom()
ax2.get_yaxis().tick_left()

ax2.get_xaxis().tick_bottom()
ax2.get_yaxis().tick_left()

# Turn off top and right splines
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)

# Plotting

# Make panel A of displacement/stiffness
ax2.text(-0.1,0.9,'B',transform = ax2.transAxes,fontsize=24)

low_color = 10./1000.
high_color = 4000./1000.
color_map = plt.get_cmap('rainbow_r')
marker_size = 40
marker_alpha=0.5
color_col=11

for key in experiment_event_data:
    event_data = experiment_event_data[key]

    k_kc_ratio = []
    for k,disp in zip(event_data[:,5],event_data[:,9]/1000.):
        k_kc_ratio.append(k/get_kc(disp))

    sc = ax2.scatter(event_data[:,9]/1000.,k_kc_ratio,c=event_data[:,color_col]/1000.,s=marker_size,alpha=marker_alpha,vmin=low_color,vmax=high_color,cmap=color_map)
    print key,np.min(event_data[:,color_col]), np.max(event_data[:,color_col])
# cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
# plt.colorbar(sc,cax=cbar_ax)
# for experiment in experiments_with_unload_reload:
#     df = pd.read_csv('/Users/jleeman/Dropbox/PennState/BiaxExperiments/%s/%s_stiffness_cycles.txt'%(experiment,experiment))
#
#     ax2.scatter(df['AvgDisp']/1000.,df['Slope'],color='g',s=50,alpha=0.6)

position=fig.add_axes([0.37,0.16,0.5,0.02])  ## the parameters are the specified position you set [left, bottom, width, height]
cb = fig.colorbar(sc,cax=position,orientation='horizontal')
cb.solids.set_edgecolor("face")
cb.set_label(r'Peak Slip Velocity [$mm/s$]',fontsize=14)

ax2.set_ylim(0,1.4)
ax2.set_xlim(8,52)

ax2.axvspan(40, 50, alpha=0.2, color='k', zorder=0)

plt.savefig('figure.png', bbox_inches="tight");
