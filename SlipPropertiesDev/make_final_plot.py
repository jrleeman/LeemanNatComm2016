import numpy as np
import matplotlib.pyplot as plt

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

# Tuple of experiments we'll consider for plotting even data from
experiments_with_event_data = ('p4342','p4343','p4344','p4345','p4346',
                               'p4347','p4348','p4350','p4351')

# Read those experiments into a dictionary of event data
experiment_event_data = dict()

for experiment in experiments_with_event_data:
    experiment_event_data[experiment] = load_events(experiment)

# Make the plot
# Setup figure and axes
# Generally plots is ~1.33x width to height (10,7.5 or 12,9)
fig = plt.figure(figsize=(12,9))
ax1 = plt.subplot(2,1,1)
ax2 = plt.subplot(2,2,3)
ax3 = plt.subplot(2,2,4)

# Set labels and tick sizes
ax1.set_xlabel(r'Load Point Displacement [$\mu m$]',fontsize=18)
ax1.set_ylabel(r'Stiffness [1/$\mu m$]',fontsize=18)
ax1.tick_params(axis='both', which='major', labelsize=16)

ax2.set_xlabel(r'Stiffness [1/$\mu m$]',fontsize=18)
ax2.set_ylabel(r'Peak Slip Velocity [$\mu m/s$]',fontsize=18)
ax2.tick_params(axis='both', which='major', labelsize=16)

ax3.set_xlabel(r'Stiffness [1/$\mu m$]',fontsize=18)
ax3.set_ylabel(r'Slip Duration [$s$]',fontsize=18)
ax3.tick_params(axis='both', which='major', labelsize=16)

# Turns off chart clutter

# Turn off top and right tick marks
ax1.get_xaxis().tick_bottom()
ax1.get_yaxis().tick_left()

ax2.get_xaxis().tick_bottom()
ax2.get_yaxis().tick_left()

ax3.get_xaxis().tick_bottom()
ax3.get_yaxis().tick_left()

# Turn off top and right splines
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)

ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)

ax3.spines["top"].set_visible(False)
ax3.spines["right"].set_visible(False)

# Plotting

# Make panel A of displacement/stiffness
low_color = 0.
high_color = 1000.
color_map = plt.get_cmap('jet')
marker_size = 40
marker_alpha=0.5
color_col=11

for key in experiment_event_data:
    event_data = experiment_event_data[key]
    sc = ax1.scatter(event_data[:,9]/1000.,event_data[:,5],c=event_data[:,color_col],s=marker_size,alpha=marker_alpha,vmin=low_color,vmax=high_color,cmap=color_map)

# cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
# plt.colorbar(sc,cax=cbar_ax)

ax1.set_ylim(0,0.0009)
ax1.set_xlim(8,52)

ax1.axvspan(40, 50, alpha=0.2, color='k', zorder=0)

# Panel B

filter_col = 9
low_val = 40000.
high_val = 50000.
y_col = 11
marker_alpha = 0.3
for key in experiment_event_data:
    event_data = experiment_event_data[key]
    event_data = filter(event_data,filter_col,low_val,high_val)

    ax2.scatter(event_data[:,5],event_data[:,y_col],color='k',alpha=marker_alpha)
    ax2.errorbar(np.mean(event_data[:,5]),np.mean(event_data[:,y_col]),fmt='ro',ecolor='k',elinewidth=2,xerr=np.std(event_data[:,5]),yerr=np.std(event_data[:,y_col]))

ax2.set_ylim(0,900)
ax2.set_xlim(0.00045,0.00085)
ax2.set_xticks([0.00045,0.00055,0.00065,0.00075,0.00085])

# Panel C

filter_col = 9
low_val = 40000.
high_val = 50000.
y_col = 4
marker_alpha = 0.3
for key in experiment_event_data:
    event_data = experiment_event_data[key]
    event_data = filter(event_data,filter_col,low_val,high_val)

    ax3.scatter(event_data[:,5],event_data[:,y_col],color='k',alpha=marker_alpha)
    ax3.errorbar(np.mean(event_data[:,5]),np.mean(event_data[:,y_col]),fmt='ro',ecolor='k',elinewidth=2,xerr=np.std(event_data[:,5]),yerr=np.std(event_data[:,y_col]))

ax3.set_ylim(0,1.2)
ax3.set_xlim(0.00045,0.00085)
ax3.set_xticks([0.00045,0.00055,0.00065,0.00075,0.00085])

plt.show()
