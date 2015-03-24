import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Rectangle

def load_event_properties(experiment):
    """
    Load event property file picks for a given experiment number and return
    that data as an array
    """
    return np.loadtxt('../Slip_Property_Data/%s_event_properties.txt'%experiment,delimiter=',',skiprows=1)

def load_blacklist(experiment):
    """
    Load event numbers from the blacklist file for each experiment and
    return them as an array
    """
    blacklist = np.loadtxt('../Slip_Property_Data/%s_blacklist.txt'%experiment)
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
experiments_with_event_data = ('p4343','p4344','p4345','p4346',
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
fig = plt.figure(figsize=(6,6))
ax1 = plt.subplot(111)


# Panel A
# Set labels and tick sizes
ax1.set_xlabel(r'$\kappa$',fontsize=18)
ax1.set_ylabel(r'Rise Time [s]',fontsize=18)
ax1.tick_params(axis='both', which='major', labelsize=16)

# Turns off chart clutter

# # Turn off top and right tick marks
# ax1.get_xaxis().tick_bottom()
# ax1.get_yaxis().tick_left()
#
# # Turn off top and right splines
# ax1.spines["top"].set_visible(False)
# ax1.spines["right"].set_visible(False)

#ax1.text(-0.2,0.95,'A',transform = ax1.transAxes,fontsize=24)

filter_col = 9
low_val = 40000.
high_val = 50000.
y_col = 11

marker_alpha = 0.3
for key in experiment_event_data:
    event_data = experiment_event_data[key]
    event_data = filter(event_data,filter_col,low_val,high_val)
    y_data = (event_data[:,2]-event_data[:,1])/1000.
    ax1.scatter(event_data[:,5]/0.0007,y_data,color='k',alpha=marker_alpha)
    ax1.errorbar(np.mean(event_data[:,5]/0.0007),np.mean(y_data),fmt='ro',ecolor='w',elinewidth=2,xerr=np.std(event_data[:,5]/0.0007),yerr=np.std(y_data))
    ax1.errorbar(np.mean(event_data[:,5]/0.0007),np.mean(y_data),fmt='ro',markeredgecolor='w',ecolor='k',elinewidth=1,xerr=np.std(event_data[:,5]/0.0007),yerr=np.std(y_data))

# Add audible/non-audible annotation
ax1.axvline(x=0.75,color='k',linestyle='--')
ax1.axvline(x=1.1,color='k',linestyle='--')
#ax1.text(0.647,4.4,'Audible',fontsize=14)
#ax1.text(0.9,4.4,'Silent',fontsize=14)
#ax1.text(1.11,2.555,'Stable',fontsize=14,rotation=90)
ax1.axvspan(1.1,1.15,color='k',alpha=0.2)

ax1.annotate("",
            xy=(0.63, 4.35), xycoords='data',
            xytext=(0.75,4.35), textcoords='data',
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3"),
            )

ax1.annotate("",
            xy=(0.75, 4.35), xycoords='data',
            xytext=(1.1,4.35), textcoords='data',
            arrowprops=dict(arrowstyle="<->",
                            connectionstyle="arc3"),
            )

ax1.annotate("",
            xy=(1.1, 4.35), xycoords='data',
            xytext=(1.15,4.35), textcoords='data',
            arrowprops=dict(arrowstyle="<-",
                            connectionstyle="arc3"),
            )

ax1.annotate("",
            xy=(1.1, 0.35), xycoords='data',
            xytext=(1.15,0.35), textcoords='data',
            arrowprops=dict(arrowstyle="<-",
                            connectionstyle="arc3"),
            )


#ax1.set_ylim(0,4.7)
ax1.set_xlim(0.63,1.15)

plt.savefig('figure.png', bbox_inches="tight");
