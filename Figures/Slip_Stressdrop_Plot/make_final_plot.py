import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from biaxread import *
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

def GetDrops(key,event_data):
    drops = []
    exp_path = '/Users/jleeman/Dropbox/PennState/BiaxExperiments/%s/%s_data.txt' %(key,key)
    exp_data = ReadAscii(exp_path)
    for event in event_data:
        start_row = event[2]
        end_row = event[3]
        friction = exp_data['mu']
        drops.append(friction[start_row]-friction[end_row])
    return drops

def GetSlips(key,event_data):
    drops = []
    exp_path = '/Users/jleeman/Dropbox/PennState/BiaxExperiments/%s/%s_data.txt' %(key,key)
    exp_data = ReadAscii(exp_path)
    for event in event_data:
        start_row = event[2]
        end_row = event[3]
        friction = exp_data['OB_Top']
        drops.append(friction[end_row]-friction[start_row])
    return drops


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
fig = plt.figure(figsize=(12,11))
ax1 = plt.subplot(111)

# Panel A
# Set labels and tick sizes
ax1.set_xlabel(r'Friction Drop',fontsize=24)
ax1.set_ylabel(r'Slip [$\mu m$]',fontsize=18)
ax1.tick_params(axis='both', which='major', labelsize=16)

filter_col = 9
low_val = 40000.
high_val = 50000.
y_col = 11
marker_alpha = 0.3
for key in experiment_event_data:
    event_data = experiment_event_data[key]
    event_data = filter(event_data,filter_col,low_val,high_val)

    stress_drops = GetDrops(key,event_data)
    slips = GetSlips(key,event_data)
    ax1.scatter(stress_drops,slips,color='k',alpha=marker_alpha)
    ax1.scatter(np.mean(stress_drops),np.mean(slips),color='r',s=50,alpha=marker_alpha)

plt.show()
