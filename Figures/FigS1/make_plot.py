import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from biaxread import *

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


# Path to folders of biax data
data_path = '/Users/jleeman/Dropbox/PennState/BiaxExperiments'

# Make Nice Plot Colors
# These are the "Tableau 20" colors as RGB.
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]

# Scale the RGB values to the [0, 1] range, which is the format
# matplotlib accepts.
for i in range(len(tableau20)):
    r, g, b = tableau20[i]
    tableau20[i] = (r / 255., g / 255., b / 255.)

# Read Data
p4316 = ReadAscii(data_path + '/p4316/p4316_data.txt')
ur_disp,ur_stiffness = np.loadtxt('p4316_stiffness_cycles.txt',unpack=True,delimiter=',',skiprows=1,usecols=[3,4])
#ev_disp,ev_stiffness = np.loadtxt('p4316_event_properties.txt',unpack=True,delimiter=',',skiprows=1,usecols=[9,5])
events = load_events('p4316')
ev_disp = events[:,9]
ev_stiffness = events[:,5]

# Setup figure and axes
# Generally plots is ~1.33x width to height (10,7.5 or 12,9)
fig = plt.figure(figsize=(12,9))
axA = plt.subplot(111)
axB = plt.axes([0.57,0.17,0.3,0.3])
#
# Plot A
#

# Label experiment
axA.text(0.01,0.96,'p4316',transform = axA.transAxes,fontsize=14,color='k')

# Set labels and tick sizes
axA.set_xlabel(r'Load Point Displacement [mm]',fontsize=18)
axA.set_ylabel(r'Friction',fontsize=18)
axA.tick_params(axis='both', which='major', labelsize=16)

# Turns off chart clutter

# Turn off top and right tick marks
#axA.get_xaxis().tick_bottom()
#axA.get_yaxis().tick_left()

# Turn off top and right splines
#axA.spines["top"].set_visible(False)
#axA.spines["right"].set_visible(False)

axA.plot(p4316['LP_Disp'][::10]/1000.,p4316['mu'][::10],color='k',linewidth=1)

# Add horizontal lines to denote where stiffness was picked
axA.axhline(y=0.3,color='k',linestyle='--')
axA.axhline(y=0.4,color='k',linestyle='--')
axA.axhspan(0.3,0.4,color='0.5')

axA.set_ylim(0,0.8)
axA.set_xlim(0,42)

rect = mpatches.Rectangle((22,0.0),19,0.39, ec="none",fc="white",zorder=10)
axA.add_patch(rect)

#
# Inset Axes
#
axB.scatter(ur_disp/1000.,ur_stiffness*10000,color='0.5',s=50)
axB.scatter(ev_disp/1000.,ev_stiffness*10000,color='k',s=15)
axB.set_xlabel(r'Load Point Displacement [mm]',fontsize=14)
axB.set_ylabel(r'Stiffnessx10000 [$\mu m^{-1}$]',fontsize=14)
axB.tick_params(axis='both', which='major', labelsize=12)

axB.set_ylim(0,6)
axB.set_xlim(0,42)

plt.savefig('figure.png', bbox_inches="tight")
