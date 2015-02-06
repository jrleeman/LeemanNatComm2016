import numpy as np
import matplotlib.pyplot as plt
from biaxread import *

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

# Setup figure and axes
# Generally plots is ~1.33x width to height (10,7.5 or 12,9)
fig = plt.figure(figsize=(12,9))
axA = plt.subplot(111)

#
# Plot A
#

# Label experiment
axA.text(0.01,0.98,'p4316',transform = axA.transAxes,fontsize=14,color=tableau20[2])

# Set labels and tick sizes
axA.set_xlabel(r'Load Point Displacement [mm]',fontsize=18)
axA.set_ylabel(r'Friction',fontsize=18)
axA.tick_params(axis='both', which='major', labelsize=16)

# Turns off chart clutter

# Turn off top and right tick marks
axA.get_xaxis().tick_bottom()
axA.get_yaxis().tick_left()

# Turn off top and right splines
axA.spines["top"].set_visible(False)
axA.spines["right"].set_visible(False)

axA.plot(p4316['LP_Disp'][::10]/1000.,p4316['mu'][::10],color=tableau20[2],linewidth=1)

axA.set_ylim(0,0.75)
axA.set_xlim(0,42)

plt.savefig('figure.png', bbox_inches="tight")
