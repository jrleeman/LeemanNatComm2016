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

#
# Read Data
#

p4309 = ReadAscii(data_path + '/p4309/p4309_data.txt')
p4311 = ReadAscii(data_path + '/p4311/p4311_data.txt')
p4316 = ReadAscii(data_path + '/p4316/p4316_data.txt')

#
# Interpolate Data to 1Hz
#
f = interpolate.interp1d(p4309['Time'],p4309['LP_Disp'])
p4309_LP_1Hz = interpolate.interp1d(p4309['Time'],p4309['LP_Disp'])


# 4 Panel figure

# A - Runplot of 3 experiments
# B - Zoom of stable velocity step
# C - Zoom of slow-slip
# D - Zoom of stick-slip

# Setup figure and axes
# Generally plots is ~1.33x width to height (10,7.5 or 12,9)
fig = plt.figure(figsize=(12,9))
axA = fig.add_subplot(2, 1, 1)
axB = fig.add_subplot(2, 3, 4)
axC = fig.add_subplot(2, 3, 5)
axD = fig.add_subplot(2, 3, 6)

#
# Plot A
#

# Label Plot
axA.text(0.01,0.9,'A',transform = axA.transAxes,fontsize=24)

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

axA.plot(p4311['LP_Disp'][::100]/1000.,p4311['mu'][::100],color=tableau20[2],linewidth=1,
        label='p4311')

axA.plot(p4316['LP_Disp'][::100]/1000.,p4316['mu'][::100],color=tableau20[4],linewidth=1,
        label='p4316',alpha=0.6)

axA.plot(p4309['LP_Disp'][::100]/1000.,p4309['mu'][::100],color=tableau20[0],linewidth=1,
        label='p4309')

axA.set_ylim(0,0.8)


#
# Plot B
#

# Label Plot
axB.text(0.01,0.9,'B',transform = axB.transAxes,fontsize=24)

# Set labels and tick sizes
axB.set_xlabel(r'',fontsize=18)
axB.set_ylabel(r'Friction',fontsize=18)
axB.tick_params(axis='both', which='major', labelsize=16)

# Turns off chart clutter

# Turn off top and right tick marks
axB.get_xaxis().tick_bottom()
axB.get_yaxis().tick_left()

# Turn off top and right splines
axB.spines["top"].set_visible(False)
axB.spines["right"].set_visible(False)

#axB.plot(p4309['LP_Disp']/1000.,p4309['mu'],color=tableau20[0],linewidth=1,
#        label='p4309')

#
# Plot C
#

# Label Plot
axC.text(0.01,0.9,'C',transform = axC.transAxes,fontsize=24)

# Set labels and tick sizes
axC.set_xlabel(r'Load Point Displacement [mm]',fontsize=18)
axC.set_ylabel(r'',fontsize=18)
axC.tick_params(axis='both', which='major', labelsize=16)

# Turns off chart clutter

# Turn off top and right tick marks
axC.get_xaxis().tick_bottom()
axC.get_yaxis().tick_left()

# Turn off top and right splines
axC.spines["top"].set_visible(False)
axC.spines["right"].set_visible(False)

#axC.plot(p4311['LP_Disp']/1000.,p4311['mu'],color=tableau20[2],linewidth=1,
#        label='p4311')

#
# Plot D
#

# Label Plot
axD.text(0.01,0.9,'D',transform = axD.transAxes,fontsize=24)

# Set labels and tick sizes
axD.set_xlabel(r'',fontsize=18)
axD.set_ylabel(r'',fontsize=18)
axD.tick_params(axis='both', which='major', labelsize=16)

# Turns off chart clutter

# Turn off top and right tick marks
axD.get_xaxis().tick_bottom()
axD.get_yaxis().tick_left()

# Turn off top and right splines
axD.spines["top"].set_visible(False)
axD.spines["right"].set_visible(False)

#axD.plot(p4316['LP_Disp']/1000.,p4316['mu'],color=tableau20[4],linewidth=1,
#        label='p4316')



plt.show()
