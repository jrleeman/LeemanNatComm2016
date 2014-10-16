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
#f = interpolate.interp1d(p4309['Time'],p4309['LP_Disp'])
#p4309_LP_1Hz = interpolate.interp1d(p4309['Time'],p4309['LP_Disp'])


# 4 Panel figure

# A - Runplot of 3 experiments
# B - Zoom of stable velocity step
# C - Zoom of slow-slip
# D - Zoom of stick-slip

# Setup figure and axes
# Generally plots is ~1.33x width to height (10,7.5 or 12,9)
fig = plt.figure(figsize=(9,9))
axA = fig.add_subplot(2, 1, 1)
axB = fig.add_subplot(2, 3, 4)
axC = fig.add_subplot(2, 3, 5)
axD = fig.add_subplot(2, 3, 6)
plt.subplots_adjust(hspace=0.35)

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

axA.plot(p4311['LP_Disp'][::10]/1000.,p4311['mu'][::10],color=tableau20[2],linewidth=1,
        label='p4311')

axA.plot(p4316['LP_Disp'][::10]/1000.,p4316['mu'][::10],color=tableau20[4],linewidth=1,
        label='p4316',alpha=0.6)

axA.plot(p4309['LP_Disp'][::10]/1000.,p4309['mu'][::10],color=tableau20[0],linewidth=1,
        label='p4309')

axA.set_ylim(0,0.8)
axA.set_xlim(0,25)

#
# Plot B
#

# Cut data down to desired snippet and plot bar on top plot
p4309 = p4309[211642:220000]

# Make second y-axis and turn off everything
axB2 = axB.twinx()
axB2.get_xaxis().set_visible(False)
axB2.get_yaxis().set_visible(False)

# Label Plot
axB.text(0.01,1.0,'B',transform = axB.transAxes,fontsize=24)

# Set labels and tick sizes
axB.set_xlabel(r'',fontsize=18)
axB.set_ylabel(r'Friction',fontsize=18)
axB.xaxis.set_ticklabels([])
axB.yaxis.set_ticklabels([])

# Turns off chart clutter

# Turn off tick marks
axB.get_xaxis().set_ticks([])
axB.get_yaxis().set_ticks([])


# Turn off top and right splines
axB.spines["top"].set_visible(False)
axB.spines["right"].set_visible(False)

# Plot
axB.plot(p4309['Time'] - p4309['Time'][0],p4309['mu'],color=tableau20[0],linewidth=1,
        label='p4309')

axB2.plot(p4309['Time'] - p4309['Time'][0],p4309['On_Board'] - p4309['On_Board'][0],color='k',linewidth=2,
        label='p4309')

# Add scale bars
mu_ref = 0.07 * (axB.get_ylim()[1] - axB.get_ylim()[0]) + axB.get_ylim()[0]
dmu = 0.1 * (axB.get_ylim()[1] - axB.get_ylim()[0])
t_ref = 0.6 * (axB.get_xlim()[1] - axB.get_xlim()[0]) + axB.get_xlim()[0]
dt = 10.
d_ref = 0.07 * (axB2.get_ylim()[1] - axB2.get_ylim()[0]) + axB2.get_ylim()[0]
ddis = 0.1 * (axB2.get_ylim()[1] - axB2.get_ylim()[0])

axB.plot([t_ref,t_ref+dt],[mu_ref,mu_ref],color='k')
axB.plot([t_ref,t_ref],[mu_ref,mu_ref+dmu],color='k')
axB2.plot([t_ref+dt,t_ref+dt],[d_ref,d_ref+ddis],color='k')
axB.text(t_ref,mu_ref-0.5*dmu,'%d s'%dt,fontsize=8)
axB.text(t_ref-2.2*dt,mu_ref+0.5*dmu,'%.4f $\mu$'%dmu,fontsize=8)
axB2.text(t_ref+1.2*dt,d_ref+0.5*ddis,'%.1f $\mu m$'%ddis,fontsize=8)

# Add starting disp text to plot
axB.text(0,-0.07,'%.1f mm' %(p4309['LP_Disp'][0]/1000.),fontsize=10,transform = axB.transAxes)

#
# Plot C
#

# Cut data down to desired snippet and plot bar on top plot
p4311 = p4311[359050:363951]

# Make second y-axis and turn off everything
axC2 = axC.twinx()
axC2.get_xaxis().set_visible(False)
axC2.get_yaxis().set_visible(False)

# Label Plot
axC.text(0.01,1.0,'C',transform = axC.transAxes,fontsize=24)

# Set labels and tick sizes
axC.set_xlabel(r'Time',fontsize=18)
axC.set_ylabel(r'',fontsize=18)
axC.xaxis.set_ticklabels([])
axC.yaxis.set_ticklabels([])

# Turns off chart clutter

# Turn off tick marks
axC.get_xaxis().set_ticks([])
axC.get_yaxis().set_ticks([])


# Turn off top and right splines
axC.spines["top"].set_visible(False)
axC.spines["right"].set_visible(False)
axC.spines["left"].set_visible(False)

# Plot
axC.plot(p4311['Time'] - p4311['Time'][0],p4311['mu'],color=tableau20[2],linewidth=1,
        label='p4311')

axC2.plot(p4311['Time'] - p4311['Time'][0],p4311['On_Board'] - p4311['On_Board'][0],color='k',linewidth=2,
        label='p4311')

# Add scale bars
mu_ref = 0.07 * (axC.get_ylim()[1] - axC.get_ylim()[0]) + axC.get_ylim()[0]
dmu = 0.1 * (axC.get_ylim()[1] - axC.get_ylim()[0])
t_ref = 0.4 * (axC.get_xlim()[1] - axC.get_xlim()[0]) + axC.get_xlim()[0]
dt = 1.
d_ref = 0.07 * (axC2.get_ylim()[1] - axC2.get_ylim()[0]) + axC2.get_ylim()[0]
ddis = 0.1 * (axC2.get_ylim()[1] - axC2.get_ylim()[0])

axC.plot([t_ref,t_ref+dt],[mu_ref,mu_ref],color='k')
axC.plot([t_ref,t_ref],[mu_ref,mu_ref+dmu],color='k')
axC2.plot([t_ref+dt,t_ref+dt],[d_ref,d_ref+ddis],color='k')
axC.text(t_ref+0.3*dt,mu_ref-0.5*dmu,'%d s'%dt,fontsize=8)
axC.text(t_ref-1.3*dt,mu_ref+0.5*dmu,'%.4f $\mu$'%dmu,fontsize=8)
axC2.text(t_ref+1.1*dt,d_ref+0.5*ddis,'%.1f $\mu m$'%ddis,fontsize=8)

# Add starting disp text to plot
axC.text(0,-0.07,'%.1f mm' %(p4311['LP_Disp'][0]/1000.),fontsize=10,transform = axC.transAxes)


#
# Plot D
#

# Cut data down to desired snippet and plot bar on top plot
p4316 = p4316[5120776:5132298]

# Make second y-axis and turn off everything
axD2 = axD.twinx()
axD.get_xaxis().set_visible(False)
axD.get_yaxis().set_visible(False)

# Label Plot
axD.text(0.01,1.0,'D',transform = axD.transAxes,fontsize=24)

# Set labels and tick sizes
axD.set_xlabel(r'',fontsize=18)
axD.set_ylabel(r'',fontsize=18)
axD2.set_ylabel(r'On-Board Displacement',fontsize=18)
axD.xaxis.set_ticklabels([])
axD.yaxis.set_ticklabels([])
axD2.xaxis.set_ticklabels([])
axD2.yaxis.set_ticklabels([])

# Turns off chart clutter

# Turn off tick marks
axD.get_xaxis().set_ticks([])
axD.get_yaxis().set_ticks([])
axD2.get_xaxis().set_ticks([])
axD2.get_yaxis().set_ticks([])


# Turn off top and right splines
axD.spines["top"].set_visible(False)
axD.spines["left"].set_visible(False)

# Plot
axD.plot(p4316['Time'] - p4316['Time'][0],p4316['mu'],color=tableau20[4],linewidth=1,
        label='p4309')

axD2.plot(p4316['Time'] - p4316['Time'][0],p4316['On_Board'] - p4316['On_Board'][0],color='k',linewidth=2,
        label='p4316')

# Add scale bars
mu_ref = 0.07 * (axD.get_ylim()[1] - axD.get_ylim()[0]) + axD.get_ylim()[0]
dmu = 0.1 * (axD.get_ylim()[1] - axD.get_ylim()[0])
t_ref = 0.6 * (axD.get_xlim()[1] - axD.get_xlim()[0]) + axD.get_xlim()[0]
dt = 1.
d_ref = 0.07 * (axD2.get_ylim()[1] - axD2.get_ylim()[0]) + axD2.get_ylim()[0]
ddis = 0.1 * (axD2.get_ylim()[1] - axD2.get_ylim()[0])

axD.plot([t_ref,t_ref+dt],[mu_ref,mu_ref],color='k')
axD.plot([t_ref,t_ref],[mu_ref,mu_ref+dmu],color='k')
axD2.plot([t_ref+dt,t_ref+dt],[d_ref,d_ref+ddis],color='k')
axD.text(t_ref+0.3*dt,mu_ref-0.5*dmu,'%d s'%dt,fontsize=8)
axD.text(t_ref-3.0*dt,mu_ref+0.5*dmu,'%.4f $\mu$'%dmu,fontsize=8)
axD2.text(t_ref+1.1*dt,d_ref+0.5*ddis,'%.1f $\mu m$'%ddis,fontsize=8)

# Add starting disp text to plot
axD.text(0,-0.07,'%.1f mm' %(p4316['LP_Disp'][0]/1000.),fontsize=10,transform = axD.transAxes)


plt.savefig('runplot.png', bbox_inches="tight")
