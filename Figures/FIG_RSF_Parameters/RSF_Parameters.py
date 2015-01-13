import matplotlib.pylab as plt
import numpy as np
import pandas as pd

# These are the "Tableau 20" colors as RGB.
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]

# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.
for i in range(len(tableau20)):
    r, g, b = tableau20[i]
    tableau20[i] = (r / 255., g / 255., b / 255.)

df = pd.read_excel('p4309_rsf_fits.xlsx')

data = df[df['Law']=='r']
data = data[data['k']==0.0055]
data =  data.query('Grade == ["A","B"]')

# Setup figure and axes
# Generally plots is ~1.33x width to height (10,7.5 or 12,9)
fig = plt.figure(figsize=(17,8))
ax3 = fig.add_subplot(1, 2, 1)
ax4 = fig.add_subplot(1, 2, 2)
#fig.subplots_adjust(hspace=0.1, wspace=0.35)


#
# Left Plot of a-b
#

# Label Plot
ax3.text(0.01,0.9,'A',transform = ax3.transAxes,fontsize=24)

# Set labels and tick sizes
ax3.set_xlabel(r'Displacement [mm]',fontsize=18)
ax3.set_ylabel(r'(a-b)',fontsize=18)
ax3.tick_params(axis='both', which='major', labelsize=16)

# Turns off chart clutter

# Turn off top and right tick marks
ax3.get_xaxis().tick_bottom()
ax3.get_yaxis().tick_left()

# Turn off top and right splines
ax3.spines["top"].set_visible(False)
ax3.spines["right"].set_visible(False)

# Plotting
up = data[data['Type']=='Up']
ax3.scatter(up['LP_Disp']/1000,(up['a']-up['b']),color=tableau20[6],
            s=50,marker='^', label='Velocity Step Up')

down = data[data['Type']=='Down']
ax3.scatter(down['LP_Disp']/1000,(down['a']-down['b']),color=tableau20[7],
            s=50,marker='v', label='Velocity Step Down')

ax3.axhline(y=0,color='k',linewidth='2',linestyle='--')

# Label velocity regions
ax3.text(15,0.0005,'Velocity Strengthening',fontsize=10)
ax3.text(16,-0.0007,'Velocity Weakening',fontsize=10)

ax3.set_ylim(-0.006 ,0.006)

#
# Right Plot of Dc
#

# Label Plot
ax4.text(0.01,0.9,'B',transform = ax4.transAxes,fontsize=24)

# Set labels and tick sizes
ax4.set_xlabel(r'Displacement [mm]',fontsize=18)
ax4.set_ylabel(r'Dc [$\mu m$]',fontsize=18)
ax4.tick_params(axis='both', which='major', labelsize=16)

# Turns off chart clutter

# Turn off top and right tick marks
ax4.get_xaxis().tick_bottom()
ax4.get_yaxis().tick_left()

# Turn off top and right splines
ax4.spines["top"].set_visible(False)
ax4.spines["right"].set_visible(False)

# Plotting
up = data[data['Type']=='Up']
ax4.scatter(up['LP_Disp']/1000,up['Dc'],color=tableau20[4],
            s=50,marker='^', label='Velocity Step Up')

down = data[data['Type']=='Down']
ax4.scatter(down['LP_Disp']/1000,down['Dc'],color=tableau20[5],
            s=50,marker='v', label='Velocity Step Down')

ax4.set_ylim(0.0,35)

plt.savefig('RSF_Parameters.png', bbox_inches="tight")
