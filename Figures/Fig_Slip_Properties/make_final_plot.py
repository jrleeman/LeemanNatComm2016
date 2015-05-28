import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap
import matplotlib
from biaxread import *

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

def get_kc(disp):
    slope = (7.e-4-2.6e-6)/(16.-6.)
    if disp >= 16:
        return 7e-4
    else:
        return slope*disp - 0.0004

def get_aminusb(experiment,steps):
    data = ReadAscii('/Users/jleeman/Dropbox/PennState/BiaxExperiments/%s/%s_data.txt'%(experiment,experiment))
    picks = np.loadtxt('%s_picks.txt'%experiment,delimiter=',')
    V0,V1,row= np.loadtxt('%s_step_rows.csv'%experiment,delimiter=',',skiprows=1,unpack=True)
    dv = V1-V0
    friction = picks[:,1].reshape((steps,2))
    temp = picks[:,0].reshape((steps,2))
    disp = temp[:,0]/1000

    d_mu = friction[:,1]-friction[:,0]
    amb = d_mu/np.log(V1/V0)

    res = np.array([disp,amb,dv])

    return np.transpose(res)


def bin_steps(steps,bin_width):
    min_disp = np.min(steps[:,0])
    max_disp = np.max(steps[:,0])
    print "min, max", min_disp,max_disp
    print np.shape(steps)

    exclude_dv = [-7]

    for dv in exclude_dv:
        steps = steps[steps[:,2]!=dv]

    disp_means = []
    amb_means = []

    for i in range(int(max_disp/bin_width)+1):
        bin_bottom = i * bin_width
        bin_top = i * bin_width + bin_width
        print "Bin bot,top", bin_bottom, bin_top
        #print steps[:,0] > bin_bottom
        #print steps[:,0] < bin_top
        bin_steps = steps[(steps[:,0] > bin_bottom)]
        bin_steps = bin_steps[(bin_steps[:,0] < bin_top)]
        print "Steps:", np.shape(bin_steps)
        if len(bin_steps)!= 0:
            disp_means.append(np.mean(bin_steps[:,0]))
            amb_means.append(np.mean(bin_steps[:,1]))
            #amb_means.append(np.median(bin_steps[:,1]))
    print bin_steps[:,2]
    return disp_means,amb_means


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

# Tuple of experiments we'll consider for plotting even data from
# Removed p4342 due to data quality issues 2/16/15
experiments_with_event_data = ('p4343','p4344','p4345','p4346',
                               'p4347','p4348','p4350','p4351')

# Tuple of experiments we'll plot unload/reload stiffness from
experiments_with_unload_reload = ('p4267','p4268','p4269','p4270','p4271',
                                  'p4272','p4273','p4309','p4310','p4311',
                                  'p4312','p4313','p4314','p4316','p4317',
                                  'p4327','p4328','p4329','p4330','p4338',
                                  'p4339')

# Read those experiments into a dictionary of event data
experiment_event_data = dict()

for experiment in experiments_with_event_data:
    experiment_event_data[experiment] = load_events(experiment)

# Make the plot
# Setup figure and axes
# Generally plots is ~1.33x width to height (10,7.5 or 12,9)
fig = plt.figure(figsize=(12,13))
#ax1 = plt.subplot(311)
#ax2 = plt.subplot(312)
#ax3 = plt.subplot(313)
ax1 = plt.subplot2grid((5,1), (0,0), rowspan=1)
ax2 = plt.subplot2grid((5,1), (1,0), rowspan=2)
ax3 = plt.subplot2grid((5,1), (3,0), rowspan=2)
ax1.set_position([0.125,0.745,0.775,0.2])
ax3.set_position([0.125,0.1,0.775,0.28])

#
# Plot A top (a-b)
#
p4309_a,p4309_b,p4309_Dc,p4309_amb,step_row = np.loadtxt('p4309_ruina_fits.csv',usecols=[0,2,4,6,9],delimiter=',',skiprows=1,unpack=True)
p4309_data = ReadAscii('/Users/jleeman/Dropbox/PennState/BiaxExperiments/p4309/p4309_data.txt')
step_row = step_row.astype(int)
step_disp = p4309_data['LP_Disp'][step_row]
p4309_step_disp = step_disp/1000.

ax1.set_ylabel(r'(a-b)',fontsize=16)
ax1.tick_params(axis='both', which='major', labelsize=14)

ax1.text(-0.1,0.9,'A',transform = ax1.transAxes,fontsize=24)

ax1.set_xticklabels([])
ax1.get_yaxis().set_ticks([-0.004,-0.002,0.,0.002,0.004])

ax1.scatter(p4309_step_disp,p4309_amb,color='k',
            s=70,marker='.',label='p4309')

ax1.axhline(y=0,color='k',linewidth='2',linestyle='--')

# Label velocity regions
ax1.text(35,0.001,'Velocity Strengthening',fontsize=12)
ax1.text(35,-0.003,'Velocity Weakening',fontsize=12)

ax1.set_xlim(0, 55)
ax1.set_ylim(-0.005 ,0.004)

#ax1.text(48,0.003,'p4309',fontsize=12)

p4381_steps = get_aminusb('p4381',83)
p4382_steps = get_aminusb('p4382',84)

p4381_d,p4381_amb = bin_steps(p4381_steps,5)
p4382_d,p4382_amb = bin_steps(p4382_steps,5)

ax1.scatter(p4381_d,p4381_amb,color='k',marker='v',s=70,label='P4381')
ax1.scatter(p4382_d,p4382_amb,color='k',marker='*',s=70,label='P4382')


handles, labels = ax1.get_legend_handles_labels()
ax1.legend(handles, labels, scatterpoints=1, frameon=False)
#
# Plot A
#

exps = ['p4267','p4268','p4269','p4270','p4271','p4272','p4273',
        'p4309','p4310','p4311','p4312','p4313','p4314','p4316','p4317',
        'p4327','p4328','p4329','p4330']

# Set labels and tick sizes
#ax2.set_xlabel(r'Average LP Displacement [mm]',fontsize=18)
ax2.set_ylabel(r'Stiffness, $k$ [1/$\mu$m]x1000',fontsize=18)
ax2.tick_params(axis='both', which='major', labelsize=16)
ax2.get_yaxis().set_ticks([0,0.5,1,1.5,2,2.5,3,3.5])

ax2.text(-0.1,0.9,'B',transform = ax2.transAxes,fontsize=24)

# Turns off chart clutter

# Turn off top and right tick marks
#ax2.get_xaxis().tick_bottom()
#ax2.get_yaxis().tick_left()

# Turn off top and right splines
#ax2.spines["top"].set_visible(False)
#ax2.spines["right"].set_visible(False)

# Plotting

for exp in experiments_with_unload_reload:

    df = pd.read_csv('/Users/jleeman/Dropbox/PennState/BiaxExperiments/%s/%s_stiffness_cycles.txt'%(exp,exp))

    temp = df[df['Behavior']=='stable']
    ax2.scatter(temp['AvgDisp']/1000.,temp['Slope']*1000,color='k',s=50,alpha=0.6,zorder=50,edgecolor='k')

    #temp = df[df['Behavior']=='slow']
    #ax2.scatter(temp['AvgDisp']/1000.,temp['Slope'],color='r',s=50,alpha=0.6)

    #temp = df[df['Behavior']=='fast']
    #ax2.scatter(temp['AvgDisp']/1000.,temp['Slope'],color='r',s=50,alpha=0.6)

# Add rectangle for where figure B comes from
# rect_x1 = 10.
# rect_x2 = 50.
# rect_y1 = 0.
# rect_y2 = 0.0009*1000
# rect_width = rect_x2-rect_x1
# rect_height = rect_y2-rect_y1
# ax2.add_patch(Rectangle((rect_x1,rect_y1),rect_width,rect_height,alpha=0.2, zorder=0,facecolor="k"))

# Set limits
ax2.set_xlim(0,52)
ax2.set_ylim(0,0.004*1000)

low_color = 10./1000.
high_color = 4600./1000.

marker_size = 40
marker_alpha=0.7
color_col=11

for key in experiment_event_data:
    event_data = experiment_event_data[key]
    sc = ax2.scatter(event_data[:,9]/1000.,event_data[:,5]*1000,s=40,alpha=marker_alpha,color='r',edgecolor='r')
    print key,np.min(event_data[:,color_col]), np.max(event_data[:,color_col])

# Plot line for kc definition
ax2.plot([6,16,52],[2.6e-6*1000,7e-4*1000,7e-4*1000],color='k',linewidth=4)

# Add text
ax2.text(35,0.95,'Stable',fontsize=22)
ax2.text(35,0.15,'Unstable',fontsize=22,color='r')
ax2.text(46,0.88,r'$k_c$',fontsize=26,color='k')



# # Plot Kc
# df = pd.read_excel('/Users/jleeman/Dropbox/PennState/BiaxExperiments/p4309/p4309_rsf_fits.xlsx')
#
#
#
#
# for i,fit in df.iterrows():
#
#     if fit['Grade'] == 'A':
#         #color='#000066'
#         #color='#FFFFFF'
#         color='#0000FF'
#     elif fit['Grade'] == 'B':
#         color='#0066CC'
#         color='#0000FF'
#         #color='#FFFFFF'
#     elif fit['Grade'] == 'C':
#         #color='#00CCFF'
#         color='#FFFFFF'
#         continue
#     elif fit['Grade'] == 'D':
#         #color='#00FFFF'
#         color='#FFFFFF'
#         continue
#
#     if fit['Type']=='Down' and fit['Law']=='r' and fit['k']==0.0055:
#         #ax2.scatter(fit['LP_Disp']/1000.,fit['Kc']*1000,c=color,s=60,marker='v',zorder=50)
#         print fit['LP_Disp']/1000.,fit['Kc']
#
#     elif fit['Type']=='Up' and fit['Law']=='r' and fit['k']==0.0055:
#         #ax2.scatter(fit['LP_Disp']/1000.,fit['Kc']*1000,c=color,s=60,marker='^',zorder=50)
#         print fit['LP_Disp']/1000.,fit['Kc']
#     else:
#         pass



#
# Plot B
#

# Set labels and tick sizes
ax3.set_xlabel(r'Load Point Displacement [mm]',fontsize=18,labelpad=15)
ax3.set_ylabel(r'$\kappa = k/k_c$',fontsize=25)
ax3.tick_params(axis='both', which='major', labelsize=16)

# Turns off chart clutter

# Turn off top and right tick marks
#ax3.get_xaxis().tick_bottom()
#ax3.get_yaxis().tick_left()
ax3.get_yaxis().set_ticks([0,0.2,0.4,0.6,0.8,1.0,1.2])
# Turn off top and right splines
#ax3.spines["top"].set_visible(False)
#ax3.spines["right"].set_visible(False)

# Plotting

# Make panel A of displacement/stiffness
ax3.text(-0.1,0.9,'C',transform = ax3.transAxes,fontsize=24)

low_color = 10./1000.
low_color = 0
high_color = 4000./1000.

cmap = plt.get_cmap('rainbow_r')
start=0.15
stop = 1.
colors = cmap(np.linspace(start, stop, cmap.N))
# Create a new colormap from those colors
color_map = LinearSegmentedColormap.from_list('Upper Half', colors)

marker_size = 40
marker_alpha=0.5
color_col=11

for key in experiment_event_data:
    event_data = experiment_event_data[key]

    k_kc_ratio = []
    for k,disp in zip(event_data[:,5],event_data[:,9]/1000.):
        k_kc_ratio.append(k/get_kc(disp))

    sc = ax3.scatter(event_data[:,9]/1000.,k_kc_ratio,c=event_data[:,color_col]/1000.,s=marker_size,alpha=marker_alpha,vmin=low_color,vmax=high_color,cmap=color_map)
    print key,np.min(event_data[:,color_col]), np.max(event_data[:,color_col])
# cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
# plt.colorbar(sc,cax=cbar_ax)
# for experiment in experiments_with_unload_reload:
#     df = pd.read_csv('/Users/jleeman/Dropbox/PennState/BiaxExperiments/%s/%s_stiffness_cycles.txt'%(experiment,experiment))
#
#     ax3.scatter(df['AvgDisp']/1000.,df['Slope'],color='g',s=50,alpha=0.6)

position=fig.add_axes([0.37,0.16,0.5,0.02])  ## the parameters are the specified position you set [left, bottom, width, height]
cb = fig.colorbar(sc,cax=position,orientation='horizontal', drawedges=False)
cb.solids.set_edgecolor("face")
cb.set_label(r'Peak Slip Velocity [$mm/s$]',fontsize=14)
cb.set_alpha(1)
cb.draw_all()
#position.set_xlim(0,4)

ax3.set_ylim(0,1.4)
ax3.set_xlim(16,50)

ax3.axvspan(40, 50, alpha=0.2, color='k', zorder=0)

# Add call out lines between plots

transFigure = fig.transFigure.inverted()
### LEFT
coord1 = transFigure.transform(ax2.transData.transform([16,0]))
coord2 = transFigure.transform(ax2.transData.transform([16,-0.3]))
line1 = matplotlib.lines.Line2D((coord1[0],coord2[0]),(coord1[1],coord2[1]),
                               transform=fig.transFigure,color='k')

coord3 = transFigure.transform(ax3.transData.transform([16,1.4]))
line2 = matplotlib.lines.Line2D((coord2[0],coord3[0]),(coord2[1],coord3[1]),
                               transform=fig.transFigure,color='k')
### RIGHT
coord1 = transFigure.transform(ax2.transData.transform([50,0]))
coord2 = transFigure.transform(ax2.transData.transform([50,-0.3]))
line3 = matplotlib.lines.Line2D((coord1[0],coord2[0]),(coord1[1],coord2[1]),
                               transform=fig.transFigure,color='k')

coord3 = transFigure.transform(ax3.transData.transform([50,1.4]))
line4 = matplotlib.lines.Line2D((coord2[0],coord3[0]),(coord2[1],coord3[1]),
                               transform=fig.transFigure,color='k')

fig.lines = line1,line2,line3,line4

# coord1 = transFigure.transform(ax2.transData.transform([16,0]))
# coord2 = transFigure.transform(ax3.transData.transform([16,1.4]))
#
# line = matplotlib.lines.Line2D((coord1[0],coord2[0]),(coord1[1],coord2[1]),
#                                transform=fig.transFigure,color='k')
# fig.lines = line,

plt.savefig('figure.png', bbox_inches="tight");
