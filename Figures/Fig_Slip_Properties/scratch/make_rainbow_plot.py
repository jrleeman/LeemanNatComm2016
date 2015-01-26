import numpy as np
import matplotlib.pyplot as plt

def load_event_properties(experiment):
    """
    Load event property file picks for a given experiment
    number and exclude events in the black list (a zero-indexed)
    list of event numbers.
    """

    blacklist = load_blacklist(experiment)
    event_properties = np.loadtxt('%s_event_properties.txt'%experiment,delimiter=',',skiprows=1)
    return np.delete(event_properties,blacklist,axis=0)

def load_blacklist(experiment):
    blacklist = np.loadtxt('%s_blacklist.txt'%experiment)
    return blacklist

p4342 = load_event_properties('p4342')
p4343 = load_event_properties('p4343')
p4344 = load_event_properties('p4344')
p4345 = load_event_properties('p4345')
p4346 = load_event_properties('p4346')
p4347 = load_event_properties('p4347')
p4348 = load_event_properties('p4348')
p4350 = load_event_properties('p4350')
p4351 = load_event_properties('p4351')


size=40
alpha=0.5
max_color = 1000.0
color_col=11
ax1 = plt.subplot(111)
ax1.scatter(p4342[:,9]/1000.,p4342[:,5],c=p4342[:,color_col],s=size,alpha=alpha,vmin=0.0,vmax=max_color,cmap=plt.get_cmap('jet'))
ax1.scatter(p4343[:,9]/1000.,p4343[:,5],c=p4343[:,color_col],s=size,alpha=alpha,vmin=0.0,vmax=max_color,cmap=plt.get_cmap('jet'))
ax1.scatter(p4344[:,9]/1000.,p4344[:,5],c=p4344[:,color_col],s=size,alpha=alpha,vmin=0.0,vmax=max_color,cmap=plt.get_cmap('jet'))
ax1.scatter(p4345[:,9]/1000.,p4345[:,5],c=p4345[:,color_col],s=size,alpha=alpha,vmin=0.0,vmax=max_color,cmap=plt.get_cmap('jet'))
ax1.scatter(p4346[:,9]/1000.,p4346[:,5],c=p4346[:,color_col],s=size,alpha=alpha,vmin=0.0,vmax=max_color,cmap=plt.get_cmap('jet'))
ax1.scatter(p4347[:,9]/1000.,p4347[:,5],c=p4347[:,color_col],s=size,alpha=alpha,vmin=0.0,vmax=max_color,cmap=plt.get_cmap('jet'))
ax1.scatter(p4348[:,9]/1000.,p4348[:,5],c=p4348[:,color_col],s=size,alpha=alpha,vmin=0.0,vmax=max_color,cmap=plt.get_cmap('jet'))
ax1.scatter(p4350[:,9]/1000.,p4350[:,5],c=p4350[:,color_col],s=size,alpha=alpha,vmin=0.0,vmax=max_color,cmap=plt.get_cmap('jet'))
ax1.scatter(p4351[:,9]/1000.,p4351[:,5],c=p4351[:,color_col],s=size,alpha=alpha,vmin=0.0,vmax=max_color,cmap=plt.get_cmap('jet'))

ax1.set_ylabel(r'Stiffness [1/$\mu$m]')
ax1.set_xlabel(r'Load Point Displacement [mm]')
#plt.colorbar()

plt.show()
