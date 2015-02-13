import numpy as np
import matplotlib.pyplot as plt

p4342 = np.loadtxt('p4342_event_properties.txt',delimiter=',',skiprows=1)
p4343 = np.loadtxt('p4343_event_properties.txt',delimiter=',',skiprows=1)
p4344 = np.loadtxt('p4344_event_properties.txt',delimiter=',',skiprows=1)
p4345 = np.loadtxt('p4345_event_properties.txt',delimiter=',',skiprows=1)
p4346 = np.loadtxt('p4346_event_properties.txt',delimiter=',',skiprows=1)
p4347 = np.loadtxt('p4347_event_properties.txt',delimiter=',',skiprows=1)
p4348 = np.loadtxt('p4348_event_properties.txt',delimiter=',',skiprows=1)
p4350 = np.loadtxt('p4350_event_properties.txt',delimiter=',',skiprows=1)
p4351 = np.loadtxt('p4351_event_properties.txt',delimiter=',',skiprows=1)

def filter(data,col,low,high):
    inds = np.where(data[:,col]>=low)
    data_trim = data[inds]
    inds = np.where(data_trim[:,col]<=high)
    data_trim = data_trim[inds]
    return data_trim

## Filter based on displacement
low = 30000.
high=32000.
p4342 = filter(p4342,-1,low,high)
p4343 = filter(p4343,-1,low,high)
p4344 = filter(p4344,-1,low,high)
p4345 = filter(p4345,-1,low,high)
p4346 = filter(p4346,-1,low,high)
p4347 = filter(p4347,-1,low,high)
p4348 = filter(p4348,-1,low,high)
p4350 = filter(p4350,-1,low,high)
p4351 = filter(p4351,-1,low,high)

## Filter based on slip time
low = 0.
high=5.
p4342 = filter(p4342,4,low,high)
p4343 = filter(p4343,4,low,high)
p4344 = filter(p4344,4,low,high)
p4345 = filter(p4345,4,low,high)
p4346 = filter(p4346,4,low,high)
p4347 = filter(p4347,4,low,high)
p4348 = filter(p4348,4,low,high)
p4350 = filter(p4350,4,low,high)
p4351 = filter(p4351,4,low,high)

print p4342[:,0]
print p4343[:,0]
print p4344[:,0]

p4346 = p4346[np.in1d(p4346[:,0],[239,240,243,244,245,246,247,248,249,250,251,252,253,254,255,256,257,258,259,260,261,262])]
p4347 = p4347[np.in1d(p4347[:,0],[211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,229,231,232])]
p4348 = p4348[np.in1d(p4348[:,0],[201,202,203,204,205,207,208,211,212,213,216,217,218,221])]

# Setup figure and axes
# Generally plots is ~1.33x width to height (10,7.5 or 12,9)
fig = plt.figure(figsize=(12,9))
ax1 = plt.subplot(111)

# Set labels and tick sizes
ax1.set_ylabel(r'Slip Duration [s]',fontsize=18)
ax1.set_xlabel(r'Stiffness [1/$\mu$m]',fontsize=18)
ax1.tick_params(axis='both', which='major', labelsize=16)

# Turn off top and right tick marks
ax1.get_xaxis().tick_bottom()
ax1.get_yaxis().tick_left()

# Turn off top and right splines
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)

ax1 = plt.subplot(111)
ax1.scatter(p4342[:,5],p4342[:,4],color='k')
ax1.scatter(p4343[:,5],p4343[:,4],color='k')
ax1.scatter(p4344[:,5],p4344[:,4],color='k')
ax1.scatter(p4345[:,5],p4345[:,4],color='k')
ax1.scatter(p4346[:,5],p4346[:,4],color='k')
ax1.scatter(p4347[:,5],p4347[:,4],color='k')
ax1.scatter(p4348[:,5],p4348[:,4],color='k')
ax1.scatter(p4350[:,5],p4350[:,4],color='k')
ax1.scatter(p4351[:,5],p4351[:,4],color='k')

ax1.errorbar(np.mean(p4342[:,5]),np.mean(p4342[:,4]),fmt='ko',ecolor='k',elinewidth=2,xerr=np.std(p4342[:,5]),yerr=np.std(p4342[:,4]))
ax1.errorbar(np.mean(p4343[:,5]),np.mean(p4343[:,4]),fmt='ko',ecolor='k',elinewidth=2,xerr=np.std(p4343[:,5]),yerr=np.std(p4343[:,4]))
ax1.errorbar(np.mean(p4344[:,5]),np.mean(p4344[:,4]),fmt='ko',ecolor='k',elinewidth=2,xerr=np.std(p4344[:,5]),yerr=np.std(p4344[:,4]))
ax1.errorbar(np.mean(p4345[:,5]),np.mean(p4345[:,4]),fmt='ko',ecolor='k',elinewidth=2,xerr=np.std(p4345[:,5]),yerr=np.std(p4345[:,4]))
ax1.errorbar(np.mean(p4346[:,5]),np.mean(p4346[:,4]),fmt='mo',ecolor='k',elinewidth=2,xerr=np.std(p4346[:,5]),yerr=np.std(p4346[:,4]))
ax1.errorbar(np.mean(p4347[:,5]),np.mean(p4347[:,4]),fmt='mo',ecolor='k',elinewidth=2,xerr=np.std(p4347[:,5]),yerr=np.std(p4347[:,4]))
ax1.errorbar(np.mean(p4348[:,5]),np.mean(p4348[:,4]),fmt='mo',ecolor='k',elinewidth=2,xerr=np.std(p4348[:,5]),yerr=np.std(p4348[:,4]))
ax1.errorbar(np.mean(p4350[:,5]),np.mean(p4350[:,4]),fmt='ko',ecolor='k',elinewidth=2,xerr=np.std(p4350[:,5]),yerr=np.std(p4350[:,4]))
ax1.errorbar(np.mean(p4351[:,5]),np.mean(p4351[:,4]),fmt='ko',ecolor='k',elinewidth=2,xerr=np.std(p4351[:,5]),yerr=np.std(p4351[:,4]))

#ax1.set_xlim(3e-4,8e-4)

plt.show()
