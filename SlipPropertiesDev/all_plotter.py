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

size=40
alpha=0.5
ax1 = plt.subplot(111)
ax1.scatter(p4342[:,9]/1000.,p4342[:,5],c=p4342[:,4],s=size,alpha=alpha)
ax1.scatter(p4343[:,9]/1000.,p4343[:,5],c=p4343[:,4],s=size,alpha=alpha)
ax1.scatter(p4344[:,9]/1000.,p4344[:,5],c=p4344[:,4],s=size,alpha=alpha)
ax1.scatter(p4345[:,9]/1000.,p4345[:,5],c=p4345[:,4],s=size,alpha=alpha)
ax1.scatter(p4346[:,9]/1000.,p4346[:,5],c=p4346[:,4],s=size,alpha=alpha)
ax1.scatter(p4347[:,9]/1000.,p4347[:,5],c=p4347[:,4],s=size,alpha=alpha)
ax1.scatter(p4348[:,9]/1000.,p4348[:,5],c=p4348[:,4],s=size,alpha=alpha)
ax1.scatter(p4350[:,9]/1000.,p4350[:,5],c=p4350[:,4],s=size,alpha=alpha)
ax1.scatter(p4351[:,9]/1000.,p4351[:,5],c=p4351[:,4],s=size,alpha=alpha)
plt.show()
