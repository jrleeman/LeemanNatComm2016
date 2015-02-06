import pylab as plt
import sys
import numpy as np
from biaxread import *


event_rows_adjusted = []

#pairs = rows.reshape(rows.size/2,2)

experiment = sys.argv[1]

# Load the biax data from the given project folder
path = '/Users/jleeman/Dropbox/PennState/BiaxExperiments/'

rows = np.loadtxt('%s_manual_picks.txt'%experiment,usecols=[0],unpack=True,skiprows=1)

# Get experiment name and read data
#experiment = raw_input('Experiment Number: ')
data = ReadAscii(path + '%s/%s_data.txt' %(experiment,experiment))


for i in range(len(rows)-1):
    pair = [0,0]
    pair[0] = int(rows[i])
    pair[1] = int(rows[i+1])
    event_friction = np.ravel(data['mu'])[pair[0]:pair[1]]
    event_rows = np.ravel(data['row_num'])[pair[0]:pair[1]]
    #print pair
    max_friction = np.max(event_friction)
    min_friction = np.min(event_friction)

    max_index = np.where(event_friction==max_friction)[0][0]
    min_index = np.where(event_friction==min_friction)[0][0]

    print max_index,min_index

    max_index = int(max_index)
    min_index = int(min_index)


    row_max_friction = event_rows[max_index]
    row_min_friction = event_rows[min_index]

    #print row_max_friction,row_min_friction
    event_rows_adjusted.append([row_max_friction,row_min_friction])

plt.plot(data['row_num'],data['mu'],color='k')

outfile = open('%s_event_rows.txt'%experiment,'w')
outfile.write('MaxFrictionRow,MinFrictionRow\n')

for event in event_rows_adjusted:
    plt.axvline(x=event[0],color='r')
    plt.axvline(x=event[1],color='b')
    outfile.write('%d,%d\n'%(event[0],event[1]))

outfile.close()
plt.show()
