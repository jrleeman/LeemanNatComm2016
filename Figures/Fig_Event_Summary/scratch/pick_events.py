import numpy as np
import matplotlib.pyplot as plt

from biaxread import *
from rslope import rslope


def zerocrossings(data):
    """
    Pick where zero crossings occur in an array
    Takes: numpy array of data
    Returns: index of sign difference
    """
    index = np.where(np.diff(np.sign(data)))[0]
    return index

def getvals(array,indexes):
    """
    Get values at given index values from the
    given array and return only those values.
    """
    vals = []
    for i in indexes:
        vals.append(array[i])
    return np.array(vals)

def pick_zero_pairs(x,y,window):
    # Take x,y data and return pairs of index rows
    # that characterize an event

    y_derivative = rslope(x,y,window)

    index_zeros = zerocrossings(y_derivative)

    # Pair data
    if index_zeros.size%2 == 0:
        return index_zeros.reshape(index_zeros.size/2,2)

    else:
        print "Odd number of zero crossings, omitting first value!"
        print "Check output as this may cause an error."
        index_zeros = index_zeros[1:]
        return index_zeros.reshape(index_zeros.size/2,2)

def threshold_events(pairs, minimum, maximum):
    # weed out events that don't meet the criteria
    filtered = []
    for i,pair in enumerate(pairs):
        stress_drop = shear_stress[pair[0]] - shear_stress[pair[1]]
        if (stress_drop < minimum or stress_drop > maximum):
            pass
        else:
            filtered.append(pair)
    return filtered

def min_max_pick(pairs,row_window):
    new_pairs = []
    for i,pair in enumerate(pairs):
        max_row = pair[0]
        min_row = pair[1]
        max_time = data['Time'][max_row]
        min_time = data['Time'][min_row]
        #print "Slip Time ", min_time-max_time
        print max_row,min_row,data['mu'][max_row],data['mu'][min_row]
        max_mu = max(data['mu'][max_row-row_window:max_row+row_window])
        min_mu = min(data['mu'][min_row-row_window:min_row+row_window])
        print max_mu,min_mu
        max_idx = np.where(data['mu'][max_row-row_window:max_row+row_window]==max_mu)[0]
        min_idx = np.where(data['mu'][min_row-row_window:min_row+row_window]==min_mu)[0]
        max_row += (max_idx-row_window)
        min_row += (min_idx-row_window)
        new_pairs.append([int(min_row[0]),int(max_row[0])])
    return new_pairs

path = '/Users/jleeman/Dropbox/PennState/BiaxExperiments/'

# Get experiment name and read data
#experiment = raw_input('Experiment Number: ')
experiment = 'p4344'
data = ReadAscii(path + '%s/%s_data.txt' %(experiment,experiment))


# Get parameters for picking
# window = input('Smoothing Window (points either side of center): ')
# start_row = input('Enter beginning row: ')
# end_row = input('Enter ending row: ')
# min_drop = input('Enter minimum friction drop couted: ')
# max_drop = input('Enter maximum friction drop couted: ')
# outfname = raw_input('Enter Output File Name: ')

window = 50
start_row = 1224502
end_row = 1455413
min_drop = 0.0008
max_drop = 0.1
outfname = 'p4344_1.txt'

data_clipped = data[start_row:end_row+1]

time = data_clipped['Time']
time = time.reshape(time.size)

shear_stress = data_clipped['Shr_stress']
shear_stress = data_clipped['mu']
shear_stress = shear_stress.reshape(shear_stress.size)

pairs = pick_zero_pairs(time,shear_stress,window)

#print "Pre threshold length: ", pairs.shape()

pairs = threshold_events(pairs,min_drop,max_drop)

pairs = min_max_pick(pairs,1)

#print "Post threshold length: ", pairs.shape()

plt.plot(data['Time'],data['mu'],color='k')

for pair in pairs:

    plt.axvline(x=data['Time'][pair[0]+start_row],color='b')
    plt.axvline(x=data['Time'][pair[1]+start_row],color='r')

plt.show()

f = open(path + experiment + "/" + outfname,'w')
f.write('# Starting at row: %d\n'%start_row)
f.write('# Ending at row: %d\n'%end_row)
f.write('# Window: %d\n'%window)
f.write('# Minimum Friction Drop: %f\n'%min_drop)
f.write('# Maximum Friction Drop: %f\n'%max_drop)
f.write('Row_min,Row_max\n')
for pair in pairs:
    f.write('%d,%d\n' %(pair[0]+start_row,pair[1]+start_row))
f.close()
