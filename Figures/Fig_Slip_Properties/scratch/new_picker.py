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
        friction_drop = friction[pair[0]] - friction[pair[1]]
        if (friction_drop < minimum or friction_drop > maximum):
            pass
        else:
            filtered.append(pair)
    return filtered

def min_max_pick(pairs,row_window):
    new_pairs = []
    for i,pair in enumerate(pairs):
        max_row = pair[0]
        min_row = pair[1]
        #max_time = data['Time'][max_row]
        #min_time = data['Time'][min_row]
        #print "Slip Time ", min_time-max_time
        print max_row,min_row,data['mu'][max_row],data['mu'][min_row]

        max_row_search_start = max_row-row_window
        max_row_search_end = max_row+row_window

        min_row_search_start = min_row-row_window
        min_row_search_end = min_row+row_window

        if max_row_search_start < 0:
            max_row_search_start = 0
        if max_row_search_end < 0:
            max_row_search_end = 0

        if min_row_search_start < 0:
            min_row_search_start = 0
        if min_row_search_end < 0:
            min_row_search_end = 0

        max_mu = max(data_clipped['mu'][max_row_search_start:max_row_search_end])
        min_mu = min(data_clipped['mu'][min_row_search_start:min_row_search_end])
        print max_mu,min_mu
        max_idx = np.where(data_clipped['mu'][max_row_search_start:max_row_search_end]==max_mu)[0]
        min_idx = np.where(data_clipped['mu'][min_row_search_start:min_row_search_end]==min_mu)[0]
        print max_idx,min_idx
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
# row_min_max = input('Enter rows for min/max or N: ')
# outfname = raw_input('Enter Output File Name: ')

window = 200

# # For first segment
# start_row = 1214480
# end_row = 1458652

# For second segment
start_row = 1535889
end_row = 2285413

# # For third segment
# start_row = 2543838
# end_row = 3952242

# # For fourth segment
# start_row = 4076047
# end_row = 5488405

min_drop = 0.0007
max_drop = 0.1
outfname = 'p4345_2.txt'

row_min_max = 1000

# Clip on the data we are going to work with, add one to end since
# row numbers are zero indexed
data_clipped = data[start_row:end_row+1]

# Since we'll just use time,friction we make them individual arrays
time = np.ravel(data_clipped['Time'])
friction = np.ravel(data_clipped['mu'])

# Take the smoothed derivative of the time,friction curve
y_derivative = rslope(time,friction,window)

# Get the index of all zero crossings of this derivative
index_zeros = zerocrossings(y_derivative)

# Pair data if possible. First should be max stress, second min stress
if index_zeros.size%2 == 0:
    pairs = index_zeros.reshape(index_zeros.size/2,2)

else:
    print "Odd number of zero crossings, omitting first value!"
    print "Check output as this may cause an error."
    index_zeros = index_zeros[1:]
    pairs =  index_zeros.reshape(index_zeros.size/2,2)

# Threshold the events to get rid of noise
pairs = threshold_events(pairs,min_drop,max_drop)

# Do min/max tweaks if needed
if row_min_max == 'N':
    pass
else:
    row_min_max = int(row_min_max)
pairs = min_max_pick(pairs,row_min_max)

#print "Post threshold length: ", pairs.shape()

plt.plot(data['Time'],data['mu'],color='k')

for pair in pairs:

    plt.axvline(x=data['Time'][pair[0]+start_row],color='b')
    plt.axvline(x=data['Time'][pair[1]+start_row],color='r')

plt.show()

f = open(outfname,'w')
f.write('# Starting at row: %d\n'%start_row)
f.write('# Ending at row: %d\n'%end_row)
f.write('# Window: %d\n'%window)
f.write('# Minimum Friction Drop: %f\n'%min_drop)
f.write('# Maximum Friction Drop: %f\n'%max_drop)
f.write('# Rows for min/max adjust: %s\n'%row_min_max)
f.write('Row_max_friction,Row_min_friction\n')
for pair in pairs:
    f.write('%d,%d\n' %(pair[1]+start_row,pair[0]+start_row))
f.close()
