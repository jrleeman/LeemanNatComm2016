import pylab as plt
import sys
import numpy as np
from biaxread import *

try:
    experiment         = sys.argv[1]
    print "File: %s" %experiment


except:
    print "Incorrect Usage!"
    sys.exit()


def onclick(event):

    if event.button == 3:
        print 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(
        event.button, event.x, event.y, event.xdata, event.ydata)

        temp.write('\n%.4f %.4f' %(event.xdata,event.ydata))

        x.append(event.xdata); y.append(event.ydata)
        ax1.scatter(x,y, color='r', s=50)
        plt.draw()

def onkey(event):
    if event.key =='z':
        print "Zooming plot"
        ax = plt.gca()
        lims = ax.get_xlim()
        ZoomAx(lims)

def PadLimits(limits,pad):
    min = limits[0] - limits[0]*pad
    max = limits[1] + limits[1]*pad
    return min,max

def ZoomAx(x_lims):
    x_min = FindNearest(biax['row_num'],x_lims[0])
    x_max = FindNearest(biax['row_num'],x_lims[1])
    ss_min = min(biax['Shr_stress'][x_min:x_max])[0]
    ss_max = max(biax['Shr_stress'][x_min:x_max])[0]

    ss_min,ss_max = PadLimits([ss_min,ss_max],0.02)

    ax1.set_ylim(ss_min,ss_max)
    plt.draw()

    return 1

def FindNearest(array,value):
    index=(np.abs(array-value)).argmin()
    return index


def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

temp = open('temp.txt','w')
temp.write('RowNumber Stress')

x=[];y=[]
# Load the biax data from the given project folder
path = '/Users/jleeman/Dropbox/PennState/BiaxExperiments/'

# Get experiment name and read data
#experiment = raw_input('Experiment Number: ')
biax = ReadAscii(path + '%s/%s_data.txt' %(experiment,experiment))

fig = plt.figure()

ax1 = plt.subplot(111)

cid = fig.canvas.mpl_connect('button_press_event', onclick)
cid = fig.canvas.mpl_connect('key_press_event', onkey)

ax1.plot(biax['row_num'],biax['Shr_stress'],color='k')
plt.show()

temp.close()
