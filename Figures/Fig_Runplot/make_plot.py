import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from biaxread import *
from scipy.signal import detrend

def ReadExp(exp,path,disp_low,disp_high):
    data = ReadAscii('%s/%s/%s_data.txt'%(path,exp,exp))
    lower_row = find_nearest(np.ravel(data['LP_Disp']),disp_low)
    upper_row = find_nearest(np.ravel(data['LP_Disp']),disp_high)
    print lower_row,upper_row
    data = data[lower_row:upper_row]
    data['mu'] = data['mu'] - data['mu'][0]
    return data

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
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

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

def rslope(x,y,window):
    """
    Takes a data vector and a window to produce a vector of the running average slope.
    The window specifies the number of points on either side of the central point, so
    the total number of points in the slope fitting is 2*window+1.  Fitting is
    done by the least squares method where the slope is defined by the equation below.
    the beginning and ends are padded with NaN, so fewer points are in those slope
    estimates.  Addition and subtraction to the totals is used so that the sum is not
    recomputed each time, speeding the process.

                    sum(x)*sum(y)
        Sum(x*y) -  -------------
                          n
    m = -------------------------
                     (sum(x))^2
        sum(x^2) - --------------
                          n
    """

    import numpy as np

    # Check that x and y are the same length
    if len(x) != len(y):
        print "Error: x and y must be the same length"
        return 0

    N = len(x) # Number of points in the dataset
    slopes = np.ones(N) # Make array for slopes

    # Pad data with window number of points NaN on either side
    x_padded = np.empty(2*window+N)
    x_padded[0:window] = 0
    x_padded[window:N+window] = x
    x_padded[N+window:2*window+N] = 0

    y_padded = np.empty(2*window+N)
    y_padded[0:window] = 0
    y_padded[window:N+window] = y
    y_padded[N+window:2*window+N] = 0

    sum_x    = np.sum(x_padded[0:2*window+1])
    sum_y    = np.sum(y_padded[0:2*window+1])
    sum_x_sq = np.sum(x_padded[0:2*window+1]*x_padded[0:2*window+1])
    sum_xy   = np.sum(x_padded[0:2*window+1]*y_padded[0:2*window+1])

    n = np.empty(N)
    n[0:window] = np.arange(window+1,2*window+1)
    n[window:N-window] = window*2+1
    n[N-window:N] = np.arange(2*window,window,-1)

    slopes[0] = (sum_xy - (sum_x*sum_y/n[0]))/(sum_x_sq - (sum_x*sum_x/n[0]))

    for i in range(1,N):
        sum_x    = sum_x - x_padded[i-1] + x_padded[2*window+i]
        sum_y    = sum_y - y_padded[i-1] + y_padded[2*window+i]
        sum_x_sq = sum_x_sq - x_padded[i-1]*x_padded[i-1] + \
            x_padded[2*window+i]*x_padded[2*window+i]
        sum_xy   = sum_xy - x_padded[i-1]*y_padded[i-1] +\
            x_padded[2*window+i]*y_padded[2*window+i]
        slopes[i] = (sum_xy - (sum_x*sum_y/n[i]))/(sum_x_sq - (sum_x*sum_x/n[i]))
    return slopes

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

# Read Data
p4347 = ReadAscii(data_path + '/p4347/p4347_data.txt')

#
# Interpolate Data to 1Hz
#
#f = interpolate.interp1d(p4347['Time'],p4347['LP_Disp'])
#p4347_LP_1Hz = interpolate.interp1d(p4347['Time'],p4347['LP_Disp'])


# 4 Panel figure

# A - Runplot of 3 experiments
# B - Zoom multiple runplots

# Setup figure and axes
# Generally plots is ~1.33x width to height (10,7.5 or 12,9)
fig = plt.figure(figsize=(12,9))
axA = plt.subplot(111)

#
# Runplot of p4347
#



# Label experiment
axA.text(0.93,0.02,'p4347',transform = axA.transAxes,fontsize=14)

# Set labels and tick sizes
axA.set_xlabel(r'Load Point Displacement [mm]',fontsize=18)
axA.set_ylabel(r'Friction',fontsize=18)
axA.tick_params(axis='both', which='major', labelsize=16)

# Smooth and plot
smoothed_mu = savitzky_golay(np.ravel(p4347['mu']), 11, 1)
axA.plot(p4347['LP_Disp'][::10]/1000.,smoothed_mu[::10],color='k',linewidth=1)

axA.set_ylim(0.4,0.73)
axA.set_xlim(1,30)

# Add rectangle for where figure B comes from
rect_x1 = 13.
rect_x2 = 13.25
rect_y1 = 0.65
rect_y2 = 0.72
rect_width = rect_x2-rect_x1
rect_height = rect_y2-rect_y1
axA.add_patch(Rectangle((rect_x1,rect_y1),rect_width,rect_height,alpha=0.3, zorder=0,facecolor="k"))



#
# Inset Plot of Stick-Slip Events
#
axB = plt.axes([.45, .47, .4, .2])

axB.set_xlabel(r'Load Point Displacement [mm]',fontsize=16)
axB.set_ylabel(r'Friction',fontsize=16)

# Turns off chart clutter

# Turn off top and right tick marks
#axB.get_xaxis().tick_bottom()
#axB.get_yaxis().tick_left()
axB.get_yaxis().set_ticks([0.67,0.68,0.69,0.7,0.71])
#axB.get_xaxis().set_ticks([])

axB.plot(p4347['LP_Disp'][::10]/1000.,p4347['mu'][::10],color='k',linewidth=2)

axB.set_ylim(0.66,0.71)
axB.set_xlim(13,13.25)

#
# Another Inset
#

inset_data = p4347[1084500:1367671]

# [left, bottom, width, height]
axC = plt.axes([.45, .19, .4, .2])
axC2 = axC.twinx()

axC.set_xlabel(r'Time [sec]',fontsize=16)
axC.set_ylabel(r'Friction',fontsize=16)
axC2.set_ylabel(r'Fault Displacement',fontsize=16)

# Turns off chart clutter

# Turn off top and right tick marks
#axB.get_xaxis().tick_bottom()
#axB.get_yaxis().tick_left()
axC.get_yaxis().set_ticks([])
axC.get_xaxis().set_ticks([])
axC2.get_yaxis().set_ticks([])
axC2.get_xaxis().set_ticks([])

OB_detrend = inset_data['OB_Top'] - inset_data['Time']*((inset_data['OB_Top'][-1] - inset_data['OB_Top'][0])/(inset_data['Time'][-1] - inset_data['Time'][0]))
axC.plot(inset_data['Time'],inset_data['mu'],color='k',linewidth=1)
axC2.plot(inset_data['Time'],OB_detrend,color='b',linewidth=1)

axC.set_ylim(0.672,0.7)
#axC2.set_ylim(-30,5)
print axC2.get_ylim()
axC.set_xlim(3750,3900)
axC2.set_xlim(3750,3900)


plt.savefig('runplot.png', bbox_inches="tight")
