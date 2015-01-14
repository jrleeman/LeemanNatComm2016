import numpy as np
import matplotlib.pyplot as plt
from biaxread import *

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
p4309 = ReadAscii(data_path + '/p4309/p4309_data.txt')
p4311 = ReadAscii(data_path + '/p4311/p4311_data.txt')
p4316 = ReadAscii(data_path + '/p4316/p4316_data.txt')

path = '/Users/jleeman/Dropbox/PennState/BiaxExperiments'
dis_low = 30*1000
dis_high = 31*1000
p4343 = ReadExp('p4343',path,dis_low,dis_high)
p4345 = ReadExp('p4345',path,dis_low,dis_high)
p4347 = ReadExp('p4347',path,dis_low,dis_high)
p4342 = ReadExp('p4342',path,dis_low,dis_high)
p4351 = ReadExp('p4351',path,dis_low,dis_high)


#
# Interpolate Data to 1Hz
#
#f = interpolate.interp1d(p4309['Time'],p4309['LP_Disp'])
#p4309_LP_1Hz = interpolate.interp1d(p4309['Time'],p4309['LP_Disp'])


# 4 Panel figure

# A - Runplot of 3 experiments
# B - Zoom multiple runplots

# Setup figure and axes
# Generally plots is ~1.33x width to height (10,7.5 or 12,9)
fig = plt.figure(figsize=(9,13))
axA = fig.add_subplot(2, 1, 1)
axB = fig.add_subplot(2, 1, 2)
plt.subplots_adjust(hspace=0.35)

#
# Plot A
#



# Label Plot
axA.text(-0.13,0.95,'A',transform = axA.transAxes,fontsize=32)

# Set labels and tick sizes
axA.set_xlabel(r'Load Point Displacement [mm]',fontsize=18)
axA.set_ylabel(r'Friction',fontsize=18)
axA.tick_params(axis='both', which='major', labelsize=16)

# Turns off chart clutter

# Turn off top and right tick marks
axA.get_xaxis().tick_bottom()
axA.get_yaxis().tick_left()

# Turn off top and right splines
axA.spines["top"].set_visible(False)
axA.spines["right"].set_visible(False)

axA.plot(p4311['LP_Disp'][::10]/1000.,p4311['mu'][::10],color=tableau20[2],linewidth=1,
        label='p4311')

axA.plot(p4316['LP_Disp'][::10]/1000.,p4316['mu'][::10],color=tableau20[4],linewidth=1,
        label='p4316',alpha=0.6)

axA.plot(p4309['LP_Disp'][::10]/1000.,p4309['mu'][::10],color=tableau20[0],linewidth=1,
        label='p4309')

axA.set_ylim(0,0.8)
axA.set_xlim(0,25)

#
# Plot B
#



# Set labels and tick sizes
axB.set_xlabel(r'Load Point Displacement [mm]',fontsize=18)
axB.set_ylabel(r'Friction',fontsize=18)
axB.tick_params(axis='both', which='major', labelsize=16)

# Label Plot
axB.text(-0.13,0.95,'B',transform = axB.transAxes,fontsize=32)

# Turns off chart clutter

# Turn off top and right tick marks
axB.get_xaxis().tick_bottom()
axB.get_yaxis().tick_left()
axB.get_yaxis().set_ticks([])

# Turn off top and right splines
axB.spines["top"].set_visible(False)
axB.spines["right"].set_visible(False)

# Mask unload in p4338
# indices_to_mask = p4338['mu'] < -0.03
# p4338['mu'][indices_to_mask] = np.nan

# Plotting
window_size = 5
order = 3
axB.plot(p4343['LP_Disp']/1000.,savitzky_golay(np.ravel(p4343['mu']), window_size, order)+0.05*0.5,label='6 MPa',color=tableau20[0])
axB.plot(p4345['LP_Disp']/1000.,savitzky_golay(np.ravel(p4345['mu']), window_size, order)+0.05*2.5,label='8 MPa',color=tableau20[4])
axB.plot(p4347['LP_Disp']/1000.,savitzky_golay(np.ravel(p4347['mu']), window_size, order)+0.05*4.5,label='10 MPa',color=tableau20[8])
axB.plot(p4342['LP_Disp']/1000.,savitzky_golay(np.ravel(p4342['mu']), window_size, order)+0.05*6.5,label='12 MPa',color=tableau20[12])
axB.plot(p4351['LP_Disp']/1000.,savitzky_golay(np.ravel(p4351['mu']), window_size, order)+0.05*8.5,label='14 MPa',color=tableau20[18])

axB.text(30.4,np.max(p4343['mu'])+0.05*0.6,r'$\sigma_n$ = 6 MPa',fontsize=12,color=tableau20[0])
axB.text(30.4,np.max(p4345['mu'])+0.05*2.6,r'$\sigma_n$ = 8 MPa',fontsize=12,color=tableau20[4])
axB.text(30.4,np.max(p4347['mu'])+0.05*4.6,r'$\sigma_n$ = 10 MPa',fontsize=12,color=tableau20[8])
axB.text(30.4,np.max(p4342['mu'])+0.05*6.6,r'$\sigma_n$ = 12 MPa',fontsize=12,color=tableau20[12])
axB.text(30.4,np.max(p4351['mu'])+0.05*8.6,r'$\sigma_n$ = 14 MPa',fontsize=12,color=tableau20[18])

axB.text(30.4,np.min(p4343['mu'])+0.05*0.3,r'p4343',fontsize=10,color=tableau20[0])
axB.text(30.4,np.min(p4345['mu'])+0.05*2.4,r'p4345',fontsize=10,color=tableau20[4])
axB.text(30.4,np.min(p4347['mu'])+0.05*4.4,r'p4347',fontsize=10,color=tableau20[8])
axB.text(30.4,np.min(p4342['mu'])+0.05*6.3,r'p4342',fontsize=10,color=tableau20[12])
axB.text(30.4,np.min(p4351['mu'])+0.05*8.2,r'p4351',fontsize=10,color=tableau20[18])

# Set limits
axB.set_xlim(30,30.5)
axB.set_ylim(0,0.45)
print axB.get_ylim()

plt.savefig('runplot.png', bbox_inches="tight")
