import numpy as np
import matplotlib.pyplot as plt
from biaxread import *

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

#
# Read Data
#

p4309 = ReadAscii(data_path + '/p4309/p4309_data.txt')
p4311 = ReadAscii(data_path + '/p4311/p4311_data.txt')
p4316 = ReadAscii(data_path + '/p4316/p4316_data.txt')
window = 20 # Window for velocity calculation

#
# Interpolate Data to 1Hz
#
#f = interpolate.interp1d(p4309['Time'],p4309['LP_Disp'])
#p4309_LP_1Hz = interpolate.interp1d(p4309['Time'],p4309['LP_Disp'])


# 4 Panel figure

# A - Runplot of 3 experiments
# B - Zoom of stable velocity step
# C - Zoom of slow-slip
# D - Zoom of stick-slip

# Setup figure and axes
# Generally plots is ~1.33x width to height (10,7.5 or 12,9)
fig = plt.figure(figsize=(9,9))
axB = fig.add_subplot(3, 3, 1)
axC = fig.add_subplot(3, 3, 2)
axD = fig.add_subplot(3, 3, 3)
axE = fig.add_subplot(3, 3, 4)
axF = fig.add_subplot(3, 3, 5)
axG = fig.add_subplot(3, 3, 6)
plt.subplots_adjust(hspace=0.35)



#
# Trim Data for zoomed in plots
#

# Cut data down to desired snippet and plot bar on top plot
p4309 = p4309[211642:220000]
p4311 = p4311[359050:363951]
p4316 = p4316[5120776:5132298]


#
# Plot B
#

# Label Plot
axB.text(0.01,1.0,'A',transform = axB.transAxes,fontsize=24)

# Set labels and tick sizes
axB.set_xlabel(r'',fontsize=18)
axB.set_ylabel(r'Friction',fontsize=18)
axB.xaxis.set_ticklabels([])
axB.yaxis.set_ticklabels([])

# Turns off chart clutter

# Turn off tick marks
axB.get_xaxis().set_ticks([])
axB.get_yaxis().set_ticks([])


# Turn off top and right splines
axB.spines["top"].set_visible(False)
axB.spines["right"].set_visible(False)

# Plot
#axB.plot(p4309['Time'] - p4309['Time'][0],p4309['mu'],color='k',linewidth=1,
#        label='p4309')

axB.plot(p4309['Time'] - p4309['Time'][0],savitzky_golay(np.ravel(p4309['mu']),201,5),color=tableau20[0],linewidth=1,
        label='p4309')

# Add scale bars
mu_ref = 0.07 * (axB.get_ylim()[1] - axB.get_ylim()[0]) + axB.get_ylim()[0]
dmu = 0.1 * (axB.get_ylim()[1] - axB.get_ylim()[0])
t_ref = 0.6 * (axB.get_xlim()[1] - axB.get_xlim()[0]) + axB.get_xlim()[0]
dt = 10.

axB.plot([t_ref,t_ref+dt],[mu_ref,mu_ref],color='k')
axB.plot([t_ref,t_ref],[mu_ref,mu_ref+dmu],color='k')
axB.text(t_ref,mu_ref-0.5*dmu,'%d s'%dt,fontsize=8)
axB.text(t_ref-2.2*dt,mu_ref+0.5*dmu,'%.4f $\mu$'%dmu,fontsize=8)

# Add starting disp text to plot
axB.text(0,-0.07,'%.1f mm' %(p4309['LP_Disp'][0]/1000.),fontsize=10,transform = axB.transAxes)

#
# Plot C
#


# Label Plot
axC.text(0.01,1.0,'B',transform = axC.transAxes,fontsize=24)

# Set labels and tick sizes
axC.set_xlabel(r'Time',fontsize=18)
axC.set_ylabel(r'',fontsize=18)
axC.xaxis.set_ticklabels([])
axC.yaxis.set_ticklabels([])

# Turns off chart clutter

# Turn off tick marks
axC.get_xaxis().set_ticks([])
axC.get_yaxis().set_ticks([])


# Turn off top and right splines
axC.spines["top"].set_visible(False)
axC.spines["right"].set_visible(False)
axC.spines["left"].set_visible(False)

# Plot
#axC.plot(p4311['Time'] - p4311['Time'][0],p4311['mu'],color='k',linewidth=1,
#        label='p4311')

axC.plot(p4311['Time'] - p4311['Time'][0],savitzky_golay(np.ravel(p4311['mu']),201,5),color=tableau20[2],linewidth=1,
        label='p4311')


# Add scale bars
mu_ref = 0.07 * (axC.get_ylim()[1] - axC.get_ylim()[0]) + axC.get_ylim()[0]
dmu = 0.1 * (axC.get_ylim()[1] - axC.get_ylim()[0])
t_ref = 0.4 * (axC.get_xlim()[1] - axC.get_xlim()[0]) + axC.get_xlim()[0]
dt = 1.

axC.plot([t_ref,t_ref+dt],[mu_ref,mu_ref],color='k')
axC.plot([t_ref,t_ref],[mu_ref,mu_ref+dmu],color='k')
axC.text(t_ref+0.3*dt,mu_ref-0.5*dmu,'%d s'%dt,fontsize=8)
axC.text(t_ref-1.3*dt,mu_ref+0.5*dmu,'%.4f $\mu$'%dmu,fontsize=8)

# Add starting disp text to plot
axC.text(0,-0.07,'%.1f mm' %(p4311['LP_Disp'][0]/1000.),fontsize=10,transform = axC.transAxes)


#
# Plot D
#

# Label Plot
axD.text(0.01,1.0,'C',transform = axD.transAxes,fontsize=24)

# Set labels and tick sizes
axD.set_xlabel(r'',fontsize=18)
axD.set_ylabel(r'',fontsize=18)
axD.xaxis.set_ticklabels([])
axD.yaxis.set_ticklabels([])


# Turns off chart clutter

# Turn off tick marks
axD.get_xaxis().set_ticks([])
axD.get_yaxis().set_ticks([])


# Turn off top and right splines
axD.spines["top"].set_visible(False)
axD.spines["left"].set_visible(False)

# Plot
axD.plot(p4316['Time'] - p4316['Time'][0],p4316['mu'],color=tableau20[4],linewidth=1,
        label='p4309')

# Add scale bars
mu_ref = 0.07 * (axD.get_ylim()[1] - axD.get_ylim()[0]) + axD.get_ylim()[0]
dmu = 0.1 * (axD.get_ylim()[1] - axD.get_ylim()[0])
t_ref = 0.6 * (axD.get_xlim()[1] - axD.get_xlim()[0]) + axD.get_xlim()[0]
dt = 1.

axD.plot([t_ref,t_ref+dt],[mu_ref,mu_ref],color='k')
axD.plot([t_ref,t_ref],[mu_ref,mu_ref+dmu],color='k')
axD.text(t_ref+0.3*dt,mu_ref-0.5*dmu,'%d s'%dt,fontsize=8)
axD.text(t_ref-3.0*dt,mu_ref+0.5*dmu,'%.4f $\mu$'%dmu,fontsize=8)

# Add starting disp text to plot
axD.text(0,-0.07,'%.1f mm' %(p4316['LP_Disp'][0]/1000.),fontsize=10,transform = axD.transAxes)

#
# Plot E
#

# Make second y-axis and turn off everything
axE2 = axE.twinx()
axE2.get_xaxis().set_visible(False)
axE2.get_yaxis().set_visible(False)

# Label Plot
#axE.text(0.01,1.0,'C',transform = axE.transAxes,fontsize=24)

# Set labels and tick sizes
axE.set_xlabel(r'',fontsize=18)
axE.set_ylabel(r'Displacement',fontsize=18)
axE2.set_ylabel(r'',fontsize=18)
axE.xaxis.set_ticklabels([])
axE.yaxis.set_ticklabels([])
axE2.xaxis.set_ticklabels([])
axE2.yaxis.set_ticklabels([])

# Turns off chart clutter

# Turn off tick marks
axE.get_xaxis().set_ticks([])
axE.get_yaxis().set_ticks([])
axE2.get_xaxis().set_ticks([])
axE2.get_yaxis().set_ticks([])


# Turn off top and right splines
axE.spines["top"].set_visible(False)
axE.spines["right"].set_visible(False)
axE2.spines["top"].set_visible(False)
axE2.spines["right"].set_visible(False)

# Plot
axE.plot(p4309['Time'] - p4309['Time'][0],p4309['On_Board'] - p4309['On_Board'][0],color='k',linewidth=2,
        label='p4311')

velocity = rslope(p4309['Time'].reshape(p4309['Time'].size),p4309['On_Board'].reshape(p4309['On_Board'].size),window)
axE2.plot(p4309['Time'] - p4309['Time'][0],velocity,color='0.7',linewidth=2,
        linestyle=':',label='p4309')

# Set Velocity Axis Limits
axE2.set_ylim(0,15)

# Add annotation for velocity
axE2.text(10,11,r'10 $\mu m/s$',fontsize=10)
axE2.text(25,4,r'1 $\mu m/s$',fontsize=10)


# Add scale bars
#mu_ref = 0.07 * (axD.get_ylim()[1] - axD.get_ylim()[0]) + axD.get_ylim()[0]
#dmu = 0.1 * (axD.get_ylim()[1] - axD.get_ylim()[0])
#t_ref = 0.6 * (axD.get_xlim()[1] - axD.get_xlim()[0]) + axD.get_xlim()[0]
#dt = 1.
#d_ref = 0.07 * (axD2.get_ylim()[1] - axD2.get_ylim()[0]) + axD2.get_ylim()[0]
#ddis = 0.1 * (axD2.get_ylim()[1] - axD2.get_ylim()[0])

#axE.plot([t_ref,t_ref+dt],[mu_ref,mu_ref],color='k')
#axE.plot([t_ref,t_ref],[mu_ref,mu_ref+dmu],color='k')
#axE2.plot([t_ref+dt,t_ref+dt],[d_ref,d_ref+ddis],color='k')
#axE.text(t_ref+0.3*dt,mu_ref-0.5*dmu,'%d s'%dt,fontsize=8)
#axE.text(t_ref-3.0*dt,mu_ref+0.5*dmu,'%.4f $\mu$'%dmu,fontsize=8)
#axE2.text(t_ref+1.1*dt,d_ref+0.5*ddis,'%.1f $\mu m$'%ddis,fontsize=8)

# Add starting disp text to plot
#axE.text(0,-0.07,'%.1f mm' %(p4316['LP_Disp'][0]/1000.),fontsize=10,transform = axD.transAxes)


#
# Plot F
#

# Make second y-axis and turn off everything
axF2 = axF.twinx()
axF2.get_xaxis().set_visible(False)
axF2.get_yaxis().set_visible(False)

# Label Plot
#axF.text(0.01,1.0,'C',transform = axF.transAxes,fontsize=24)

# Set labels and tick sizes
axF.set_xlabel(r'Time',fontsize=18)
axF.set_ylabel(r'',fontsize=18)
axF.xaxis.set_ticklabels([])
axF.yaxis.set_ticklabels([])

# Turns off chart clutter

# Turn off tick marks
axF.get_xaxis().set_ticks([])
axF.get_yaxis().set_ticks([])


# Turn off top and right splines
axF.spines["top"].set_visible(False)
axF.spines["right"].set_visible(False)
axF.spines["left"].set_visible(False)
axF2.spines["top"].set_visible(False)
axF2.spines["right"].set_visible(False)
axF2.spines["left"].set_visible(False)

# Plot
#axF.plot(p4311['Time'] - p4311['Time'][0],p4311['mu'],color=tableau20[2],linewidth=1,
#        label='p4311')

axF.plot(p4311['Time'] - p4311['Time'][0],p4311['On_Board'] - p4311['On_Board'][0],color='k',linewidth=2,
        label='p4311')

velocity = rslope(p4311['Time'].reshape(p4311['Time'].size),p4311['On_Board'].reshape(p4311['On_Board'].size),window)
axF2.plot(p4311['Time'] - p4311['Time'][0],velocity,color='0.7',linewidth=2,
        linestyle=':',label='p4311')

# Set Velocity Axis Limits
axF2.set_ylim(0,50)

# Add annotation for velocity
axF2.text(1.5,1.5,r'10 $\mu m/s$',fontsize=10)
axF2.text(2.6,45,r'45 $\mu m/s$',fontsize=10)

# # Add scale bars
# mu_ref = 0.07 * (axF.get_ylim()[1] - axF.get_ylim()[0]) + axF.get_ylim()[0]
# dmu = 0.1 * (axF.get_ylim()[1] - axF.get_ylim()[0])
# t_ref = 0.4 * (axF.get_xlim()[1] - axF.get_xlim()[0]) + axF.get_xlim()[0]
# dt = 1.
# d_ref = 0.07 * (axF2.get_ylim()[1] - axF2.get_ylim()[0]) + axF2.get_ylim()[0]
# ddis = 0.1 * (axF2.get_ylim()[1] - axF2.get_ylim()[0])
#
# axF.plot([t_ref,t_ref+dt],[mu_ref,mu_ref],color='k')
# axF.plot([t_ref,t_ref],[mu_ref,mu_ref+dmu],color='k')
# axF2.plot([t_ref+dt,t_ref+dt],[d_ref,d_ref+ddis],color='k')
# axF.text(t_ref+0.3*dt,mu_ref-0.5*dmu,'%d s'%dt,fontsize=8)
# axF.text(t_ref-1.3*dt,mu_ref+0.5*dmu,'%.4f $\mu$'%dmu,fontsize=8)
# axF2.text(t_ref+1.1*dt,d_ref+0.5*ddis,'%.1f $\mu m$'%ddis,fontsize=8)

# Add starting disp text to plot
#axF.text(0,-0.07,'%.1f mm' %(p4311['LP_Disp'][0]/1000.),fontsize=10,transform = axF.transAxes)

#
# Plot G
#

# Make second y-axis and turn off everything
axG2 = axG.twinx()
axG.get_xaxis().set_visible(False)
axG.get_yaxis().set_visible(False)

# Label Plot
#axG.text(0.01,1.0,'D',transform = axG.transAxes,fontsize=24)

# Set labels and tick sizes
axG.set_xlabel(r'',fontsize=18)
axG.set_ylabel(r'',fontsize=18)
axG2.set_ylabel(r'Velocity',fontsize=18)
axG.xaxis.set_ticklabels([])
axG.yaxis.set_ticklabels([])
axG2.xaxis.set_ticklabels([])
axG2.yaxis.set_ticklabels([])

# Turns off chart clutter

# Turn off tick marks
axG.get_xaxis().set_ticks([])
axG.get_yaxis().set_ticks([])
axG2.get_xaxis().set_ticks([])
axG2.get_yaxis().set_ticks([])


# Turn off top and right splines
axG.spines["top"].set_visible(False)
axG.spines["left"].set_visible(False)
axG2.spines["top"].set_visible(False)
axG2.spines["left"].set_visible(False)

# Plot

axG.plot(p4316['Time'] - p4316['Time'][0],p4316['On_Board'] - p4316['On_Board'][0],color='k',linewidth=2,
        label='p4316')

velocity = rslope(p4316['Time'].reshape(p4316['Time'].size),p4316['On_Board'].reshape(p4316['On_Board'].size),window)
axG2.plot(p4316['Time'] - p4316['Time'][0],velocity,color='0.7',linewidth=2,
        linestyle=':',label='p4316')

# Set Velocity Axis Limits
axG2.set_ylim(0,2100)

# Add annotation for velocity
axG2.text(3,1900,r'1900 $\mu m/s$',fontsize=10)

# Add scale bars
# mu_ref = 0.07 * (axG.get_ylim()[1] - axG.get_ylim()[0]) + axG.get_ylim()[0]
# dmu = 0.1 * (axG.get_ylim()[1] - axG.get_ylim()[0])
# t_ref = 0.6 * (axG.get_xlim()[1] - axG.get_xlim()[0]) + axG.get_xlim()[0]
# dt = 1.
# d_ref = 0.07 * (axG2.get_ylim()[1] - axG2.get_ylim()[0]) + axG2.get_ylim()[0]
# ddis = 0.1 * (axG2.get_ylim()[1] - axG2.get_ylim()[0])
#
# axG.plot([t_ref,t_ref+dt],[mu_ref,mu_ref],color='k')
# axG.plot([t_ref,t_ref],[mu_ref,mu_ref+dmu],color='k')
# axG2.plot([t_ref+dt,t_ref+dt],[d_ref,d_ref+ddis],color='k')
# axG.text(t_ref+0.3*dt,mu_ref-0.5*dmu,'%d s'%dt,fontsize=8)
# axG.text(t_ref-3.0*dt,mu_ref+0.5*dmu,'%.4f $\mu$'%dmu,fontsize=8)
# axG2.text(t_ref+1.1*dt,d_ref+0.5*ddis,'%.1f $\mu m$'%ddis,fontsize=8)

# Add starting disp text to plot
#axG.text(0,-0.07,'%.1f mm' %(p4316['LP_Disp'][0]/1000.),fontsize=10,transform = axG.transAxes)


plt.savefig('event_zooms.png', bbox_inches="tight")
