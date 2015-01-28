import matplotlib.pyplot as plt
import numpy as np
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

path = '/Users/jleeman/Dropbox/PennState/BiaxExperiments'
dis_low = 9*1000
dis_high = 17*1000
p4343 = ReadExp('p4343',path,dis_low,dis_high)
p4344 = ReadExp('p4344',path,dis_low,dis_high)
p4345 = ReadExp('p4345',path,dis_low,dis_high)
p4346 = ReadExp('p4346',path,dis_low,dis_high)
p4347 = ReadExp('p4347',path,dis_low,dis_high)
p4348 = ReadExp('p4348',path,dis_low,dis_high)
p4342 = ReadExp('p4342',path,dis_low,dis_high)
p4350 = ReadExp('p4350',path,dis_low,dis_high)
p4351 = ReadExp('p4351',path,dis_low,dis_high)

# These are the "Tableau 20" colors as RGB.
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]

# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.
for i in range(len(tableau20)):
    r, g, b = tableau20[i]
    tableau20[i] = (r / 255., g / 255., b / 255.)

# Setup figure and axes
# Generally plots is ~1.33x width to height (10,7.5 or 12,9)
fig = plt.figure(figsize=(24,18))
ax1 = plt.subplot(111)

# Set labels and tick sizes
ax1.set_xlabel(r'Load Point Displacement [mm]',fontsize=32)
ax1.set_ylabel(r'Friction',fontsize=32)
ax1.tick_params(axis='both', which='major', labelsize=24)

# Turns off chart clutter

# Turn off top and right tick marks
ax1.get_xaxis().tick_bottom()
ax1.get_yaxis().tick_left()
ax1.get_yaxis().set_ticks([])

# Turn off top and right splines
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)

# Plotting
window_size = 51
order = 5
ax1.plot(p4343['LP_Disp']/1000.,savitzky_golay(np.ravel(p4343['mu']), window_size, order),label='6 MPa',color=tableau20[0])
ax1.plot(p4344['LP_Disp']/1000.,savitzky_golay(np.ravel(p4344['mu']), window_size, order)+0.05,label='7 MPa',color=tableau20[2])
ax1.plot(p4345['LP_Disp']/1000.,savitzky_golay(np.ravel(p4345['mu']), window_size, order)+0.05*2,label='8 MPa',color=tableau20[4])
ax1.plot(p4346['LP_Disp']/1000.,savitzky_golay(np.ravel(p4346['mu']), window_size, order)+0.05*3,label='9 MPa',color=tableau20[6])
ax1.plot(p4347['LP_Disp']/1000.,savitzky_golay(np.ravel(p4347['mu']), window_size, order)+0.05*4,label='10 MPa',color=tableau20[8])
ax1.plot(p4348['LP_Disp']/1000.,savitzky_golay(np.ravel(p4348['mu']), window_size, order)+0.05*5,label='11 MPa',color=tableau20[10])
ax1.plot(p4342['LP_Disp']/1000.,savitzky_golay(np.ravel(p4342['mu']), window_size, order)+0.05*6,label='12 MPa',color=tableau20[12])
ax1.plot(p4350['LP_Disp']/1000.,savitzky_golay(np.ravel(p4350['mu']), window_size, order)+0.05*7,label='13 MPa',color=tableau20[14])
ax1.plot(p4351['LP_Disp']/1000.,savitzky_golay(np.ravel(p4351['mu']), window_size, order)+0.05*8,label='14 MPa',color=tableau20[18])

#ax1.text(9.1,np.mean(p4343['mu'])+0.05*0.25,r'$\mu = $%0.2f' %np.mean(p4343['mu']),fontsize=16,color=tableau20[0])
#ax1.text(9.1,np.mean(p4344['mu'])+0.05*1.25,r'$\mu = $%0.2f' %np.mean(p4344['mu']),fontsize=16,color=tableau20[2])
#ax1.text(9.1,np.mean(p4345['mu'])+0.05*2.25,r'$\mu = $%0.2f' %np.mean(p4345['mu']),fontsize=16,color=tableau20[4])
#ax1.text(9.1,np.mean(p4346['mu'])+0.05*3.25,r'$\mu = $%0.2f' %np.mean(p4346['mu']),fontsize=16,color=tableau20[6])
#ax1.text(9.1,np.mean(p4347['mu'])+0.05*4.25,r'$\mu = $%0.2f' %np.mean(p4347['mu']),fontsize=16,color=tableau20[8])
#ax1.text(9.1,np.mean(p4348['mu'])+0.05*5.25,r'$\mu = $%0.2f' %np.mean(p4348['mu']),fontsize=16,color=tableau20[10])
#ax1.text(9.1,np.mean(p4342['mu'])+0.05*6.25,r'$\mu = $%0.2f' %np.mean(p4342['mu']),fontsize=16,color=tableau20[12])
#ax1.text(9.1,np.mean(p4350['mu'])+0.05*7.25,r'$\mu = $%0.2f' %np.mean(p4350['mu']),fontsize=16,color=tableau20[14])
#ax1.text(9.1,np.mean(p4351['mu'])+0.05*8.25,r'$\mu = $%0.2f' %np.mean(p4351['mu']),fontsize=16,color=tableau20[18])

ax1.text(17.1,np.mean(p4343['mu'])+0.05*0,r'$\sigma_n$ = 6 MPa',fontsize=22,color=tableau20[0])
ax1.text(17.1,np.mean(p4344['mu'])+0.05*1,r'$\sigma_n$ = 7 MPa',fontsize=22,color=tableau20[2])
ax1.text(17.1,np.mean(p4345['mu'])+0.05*2,r'$\sigma_n$ = 8 MPa',fontsize=22,color=tableau20[4])
ax1.text(17.1,np.mean(p4346['mu'])+0.05*3,r'$\sigma_n$ = 9 MPa',fontsize=22,color=tableau20[6])
ax1.text(17.1,np.mean(p4347['mu'])+0.05*4,r'$\sigma_n$ = 10 MPa',fontsize=22,color=tableau20[8])
ax1.text(17.1,np.mean(p4348['mu'])+0.05*5,r'$\sigma_n$ = 11 MPa',fontsize=22,color=tableau20[10])
ax1.text(17.1,np.mean(p4342['mu'])+0.05*6,r'$\sigma_n$ = 12 MPa',fontsize=22,color=tableau20[12])
ax1.text(17.1,np.mean(p4350['mu'])+0.05*7,r'$\sigma_n$ = 13 MPa',fontsize=22,color=tableau20[14])
ax1.text(17.1,np.mean(p4351['mu'])+0.05*8,r'$\sigma_n$ = 14 MPa',fontsize=22,color=tableau20[18])

ax1.text(17.1,np.mean(p4343['mu'])+0.05*0-0.01,r'p4343',fontsize=14,color=tableau20[0])
ax1.text(17.1,np.mean(p4344['mu'])+0.05*1-0.01,r'p4344',fontsize=14,color=tableau20[2])
ax1.text(17.1,np.mean(p4345['mu'])+0.05*2-0.01,r'p4345',fontsize=14,color=tableau20[4])
ax1.text(17.1,np.mean(p4346['mu'])+0.05*3-0.01,r'p4346',fontsize=14,color=tableau20[6])
ax1.text(17.1,np.mean(p4347['mu'])+0.05*4-0.01,r'p4347',fontsize=14,color=tableau20[8])
ax1.text(17.1,np.mean(p4348['mu'])+0.05*5-0.01,r'p4348',fontsize=14,color=tableau20[10])
ax1.text(17.1,np.mean(p4342['mu'])+0.05*6-0.01,r'p4342',fontsize=14,color=tableau20[12])
ax1.text(17.1,np.mean(p4350['mu'])+0.05*7-0.01,r'p4350',fontsize=14,color=tableau20[14])
ax1.text(17.1,np.mean(p4351['mu'])+0.05*8-0.01,r'p4351',fontsize=14,color=tableau20[18])




# Set limits
ax1.set_xlim(9,17)
#ax1.set_ylim(0.65,1.12)

plt.savefig('AllNormal_Runplots.pdf', bbox_inches="tight")
