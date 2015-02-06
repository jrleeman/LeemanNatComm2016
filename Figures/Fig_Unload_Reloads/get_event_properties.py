import pylab as plt
import sys
import numpy as np
from biaxread import *
from math import sqrt

def CalcCorrelation(x,y):

    # Check that both arrays are the same length
    if len(x) != len(y):
        print "Error! Displacement and stress must be the same length"
        return 0


    npts = len(x)
    n_min = 50 # Require MORE THAN n_min points to fit

    # Make the sums we will need
    sum_x = sum(x[0:n_min])
    sum_y = sum(y[0:n_min])
    sum_x_sq = sum(x[0:n_min]*x[0:n_min])
    sum_y_sq = sum(y[0:n_min]*y[0:n_min])
    sum_xy = sum(x[0:n_min]*y[0:n_min])

    correlations = []
    for i in range(n_min+1):
        correlations.append(0)

    for i in range(n_min,npts):

        n = i+1

        # Advance sums by one point
        sum_x = sum_x + x[i]
        sum_y = sum_y + y[i]
        sum_x_sq = sum_x_sq + x[i]**2
        sum_y_sq = sum_y_sq + y[i]**2
        sum_xy = sum_xy + x[i]*y[i]

        correlation = (n * sum_xy - sum_x * sum_y) / (sqrt(n*sum_x_sq - sum_x*sum_x) * sqrt(n * sum_y_sq - sum_y*sum_y))
        correlations.append(correlation)

    return correlations

def FitBestCorrelation(x,y,correlations):

    max_correlation = max(correlations)
    npts_bestfit = np.where(correlations==max_correlation)[0][0]

    #print npts_bestfit

    x = x[0:npts_bestfit]
    y = y[0:npts_bestfit]

    #print len(x),len(y)

    coeffs = np.polyfit(x,y,1)

    return npts_bestfit, coeffs

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

class SlipEvent():

    def __init__(self):
        self.start_row = None
        self.failure_row = None
        self.end_row = None
        self.stiffness = None
        self.stiffness_intercept = None
        self.friction_drop = None
        self.slip_duration = None
        self.failure_time = None
        self.failure_displacement = None
        self.npts_stiffness_fit = None
        self.friction_25_percent = None
        self.friction_75_percent = None
        self.friction_25_percent_idx = None
        self.friction_75_percent_idx = None
        self.velocity = None
        self.mean_middle_50_percent_velocity =None
        self.stiffness_correlation = None
        self.peak_velocity = None

    def plot(self,data,plot_name):
        ax1 = plt.subplot(111)
        ax2 = ax1.twinx()

        # Plot friction Data
        ax1.plot(data['LP_Disp'][self.start_row:self.end_row],data['mu'][self.start_row:self.end_row],color='k')

        # Plot velocity Data
        print np.shape(data['LP_Disp'][self.start_row:self.end_row]),np.shape(self.velocity)
        ax2.plot(data['LP_Disp'][self.start_row:self.end_row],self.velocity,color='b')

        # Plot stiffness fit and fit limit
        ax1.axvline(x=data['LP_Disp'][self.start_row+self.npts_stiffness_fit],color='r')

        x_stiffness = [data['LP_Disp'][self.start_row],data['LP_Disp'][self.end_row]]
        y_stiffness = np.polyval([self.stiffness,self.stiffness_intercept],x_stiffness)
        ax1.plot(x_stiffness,y_stiffness,color='r',linewidth=2)

        # Plot max, min, and quartile friction values
        ax1.scatter(data['LP_Disp'][self.start_row],data['mu'][self.start_row],color='g')
        ax1.scatter(data['LP_Disp'][self.failure_row],data['mu'][self.failure_row],color='r')
        ax1.scatter(data['LP_Disp'][self.end_row],data['mu'][self.end_row],color='g')

        #ax1.scatter(data['LP_Disp'][self.friction_75_percent_idx],data['mu'][self.friction_75_percent_idx],color='m')
        #ax1.scatter(data['LP_Disp'][self.friction_25_percent_idx],data['mu'][self.friction_25_percent_idx],color='m')

        ax1.axhline(y=self.friction_25_percent,color='m')
        ax1.axhline(y=self.friction_75_percent,color='m')

        ax1.set_xlabel(r'LP Displacement [$\mu m$]')
        ax1.set_ylabel(r'Friction')

        plt.title('%s'%plot_name)
        plt.savefig('%s.png'%plot_name)
        plt.clf()

    def velocity_analysis(self,data):
        x = data['Time'][self.start_row:self.end_row]
        y = data['OB_Top'][self.start_row:self.end_row]
        mu = data['mu'][self.start_row:self.end_row]
        x = x.flatten()
        y = y.flatten()

        window_size = 11

        #print failure_row,ending_row
        #print np.shape(x),np.shape(y)

        self.velocity = rslope(x,y,window_size)

        minimum_friction = data['mu'][self.end_row]
        maximum_friction = data['mu'][self.failure_row]

        self.friction_75_percent = 0.75*(maximum_friction - minimum_friction) + minimum_friction
        self.friction_25_percent = 0.25*(maximum_friction - minimum_friction) + minimum_friction

        self.friction_75_percent_idx = find_nearest(data['mu'][self.failure_row:self.end_row],self.friction_75_percent)
        self.friction_25_percent_idx = find_nearest(data['mu'][self.failure_row:self.end_row],self.friction_25_percent)

        self.mean_middle_50_percent_velocity = np.mean(self.velocity[self.friction_75_percent_idx:self.friction_25_percent_idx])

        self.peak_velocity = max(self.velocity)

experiment = sys.argv[1]

# Load the biax data from the given project folder
path = '/Users/jleeman/Dropbox/PennState/BiaxExperiments/'

events = np.loadtxt('%s_event_rows.txt'%experiment,skiprows=1,delimiter=',')

# Get experiment name and read data
data = ReadAscii(path + '%s/%s_data.txt' %(experiment,experiment))

# We can't use the first event, because we don't have an end point
# of an earlier event for stiffness calculation. So, let's make a
# list of events.
event_list = []
for i in np.arange(1,len(events)):
    slip_event = SlipEvent()

    # Assign rows for event parts
    slip_event.start_row = events[i-1][1]
    slip_event.failure_row= events[i][0]
    slip_event.end_row = events[i][1]

    # Assign time and displacement to failure
    slip_event.failure_time = data['Time'][slip_event.failure_row]
    slip_event.failure_displacement = data['LP_Disp'][slip_event.failure_row]

    # Assign basic mechanical properties
    slip_event.friction_drop = data['mu'][slip_event.failure_row] - data['mu'][slip_event.end_row]
    slip_event.slip_duration = data['Time'][slip_event.end_row] - data['Time'][slip_event.failure_row]

    # Assign quartile portions, velocity, etc
    try:
        slip_event.velocity_analysis(data)
        print i,slip_event.mean_middle_50_percent_velocity
    except:
        pass


    # Calculate the stiffness
    try:
        correlation_results = CalcCorrelation(np.ravel(data['LP_Disp'][slip_event.start_row:slip_event.failure_row]),np.ravel(data['mu'][slip_event.start_row:slip_event.failure_row]))
        npts_bestfit,coeffs = FitBestCorrelation(np.ravel(data['LP_Disp'][slip_event.start_row:slip_event.failure_row]),np.ravel(data['mu'][slip_event.start_row:slip_event.failure_row]),correlation_results)
    except:
        print "Stiffness calculation bombed!"
        coeffs = [0,0]
        npts_bestfit = 0

    slip_event.stiffness = coeffs[0]
    slip_event.stiffness_intercept = coeffs[1]
    slip_event.npts_stiffness_fit = npts_bestfit

    if slip_event.failure_row < slip_event.end_row:
        event_list.append(slip_event)
    else:
        pass

outfile = open('%s_event_properties.txt'%experiment,'w')
outfile.write('Event,StartRow,FailRow,EndRow,SlipDuration,Stiffness,NptsStiffness,StiffnessIntercept,FailTime,FailDisplacement,50_Mean_Velocity,Peak_Velocity\n')
for i,event in enumerate(event_list):
    try:
        outfile.write('%d,%d,%d,%d,%f,%f,%d,%f,%f,%f,%f,%f\n'%(i,event.start_row,event.failure_row,event.end_row,event.slip_duration,event.stiffness,event.npts_stiffness_fit,event.stiffness_intercept,event.failure_time,event.failure_displacement,event.mean_middle_50_percent_velocity,event.peak_velocity))
    except:
        print i,event.start_row,event.failure_row,event.end_row,event.slip_duration,event.stiffness,event.npts_stiffness_fit,event.stiffness_intercept,event.failure_time,event.failure_displacement,event.mean_middle_50_percent_velocity,event.peak_velocity
        break
    event.plot(data,'%s_plots/%s_%d'%(experiment,experiment,i))
outfile.close()
