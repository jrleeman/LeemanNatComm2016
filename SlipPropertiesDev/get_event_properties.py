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

    def plot(self,data,plot_name):
        ax1 = plt.subplot(111)
        ax1.plot(data['LP_Disp'][self.start_row:self.end_row],data['mu'][self.start_row:self.end_row],color='k')
        ax1.axvline(x=data['LP_Disp'][self.start_row+self.npts_stiffness_fit],color='b')
        x_stiffness = [data['LP_Disp'][self.start_row],data['LP_Disp'][self.end_row]]
        y_stiffness = np.polyval([self.stiffness,self.stiffness_intercept],x_stiffness)
        ax1.plot(x_stiffness,y_stiffness,color='r',linewidth=2)
        ax1.set_xlabel(r'LP Displacement [$\mu m$]')
        ax1.set_ylabel(r'Friction')
        plt.title('%s'%plot_name)
        plt.savefig('%s.png'%plot_name)
        plt.clf()

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

    event_list.append(slip_event)

outfile = open('%s_event_properties.txt'%experiment,'w')
outfile.write('Event,StartRow,FailRow,EndRow,SlipDuration,Stiffness,NptsStiffness,StiffnessIntercept,FailTime,FailDisplacement\n')
for i,event in enumerate(event_list):
    try:
        outfile.write('%d,%d,%d,%d,%f,%f,%d,%f,%f,%f\n'%(i,event.start_row,event.failure_row,event.end_row,event.slip_duration,event.stiffness,event.npts_stiffness_fit,event.stiffness_intercept,event.failure_time,event.failure_displacement))
    except:
        print i,event.start_row,event.failure_row,event.end_row,event.slip_duration,event.stiffness,event.npts_stiffness_fit,event.stiffness_intercept,event.failure_time,event.failure_displacement
        break
    event.plot(data,'%s_%d'%(experiment,i))
outfile.close()
