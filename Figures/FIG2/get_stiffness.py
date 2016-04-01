import numpy as np
import matplotlib.pyplot as plt
import sys

def load_event_properties(experiment):
    """
    Load event property file picks for a given experiment number and return
    that data as an array
    """
    return np.loadtxt('%s_event_properties.txt'%experiment,delimiter=',',skiprows=1)

def load_blacklist(experiment):
    """
    Load event numbers from the blacklist file for each experiment and
    return them as an array
    """
    blacklist = np.loadtxt('%s_blacklist.txt'%experiment)
    return blacklist

def load_events(experiment):
    """
    Loads all events from a given experiment that are not on the blacklist
    file for that experiment. Returns array of event properties.
    """
    event_properties = load_event_properties(experiment)
    blacklist = load_blacklist(experiment)
    return np.delete(event_properties,blacklist,axis=0)

#Event,StartRow,FailRow,EndRow,SlipDuration,Stiffness,NptsStiffness,StiffnessIntercept,FailTime,FailDisplacement,50_Mean_Velocity,Peak_Velocity

normal = []
stiffness = []

experiments = ['p4351','p4350','p4342','p4348','p4347','p4346','p4345','p4344','p4343']
sigma_n=[14,13,12,11,10,9,8,7,6]

for exp,stress in zip(experiments,sigma_n):
    print "Experiment: ",exp
    #data = np.loadtxt('%s_event_properties.txt'%exp,unpack=True,delimiter=',',skiprows=1,usecols=[9,5])
    data = load_events(exp)
    data = data[data[:,9]>=40000]
    data = data[data[:,9]<=50000]

    normal.append(stress)
    stiffness.append(np.mean(data[:,5]))
    plt.axhline(y=np.mean(data[:,5]))
    plt.scatter(data[:,0],data[:,5])
    plt.ylim(0,0.0007)
    #plt.show()

k_kc = np.array(stiffness)/7e-4
plt.clf()
print normal
print k_kc
plt.scatter(normal,k_kc,s=50)
plt.xlabel('Normal Stress')
plt.ylabel('k/kc')
plt.show()
