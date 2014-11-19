import matplotlib.pyplot as plt
import pandas as pd

exps = ['p4267','p4268','p4269','p4270','p4271','p4272','p4273',
        'p4309','p4310','p4311','p4312','p4313','p4314','p4316','p4317',
        'p4327','p4328','p4329','p4330']

exp_list = pd.read_excel('../../experiment_list.md')

# Setup figure and axes
fig = plt.figure(figsize=(12,9))
ax1 = plt.subplot(111)

# Set labels and tick sizes
ax1.set_xlabel(r'Average LP Displacement [mm]',fontsize=18)
ax1.set_ylabel(r'Stiffness [1/um]',fontsize=18)
ax1.tick_params(axis='both', which='major', labelsize=16)

# Turns off chart clutter

# Turn off top and right tick marks
ax1.get_xaxis().tick_bottom()
ax1.get_yaxis().tick_left()

# Turn off top and right splines
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)

# Plotting

for exp in exps:

    details = exp_list[exp_list['Experiment']==exp]
    nor_stress = details['Normal Stress']

    df = pd.read_csv('/Users/jleeman/Dropbox/PennState/BiaxExperiments/%s/%s_stiffness_cycles.txt'%(exp,exp))

    temp = df[df['Behavior']=='stable']
    ax1.scatter(temp['AvgDisp']/1000.,temp['Slope'],color='g',s=(nor_stress*1)**2,alpha=0.6)

    temp = df[df['Behavior']=='slow']
    ax1.scatter(temp['AvgDisp']/1000.,temp['Slope'],color='y',s=(nor_stress*1)**2,alpha=0.6)

    temp = df[df['Behavior']=='fast']
    ax1.scatter(temp['AvgDisp']/1000.,temp['Slope'],color='r',s=(nor_stress*1)**2,alpha=0.6)

# Set limits
ax1.set_xlim(0,50)
ax1.set_ylim(0,0.006)

# Plot Kc
df = pd.read_excel('/Users/jleeman/Dropbox/PennState/BiaxExperiments/p4309/p4309_rsf_fits.xlsx')




for i,fit in df.iterrows():

    if fit['Grade'] == 'A':
        color='#000066'
    elif fit['Grade'] == 'B':
        color='#0066CC'
    elif fit['Grade'] == 'C':
        color='#00CCFF'
        color='#000000'
    elif fit['Grade'] == 'D':
        color='#00FFFF'
        color='#000000'

    if fit['Type']=='Down' and fit['Law']=='r' and fit['k']==0.0055:
        ax1.scatter(fit['LP_Disp']/1000.,fit['Kc'],c=color,s=60,marker='v')

    elif fit['Type']=='Up' and fit['Law']=='r' and fit['k']==0.0055:
        ax1.scatter(fit['LP_Disp']/1000.,fit['Kc'],c=color,s=60,marker='^')

    else:
        pass

plt.savefig('Stiffness_Evolution.pdf',bbox_inches='tight')
