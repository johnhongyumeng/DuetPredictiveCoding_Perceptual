# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 15:39:04 2024

@author: yawnt
"""
import numpy as np

import matplotlib.pyplot as plt  # Import matplotlib for plotting

from class_main_oddball import oddball_model 


plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False
#dev_degree_array=np.array([0,45,90,135,180,225,270,315])
#shift_axis=np.append(np.arange(int(2),int(8)),np.arange(int(2)))

#dev_degree_array=np.array([0])
#dev_degree_array=np.arange(0,360,22.5)
#shift_axis=np.append(np.arange(int(4),int(16)),np.arange(int(4)))
#interval_array=np.arange(0.5,3.1,0.5)
#interval_array=np.arange(0.1,0.11,0.1)   # Testing omission
#interval_array=np.array([0.1,0.125,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8])
interval_array=np.array([0.2])

output_dic_list=[]
DD_list=[]
nOmiss_list=[]

for time_interval in interval_array:
    main_model=oddball_model()
    output_dic=main_model.main(Tinter=time_interval,Tstim=0.1,
                               omission_flag=True,saving_tag=f"oddball_omission_int{time_interval:.2f}",
                               Tau_int=0.1,g_IntE=6,Stim_g=1.5,
                               Trained_Pathway='\\net\\Output_noHomeoN200.pkl', flag_homeo=False)
    DD_value,nOmiss_ite=main_model.analyze(Plot_flag=True)
    output_dic_list.append(output_dic)
    DD_list.append(DD_value)
    nOmiss_list.append(nOmiss_ite)
    
DD_array=np.array(DD_list, dtype=float)
fig, axes=plt.subplots(1,1,figsize=(2.7,1.8))
axes.plot(interval_array, DD_array, marker='o',  color='b')
axes.plot(interval_array, DD_array, linestyle='-', color='k')

axes.set_xlabel('Temporal Interval')
axes.set_ylabel('DD value')
fig.show()
#fig.savefig('../figs/looped/'+'PopDDoverInterval'+ '.jpg',bbox_inches='tight')
#fig.savefig('../figs/looped/'+'PopDDoverInterval'+ '.svg',bbox_inches='tight')

fig, axes=plt.subplots(1,1,figsize=(2.7,1.8))
nOmiss_array=np.array(nOmiss_list, dtype=float)
axes.plot(interval_array, nOmiss_array, marker='o', linestyle='-', color='b')
axes.plot(interval_array, nOmiss_array, linestyle='-', color='k')

axes.set_xlabel('Temporal Interval')
axes.set_ylabel('Num of Omiss Neu')
#fig.savefig('../figs/looped/'+'nOmissoverInterval'+ '.jpg',bbox_inches='tight')
#fig.savefig('../figs/looped/'+'nOmissoverInterval'+ '.svg',bbox_inches='tight')



# Here are sth about reploting the data.
# Data for ISI values
isi_values = [0.125, 0.25, 0.5]
# Experimental measurements (%f values)
exp_values = [10.94, 5.74, 2.15]

# Modeling results (r values)
model_values = [0.925, 0.403, 0.090]
# Create a figure and axis

fig, ax1 = plt.subplots(figsize=(2.7, 1.8))

# Define bar width and positions
bar_width = 0.3
x_positions = np.arange(len(isi_values))

# Plot experimental data (left y-axis)
ax1.bar(x_positions - bar_width/2, exp_values, width=bar_width, color='black', label='Experimental (%f)')
ax1.set_ylabel('Exp. (% Seqs)', color='black')
ax1.tick_params(axis='y', labelcolor='black')

# Create second y-axis for modeling results
ax2 = ax1.twinx()
ax2.bar(x_positions + bar_width/2, model_values, width=bar_width, color='#E34A33', label='Modeling (r)')
ax2.set_ylabel('Model (rate)', color='#E34A33')
ax2.tick_params(axis='y', labelcolor='#E34A33')
ax2.spines['right'].set_visible(True)

# Set x-axis labels
ax1.set_xticks(x_positions)
ax1.set_xticklabels(isi_values)
#ax1.set_xlabel('ISI')

# Title and layout adjustments
#plt.title('Comparison of Experimental and Modeling Results')
fig.tight_layout()
plt.show()
fig.savefig('../figs/looped/'+'PopDDoverISIModelwExp'+ '.jpg',bbox_inches='tight')
fig.savefig('../figs/looped/'+'PopDDoverISIModelwExp'+ '.svg',bbox_inches='tight')


# Now let's redo the ploting. Based on PlotPyrTimes_col

Exp_iOmi = np.array([1, 1, 7, 1, 2, 2, 6, 15, 10, 11, 19, 18, 10, 7, 5, 2, 1, 1])
sum_count= np.sum(Exp_iOmi)
bin_edges = np.arange(-1, 1.01, 1.0/9)

fig, axes=plt.subplots(1,1,figsize=(2,2))
# Compute the histogram data
#counts, edges = np.histogram(Exp_iOmi, bins=bin_edges, density=True)

# Plot the histogram manually
for i in range(len(bin_edges) - 1):
    # Get the bin's start, end, and count
    left = bin_edges[i]
    right = bin_edges[i + 1]
    height = Exp_iOmi[i]/sum_count*9  # normalized by sum count and bin width

    # Set the color based on the bin height
    facecolor = 'red' if left > 0.3 else 'grey'
    # Plot the bar
    axes.bar((left+right)/2, height, width=right - left, color=facecolor, alpha=0.5, edgecolor='black', linewidth=0.5)


axes.axvline(0.0771, color='grey', linestyle='--', linewidth=1.5)
plt.show()

fig.savefig('../figs/oddball_omission/'+'iOmidistExp'+ '.jpg',bbox_inches='tight')
fig.savefig('../figs/oddball_omission/'+'iOmidistExp'+ '.svg',bbox_inches='tight')
axes.set_title('iOmidist')
