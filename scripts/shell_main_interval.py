# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 15:39:04 2024

@author: yawnt
"""
import numpy as np

from class_main_oddball import oddball_model 



#dev_degree_array=np.array([0,45,90,135,180,225,270,315])
#shift_axis=np.append(np.arange(int(2),int(8)),np.arange(int(2)))

#dev_degree_array=np.array([0])
dev_degree_array=np.arange(0,360,22.5)
shift_axis=np.append(np.arange(int(4),int(16)),np.arange(int(4)))

interval_array=np.arange(0.25,2.1,0.25)
Tau_int=1

output_dic_list=[]
DD_list=[]

for time_interval in interval_array:
#for dev_degree in dev_degree_array:

    main_model=oddball_model()
    # output_dic=main_model.main(dev_degree=dev_degree)
    output_dic=main_model.main(Tinter=time_interval, Tau_int=Tau_int)
    PltTool=main_model.selectAnalyz(Plot_flag=False)
    DD_value=PltTool.PlotPyrTimes_pop_all()

#    output_dic_list.append(output_dic)
    DD_list.append(DD_value)

import matplotlib.pyplot as plt  # Import matplotlib for plotting
DD_array=np.array(DD_list, dtype=float)

fig, axes=plt.subplots(1,1,figsize=(2.7,1.8))
axes.plot(interval_array, DD_array, marker='o',  color='b')
axes.plot(interval_array, DD_array, linestyle='-', color='black')

axes.set_xlabel('Temporal Interval')
axes.set_xlim(left=0)
axes.set_ylim(top=.41)

axes.set_ylabel('DD value')
fig.savefig('../figs/loopedOddball/'+'PopDDoverTemporalDif'+f'tau{Tau_int}' +'.jpg',bbox_inches='tight')
fig.savefig('../figs/loopedOddball/'+'PopDDoverTemporalDif'+f'tau{Tau_int}'+ '.svg',bbox_inches='tight')
plt.show()
'''
fig, axes=plt.subplots(1,1,figsize=(2.7,1.8))
axes.plot((dev_degree_array-180)/2, DD_array[shift_axis], marker='o', linestyle='-', color='b')
axes.plot((dev_degree_array-180)/2, DD_array[shift_axis], linestyle='-', color='black')

axes.set_xlabel('Tuning Degree Dif')
axes.set_ylabel('DD value')
fig.savefig('../figs/loopedOddball/'+'PopDDoverDegreeDif'+ '.jpg',bbox_inches='tight')
fig.savefig('../figs/loopedOddball/'+'PopDDoverDegreeDif'+ '.svg',bbox_inches='tight')
'''
