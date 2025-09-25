# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 11:15:10 2024
This script is used to plot the omiss ctrl vs the omission paradigm.
@author: yawnt
"""
import pickle
import os
import sys
import numpy as np
from class_main_oddball_prob_manyctrl_normset import oddball_model 
#import class_main_oddball_prob_manyctrl_normset_realistic
current_dir = os.getcwd()
parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))

lib_dir = os.path.join(parent_dir, 'lib')
sys.path.append(lib_dir)
        
        

# Specify the file path
current_dir = os.getcwd()
parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))



file_path = f"{parent_dir}\\net\\realistic\\Output_oddball_normset_ctrl.pkl"
with open(file_path, 'rb') as f:
    main_ctrl = pickle.load(f)
    Output_ctrl = pickle.load(f)
    plt_tool_ctrl = pickle.load(f)
    print(f"Successfully loaded data fornormset_ctrl")

file_path = f"{parent_dir}\\net\\realistic\\Output_oddball_normset_realistic.pkl"



with open(file_path, 'rb') as f:
    model_real = pickle.load(f)
    output_real= pickle.load(f)
    plt_tool_real = pickle.load(f)
    print(f"Successfully loaded data fornormset_real")


# Now combining the two. 
plt_tool_ctrl.Plot_nRep(left=50, right=60)

plt_tool_real.Plot_real_odd(left=50, right=60)

plt_tool_real.Plot_real_odd_combine(plt_tool_ctrl)
plt_tool_real.std+plt_tool_real.N_eachC//2-1



plt_tool_ctrl.PlotPyrPop(left=50,right=60)
plt_tool_ctrl.PlotInt(left=50,right=60,labels=['x','y'])
plt_tool_ctrl.PlotPyr_pPExnPEy(left=50,right=60,labels=['pPE(x)','nPE(y)'])

plt_tool_real.PlotPyrPop(left=50,right=60)
plt_tool_real.PlotInt(left=50,right=60,labels=['x','y'])
plt_tool_real.PlotPyr_pPExnPEy(left=50,right=60,labels=['pPE(x)','nPE(y)'])


plt_tool_real.PlotInt(left=-3,right=10,labels=['x','y'])

'''
Rsquare_vec=[]
decay_rate_vec = np.arange(0.4, 0.65, 0.01)
for decay_rate in decay_rate_vec:
    Rtemp=PltTool[0].PlotPyrTimes_pop_all_prob(Output_dic,PltTool,decay_rate=decay_rate)   # Requires param_stim['choices']
    Rsquare_vec.append(Rtemp)

import matplotlib.pyplot as plt
fig, axes=plt.subplots(1,1,figsize=(2.7,1.8))  # sanity test on the firing rate
axes.plot(decay_rate_vec,Rsquare_vec,'k')
axes.plot(decay_rate_vec,Rsquare_vec,'.b')
'''

# Testing the parres et al. 2017 Not good. Maybe due to the tuning curve differences.
#PltTool_odd.PlotPyrTimes_all_parras(thr=1.)   # New analysis on generating the distribution.
#PltTool_odd.PlotPyrTimes_col()
