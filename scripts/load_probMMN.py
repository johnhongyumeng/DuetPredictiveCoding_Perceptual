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
# Specify the file path
current_dir = os.getcwd()
parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
temp_probs = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

main_model = []
Output_dic = []
PltTool = []

#oddball_path =  '\\net\\Output_oddball_noHomeo.pkl'
#omissctrl_path='\\net\\Output_oddball_noHomeo_omissctrl.pkl'

# Loop through the temp_prob values to load each file into the lists
for idx, temp_prob in enumerate(temp_probs):
    file_path = f"{parent_dir}\\net\\prop\\Output_oddball.pklprob{temp_prob:.1f}"
    with open(file_path, 'rb') as f:
        model = pickle.load(f)
        output = pickle.load(f)
        plt_tool = pickle.load(f)
        main_model.append(model)
        Output_dic.append(output)
        PltTool.append(plt_tool)
        print(f"Successfully loaded data for temp_prob = {temp_prob:.1f} (Index: {idx})")

    
# Test # So actually this is right. 
PltTool[2].PlotPyrTimes_pop_all_long(left=50, right=59.9)   # Requires param_stim['choices']
#PltTool[2].PlotPyrPop_prob(left=-4,right=5) debug
PltTool[2].time_draw_start=1
PltTool[2].PlotPyrPop_prob(left=50, right=59.9) 

PltTool[2].PlotPyr_conductance(shift_id=-model.N_eachC//2)

PltTool[2].PlotInt_prob(left=50, right=59.9)  


PltTool[0].PlotPyrTimes_pop_all_prob(Output_dic,PltTool,decay_rate=0.600)   # Requires param_stim['choices']
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
