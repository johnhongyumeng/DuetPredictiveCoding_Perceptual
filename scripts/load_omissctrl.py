# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 11:15:10 2024
This script is used to plot the omiss ctrl vs the omission paradigm.
@author: yawnt
"""
import pickle
import os
import sys
from class_main_oddball import oddball_model 

# Specify the file path
current_dir = os.getcwd()
parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))

oddball_path =  '\\net\\Output_oddball.pkl'
omissctrl_path='\\net\\Output_oddball_omissctrl.pkl'

#with open(parent_dir+'\\net\\Output_oddball_noHomeo.pkl', 'wb') as f:
#with open(parent_dir+'\\net\\Output_oddball_noHomeo_omissctrl.pkl', 'wb') as f:

#oddball_path =  '\\net\\Output_oddball_noHomeo.pkl'
#omissctrl_path='\\net\\Output_oddball_noHomeo_omissctrl.pkl'

# Open the file and load the objects
with open(parent_dir+oddball_path, 'rb') as f:
    main_model_odd= pickle.load(f)
    Output_dic_odd = pickle.load(f)
    PltTool_odd = pickle.load(f)
    
with open(parent_dir+omissctrl_path, 'rb') as f:
    main_model_ctrl= pickle.load(f)
    Output_dic_ctrl = pickle.load(f)
    PltTool_ctrl = pickle.load(f)
    
#PltTool_odd.PlotPyrDDInd_all()   sanity test
PltTool_odd.PlotPyrDDInd_all() 
PltTool_odd.PlotPyrDDOmiss_all(PltTool_ctrl.Odd_rate_vec,linewidth=1,top=10)
#Plot_nPE_barplot


# Testing the parres et al. 2017 Not good. Maybe due to the tuning curve differences.
#PltTool_odd.PlotPyrTimes_all_parras(thr=1.)   # New analysis on generating the distribution.
#PltTool_odd.PlotPyrTimes_col()
