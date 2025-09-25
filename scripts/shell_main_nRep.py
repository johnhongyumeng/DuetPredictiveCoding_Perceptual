# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 15:39:04 2024

@author: yawnt
"""
import numpy as np

import matplotlib.pyplot as plt  # Import matplotlib for plotting

from class_main_oddball_prob_nRep_dev import oddball_model 

Tinter=0.8 
Tstim=0.2
Tau_int=0.5
g_IntE=15         # first bunch 15
g_back_int=44        # first bunch 44

nbr_rep_dev=1
prob_flag= False

plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False
#dev_degree_array=np.array([0,45,90,135,180,225,270,315])
nRep_array=np.arange(1,11)

MMN_pop_vec=[]
MMN_int_vec=[]

for nbr_rep_std in nRep_array:
    Input_blocker=nbr_rep_std
    main_model=oddball_model()
    Output_dic=main_model.main(omission_flag=False,saving_tag=f'oddball_nRep{nbr_rep_std}Tau_int{Tau_int}',adap_amount=0,
                               Tau_int=Tau_int, Tinter=Tinter,Tstim=Tstim, 
                               nbr_rep_std=nbr_rep_std,nbr_rep_dev=nbr_rep_dev,
                               g_IntE=g_IntE,g_back_int=g_back_int,prob_flag=prob_flag,
                               Input_blocker=Input_blocker)
    PltTool=main_model.selectAnalyz(Plot_flag=False)
    PltTool.PlotPyrPop(left=-1)
    MMN_Int=PltTool.PlotInt_nrep(nrep=Input_blocker)
    MMN_pop=np.mean(PltTool.Odd_rate_vec)-np.mean(PltTool.Last_rate_vec)

    MMN_pop_vec.append(MMN_pop)
    MMN_int_vec.append(MMN_Int)
    
MMN_pop_array=np.array(MMN_pop_vec, dtype=float)
MMN_int_array=np.array(MMN_int_vec, dtype=float)

fig, axes=plt.subplots(1,1,figsize=(2.7,1.8))
axes.plot(nRep_array, MMN_pop_array, marker='o',  color='b')
axes.plot(nRep_array, MMN_pop_array, linestyle='-', color='k')


axes.set_xlabel('Rep Num')
axes.set_ylabel('MMN')
axes.set_ylim(bottom=0)
axes.set_ylim(top=0.11)

fig.savefig('../figs/looped/'+f'MMNoverNumRepTau_int{Tau_int}'+ '.jpg',bbox_inches='tight')
fig.savefig('../figs/looped/'+f'MMNoverNumRepTau_int{Tau_int}'+ '.svg',bbox_inches='tight')

fig, axes=plt.subplots(1,1,figsize=(2.7,1.8))
axes.plot(nRep_array, MMN_int_array, marker='o',  color='b')
axes.plot(nRep_array, MMN_int_array, linestyle='-', color='k')

axes.set_xlabel('Rep Num')
axes.set_ylabel('Integrator Dif')
axes.set_ylim(top=0)
axes.set_ylim(bottom=-2.3)

fig.show()
fig.savefig('../figs/looped/'+f'MMNIntoverNumRepTau_int{Tau_int}'+ '.jpg',bbox_inches='tight')
fig.savefig('../figs/looped/'+f'MMNIntoverNumRepTau_int{Tau_int}'+ '.svg',bbox_inches='tight')
