# -*- coding: utf-8 -*-
"""
Modified on 10042023. So this is a totally different thing. Hope this would work.
Created on 03312023 15:33:15 2023
Current version, suggest at a specific tuning curve point, different E-cells
@author: John Meng
"""
# Certain control parameters. Should set up better later.

flag_StoE=1  # Inhibition control from SST to pyramidal cells.
flag_EtoInt=1  # Bottom to top control 
adap_amount=0.3
omission_flag=0  # 0: turning down the deviation stimuli
Tinter=0.4
#saving_tag='RingGate/'   # Additonal saving tags. Can be adap, adaponly, Ring
#saving_tag='NoIntdown/'   # Additonal saving tags. Can be adap, adaponly, Ring


#saving_tag='mainOddball_omission/'   # Additonal saving tags. Can be adap, adaponly, Ring

#saving_tag='adap/' 
#='adaponly/' 
kappa=0.3             # default 0.5. Requires a fine tuning.

flag_basebooster=1   #Cheating one. Why the line is not right. Not sure.
Input_S_ratio=1

####################
# Same trivial setup for input and output function

import os
import sys
import importlib
from datetime import date,datetime
from tqdm import tqdm

current_dir = os.getcwd()
parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))

lib_dir = os.path.join(parent_dir, 'lib')
sys.path.append(lib_dir)


sav_dir = os.path.join(parent_dir, 'runs')

now=datetime.now()
today=date.today()
current_time=now.strftime("%H%M")
current_day=today.strftime("%m%d_")
timestamp=current_day+current_time

TempPath=sav_dir+'\\'+timestamp

if os.path.isdir(TempPath)==False:
    os.makedirs(TempPath)
    
FigPath=os.path.join(parent_dir+'\\figs\\mainOddball_omission')
if os.path.isdir(FigPath)==False:
    os.makedirs(FigPath)

filename = os.path.basename(__file__)
from shutil import copyfile      # Save a copy of the code to the saving folder.
copyfile(filename,TempPath+'\\'+filename)
para_loc='..\lib\\'

script_name=os.path.basename(__file__)
copyfile(para_loc+'parameters_tuningcurve.py',TempPath+'\parameters_tuningcurve.py')
copyfile(script_name, TempPath+script_name)

import basic_functions as bf
import parameters_tuningcurve as params
import data_structure
import update_functions
import numpy as np
from collections import OrderedDict
import matplotlib.pyplot as plt
import pickle

# All right, I really hate this. The work is not that transparaent. But it is one way
PARAS=params.PARAMS_ALL.copy()
PARAS['PARAMS_Pyr']['kappa_NMDA']=0
PARAS['PARAMS_INT']['kappa_NMDA']=0
params_E = PARAS['PARAMS_Pyr'].copy()
params_Int = PARAS['PARAMS_INT'].copy()
params_P = PARAS['PARAMS_PV'].copy()
params_S = PARAS['PARAMS_SST'].copy()

# list_param= np.array([-0.01,2,0.5,8,50,1,1,1,1,0.1,2])   

# adaptation, top-down feedback, Tinter, dev_iud, delay, 
# ndf_plasticity, int_plasticity. bool_sigma_ndf, bool_sigma_int, adaptation_timeconstant
# nbr_rep. This description is based on the dict_param in Script7

params_stim=PARAS['Stimulus'].copy()
params_stim['prob_std'] = 0.8


params_sim=PARAS['Simulation'].copy()
    
# The following par values are from the 1st parameter set in script7.    
delay=50*0.001   # The unit should be second
params_stim['Tinter'] = Tinter
params_stim['nbr_rep_std']=6

# Generate a stimulus. This may change into a class.

#nbr_stim_std = 0
#nbr_stim_dev = 0   # default is 2. number of repetited deviation
n_stim_rep=params_stim['nbr_rep_std']+params_stim['nbr_rep_dev']+2
params_sim['t_total'] = params_stim['Tresting'] + n_stim_rep*(params_stim['Tinter']+params_stim['Tstim'])   # Why is adding 4? 
params_stim['t_total']=params_sim['t_total']

ndt=round(params_sim['t_total']/params_sim['dt'])
ndt_per_trial=round((params_stim['Tinter']+params_stim['Tstim'])/params_sim['dt'])
ndt_rest= round(params_stim['Tresting']/params_sim['dt'])

N_column=8
N_eachC=40
N_neu=N_column*N_eachC
params_stim['num_perG'] = N_neu   # Number of cells in each group. 

N_column=int(round(N_column))
N_eachC=int(round(N_eachC))
N_neu=int(round(N_neu))


with open(parent_dir+'\\net\\Output.pkl', 'rb') as f:
    my_dict_loaded = pickle.load(f)  
    my_para_loaded= pickle.load(f)

Trained_W_EdS=my_dict_loaded['W_EdS'][:, my_para_loaded['n_frames']-1];
Trained_W_EP=my_dict_loaded['W_EP'][:, my_para_loaded['n_frames']-1];
stimuli_inc_base_eachC = np.linspace(
    params_stim['str_g_back']*params_stim['ratio_min'], params_stim['str_g_back']*params_stim['ratio_max'], N_eachC)
stimuli_inc_base=np.tile(stimuli_inc_base_eachC,N_column)

params_stim['N_eachC'] = N_eachC   # Number of cells in each group. 
params_stim['N_column'] = N_column   # Number of cells in each group. 
params_stim['N_neu'] = N_neu   # Number of cells in each group. 

Ind_vec_rest=np.zeros(ndt_rest)
Ind_vec_per_trial= np.arange(0,ndt_per_trial)
Ind_vec_control= np.zeros(ndt_per_trial)

Ind_vec_std_temp= np.concatenate((Ind_vec_rest,
                               np.tile(Ind_vec_per_trial,params_stim['nbr_rep_std']),
                               np.tile(Ind_vec_control,params_stim['nbr_rep_dev']),
                               np.tile(Ind_vec_per_trial,2)))
Ind_vec_dev_temp= np.concatenate((Ind_vec_rest,
                               np.tile(Ind_vec_control,params_stim['nbr_rep_std']),
                               np.tile(Ind_vec_per_trial,params_stim['nbr_rep_dev']),
                               np.tile(Ind_vec_control,2)))

fig, ax=plt.subplots(dpi=250,figsize=(2,1.7))
ax.plot(Ind_vec_std_temp)
plt.show()
fig, ax=plt.subplots(dpi=250,figsize=(2,1.7))
ax.plot(stimuli_inc_base)
plt.show()

Ind_vec_std=Ind_vec_std_temp>(params_stim['Tinter']/params_sim['dt'])
Ind_vec_dev=Ind_vec_dev_temp>(params_stim['Tinter']/params_sim['dt'])
Ind_Gate= Ind_vec_std*1+Ind_vec_dev*1;

adaptive_vector=Ind_vec_rest
adap_ratio=0.5
for i_sti in range(params_stim['nbr_rep_std']):
    append_vec= np.ones(ndt_per_trial)*adap_amount*(adap_ratio**i_sti)+(1-adap_amount)
    adaptive_vector=np.append(adaptive_vector, append_vec)
    
for i_sti in range(params_stim['nbr_rep_dev']):
    append_vec= np.ones(ndt_per_trial)*adap_amount*(adap_ratio**i_sti)+(1-adap_amount)
    adaptive_vector=np.append(adaptive_vector, append_vec)
for i_sti in range(2):
    append_vec= np.ones(ndt_per_trial)*adap_amount*(adap_ratio**i_sti)+(1-adap_amount)
    adaptive_vector=np.append(adaptive_vector, append_vec)

adaptive_vector_reshaped = adaptive_vector.reshape(1, -1)

#index_1d = j * N_col + i   for stim_2d[i,j]
# Here stim2d for debugging.

stimulus_std_inc,stim2d =bf.generate_ringblock_stim(params_stim['std_id'], N_column, N_eachC,
                                            params_stim['str_g']*params_stim['ratio_min'],
                                            params_stim['str_g']*params_stim['ratio_max'],
                                            params_stim['sigma'])
fig, ax=plt.subplots(dpi=250,figsize=(2,1.7))
ax.plot(stimulus_std_inc)
plt.show()

temp_vec=  np.linspace(1.0,1.1, N_eachC)
background_inc_ratio=np.tile(temp_vec, N_column)

stimulus_dev_inc,_=bf.generate_ringblock_stim(params_stim['dev_id'], N_column, N_eachC,
                                            omission_flag*params_stim['str_g']*params_stim['ratio_min'],
                                            omission_flag*params_stim['str_g']*params_stim['ratio_max'],
                                            params_stim['sigma'])

stimulus_std_ave_sst,_=bf.generate_ringblock_stim(params_stim['std_id'], N_column, N_eachC,
                                            params_stim['str_g']*params_stim['SST_ratio'],
                                            params_stim['str_g']*params_stim['SST_ratio'],
                                            params_stim['sigma'])

stimulus_dev_ave_sst,_=bf.generate_ringblock_stim(params_stim['dev_id'], N_column, N_eachC,
                                            omission_flag*params_stim['str_g']*params_stim['SST_ratio'],
                                            omission_flag*params_stim['str_g']*params_stim['SST_ratio'],
                                            params_stim['sigma'])


fig, ax=plt.subplots(dpi=250,figsize=(2,1.7))
ax.plot(stimulus_std_ave_sst)
plt.show()

stimuli_inc_back = np.linspace(
    params_stim['str_g_back']*params_stim['ratio_min'], params_stim['str_g_back']*params_stim['ratio_max'], N_eachC)


# Now has two different stimuli, for pyramidal and SST separately. 
stimuli_inc = np.zeros((params_stim['N_neu'] , ndt)) #full matrix of stimuli for every neurons
stimuli_ave_SST = np.zeros((params_stim['N_neu'] , ndt)) #full matrix of stimuli for every neurons

stimuli_inc[:,Ind_vec_std] = np.transpose(np.tile(stimulus_std_inc, (np.sum(Ind_vec_std),1) ))
stimuli_inc[:,Ind_vec_dev] = np.transpose(np.tile(stimulus_dev_inc, (np.sum(Ind_vec_dev),1) ))

stimuli_ave_SST[:,Ind_vec_std] = np.transpose(np.tile(stimulus_std_ave_sst, (np.sum(Ind_vec_std),1) ))
stimuli_ave_SST[:,Ind_vec_dev] = np.transpose(np.tile(stimulus_dev_ave_sst, (np.sum(Ind_vec_dev),1) ))

stimuli_ave_SST_adap=stimuli_ave_SST*adaptive_vector_reshaped*Input_S_ratio
stimuli_inc_adap=stimuli_inc*adaptive_vector_reshaped
# Test is done in main_stitest. How about let me do it here before moving on. Should I write this
# into a seperate function? Not sure yet.
# Testing the stimuli:
fig, ax =  plt.subplots(dpi=250,figsize=(4,3))
im=ax.imshow(stimuli_inc_adap,interpolation='nearest',aspect='auto',extent=[-2.5,10,N_neu-1,0]) 
cbar = ax.figure.colorbar(im, ax=ax)
ax.set_xlim(left=-1)  # This changes the view to start from -1 on the x-axis
#fig.savefig(PltTool.sdir+'AdapStimuli'+'.jpg',bbox_inches='tight')
#fig.savefig(PltTool.sdir+'AdapStimuli'+'.svg',bbox_inches='tight')
ax.set_title('Adaptive Stimuli')
plt.show()

fig, ax =  plt.subplots(dpi=250,figsize=(4,3))
im=ax.imshow(stimuli_ave_SST_adap,interpolation='nearest',aspect='auto') 
cbar = ax.figure.colorbar(im, ax=ax)
plt.show()


#### Now setting up the single cells. Using the structure set up in data_structure.py.
# Notice the STP related features are saved heretoo
E_pop= data_structure.Pyramidal(params_stim['N_neu'])  # Initial pyramidal cells data strcuture. Including rate(some), vdend
Int_pop = data_structure.Integrator(params_stim['N_column'])
P_pop= data_structure.Interneuron(params_stim['N_neu'])
S_pop= data_structure.Interneuron(params_stim['N_neu'])
P0_pop= data_structure.Interneuron(params_stim['N_neu']) # New one that generate a global control.

### Initiate the connectivity matrix. Using the create_ring_connectivity() from Basic_functions.
# This should be constrained by the Allen's data Compagnota et al. 2022. In the sense of relative width. 
# Potential modulation of existing parameters here. 



W_dict = OrderedDict()   # 10052023 Now here would be the major changes. 

W_dict['EE']= 0 # Hope this work. Double check if error returns.
W_dict['PE']= 0
#'''
W_dict['P0E']= bf.create_ringbulk_connectivity(
    params_E['weight_to_pv0'],params_E['weight_to_pv0'],
    N_col=N_column, N_eachc=N_eachC, type='BulktoBulk')
#'''


W_dict['SE']= 0

#W_dict['EE']= bf.create_ring_connectivity(N_neu,N_neu,params_E['weight_to_soma']*Wscale*WeeScale,params_E['sigma_to_pc'])  # Jmax, sigma=43.2, Jmin=0
#W_dict['PE']= bf.create_ring_connectivity(N_neu,N_neu,params_E['weight_to_pv']*Wscale,params_E['sigma_to_pv'])    # sigma=43.2
#W_dict['SE']= bf.create_ring_connectivity(N_neu,N_neu,params_E['weight_to_sst']*Wscale,params_E['sigma_to_sst'])  # Uniform excitation. I 
# I need a bulk connectivity that do not specify the source of Top area.
W_dict['IntE']=bf.create_ringbulk_connectivity(
    params_E['weight_to_integrator'],params_E['weight_to_integrator'],
    N_col=N_column, N_eachc=N_eachC, type='BulktoInt')*flag_EtoInt


#W_dict['IntE']= bf.create_ring_connectivity(N_neu,N_neu,params_E['weight_to_integrator']*Wscale*WbottopupScale,params_E['sigma_to_integrator'])  #sigma=1
'''
fig, ax =  plt.subplots(dpi=250,figsize=(2,1.7))
im=ax.imshow(W_dict['IntE'])    
cbar = ax.figure.colorbar(im, ax=ax)
plt.show()
fig, ax =  plt.subplots(dpi=250,figsize=(2,1.7))
im=ax.imshow(W_dict['PE'])    
cbar = ax.figure.colorbar(im, ax=ax)
plt.show()
'''
W_dict['EP'] = bf.load_ringbulk_connectivity(
    Trained_W_EP,N_column,type='BulktoBulk_diag')

W_dict['EP0'] = bf.create_delta_connectivity(
    params_P['weight_to_soma0'],N_neu)

W_dict['PP'] =0 
W_dict['SP'] =0

#W_dict['EP'] = bf.create_ring_connectivity(N_neu,N_neu, params_P['weight_to_soma']*Wscale*WPVInhScale,params_P['sigma_to_pc'])  # Sigma=20. Need to study later. Diagonal in Kevin code
#W_dict['PP'] = bf.create_ring_connectivity(N_neu,N_neu, params_P['weight_to_pv']*Wscale,params_P['sigma_to_pv'])  #  Sigma=20
#W_dict['SP'] = bf.create_ring_connectivity(N_neu,N_neu, params_P['weight_to_sst']*Wscale,params_P['sigma_to_sst'])

fig, ax =  plt.subplots(dpi=250,figsize=(2,1.7))
ax.plot(np.diag(W_dict['EP']))   
ax.set_title('From P to E')   
plt.show()



'''
fig, ax =  plt.subplots(dpi=250,figsize=(2,1.7))
im=ax.imshow(W_dict['PP'])  
cbar = ax.figure.colorbar(im, ax=ax)  
plt.show()
'''
W_dict['EdS']=bf.load_ringbulk_connectivity(Trained_W_EdS,N_column,type='BulktoBulk_diag')
#W_dict['EdS']= bf.create_delta_connectivity(params_S['weight_to_dend'],N_neu)*flag_StoE
#W_dict['EdS']= bf.create_ring_connectivity(N_neu,N_neu,params_S['weight_to_dend']*WSSTInhScale,params_S['sigma_to_pc'])

fig, ax =  plt.subplots(dpi=250,figsize=(2,1.7))
ax.plot(np.diag(W_dict['EdS']))   
ax.set_title('From S to Ed')   
plt.show()

W_dict['PS']= 0
#W_dict['PS']= bf.create_ring_connectivity(N_neu,N_neu,params_S['weight_to_pv'],params_S['sigma_to_pv'])

# Here I should do more systematically, so far just directly change it.
W_dict['EdInt']= bf.create_ringbulk_connectivity(
    params_Int['weight_to_dend']*params_Int['scaleMax'],params_Int['weight_to_dend']*params_Int['scaleMin'],
    N_col=N_column, N_eachc=N_eachC, type='InttoBulk')

'''
W_dict['SInt']= bf.create_ringbulk_connectivity(
    params_Int['weight_to_sst']*WsTopSMax,params_Int['weight_to_sst']*WsTopSmin,
    N_col=N_column, N_eachc=N_eachC, type='InttoBulk')*WsTopscale*WSIntscale
'''

W_dict['SInt']= np.zeros((N_neu,N_column))
# Using the trained connectivity to mimick the learned results. 

#W_dict['SInt']= bf.create_ring_connectivity(N_neu,N_neu,params_Int['weight_to_sst']*WsTopscale,params_Int['sigma_to_sst'])   # Wasn't existed.

W_dict['PInt']= bf.create_ringbulk_connectivity(
    params_Int['weight_to_pv'],params_Int['weight_to_pv'],
    N_col=N_column, N_eachc=N_eachC, type='InttoBulk')
#W_dict['PInt']= bf.create_ring_connectivity(N_neu,N_neu,params_Int['weight_to_pv']*WsTopscale,params_Int['sigma_to_pv'])   # Wasn't existed.

W_dict['IntInt']= bf.create_ring_connectivity(
    N_column,N_column,params_Int['Jmax'],params_Int['sigma_to_int'], 
    N_column= N_column, Jmin=params_Int['Jmin'])

fig, ax =  plt.subplots(dpi=250,figsize=(4,3))
im=ax.imshow(W_dict['PInt'],aspect='auto') 
#im=ax.imshow(W_dict['EdInt'],aspect='auto') 

cbar = ax.figure.colorbar(im, ax=ax)
plt.title('Connectivity int to P')
#plt.title('Connectivity Int to Ed')
plt.show()


'''
fig, ax =  plt.subplots(dpi=200,figsize=(2,1.7))
im=ax.imshow(W_dict['SInt'])   
cbar = ax.figure.colorbar(im, ax=ax)
ax.set_title('From Int to Int')
plt.show()
'''
g_back = dict(
    Es=params_E['gbg_soma'],
    Ed=params_E['gbg_dend'],
    Int=params_Int['gbg'],

    P=params_P['gbg'],
    P0=params_P['gbg'],

    S=params_S['gbg'],
    
#    Ekappa=params_E['kappa_NMDA'],  # How much excitation can get from NMDA.
#    Intkappa=params_Int['kappa_NMDA'],  # default 0.5 for both cases

    Ekappa=kappa,  # How much excitation can get from NMDA.
    Intkappa=kappa,  # default 0.5 for both cases
    IntIntkappa=1, # at top circuit, only NMDA works. 

)


### Initialize a dataframe to save desired output. 
n_frames=5000   # Number of time frames I save. 
save_interval = ndt // n_frames  # Calculate interval at which to save frames
i_frames=0
Output_dic= dict(
    E_pop_rate= np.full((N_neu,n_frames),np.nan, dtype=np.float32) ,
    E_pop_Vdend= np.full((N_neu,n_frames),np.nan, dtype=np.float32) ,
    E_pop_Isd= np.full((N_neu,n_frames),np.nan, dtype=np.float32) ,
    E_pop_Isum= np.full((N_neu,n_frames),np.nan, dtype=np.float32) ,
    P_pop_rate= np.full((N_neu,n_frames),np.nan, dtype=np.float32) ,
    P0_pop_rate= np.full((N_neu,n_frames),np.nan, dtype=np.float32) ,

    P_pop_sumI= np.full((N_neu,n_frames),np.nan, dtype=np.float32) ,
    P0_pop_sumI= np.full((N_neu,n_frames),np.nan, dtype=np.float32) ,
    Int_pop_rate= np.full((N_column,n_frames),np.nan, dtype=np.float32) ,
    Int_pop_sumI= np.full((N_column,n_frames),np.nan, dtype=np.float32) ,
    S_pop_rate= np.full((N_neu,n_frames),np.nan, dtype=np.float32) ,
    S_pop_sumI= np.full((N_neu,n_frames),np.nan, dtype=np.float32) ,
    
    g_E_sE=np.full((N_neu,n_frames),np.nan, dtype=np.float32) ,
    g_E_sI=np.full((N_neu,n_frames),np.nan, dtype=np.float32) ,
    g_E_dE=np.full((N_neu,n_frames),np.nan, dtype=np.float32) ,
    g_E_dI=np.full((N_neu,n_frames),np.nan, dtype=np.float32) ,

    g_E_sEN=np.full((N_neu,n_frames),np.nan, dtype=np.float32) ,
    g_E_dEN=np.full((N_neu,n_frames),np.nan, dtype=np.float32) ,

    g_Int=np.full((N_column,n_frames),np.nan, dtype=np.float32) ,
    g_Int_N=np.full((N_column,n_frames),np.nan, dtype=np.float32) ,

    g_P_E=np.full((N_neu,n_frames),np.nan, dtype=np.float32) ,
    g_P_EN=np.full((N_neu,n_frames),np.nan, dtype=np.float32) ,
    g_P_I=np.full((N_neu,n_frames),np.nan, dtype=np.float32) ,

    g_S_E=np.full((N_neu,n_frames),np.nan, dtype=np.float32) ,
    g_S_EN=np.full((N_neu,n_frames),np.nan, dtype=np.float32) ,
    g_S_I=np.full((N_neu,n_frames),np.nan, dtype=np.float32) ,
    
    h_E= np.full((N_neu,n_frames),np.nan, dtype=np.float32) ,
    h_EN= np.full((N_neu,n_frames),np.nan, dtype=np.float32) ,
    h_Int= np.full((N_column,n_frames),np.nan, dtype=np.float32) ,
    h_IntN= np.full((N_column,n_frames),np.nan, dtype=np.float32) ,
    )


### Initiate the updating functions, using update_functions.py
RHS_update= update_functions.RHS_update(W_dict, g_back)
#RHS_update.g_back['Es'] +=stimuli_inc_back

LHS_update_E= update_functions.LHS_update.Pyr(params_E,params_sim['dt'])
LHS_update_Int= update_functions.LHS_update.Int(params_Int,params_sim['dt'])

LHS_update_P= update_functions.LHS_update.P(params_P,params_sim['dt'])
LHS_update_P0= update_functions.LHS_update.P(params_P,params_sim['dt'])

LHS_update_S= update_functions.LHS_update.S(params_S,params_sim['dt'])


### main loop. Noticing currently, I don't include NMDA, STP, LTP, Adaptation.
#t_PVadvance=int(0.05//params_sim['dt']) 

#Int_hgate=0.01
#Int_hNgate=0.16128

for i_time in tqdm(range(ndt)):
#for i_time in tqdm(range(int(10/params_sim['dt']))):    

    if i_time == ndt_rest-1:
        Int_hgate=Int_pop.h
        Int_hNgate=Int_pop.hN
        Int_rategate=Int_pop.rate

    # In each timestep, rate, conductance, and current are updated.
    # tau_I=0, so it's calculated as RHS. 
    # conductance (pre-synaptic) depends on the firing rate of the same cell. Update second
    # rate, doesn't depend on conductance. So update third.
    
    # Input
    RHS_update.g_back['Es'] = stimuli_inc_adap[:,i_time]+g_back['Es']+stimuli_inc_base
    #RHS_update.g_back['Ed'] = g_back['Ed']
    RHS_update.g_back['S'] = stimuli_ave_SST_adap[:,i_time]+g_back['S']

    # Gate
    if Ind_Gate[i_time]>0.5 or i_time<ndt_rest:
        Int_pop_h_iter=Int_pop.h
        Int_pop_hN_iter=Int_pop.hN
    else:   # Calculate with Int_pop.rate=5 Can change to the baseline later
        Int_pop_h_iter= Int_hgate            #Int_pop.h*0+0.01    
        Int_pop_hN_iter= Int_hNgate          #Int_pop.hN*0+0.16128

    
    # RHS update. Calculate the post-synaptic conductance based on connectivity.
    E_pop.gE_soma= RHS_update.g_pyr_somaE(E_pop.h)
    E_pop.gEN_soma= RHS_update.g_pyr_somaEN(E_pop.h)
    E_pop.gI_soma= RHS_update.g_pyr_somaI(P_pop.h,P0_pop.h)


    E_pop.gE_dend= RHS_update.g_pyr_dendE(Int_pop_h_iter)
    E_pop.gEN_dend= RHS_update.g_pyr_dendEN(Int_pop_hN_iter)
    E_pop.gI_dend= RHS_update.g_pyr_dendI(S_pop.h)    
    
    Int_pop.g= RHS_update.g_Int(E_pop.h, Int_pop.h)
    Int_pop.gN= RHS_update.g_IntN(E_pop.hN, Int_pop.hN)
    #Int_pop.input=stimuli[:,i_time]
    '''
    Int_pop.g= RHS_update.g_Int(E_pop.h, Int_pop.h)
    Int_pop.g_Iglobal=RHS_update.g_Int_Iglobal(Int_pop.h)   # Added 0913. ??? Should be gaba
    '''    
    P_pop.gE= RHS_update.g_pv_E(E_pop.h, Int_pop_h_iter)    
    P_pop.gEN= RHS_update.g_pv_EN(E_pop.hN, Int_pop_hN_iter)    
    P_pop.gI= RHS_update.g_pv_I(P_pop.h, S_pop.h) 
    
    P0_pop.gE= RHS_update.g_pv0_E(E_pop.h)    
    P0_pop.gEN= RHS_update.g_pv0_EN(E_pop.hN)    
    P0_pop.gI= 0   #RHS_update.g_pv_I(P_pop.h, S_pop.h) 
    
    S_pop.gE= RHS_update.g_sst_E(E_pop.h, Int_pop_h_iter)   
    S_pop.gEN= RHS_update.g_sst_EN(E_pop.hN, Int_pop_hN_iter)    
    S_pop.gI= RHS_update.g_sst_I(P_pop.h) 

    # LHS update all the conductance.
    E_pop.h= LHS_update_E.h(E_pop.h,E_pop.rate)
    E_pop.hN= LHS_update_E.hNMDA(E_pop.h,E_pop.rate)
    
    
    Int_pop.h= LHS_update_Int.h(Int_pop.h,Int_pop.rate)
    Int_pop.hN= LHS_update_Int.hNMDA(Int_pop.hN,Int_pop.rate)

    P_pop.h= LHS_update_P.h(P_pop.h,P_pop.rate)
    P0_pop.h= LHS_update_P0.h(P0_pop.h,P0_pop.rate)

    S_pop.h= LHS_update_S.h(S_pop.h,S_pop.rate)

    # LHS update all the firing rates.
    E_pop.rate, E_pop.Vdend, Isd, Isum= LHS_update_E.rate(E_pop )
    P_pop.rate, p_total_current=LHS_update_P.rate(P_pop )
    P0_pop.rate, p0_total_current=LHS_update_P.rate(P0_pop )
    
    Int_pop.rate, Int_total_current=LHS_update_Int.rate(Int_pop )
    S_pop.rate, S_total_current=LHS_update_S.rate(S_pop )

    
    # saving. 
    if i_time % save_interval == 0 and i_frames < n_frames:
        # Replace 'your_data' with the actual data you want to save
        
        Output_dic['E_pop_rate'][:,i_frames]=E_pop.rate
        Output_dic['E_pop_Vdend'][:,i_frames]=E_pop.Vdend
        Output_dic['E_pop_Isd'][:,i_frames] = Isd
        Output_dic['E_pop_Isum'][:,i_frames] = Isum
    
        Output_dic['P_pop_rate'][:,i_frames]=P_pop.rate
        Output_dic['P_pop_sumI'][:,i_frames]=p_total_current
        Output_dic['P0_pop_rate'][:,i_frames]=P0_pop.rate
        Output_dic['P0_pop_sumI'][:,i_frames]=p_total_current

        Output_dic['Int_pop_rate'][:,i_frames]=Int_pop.rate
        Output_dic['Int_pop_sumI'][:,i_frames]=Int_total_current
        
        Output_dic['S_pop_rate'][:,i_frames]=S_pop.rate
        Output_dic['S_pop_sumI'][:,i_frames]=S_total_current
    
        Output_dic['g_E_sE'][:,i_frames]=E_pop.gE_soma
        Output_dic['g_E_sI'][:,i_frames]=E_pop.gI_soma
        Output_dic['g_E_dE'][:,i_frames]=E_pop.gE_dend
        Output_dic['g_E_dI'][:,i_frames]=E_pop.gI_dend
    
        Output_dic['g_E_sEN'][:,i_frames]=E_pop.gEN_soma
        Output_dic['g_E_dEN'][:,i_frames]=E_pop.gEN_dend
    
        Output_dic['g_Int'][:,i_frames]=Int_pop.g
        Output_dic['g_Int_N'][:,i_frames]=Int_pop.gN

    
        Output_dic['g_P_E'][:,i_frames]=P_pop.gE
        Output_dic['g_P_EN'][:,i_frames]=P_pop.gEN
        Output_dic['g_P_I'][:,i_frames]=P_pop.gI
        Output_dic['g_S_E'][:,i_frames]=S_pop.gE
        Output_dic['g_S_EN'][:,i_frames]=S_pop.gEN
        Output_dic['g_S_I'][:,i_frames]=S_pop.gI
        
        Output_dic['h_E'][:,i_frames]=E_pop.h
        Output_dic['h_EN'][:,i_frames]=E_pop.hN
        Output_dic['h_Int'][:,i_frames]=Int_pop.h
        Output_dic['h_IntN'][:,i_frames]=Int_pop.hN

        i_frames+=1
    

### output, maybe later analyze in a separate code. 
checkind_0=20                         
t_axis=np.arange(ndt)*params_sim['dt']


'''
xyscale=ndt/100/2

fig, axes=plt.subplots(2,1,figsize=(4,3))

# t_axis
im=axes[0].imshow(Output_dic['E_pop_rate'][:,int(1/params_sim['dt']):],aspect=xyscale)
axes[0].set_title('Ending E rate'+ f"{Output_dic['E_pop_rate'][checkind_0,-1]:.4f}")
cbar = ax.figure.colorbar(im, ax=axes[0])
im=axes[1].imshow(Output_dic['E_pop_Vdend'][:,int(1/params_sim['dt']):],aspect=xyscale)
axes[1].set_title('Ending E dendV'+ f"{Output_dic['E_pop_Vdend'][checkind_0,-1]:.4f}")
fig.subplots_adjust(hspace=0.5)
cbar = ax.figure.colorbar(im, ax=axes[1])

plt.show()

fig, axes=plt.subplots(2,1,figsize=(4,3))

# t_axis
im=axes[0].imshow(Output_dic['P_pop_rate'][:,int(1/params_sim['dt']):],aspect=xyscale)
axes[0].set_title('Ending P rate'+ f"{Output_dic['P_pop_rate'][checkind_0,-1]:.4f}")
cbar = ax.figure.colorbar(im, ax=axes[0])
im=axes[1].imshow(Output_dic['P_pop_sumI'][:,int(1/params_sim['dt']):],aspect=xyscale)
axes[1].set_title('Ending PsumI'+ f"{Output_dic['P_pop_sumI'][checkind_0,-1]:.4f}")
cbar = ax.figure.colorbar(im, ax=axes[1])

fig.subplots_adjust(hspace=0.5)
plt.show()

fig, axes=plt.subplots(2,1,figsize=(4,3))

# t_axis
im=axes[0].imshow(Output_dic['Int_pop_rate'][:,int(1/params_sim['dt']):],aspect=xyscale)
axes[0].set_title('Ending rate Int'+ f"{Output_dic['Int_pop_rate'][checkind_0,-1]:.4f}")
cbar = ax.figure.colorbar(im, ax=axes[0])


im=axes[1].imshow(Output_dic['S_pop_rate'][:,int(1/params_sim['dt']):],aspect=xyscale)
axes[1].set_title('Ending S rate'+ f"{Output_dic['S_pop_rate'][checkind_0,-1]:.4f}")
cbar = ax.figure.colorbar(im, ax=axes[1])

fig.subplots_adjust(hspace=0.5)
plt.show()
'''
''' #A few variables to check
Output_dic['E_pop_rate'][20,:100]
Output_dic['E_pop_Vdend'][20,:100]
Output_dic['P_pop_rate'][20,:100]
Output_dic['P_pop_sumI'][20,:100]

Output_dic['S_pop_rate'][20,:100]
Output_dic['Int_pop_rate'][20,:100]

Output_dic['g_E_sE'][20,:100]
Output_dic['g_E_sI'][20,:100]
Output_dic['g_E_dE'][20,:100]
Output_dic['g_E_dI'][20,:100]

Output_dic['g_Int'][20,:100]
Output_dic['g_P_E'][20,:100]
Output_dic['g_P_I'][20,:100]
Output_dic['g_S_E'][20,:100]
Output_dic['g_S_I'][20,:100]

'''


import AnalysisTool as AT
PltTool=AT.Tools(Output_dic, ndt,n_frames, params_sim['dt'], params_stim, 
                 Plot_flag=True, omission_flag=True, save_dir=FigPath+'\\')  # Can add a output directory at the end.
PltTool.PlotPyrPop()
PltTool.PlotPVPop()
PltTool.PlotSSTPop()
PltTool.PlotIntPop()

PltTool.PlotPyrTimes()
PltTool.PlotPyrTimes_all()

PltTool.PlotPyrDDInd()
PltTool.PlotPyrDDInd_all()

PltTool.PlotPyrTimes_sample(shift_id=-N_eachC//2)
PltTool.PlotPyrTimes_sample(shift_id=N_eachC//2-1)
MMNshift=PltTool.PlotPyr(shift_id=-N_eachC//2)
MMNshift=PltTool.PlotPyr(shift_id=N_eachC//2-1)

PltTool.PlotPyrTimes_sample(shift_id=N_neu//2-N_eachC//2)
PltTool.PlotPyrTimes_sample(shift_id=N_neu//2+N_eachC//2-1)
MMNshift=PltTool.PlotPyr(shift_id=N_neu//2-N_eachC//2)
MMNshift=PltTool.PlotPyr(shift_id=N_neu//2+N_eachC//2-1)



PltTool.PlotPyrTimes_pop()
PltTool.PlotPyrTimes_pop_all()


PltTool.PlotInhPopITimes()
PltTool.PlotPyrTimesDDSSA()
PltTool.PlotPyrTimesDDSSA_all()




MMN=PltTool.PlotPyr()


#MMNshift=PltTool.PlotPyr(shift_id=-N_eachC//4)
#MMNshift=PltTool.PlotPyr(shift_id=N_eachC//2-1)

MMNInt=PltTool.PlotInt() 
MMNSST=PltTool.PlotSST()
MMNP=PltTool.PlotPV()
PltTool.PlotPyr_current(shift_id=-N_eachC//2)
#PltTool.PlotPyr_current()
#PltTool.PlotPyr_current(shift_id=N_eachC//2-1)


PltTool.PlotInt_current()  
PltTool.PlotP_current(shift_id=-N_eachC//2)
PltTool.PlotS_current(shift_id=-N_eachC//4)
PltTool.PlotS_current(shift_id=N_eachC//4)

#PltTool.Plot_Dict(Output_dic['g_Int_N'][1,:])
'''
#PltTool.std
#PltTool.dev
PltTool.Plot_Times_Dict(Output_dic['S_pop_rate'],name='SSTrate') 
PltTool.Plot_Times_Dict(Output_dic['P_pop_rate'],name='PVrate') 


#PltTool.Plot_Snapshot(Output_dic['E_pop_rate'][:,n_frames-1],'E, last frame')

PltTool.Plot_Snapshot(PltTool.FirstInd+5,tag='First')  # Which frame and tag
PltTool.Plot_Snapshot(PltTool.LastInd+5,tag='Last')  # Which frame and tag
PltTool.Plot_Snapshot(PltTool.OddInd+5,tag='Odd')  # Which frame and tag

PltTool.Plot_Snapshot(PltTool.FirstInd-5,tag='Baseline')  # Which frame and tag
'''


np.mean(PltTool.First_rate_vec)
np.mean(PltTool.Odd_rate_vec)
np.mean(PltTool.baseline_rate_vec)
# The following is calculate how many neurons show significant omission response
res_vec= PltTool.Odd_rate_vec-PltTool.baseline_rate_vec*(2)

n_omiss=np.sum(res_vec > 0)



'''
# Random control test.
with open(parent_dir+'\\net\\Output_rand\\Output0.pkl', 'rb') as f:
    dict_rand = pickle.load(f)  
    para_rand= pickle.load(f)

PltTool.PlotPyrTimes_randomcontrol(dict_rand)
PltTool.PlotPyrTimesDif_randomcontrol(dict_rand)

PltTool.Plot_Dict(dict_rand['Int_pop_rate'][6,:],name='IntRandom')
'''

#Output_dic['g_Int_N'][1,:]
#plt.close('all')
sys.path.remove(lib_dir)
