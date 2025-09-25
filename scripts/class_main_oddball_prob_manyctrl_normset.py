# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 14:27:37 2024
Modified based on previous oddball paradigm code. This one is specific to run the protocol 
with stimuli determined by a probability.
Please notice lots of analysis are different comparing to older ones. 
@author: John Meng
"""
import numpy as np

class oddball_model:
    def main(self, dev_degree=270, Tinter=0.5, Tstim=0.5, adap_amount=0.2,
             omission_flag=False, saving_tag='oddball', Stim_g=None,
             g_back_int=53,Tau_int=None, g_IntE=None, nbr_rep_std=6,nbr_rep_dev=2,
             Trained_Pathway='\\net\\Output.pkl', flag_homeo=True,Input_blocker=None,
             prob_flag=True):
        # Certain control parameters. Should set up better later.
        # Input blocker indicate starting which trial, the input is 0
        flag_StoE=1  # Inhibition control from SST to pyramidal cells.
        flag_EtoInt=1  # Bottom to top control 
        #adap_amount=0.2
        kappa=0.3             # default 0.5. Requires a fine tuning.
        
        Input_S_ratio=1
        self.omission_flag=omission_flag
        if omission_flag:
            omission_value=0.0
        else:
            omission_value=1.0
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
        if Input_blocker is not None:
            FigPath=os.path.join(parent_dir+'\\figs\\blocker\\'+saving_tag)
            if os.path.isdir(FigPath)==False:
                os.makedirs(FigPath)
        else:
            FigPath=os.path.join(parent_dir+'\\figs\\'+saving_tag)
            if os.path.isdir(FigPath)==False:
                os.makedirs(FigPath)

        filename = os.path.basename(__file__)
        
        from shutil import copyfile      # Save a copy of the code to the saving folder.
        copyfile(filename,TempPath+'\\'+filename)
        para_loc='..\lib\\'
        
        copyfile(para_loc+'parameters_tuningcurve.py',TempPath+'\parameters_tuningcurve.py')
        copyfile(filename,TempPath+'\\'+filename)
        
        import basic_functions as bf
        import parameters_tuningcurve as params
        import data_structure
        import update_functions
        import numpy as np
        from collections import OrderedDict
        import matplotlib.pyplot as plt
        import pickle
        np.random.seed(42)

        # All right, I really hate this. The work is not that transparaent. But it is one way
        PARAS=params.PARAMS_ALL.copy()
        PARAS['PARAMS_Pyr']['kappa_NMDA']=0
        PARAS['PARAMS_INT']['kappa_NMDA']=0
        params_E = PARAS['PARAMS_Pyr'].copy()
        params_Int = PARAS['PARAMS_INT'].copy()
        if Tau_int is not None:
            params_Int['tau']=Tau_int
            
            
        params_P = PARAS['PARAMS_PV'].copy()
        params_S = PARAS['PARAMS_SST'].copy()
        
        # list_param= np.array([-0.01,2,0.5,8,50,1,1,1,1,0.1,2])   
        
        # adaptation, top-down feedback, Tinter, dev_iud, delay, 
        # ndf_plasticity, int_plasticity. bool_sigma_ndf, bool_sigma_int, adaptation_timeconstant
        # nbr_rep. This description is based on the dict_param in Script7
        
        params_stim=PARAS['Stimulus'].copy()
        params_stim['prob_std'] = 0.8
        
        params_stim['dev_id']= dev_degree
        params_stim['Tinter']= Tinter
        params_stim['Tstim']= Tstim
        if Stim_g is not None:
            params_stim['str_g']=Stim_g
        
        params_sim=PARAS['Simulation'].copy()
            
        # The following par values are from the 1st parameter set in script7.    
        delay=50*0.001   # The unit should be second
        params_stim['nbr_rep_std']=nbr_rep_std
        params_stim['nbr_rep_dev']=nbr_rep_dev
        std_prob= nbr_rep_std/(nbr_rep_std+nbr_rep_dev)
        params_stim['std_prob']= std_prob      # Only specific analysis use this. Required.
        params_stim['Tresting']=3      # John, add here for longer tau.
        # Generate a stimulus. This may change into a class.
        
        #nbr_stim_std = 0
        #nbr_stim_dev = 0   # default is 2. number of repetited deviation
        n_stim_rep=params_stim['nbr_rep_std']+params_stim['nbr_rep_dev']
        params_sim['t_total'] = params_stim['Tresting'] + n_stim_rep*(params_stim['Tinter']+params_stim['Tstim'])   

        if prob_flag is False:
            n_stim_rep=12
            params_sim['t_total'] = params_stim['Tresting'] + 12*(params_stim['Tinter']+params_stim['Tstim'])
        ndt=round(params_sim['t_total']/params_sim['dt'])
        ndt_per_trial=round((params_stim['Tinter']+params_stim['Tstim'])/params_sim['dt'])
        ndt_rest= round(params_stim['Tresting']/params_sim['dt'])
        params_stim['t_total']= params_sim['t_total']

        if Input_blocker is not None:
            time_block=  params_stim['Tresting'] + (Input_blocker+1)*(params_stim['Tinter']+params_stim['Tstim']) 
            n_block=  round(time_block/params_sim['dt'])
            print(ndt)
            print(n_block)
        else:
            n_block=ndt+1

        with open(parent_dir+Trained_Pathway, 'rb') as f:
            my_dict_loaded = pickle.load(f)  
            my_para_loaded= pickle.load(f)

        
        N_column=8
        N_eachC=my_para_loaded['N_neu']
        N_neu=N_column*N_eachC
        params_stim['num_perG'] = N_neu   # Number of cells in each group. 
        
        N_column=int(round(N_column))
        N_eachC=int(round(N_eachC))
        N_neu=int(round(N_neu))
        
        

        
        Trained_W_EdS=my_dict_loaded['W_EdS'][:, my_para_loaded['n_frames']-1];
        Trained_W_EP=my_dict_loaded['W_EP'][:, my_para_loaded['n_frames']-1];


        stimuli_inc_back = np.linspace(
            params_stim['str_g_back']*params_stim['ratio_min'], params_stim['str_g_back']*params_stim['ratio_max'], N_eachC)

        np.random.seed(5)
        if not flag_homeo:
            shuffle_input_indices = np.random.permutation(N_eachC)
            stimuli_inc_back=stimuli_inc_back[shuffle_input_indices]
        
        stimuli_inc_base=np.tile(stimuli_inc_back,N_column)
        
        params_stim['N_eachC'] = N_eachC   # Number of cells in each group. 
        params_stim['N_column'] = N_column   # Number of cells in each group. 
        params_stim['N_neu'] = N_neu   # Number of cells in each group. 
        
        Ind_vec_rest=np.zeros(ndt_rest)
        Ind_vec_per_trial= np.arange(0,ndt_per_trial)
        Ind_vec_control= np.zeros(ndt_per_trial)
        Ind_vec_all=np.zeros(ndt)

        # n_stim_rep=params_stim['nbr_rep_std']+params_stim['nbr_rep_dev']

        # Randomly choose between 8 different angels while first choise 1
        if prob_flag is True:
            #choices = np.random.choice([1, 2], size=n_stim_rep, p=[std_prob, 1-std_prob])
            choices = np.random.choice(range(8), size=n_stim_rep)
            choices[0] = 2   # Guarentee the first one is standard to compare.
        else:
            choices = ([1] * params_stim['nbr_rep_std'] + [2] * params_stim['nbr_rep_dev'] +
                       [1] * (12-1-params_stim['nbr_rep_std'])        ) 
        params_stim['choices']=choices
        
        '''
        # Build the resulting vector by sampling from vect1 or vect2 based on the choices
        Ind_vec_std_2d = np.array([Ind_vec_per_trial if c == 1 else Ind_vec_control for c in choices])
        Ind_vec_std_temp = Ind_vec_std_2d.flatten()  # Or use Ind_vec_std_temp.ravel()
        Ind_vec_std_temp =  np.concatenate((Ind_vec_rest,Ind_vec_std_temp))

        Ind_vec_dev_2d = np.array([Ind_vec_control if c == 1 else Ind_vec_per_trial for c in choices])
        Ind_vec_dev_temp = Ind_vec_dev_2d.flatten()  # Or use Ind_vec_std_temp.ravel()
        Ind_vec_dev_temp =  np.concatenate((Ind_vec_rest,Ind_vec_dev_temp))
        
        # I will leave the gating part here 
        Ind_vec_std=Ind_vec_std_temp>=(params_stim['Tinter']/params_sim['dt'])
        Ind_vec_dev=Ind_vec_dev_temp>=(params_stim['Tinter']/params_sim['dt'])
        Ind_Gate= Ind_vec_std*1+Ind_vec_dev*1;
        '''
        # Rewrite the gating part 01082024
        Ind_vec_2d = np.array([Ind_vec_per_trial if c is not None else Ind_vec_per_trial for c in choices])
        Ind_vec_temp = Ind_vec_2d.flatten()  # Or use Ind_vec_std_temp.ravel()
        Ind_vec_temp =  np.concatenate((Ind_vec_rest,Ind_vec_temp))
        Ind_vec=Ind_vec_temp>=(params_stim['Tinter']/params_sim['dt'])
        Ind_Gate= Ind_vec

        '''
        Ind_vec_std_temp= np.concatenate((Ind_vec_rest,
                                       np.tile(Ind_vec_per_trial,params_stim['nbr_rep_std']),
                                       np.tile(Ind_vec_control,params_stim['nbr_rep_dev']),
                                       np.tile(Ind_vec_per_trial,2)))
        Ind_vec_dev_temp= np.concatenate((Ind_vec_rest,
                                       np.tile(Ind_vec_control,params_stim['nbr_rep_std']),
                                       np.tile(Ind_vec_per_trial,params_stim['nbr_rep_dev']),
                                       np.tile(Ind_vec_control,2)))
        '''
        fig, ax=plt.subplots(dpi=250,figsize=(2,1.7))
        ax.plot(Ind_Gate)
        plt.show()
        fig, ax=plt.subplots(dpi=250,figsize=(2,1.7))
        ax.plot(stimuli_inc_base)
        plt.show()
        
        '''
        adaptive_vector=Ind_vec_rest
        adap_ratio=0.5
        for i_sti in range(params_stim['nbr_rep_std']):
            append_vec= np.ones(ndt_per_trial)*adap_amount*(adap_ratio**i_sti)+(1-adap_amount)
            if (Input_blocker is not None) and (i_sti>=Input_blocker):
                append_vec=0*append_vec
            adaptive_vector=np.append(adaptive_vector, append_vec)
            
        for i_sti in range(params_stim['nbr_rep_dev']):
            append_vec= np.ones(ndt_per_trial)*adap_amount*(adap_ratio**i_sti)+(1-adap_amount)
            #if Input_blocker is not None:
            #    append_vec=0*append_vec
            adaptive_vector=np.append(adaptive_vector, append_vec)
        print(Input_blocker)    
        if prob_flag is False:
            for i_sti in range(12-params_stim['nbr_rep_dev']-params_stim['nbr_rep_std']):
                append_vec= np.ones(ndt_per_trial)*adap_amount*(adap_ratio**i_sti)+(1-adap_amount)
                if Input_blocker is not None:
                    append_vec=0*append_vec
                    print(f'i_sti{i_sti}')
                adaptive_vector=np.append(adaptive_vector, append_vec)
        
        
        adaptive_vector_reshaped = adaptive_vector.reshape(1, -1)
        '''
        
        #index_1d = j * N_col + i   for stim_2d[i,j]
        # Here stim2d for debugging.
        degree_vec=np.array([0,45,90,135,180,225,270,315])

        stimuli_inc_dic=[]
        stimuli_ave_sst_dic=[]
        for (i_deg,degree) in enumerate(degree_vec):
            stim_inc_ite,_=bf.generate_ringblock_stim(degree, N_column, N_eachC,
                                                        params_stim['str_g']*params_stim['ratio_min'],
                                                        params_stim['str_g']*params_stim['ratio_max'],
                                                        params_stim['sigma'])
            
            stim_ave_sst_ite,_=bf.generate_ringblock_stim(degree, N_column, N_eachC,
                                                        params_stim['str_g']*params_stim['SST_ratio'],
                                                        params_stim['str_g']*params_stim['SST_ratio'],
                                                        params_stim['sigma'])
            stimuli_inc_dic.append(stim_inc_ite)
            stimuli_ave_sst_dic.append(stim_ave_sst_ite)
            
        # sanity test     Maybe there are indexing issue
        fig, ax=plt.subplots(dpi=250,figsize=(2,1.7))
        ax.plot(stimuli_inc_dic[2])
        plt.show()
       
        fig, ax=plt.subplots(dpi=250,figsize=(2,1.7))
        ax.plot(stimuli_ave_sst_dic[2])
        plt.show()
        

        # Now has two different stimuli, for pyramidal and SST separately. 
        stimuli_inc = np.zeros((params_stim['N_neu'] , ndt)) #full matrix of stimuli for every neurons
        stimuli_ave_SST = np.zeros((params_stim['N_neu'] , ndt)) #full matrix of stimuli for every neurons
        

        # Updating here: NOT FINISHED
        for (i_rep,i_choice) in enumerate(choices):
            ind_start=ndt_rest+ ndt_per_trial*i_rep
            #print(f'i_rep{i_rep},ind_start{ind_start} ')
            stimulus_inc_rep= stimuli_inc_dic[i_choice]
            stimulus_dec_SST_rep=stimuli_ave_sst_dic[i_choice]
            
            Ind_vec_this_rep= Ind_vec_all.copy()
            Ind_vec_this_rep[ind_start:(ind_start+ndt_per_trial)]=Ind_vec_per_trial
            flag_this_rep= Ind_vec_this_rep>(params_stim['Tinter']/params_sim['dt'])
            
            stimuli_inc[:,flag_this_rep] = np.transpose(np.tile(stimulus_inc_rep, (np.sum(flag_this_rep),1) ))
            stimuli_ave_SST[:,flag_this_rep] = np.transpose(np.tile(stimulus_dec_SST_rep, (np.sum(flag_this_rep),1) ))
            
                
        
      
        #stimuli_ave_SST_adap=stimuli_ave_SST*adaptive_vector_reshaped*Input_S_ratio
        #stimuli_inc_adap=stimuli_inc*adaptive_vector_reshaped
        stimuli_ave_SST_adap=stimuli_ave_SST*Input_S_ratio
        stimuli_inc_adap=stimuli_inc
        
        
        
        # Test is done in main_stitest. How about let me do it here before moving on. Should I write this
        # into a seperate function? Not sure yet.
        # Testing the stimuli:
        fig, ax =  plt.subplots(dpi=250,figsize=(4,3))
        im=ax.imshow(stimuli_inc,interpolation='nearest',aspect='auto',extent=[-2.5,10,N_neu-1,0]) 
        cbar = ax.figure.colorbar(im, ax=ax)
        ax.set_xlim(left=-1)  # This changes the view to start from -1 on the x-axis
        #fig.savefig(PltTool.sdir+'AdapStimuli'+'.jpg',bbox_inches='tight')
        #fig.savefig(PltTool.sdir+'AdapStimuli'+'.svg',bbox_inches='tight')
        ax.set_title('Raw Input vector')
        plt.show()    
            
 
            
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
        
        ### Initiate the connectivity matrix. Using the create_ring_connectivity() from Basic_functions.
        # This should be constrained by the Allen's data Compagnota et al. 2022. In the sense of relative width. 
        # Potential modulation of existing parameters here. 
        
        
        
        W_dict = OrderedDict()   # 10052023 Now here would be the major changes. 
        
        W_dict['EE']= 0 # Hope this work. Double check if error returns.
        W_dict['PE']= 0

        
        
        W_dict['SE']= 0
        
        #W_dict['EE']= bf.create_ring_connectivity(N_neu,N_neu,params_E['weight_to_soma']*Wscale*WeeScale,params_E['sigma_to_pc'])  # Jmax, sigma=43.2, Jmin=0
        #W_dict['PE']= bf.create_ring_connectivity(N_neu,N_neu,params_E['weight_to_pv']*Wscale,params_E['sigma_to_pv'])    # sigma=43.2
        #W_dict['SE']= bf.create_ring_connectivity(N_neu,N_neu,params_E['weight_to_sst']*Wscale,params_E['sigma_to_sst'])  # Uniform excitation. I 
        # I need a bulk connectivity that do not specify the source of Top area.
        if g_IntE is not None:
            params_E['weight_to_integrator']=g_IntE   # 6 by default
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
        if not flag_homeo:
            W_dict['EdInt']= bf.create_ringbulk_connectivity(
            params_Int['weight_to_dend']*params_Int['scaleMax'],params_Int['weight_to_dend']*params_Int['scaleMin'],
            N_col=N_column, N_eachc=N_eachC, type='InttoBulk',shuffle_inds=shuffle_input_indices)
        
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
        im=ax.imshow(W_dict['IntInt'],aspect='auto') 
        #im=ax.imshow(W_dict['EdInt'],aspect='auto') 
        
        cbar = ax.figure.colorbar(im, ax=ax)
        plt.title('Integrator Connectivity')
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
        
            S=params_S['gbg'],
            
        #    Ekappa=params_E['kappa_NMDA'],  # How much excitation can get from NMDA.
        #    Intkappa=params_Int['kappa_NMDA'],  # default 0.5 for both cases
        
            Ekappa=kappa,  # How much excitation can get from NMDA.
            Intkappa=kappa,  # default 0.5 for both cases
            IntIntkappa=1, # at top circuit, only NMDA works. 
        
        )
        if g_back_int is not None:
            g_back['Int']=g_back_int
        
        ### Initialize a dataframe to save desired output. This varies due to the round-up problem in the vaires interval case. Let's change it.
        save_interval_time= 0.025
        save_interval=round(save_interval_time /params_sim['dt'])

        #n_frames= round(ndt//save_interval)
        n_frames= round(params_sim['t_total']/save_interval_time)
        temp_t_total=params_sim['t_total']

        '''
        ndt=5000   # Number of time frames I save. Uf totak
        save_interval = ndt // n_frames  # Calculate interval at which to save frames
        '''
        i_frames=0
        Output_dic= dict(
            E_pop_rate= np.full((N_neu,n_frames),np.nan, dtype=np.float32) ,
            E_pop_Vdend= np.full((N_neu,n_frames),np.nan, dtype=np.float32) ,
            E_pop_Isd= np.full((N_neu,n_frames),np.nan, dtype=np.float32) ,
            E_pop_Isum= np.full((N_neu,n_frames),np.nan, dtype=np.float32) ,
            P_pop_rate= np.full((N_neu,n_frames),np.nan, dtype=np.float32) ,
        
            P_pop_sumI= np.full((N_neu,n_frames),np.nan, dtype=np.float32) ,
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
            if (Ind_Gate[i_time]>0.5 or i_time<ndt_rest) and i_time<n_block:   # Last one is for blocking.
                Int_pop_h_iter=Int_pop.h
                Int_pop_hN_iter=Int_pop.hN
            else:   # Calculate with Int_pop.rate=5 Can change to the baseline later
                Int_pop_h_iter= Int_hgate            #Int_pop.h*0+0.01    
                Int_pop_hN_iter= Int_hNgate          #Int_pop.hN*0+0.16128
        
            
            # RHS update. Calculate the post-synaptic conductance based on connectivity.
            E_pop.gE_soma= RHS_update.g_pyr_somaE(E_pop.h)
            E_pop.gEN_soma= RHS_update.g_pyr_somaEN(E_pop.h)
            E_pop.gI_soma= RHS_update.g_pyr_somaI(P_pop.h)
        
        
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
        
            S_pop.gE= RHS_update.g_sst_E(E_pop.h, Int_pop_h_iter)   
            S_pop.gEN= RHS_update.g_sst_EN(E_pop.hN, Int_pop_hN_iter)    
            S_pop.gI= RHS_update.g_sst_I(P_pop.h) 
        
            # LHS update all the conductance.
            E_pop.h= LHS_update_E.h(E_pop.h,E_pop.rate)
            E_pop.hN= LHS_update_E.hNMDA(E_pop.h,E_pop.rate)
            
            
            Int_pop.h= LHS_update_Int.h(Int_pop.h,Int_pop.rate)
            Int_pop.hN= LHS_update_Int.hNMDA(Int_pop.hN,Int_pop.rate)
        
            P_pop.h= LHS_update_P.h(P_pop.h,P_pop.rate)
        
            S_pop.h= LHS_update_S.h(S_pop.h,S_pop.rate)
        
            # LHS update all the firing rates.
            E_pop.rate, E_pop.Vdend, Isd, Isum= LHS_update_E.rate(E_pop )
            P_pop.rate, p_total_current=LHS_update_P.rate(P_pop )
            
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
        # adding the follwoing to be used in the analyze function
        self.Output_dic=Output_dic
        self.params_sim=params_sim
        self.params_stim=params_stim
        self.ndt=ndt
        self.n_frames=n_frames
        self.saving_tag=saving_tag
        self.N_eachC=N_eachC
        self.N_column=N_column
        self.Input_blocker=Input_blocker
        return Output_dic
    def analyze(self,Plot_flag=True): 
        import AnalysisTool as AT
        if self.Input_blocker is None:
            save_dir= ('../figs/'+self.saving_tag+'/')
        else:
            save_dir= ('../figs/blocker/'+self.saving_tag+'/')        
        PltTool=AT.Tools(self.Output_dic, self.ndt,self.n_frames, self.params_sim['dt'], self.params_stim, Plot_flag=Plot_flag,omission_flag=self.omission_flag, save_dir=save_dir)  # Can add a output directory at the end.
        PltTool.PlotPyrPop(left=-1)
        PltTool.PlotIntPop()
    
        PltTool.PlotPyrTimes_all()
        PltTool.PlotPyrDDInd_all()
        PltTool.PlotPyrTimes_sample(shift_id=-self.N_eachC//2)
        PltTool.PlotPyrTimes_sample(shift_id=self.N_eachC//2-1)
        MMNshift=PltTool.PlotPyr(shift_id=-self.N_eachC//2)
        MMNshift=PltTool.PlotPyr(shift_id=self.N_eachC//2-1)
        MMN=PltTool.PlotPyr()
        
        PltTool.PlotPyrTimes_pop_all()

    
   
        #MMNshift=PltTool.PlotPyr(shift_id=-N_eachC//4)
        #MMNshift=PltTool.PlotPyr(shift_id=N_eachC//2-1)
    
        MMNInt=PltTool.PlotInt() 

        dtheta=360.0/self.N_column
        dev_shift= (int(self.params_stim['dev_id']/dtheta)-
                    int(self.params_stim['std_id']/dtheta))*self.N_eachC

        PltTool.PlotPyr_current(shift_id=-self.N_eachC//2)
        PltTool.PlotPyr_current(shift_id=self.N_eachC//2-1)
        PltTool.PlotPyr_conductance(shift_id=-self.N_eachC//2)
        PltTool.PlotPyr_conductance(shift_id=self.N_eachC//2-1)
        PltTool.PlotPyr_conductance(shift_id=dev_shift+self.N_eachC//2-1)
        PltTool.PlotPyr_conductance(shift_id=dev_shift-self.N_eachC//2)

    
        #PltTool.PlotInt_current()  
        #PltTool.PlotP_current(shift_id=-self.N_eachC//2)
        #PltTool.PlotS_current(shift_id=-self.N_eachC//4)
    
        PltTool.PlotPyrTimes_col()



        self.PltTool=PltTool
        '''
        if not self.omission_flag:
            return DD_value
        else:
            return DD_value, n_omiss
        '''
        
    def selectAnalyz(self,Plot_flag=True): 
        import AnalysisTool as AT
        if self.Input_blocker is None:
            save_dir= ('../figs/'+self.saving_tag+'/')
        else:
            save_dir= ('../figs/blocker/'+self.saving_tag+'/')        
        PltTool=AT.Tools(self.Output_dic, self.ndt,self.n_frames, self.params_sim['dt'], self.params_stim, Plot_flag=Plot_flag,omission_flag=self.omission_flag, save_dir=save_dir)  # Can add a output directory at the end.
        PltTool.PlotPyrPop(left=-1)

        #PltTool=AT.Tools(self.Output_dic, self.ndt,self.n_frames, self.params_sim['dt'], self.params_stim, Plot_flag=Plot_flag,omission_flag=self.omission_flag, save_dir= ('../figs/blocker/'+self.saving_tag+'/'))  # Can add a output directory at the end.
        return PltTool
if __name__=="__main__":
    # some parameters can be changed here. The parameters are changed to oddball parameter set. 
    main_model=oddball_model()
    Tinter=0.5 
    Tstim=0.5

#    Tau_int=1
#    g_IntE=15         # first bunch 15
#    g_back_int=44        # first bunch 44

#    nbr_rep_std=70
#    nbr_rep_dev=30
#    prob_flag= True
    nbr_rep_std=100   #100
    nbr_rep_dev=100  #100
    Input_blocker=None
    prob_flag= True

    '''
    # g_IntE=15, g_back_int=44 seem better.
    # For omission, default values are:  Tinter=0.2   Tstim=0.1  Tau_int=0.1
    # Stim_g=1.5 , ,g_IntE=6,  g_back_int=53     
    if Input_blocker is None:
        saving_tag='omission_noHomeoN200'+f'Tinter{Tinter:.1f}'
    else:
        saving_tag='omission_noHomeoN200'+f'Input_blocker{Input_blocker:.1f}'
    saving_tag=saving_tag,
    Tau_int=Tau_int,g_IntE=6,Stim_g=1.5,
    '''    
    #Output_dic=main_model.main(omission_flag=False,saving_tag=f'oddball_manyctrl',adap_amount=0,
    #                           Tau_int=Tau_int, Tinter=Tinter,Tstim=Tstim, 
    #                           nbr_rep_std=nbr_rep_std,nbr_rep_dev=nbr_rep_dev,
    #                           g_IntE=g_IntE,g_back_int=g_back_int,prob_flag=prob_flag,
    #                           Input_blocker=Input_blocker)
    Output_dic=main_model.main(omission_flag=False,saving_tag=f'oddball_manyctrl_normset',adap_amount=0,
                               Tinter=Tinter,Tstim=Tstim, 
                               nbr_rep_std=nbr_rep_std,nbr_rep_dev=nbr_rep_dev,
                               prob_flag=prob_flag,
                               Input_blocker=Input_blocker)    
    #Output_dic=main_model.main(omission_flag=False,saving_tag='oddball')
    #Output_dic=main_model.main(omission_flag=True,saving_tag='oddball_omissctrl')
    #Output_dic=main_model.main(Trained_Pathway='\\net\\Output_noHomeoN200.pkl', flag_homeo=False,saving_tag='oddballnoHomeoN200')
  
    main_model.analyze(Plot_flag=True)
    PltTool=main_model.selectAnalyz()
    PltTool.PlotPyrPop(left=50, right=60)
    PltTool.Plot_nRep(left=50, right=60)

    # Ploting goal,first vs. Many ctrl of Integrator, Pop average, Res. Pop average. 
    PltTool.PlotPyrPop(left=-1)
    MMN_Int=PltTool.PlotInt_nrep(nrep=Input_blocker)
    MMN_pop=np.mean(PltTool.Odd_rate_vec)-np.mean(PltTool.Last_rate_vec)
    PltTool.PlotPyrTimes_pop_all_long(left=0, right=10,nRep=Input_blocker)   # Requires param_stim['choices']



    PltTool.PlotPyr_conductance(shift_id=-PltTool.N_eachC//2)
    #PltTool.PlotPyr_conductance(shift_id=PltTool.N_eachC//2-1)
    PltTool.PlotPyrTimes_col()
    print(f'MMN_Int{MMN_Int:2f},MMN_pop{MMN_pop:2f}')

    PltTool.PlotInt(left=0,right=10)
    #PltTool.PlotPyrTimes_pop_all_long(left=50, right=60)   # Requires param_stim['choices']
    #PltTool.PlotInt(left=50,right=60)
    
    
    # Count occurrences of 2 and 3
    choices= PltTool.Paras['choices']    # 1 for std, 2 for dev

    
''' 
import os
import pickle

current_dir = os.getcwd()
parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))

with open(parent_dir+'\\net\\Output_oddball_omissctrl.pkl', 'wb') as f:
with open(parent_dir+'\\net\\Output_oddball.pkl', 'wb') as f:

temp_prob=nbr_rep_std/(nbr_rep_std+nbr_rep_dev)
with open(parent_dir+'\\net\\prop\\Output_oddball.pkl'+f'prob{temp_prob:.1f}', 'wb') as f:

with open(parent_dir+'\\net\\realistic\\Output_oddball_normset_ctrl.pkl', 'wb') as f:

    pickle.dump(main_model, f)
    pickle.dump(Output_dic, f)
    pickle.dump(PltTool, f)
    
    
    
'''
