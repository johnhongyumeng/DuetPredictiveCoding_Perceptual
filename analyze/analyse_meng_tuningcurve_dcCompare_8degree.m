% This is preprocessing the data. Change the key work regular_oddball from
% all this to analyze omission data. 
% 02042025 Include the data to calculate the turning curve of every cell.  
% the manystandards control segment, the following numbers refer to the orientations presented 
% 1=90 degrees; 2=0 degrees; 3=22.5 deg; 4=67.5deg; 5=45deg;  6=135deg; 7=157.5; 8=112.5deg.
% Notice begin at 15.

% I have an old one that is:
% 1=90 degrees; 2=0 degrees; 3=22.5 deg; 4=112.5deg; 5=45deg; 6=135deg; 7=157.5; 8=67.5deg.
% And if that is the one i used, then the original 3 and 4 codes were correct.
% 
% But the newer one is 
% 1=90 degrees; 2=0 degrees; 3=22.5 deg; 4=67.5deg; 5=45deg; 6=135deg; 7=157.5; 8=112.5deg.
% 

clear
% load('omission_and_oddball_workspace1.mat')
% file='regular_oddball_data_JMENG1';
% 



tagaddon='Unpublished8degree';
files = {'CA_workspace_SORTED_regular_oddball_data_JMENG1'};
flag_degree=1;

% tagaddon='PYR_maus18degree';
% files = {'CA_workspace_SORTEDrescore_031521_vmmn_PYR_maus1_001'};
% flag_degree=2;

% tagaddon='PYR_maus28degree';
% files = {'CA_workspace_SORTEDrescore_032321_vmmn_PYR_maus2_002'};
% flag_degree=3;


file_tag=['./FigsTest/' tagaddon];
load(files{1},'vectrials','framerate','phzstd','phases','dfofanalyse','stimvisMAIN','stimvistypeMAIN');


data_analyze=dfofanalyse; 
stimvistypeall=stimvistypeMAIN;
% phases=regular_oddball_runs;
n_frame=length(phases);

Time_axis= 1/28*(1:n_frame);

figure(50)
clf
plot(Time_axis,stimvistypeall,'k')
title('Stimulus')


if flag_degree==1
    stimvistypeall(1,stimvistypeall==20 & phases==1 )=5;
    stimvistypeall(1,stimvistypeall==21 & phases==1 )=6;
elseif flag_degree==2
    stimvistypeall(1,stimvistypeall==20 & phases==1 )=1;
    stimvistypeall(1,stimvistypeall==21 & phases==1 )=2;    
elseif flag_degree==3
    stimvistypeall(1,stimvistypeall==20 & phases==1 )=3;
    stimvistypeall(1,stimvistypeall==21 & phases==1 )=4;    
end


figure(51)
clf
plot(Time_axis,stimvistypeall,'k')
title('Stimulus')



framerate=28;
n_neu=size(data_analyze,1);

% starting Index:

flag_startind= (stimvistypeall(2:end)-stimvistypeall(1:end-1)>0.5)...
              & phases(2:end)==1;   
ind_all_array=2: n_frame;
ind_start=ind_all_array(flag_startind);

% starting Index for deviant:

flag_startind_dev= (stimvistypeall(2:end)==15) & (stimvistypeall(1:end-1)==0);
ind_start_dev=ind_all_array(flag_startind_dev);

hold on
plot(Time_axis(ind_start),0*Time_axis(ind_start)+5,'r.')
plot(Time_axis(ind_start_dev),0*Time_axis(ind_start_dev)+5,'b.')

title('Stimulus')
hold off


n_stim= length(ind_start);
vecdat_dev=nan(n_neu, framerate*1.5+1,n_stim);
vecdat_meanNorm=nan(n_neu, framerate*1.5+1,n_stim);

vecstim=nan(1,n_stim);

% Need some normalization. Slide data is done here. Now let's do the
% normalization. When averaging the mean, it can be done by baseline per
% stim, or per cell. Here my hypothesis to test is wheter the stimulus
% trigger a sig. response, so the data should be normalized per trial. 
% So here, the data is normalized by the baseline before the stimulus. 
baseBEGIN=7; baseEND=13; %what you call the baseline
BEGINo=15; ENDo=28; %these are the analysis windows. note: stimulus starts at 15 and finishes at 28

tmpbasos_all=data_analyze(:,phases==1);
norm_vec= nan(n_neu,2);  % Baseline and Std.
tempbaseline=mean(mean(vecdat_dev(:,baseBEGIN:baseEND,:),3),2);
norm_vec(:,1)=tempbaseline;

% Mean is normalized by the baseline just before the response
% i.e., the response here is stimuli triggered response

for i_stim=1 : n_stim
    t_ind=ind_start(i_stim);
    vecdat_dev(:,:,i_stim)=data_analyze(:,t_ind-(framerate*.5):t_ind+framerate);
    vecdat_meanNorm(:,:,i_stim)=vecdat_dev(:,:,i_stim)-mean(vecdat_dev(:,baseBEGIN:baseEND,i_stim),2);
    vecstim(i_stim)=stimvistypeall(t_ind);
end
% The z score is calculated by normalizing over the noise of the whole cell
vecdat_Norm=nan(n_neu, framerate*1.5+1,n_stim); % mean and std norm

for i_cell=1:n_neu
    tmpbasos=tmpbasos_all(i_cell,:);
    tmpbasos=tmpbasos(tmpbasos>0);
    tmpbasos=tmpbasos(tmpbasos<prctile(tmpbasos,92));
    norm_vec(i_cell,2)=std(tmpbasos);
    vecdat_Norm(i_cell,:,:)=vecdat_meanNorm(i_cell,:,:)/norm_vec(i_cell,2);
end


orient_ave_data=nan(n_neu, framerate*1.5+1,8);
orient_num_sample=nan(1,8);


% Plot the sample data
% For the unpublished data 
% 1=90 degrees; 2=0 degrees; 3=22.5 deg; 4=67.5deg; 5=45deg;  6=135deg; 7=157.5; 8=112.5deg.
% mapped to
%   5;          1;             2;         4;         3;          7;       8;        6       
% For the maus1 and maus2: Maybe I want to check the tuning curve using the
% later one. 
% 1=90 degrees; 2=0 degrees; 3=22.5 deg; 4=112.5deg; 5=45deg;  6=135deg; 7=157.5; 8=67.5deg.
% mapped to
%   5;          1;             2;         6;         3;          7;       8;        4;       
% Consider align the small degree to poistion 3, then
% In the case of flag_degree=2, norm: 1 =90 degree maped to 3, dev 2= 0 to 7
% Then 1=90 degrees; 2=0 degrees; 3=22.5 deg; 4=112.5deg; 5=45deg;  6=135deg; 7=157.5; 8=67.5deg.
% mapped to
%   3;          7;             8;         4;         1;          5;       6;        2;       

% In the case of flag_degree=3, norm: 3 =22.5 degree maped to 3, dev 4= 112.5 to 7
% 1=90 degrees; 2=0 degrees; 3=22.5 deg; 4=112.5deg; 5=45deg;  6=135deg; 7=157.5; 8=67.5deg.
% mapped to
%   6;          2;             3;         7;         4;          8;       1;        5;       





% This is mapped to 0 to 135.5  
if flag_degree==1 
    map_vector= [5,1,2,4,3,7,8,6];   
elseif flag_degree==2   % Switched to map the norm at 2.
    % map_vector= [5,1,2,6,3,7,8,4];
    map_vector= [3,7,8,4,1,5,6,2];   
elseif flag_degree==3   % Switched to map the norm at 2.
    % map_vector= [5,1,2,6,3,7,8,4];
    map_vector= [6,2,3,7,4,8,1,5];   
end

% mapped to stimulus starting from smaller  

mapped_stim= map_vector(vecstim);

for i_degree= 1: 8
    temp_data= vecdat_Norm(:,:, mapped_stim==i_degree);
    orient_num_sample(i_degree)= size(temp_data,3);
    orient_ave_data(:,:,i_degree)=mean(temp_data,3);
end


% Now time to calculate the significance of the response.

Tuning_response_pre_mat=nan(n_neu,8); % Activity 
Tuning_response_post_mat=nan(n_neu,8); % Activity 

Tuning_res_sig_mat=nan(n_neu,8);  % significance
Tuning_res_pvalue_mat=nan(n_neu,8);  % significance
Tuning_res_BHpvalue_mat=nan(n_neu,8);  % significance

for i_cell=1:n_neu
    for i_degree=1:8
        temp_data= vecdat_Norm(i_cell,:, mapped_stim==i_degree);
        y_pre=squeeze(mean(temp_data(1,baseBEGIN:baseEND,:),2));
        y_post=squeeze(mean(temp_data(1,BEGINo:ENDo,:),2));
        [h, p, ci, stats] = ttest(y_post, y_pre);
        Tuning_response_pre_mat(i_cell,i_degree)=mean(y_pre);
        Tuning_response_post_mat(i_cell,i_degree)=mean(y_post);
        Tuning_res_pvalue_mat(i_cell,i_degree)=p;

    end
end

% BH correction for each degree.
for i_degree=1:8
    p_values=squeeze(Tuning_res_pvalue_mat(:,i_degree));
    adj_p_values = mafdr(p_values, 'BHFDR', true); % Apply Benjamini-Hochberg correction
    Tuning_res_BHpvalue_mat(:,i_degree)=adj_p_values;  % significance
    
    for i_cell=1:n_neu
        p_adj = Tuning_res_BHpvalue_mat(i_cell, i_degree); % Get adjusted p-value
        % Apply significance criteria
        if p_adj >= 0.05
            Tuning_res_sig_mat(i_cell, i_degree) = 0;
        elseif Tuning_response_pre_mat(i_cell, i_degree) < Tuning_response_post_mat(i_cell, i_degree)
            Tuning_res_sig_mat(i_cell, i_degree) = 1;
        else
            Tuning_res_sig_mat(i_cell, i_degree) = -1;
        end    
    end
end



open('Tuning_res_sig_mat')
open('Tuning_res_BHpvalue_mat')
            




%% Dev analysis starting here
% Including the data during deviant for analyze. 

tmpbasos_all_dev=data_analyze(:,phases==2 | phases==3);

norm_dev_vec= nan(n_neu,2);  % Baseline and Std.

n_stim_dev= length(ind_start_dev);
vecdat_dev=nan(n_neu, framerate*1.5+1,n_stim_dev);
vecdat_meanNorm_dev=nan(n_neu, framerate*1.5+1,n_stim_dev);
vecphase_dev=nan(1,n_stim_dev);


for i_stim=1 : n_stim_dev
    t_ind=ind_start_dev(i_stim);
    vecdat_dev(:,:,i_stim)=data_analyze(:,t_ind-(framerate*.5):t_ind+framerate);
    vecdat_meanNorm_dev(:,:,i_stim)=vecdat_dev(:,:,i_stim)-mean(vecdat_dev(:,baseBEGIN:baseEND,i_stim),2);
    vecphase_dev(i_stim)=phases(t_ind);
end


tempbaseline_dev=mean(mean(vecdat_dev(:,baseBEGIN:baseEND,:),3),2); % This is averaged baseline for each cell
norm_vec_dev(:,1)=tempbaseline;

% The z score is calculated by normalizing over the noise of the whole cell
vecdat_Norm_dev=nan(n_neu, framerate*1.5+1,n_stim_dev); % mean and std norm

for i_cell=1:n_neu
    tmpbasos_dev=tmpbasos_all_dev(i_cell,:);
    tmpbasos_dev=tmpbasos_dev(tmpbasos_dev>0);
    tmpbasos_dev=tmpbasos_dev(tmpbasos_dev<prctile(tmpbasos_dev,92));
    norm_vec_dev(i_cell,2)=std(tmpbasos);
    vecdat_Norm_dev(i_cell,:,:)=vecdat_meanNorm_dev(i_cell,:,:)/norm_vec_dev(i_cell,2);
end

% No need for orientation here. Just for phase phase=2, stim=7. phase=3, stim=3 
% Still need somehting to do sanity check. 

% phase=3, stim=3 
Dev_ave_data=nan(n_neu, framerate*1.5+1,2);

temp_data= vecdat_Norm_dev(:,:, vecphase_dev==3);
Dev_num_sample=[0,0];
Dev_num_sample(2)= size(temp_data,3);
Dev_ave_data(:,:,2)=mean(temp_data,3);

temp_data= vecdat_Norm_dev(:,:, vecphase_dev==2);
Dev_num_sample(1)= size(temp_data,3);
Dev_ave_data(:,:,1)=mean(temp_data,3);


Dev_res_pre_mat=nan(n_neu,2);  % Dev activity
Dev_res_post_mat=nan(n_neu,2);  % Dev activity
Dev_ctrl_pre_mat=nan(n_neu,2);  % Same as many ctrl. Copied here for readibility.

Dev_res_pvalue_mat=nan(n_neu,2);  % Dev vs pre-stim p_value
Dev_ctrl_pvalue_mat=nan(n_neu,2);  % Dev vs baseline p_value
Dev_res_BHpvalue_mat=nan(n_neu,2);  % BH correction  p_value
Dev_ctrl_BHpvalue_mat=nan(n_neu,2);  % 


Dev_res_sig_mat=nan(n_neu,2);  % significance
Dev_ctrl_sig_mat=nan(n_neu,2);  % significance

% Loop to test Dev vs Ctrl and Dev vs Pre-stim. Baseline
for i_cell=1:n_neu
    for i_phase=2:3
        %  phase phase=2, stim=7. 
        if i_phase==2
            i_degree=7;
        elseif i_phase==3
            i_degree=3;
        end
        temp_data= vecdat_Norm_dev(i_cell,:, vecphase_dev==i_phase);
        y_pre=squeeze(mean(temp_data(1,baseBEGIN:baseEND,:),2));
        y_post=squeeze(mean(temp_data(1,BEGINo:ENDo,:),2));
        [h, p, ci, stats] = ttest(y_post, y_pre);
        Dev_res_pre_mat(i_cell,i_phase-1)=mean(y_pre);
        Dev_res_post_mat(i_cell,i_phase-1)=mean(y_post);
        Dev_res_pvalue_mat(i_cell,i_phase-1)=p;

        temp_data2=vecdat_Norm(i_cell,:, mapped_stim==i_degree);
        y_pre2=squeeze(mean(temp_data2(1,BEGINo:ENDo,:),2));
        Dev_ctrl_pre_mat(i_cell,i_phase-1)=mean(y_pre2);

        [h2, p2, ci2, stats2] = ttest2(y_pre2, y_post, 'Vartype', 'unequal');
        Dev_ctrl_pvalue_mat(i_cell,i_phase-1)=p2;
    end
end


% BH correction for each degree. This is wrong. I need to cross validate
% with other forms. 
for i_phase=2:3

    p_values=squeeze(Dev_res_pvalue_mat(:,i_phase-1));
    adj_p_values = mafdr(p_values, 'BHFDR', true); % Apply Benjamini-Hochberg correction
    Dev_res_BHpvalue_mat(:,i_phase-1)=adj_p_values;  % significance

    p_values=squeeze(Dev_ctrl_pvalue_mat(:,i_phase-1));
    adj_p_values = mafdr(p_values, 'BHFDR', true); % Apply Benjamini-Hochberg correction
    Dev_ctrl_BHpvalue_mat(:,i_phase-1)=adj_p_values;  % significance


    for i_cell=1:n_neu
        % Sig for Dev vs baseline
        p_adj = Dev_res_BHpvalue_mat(i_cell, i_phase-1); % Get adjusted p-value
        % Apply significance criteria
        if p_adj >= 0.05
            Dev_res_sig_mat(i_cell, i_phase-1) = 0;
        elseif Dev_res_pre_mat(i_cell, i_phase-1) < Dev_res_post_mat(i_cell, i_phase-1)
            Dev_res_sig_mat(i_cell, i_phase-1) = 1;
        else
            Dev_res_sig_mat(i_cell, i_phase-1) = -1;
        end    

        % Sig for Dev vs Ctrl
        p_adj = Dev_ctrl_BHpvalue_mat(i_cell, i_phase-1); % Get adjusted p-value
        % Apply significance criteria
        if p_adj >= 0.05
            Dev_ctrl_sig_mat(i_cell, i_phase-1) = 0;
        elseif Dev_ctrl_pre_mat(i_cell, i_phase-1) < Dev_res_post_mat(i_cell, i_phase-1)
            Dev_ctrl_sig_mat(i_cell, i_phase-1) = 1;
        else
            Dev_ctrl_sig_mat(i_cell, i_phase-1) = -1;
        end    


    end
end

open('Dev_ctrl_sig_mat')
open('Dev_res_sig_mat')
open('Tuning_response_post_mat')



check_cell=5;   % The lower the baseline firing rate, the stronger the responsive can be, which make sense. 
time_axis_trial= (-(framerate*.5):framerate)/framerate;
figure(700+check_cell)
for i_degree=1:8
    subplot(1,8,i_degree)
    plot(time_axis_trial,orient_ave_data(check_cell,:,i_degree))
end


figure(800+check_cell)
for i_dev=1:2
    if i_dev==1
        i_degree=7;
    elseif i_dev==2
        i_degree=3;
    end
    subplot(2,2,i_dev)
    plot(time_axis_trial,orient_ave_data(check_cell,:,i_degree))
    title(['control response',num2str(i_degree)])
    subplot(2,2,i_dev+2)
    plot(time_axis_trial,Dev_ave_data(check_cell,:,i_dev))
    title(['dev response',num2str(i_degree)])

end


% All right, now I need to think a way combine all these analysis, and
% report them. First save into a matrix, then export into cvs. 

Results_table=nan(n_neu,4+3+3,2); % Orientation 1,2 (pre, ctrl, pre, dev), P1(ctrl), P2(dev), P3 (dev-ctrl), Sig1,S2,S3.

column_names = {'Pre1', 'Ctrl', 'Pre2', 'Dev', ...
                'P1: Ctrl', 'P2: Dev', 'P3: Dev-Ctrl', ...
                'Sig1: Ctrl', 'Sig2: Dev', 'Sig3: Dev-Ctrl'};


Results_table(:,1,1)=Tuning_response_pre_mat(:,3);
Results_table(:,2,1)=Tuning_response_post_mat(:,3);
Results_table(:,3,1)=Dev_res_pre_mat(:,2);  
Results_table(:,4,1)=Dev_res_post_mat(:,2);

Results_table(:,5,1)=Tuning_res_BHpvalue_mat(:,3);  
Results_table(:,6,1)=Dev_res_BHpvalue_mat(:,2);
Results_table(:,7,1)=Dev_ctrl_BHpvalue_mat(:,2);

Results_table(:,8,1)=Tuning_res_sig_mat(:,3);  
Results_table(:,9,1)=Dev_res_sig_mat(:,2);
Results_table(:,10,1)=Dev_ctrl_sig_mat(:,2);

Results_table(:,1,2)=Tuning_response_pre_mat(:,7);
Results_table(:,2,2)=Tuning_response_post_mat(:,7);
Results_table(:,3,2)=Dev_res_pre_mat(:,1);  
Results_table(:,4,2)=Dev_res_post_mat(:,1);

Results_table(:,5,2)=Tuning_res_BHpvalue_mat(:,7);  
Results_table(:,6,2)=Dev_res_BHpvalue_mat(:,1);
Results_table(:,7,2)=Dev_ctrl_BHpvalue_mat(:,1);

Results_table(:,8,2)=Tuning_res_sig_mat(:,7);  
Results_table(:,9,2)=Dev_res_sig_mat(:,1);
Results_table(:,10,2)=Dev_ctrl_sig_mat(:,1);

% Extract slices from the third dimension
degree_names={'0 degree', '22.5 degree', '45 degree', '77.5 degree', ...
                '90 degree', '112.5 degree', '135 degree', '147.5 degree'};
data1 = array2table(Results_table(:,:,1), 'VariableNames', column_names);
data2 = array2table(Results_table(:,:,2), 'VariableNames', column_names);
data3 = array2table(Tuning_response_post_mat, 'VariableNames', degree_names);


filename = ['./Figs/' tagaddon 'Results.xlsx']; % Define output file name
% Only include the 8 degrees incase I made a mistake.

% Write the first dataset to "Sheet1"
% writetable(data1, filename, 'Sheet', '45 degree');

% Write the second dataset to "Sheet2"
% writetable(data2, filename, 'Sheet', '135 degree');
writetable(data3 ,filename, 'Sheet','Eight Degrees')

disp('Excel file successfully created with separate sheets!');


save(['./Figs/' tagaddon 'TuningCurve.mat'],...
    'Tuning_response_pre_mat','Tuning_response_post_mat',...
    'Tuning_res_sig_mat','Tuning_res_BHpvalue_mat', ...
    'Results_table');        



