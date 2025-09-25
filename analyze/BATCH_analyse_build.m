%%% Modified by JOhn on 020192025. Only include the part for building the
%%% matrix. Run this before specific contingences.
%%% do stringent test and pre/post for supplemental figure, to show
%%% stability of conclusions. VIPs=no mod. SSTs=slight DD, strong SSA.
%%% PYRs=strong DD, strong SSA jojo
% clcl
close all
clear 
for file_load=1

% files = {'regular_oddball_data_JMENG1'};
% files = {'omission_oddball_data_JMENG1'};


% files = {'CA_workspace_SORTED_regular_oddball_data_JMENG1'};

files = {'CA_workspace_SORTEDrescore_031521_vmmn_PYR_maus1_001'};
% files = {'CA_workspace_SORTEDrescore_032321_vmmn_PYR_maus2_002'};

% files = {'CA_workspace_SORTEDrescore_inclloco_102020_vmm_vglut_001'};
% files = {'CA_workspace_SORTEDrescore_inclloco_102020_vmm_vglut_002'};
% files = {'CA_workspace_SORTEDrescore_inclloco_102020_vmm_vglut_003'};
% files = {'CA_workspace_SORTEDrescore_inclloco_102020_vmm_vglut_004'};
end

limtrials1=99; %limit the number of trials analysed for each stimulus type(control/redundandt/deviant). IF set high, then it normalizes to the smallest number of trials in any condition (cont, rdnt, dev)
baseBEGIN=7; baseEND=13; %what you call the baseline
base_vec=baseBEGIN:baseEND;
BEGINo=15; ENDo=28; %these are the analysis windows. note: stimulus starts at 14 and finishes at 28
onestim=0; %if 0, then plot responses to both stimuli (e.g. 45 and 135 deg). if 1, then plot only to first orientation. if 2, then plot to second
dobeforeafter=0; %if this is "1", then it compares the first trials to the last for the test. if it's zero, it does even/odd
stdfromaverage=0; %normalize based on average, rather than background Df. not advisable, but wouth a look if activity is very noisy at rest.

%%%%%%%this compiles all the files into three matrices
for buildingmatrices=1;
    effectsallBUILD=[]; %even
    effectsallTEST=[]; %odd
    singtrialsall=[];   % Not used here
    count=0; counttest=0;
    effectsall=[];  %this is the main varable
    effectsallSTD=[];  %this is for limiting trials
    actcont=[];
    mousevec=[];   % Not used here
    countall=0;
end
for file=1:size(files,1)   
        load(files{file},'vectrials','framerate','phzstd','phases','dfofanalyse','stimvisMAIN','stimvistypeMAIN');
      
    % 
    figure(1999);  % Get a sense of the z-score distribution. Need to know the z-score definition. 
    subplot(1,2,1)
    data = vectrials{1,1}(:);
    histogram(data);
    title('Vectrials stim 1')
    subplot(1,2,2)
    data = vectrials{2,1}(:);
    histogram(data);
    title('Vectrials stim 1')



    %figure; plot(stimvistypeMAIN); title(num2str(file));
    disp(files{file});
    for remove_first_trial_of_each_type=0
        if remove_first_trial_of_each_type==1
            if limtrials1>0
                limtrials=limtrials1;
                for ss1=1:2;
                    vecstim=vectrials{ss1,2};
                    tmps1=[];
                    for tr1=1:size(vecstim,2);
                        
                        if and(sum(vecstim(tr1)==tmps1)==0,vecstim(tr1)~=15)
                            tmps1=horzcat(tmps1,vecstim(tr1));
                            vecstim(tr1)=99;
                        end
                        
                    end
                    for additionally_remove_first_few_controls=1;
                        tmps1=[]; gh=size(vecstim,2)+1;
                        disp(sum(vecstim==(20+(ss1-1)))); 
                        for tr1=1:size(vecstim,2);
                            if and(vecstim(gh-tr1)>15,vecstim(gh-tr1)<22)
                                if sum(tmps1==1)<(limtrials)
                                    tmps1=horzcat(tmps1,1);
                                else
                                    vecstim(gh-tr1)=99;
                                end
                            end
                        end
                    end
                    vectrials{ss1,2}=vecstim;
                end
            else
                for ss1=1:2;
                    vecstim=vectrials{ss1,2};
                    tmps1=[];
                    for tr1=1:size(vecstim,2);
                        
                        if and(sum(vecstim(tr1)==tmps1)==0,vecstim(tr1)~=15)
                            tmps1=horzcat(tmps1,vecstim(tr1));
                            vecstim(tr1)=99;
                        end
                        
                    end
                    vectrials{ss1,2}=vecstim;
                end
            end
        end
    end
    
  
    for analysis=1
        
        samp=1/framerate;
        timeaxis=-(samp*(framerate/2)):samp:(samp*framerate); timeaxis=timeaxis+.088;
        
        for stim=1:2;
            limtrials=limtrials1;
            vecdat=vectrials{stim,1};
            mousetmp=zeros(size(vecdat,1),1)+file; 
            vecstim=vectrials{stim,2};
            if limtrials1>0; 
            for checkingfortoo_few_trials=1;
                ss=[1 2 3 4 15 phzstd(stim)];
                for gss=1:6; 
                    sss(gss)=sum(vecstim==ss(gss));
                end
                if min(sss)<limtrials1;
                   limtrials=min(sss); 
                end               
                disp(strcat('stim',num2str(stim),'min number of trials',num2str(limtrials)));     % 16 is the limitation
            end
            
            effects=[];
            for z1=1:8
                tmp1=vecdat(:,:,vecstim==z1);
                if size(tmp1,3)>(limtrials-1)
                    effects(:,:,z1+1)=mean(tmp1(:,:,1:limtrials),3);
                else
                    effects(:,:,z1+1)=mean(tmp1(:,:,1:end),3);
                end
            end
            tmp1=vecdat(:,:,vecstim==15);
            if size(tmp1,3)>(limtrials-1)
                effects(:,:,10)=mean(tmp1(:,:,1:limtrials),3);
            else
                effects(:,:,10)=mean(tmp1(:,:,1:end),3);
            end
            tmp1=vecdat(:,:,vecstim==phzstd(stim));
            
                if size(tmp1,3)>(limtrials-1)
                    effects(:,:,1)=mean(tmp1(:,:,1:limtrials),3);
                     else
                    effects(:,:,1)=mean(tmp1(:,:,1:end),3);
                end
                
            else
                 vecdat=vectrials{stim,1};
                vecstim=vectrials{stim,2};
                effects=[];
                stimax=[phzstd(stim) 1 2 3 4 5 6 7 8 15];
                for z1=1:10
                    effects(:,:,z1)=mean(vecdat(:,:,vecstim==stimax(z1)),3);
                end
            end
                
            
            for c1=1:size(effects,1) %baseline correction
                basos=[1 2 2 2 2 2 2 2 2 3; 1 3 3 3 3 3 3 3 3 2];
                at=squeeze(effects(c1,baseBEGIN:baseEND,:));
                
                for z1=1:size(effects,3)
                    % mnn=min(effects(c1,baseBEGIN:baseEND,z1),[],2); 
                    mnn=mean(effects(c1,baseBEGIN:baseEND,z1),2); 
                    tmpbasos=dfofanalyse(c1,phases==basos(stim,z1));
                    tmpbasos=tmpbasos(tmpbasos>0); tmpbasos=tmpbasos(tmpbasos<prctile(tmpbasos,92));
                    stt=std(tmpbasos); 
                    if stdfromaverage==1; stt=std(at(:)); end
                    
                    if stt<.001; stt=std(effects(c1,ENDo:end,z1)); end
                    if stt<.001; stt=10000; end
                    for t1=1:size(effects,2)
                        effectsSTD(c1,t1,z1)=(effects(c1,t1,z1)-mnn)./stt;
                        effects(c1,t1,z1)=(effects(c1,t1,z1)-mnn);
                    end
                    mnns(file,stim,c1,z1)=mnn;
                    stts(file,stim,c1,z1)=stt;
                end
                
            end
            % Why calculated mean and stt here. So untuituive.
            
            for getspontaneousactivity=1
                if stim==1
                    for c1=1:size(effects,1)
                        tmp=[];
                        for ph=1:3
                            
                            tmp2=dfofanalyse(c1,phases==ph);
                            tmp3=stimvisMAIN(1,phases==ph);
                            t1=0; t2=0; %t1 becomes teh beginning of stims. t2 becomes the end.
                            for t=100:size(tmp3,2)-100
                                if and(t1==0,tmp3(1,t)>0)
                                    t1=t;
                                end
                                if and(t2==0,t1>0)
                                    if max(tmp3(1,t-99:t+99))==0
                                        t2=t-99;
                                    end
                                end
                            end
                            if t2==0; t2=size(tmp3,2)-99; end
                            
                          tmp22=tmp2(1,tmp2(1,:)>0); pr=prctile(tmp22,95); tmp2=tmp2/std(tmp22<pr);
                            
                            for alltimepoints=0 %if you want to just analyse interstimulus activity, make this =0
                                if alltimepoints==1
                                    bin1=floor((t2-t1)/10);
                                    for b=1:10;
                                        tmp=horzcat(tmp,mean(tmp2(1,((b-1)*bin1)+t1+1:(b*bin1)+t1)));
                                    end
                                else
                                    tmp2=tmp2(1,t1:t2); tmp3=tmp3(1,t1:t2);
                                    tmp2=tmp2(1,tmp3==0);
                                    bin1=floor(size(tmp2,2)/10);
                                    for b=1:10
                                        tmp=horzcat(tmp,mean(tmp2(1,((b-1)*bin1)+1:(b*bin1))));
                                    end
                                end
                            end
                            
                        end
                        actcont=vertcat(actcont,tmp);
                    end
                end
            end
            
            %for getting single trial data. make sure to comment out the
            %"exclude first few trials" section above if you want to really see the truth!
            tmpcells=[];
            for cc=1:size(vecdat,1)
                
                tmp00=dfofanalyse(cc,30:449)./stts(file,stim,cc,1); 
                tmp1=squeeze(vecdat(cc,:,vecstim==phzstd(stim)))./stts(file,stim,cc,1);
                tmp0=[]; stim0=[];
                for g=1:limtrials;
                    tmp0=vertcat(tmp0,tmp1(:,g));
                    stim0=horzcat(stim0,zeros(1,13),ones(1,14),zeros(1,16));
                end
                stim0=horzcat(zeros(1,420),stim0); tmp0=vertcat(tmp00',tmp0); 
                tmp1=squeeze(vecdat(cc,:,vecstim==1))./stts(file,stim,cc,2);  tmp2=squeeze(vecdat(cc,:,vecstim==2))./stts(file,stim,cc,3);
                tmp3=squeeze(vecdat(cc,:,vecstim==3))./stts(file,stim,cc,4);  tmp4=squeeze(vecdat(cc,:,vecstim==4))./stts(file,stim,cc,5); 
                tmp5=squeeze(vecdat(cc,:,vecstim==5))./stts(file,stim,cc,6);  tmp6=squeeze(vecdat(cc,:,vecstim==6))./stts(file,stim,cc,7); 
                tmp7=squeeze(vecdat(cc,:,vecstim==7))./stts(file,stim,cc,8);  tmp8=squeeze(vecdat(cc,:,vecstim==8))./stts(file,stim,cc,9); 
                tmp15=squeeze(vecdat(cc,:,vecstim==15))./stts(file,stim,cc,10);  ggg=0;
                for g=1:5
                    try
                        tmp0=vertcat(tmp0,tmp1(:,g),tmp2(:,g),tmp3(:,g),tmp4(:,g),tmp5(:,g),tmp6(:,g),tmp7(:,g),tmp8(:,g),tmp15(:,g));
                        
                        stim0=horzcat(stim0,zeros(1,13),ones(1,14)-.5,zeros(1,16),zeros(1,13),ones(1,14)-.5,zeros(1,16),zeros(1,13),ones(1,14)-.5,zeros(1,16),zeros(1,13),ones(1,14)-.5,zeros(1,16),zeros(1,13),ones(1,14)-.5,zeros(1,16),zeros(1,13),ones(1,14)-.5,zeros(1,16),zeros(1,13),ones(1,14)-.5,zeros(1,16),zeros(1,13),ones(1,14)-.5,zeros(1,16),zeros(1,13),ones(1,14)+.5,zeros(1,16));
                        ggg=g;
                    end
                end
                if ggg==5;
                tmpcells=vertcat(tmpcells,tmp0');
                stim10=stim0;
                end
            end
            
            
            
            if onestim==0;
                effectsall((1+countall):(countall+size(effects,1)),1:size(effects,2),2:10,stim)=effects(:,1:end,2:10);
                effectsall((1+countall):(countall+size(effects,1)),1:size(effects,2),1,stim)=effects(:,1:end,1);
                
                singtrialsall((1+countall):(countall+size(effects,1)),1:size(tmpcells,2),stim)=tmpcells;
                mousevec((1+countall):(countall+size(effects,1)),stim)=mousetmp;
            elseif onestim==1;
                if stim==1;
                    effectsall((1+countall):(countall+size(effects,1)),1:size(effects,2),2:10,stim)=effects(:,1:end,2:10);
                    effectsall((1+countall):(countall+size(effects,1)),1:size(effects,2),1,stim)=effects(:,1:end,1);
                    
                    singtrialsall((1+countall):(countall+size(effects,1)),1:size(tmpcells,2),stim)=tmpcells;
                mousevec((1+countall):(countall+size(effects,1)),stim)=mousetmp;
                else
                    effectsall((1+countall):(countall+size(effects,1)),1:size(effects,2),2:10,stim)=effects(:,1:end,2:10)-effects(:,1:end,2:10);
                    effectsall((1+countall):(countall+size(effects,1)),1:size(effects,2),1,stim)=effects(:,1:end,1)-effects(:,1:end,1);
                    singtrialsall((1+countall):(countall+size(effects,1)),1:size(tmpcells,2),stim)=tmpcells-tmpcells;
                    mousevec((1+countall):(countall+size(effects,1)),stim)=mousetmp-mousetmp;
                end
            elseif onestim==2;
                if stim==2;
                    effectsall((1+countall):(countall+size(effects,1)),1:size(effects,2),2:10,stim)=effects(:,1:end,2:10);
                    effectsall((1+countall):(countall+size(effects,1)),1:size(effects,2),1,stim)=effects(:,1:end,1);
                    singtrialsall((1+countall):(countall+size(effects,1)),1:size(tmpcells,2),stim)=tmpcells;
                mousevec((1+countall):(countall+size(effects,1)),stim)=mousetmp;
                else
                    effectsall((1+countall):(countall+size(effects,1)),1:size(effects,2),2:10,stim)=effects(:,1:end,2:10)-effects(:,1:end,2:10);
                    effectsall((1+countall):(countall+size(effects,1)),1:size(effects,2),1,stim)=effects(:,1:end,1)-effects(:,1:end,1);
                    singtrialsall((1+countall):(countall+size(effects,1)),1:size(tmpcells,2),stim)=tmpcells-tmpcells;
                mousevec((1+countall):(countall+size(effects,1)),stim)=mousetmp-mousetmp;
                end
            end

            if onestim==0;
                effectsallSTD((1+countall):(countall+size(effects,1)),1:size(effects,2),2:10,stim)=effectsSTD(:,1:end,2:10);
                effectsallSTD((1+countall):(countall+size(effects,1)),1:size(effects,2),1,stim)=effectsSTD(:,1:end,1);
                
            elseif onestim==1;
                if stim==1;
                    effectsallSTD((1+countall):(countall+size(effects,1)),1:size(effects,2),2:10,stim)=effectsSTD(:,1:end,2:10);
                    effectsallSTD((1+countall):(countall+size(effects,1)),1:size(effects,2),1,stim)=effectsSTD(:,1:end,1);
                    
                else
                    effectsallSTD((1+countall):(countall+size(effects,1)),1:size(effects,2),2:10,stim)=effectsSTD(:,1:end,2:10)-effectsSTD(:,1:end,2:10);
                    effectsallSTD((1+countall):(countall+size(effects,1)),1:size(effects,2),1,stim)=effectsSTD(:,1:end,1)-effectsSTD(:,1:end,1);
                    
                end
            elseif onestim==2;
                if stim==2;
                    effectsallSTD((1+countall):(countall+size(effects,1)),1:size(effects,2),2:10,stim)=effectsSTD(:,1:end,2:10);
                    effectsallSTD((1+countall):(countall+size(effects,1)),1:size(effects,2),1,stim)=effectsSTD(:,1:end,1);
                  
                else
                    effectsallSTD((1+countall):(countall+size(effects,1)),1:size(effects,2),2:10,stim)=effectsSTD(:,1:end,2:10)-effectsSTD(:,1:end,2:10);
                    effectsallSTD((1+countall):(countall+size(effects,1)),1:size(effects,2),1,stim)=effectsSTD(:,1:end,1)-effectsSTD(:,1:end,1);
                    
                end
            end
            
        end
        
        countall=countall+size(effects,1);
    end
    
    for build=1
        
        samp=1/framerate;
        timeaxis=-0.4286:samp:1.0714; timeaxis=timeaxis;
        for stim=1:2
            limtrials=limtrials1;
            if limtrials1>0;
                vecdat=vectrials{stim,1};
                vecstim=vectrials{stim,2};
                for checkingfortoo_few_trials=1;
                    ss=[1 2 3 4 15 phzstd(stim)];
                    for gss=1:6;
                        sss(gss)=sum(vecstim==ss(gss));
                    end
                    if min(sss)<limtrials1;
                        limtrials=min(sss);
                    end
                  
                end
                
                effects=[];
                for z1=1:8
                    tmp1=vecdat(:,:,vecstim==z1);
                    ex=[];
                    for t1=1:size(tmp1,3)
                        if rem(t1,2)==1;
                            ex(t1)=1;
                        else
                            ex(t1)=0;
                        end
                    end
                    if dobeforeafter==1; pt=ceil(limtrials/2); clear ex; ex(1,1:pt)=ones(1,pt); try ex(1,(pt+1):limtrials)=zeros(1,pt); catch ex(1,(pt+1):limtrials)=zeros(1,pt-1); end; end
                    
                    if size(tmp1,3)>(limtrials-1)
                        atmp=(tmp1(:,:,1:limtrials));
                        ex=ex(1:limtrials);
                        effects(:,:,z1+1)=mean(atmp(:,:,ex==1),3);
                    else
                        atmp=(tmp1(:,:,1:end));
                        ex=ex(1:size(atmp,3));
                        effects(:,:,z1+1)=mean(atmp(:,:,ex==1),3);
                    end
                end
                
                tmp1=vecdat(:,:,vecstim==15);
                ex=[];
                for t1=1:size(tmp1,3)
                    if rem(t1,2)==1
                        ex(t1)=1;
                    else
                        ex(t1)=0;
                    end
                end
                if dobeforeafter==1; pt=ceil(limtrials/2); clear ex; ex(1,1:pt)=ones(1,pt); try ex(1,(pt+1):limtrials)=zeros(1,pt); catch ex(1,(pt+1):limtrials)=zeros(1,pt-1); end; end
                
                if size(tmp1,3)>(limtrials-1)
                    atmp=(tmp1(:,:,1:limtrials));
                    ex=ex(1:limtrials);
                    effects(:,:,10)=mean(atmp(:,:,ex==1),3);
                else
                    atmp=(tmp1(:,:,1:end));
                    ex=ex(1:size(atmp,3));
                    effects(:,:,10)=mean(atmp(:,:,ex==1),3);
                end
                tmp1=vecdat(:,:,vecstim==phzstd(stim));
                ex=[];
                for t1=1:size(tmp1,3)
                    if rem(t1,2)==1;
                        ex(t1)=1;
                    else
                        ex(t1)=0;
                    end
                end
                if dobeforeafter==1; pt=ceil(limtrials/2); clear ex; ex(1,1:pt)=ones(1,pt); try ex(1,(pt+1):limtrials)=zeros(1,pt); catch ex(1,(pt+1):limtrials)=zeros(1,pt-1); end; end
                
                atmp=(tmp1(:,:,1:limtrials));
                ex=ex(1:limtrials);
                effects(:,:,1)=mean(atmp(:,:,ex==1),3);
            else
                vecdat=vectrials{stim,1};
                vecstim=vectrials{stim,2};
                effects=[];
                stimax=[phzstd(stim) 1 2 3 4 5 6 7 8 15];
                for z1=1:10;
                    tmp1=vecdat(:,:,vecstim==stimax(z1));
                    ex=[];
                    for t1=1:size(tmp1,3)
                        if rem(t1,2)==1;
                            ex(t1)=1;
                        else
                            ex(t1)=0;
                        end
                    end
                    tmp1=tmp1(:,:,ex==1);
                    effects(:,:,z1)=mean(tmp1,3);
                end
            end
                
                
            for c1=1:size(effects,1) %baseline correction
                
                
                for z1=1:size(effects,3)
                    %                     mnn=min(effects(c1,baseBEGIN:baseEND,z1),[],2);
                    %                     stt=std(effects(c1,baseBEGIN:baseEND,z1));
                    
                    mnn=mnns(file,stim,c1,z1);
                    stt=stts(file,stim,c1,z1);
                    for t1=1:size(effects,2)
                        effects(c1,t1,z1)=(effects(c1,t1,z1)-mnn)./stt;
                    end
                end
                
            end
            
            if onestim==0;
                effectsallBUILD((1+count):(count+size(effects,1)),1:size(effects,2),2:10,stim)=effects(:,1:end,2:10);
                effectsallBUILD((1+count):(count+size(effects,1)),1:size(effects,2),1,stim)=effects(:,1:end,1);
            elseif onestim==1;
                if stim==1;
                    effectsallBUILD((1+count):(count+size(effects,1)),1:size(effects,2),2:10,stim)=effects(:,1:end,2:10);
                    effectsallBUILD((1+count):(count+size(effects,1)),1:size(effects,2),1,stim)=effects(:,1:end,1);
                else
                    effectsallBUILD((1+count):(count+size(effects,1)),1:size(effects,2),2:10,stim)=effects(:,1:end,2:10)-effects(:,1:end,2:10);
                    effectsallBUILD((1+count):(count+size(effects,1)),1:size(effects,2),1,stim)=effects(:,1:end,1)-effects(:,1:end,1);
                end
            elseif onestim==2;
                if stim==2;
                    effectsallBUILD((1+count):(count+size(effects,1)),1:size(effects,2),2:10,stim)=effects(:,1:end,2:10);
                    effectsallBUILD((1+count):(count+size(effects,1)),1:size(effects,2),1,stim)=effects(:,1:end,1);
                else
                    effectsallBUILD((1+count):(count+size(effects,1)),1:size(effects,2),2:10,stim)=effects(:,1:end,2:10)-effects(:,1:end,2:10);
                    effectsallBUILD((1+count):(count+size(effects,1)),1:size(effects,2),1,stim)=effects(:,1:end,1)-effects(:,1:end,1);
                end
            end
            
        end
        
        count=count+size(effects,1);
    end
    
    
    for test=1
        samp=1/framerate;
        for stim=1:2
            limtrials=limtrials1;
            if limtrials1>0;
                vecdat=vectrials{stim,1};
                vecstim=vectrials{stim,2};
                for checkingfortoo_few_trials=1;
                    ss=[1 2 3 4 15 phzstd(stim)];
                    for gss=1:6;
                        sss(gss)=sum(vecstim==ss(gss));
                    end
                    if min(sss)<limtrials1;
                        limtrials=min(sss);
                    end
                   
                end
                
                effects=[];
                for z1=1:8
                    tmp1=vecdat(:,:,vecstim==z1);
                    ex=[];
                    for t1=1:size(tmp1,3)
                        if rem(t1,2)==1;
                            ex(t1)=0;
                        else
                            ex(t1)=1;
                        end
                    end
                    if dobeforeafter==1; pt=ceil(limtrials/2); clear ex; ex(1,1:pt)=zeros(1,pt); try ex(1,(pt+1):limtrials)=ones(1,pt); catch ex(1,(pt+1):limtrials)=ones(1,pt-1); end; end
                    
                    if size(tmp1,3)>(limtrials-1)
                        atmp=(tmp1(:,:,1:limtrials));
                        ex=ex(1:limtrials);
                        effects(:,:,z1+1)=mean(atmp(:,:,ex==1),3);
                    else
                        atmp=(tmp1(:,:,1:end));
                        ex=ex(1:size(atmp,3));
                        effects(:,:,z1+1)=mean(atmp(:,:,ex==1),3);
                    end
                end
                
                tmp1=vecdat(:,:,vecstim==15);
                ex=[];
                for t1=1:size(tmp1,3)
                    if rem(t1,2)==1;
                        ex(t1)=0;
                    else
                        ex(t1)=1;
                    end
                end
                if dobeforeafter==1; pt=ceil(limtrials/2); clear ex; ex(1,1:pt)=zeros(1,pt); try ex(1,(pt+1):limtrials)=ones(1,pt); catch ex(1,(pt+1):limtrials)=ones(1,pt-1); end; end
                
                if size(tmp1,3)>(limtrials-1)
                    atmp=(tmp1(:,:,1:limtrials));
                    ex=ex(1:limtrials);
                    effects(:,:,10)=mean(atmp(:,:,ex==1),3);
                else
                    atmp=(tmp1(:,:,1:end));
                    ex=ex(1:size(atmp,3));
                    effects(:,:,10)=mean(atmp(:,:,ex==1),3);
                end
                tmp1=vecdat(:,:,vecstim==phzstd(stim));
                ex=[];
                for t1=1:size(tmp1,3)
                    if rem(t1,2)==1;
                        ex(t1)=0;
                    else
                        ex(t1)=1;
                    end
                end
                if dobeforeafter==1; pt=ceil(limtrials/2); clear ex; ex(1,1:pt)=zeros(1,pt); try ex(1,(pt+1):limtrials)=ones(1,pt); catch ex(1,(pt+1):limtrials)=ones(1,pt-1); end; end
                
                atmp=(tmp1(:,:,1:limtrials));
                ex=ex(1:limtrials);
                effects(:,:,1)=mean(atmp(:,:,ex==1),3);
            else
                vecdat=vectrials{stim,1};
                vecstim=vectrials{stim,2};
                effects=[];
                stimax=[phzstd(stim) 1 2 3 4 5 6 7 8 15];
                for z1=1:10;
                    tmp1=vecdat(:,:,vecstim==stimax(z1));
                    ex=[];
                    for t1=1:size(tmp1,3)
                        if rem(t1,2)==1;
                            ex(t1)=0;
                        else
                            ex(t1)=1;
                        end
                    end
                    tmp1=tmp1(:,:,ex==1);
                    effects(:,:,z1)=mean(tmp1,3);
                end
            end
            
            for c1=1:size(effects,1) %baseline correction
                
                
                for z1=1:size(effects,3)
                    %                     mnn=min(effects(c1,baseBEGIN:baseEND,z1),[],2);
                    %                     stt=std(effects(c1,baseBEGIN:baseEND,z1));
                    mnn=mnns(file,stim,c1,z1);
                    stt=stts(file,stim,c1,z1);
                    for t1=1:size(effects,2)
                        effects(c1,t1,z1)=(effects(c1,t1,z1)-mnn)./stt;
                    end
                end
                
            end
            
            
            if onestim==0
                effectsallTEST((1+counttest):(counttest+size(effects,1)),1:size(effects,2),2:10,stim)=effects(:,1:end,2:10);
                effectsallTEST((1+counttest):(counttest+size(effects,1)),1:size(effects,2),1,stim)=effects(:,1:end,1);
            elseif onestim==1;
                if stim==1;
                    effectsallTEST((1+counttest):(counttest+size(effects,1)),1:size(effects,2),2:10,stim)=effects(:,1:end,2:10);
                    effectsallTEST((1+counttest):(counttest+size(effects,1)),1:size(effects,2),1,stim)=effects(:,1:end,1);
                else
                    % Just 0 here. 
                    effectsallTEST((1+counttest):(counttest+size(effects,1)),1:size(effects,2),2:10,stim)=effects(:,1:end,2:10)-effects(:,1:end,2:10);
                    effectsallTEST((1+counttest):(counttest+size(effects,1)),1:size(effects,2),1,stim)=effects(:,1:end,1)-effects(:,1:end,1);
                end
            elseif onestim==2;
                if stim==2;
                    effectsallTEST((1+counttest):(counttest+size(effects,1)),1:size(effects,2),2:10,stim)=effects(:,1:end,2:10);
                    effectsallTEST((1+counttest):(counttest+size(effects,1)),1:size(effects,2),1,stim)=effects(:,1:end,1);
                else
                    % Just 0 here. 
                    effectsallTEST((1+counttest):(counttest+size(effects,1)),1:size(effects,2),2:10,stim)=effects(:,1:end,2:10)-effects(:,1:end,2:10);
                    effectsallTEST((1+counttest):(counttest+size(effects,1)),1:size(effects,2),1,stim)=effects(:,1:end,1)-effects(:,1:end,1);
                end
            end
            
        end
        
        counttest=counttest+size(effects,1);
    end
    
    
end

figure(2000);  % Get a sense of the z-score distribution. Need to know the z-score definition. 
subplot(1,4,1)
data = effectsall(:);
histogram(data);
title('Effects All')
subplot(1,4,2)
data = effectsallSTD(:);
histogram(data);
title('Effects STD')

subplot(1,4,3)
data = effectsallBUILD(:);
histogram(data);
title('Effects BUILD(even)')

subplot(1,4,4)
data = effectsallTEST(:);
histogram(data);
title('Effects TEST (odd)')



% This one generate a sanity test. 
% Notice the order changed. ctrl 1, rdnt 2 to 9, dev 10
% 1st row stim 1, 2nd rwo stim 2. So the neuron sorted based
% So the neuron is sorted by the top-right panel. Let's including
% the sorting part.
% colormap 1 to 4
figure(2001)
ak=mean(effectsallSTD(:,16:35,10,1),2); 
[ix,b]=sort(ak,'descend');



for stim=1:2
 
    for z=1:size(effects,3)
        subplot(2,size(effects,3),z+((stim-1)*size(effects,3))); 
        imagesc(effectsallSTD(b,:,z,stim)); 
        % imagesc(effectsallBUILD(:,:,z,stim)); 
        % imagesc(effectsallTEST(:,:,z,stim)); 
        caxis([-2 2]); 
        
        %colormap('Gray')
        
    end
    colorbar
end
      

titles={'Control';' ';' ';' ';'redundants';' ';' ';' ';' ';'DEVIANT'};
redundantnumberA=0; % this is the redundant in the chain that you will plot/analyse; if it's 0, then average over all redundants
include_cutoff=1; %number of standard deviations above baseline for inclusion Using this, then there won't be any inhibited cells. 
exclude_outlier=100; %get rid of the cell if any responses are this big

%if you want to just plot averages accross all trials, make stringent test
%%% equal to 0. if you want to do a hard test of clusters,make it 1... use
%%% 1 only if you are actually testing teh reality of a subcluster


