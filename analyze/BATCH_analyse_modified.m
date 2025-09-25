%%% Working memo for John: 
% There are inhibited part and excited part. Let's get those out, without
% touching the normal part. I still don't understand the build and test.
% And I don't need to work on that, too.
% Do I have a new story on the rndt? Yes, for the nPE. Let's start. 
% Identify what needs to be done
% Scatter not found. Aggregated found at fig 1300
% Imagesc found at 1600, but only with responsive cells. Maybe it's
% overwhittened. 
% Notice that the normalization is fully changed in this code
% mnn is not minimum of the baseline, but average.


%%% do stringent test and pre/post for supplemental figure, to show
%%% stability of conclusions. VIPs=no mod. SSTs=slight DD, strong SSA.
%%% PYRs=strong DD, strong SSA jojo
% clcl
close all
clear 
for file_load=1

% files = {'regular_oddball_data_JMENG1'};
%files = {'omission_oddball_data_JMENG1'};
files = {'CA_workspace_SORTEDrescore_031521_vmmn_PYR_maus1_001'};

% files = {'CA_workspace_SORTEDrescore_032321_vmmn_PYR_maus2_002'};
% files = {'CA_workspace_SORTEDrescore_inclloco_102020_vmm_vglut_001'};
% files = {'CA_workspace_SORTEDrescore_inclloco_102020_vmm_vglut_002'};
% files = {'CA_workspace_SORTEDrescore_inclloco_102020_vmm_vglut_003'};
% files = {'CA_workspace_SORTEDrescore_inclloco_102020_vmm_vglut_004'};
end

limtrials1=99; %limit the number of trials analysed for each stimulus type(control/redundandt/deviant). IF set high, then it normalizes to the smallest number of trials in any condition (cont, rdnt, dev)
baseBEGIN=1; baseEND=12; %what you call the baseline
base_vec=baseBEGIN:baseEND;
% BEGINo=14;
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


for stringentTEST=0;
    if stringentTEST==0;
        %%%%%%%%%figures out relative responses to different stimuli
        simpleeffects=[]; simpleeffectsPRE=[];simpleeffectsPOST=[];c1=1; c2=1; cols=[]; newindices=[];notindices=[];
        for stim=1:2
            for c=1:size(effectsallSTD,1)
                if redundantnumberA==0; 
                	r=mean(mean(effectsallSTD(c,BEGINo:ENDo,2:9,stim)));
                else
                    r=mean(effectsallSTD(c,BEGINo:ENDo,redundantnumberA+1,stim));
                end
                cn=mean(effectsallSTD(c,BEGINo:ENDo,1,stim));
                d=mean(effectsallSTD(c,BEGINo:ENDo,10,stim));
                if redundantnumberA==0 
                	rP=mean(mean(effectsallBUILD(c,BEGINo:ENDo,2:9,stim)));
                else
                    rP=mean(effectsallBUILD(c,BEGINo:ENDo,redundantnumberA+1,stim));
                end
                cnP=mean(effectsallBUILD(c,BEGINo:ENDo,1,stim));
                dP=mean(effectsallBUILD(c,BEGINo:ENDo,10,stim));
                if redundantnumberA==0 
                	rT=mean(mean(effectsallTEST(c,BEGINo:ENDo,2:9,stim)));
                else
                    rT=mean(effectsallTEST(c,BEGINo:ENDo,redundantnumberA+1,stim));
                end
                cnT=mean(effectsallTEST(c,BEGINo:ENDo,1,stim));
                dT=mean(effectsallTEST(c,BEGINo:ENDo,10,stim));
                % Get the responsive cells
                % if and(or(or(r>include_cutoff,cn>include_cutoff),d>include_cutoff),and(and(r<exclude_outlier,cn<exclude_outlier),d<exclude_outlier))    %??
                if 1==1
                    simpleeffects(c1,1)=(cn);
                    simpleeffects(c1,2)=(r);
                    simpleeffects(c1,3)=(d);
                    simpleeffectsPRE(c1,1)=(cnP);
                    simpleeffectsPRE(c1,2)=(rP);
                    simpleeffectsPRE(c1,3)=(dP);
                    simpleeffectsPOST(c1,1)=(cnT);
                    simpleeffectsPOST(c1,2)=(rT);
                    simpleeffectsPOST(c1,3)=(dT);
                    newindices(c1,1)=c;
                    newindices(c1,2)=stim;
                    c1=1+c1;
                else
                    notindices(c2,1)=c;
                    notindices(c2,2)=stim;
                    c2=1+c2;
                    
                end
                
            end
        end
        for makingcolors=1;
            for c=1:(c1-1)
                s=simpleeffects(c,1)./max(simpleeffects(c,:));
                cols(c,2)=s;
                s=simpleeffects(c,3)./max(simpleeffects(c,:));
                cols(c,1)=s;
                s=simpleeffects(c,2)./max(simpleeffects(c,:));
                cols(c,3)=s;
            end
            for c=1:(c1-1);
                a=cols(c,:)-min(cols(c,:));
                cols(c,:)=cols(c,:)./max(cols(c,:));
            end
        end
        
    else
        simpleeffects=[]; simpleeffectsPRE=[];simpleeffectsPOST=[];c1=1; c2=1; cols=[]; newindices=[];notindices=[];
        for stim=1:2
            for c=1:size(effectsallBUILD,1)
                if redundantnumberA==0; 
                	r=mean(mean(effectsallBUILD(c,BEGINo:ENDo,2:9,stim)));
                else
                    r=mean(effectsallBUILD(c,BEGINo:ENDo,redundantnumberA+1,stim));
                end
                cn=mean(effectsallBUILD(c,BEGINo:ENDo,1,stim));
                d=mean(effectsallBUILD(c,BEGINo:ENDo,10,stim));
                if redundantnumberA==0; 
                	rT=mean(mean(effectsallTEST(c,BEGINo:ENDo,2:9,stim)));
                else
                    rT=mean(effectsallTEST(c,BEGINo:ENDo,redundantnumberA+1,stim));
                end
                cnT=mean(effectsallTEST(c,BEGINo:ENDo,1,stim));
                dT=mean(effectsallTEST(c,BEGINo:ENDo,10,stim));
                if and(or(or(r>include_cutoff,cn>include_cutoff),d>include_cutoff),and(and(r<exclude_outlier,cn<exclude_outlier),d<exclude_outlier))
                    simpleeffects(c1,1)=(cn);
                    simpleeffects(c1,2)=(r);
                    simpleeffects(c1,3)=(d);
                    simpleeffectsPRE(c1,1)=(cn);
                    simpleeffectsPRE(c1,2)=(r);
                    simpleeffectsPRE(c1,3)=(d);
                    simpleeffectsPOST(c1,1)=(cnT);
                    simpleeffectsPOST(c1,2)=(rT);
                    simpleeffectsPOST(c1,3)=(dT);
                    newindices(c1,1)=c;
                    newindices(c1,2)=stim;
                    c1=1+c1;
                else
                    notindices(c2,1)=c;
                    notindices(c2,2)=stim;
                    c2=1+c2;
                    
                end
                
            end
        end
        for makingcolors=1;
            for c=1:(c1-1)
                s=simpleeffects(c,1)./max(simpleeffects(c,:));
                cols(c,2)=s;
                s=simpleeffects(c,3)./max(simpleeffects(c,:));
                cols(c,1)=s;
                s=simpleeffects(c,2)./max(simpleeffects(c,:));
                cols(c,3)=s;
            end
            for c=1:(c1-1);
                a=cols(c,:)-min(cols(c,:));
                cols(c,:)=cols(c,:)./max(cols(c,:));
            end
        end
    end
end
% Till now are all necessary. 
countcellstotal=0;
% Count the responsive cells
if countcellstotal==1

% for countcellstotal=1  % for only focusing on one response per cell
    numcells=size(unique(newindices(:,1)),1);
    disp(strcat('all cells=',num2str(size(effectsall,1)),'; responsive=',num2str(numcells),'; percent=',num2str(numcells/size(effectsall,1))));
    ex1=zeros(size(simpleeffects,1),1);
    for c=1:size(effectsall,1)
        if sum(newindices(:,1)==c)>1;
            clear tmp1 tmp2; cnt1=1;
            for c1=1:size(newindices,1)
                if newindices(c1,1)==c;
                    tmp1(cnt1,1)=mean(simpleeffects(c1,[1 2 3]));
                    tmp2(cnt1,1)=c1;
                    cnt1=1+cnt1;
                end
            end
            if tmp1(1,1)>tmp1(2,1);
                ex1(tmp2(2,1),1)=1;
            else
                ex1(tmp2(1,1),1)=1;
            end
        end
    end
    for onlylookatonePERcell=1;
        if onlylookatonePERcell==1;
            simpleeffectsPRE=simpleeffectsPRE(ex1==0,:);
            simpleeffects=simpleeffects(ex1==0,:);
            simpleeffectsPOST=simpleeffectsPOST(ex1==0,:);
            newindices=newindices(ex1==0,:);
        end
    end
    for makingcolors=1;
    c1=size(simpleeffects,1)+1; clear cols
    for c=1:(c1-1)
        s=simpleeffects(c,1)./max(simpleeffects(c,:));
        cols(c,2)=s;
        s=simpleeffects(c,3)./max(simpleeffects(c,:));
        cols(c,1)=s;
        s=simpleeffects(c,2)./max(simpleeffects(c,:));
        cols(c,3)=s;
        
    end
    for c=1:(c1-1);
        a=cols(c,:)-min(cols(c,:));
        cols(c,:)=cols(c,:)./max(cols(c,:));
        cols(c,2)=0;cols(c,3)=0;%cols(c,1)=0;
    end
end

end

% John. Select what is called redundant.

redundantnumber=0; %if you want to plot (but not "select by") a redudnant later in the stream, use this line
for plottingeverything=1; 
%%%%for selecting all or subset of cells for subesequent plots
%%%for selecting all or subset of cells for subesequent plots
plotanti=0; %if =1, then plot everything outside of the selected box
selections=1; %if the cells you want to select are scattered more than a single box can capture, make this 2 or 3 and select multiple regions.
for plotpre=1;
    %close all force
    xcut=std((.5*simpleeffects(:,1)+.5*simpleeffects(:,3))-simpleeffects(:,2))*4; ycut=std(simpleeffects(:,3)-simpleeffects(:,1))*4;
    xmn=mean((.5*simpleeffects(:,1)+.5*simpleeffects(:,3))-simpleeffects(:,2)); ymn=mean(simpleeffects(:,3)-simpleeffects(:,1));
    
    % pind=[];
    % for selecto=1:selections
    %     figure(100+selecto);   % I only see selecto 1.  
    %     scatter((.5*simpleeffects(:,1)+.5*simpleeffects(:,3))-simpleeffects(:,2),simpleeffects(:,3)-simpleeffects(:,1),50,cols,'fill');
    % 
    %     if selecto>1;
    %         hold on;
    %         scatter([x11(1) x11(1) x11(2) x11(2)],[y11(1) y11(2) y11(1) y11(2)]);
    %     end
    %     ylabel('Deviance Detection','FontSize',12,'FontWeight','bold'); xlabel('Stim specific adaptation','FontSize',12,'FontWeight','bold');
    %     title('pick region to plot; upper left, then lower right');
    %     [x,y]=ginput(2); close
    %     if selecto==1;
    %         x11=x; y11=y;
    %     end
    %     if selecto==2;
    %         x22=x; y22=y;
    %     end
    % 
    %     for c=1:size(simpleeffects,1);
    %         if and(((.5*simpleeffects(c,1)+.5*simpleeffects(c,3))-simpleeffects(c,2))>x(1),((.5*simpleeffects(c,1)+.5*simpleeffects(c,3))-simpleeffects(c,2))<x(2));
    %             if and((simpleeffects(c,3)-simpleeffects(c,1))<y(1),(simpleeffects(c,3)-simpleeffects(c,1))>y(2))
    %                 if sum(c==pind)==0;
    %                     pind=vertcat(pind,c);
    %                 end
    %             end
    %         end
    %     end
    % end
    pind=1:size(simpleeffects,1);
    pind=pind';
    
    if plotanti==1;
        ppp=[];
        pp=1:size(simpleeffects,1); pp=pp';
        for c=1:size(simpleeffects,1);
            if sum(pp(c)==pind)==0;
                ppp=vertcat(ppp,c);
            end
        end
        pind=ppp;
    end
    numcells=size(pind,1);
    tmpeff=[];
    for c=1:numcells;
        tmpeff(c,:,:)=effectsallBUILD(newindices(pind(c),1),:,:,(newindices(pind(c),2)));
        
        for t=1:size(effectsallBUILD,2)
            if redundantnumber>0;
                mn=mean(mean(effectsallBUILD(newindices(pind(c),1),t,[redundantnumber+1 1 10],newindices(pind(c),2)),3),4);
                a1(c,t)=mean(effectsallBUILD(newindices(pind(c),1),t,redundantnumber+1,newindices(pind(c),2)),4)-mn;
            else
                a1tmp=mean(mean(effectsallBUILD(newindices(pind(c),1),t,2:9,newindices(pind(c),2)),3),4);
                a2tmp=mean(effectsallBUILD(newindices(pind(c),1),t,1,newindices(pind(c),2)),4);
                a3tmp=mean(effectsallBUILD(newindices(pind(c),1),t,10,newindices(pind(c),2)),4);
                mn=(a1tmp+a2tmp+a3tmp)/3;
                a1(c,t)=mean(mean(effectsallBUILD(newindices(pind(c),1),t,2:9,newindices(pind(c),2)),3),4)-mn;
            end
            a2(c,t)=mean(effectsallBUILD(newindices(pind(c),1),t,1,newindices(pind(c),2)),4)-mn;
            a3(c,t)=mean(effectsallBUILD(newindices(pind(c),1),t,10,newindices(pind(c),2)),4)-mn;
        end
    end
    x=timeaxis;
    y1=mean(mean(mean(tmpeff(:,:,redundantnumber+1),3),4),1); 
    if redundantnumber==0;y1=mean(mean(mean(tmpeff(:,:,2:9),3),4),1); end  
    y1a=y1-min(y1(1,1:BEGINo-1));
    err1=std(a1,0,1)./sqrt(numcells);
    
    y2=mean(mean(mean(tmpeff(:,:,1),3),4),1);y2a=y2-min(y2(1,1:BEGINo-1));
    err2=std(a2,0,1)./sqrt(numcells);
    
    y3=mean(mean(mean(tmpeff(:,:,10),3),4),1); y3a=y3-min(y3(1,1:BEGINo-1));
    err3=std(a3,0,1)./sqrt(numcells);
    
    
    if numcells>1;
        figure(200); 
        shadedErrorBar(x,y1a,err1,'lineProps','b'); 
        hold on;shadedErrorBar(x,y2a,err1,'lineProps','k');
        hold on;shadedErrorBar(x,y3a,err1,'lineProps','r');
        xlim([-.2 1]); 
        title ('Pre trials', 'FontSize', 16, 'FontWeight', 'bold');
        ylabel('z-scores','FontSize',16,'FontWeight','bold'); 
        xlabel('time (sec)','FontSize',16,'FontWeight','bold'); 
        set(gca,'FontSize',16,'FontWeight','bold');
        set(gcf,'Color','w');
        
        yy=ylim;
        % figure(200+1);    % The errorbar_groups contains the figID.
        errorbar_groups([mean(mean(y2a(:,BEGINo:ENDo))) mean(mean(y1a(:,BEGINo:ENDo))) mean(mean(y3a(:,BEGINo:ENDo)))]',[mean(mean(err2(:,BEGINo:ENDo))) mean(mean(err1(:,BEGINo:ENDo))) mean(mean(err3(:,BEGINo:ENDo)))]', ...
            'bar_colors',[.5 .5 .5; 0 0 1; 1 0 0],'FigID',201);
         % john modified
         % figure;

         
        ylabel('z-scores','FontSize',16,'FontWeight','bold'); xlabel('time (sec)','FontSize',16,'FontWeight','bold'); set(gca,'FontSize',16,'FontWeight','bold');
        set(gcf,'Color','w');
        for statstat=1;
            r1=mean(mean(tmpeff(:,BEGINo:ENDo,redundantnumber+1),3),2)-mean(mean(tmpeff(:,BEGINo-9:BEGINo-1,redundantnumber+1),3),2);
            if redundantnumber==0;
                 r1=mean(mean(tmpeff(:,BEGINo:ENDo,2:9),3),2)-mean(mean(tmpeff(:,BEGINo-9:BEGINo-1,2:9),3),2);
            end
            c1=mean(mean(tmpeff(:,BEGINo:ENDo,1),3),2)-mean(mean(tmpeff(:,BEGINo-9:BEGINo-1,1),3),2);
            d1=mean(mean(tmpeff(:,BEGINo:ENDo,10),3),2)-mean(mean(tmpeff(:,BEGINo-9:BEGINo-1,10),3),2);
            
            [h,p1,ci,stats1]=ttest(r1,c1);
            [h,p2,ci,stats2]=ttest(d1,c1);
        end
        title(strcat('T_S_S_A','=',num2str(round(stats1.tstat*100)/100),', p=',num2str(round(p1*1000)/1000),'; T_D_D','=',num2str(round(stats2.tstat*100)/100),', p=',num2str(round(p2*1000)/1000)));
        title('pre trials');
    else
        figure(400);
        plot(x,y1,'b','LineWidth',3); 
        hold on;plot(x,y2,'k','LineWidth',3); 
        hold on;plot(x,y3,'r','LineWidth',3);
        yy=ylim;
        xlim([-.2 1]); 
        title ('Average of cell responses: pre trials', 'FontSize', 16, 'FontWeight', 'bold');
    end
    
    for manyplotthing=1;
        figure(500);
        for cnd=1:10;
            y1=mean(tmpeff(:,:,cnd),1);
            y1=y1-min(y1(1,baseBEGIN:baseEND));
            if cnd==1; err1=std(a1,0,1)./sqrt(numcells);
                err=std(a1,0,1)./sqrt(numcells);
                
            elseif cnd==10;
                err=std(a3,0,1)./sqrt(numcells);
            else
                err=std(a2,0,1)./sqrt(numcells);
            end
            subplot(1,10,cnd); shadedErrorBar(x,y1,err,'lineProps','k'); set(gcf,'Color','w');
            if cnd==1;
                ylabel('z-scores','FontSize',16,'FontWeight','bold'); xlabel('time (sec)','FontSize',16,'FontWeight','bold'); set(gca,'FontSize',16,'FontWeight','bold');
                title(titles{cnd}, 'FontSize', 18, 'FontWeight', 'bold');
            else
                title(titles{cnd}, 'FontSize', 18, 'FontWeight', 'bold'); yticks([]);set(gca,'FontSize',16,'FontWeight','bold');
            end
        end
        for fixaxes=1
            
            ymin1=1; ymax1=-1;
            for cnd=1:10
                subplot(1,10,cnd); yy=ylim;
                if yy(1)<ymin1; ymin1=yy(1); end
                if yy(2)>ymax1; ymax1=yy(2); end
                
            end
            for cnd=1:10
                subplot(1,10,cnd); ylim([ymin1 ymax1]);
            end
            
        end
    end
    
    for testoverlapppingcelss=1;
        for pp=1:size(pind,1);
            aaa(pp)=newindices(pind(pp),1);
            bbb(pp)=newindices(pind(pp),2);
        end
        zerp=zeros(size(effectsall,1),2);
        for pp=1:size(aaa,2);
            zerp(aaa(pp),1)=1+zerp(aaa(pp));
            
            zerp(aaa(pp),2)=bbb(pp);
        end
        onestimDD=sum(zerp(:,1)==1);
        twostimDD=sum(zerp(:,1)==2);
        weird=zeros(size(zerp,1),1);ffllpp=[2 1];
        for ppp=1:size(zerp,1);
            if zerp(ppp,1)==1;
                if effectsall(ppp,BEGINo:ENDo,1,ffllpp(zerp(ppp,2)))>include_cutoff;
                    weird(ppp)=1;
                end
            end
        end
        figure(600); pie([onestimDD twostimDD]);   title('proportion of DDs to one stim vs both');
    end
end

for PLOTpost=1;
    
    tmpeff=[];
    for c=1:numcells;
        tmpeff(c,:,:)=effectsallTEST(newindices(pind(c),1),:,:,(newindices(pind(c),2)));
        for t=1:size(effectsall,2)
            if redundantnumber>0;
                mn=mean(mean(effectsallTEST(newindices(pind(c),1),t,[redundantnumber+1 1 10],newindices(pind(c),2)),3),4);
                a1(c,t)=mean(effectsallTEST(newindices(pind(c),1),t,redundantnumber+1,newindices(pind(c),2)),4)-mn;
            else
                 a1tmp=mean(mean(effectsallTEST(newindices(pind(c),1),t,2:9,newindices(pind(c),2)),3),4);
                a2tmp=mean(effectsallTEST(newindices(pind(c),1),t,1,newindices(pind(c),2)),4);
                a3tmp=mean(effectsallTEST(newindices(pind(c),1),t,10,newindices(pind(c),2)),4);
                mn=(a1tmp+a2tmp+a3tmp)/3;
                a1(c,t)=mean(mean(effectsallTEST(newindices(pind(c),1),t,2:9,newindices(pind(c),2)),3),4)-mn;
            end
            a2(c,t)=mean(effectsallTEST(newindices(pind(c),1),t,1,newindices(pind(c),2)),4)-mn;
            a3(c,t)=mean(effectsallTEST(newindices(pind(c),1),t,10,newindices(pind(c),2)),4)-mn;
        end
    end
    
    x=timeaxis;
    y1=mean(mean(mean(tmpeff(:,:,redundantnumber+1),3),4),1); 
    
    if redundantnumber==0; y1=mean(mean(mean(tmpeff(:,:,2:9),3),4),1); end; 
    y1a=y1-min(y1(1,1:BEGINo-1));
    err1=std(a1,0,1)./sqrt(numcells);
    
    y2=mean(mean(mean(tmpeff(:,:,1),3),4),1);y2a=y2-min(y2(1,1:BEGINo-1));
    err2=std(a2,0,1)./sqrt(numcells);
    
    y3=mean(mean(mean(tmpeff(:,:,10),3),4),1); y3a=y3-min(y3(1,1:BEGINo-1));
    err3=std(a3,0,1)./sqrt(numcells);
    
    
    if numcells>1;
        figure(700); 
        shadedErrorBar(x,y1a,err1,'lineProps','b'); 
        hold on;shadedErrorBar(x,y2a,err1,'lineProps','k');hold on;shadedErrorBar(x,y3a,err1,'lineProps','r');
        xlim([-.2 1]); title ('Post trials', 'FontSize', 16, 'FontWeight', 'bold');ylabel('z-scores','FontSize',16,'FontWeight','bold'); xlabel('time (sec)','FontSize',16,'FontWeight','bold'); set(gca,'FontSize',16,'FontWeight','bold');
        set(gcf,'Color','w');
        
        yy=ylim;
        % figure(701);
        errorbar_groups([mean(mean(y2a(:,BEGINo:ENDo))) mean(mean(y1a(:,BEGINo:ENDo))) mean(mean(y3a(:,BEGINo:ENDo)))]',[mean(mean(err2(:,BEGINo:ENDo))) mean(mean(err1(:,BEGINo:ENDo))) mean(mean(err3(:,BEGINo:ENDo)))]', ...
            'bar_colors',[.5 .5 .5; 0 0 1; 1 0 0],'FigID',701);
        
        ylabel('z-scores','FontSize',16,'FontWeight','bold'); xlabel('time (sec)','FontSize',16,'FontWeight','bold'); set(gca,'FontSize',16,'FontWeight','bold');
        set(gcf,'Color','w');
        for statstat=1;
            r1=mean(mean(tmpeff(:,BEGINo:ENDo,redundantnumber+1),3),2);
            if redundantnumber==0;
                 r1=mean(mean(tmpeff(:,BEGINo:ENDo,2:9),3),2)-mean(mean(tmpeff(:,BEGINo-9:BEGINo-1,2:9),3),2);
            end
            c1=mean(mean(tmpeff(:,BEGINo:ENDo,1),3),2);
            d1=mean(mean(tmpeff(:,BEGINo:ENDo,10),3),2);
            
            [h,p1,ci,stats1]=ttest(r1,c1);
            [h,p2,ci,stats2]=ttest(d1,c1);
            clear simpleeffectsTEST
            simpleeffectsTEST(:,1)=c1;simpleeffectsTEST(:,2)=r1;simpleeffectsTEST(:,3)=d1;
        end
        title(strcat('T_S_S_A','=',num2str(round(stats1.tstat*100)/100),', p=',num2str(round(p1*1000)/1000),'; T_D_D','=',num2str(round(stats2.tstat*100)/100),', p=',num2str(round(p2*1000)/1000)));
        
    else
        figure(702); plot(x,y1,'b','LineWidth',3); hold on;plot(x,y2,'k','LineWidth',3); hold on;plot(x,y3,'r','LineWidth',3);
        xlim([-.2 1]); title ('Mismatch Z Scores', 'FontSize', 16, 'FontWeight', 'bold');
        
        
    end
    tmpeff_post=tmpeff; 
    for manyplotthing=1;
        figure(800);
        for cnd=1:10;
            y1=mean(tmpeff(:,:,cnd),1);
            y1=y1-min(y1(1,baseBEGIN:baseEND));
            if cnd==1; err1=std(a1,0,1)./sqrt(numcells);
                err=std(a1,0,1)./sqrt(numcells);
                
            elseif cnd==10;
                err=std(a3,0,1)./sqrt(numcells);
            else
                err=std(a2,0,1)./sqrt(numcells);
            end
            subplot(1,10,cnd); shadedErrorBar(x,y1,err,'lineProps','k'); set(gcf,'Color','w');
            if cnd==1;
                ylabel('z-scores','FontSize',16,'FontWeight','bold'); xlabel('time (sec)','FontSize',16,'FontWeight','bold'); set(gca,'FontSize',16,'FontWeight','bold');
                title(titles{cnd}, 'FontSize', 18, 'FontWeight', 'bold');
            else
                title(titles{cnd}, 'FontSize', 18, 'FontWeight', 'bold'); yticks([]);set(gca,'FontSize',16,'FontWeight','bold');
            end
        end
        for fixaxes=1
            
            ymin1=1; ymax1=-1;
            for cnd=1:10
                subplot(1,10,cnd); yy=ylim;
                if yy(1)<ymin1; ymin1=yy(1); end
                if yy(2)>ymax1; ymax1=yy(2); end
                
            end
            for cnd=1:10
                subplot(1,10,cnd); ylim([ymin1 ymax1]);
            end
            
        end
    end
    
     
   
    
end

%%%plot average of both


%%%%fo plotting ongoing activity, i.e. not stimulus locked
for plot_ongoing=1;
    dosubsetting=0;
    if dosubsetting==1;
        
    tmpeff2=[];
    actcont2=actcont;
    actcont2(:,:,2)=actcont;
    
    for c=1:numcells
        tmpeff2(c,:)=actcont2(newindices(pind(c),1),:,(newindices(pind(c),2)));
                   tmpeff2(c,1:10)=tmpeff2(c,1:10)-tmpeff2(c,1);
                   tmpeff2(c,11:20)=tmpeff2(c,11:20)-tmpeff2(c,1);
                  tmpeff2(c,21:30)=tmpeff2(c,21:30)-tmpeff2(c,1);
    end
    
    x=mean(tmpeff2); y=std(tmpeff2/sqrt(size(tmpeff2,1)));
    
    figure(900);
    subplot(1,2,1); shadedErrorBar(1:10,x(1:10),y(1:10)); title({'Ongoing activity of cells in cluster';'CONTROL                        ODDBALL1                        ODDBALL2'}); hold on
    subplot(1,2,1); shadedErrorBar(12:21,x(11:20),y(11:20)); hold on
    subplot(1,2,1); shadedErrorBar(23:32,x(21:30),y(21:30)); hold on
    
    tmpeff2=[];
    actcont2=actcont;
    actcont2(:,:,2)=actcont;
    
    for c=1:size(notindices,1)
        tmpeff2(c,:)=actcont2(notindices(c,1),:,(notindices(c,2)));
                     tmpeff2(c,1:10)=tmpeff2(c,1:10)-tmpeff2(c,1);
                   tmpeff2(c,11:20)=tmpeff2(c,11:20)-tmpeff2(c,1);
                   tmpeff2(c,21:30)=tmpeff2(c,21:30)-tmpeff2(c,2);
        
    end
    
    x=mean(tmpeff2); y=std(tmpeff2/sqrt(size(tmpeff2,1)));
    
    subplot(1,2,2); shadedErrorBar(1:10,x(1:10),y(1:10)); title({'Ongoing activity non stim-driven cells';'CONTROL                        ODDBALL1                        ODDBALL2'}); hold on
    subplot(1,2,2); shadedErrorBar(12:21,x(11:20),y(11:20)); hold on
    subplot(1,2,2); shadedErrorBar(23:32,x(21:30),y(21:30)); hold on
    
    y1=axis(subplot(1,2,1)); y2=axis(subplot(1,2,2));
    yn=min([y1(3) y2(3)]); yx=max([y1(4) y2(4)]);
    for z=1:2; subplot(1,2,z); ylim([yn yx]); end
    else
        tmpeff2=[];
    actcont2=actcont;
    actcont2(:,:,2)=actcont;
    c1=0;
    sponcnt=[]; sponMMN=[];
    for c=1:size(actcont2,1);
        if sum(isnan(actcont2(c,:,1)))==0
            c1=1+c1; 
                tmpeff2(c1,:)=actcont2(c,:,1);
                kb=tmpeff2(c1,1);
                   tmpeff2(c1,1:10)=tmpeff2(c1,1:10)-kb;
                   tmpeff2(c1,11:20)=tmpeff2(c1,11:20)-kb;
                  tmpeff2(c1,21:30)=tmpeff2(c1,21:30)-kb;
                  sponcnt(c1,1)=mean(tmpeff2(c1,1:10),2);
                  sponMMN(c1,1)=mean(tmpeff2(c1,11:end),2);
        end
    end
    
    x=mean(tmpeff2); y=std(tmpeff2/sqrt(size(tmpeff2,1)));

    % John: Somehow the following is commented out
% figure;  
%      shadedErrorBar(1:10,x(1:10),y(1:10)); title({'Ongoing activity of cells';'CONTROL                        ODDBALL1                        ODDBALL2'}); hold on
%     shadedErrorBar(12:21,x(11:20),y(11:20)); hold on
%      shadedErrorBar(23:32,x(21:30),y(21:30)); hold on
%      
%     figure;scatter(sponcnt,sponMMN,50,'k','fill'); 
%     title ('average spontaneous activity', 'FontSize', 20, 'FontWeight', 'bold');
%     ylabel('during mismatch','FontSize',18,'FontWeight','bold'); xlabel('during control run','FontSize',18,'FontWeight','bold'); set(gca,'FontSize',16,'FontWeight','bold');set(gcf,'Color','w'); hold on
%    plot(min(min(sponcnt),min(sponMMN)):.1:max(max(sponcnt),max(sponMMN)),min(min(sponcnt),min(sponMMN)):.1:max(max(sponcnt),max(sponMMN)),'k');
%     axis([-5 5 -5 5]); 
    end
end

for plotbOLTH=1;
    effectsallBOTH=effectsallSTD; %(effectsallBUILD+effectsallTEST)/2;  ???

    DEBUGeffect1=effectsallBOTH;
    tmpeff=[];
    for c=1:numcells;
        tmpeff(c,:,:)=effectsallBOTH(newindices(pind(c),1),:,:,(newindices(pind(c),2)));
        for t=1:size(effectsallBOTH,2)
            if redundantnumber>0;
                mn=mean(mean(effectsallBOTH(newindices(pind(c),1),t,[redundantnumber+1 1 10],newindices(pind(c),2)),3),4);
                a1(c,t)=mean(effectsallBOTH(newindices(pind(c),1),t,redundantnumber+1,newindices(pind(c),2)),4)-mn;
            else
                
                a1tmp=mean(mean(effectsallBOTH(newindices(pind(c),1),t,2:9,newindices(pind(c),2)),3),4);
                a2tmp=mean(effectsallBOTH(newindices(pind(c),1),t,1,newindices(pind(c),2)),4);
                a3tmp=mean(effectsallBOTH(newindices(pind(c),1),t,10,newindices(pind(c),2)),4);
                mn=(a1tmp+a2tmp+a3tmp)/3;
                a1(c,t)=mean(mean(effectsallBOTH(newindices(pind(c),1),t,2:9,newindices(pind(c),2)),3),4)-mn;
            end
            a2(c,t)=mean(effectsallBOTH(newindices(pind(c),1),t,1,newindices(pind(c),2)),4)-mn;
            a3(c,t)=mean(effectsallBOTH(newindices(pind(c),1),t,10,newindices(pind(c),2)),4)-mn;
        end
    end
   

    x=timeaxis;
    y1=mean(mean(mean(tmpeff(:,:,redundantnumber+1),3),4),1);
    if redundantnumber==0; y1=mean(mean(mean(tmpeff(:,:,2:9),3),4),1); end;
    y1a=y1-mean(y1(1,BEGINo-9:BEGINo-1));
    err1=std(a1,0,1)./sqrt(numcells);
    
    y2=mean(mean(mean(tmpeff(:,:,1),3),4),1);y2a=y2-mean(y2(1,BEGINo-9:BEGINo-1));
    err2=std(a2,0,1)./sqrt(numcells);
    
    y3=mean(mean(mean(tmpeff(:,:,10),3),4),1); y3a=y3-mean(y3(1,BEGINo-9:BEGINo-1));
    err3=std(a3,0,1)./sqrt(numcells);
    
       % mm=min(horzcat(y1a(1,1:BEGINo-1),y2a(1,1:BEGINo-1),y3a(1,1:BEGINo-1)));
    mm=0;
     y1a=y1a-mm;y2a=y2a-mm;y3a=y3a-mm;
    
    if numcells>1;
        figure(1000); shadedErrorBar(x,y1a,err1,'lineProps','b'); hold on;shadedErrorBar(x,y2a,err2,'lineProps','k');hold on;shadedErrorBar(x,y3a,err3,'lineProps','r');
        xlim([-.2 1]); title ('Average of cell responses', 'FontSize', 16, 'FontWeight', 'bold');ylabel('z-scores','FontSize',16,'FontWeight','bold'); xlabel('time (sec)','FontSize',16,'FontWeight','bold'); set(gca,'FontSize',16,'FontWeight','bold');
        set(gcf,'Color','w');
        
        yy=ylim;
        % figure(1001); 
        errorbar_groups([mean(mean(y2a(:,BEGINo:ENDo))) mean(mean(y1a(:,BEGINo:ENDo))) mean(mean(y3a(:,BEGINo:ENDo)))]',[mean(mean(err2(:,BEGINo:ENDo))) mean(mean(err1(:,BEGINo:ENDo))) mean(mean(err3(:,BEGINo:ENDo)))]', ...
            'bar_colors',[.5 .5 .5; 0 0 1; 1 0 0],'FigID',1001);
        ylabel('z-scores','FontSize',16,'FontWeight','bold'); xlabel('time (sec)','FontSize',16,'FontWeight','bold'); set(gca,'FontSize',16,'FontWeight','bold');
        set(gcf,'Color','w');
        for statstat=1;
            r1=mean(mean(tmpeff(:,BEGINo:ENDo,redundantnumber+1),3),2)-mean(mean(tmpeff(:,BEGINo-9:BEGINo-1,redundantnumber+1),3),2);
            if redundantnumber==0;
                 r1=mean(mean(tmpeff(:,BEGINo:ENDo,2:9),3),2)-mean(mean(tmpeff(:,BEGINo-9:BEGINo-1,2:9),3),2);
            end
            c1=mean(mean(tmpeff(:,BEGINo:ENDo,1),3),2)-mean(mean(tmpeff(:,BEGINo-9:BEGINo-1,1),3),2);
            d1=mean(mean(tmpeff(:,BEGINo:ENDo,10),3),2)-mean(mean(tmpeff(:,BEGINo-9:BEGINo-1,10),3),2);
            
            [h,p1,ci,stats1]=ttest(r1,c1);
            [h,p2,ci,stats2]=ttest(d1,c1);
        end
        title(strcat('T_S_S_A','=',num2str(round(stats1.tstat*100)/100),', p=',num2str(round(p1*1000)/1000),'; T_D_D','=',num2str(round(stats2.tstat*100)/100),', p=',num2str(round(p2*1000)/1000)));
    else
        figure(1002); plot(x,y1,'b','LineWidth',3); hold on;plot(x,y2,'k','LineWidth',3); hold on;plot(x,y3,'r','LineWidth',3);yy=ylim;
        xlim([-.2 1]); title ('Average of cell responses', 'FontSize', 16, 'FontWeight', 'bold');
    end
    
    for manyplotthing=1;
        figure(1100);
        for cnd=1:10;
            y1=mean(tmpeff(:,:,cnd),1);
            y1=y1-mean(y1(1,1:BEGINo-1));
            if cnd==1; err1=std(a2,0,1)./sqrt(numcells);
                err=std(a2,0,1)./sqrt(numcells);
                
            elseif cnd==10
                err=std(a3,0,1)./sqrt(numcells);
            else
                err=std(a1,0,1)./sqrt(numcells);
            end
            subplot(1,10,cnd); shadedErrorBar(x,y1,err,'lineProps','k'); set(gcf,'Color','w'); xlim([-.2 .8]);
            if cnd==1
                ylabel('z-scores','FontSize',16,'FontWeight','bold'); xlabel('time (sec)','FontSize',16,'FontWeight','bold'); set(gca,'FontSize',16,'FontWeight','bold');
                title(titles{cnd}, 'FontSize', 18, 'FontWeight', 'bold');
            else
                title(titles{cnd}, 'FontSize', 18, 'FontWeight', 'bold'); yticks([]);set(gca,'FontSize',16,'FontWeight','bold');
            end
        end
        for fixaxes=1
            
            ymin1=1; ymax1=-1;
            for cnd=1:10
                subplot(1,10,cnd); yy=ylim;
                if yy(1)<ymin1; ymin1=yy(1); end
                if yy(2)>ymax1; ymax1=yy(2); end
                
            end
            for cnd=1:10
                subplot(1,10,cnd); ylim([ymin1 ymax1]);
            end
            
        end
    end
    
    for manyplotthing=2;
        
        if redundantnumber==0;
            figure(1200);
            ccnds={1;2:9;10};
            cccf=[1 5 10];
            
            colars={'k' 'b' 'r'};
            for cnd=1:3;
                y1=mean(mean(tmpeff(:,:,ccnds{cnd}),3),1);
                y1=y1-mean(y1(1,BEGINo-9:BEGINo-1));
                if cnd==1; err1=std(a2,0,1)./sqrt(numcells/2);
                    err=std(a2,0,1)./sqrt(numcells/2);
                    
                elseif cnd==3;
                    err=std(a3,0,1)./sqrt(numcells/2);
                else
                    err=std(a1,0,1)./sqrt(numcells/2);
                end
                subplot(1,3,cnd); shadedErrorBar(x,y1,err,'lineProps',colars{cnd}); set(gcf,'Color','w'); xlim([-.2 .8]);
                if cnd==1;
                    ylabel('z-scores','FontSize',16,'FontWeight','bold'); xlabel('time (sec)','FontSize',16,'FontWeight','bold'); set(gca,'FontSize',16,'FontWeight','bold');
                    title(titles{cccf(cnd)}, 'FontSize', 18, 'FontWeight', 'bold');
                    cont1=y1./max(y1); cont1err=err; y11=y1;
                else
                    if cnd==2;
                        title(titles{5}, 'FontSize', 18, 'FontWeight', 'bold'); yticks([]);set(gca,'FontSize',16,'FontWeight','bold');
                        rdnt1=y1./max(y11); rdnt1err=err;
                    else
                        title(titles{cccf(cnd)}, 'FontSize', 18, 'FontWeight', 'bold'); yticks([]);set(gca,'FontSize',16,'FontWeight','bold');
                        dev1=y1./max(y11); dev1err=err;
                        
                    end
                end
            end
            for fixaxes=1
                
                ymin1=1; ymax1=-1;
                for cnd=1:3
                    subplot(1,3,cnd); yy=ylim;
                    if yy(1)<ymin1; ymin1=yy(1); end
                    if yy(2)>ymax1; ymax1=yy(2); end
                    
                end
                for cnd=1:3
                    subplot(1,3,cnd); ylim([ymin1 ymax1]);
                end
                
            end
        else
            figure(1300);
            ccnds=[1 redundantnumber+1 10];
            % redundantnumber: what to plot as rdnt. line 770
            colars={'k' 'b' 'r'};
            for cnd=1:3
                y1=mean(mean(tmpeff(:,:,ccnds(cnd)),3),1);
                % line 831
                y1=y1-mean(y1(1,BEGINo-9:BEGINo-1));
                if cnd==1; err1=std(a1,0,1)./sqrt(numcells/2);
                    err=std(a1,0,1)./sqrt(numcells/2);
                    
                elseif cnd==3
                    err=std(a3,0,1)./sqrt(numcells/2);
                else  % cnd=2
                    err=std(a2,0,1)./sqrt(numcells/2);
                end
                subplot(1,3,cnd); shadedErrorBar(x,y1,err,'lineProps',colars{cnd}); set(gcf,'Color','w'); xlim([-.2 .8]);
                if cnd==1;
                    ylabel('z-scores','FontSize',16,'FontWeight','bold'); xlabel('time (sec)','FontSize',16,'FontWeight','bold'); set(gca,'FontSize',16,'FontWeight','bold');
                    title(titles{ccnds(cnd)}, 'FontSize', 18, 'FontWeight', 'bold');
                    cont1=y1./max(y1); cont1err=err; y11=y1;
                else
                    if cnd==2;
                        title(titles{5}, 'FontSize', 18, 'FontWeight', 'bold'); yticks([]);set(gca,'FontSize',16,'FontWeight','bold');
                        rdnt1=y1./max(y11); rdnt1err=err;
                    else
                        title(titles{ccnds(cnd)}, 'FontSize', 18, 'FontWeight', 'bold'); yticks([]);set(gca,'FontSize',16,'FontWeight','bold');
                        dev1=y1./max(y11); dev1err=err;
                    end
                end
            end
            for fixaxes=1
                
                ymin1=1; ymax1=-1;
                for cnd=1:3
                    subplot(1,3,cnd); yy=ylim;
                    if yy(1)<ymin1; ymin1=yy(1); end
                    if yy(2)>ymax1; ymax1=yy(2); end
                    
                end
                for cnd=1:3
                    subplot(1,3,cnd); ylim([ymin1 ymax1]);
                end
                
            end
        end
    end
    figure(1400); shadedErrorBar(x,rdnt1,rdnt1err,'lineProps','b'); hold on;shadedErrorBar(x,cont1,cont1err,'lineProps','k');hold on;shadedErrorBar(x,dev1,dev1err,'lineProps','r');
    xlim([-.2 1]); title ('Average of cell responses; cont normalized', 'FontSize', 16, 'FontWeight', 'bold');ylabel('z-scores','FontSize',16,'FontWeight','bold'); xlabel('time (sec)','FontSize',16,'FontWeight','bold'); set(gca,'FontSize',16,'FontWeight','bold');
    set(gcf,'Color','w');
    
    
    
    for PLOTsingtrials=1;
        
        tmpeff1=[];
        for c=1:numcells;
            tmpeff1(c,:)=singtrialsall(newindices(pind(c),1),:,(newindices(pind(c),2)));
        end
        %sz=max([size(stim0,2),size(tmpeff1,2)])
        figure(1401);
        subplot(2,1,1); plot(stim10); %xlim([1 sz]);
        subplot(2,1,2); plot(mean(tmpeff1,1)); %xlim([1 sz]);
    end
end



cols3=cols-cols;
% try 
DDs1=simpleeffectsPRE(:,3)-(mean(simpleeffectsPRE(:,1),2));DDs2=simpleeffectsPOST(:,3)-(mean(simpleeffectsPOST(:,1),2));
DDs3=simpleeffectsPRE(:,3)-(mean(simpleeffectsPRE(:,2),2));DDs4=simpleeffectsPOST(:,3)-(mean(simpleeffectsPOST(:,2),2));
sm=(DDs1>1.67)+(DDs2>1.67)+(DDs3>.1)+(DDs4>.1);%
figure(1500);  % Seems some kind of the significant test. 
ax = gca(); 
pieData = [mean(sm==4) mean(sm~=4)]; 
h = pie(ax, pieData); 
newColors = [...
1,       0, 0;   
.5,       .5,       0.5];  
ax.Colormap = newColors; 
colDD=cols-cols;
for c=1:size(sm,1)
    if sm(c,1)==4;
        colDD(c,1)=1;
    end
end
title(num2str(mean(sm~=4)));
%     for determiningcombinedPvals=1;
%         pp=1;
%         for p1=.01:.01:1;
%             pvals(pp)=chi2pdf((-2)*((log(p1)+log(p1))),4);pp=1+pp;
%         end
%         figure; plot(.01:.01:1,pvals);
%     end
% end
%axis(sctax)

figure(1600);
scatter((.5*simpleeffects(:,1)+.5*simpleeffects(:,3))-simpleeffects(:,2),simpleeffects(:,3)-simpleeffects(:,1),50,colDD,'fill'); title ('Cell responses', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Deviance Detection (Z)','FontSize',18,'FontWeight','bold'); xlabel('Stim specific adaptation (Z)','FontSize',18,'FontWeight','bold'); set(gca,'FontSize',16,'FontWeight','bold');set(gcf,'Color','w');
SSAs=(.5*simpleeffects(:,1)+.5*simpleeffects(:,3))-simpleeffects(:,2);DDs=(simpleeffects(:,3)-simpleeffects(:,1));
axis([-15 30 -15 30]);

% 
figure(1601); scatter((.5*simpleeffectsTEST(:,1)+.5*simpleeffectsTEST(:,3))-simpleeffectsTEST(:,2),simpleeffectsTEST(:,3)-simpleeffectsTEST(:,1),50,cols,'fill'); title ('Cell responses', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Deviance Detection (Z)','FontSize',18,'FontWeight','bold'); xlabel('Stim specific adaptation (Z)','FontSize',18,'FontWeight','bold'); set(gca,'FontSize',16,'FontWeight','bold');set(gcf,'Color','w');
figure(1602); scatter((.5*simpleeffects(:,1)+.5*simpleeffects(:,3))-simpleeffects(:,2),simpleeffects(:,3)-simpleeffects(:,1),50,cols,'fill'); title ('Cell responses', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Deviance Detection (Z)','FontSize',18,'FontWeight','bold'); xlabel('Stim specific adaptation (Z)','FontSize',18,'FontWeight','bold'); set(gca,'FontSize',16,'FontWeight','bold');set(gcf,'Color','w');


   sortby=[1 redundantnumber+1 10];
   BOOSTscalefactor=50; %between 0 and 100. higher is brighter
   figure(1700); 
   
   for rasters=1
       
       effects1=tmpeff;
       for c=1:size(tmpeff,1)
           for cnd=1:10
               m=mean(effects1(c,BEGINo-9:BEGINo-1,cnd));
               effects1(c,:,cnd)=effects1(c,:,cnd)-m;
           end
       end
       ak1=mean(effects1(:,:,sortby),3); [ix,b]=sort(mean(ak1(:,BEGINo:ENDo),2),'descend');
       effects2=effects1(b,:,:);
       
       for z1=1:size(effects2,3)
           subplot(1,size(effects2,3),z1); imagesc(timeaxis,1:size(effects2,1),effects2(:,:,z1)); colormap('Gray'); set(gca,'FontSize',14,'FontWeight','bold');
           xlim([-.2 1]);
           title(titles{z1});
           
           if z1==1; xlabel('sec','FontSize',14,'FontWeight','bold');
               ylabel('neurons','FontSize',14,'FontWeight','bold');
           end
       end
       for settingcolorscale=1;
           scalefactor=(BOOSTscalefactor/100); if scalefactor==1; scalefactor=.99; end
           tmpmax=0; tmpmin=0;
           for z1=1:size(effects2,3)
               subplot(1,size(effects2,3),z1);
               cb=caxis; if cb(1)<tmpmin; tmpmin=cb(1); end
               if cb(2)>tmpmax; tmpmax=cb(2); end
           end
           for z1=1:size(effects2,3)
               subplot(1,size(effects2,3),z1);
               cb=[tmpmin tmpmax];
               caxis(cb*(1-scalefactor))
           end
       end
   end


end
% 
