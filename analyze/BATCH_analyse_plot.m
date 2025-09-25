% add on part after running Batch_analysis_modified.m
% Focusing on the inhibitory part.
close all
include_cutoff=1; %number of standard deviations above baseline for inclusion Using this, then there won't be any inhibited cells. 
rdnt_vec=4:9;

% tagaddon='Unpublished';

% tagaddon='PYR_maus1';
% tagaddon='PYR_maus2';


tagaddon='CombinedBH';

% file_tag='./Figs/dlarger1';
% file_tag='./Figs/dlarger1ANDctrllarger1';
%file_tag='./Figs/dlarger1ANDctrlsmaller1';
% file_tag='./Figs/rdntlargerctrl';
% file_tag='./Figs/everything';
% file_tag='./Figs/ctrllarger1';
% file_tag='./Figs/' tagaddon 'ctrlNosigExcitDlarger1';
% file_tag=['./Figs/' tagaddon 'ctrlsigDvsCNosig'];

file_tag=['./Figs/' tagaddon 'everything'];
% file_tag=['./Figs/' tagaddon 'DvsCsigExcit'];
% file_tag=['./Figs/' tagaddon 'ReverseTest_DvsCsigExcit'];

% file_tag=['./Figs/' tagaddon 'ctrlsigDvsCsig'];  % pPE(x) to xxxy
% file_tag=['./Figs/' tagaddon 'ctrlsigDvsCsig_Reverse'];

% file_tag=['./Figs/' tagaddon 'CnoSig_DvsCsig'];   % nPE(y)
% file_tag=['./Figs/' tagaddon 'CnoSig_DvsCsig_Reverse'];

% file_tag=['./Figs/' tagaddon 'ReverseTest_CnoSig_DvsCsig'];


% file_tag=['./Figs/' tagaddon 'ctrlsigExcit'];
% file_tag=['./Figs/' tagaddon 'ctrlNosigExcit'];
% file_tag=['./Figs/' tagaddon 'ctrlsigInhibit'];

% file_tag=['./Figs/' tagaddon 'DvsCsiginhibit'];

% file_tag=['./Figs/' tagaddon 'CnoSig_DvsCsigInhib'];   % nPE(x)


load(['./Figs/' tagaddon 'TuningCurve.mat']);   % stim 1: 45 degree no. 3, stim 2: 135 degree no. 7 

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

for stringentTEST=0;
    if stringentTEST==0;
        %%%%%%%%%figures out relative responses to different stimuli
        simpleeffects=[]; simpleeffectsPRE=[];simpleeffectsPOST=[];c1=1; c2=1; cols=[]; newindices=[];notindices=[];
        for stim=1:2
            for c=1:size(effectsallSTD,1)
                if redundantnumberA==0; 
                	r=mean(mean(effectsallSTD(c,BEGINo:ENDo,rdnt_vec,stim)));
                else
                    r=mean(effectsallSTD(c,BEGINo:ENDo,redundantnumberA+1,stim));
                end
                cn=mean(effectsallSTD(c,BEGINo:ENDo,1,stim));
                d=mean(effectsallSTD(c,BEGINo:ENDo,10,stim));
                if redundantnumberA==0 
                	rP=mean(mean(effectsallBUILD(c,BEGINo:ENDo,rdnt_vec,stim)));
                else
                    rP=mean(effectsallBUILD(c,BEGINo:ENDo,redundantnumberA+1,stim));
                end
                cnP=mean(effectsallBUILD(c,BEGINo:ENDo,1,stim));
                dP=mean(effectsallBUILD(c,BEGINo:ENDo,10,stim));
                if redundantnumberA==0 
                	rT=mean(mean(effectsallTEST(c,BEGINo:ENDo,rdnt_vec,stim)));
                else
                    rT=mean(effectsallTEST(c,BEGINo:ENDo,redundantnumberA+1,stim));
                end
                cnT=mean(effectsallTEST(c,BEGINo:ENDo,1,stim));
                dT=mean(effectsallTEST(c,BEGINo:ENDo,10,stim));
                % Get the responsive cells
                % if and(or(or(r<include_cutoff,cn<include_cutoff),d<include_cutoff),and(and(r<exclude_outlier,cn<exclude_outlier),d<exclude_outlier))
                % if and(or(or(r>include_cutoff,cn>include_cutoff),d>include_cutoff),and(and(r<exclude_outlier,cn<exclude_outlier),d<exclude_outlier))
                % if and(and(d>include_cutoff,cn>include_cutoff),and(and(r<exclude_outlier,cn<exclude_outlier),d<exclude_outlier))

                %if and(and(d>include_cutoff,cn<include_cutoff),and(and(r<exclude_outlier,cn<exclude_outlier),d<exclude_outlier))
                
                % if (stim==1 && Tuning_res_sig_mat(c,3)>0.5) ||  (stim==2 && Tuning_res_sig_mat(c,7)>0.5)
                % if (stim==1 && Tuning_res_sig_mat(c,3)<0.5) ||  (stim==2 && Tuning_res_sig_mat(c,7)<0.5)
                % if (stim==1 && Tuning_res_sig_mat(c,3)>0.5 && d>include_cutoff) ||  (stim==2 && Tuning_res_sig_mat(c,7)>0.5 && d>include_cutoff)
                % if (stim==1 && Tuning_res_sig_mat(c,3)<0.5)&& d>include_cutoff ||  (stim==2 && Tuning_res_sig_mat(c,7)<0.5 && d>include_cutoff)
                % if and(1>0,and(and(r<exclude_outlier,cn<exclude_outlier),d<exclude_outlier))
                % if and(d>include_cutoff,and(and(r<exclude_outlier,cn<exclude_outlier),d<exclude_outlier))
                
                % New way using dev vs cntrl sig 

                if 1==1
                % if Results_table(c,10,stim)>0.5 % d vs C sig
                % if Results_table(c,10,3-stim)>0.5   % Reverse test
                % if Results_table(c,8,stim)>0.5 && Results_table(c,10,stim)>0.5    % pPE
                % if Results_table(c,8,3-stim)>0.5 && Results_table(c,10,3-stim)>0.5    % Reverse pPE
                % if Results_table(c,8,stim)<0.5 && Results_table(c,10,stim)>0.5    % nPE
                % if Results_table(c,8,3-stim)<0.5 && Results_table(c,10,3-stim)>0.5    % Reverse nPE


                % if (stim==1 && Results_table(c,10,1)>0.5 ) ||  (stim==2 &&  Results_table(c,10,2)>0.5)
                % if (stim==1 && Results_table(c,8,1)>0.5 && Results_table(c,10,1)>0.5 ) ||  (stim==2 && Results_table(c,8,2)>0.5 && Results_table(c,10,2)>0.5)
                % if (stim==1 && Results_table(c,8,1)<0.5 && Results_table(c,10,1)>0.5 ) ||  (stim==2 && Results_table(c,8,2)<0.5 && Results_table(c,10,2)>0.5)

                % if (stim==1 && Results_table(c,8,1)>0.5 ) ||  (stim==2 &&  Results_table(c,8,2)>0.5)
                % if (stim==1 && Results_table(c,8,1)<0.5 ) ||  (stim==2 &&  Results_table(c,8,2)<0.5)
                % if (stim==1 && Results_table(c,8,1)<-0.5 ) ||  (stim==2 &&  Results_table(c,8,2)<-0.5)
    
                % Why I just use Results_table(c,10,stim)  ???
                % if (stim==1 && Results_table(c,10,1)<-0.5 ) ||  (stim==2 &&  Results_table(c,10,2)<-0.5)

                % if (stim==1 && Results_table(c,8,1)<0.5 && Results_table(c,10,1)<-0.5 ) ||  (stim==2 && Results_table(c,8,2)<0.5 && Results_table(c,10,2)<-0.5)

                    
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
                    newindices(c1,3)=  Results_table(c,5,stim);
                    newindices(c1,4)=  Results_table(c,8,stim);

                    c1=1+c1;
                else
                    notindices(c2,1)=c;
                    notindices(c2,2)=stim;
                    c2=1+c2;
                    
                end
                
            end
        end
        save(strcat(file_tag, '_NewIndices.mat'), 'newindices');
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
        
    else  % Shut down here
        simpleeffects=[]; simpleeffectsPRE=[];simpleeffectsPOST=[];c1=1; c2=1; cols=[]; newindices=[];notindices=[];
        for stim=1:2
            for c=1:size(effectsallBUILD,1)
                if redundantnumberA==0; 
                	r=mean(mean(effectsallBUILD(c,BEGINo:ENDo,rdnt_vec,stim)));
                else
                    r=mean(effectsallBUILD(c,BEGINo:ENDo,redundantnumberA+1,stim));
                end
                cn=mean(effectsallBUILD(c,BEGINo:ENDo,1,stim));
                d=mean(effectsallBUILD(c,BEGINo:ENDo,10,stim));
                if redundantnumberA==0; 
                	rT=mean(mean(effectsallTEST(c,BEGINo:ENDo,rdnt_vec,stim)));
                else
                    rT=mean(effectsallTEST(c,BEGINo:ENDo,redundantnumberA+1,stim));
                end
                cnT=mean(effectsallTEST(c,BEGINo:ENDo,1,stim));
                dT=mean(effectsallTEST(c,BEGINo:ENDo,10,stim));
                if and(or(or(r<include_cutoff,cn<include_cutoff),d<include_cutoff),and(and(r<exclude_outlier,cn<exclude_outlier),d<exclude_outlier))
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
DDind_hist_data= (simpleeffects(:,3)-simpleeffects(:,1))./(simpleeffects(:,3)+simpleeffects(:,1)) ;  % This is normalized in Calcium data. 

figure(2002)
histogram(DDind_hist_data)
title('DDind hist')

figure(2004)
plot(simpleeffects(:,1),simpleeffects(:,3),'.',markersiz=20)
title('Deviant over stimulus')
figure(2005)
plot(simpleeffects(:,1),simpleeffects(:,3)-simpleeffects(:,1),'.',markersiz=20)
title('DD over stimulus')


figure(2006)
plot(simpleeffects(:,1),DDind_hist_data,'.',markersiz=20)
title('DDind over cntrl')






% Count the responsive cells
countcellstotal=0;
% Count the responsive cells
numcells=size(newindices(:,1),1);
disp(strcat('all cells=',num2str(size(effectsall,1)),'; responsive=',num2str(numcells),'; percent=',num2str(numcells/size(effectsall,1))));
uniq_numcells=size(unique(newindices(:,1)),1);
disp(strcat('all cells=',num2str(size(effectsall,1)),'; Unique responsive=',num2str(uniq_numcells),'; percent=',num2str(uniq_numcells/size(effectsall,1))));



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
                a1(c,t)=mean(effectsallBUILD(newindices(pind(c),1),t,redundantnumber+1,newindices(pind(c),2)),4);

            else
                a1tmp=mean(mean(effectsallBUILD(newindices(pind(c),1),t,rdnt_vec,newindices(pind(c),2)),3),4);
                a2tmp=mean(effectsallBUILD(newindices(pind(c),1),t,1,newindices(pind(c),2)),4);
                a3tmp=mean(effectsallBUILD(newindices(pind(c),1),t,10,newindices(pind(c),2)),4);
                mn=(a1tmp+a2tmp+a3tmp)/3;
                a1(c,t)=mean(mean(effectsallBUILD(newindices(pind(c),1),t,rdnt_vec,newindices(pind(c),2)),3),4);
            end
            a2(c,t)=mean(effectsallBUILD(newindices(pind(c),1),t,1,newindices(pind(c),2)),4);
            a3(c,t)=mean(effectsallBUILD(newindices(pind(c),1),t,10,newindices(pind(c),2)),4);
        end
    end
    x=timeaxis;
    y1=mean(mean(mean(tmpeff(:,:,redundantnumber+1),3),4),1); 
    if redundantnumber==0;y1=mean(mean(mean(tmpeff(:,:,rdnt_vec),3),4),1); end  
    y1a=y1-mean(y1(1,base_vec));
    err1=std(a1,0,1)./sqrt(numcells);
    
    y2=mean(mean(mean(tmpeff(:,:,1),3),4),1);
    y2a=y2-mean(y2(1,base_vec));
    err2=std(a2,0,1)./sqrt(numcells);
    
    y3=mean(mean(mean(tmpeff(:,:,10),3),4),1); 
    y3a=y3-mean(y3(1,base_vec));
    err3=std(a3,0,1)./sqrt(numcells);
    
    
    if numcells>1;
        figure(200); 
        shadedErrorBar(x,y1a,err1,'lineProps','b'); 
        hold on;shadedErrorBar(x,y2a,err2,'lineProps','k');
        hold on;shadedErrorBar(x,y3a,err3,'lineProps','r');
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
            r1=mean(mean(tmpeff(:,BEGINo:ENDo,redundantnumber+1),3),2)-mean(mean(tmpeff(:,base_vec,redundantnumber+1),3),2);
            if redundantnumber==0;
                 r1=mean(mean(tmpeff(:,BEGINo:ENDo,rdnt_vec),3),2)-mean(mean(tmpeff(:,base_vec,rdnt_vec),3),2);
            end
            c1=mean(mean(tmpeff(:,BEGINo:ENDo,1),3),2)-mean(mean(tmpeff(:,base_vec,1),3),2);
            d1=mean(mean(tmpeff(:,BEGINo:ENDo,10),3),2)-mean(mean(tmpeff(:,base_vec,10),3),2);
            
            [h,p1,ci,stats1]=ttest(r1,c1);
            [h,p2,ci,stats2]=ttest(d1,c1);
        end
        title(strcat('pre trials_T_S_S_A','=',num2str(round(stats1.tstat*100)/100),', p=',num2str(round(p1*1000)/1000),'; T_D_D','=',num2str(round(stats2.tstat*100)/100),', p=',num2str(round(p2*1000)/1000)));
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
            y1=y1-mean(y1(1,baseBEGIN:baseEND));
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
                a1(c,t)=mean(effectsallTEST(newindices(pind(c),1),t,redundantnumber+1,newindices(pind(c),2)),4);
            else
                 a1tmp=mean(mean(effectsallTEST(newindices(pind(c),1),t,rdnt_vec,newindices(pind(c),2)),3),4);
                a2tmp=mean(effectsallTEST(newindices(pind(c),1),t,1,newindices(pind(c),2)),4);
                a3tmp=mean(effectsallTEST(newindices(pind(c),1),t,10,newindices(pind(c),2)),4);
                mn=(a1tmp+a2tmp+a3tmp)/3;
                a1(c,t)=mean(mean(effectsallTEST(newindices(pind(c),1),t,rdnt_vec,newindices(pind(c),2)),3),4);
            end
            a2(c,t)=mean(effectsallTEST(newindices(pind(c),1),t,1,newindices(pind(c),2)),4);
            a3(c,t)=mean(effectsallTEST(newindices(pind(c),1),t,10,newindices(pind(c),2)),4);
        end
    end
    
    x=timeaxis;
    y1=mean(mean(mean(tmpeff(:,:,redundantnumber+1),3),4),1); 
    
    if redundantnumber==0; y1=mean(mean(mean(tmpeff(:,:,rdnt_vec),3),4),1); end; 
    y1a=y1-mean(y1(1,base_vec));
    err1=std(a1,0,1)./sqrt(numcells);
    
    y2=mean(mean(mean(tmpeff(:,:,1),3),4),1);
    y2a=y2-mean(y2(1,base_vec));
    err2=std(a2,0,1)./sqrt(numcells);
    
    y3=mean(mean(mean(tmpeff(:,:,10),3),4),1); 
    y3a=y3-mean(y3(1,base_vec));
    err3=std(a3,0,1)./sqrt(numcells);
    
    
    if numcells>1;
        figure(700); 
        shadedErrorBar(x,y1a,err1,'lineProps','b'); 
        hold on;shadedErrorBar(x,y2a,err2,'lineProps','k');hold on;shadedErrorBar(x,y3a,err3,'lineProps','r');
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
                 r1=mean(mean(tmpeff(:,BEGINo:ENDo,rdnt_vec),3),2)-mean(mean(tmpeff(:,base_vec,rdnt_vec),3),2);
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
            y1=y1-mean(y1(1,baseBEGIN:baseEND));
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

for plotbOLTH=1
    effectsallBOTH= effectsallSTD; %  Here should be effectsallSTD. Original effectsall, another option  (effectsallBUILD+effectsallTEST)/2;   ???
    DEBUGeffect2=effectsallBOTH;
    tmpeff=[];
    for c=1:numcells;
        tmpeff(c,:,:)=effectsallBOTH(newindices(pind(c),1),:,:,(newindices(pind(c),2)));
        for t=1:size(effectsallBOTH,2)
            if redundantnumber>0;
                mn=mean(mean(effectsallBOTH(newindices(pind(c),1),t,[redundantnumber+1 1 10],newindices(pind(c),2)),3),4);
                a1(c,t)=mean(effectsallBOTH(newindices(pind(c),1),t,redundantnumber+1,newindices(pind(c),2)),4);
            else
                
                a1tmp=mean(mean(effectsallBOTH(newindices(pind(c),1),t,rdnt_vec,newindices(pind(c),2)),3),4);
                a2tmp=mean(effectsallBOTH(newindices(pind(c),1),t,1,newindices(pind(c),2)),4);
                a3tmp=mean(effectsallBOTH(newindices(pind(c),1),t,10,newindices(pind(c),2)),4);
                mn=(a1tmp+a2tmp+a3tmp)/3;
                a1(c,t)=mean(mean(effectsallBOTH(newindices(pind(c),1),t,rdnt_vec,newindices(pind(c),2)),3),4);
            end
            a2(c,t)=mean(effectsallBOTH(newindices(pind(c),1),t,1,newindices(pind(c),2)),4);
            a3(c,t)=mean(effectsallBOTH(newindices(pind(c),1),t,10,newindices(pind(c),2)),4);
            afirst(c,t)=mean(effectsallBOTH(newindices(pind(c),1),t,2,newindices(pind(c),2)),4);

        end
    end
   

    x=timeaxis;
    y1=mean(mean(mean(tmpeff(:,:,redundantnumber+1),3),4),1);
    if redundantnumber==0; y1=mean(mean(mean(tmpeff(:,:,rdnt_vec),3),4),1); end;
    y1a=y1-mean(y1(1,base_vec));
    % y1a=y1-min(y1(1,base_vec));
    err1=std(a1,0,1)./sqrt(numcells);
    
    y2=mean(mean(mean(tmpeff(:,:,1),3),4),1);
    y2a=y2-mean(y2(1,base_vec));
    % y2a=y2-min(y2(1,base_vec));
    err2=std(a2,0,1)./sqrt(numcells);
    
    y3=mean(mean(mean(tmpeff(:,:,10),3),4),1); 
    y3a=y3-mean(y3(1,base_vec));
    % y3a=y3-min(y3(1,base_vec));
    err3=std(a3,0,1)./sqrt(numcells);

    yfirst=mean(mean(mean(tmpeff(:,:,2),3),4),1); 
    yfirsta=yfirst-mean(yfirst(1,base_vec));
    % yfirsta=yfirst-min(yfirst(1,base_vec));
    errfirst=std(afirst,0,1)./sqrt(numcells);

    % mm=min(horzcat(y1a(1,1:BEGINo-1),y2a(1,1:BEGINo-1),y3a(1,1:BEGINo-1))); % Normalized by the mm of the baseline.
    mm=0;
    y1a=y1a-mm;y2a=y2a-mm;y3a=y3a-mm;
    
    if numcells>0
        fig=figure(1000); shadedErrorBar(x,y1a,err1,'lineProps','b'); 
        hold on;
        shadedErrorBar(x,y2a,err2,'lineProps','k');
        shadedErrorBar(x,y3a,err3,'lineProps','r');
        shadedErrorBar(x,yfirsta,errfirst,'lineProps',{'-','color','#FF8C00'});


        xlim([-.2 1]); title ('Average of cell responses', 'FontSize', 16, 'FontWeight', 'bold');ylabel('z-scores','FontSize',16,'FontWeight','bold'); xlabel('time (sec)','FontSize',16,'FontWeight','bold'); set(gca,'FontSize',16,'FontWeight','bold');
        set(gcf,'Color','w');
        figname= 'Pop Ave';
        file_name= strcat(file_tag,'_',figname,'.jpg');
        saveas(fig,file_name)        
        file_name= strcat(file_tag,'_',figname,'.svg');
        saveas(fig,file_name)        


        yy=ylim;
        % figure(1001); 
        errorbar_groups([mean(mean(y2a(:,BEGINo:ENDo)))  mean(mean(yfirsta(:,BEGINo:ENDo)))    mean(mean(y1a(:,BEGINo:ENDo))) mean(mean(y3a(:,BEGINo:ENDo)))]', ...
            [mean(mean(err2(:,BEGINo:ENDo))) mean(mean(errfirst(:,BEGINo:ENDo)))   mean(mean(err1(:,BEGINo:ENDo))) mean(mean(err3(:,BEGINo:ENDo)))]', ...
            'bar_colors',[.5 .5 .5; .8 .3 0; 0 0 1; 1 0 0],'bar_width', 1.3,'FigID',1001);
        % ylim([0, 7])   % pPE

        % ylim([-1, 1.5]) % nPE
        ylabel('z-score','FontSize',16,'FontWeight','bold'); set(gca,'FontSize',16,'FontWeight','bold');
        set(gcf,'Color','w');
        for statstat=1;
            r1=mean(mean(tmpeff(:,BEGINo:ENDo,redundantnumber+1),3),2)-mean(mean(tmpeff(:,base_vec,redundantnumber+1),3),2);
            if redundantnumber==0;
                 r1=mean(mean(tmpeff(:,BEGINo:ENDo,rdnt_vec),3),2)-mean(mean(tmpeff(:,base_vec,rdnt_vec),3),2);
            end
            c1=mean(mean(tmpeff(:,BEGINo:ENDo,1),3),2)-mean(mean(tmpeff(:,base_vec,1),3),2);
            d1=mean(mean(tmpeff(:,BEGINo:ENDo,10),3),2)-mean(mean(tmpeff(:,base_vec,10),3),2);
            first1=mean(mean(tmpeff(:,BEGINo:ENDo,2),3),2)-mean(mean(tmpeff(:,base_vec,2),3),2);
            [h,p1,ci,stats1]=ttest(r1,c1);
            [h,p2,ci,stats2]=ttest(d1,c1);
            [h,p3,ci,stats3]=ttest(first1,c1);

            [h,p4,ci,stats4]=ttest(r1,first1);
            [h,p5,ci,stats5]=ttest(d1,first1);

        end
        set(gcf,'units','points','position',[200,200,400,300])

        title({strcat('T_S_S_A','=',num2str(round(stats1.tstat*100)/100),', p=',num2str(round(p1*1000)/1000), ...
            '; T_D_D','=',num2str(round(stats2.tstat*100)/100),', p=',num2str(round(p2*1000)/1000), ...
            '; T_F','=',num2str(round(stats3.tstat*100)/100),', p=',num2str(round(p3*1000)/1000))...
            ,strcat('p_RvsF=',num2str(round(p4*1000)/1000), 'p_DvsF=',num2str(round(p5*1000)/1000)) });
        figname= 'Ttest';
        file_name= strcat(file_tag,'_',figname,'.jpg');
        saveas(gcf,file_name)        
        file_name= strcat(file_tag,'_',figname,'.svg');
        saveas(gcf,file_name)        


    else
        figure(1002); plot(x,y1,'b','LineWidth',3); hold on;plot(x,y2,'k','LineWidth',3); hold on;plot(x,y3,'r','LineWidth',3);yy=ylim;
        xlim([-.2 1]); title ('Average of cell responses', 'FontSize', 16, 'FontWeight', 'bold');
    end
    
    for manyplotthing=1;
        figure(1100);
        for cnd=1:10;
            y1=mean(tmpeff(:,:,cnd),1);
            y1=y1-mean(y1(1,base_vec));
            % y1=y1-min(y1(1,base_vec));
            if cnd==1; err1=std(a2,0,1)./sqrt(numcells);
                err=std(a2,0,1)./sqrt(numcells);
                
            elseif cnd==10
                err=std(a3,0,1)./sqrt(numcells);
            else
                err=std(a1,0,1)./sqrt(numcells);
            end
            subplot(1,10,cnd); 
            hold on
            shadedErrorBar(x,y1,err,'lineProps','k'); 
            plot([-.2 .8],[1,1]*mean(y1(1,BEGINo:ENDo)),  '--k')
            set(gcf,'Color','w'); xlim([-.2 1]);
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
        set(gcf,'units','points','position',[200,400,1000,250])
        figname= 'Ind Trace';
        file_name= strcat(file_tag,'_',figname,'.jpg');
        saveas(gcf,file_name)        
        file_name= strcat(file_tag,'_',figname,'.svg');
        saveas(gcf,file_name)        


    end
    
    for manyplotthing=2;
        
        if redundantnumber==0;
            figure(1200);
            ccnds={1;rdnt_vec;10};
            cccf=[1 5 10];
            
            colars={'k' 'b' 'r'};
            for cnd=1:3;
                y1=mean(mean(tmpeff(:,:,ccnds{cnd}),3),1);
                y1=y1-mean(y1(1,base_vec));
                % y1=y1-min(y1(1,base_vec));
                if cnd==1; err1=std(a2,0,1)./sqrt(numcells);
                    err=std(a2,0,1)./sqrt(numcells);
                    
                elseif cnd==3;
                    err=std(a3,0,1)./sqrt(numcells);
                else
                    err=std(a1,0,1)./sqrt(numcells);
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
                y1=y1-mean(y1(1,base_vec));
                if cnd==1; err1=std(a1,0,1)./sqrt(numcells);
                    err=std(a1,0,1)./sqrt(numcells);
                    
                elseif cnd==3
                    err=std(a3,0,1)./sqrt(numcells);
                else  % cnd=2
                    err=std(a2,0,1)./sqrt(numcells);
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
    figure(1400); 
    shadedErrorBar(x,rdnt1,rdnt1err,'lineProps','b'); hold on;
    shadedErrorBar(x,cont1,cont1err,'lineProps','k');hold on;
    shadedErrorBar(x,dev1,dev1err,'lineProps','r');
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


    % sortby=[1 redundantnumber+1 10];
    sortby= [1];
    BOOSTscalefactor=50; %between 0 and 100. higher is brighter
    % Define RGB values
    start_color = [0, 0.5, 0];      % Green
    middle_color = [1, 1, 1];       % White
    end_color = [0.9686, 0.6471, 0.2549]; % #F7A541 (Orange)
    
    % Number of colors
    n = 256;  % Adjust resolution if needed
    
    % Interpolate colors across 256 levels
    x = [1, n/2, n]; % Position of colors in colormap
    r = interp1(x, [start_color(1), middle_color(1), end_color(1)], 1:n);
    g = interp1(x, [start_color(2), middle_color(2), end_color(2)], 1:n);
    b = interp1(x, [start_color(3), middle_color(3), end_color(3)], 1:n);
    
    % Combine into colormap matrix
    custom_colormap = [r(:), g(:), b(:)];
        
    % Apply the colormap
    figure(1700); 
    clf
   for rasters=1  
       effects1=tmpeff;
       for c=1:size(tmpeff,1)
           for cnd=1:10
               m=mean(effects1(c,base_vec,cnd));
               effects1(c,:,cnd)=effects1(c,:,cnd)-m;
           end
       end
       ak1=mean(effects1(:,:,sortby),3); 
       [ix,rev_ind]=sort(mean(ak1(:,BEGINo:ENDo),2),'descend');
       effects2=effects1(rev_ind,:,:);
       linthresh = 1;
       log_norm_data = sign(effects2(:,:,:)) .* log10(1+abs(effects2(:,:,:)) / linthresh);
       max_df= max(log_norm_data(:))*0.8;

       for z1=1:size(effects2,3)
           subplot(1,size(effects2,3),z1); 

           imagesc(timeaxis,1:size(effects2,1),log_norm_data(:,:,z1)); 
           % colormap('Gray'); 
           colormap(custom_colormap);
           % c =colorbar;
           % c.Ticks = [-log10(1+10/linthresh), -log10(1+1/linthresh), 0, log10(1+1/linthresh), log10(1+10/linthresh)];  % Set tick positions
           % c.TickLabels = {'-10', '-1', '0', '1', '10'}; % Set tick labels 
           set(gca,'FontSize',14,'FontWeight','bold');
           xlim([-.2 1]);
           clim([-max_df, max_df]); % Set limits

           title(titles{z1});
           
           if z1==1; xlabel('sec','FontSize',14,'FontWeight','bold');
               ylabel('neurons','FontSize',14,'FontWeight','bold');
               yticks([0, 10, 20]); % Custom y-tick positions
           else
               yticklabels([]); % Hide y-axis labels but keep ticks
           end
       end
      set(gcf,'units','points','position',[200,400,1000,250])

       figname= 'Individuals';
        file_name= strcat(file_tag,'_',figname,'.jpg');
        saveas(gcf,file_name)        
        file_name= strcat(file_tag,'_',figname,'.svg');
        saveas(gcf,file_name)        

   end
       % for settingcolorscale=1
       %     scalefactor=(BOOSTscalefactor/100); if scalefactor==1; scalefactor=.99; end
       %     tmpmax=0; tmpmin=0;
       %     for z1=1:size(effects2,3)
       %         subplot(1,size(effects2,3),z1);
       %         cb=caxis; 
       %         if cb(1)<tmpmin; tmpmin=cb(1); end
       %         if cb(2)>tmpmax; tmpmax=cb(2); end
       %     end
       %     for z1=1:size(effects2,3)
       %         subplot(1,size(effects2,3),z1);
       %         cb=[tmpmin tmpmax];
       %         caxis(cb*(1-scalefactor))
       %     end
       % 
       % end

end
%  Check the activity and p_value
% shuffled_indices=newindices(rev_ind,:,:)