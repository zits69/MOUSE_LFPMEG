addpath '/home/ezp/Documents'/MATLAB/fieldtrip-20200914/%need Fieldtrip toolbox
addpath '/home/ezp/Documents/scripts_and_stats/'%need plotpatch function

%%
[Bsame X fullGLMdata_LFP_NORM fullGLMdata_REGS_NORM BMwinner] = lfp_regression_allcontacts_simple(1,0,9,1,0);%the normalization model tailored
[X1 BEV fullGLMdata_LFP_BEV fullGLMdata_REGS_BEV] = lfp_regression_allcontacts_simple(1,0,4,0,0);%the basic DVs (absolute value of integrated evidece (confidence))
[X2 BU fullGLMdata_LFP_BU fullGLMdata_REGS_BU] = lfp_regression_allcontacts_simple(1,0,5,0,0);%cue number (urgency)

%the LFP output is the same, just extracting for sanity check

%%%plotting separate GLMs (pairwise against same)
time=[-0.5:0.05:0.8];

figure;set(gcf,'color','w');
title({'Evidence over all cues and STN contacts','(pairwise GLM: cue ID vs other)'}, 'FontSize',15);

p1=plotpatch (Bsame, time,'k');
hold on
p3=plotpatch (BMwinner, time,'b');
p4=plotpatch (BEV, time,'g');
p5=plotpatch (BU, time,'m');

%stats on GLM - FDR corrected (rather than cluster corerction within
%lfp_regression
for i=11:27;[h(i) p(i)]=ttest(Bsame(:,i));end
p(1:10)=[];%the precue period
[hc, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p,0.05,'pdep','yes')
sig = find (adj_p<0.05);
sig=sig+10;
plot (time(sig), ones(1, length(sig))*0.017, 'k.','MarkerSize',15);

clear p sig

for i=11:27;[h(i) p(i)]=ttest(BMwinner(:,i));end
p(1:10)=[];%the precue period
[hc, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p,0.05,'pdep','yes')
sig = find (adj_p<0.05);
sig=sig+10;
plot (time(sig), ones(1, length(sig))*0.017, 'b.','MarkerSize',15);

clear p sig

for i=11:27;[h(i) p(i)]=ttest(BEV(:,i));end
p(1:10)=[];%the precue period
[hc, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p,0.05,'pdep','yes')
sig = find (adj_p<0.05);
sig=sig+10;
plot (time(sig), ones(1, length(sig))*0.017, 'g.','MarkerSize',15);

clear p sig

for i=11:27;[h(i) p(i)]=ttest(BU(:,i));end
p(1:10)=[];%the precue period
[hc, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p,0.05,'pdep','yes')
sig = find (adj_p<0.05);
sig=sig+10;
plot (time(sig), ones(1, length(sig))*0.017, 'm.','MarkerSize',15);

clear p sig

hold off

yline(0,'k');ylim([-0.025,0.025]);
lgd=legend ([p1(1) p3(1) p4(1) p5(1) ],'Box', 'off');
lgd.String={'cue identity (same or different)', 'normalization model/global conflict', 'absolute evidence','urgency/cue number/WM'};
xlabel ('Time from cue onset', 'FontSize',15);
ylabel ('BETA Power Regression coefficient', 'FontSize',15);
lgd.Location ='northeast';
lgd.FontSize=10;

%% %do full GLM with all 4 predictors combined
for i=1:13
    
    regs(:,1)=fullGLMdata_REGS_NORM{i}{1};
    regs(:,2)=fullGLMdata_REGS_NORM{i}{2};
    regs(:,3)=fullGLMdata_REGS_BEV{i}{2};
    regs(:,4)=fullGLMdata_REGS_BU{i}{2};
    lfp=fullGLMdata_LFP_NORM{i};
    
    for d = 1:length(lfp)
        
        y = lfp{d};
        
        %GLM with glmfit
        B = glmfit(regs,y);
        Bsame_fullglm(i,d) = B(2);
        BNORM_fullglm(i,d) = B(3);
        BEV_fullglm(i,d) = B(4);
        BU_fullglm(i,d) = B(5);
        
        %GLM with fitglm
        table_data=table; %create empty table
        table_data.same=regs(:,1);
        table_data.norm=regs(:,2);
        table_data.evidence=regs(:,3);
        table_data.urgency=regs(:,4);
        table_data.lfp=y;
        mdl=fitglm(table_data,'lfp ~ same*(norm+evidence+urgency)');% linear + interaction with same
        %            mdl=fitglm(table_data,'lfp ~ same+norm+evidence+urgency');%  only linear and equiv to glmfit
        Bsame_fullglm2(i,d) = mdl.Coefficients.Estimate(2);
        BNORM_fullglm2(i,d) = mdl.Coefficients.Estimate(3);
        BEV_fullglm2(i,d) = mdl.Coefficients.Estimate(4);
        BU_fullglm2(i,d) = mdl.Coefficients.Estimate(5);
        BSameNORM_fullglm2(i,d) = mdl.Coefficients.Estimate(6);
        BSameEV_fullglm2(i,d) = mdl.Coefficients.Estimate(7);
        BSameU_fullglm2(i,d) = mdl.Coefficients.Estimate(8);
        
    end
    
    
    clear regs lfp table_data
end

%%
stats=2;%1=FDR;2=cluster
CI=0;%with confidence interval(CI) or just SEM
%%
figure;set(gcf,'color','w');
if stats==1
    title('Evidence over all cues and STN contacts (full GLM) FDR', 'FontSize',15);
else
    title({'Evidence over all cues and STN contacts', '(full GLM) cluster'}, 'FontSize',15);
end

if CI==0 %this will plot with s.e.m
    
p1=plotpatch (Bsame_fullglm, time,'k');
p3=plotpatch (BNORM_fullglm, time,'g');
p4=plotpatch (BEV_fullglm, time,'b');
p5=plotpatch (BU_fullglm, time,'m');

elseif CI==1 %this will plot with confidence interval which is more appropriate for comparisons against zero, as here

    for p=1:4
        if p==1
        y = Bsame_fullglm;                             % Create Dependent Variable Experiments Data
        elseif p==2
        y = BNORM_fullglm;                             % Create Dependent Variable Experiments Data
        elseif p==3
        y = BEV_fullglm;                             % Create Dependent Variable Experiments Data
        elseif p==4
        y = BU_fullglm;                             % Create Dependent Variable Experiments Data
        end
        
        x = time;                                          % Create Independent Variable
        N = size(y,1);                                     % Number of Experiments In Data Set
        yMean = mean(y,1);                                 % Mean Of All Experiments At Each Value Of x
        ySEM = nanstd(y,'',1)/sqrt(N);                     % Compute Standard Error Of The Mean Of All Experiments At Each Value Of x
        CI95 = tinv([0.025 0.975], N-1);                   % Calculate 95% Probability Intervals Of t-Distribution
        yCI95 = bsxfun(@times, ySEM, CI95(:));             % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of x
        
        if p==1
            c='k';c2='k--';
        elseif p==2
            c='g'; c2='g--';
        elseif p==3
            c='b';c2='b--';
        elseif p==4
            c='m';c2='m--';
        end
        
        xx=[x fliplr(x)];
        yp  = [yCI95(2,:)+yMean fliplr(yMean-yCI95(2,:))];        
        patch(xx, yp,c,'EdgeColor','none','FaceAlpha',0.2);
        
        hold on
        
        plot(x, yMean,c,'LineWidth',2)                                      % Plot Mean Of All Experiments
%         plot(x, yCI95+yMean,c2,'LineWidth',2)                               % Plot 95% Confidence Intervals Of All Experiments
%         hold off
        xlabel ('Time from cue onset', 'FontSize',15);
        xlim([-0.5 0.8])
        % xline(0,'k');xline(0.8,'k');xline(1.6,'k');
        set(gca,'FontSize',15)
    end

end

if stats==1
%%FDR
%stats on GLM - FDR corrected (rather than cluster corerction within
%lfp_regression
        for i=11:27;[h(i) p(i)]=ttest(Bsame_fullglm(:,i));end
        p(1:10)=[];%the precue period
        [hc, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p,0.05,'pdep','yes')
        sig = find (adj_p<0.05);
        if isempty(sig)==0
        sig=sig+10;
        plot (time(sig), ones(1, length(sig))*0.015, 'k.','MarkerSize',15);
        else
        sig = find (adj_p<0.1);
        sig=sig+10;
        plot (time(sig), ones(1, length(sig))*0.015, 'k+','MarkerSize',5);
        end
        clear p sig

        for i=11:27;[h(i) p(i)]=ttest(BNORM_fullglm(:,i));end
        p(1:10)=[];%the precue period
        [hc, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p,0.05,'pdep','yes')
        sig = find (adj_p<0.05);
        if isempty(sig)==0
        sig=sig+10;
        plot (time(sig), ones(1, length(sig))*0.017, 'g.','MarkerSize',15);
        else
        sig = find (adj_p<0.1);
        sig=sig+10;
        plot (time(sig), ones(1, length(sig))*0.017, 'g+','MarkerSize',5);
        end
        clear p sig

        for i=11:27;[h(i) p(i)]=ttest(BEV_fullglm(:,i));end
        p(1:10)=[];%the precue period
        [hc, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p,0.05,'pdep','yes')
        sig = find (adj_p<0.05);
        sig=sig+10;
        plot (time(sig), ones(1, length(sig))*0.019, 'b.','MarkerSize',15);

        clear p sig

        for i=11:27;[h(i) p(i)]=ttest(BU_fullglm(:,i));end
        p(1:10)=[];%the precue period
        [hc, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p,0.05,'pdep','yes')
        sig = find (adj_p<0.05/4);
        sig=sig+10;
        plot (time(sig), ones(1, length(sig))*0.021, 'm.','MarkerSize',15);

        clear p sig

        hold off

elseif stats==2
%%CLUSTER

        cfg=[];
        cfg.statistic = 'ft_statfun_depsamplesT';
        cfg.method ='montecarlo';
        cfg.correctm = 'cluster';
        cfg.clusteralpha  = 0.05;%cluster forming thhreshold...it's complicated and FSL uses the same, that first you create some clusters, which need to be at a certain threshold, then you do stats on those, which is the values below
        cfg.alpha       = 0.05;
        cfg.tail        = 0; % two-sided test
        cfg.correcttail = 'prob';
        cfg.numrandomization = 1000;%the higher the better, 1000 is enough to ensur
        cfg.latency = [0 0.8];
        cfg.avgovertime = 'no';
        cfg.neighbours    =[];

        Nsub=size(Bsame_fullglm,1);

        subj = Nsub;
        subj2 = Nsub;
        design = zeros(2,subj+subj2);
        for i = 1:subj
            design(1,i) = i;
        end
        for i = 1:subj2
            design(1,subj+i) = i;
        end
        design(2,1:subj)        = 1;
        design(2,subj+1:subj+subj2) = 2;

        cfg.design   = design;
        cfg.uvar     = 1;
        cfg.ivar     = 2;
        
        avsame_group.avg=reshape(Bsame_fullglm,[Nsub, 1, size(time,2)]);
        dummy.avg=zeros(size(avsame_group.avg));
        avsame_group.dimord = 'subj_chan_time';dummy.dimord = 'subj_chan_time';
        avsame_group.label={'LFP'};dummy.label={'LFP'};
        avsame_group.time=time;dummy.time=time;

        [stat] = ft_timelockstatistics(cfg, avsame_group, dummy);%
        
%         if (stat.posclusters.prob < (0.05/4)) || (stat.negclusters(1).prob < (0.05/4))
        sig = find (stat.mask==1);%<0.05
%         sig = find (stat.prob<0.05/4);%correct for MC because 4 regs against zero
        sig=sig+10;
        plot (time(sig), ones(1, length(sig))*0.019, 'k.','MarkerSize',15);
%         else
%         end
        clear sig stat
        
        avsame_group.avg=reshape(BNORM_fullglm,[Nsub, 1, size(time,2)]);
        
        [stat] = ft_timelockstatistics(cfg, avsame_group, dummy);%
        
        sig = find (stat.mask==1);%<0.05
%         sig = find (stat.prob<0.05/4);%correct for MC because 4 regs against zero
        sig=sig+10;
        plot (time(sig), ones(1, length(sig))*0.021, 'g.','MarkerSize',15);

        clear sig stat
        
        avsame_group.avg=reshape(BEV_fullglm,[Nsub, 1, size(time,2)]);
        
        [stat] = ft_timelockstatistics(cfg, avsame_group, dummy);%
        
        sig = find (stat.mask==1);%<0.05
%         sig = find (stat.prob<0.05/4);%correct for MC because 4 regs against zero
        sig=sig+10;
        plot (time(sig), ones(1, length(sig))*0.023, 'b.','MarkerSize',15);
        
        clear sig stat
        
        avsame_group.avg=reshape(BU_fullglm,[Nsub, 1, size(time,2)]);
        
        [stat] = ft_timelockstatistics(cfg, avsame_group, dummy);%
        
        sig = find (stat.mask==1);%<0.05
%         sig = find (stat.prob<0.05/4);%correct for MC because 4 regs against zero
        sig=sig+10;
        plot (time(sig), ones(1, length(sig))*0.025, 'm.','MarkerSize',15);
        
end
% figure;set(gcf,'color','w');
% title('Evidence over all cues and STN contacts (full GLM with same interaction)', 'FontSize',15);
% 
% p1=plotpatch (Bsame_fullglm2, time,'k');
% p3=plotpatch (BNORM_fullglm2, time,'b');
% p4=plotpatch (BEV_fullglm2, time,'g');
% p5=plotpatch (BU_fullglm2, time,'m');
% 
% %stats on GLM - FDR corrected (rather than cluster corerction within
% %lfp_regression
% for i=1:27;[h(i) p(i)]=ttest(Bsame_fullglm2(:,i));end
% [hc, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p,0.05,'pdep','yes')
% 
% sig = find (adj_p<0.05);
% plot (time(sig), ones(1, length(sig))*0.015, 'k.','MarkerSize',15);
% 
% clear p sig
% 
% for i=1:27;[h(i) p(i)]=ttest(BNORM_fullglm2(:,i));end
% [hc, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p,0.05,'pdep','yes')
% 
% sig = find (adj_p<0.05);
% plot (time(sig), ones(1, length(sig))*0.015, 'b.','MarkerSize',15);
% 
% clear p sig
% 
% for i=1:27;[h(i) p(i)]=ttest(BEV_fullglm2(:,i));end
% [hc, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p,0.05,'pdep','yes')
% 
% sig = find (adj_p<0.05);
% plot (time(sig), ones(1, length(sig))*0.015, 'g.','MarkerSize',15);
% 
% clear p sig
% 
% for i=1:27;[h(i) p(i)]=ttest(BU_fullglm2(:,i));end
% [hc, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p,0.05,'pdep','yes')
% 
% sig = find (adj_p<0.05);
% plot (time(sig), ones(1, length(sig))*0.017, 'm.','MarkerSize',15);
% 
% clear p sig
% 
% hold off

yline(0,'k');ylim([-0.04,0.04]);xline(0,'k');
lgd=legend ([p1(1) p3(1) p4(1) p5(1) ],'Box', 'off');
lgd.String={'cue identity (different - same)', 'normalization model / global conflict', 'absolute evidence','urgency / cue number / WM'};
xlabel ('Time from cue onset', 'FontSize',15);
if CI==0
ylabel ('BETA Power Regression coefficient', 'FontSize',15);
% ylabel ('THETA Power Regression coefficient', 'FontSize',15);
else
ylabel ('BETA Power Regression coefficient (+95%CI)', 'FontSize',15);
end
lgd.Location ='northeast';
lgd.FontSize=10;


%%
function [Bsame Badd fullGLMdata_LFP fullGLMdata_REGS BMwinner time BM1 BM2 BM3 BM4 ] = lfp_regression_allcontacts_simple (freq, last, additional, tailored,excl)

% function regression {(last, additional)}
%
% Function performs regression of LFP activity on the parameters on
% various characteristic of stimuli
% Input:
% % freq: beta (1) or theta (2)
% % tailored: should model fitting (if additional ==9) but per within subject winning
% model or not)
% Optional Inputs:
%  last - if =1, the analysis is performed on the last stimulus in the
%         sequence, if=0, it is perfroming on all stimuli excluding the first one
%         and the last one (default).
%  additional - an additional regressor to inlcude: 0 = none (default)
%               1 = does the stimulus match the average of previous stim
%               2 = evidence for the current stimulus
%               3 = integrated evidence for the ipsilateral side - not possible when averaging over contacts
%               4 = absolute value of integrated evidece (confidence)
%               5 = index of the stimulus in the sequence (urgency signal)
%               6 = match of previous stim but simple binary (very similar to additional #1)
%               7 = match to response
%               8 = match to accuracy
%               9 = models est from behaviour: forgetting, bonuses and various combos
%%M1: DVt = DVt-1 + Xt
%%M2: DVt = (1-LAMBDA)DVt-1 + Xt
%%M3: DVt = DVt-1 + WtXt
%%M4: DVt = (1-LAMBDA)DVt-1 + WtXt
%               10 = same as 9 but only for the 2nd stimulus rather than same reg for all stim


if nargin < 2
    last = 0;
end
if nargin < 3
    additional = 0;
end

if freq==1
files = dir('/home/ezp/Documents/STN/GitHub/BETA_LFP/*_cue.mat');
else
files = dir('/home/ezp/Documents/STN/GitHub/THETA_LFP/*_cue.mat');
end    

n = numel(files);

duration = 1:27;

for i = 1:n
    if freq==1
    load (['GitHub/BETA_LFP/' files(i).name]);
    else
    load (['GitHub/THETA_LFP/' files(i).name]);
    end
    
    side(i) = (files(i).name(12) == 'L');%obsolete as all contacts are being averaged
    lfp = squeeze (nanmean(data.trial,2));
    include = [];
    for t = 1:length(seq)
        index = seq{t}.ind; %index of cue within a trial
        if last && index == length (seq{t}.seq)
            include = [include, t];
        end
        if additional <10 && ~last && index > 1 && index < length (seq{t}.seq) && isempty(find(bad==t))  %ignore the first and the last trial
            include = [include, t];
        end
        if additional == 10 && ~last && index ==2  && isempty(find(bad==t))  %only take 2nd stim
            include = [include, t];
        end
        %gather stimuli and responses
        stim(t) = seq{t}.seq(index);
        temp = seq{t}.resp(1);if strcmp(temp,'L')==1; choice(t)=1;else choice(t)=2;end
        temp = seq{t}.resp;if strcmp(temp,'Left_correct')==1 || strcmp(temp,'Right_correct')==1 ; acc(t)=1;else acc(t)=0;end
        if index == 1
            evidence(t) = 2*stim(t)-3;
            same(t) = 0;
            evstim(t) = 0;
            match(t) = 0;
        else
            evidence(t) = evidence(t-1) + 2*stim(t)-3;
            same(t) = (stim(t) == stim(t-1));
            if stim(t) == 2
                evstim(t) = evidence(t-1);%RB:evstim(t) = evidence(t-1);
                match(t) = sign(evidence(t-1));
            else
                evstim(t) = -evidence(t-1);%RB:evstim(t) = -evidence(t-1);
                match(t) = -sign(evidence(t-1));
            end
        end
        if seq{t}.ind > 2
            samepr(t) = (stim(t-2) == stim(t-1));
        else
            samepr(t) = 0;
        end
        ind(t) = index;
        if side(i)
            ipsi(t) = (stim(t)==1);
            evipsi(t) = -evidence(t);
        else
            ipsi(t) = (stim(t)==2);
            evipsi(t) = evidence(t);
        end
    end
    
    if additional ==9 || additional ==10  %this is when we have the 4 models as estimated from behaviour
        x = [-1*same(include)'];%*-1 is becasuse it should be DIFF-SAME;
    else
        x = [-1*same(include)'];%*-1 is becasuse it should be DIFF-SAME; s]; % previously included:  samepr(include)', ipsi(include)'
    end
    
    if additional == 0
        %
    elseif additional == 1
        x = [x, match(include)'];
    elseif additional == 2
        x = [x, evstim(include)'];
    elseif additional == 3
        x = [x, evipsi(include)'];
    elseif additional == 4
        x = [x, abs(evidence(include)')];
    elseif additional == 5
        x = [x, ind(include)'];
    elseif additional == 6
        matchbinary=match;matchbinary(match==-1)=0;
        x = [x, matchbinary(include)'];
    elseif additional == 7
        matchtoresp = stim==choice;
        x = [x, matchtoresp(include)'];
    elseif additional == 8
        matchtoacc = acc;
        x = [x, matchtoacc(include)'];
        
    elseif additional == 9
        
        %                         load('STN_4modelregs.mat');
        load('/home/ezp/Documents/STN/STN_4modelregsNORM.mat');
        if excl
            %remove C43/C46
            STN_M1([6,7],:,:)=[];STN_M2([6,7],:,:)=[];STN_M3([6,7],:,:)=[];STN_M4([6,7],:,:)=[];
        end
        
        for m=1:4
            
            eval(['STN_M = STN_M',num2str(m),';']);
            tempSTN=[];
            for j=1:size(STN_M,2);for jj=1:size(STN_M,3);if isempty(STN_M{i,j,jj})==1;continue;else temp=STN_M{i,j,jj};tempSTN=[tempSTN,temp];end;end;end
            
            if size(tempSTN,2)~=size(seq,2) == 1
                
                fprintf('reducing tempSTN');disp(i);
                
                %SEE fixing_codes.m for details!!!!
                if i==7 && freq~=2 %(LN_C56 trial 1 in last block missing first and last trigs)
                    tempSTN=tempSTN([1:784,786:789,791:end]);
                elseif i==13 %(LN_C62 end missing from D)
                    tempSTN=tempSTN(1:size(seq,2));
                end
                
            else
            end
            
            STN_Mreg=tempSTN(include);
            
            x = [x, STN_Mreg'];
            
        end
        
    end
    
    %zscore regressors
    x=zscore(x);
    
    
    clear same ipsi evstim samepr ind match evidence evipsi matchbinary matchtoresp matchtoacc
    
    if additional ==9
        
        if tailored == 1
            load('behavioral_models.mat','winning_model');winning_model=winning_model;
        else
            winning_model(i)=1;
        end
        
        %1 model at a time
        for m=1:4
            
            xm=x(:,[1,m+1]);
            %                 xm=x(:,[m]);%remove same
            
            for d = 1:length(duration)
                y = lfp(include, duration(d));
                fullGLMdata_LFP{i}{1,d}=y;
                if winning_model(i)==m
                    fullGLMdata_REGS{i}{1}=xm(:,1);
                    fullGLMdata_REGS{i}{2}=xm(:,2);
                end
                B = glmfit(xm,y);
                Bsame(i,d) = B(2);%remove same
                eval(['BM',num2str(m),'(i,d) = B(3);']);
                %                 BM2(i,d) = B(4);
                %                 BM3(i,d) = B(5);
                %                 BM4(i,d) = B(6);
            end
        end
        
    elseif additional ==10
        
        
        for d = 1:length(duration)
            y = lfp(include, duration(d));
            B = glmfit(x,y);
            Bsame(i,d) = B(2);
        end
        
        
    else
        
        for d = 1:length(duration)
            y = lfp(include, duration(d));
            fullGLMdata_LFP{i}{1,d}=y;
            fullGLMdata_REGS{i}{1}=x(:,1);
            fullGLMdata_REGS{i}{2}=x(:,2);
            B = glmfit(x,y);
            Bsame(i,d) = B(2);
            %                 Bsamepr(i,d) = B(3);
            %                 Bipsi(i,d) = B(4);
            if additional > 0
                Badd(i,d) = B(3);
            end
        end
    end
    
    nsame(i) = choice_after_same (seq, 1);
    
    if  additional ==9 && tailored == 1
        
        %            load('behavioral_models.mat','winning_model_no5');winning_model=winning_model_no5;
        load('behavioral_models.mat','winning_model');winning_model=winning_model;
        
        eval(['BMwinner(i,:)= BM',num2str(winning_model(i)),'(i,:);']);
    
        Badd=[];
        
    elseif  additional ==9 && tailored == 0       
        BMwinner(i,:)= BM1(i,:);
        Badd=[];
    end
    
end



%%%PLOTTING
% %     errorbar (data.time, mean(Bsame), std(Bsame)/sqrt(n));
% %     hold on
% %     errorbar (data.time, mean(Bsamepr), std(Bsamepr)/sqrt(n));
% %     errorbar (data.time, mean(Bipsi), std(Bipsi)/sqrt(n));
% %     if additional > 0
% %         errorbar (data.time, mean(Badd), std(Badd)/sqrt(n));
% %     end
% %     legend ('same', 'prevous same', 'ipsi');
% %     xlabel ('Time from cue onset');
% %     ylabel ('Regression coefficient');
% %
% %
%     if last == 0
%         regressor = Bsame;
%         same_time = 15:27;
%     else
%         regressor = Bsamepr;
%         same_time = 10:13;
%     end
% %
% %     figure
%     nsame_pat = nsame(1);
%     bsame_pat = mean(abs(regressor(1,same_time)));
%     patnum = files(1).name(6);
%
%     for i = 2:n
%         if files(i).name(6) == patnum
%             bsame_pat(end) = (bsame_pat(end) + mean(abs(regressor(i,same_time)))) / 2;
%         else
%             nsame_pat = [nsame_pat, nsame(i)];
%             bsame_pat = [bsame_pat, mean(abs(regressor(i,same_time)))];
%             patnum = files(i).name(6);
%         end
%     end
% %
% %     plot (bsame_pat, nsame_pat, 'o');
% %     xlabel ('Regression weight of same stimulus');
% %     ylabel ('Fraction of choices after same stimuli');
%     [r, p] = corr (bsame_pat', nsame_pat');
%     fprintf ('Correlation between the regression weight of same, and fraction of choices of the same stimulus r=%f, p=%f\n', r, p);


%stats setup before plotting
%%%% stats

cfg=[];
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.method ='montecarlo';
cfg.correctm = 'cluster';
cfg.clusteralpha  = 0.05;%cluster forming thhreshold...it's complicated and FSL uses the same, that first you create some clusters, which need to be at a certain threshold, then you do stats on those, which is the values below
cfg.alpha       = 0.025;
cfg.tail        = 0; % two-sided test
cfg.correcttail = 'prob';
cfg.numrandomization = 1000;%the higher the better, 1000 is enough to ensur
cfg.latency = [0 0.8];
cfg.avgovertime = 'no';
cfg.neighbours    =[];

Nsub=size(Bsame,1);

subj = Nsub;
subj2 = Nsub;
design = zeros(2,subj+subj2);
for i = 1:subj
    design(1,i) = i;
end
for i = 1:subj2
    design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:subj+subj2) = 2;

cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;
%

%%%EZP PLOTTING
figure;set(gcf,'color','w');

if last ==0
    title('evidence over all stim', 'FontSize',20);
else
    title('evidence at last stim', 'FontSize',20);
end

if additional == 9
    
    if tailored == 1
        
        p1=plotpatch (Bsame, data.time,[0.5 0.5 0.5]);
        hold on
        p2=plotpatch (BMwinner, data.time,'k');
        
        %%make FT data
        avsame_group.avg=reshape(Bsame,[Nsub, 1, size(data.time,2)]);%avsame_group.avg=Bsame;
        avdiff_group.avg=reshape(BMwinner,[Nsub, 1, size(data.time,2)]);%BMwinner;
        dummy.avg=zeros(size(avsame_group.avg));
        avsame_group.dimord = 'subj_chan_time';avdiff_group.dimord = 'subj_chan_time';dummy.dimord = 'subj_chan_time';
        avsame_group.label={'LFP'};avdiff_group.label={'LFP'};dummy.label={'LFP'};
        avsame_group.time=data.time;avdiff_group.time=data.time;dummy.time=data.time;
        
        
        [stat] = ft_timelockstatistics(cfg, avsame_group, avdiff_group);%
        
%           sig = find (stat.mask==1);%<0.05
        sig = find (stat.prob<0.05/3);%correct for MC
        sig=sig+10;
        plot (data.time(sig), ones(1, length(sig))*0.1, 'r.','MarkerSize',15);
        
        
        %         [h,p] = ttest (Bsame);
        %         % %     [p,praw] = ClusterCorrection2 (Bsame - BMwinner, 500, 0.05);
        %         h(p<0.05)=1;
        %         sigS = find (h);dotsloc=ylim;clear p praw h
        %         plot (data.time(sigS), ones(1, length(sigS))*dotsloc(1), 'k*');
        
        %     [h1] = ttest (Bsame,0);
        [statS] = ft_timelockstatistics(cfg, avsame_group, dummy);%
            sig1 = find (statS.mask==1);sig1=sig1+10;
%         sig1 = find (statS.prob<0.05/3);%correct for MC
        plot (data.time(sig1), ones(1, length(sig1))*0.11, 'k+');
        %     [h2] = ttest (BMwinner,0);%plain ttest against 0
        [statW] = ft_timelockstatistics(cfg, avdiff_group, dummy);%cluster stat against 0
            sig2 = find (statW.mask==1);sig2=sig2+10;
%         sig2 = find (statW.prob<0.05/3);%correct for MC
        plot (data.time(sig2), ones(1, length(sig2))*0.12, 'k.');
        
        lgd=legend ([p1(1) p2(1)],'Box', 'off');
        lgd.String={'same','tailored normalization model'};
        
    else
        
        p1=plotpatch (Bsame, data.time,[0.5 0.5 0.5]);
        p2=plotpatch (BM1, data.time,'k');
        p3=plotpatch (BM2, data.time,'r');
        p4=plotpatch (BM3, data.time,'b');
        p5=plotpatch (BM4, data.time,'c');
        lgd=legend ([p1(1) p2(1) p3(1) p4(1) p5(1) ],'Box', 'off');
        lgd.String={'same','M1','M2','M3','M4'};
        
    end
    
elseif additional ==10
    
    p1=plotpatch (Bsame, data.time,'k');
    [h1] = ttest (Bsame,0);
    sig1 = find (h1);
    plot (data.time(sig1), ones(1, length(sig1))*0.1, 'k+');
    
else
    
    p1=plotpatch (Bsame, data.time,'k');
    hold on
    % p2=plotpatch (Bsamepr, data.time,[0.5 0.5 0.5]);
    % p3=plotpatch (Bipsi, data.time,'b');
    p4=plotpatch (Badd, data.time,'c');
    
    %%make FT data
    avsame_group.avg=reshape(Bsame,[Nsub, 1, size(data.time,2)]);%avsame_group.avg=Bsame;
    avdiff_group.avg=reshape(Badd,[Nsub, 1, size(data.time,2)]);%BMwinner;
    dummy.avg=zeros(size(avsame_group.avg));
    avsame_group.dimord = 'subj_chan_time';avdiff_group.dimord = 'subj_chan_time';dummy.dimord = 'subj_chan_time';
    avsame_group.label={'LFP'};avdiff_group.label={'LFP'};dummy.label={'LFP'};
    avsame_group.time=data.time;avdiff_group.time=data.time;dummy.time=data.time;
    
    [stat] = ft_timelockstatistics(cfg, avsame_group, avdiff_group);%
    
    %         sig = find (stat.mask==1);
    sig = find (stat.prob<0.05/3);%correct for MC
    sig=sig+10;
    plot (data.time(sig), ones(1, length(sig))*0.1, 'r.','MarkerSize',15);
    
    if additional > 0
        % [h] = ttest (Badd);
        % sig = find (h);
        % plot (data.time(sig), ones(1, length(sig))*0.1, 'k.');
        [statS] = ft_timelockstatistics(cfg, avsame_group, dummy);%
            sig1 = find (statS.mask==1);sig1=sig1+10;
%         sig1 = find (statS.prob<0.05/3);%correct for MC
        plot (data.time(sig1), ones(1, length(sig1))*0.11, 'k+');
        [statA] = ft_timelockstatistics(cfg, avdiff_group, dummy);%cluster stat against 0
            sig2 = find (statA.mask==1);sig2=sig2+10;
%         sig2 = find (statA.prob<0.05/3);%correct for MC
        plot (data.time(sig2), ones(1, length(sig2))*0.12, 'k.');
        
    end
    
    lgd=legend ([p1(1)  ],'Box', 'off');%p2(1) p3(1)
    lgd.String={'same'};% 'prevous same', 'ipsi'
    
    if additional > 0 && additional <9
        %     p4=plotpatch (Badd, data.time,'c');
        lgd=legend ([p1(1) p4(1) ],'Box', 'off');% p2(1) p3(1)
        if additional == 1
            bla='match to average (ongoing)';
        elseif additional == 2
            bla='evidence for current';
        elseif additional == 3
            bla='integrated evidence for ipsi';
        elseif additional == 4
            bla='absolute integrated evidence (confidence)';
        elseif additional == 5
            bla='index of stimuli (urgency)';
        elseif additional == 6
            bla='match to average (ongoing) binary';
        elseif additional == 7
            bla='match to final response';
        elseif additional == 8
            bla='accurate?';
        end
        lgd.String={'same', bla};%
    end
    
end

xlabel ('Time from cue onset', 'FontSize',15);
ylabel ('Regression coefficient', 'FontSize',15);
lgd.Location ='northeast';
lgd.FontSize=20;

time=data.time;

end

function nsame = choice_after_same(seq, last)

nsame = 0;
ntrials = 0;
for t = 1:length(seq)
    if seq{t}.ind == 1
        ntrials = ntrials + 1;
        if last
            if length (seq{t}.seq) >= 2
                if seq{t}.seq(end-1) == seq{t}.seq(end)
                    nsame = nsame + 1;
                end
            end
        else
            if length (seq{t}.seq) >= 3
                if seq{t}.seq(end-1) == seq{t}.seq(end-2)
                    nsame = nsame + 1;
                end
            end
        end
    end
end
nsame = nsame / ntrials;
end