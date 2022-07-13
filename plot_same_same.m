
function [ avsame  avdiff betalfp cueind evidence_index avsame2 avsame3  avdiff2 avdiff3  ,sigA, statA avsame_group, avdiff_group, avsame_same,avdiff_same,avsame_diff,avdiff_diff,...
    trlwsesame trlwsediff trlwsesame2 trlwsediff2 trlwsesame3 trlwsediff3 trlnumsame trlnumdiff] = plot_same_same(freq,firsts,accuracy)


% The function compares the beta power after cue depending on whether the
% current and previous stimuli where matchig the preceeding stimuli
%%freq:
% % % 1: beta
% % % 2: theta
%%firsts:
% % % 1: (first few 3 stim only),
% % % 0: (all stim 'i'), -1 (same/diff based on 2nd cue only but plotted for all);...
% % %-1: plot from the penultimate cue backwards
%%accruacy: whether to do on 
% % 1: accurate 
% % 0: inaccuracte
% % 2: alltrials


if freq==1
    files = dir('/home/ezp/Documents/STN/GitHub/BETA_LFP/P*allchans_cue.mat');%all channels with beat a-priori
elseif freq==2
    files = dir('/home/ezp/Documents/STN/GitHub/THETA_LFP/P*allchans_cue.mat');%all channels
end

n = numel(files);
for i = 1:n
    if freq==1
        load (['/home/ezp/Documents/STN/GitHub/BETA_LFP/' files(i).name]);
    elseif freq==2
        load (['/home/ezp/Documents/STN/GitHub/THETA_LFP/' files(i).name]);
    end
    
    same_same = [];
    diff_same = [];
    same_diff = [];
    diff_diff = [];
    same_same_prev = [];
    diff_same_prev = [];
    same_diff_prev = [];
    diff_diff_prev = [];
    include=[];
    for t = 1:length(seq)
        index = seq{t}.ind; %index of cue within a trial
        
        if index > 2 && index < length (seq{t}.seq)   %ignore the first and the last trial
            stim = seq{t}.seq(index) * 2 - 3;
            prev = seq{t}.seq(index-1) * 2 - 3;
            grand = seq{t}.seq(index-2) * 2 - 3;
            if prev == grand
                if stim == prev
                    if isempty(find(bad==t))
                        same_same = [same_same, t];
                        same_same_prev = [same_same_prev, t-1];
                    end
                else
                    if isempty(find(bad==t))
                        same_diff = [same_diff, t];
                        same_diff_prev = [same_diff_prev, t-1];
                        
                    end
                end
            else
                if stim == prev
                    if isempty(find(bad==t))
                        diff_same = [diff_same, t];
                        diff_same_prev = [diff_same_prev, t-1];
                    end
                else
                    if isempty(find(bad==t))
                        diff_diff = [diff_diff, t];
                        diff_diff_prev = [diff_diff_prev, t-1];
                    end
                end
            end
        end
        
        %cues to include in LFP and urgency/evidence (these are for
        %plotting raw data for these regs from GLM)
        if  index > 1 && index < length (seq{t}.seq) && isempty(find(bad==t))  %ignore the first and the last trial
            include = [include, t];
        end
        %for urgency           
        cueindtemp(t)=index;
        %for evidence
        stim = seq{t}.seq(index)*2 - 3;
        if index == 1
            evidence(t) = stim;        
        else
            evidence(t) = evidence(t-1) + stim;
        end
        
            
    end
    
    %for raw data plotting (GLM control)
    cueind{i}= cueindtemp(include);
    evidence_index{i}=evidence(include); 
    
    lfp=squeeze(nanmean(data.trial,2));
    betalfp{i}=lfp(include,:);
    
    same = [];
    diff = [];
    same2 = [];
    diff2 = [];
    same3 = [];
    diff3 = [];
    same4=[];same5=[];same6=[];diff4=[];diff5=[];diff6=[];
    sametrlnum=[];
    difftrlnum=[];

    for t = 1:length(seq)
        index = seq{t}.ind; %index of cue within a trial
        
        %by accuracy
        temp = seq{t}.resp;if strcmp(temp,'Left_correct')==1 || strcmp(temp,'Right_correct')==1 ; acc(t)=1;else acc(t)=0;end
        if accuracy == 1
            accon=1;
        elseif accuracy ==0
            accon=0;
        else
            acc(t)=1;
            accon=1;
        end
        
        if firsts==0
            
            if index > 2 && index < length (seq{t}.seq) %ignore the first and the last trial
                stim = seq{t}.seq(index) * 2 - 3;
                prev = seq{t}.seq(index-1) * 2 - 3;
                grand = seq{t}.seq(index-2) * 2 - 3;
                if prev == grand
                    if isempty(find(bad==t)) && acc(t)==accon
                        same= [same, t-1];
                        same2 = [same2, t];
                        same3 = [same3, t+1];
                        sametrlnum=[sametrlnum, index];
                    end
                else
                    if isempty(find(bad==t, 1)) && acc(t)==accon
                        diff = [diff, t-1];
                        diff2 = [diff2, t];
                        diff3 = [diff3, t+1];
                        difftrlnum=[difftrlnum, index];
                    end
                end
            end
            



        elseif firsts==1
            
            if index == 3  && length(seq{t}.seq)>3 %only take first 4 stimulus as sequnce
                stim = seq{t}.seq(index) * 2 - 3;
                prev = seq{t}.seq(index-1) * 2 - 3;
                grand = seq{t}.seq(index-2) * 2 - 3;
                if prev == grand
                    if isempty(find(bad==t))
                        same= [same, t-2];
                        same2 = [same2, t-1];
                        same3 = [same3, t];
                        sametrlnum=[sametrlnum, seq{t}.trial];
                    end
                else
                    if isempty(find(bad==t))
                        diff = [diff, t-2];
                        diff2 = [diff2, t-1];
                        diff3 = [diff3, t];
                        difftrlnum=[difftrlnum,seq{t}.trial];
                    end
                end
            end
            
            
            
        elseif firsts==-1
            if index == length (seq{t}.seq)-1 && index>2 %penultimate trial only but only for longer than 1 cue trials
                
                stim = seq{t}.seq(index) * 2 - 3;
                prev = seq{t}.seq(index-1) * 2 - 3;
                grand = seq{t}.seq(index-2) * 2 - 3;
                if prev == grand
                    if isempty(find(bad==t))&& acc(t)==accon
                        same= [same, t-1];
                        same2 = [same2, t];
                        same3 = [same3, t+1];
                        sametrlnum=[sametrlnum, index];
                    end
                else
                    if isempty(find(bad==t))&& acc(t)==accon
                        diff = [diff, t-1];
                        diff2 = [diff2, t];
                        diff3 = [diff3, t+1];
                        difftrlnum=[difftrlnum, index];
                    end
                end
            end
            
            
            
        end
    end
    
    missingchan_subj=[8,11];
    
    if i==missingchan_subj(1) %no R01 channel
        lfp1 = nan(size(squeeze (data.trial(:,1,:))));
        lfp2 = squeeze (data.trial(:,1,:));
        lfp3 = squeeze (data.trial(:,2,:));
        lfp4 = squeeze (data.trial(:,3,:));
        lfp5 = squeeze (data.trial(:,4,:));
        lfp6 = squeeze (data.trial(:,5,:));
    elseif i==missingchan_subj(2) %no left channels
        lfp1 = squeeze (data.trial(:,1,:));
        lfp2 = squeeze (data.trial(:,2,:));
        lfp3 = squeeze (data.trial(:,3,:));
        lfp4 = nan(size(lfp1));
        lfp5 = nan(size(lfp1));
        lfp6 = nan(size(lfp1));
    else
        lfp1 = squeeze (data.trial(:,1,:));
        lfp2 = squeeze (data.trial(:,2,:));
        lfp3 = squeeze (data.trial(:,3,:));
        lfp4 = squeeze (data.trial(:,4,:));
        lfp5 = squeeze (data.trial(:,5,:));
        lfp6 = squeeze (data.trial(:,6,:));
    end
    

    for c=1:6
        eval(['avsame (i,c,:) = nanmean(lfp',num2str(c),'(same,:));']);
        eval(['avdiff (i,c,:) = nanmean(lfp',num2str(c),'(diff,:));']);
        eval(['avsame2 (i,c,:) = nanmean(lfp',num2str(c),'(same2,:));']);
        eval(['avdiff2 (i,c,:) = nanmean(lfp',num2str(c),'(diff2,:));']);
        eval(['avsame3 (i,c,:) = nanmean(lfp',num2str(c),'(same3,:));']);
        eval(['avdiff3 (i,c,:) = nanmean(lfp',num2str(c),'(diff3,:));']);
        
        eval(['avsame_same (i,c,:) = nanmean(lfp',num2str(c),'(same_same,:));']);
        eval(['avdiff_same (i,c,:) = nanmean(lfp',num2str(c),'(diff_same,:));']);
        eval(['avsame_diff (i,c,:) = nanmean(lfp',num2str(c),'(same_diff,:));']);
        eval(['avdiff_diff (i,c,:) = nanmean(lfp',num2str(c),'(diff_diff,:));']);
    end
    
    lfp=squeeze(nanmean(data.trial,2));
    trlwsesame{i} = lfp(same,:);
    trlwsediff{i} = lfp(diff,:);
    trlwsesame2{i} = lfp(same2,:);
    trlwsediff2{i} = lfp(diff2,:);
    trlwsesame3{i} = lfp(same3,:);
    trlwsediff3{i} = lfp(diff3,:);
    trlnumsame{i} = sametrlnum;
    trlnumdiff{i} = difftrlnum;
    

    
    
end

%convert the data to dB
avsame=avsame*10/log(10);
avdiff=avdiff*10/log(10);
avsame2=avsame2*10/log(10);
avdiff2=avdiff2*10/log(10);
avsame3=avsame3*10/log(10);
avdiff3=avdiff3*10/log(10);
%

figure;set(gcf,'color','w');
chans={'STN R01';'STN R12';'STN R23';'STN L01';'STN L12';'STN L23'};
combotime=-0.5:0.05:2.4;

for c=1:6
    
   
        subplot(2,3,c)
        p1=plotpatch (cat(2,squeeze(avsame(:,c,1:end-1)),squeeze(avsame2(:,c,11:end-1)),squeeze(avsame3(:,c,11:end))), combotime,'k');
        hold on
        p2=plotpatch (cat(2,squeeze(avdiff(:,c,1:end-1)),squeeze(avdiff2(:,c,11:end-1)),squeeze(avdiff3(:,c,11:end))), combotime,'b');
        title(chans(c));
        if freq==1;ylim([-.3 .3]);else ylim([-.6 .6]);end
        if c==1;ylabel ('dB', 'FontSize',15);end
        set(gca,'FontSize',15);

        %stats settings
        cfg=[];
        cfg.statistic = 'ft_statfun_depsamplesT';
        cfg.method ='montecarlo';
        cfg.correctm = 'cluster';
        cfg.clusterthreshold = 'parametric'; 
        cfg.clusteralpha  = 0.05;
        cfg.alpha       = 0.05;
        cfg.tail        = 0; 
        cfg.correcttail = 'prob';
        cfg.numrandomization = 1000;
        cfg.latency = [];
        cfg.avgovertime = 'no';
        cfg.neighbours    =[];
        cfg.minnbchan =0;
        cfg.spmversion = 'spm12' ;

        Nsub=size(avsame,1);
        
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
        
        
        avsame_group.avg(c,:,:)=cat(2,squeeze(avsame(:,c,1:end-1)),squeeze(avsame2(:,c,11:end-1)),squeeze(avsame3(:,c,11:end)));
        avdiff_group.avg(c,:,:)=cat(2,squeeze(avdiff(:,c,1:end-1)),squeeze(avdiff2(:,c,11:end-1)),squeeze(avdiff3(:,c,11:end)));

end

if freq==1;sgtitle ('Beta power');elseif freq==2; suptitle ('Theta power');end      
set(gca,'FontSize',15);

%make data struc for stats in Fieldtrip
siz = size(avsame_group.avg);
avsame_group.label = {'STN_R01';'STN_R12';'STN_R23';'STN_L01';'STN_L12';'STN_L23'};avdiff_group.label = {'STN_R01';'STN_R12';'STN_R23';'STN_L01';'STN_L12';'STN_L23'};
avsame_group.dimord = 'rpt_chan_time';avdiff_group.dimord = 'rpt_chan_time';
avsame_group.avg = reshape(avsame_group.avg, [siz(2), siz(1), siz(3)]);
avdiff_group.avg = reshape(avdiff_group.avg, [siz(2), siz(1), siz(3)]);
avsame_group.time=combotime;avdiff_group.time=combotime;

[stat] = ft_timelockstatistics(cfg, avsame_group, avdiff_group);%
figure;imagesc(stat.mask); xlabel ('Time from cue onset', 'FontSize',15);ylabel ('channel', 'FontSize',15);
set(gca, 'XTick', combotime,'XTicklabel', mat2str(combotime), 'Yticklabel',chans);


figure;set(gcf,'color','w');
p1=plotpatch (cat(2,squeeze(mean(avsame(:,:,1:end-1),2)),squeeze(mean(avsame2(:,:,11:end-1),2)),squeeze(mean(avsame3(:,:,11:end),2))), combotime,'k');
hold on
p2=plotpatch (cat(2,squeeze(mean(avdiff(:,:,1:end-1),2)),squeeze(mean(avdiff2(:,:,11:end-1),2)),squeeze(mean(avdiff3(:,:,11:end),2))), combotime,'b');
if freq==1;ylim([-.2 .2]);else ylim([-.6 .6]);end
xlim([-0.5 2.4])
% xline(0,'k');xline(0.8,'k');xline(1.6,'k');
if freq==1;ylabel ('Beta power (dB)', 'FontSize',15);elseif freq==2; ylabel ('THETA  power (dB)', 'FontSize',15);end
set(gca,'FontSize',15);

cfg.avgoverchan = 'yes';
[statA] = ft_timelockstatistics(cfg, avsame_group, avdiff_group);%
sigA = find (statA.mask==1);
disp('significant times')
disp(combotime(sigA))
plot (combotime(sigA), ones(1, length(sigA))*max(ylim), 'k.');

%% cat so can plot difference wave with errors bars

SAME=cat(3,avsame(:,:,1:end-1),avsame2(:,:,11:end-1),avsame3(:,:,11:end));
DIFF=cat(3,avdiff(:,:,1:end-1),avdiff2(:,:,11:end-1),avdiff3(:,:,11:end));

SvsD_timecourse=DIFF-SAME;

figure;set(gcf,'color','w');
combotime=-0.5:0.05:2.4;
p1=plotpatch (squeeze(nanmean(SvsD_timecourse(:,:,:),2)), combotime,'k');

xlabel ('Time from cue onset', 'FontSize',15);
ylabel ('Different-Same (+sem)', 'FontSize',15);
xlim([-0.5 2.4])

% xline(0,'k');xline(0.8,'k');xline(1.6,'k');


x = combotime;                                          % Create Independent Variable
y = squeeze(nanmean(SvsD_timecourse(:,:,:),2));                                  % Create Dependent Variable ‘Experiments’ Data
N = size(y,1);                                      % Number of ‘Experiments’ In Data Set
yMean = mean(y,1);                                    % Mean Of All Experiments At Each Value Of ‘x’
ySEM = nanstd(y,'',1)/sqrt(N);                              % Compute ‘Standard Error Of The Mean’ Of All Experiments At Each Value Of ‘x’
CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
yCI95 = bsxfun(@times, ySEM, CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’

figure;set(gcf,'color','w');
plot(x, yMean,'k','LineWidth',2)                                      % Plot Mean Of All Experiments
hold on
plot(x, yCI95+yMean,'k--','LineWidth',2)                                % Plot 95% Confidence Intervals Of All Experiments
hold off
xlabel ('Time from cue onset', 'FontSize',15);
ylabel ('BETA Diff-Same (+95%CI)', 'FontSize',15);
xlim([-0.5 2.4])
% xline(0,'k');xline(0.8,'k');xline(1.6,'k');
set(gca,'FontSize',15)

%%AT CUE i+1
figure;set(gcf,'color','w');

p1=plotpatch (squeeze(nanmean(avsame_same,2)), data.time,'k');
hold on
p2=plotpatch (squeeze(nanmean(avsame_diff,2)), data.time,[0.5 0.5 0.5]);
p3=plotpatch (squeeze(nanmean(avdiff_same,2)), data.time,'b');
p4=plotpatch (squeeze(nanmean(avdiff_diff,2)), data.time,'c');
xlabel ('Time from cue onset', 'FontSize',15);
if freq==1;ylabel ('Beta power', 'FontSize',15);elseif freq==2; ylabel ('Theta power', 'FontSize',15);end
set(gca,'FontSize',15)
if freq==1;ylim([-.1 .1]);else ylim([-.15 .15]);end
% xline(0,'k');
lgd=legend ([p1(1) p2(1) p3(1) p4(1)]);
lgd.String={'Same-Same', 'Same-Different', 'Different-Same', 'Different-Different'};
lgd.Location ='northeast';
lgd.FontSize=15;
legend boxoff  

end
