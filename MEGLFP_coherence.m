%% coherence b/w MEG SOURCE & STN
clear
%%input type
TF=1;%compute time-frequency coherence
%%
cd D:\GitHub\SOURCE\

TFmeth=1;%1:wavelet (matches Fig4), 2:fixed window multitaper (as used for LFP analysis)
%%
for s = 1:13
    
    eval(['load(''P',num2str(s,'%02d'),'_raw_source_timecourse.mat'')']);
    
    %%%%conditions
    same = [];same_1=[];same_0=[];same_pen=[];
    diff = [];diff_1 = [];diff_0 = [];diff_pen=[];
    bad = [];    % already removed
    for t = 1:length(seq)
        
        %by accuracy
        temp = seq{t}.resp;if strcmp(temp,'Left_correct')==1 || strcmp(temp,'Right_correct')==1 ; acc(t)=1;else acc(t)=0;end
        
        index = seq{t}.ind; %index of cue within a trial
        if index > 1 && index < length (seq{t}.seq)   %ignore the first and the last trial
            if seq{t}.seq(index) == seq{t}.seq(index - 1)
                if isempty(find(bad==t))
                    same = [same, t];
                end
                if isempty(find(bad==t)) && acc(t)==1
                    same_1 = [same_1, t];
                end
                if isempty(find(bad==t)) && acc(t)==0
                    same_0 = [same_0, t];
                end
            else
                if isempty(find(bad==t))
                    diff = [diff, t];
                end
                if isempty(find(bad==t)) && acc(t)==1
                    diff_1 = [diff_1, t];
                end
                if isempty(find(bad==t)) && acc(t)==0
                    diff_0 = [diff_0, t];
                end
            end
        end
        
    end
    %%%%
    
    %get only some segments of the data
    cfg=[];
    if TF==0;cfg.toilim = [0 0.2];end%the sig coh segment 0-200ms
    cfg.trials =same;
    datasame=ft_redefinetrial(cfg,data);
    cfg=[];
    if TF==0;cfg.toilim = [0 0.2];end%the sig coh segment 0-200ms
    cfg.trials =diff;
    datadiff=ft_redefinetrial(cfg,data);
    
    %frequency decomp with phase
    fsample=1/(data.time(2)-data.time(1));
    cfg = [];
    cfg.output ='powandcsd';
    cfg.keeptrials = 'no';
    cfg.keeptapers='no';
    cfg.taper = 'dpss';
    cfg.method  ='mtmconvol';
    
    if TFmeth==1
        
        cfg.foi        = 1:75;
        cfg.tapsmofrq  = 0.2*cfg.foi; %smoothing
        cfg.toi        = -0.5:1/fsample:0.8; %time window
        cfg.t_ftimwin  =  5./cfg.foi;   % number of cycles per timewindow
        
    elseif TFmeth==2
        cfg.foi        = 1:75;
        res  = 2.5*ones(size(cfg.foi));
        res(cfg.foi>25) = 0.1*cfg.foi(cfg.foi>25);
        res(cfg.foi>50) = 5;
        cfg.tapsmofrq  = res; %smoothing
        cfg.toi        = -0.5:1/fsample:0.8; %time window
        cfg.t_ftimwin  =  0.4*ones(size(cfg.foi));
    end
    
    cfg.channel      = 'all';
    cfg.channelcmb      =  'all';
    cfg.pad = 2;
    freqS = ft_freqanalysis(cfg, datasame);
    freqD = ft_freqanalysis(cfg, datadiff);
    
    
    %connectivity
    cfg = [];
    cfg.method  = 'coh';
    cohS = ft_connectivityanalysis(cfg, freqS);
    cohD = ft_connectivityanalysis(cfg, freqD);
    
    group_cohS{s}=cohS;
    group_cohD{s}=cohD;
    
    
    group_freqS{s}=freqS;
    group_freqD{s}=freqD;
    
    
    clearvars -except s subjects group* seqinfofiles incl TFmeth TF
    
end

%% group (specify which STN channels to use and combine double vertices)

TF=1;
granger=0;

for s=1:13
    
    ref = {'STN_R01','STN_R12','STN_R23'};%'STN_L01','STN_L12','STN_L23'};
    %         if s==11;continue;end%for L chans need to remove this subject
    
    for cc=1:size(ref,2)
        temp=find(strcmp(group_cohS{s}.labelcmb(:,1),ref(cc)));
        if isempty(temp)==0;ch(cc,:)=temp(1:6);else;end%only get ref to cortex
    end
    if s==8 && granger ==0;ch(1,:)=[];end%for R chans
    
    %SAME
    group_cohS_red{s}=group_cohS{s};
    %DIFF
    group_cohD_red{s}=group_cohD{s};
    
    group_cohS_red{s}.dimord='chan_freq_time';
    group_cohS_red{s}.label={'occ','temp','frontal'};
    group_cohD_red{s}.dimord='chan_freq_time';
    group_cohD_red{s}.label={'occ','temp','frontal'};
    
    
    if size(ref,2)>1
        siz=size(group_cohS_red{s}.cohspctrm);
        group_cohS_red{s}.cohspctrm=nan(3,siz(2),siz(3));group_cohD_red{s}.cohspctrm=nan(3,siz(2),siz(3));%avg over sources
        k=1;
        for cc=1:3 %avg over refs and double-vertices
            group_cohS_red{s}.cohspctrm(cc,:,:)=nanmean(group_cohS{s}.cohspctrm(ch(:,k:k+1),:,:));
            group_cohD_red{s}.cohspctrm(cc,:,:)=nanmean(group_cohD{s}.cohspctrm(ch(:,k:k+1),:,:));k=k+2;
        end
    else
        group_cohS_red{s}.cohspctrm=group_cohS{s}.cohspctrm(ch,:,:);
        group_cohD_red{s}.cohspctrm=group_cohD{s}.cohspctrm(ch,:,:);
    end
end

clear ch temp*

group_cohS_red=group_cohS_red(~cellfun('isempty',group_cohS_red));
group_cohD_red=group_cohD_red(~cellfun('isempty',group_cohD_red));


%% data for stats
group_cohSvsD=group_cohS_red;group_cohDum=group_cohD_red;
for s=1:size(group_cohSvsD,2)
    
    group_cohSvsD{s}.cohspctrm = atanh(group_cohD_red{s}.cohspctrm) - atanh(group_cohS_red{s}.cohspctrm);%%DvsS
    group_cohDum{s}.cohspctrm(:)=0;
    
end
%% stats
cfg=[];
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.method ='montecarlo';
cfg.correctm = 'cluster';
cfg.clusteralpha  = 0.05;
cfg.alpha       = 0.05;
cfg.parameter = 'cohspctrm';
cfg.numrandomization = 1000;
cfg.latency = [0 0.8];
cfg.frequency = [1 30];
cfg.neighbours       = [];
cfg.avgovertime='no';
cfg.avgoverfreq='no';

Nsub=length(group_cohD_red);

%design matrix
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

[stat] = ft_freqstatistics(cfg, group_cohSvsD{:}, group_cohDum{:});%zscore diff vs dummy

% % %find significant times/freqs
k=1;
for i=1:size(stat.mask,2)
    if isempty(find(squeeze(stat.mask(:,i,:))==1))==0
        sigfreqCOH(k)=stat.freq(i);sigfreqIND(k)=i;k=k+1;
    else
    end
end
k=1;
for i=1:size(stat.mask,1)
    if isempty(find(squeeze(stat.mask(i,:,:))==1))==0
        sigchCOH(k)=stat.label(i);sigchIND(k)=i;k=k+1;
    else
    end
end
k=1;
for i=1:size(stat.mask,3)
    if isempty(find(squeeze(stat.mask(:,:,i))==1))==0
        sigtimeCOH(k)=stat.time(i);sigtimeIND(k)=i;k=k+1;
    else
    end
end


cfg=[];
cfg.parameter       = 'cohspctrm';
group_cohSvsD_avg=ft_freqgrandaverage(cfg,group_cohSvsD{:});
group_cohS_avg=ft_freqgrandaverage(cfg,group_cohS_red{:});
group_cohD_avg=ft_freqgrandaverage(cfg,group_cohD_red{:});

%plot COH TF with stat mask (rather than plotting stat as above)
cfg = [];
cfg.zlim=[-0.05 0.05];
cfg.maskparameter = 'mask';
cfg.maskstyle = 'outline';
cfg.parameter = 'cohspctrm';

siz=size(group_cohSvsD_avg.cohspctrm);
tempmask=padarray(stat.mask,[0,stat.freq(1)-1,siz(3)-size(stat.mask,3)],'pre');
tempmask=padarray(tempmask,[0,siz(2)-size(stat.mask,2),0],'post');
group_cohSvsD_avg.mask=logical(tempmask>0.5);

figure;set(gcf,'color','w');ylabel('frequency');xlabel('time');set(gca,'FontSize',15);
labels={'occipital ','temporal','frontal'};
for i=1:3
    cfg.channel=i;
    subplot(3,1,i);
    ft_singleplotTFR(cfg, group_cohSvsD_avg);
    ylabel ('coherence frequency', 'FontSize',15);xlabel('time','FontSize',15);
    if length(ref)==1
        title([group_cohSvsD_avg.label(i),' <> ',ref],'FontSize',15);
    else
        title([labels{i},'<>right STN'],'FontSize',15);
    end
end

%% plot TF effect over different conditions
% load('sigcoh.mat')
n=13;

sigtimeIND_adj=sigtimeIND+76;%stats not run on precue timepoints so starting from time zero;

for s=1:n
              
        S_spectra(s,:)=nanmean(squeeze(nanmean(group_cohS_red{s}.cohspctrm(sigchIND,:,sigtimeIND_adj),1)),2);
        D_spectra(s,:)=nanmean(squeeze(nanmean(group_cohD_red{s}.cohspctrm(sigchIND,:,sigtimeIND_adj),1)),2);
        
        S_betacoh_timecourse(s,:)=nanmean(squeeze(nanmean(group_cohS_red{s}.cohspctrm(sigchIND,sigfreqCOH,:),1)),1);
        D_betacoh_timecourse(s,:)=nanmean(squeeze(nanmean(group_cohD_red{s}.cohspctrm(sigchIND,sigfreqCOH,:),1)),1)        
  
    clear temp*
end

figure;set(gcf,'color','w');
subplot(1,2,1);
plot(nanmean(S_spectra),'k','linewidth',5);hold on;plot(nanmean(D_spectra),'b','linewidth',5);ylim([0 0.2])%xticks(1:5:75);
ylabel('coherence');xlabel('frequency');set(gca,'Fontsize',15);
subplot(1,2,2);
plotpatch(S_betacoh_timecourse,group_cohD{1}.time,'k');hold on;plotpatch(D_betacoh_timecourse,group_cohD{1}.time,'b');
ylabel('coherence');xlabel('time');set(gca,'Fontsize',15);
suptitle('right STN <> right frontal');set(gca,'Fontsize',15);

%  with confidence intervals
x = group_cohD{1}.time;
y = D_betacoh_timecourse-S_betacoh_timecourse;
N = size(y,1);
yMean = mean(y,1);
ySEM = nanstd(y,'',1)/sqrt(N);
CI95 = tinv([0.025 0.975], N-1);
yCI95 = bsxfun(@times, ySEM, CI95(:));

figure;set(gcf,'color','w');
plot(x, yMean,'k','LineWidth',2)
hold on
plot(x, yCI95+yMean,'k--','LineWidth',2)
hold off
xlabel ('Time from cue onset', 'FontSize',15);
ylabel ('BETA coherence Diff-Same (+95%CI)', 'FontSize',15);
xlim([-0.5 0.8])
xline(0,'k');xline(0.8,'k');yline(0,'k');


% % %individual
figure;set(gcf,'color','w');
for s=1:13
    subplot(4,4,s);
    plot(S_spectra(s,:),'k','linewidth',3);hold on;plot(D_spectra(s,:),'b','linewidth',3);xticks(2:10:75);ylim([0 0.2])
    ylabel('coherence');xlabel('frequency');
end
suptitle('right STN <> frontal');
