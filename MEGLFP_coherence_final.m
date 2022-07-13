%% coherence b/w MEG SOURCE & STN
clear
 
addpath('/home/ezp/Documents/MATLAB/fieldtrip-20200914/');
%%
cd /home/ezp/Documents/STN/GitHub/SOURCE

TFmeth=2;

P=1;
%1:wavelet 
%2:fixed window multitaper (as used for LFP/MEG analysis); Figure 4 inpaper

COH=2;%1:coh;2:dwpli
%%
if P==0
subjects = {'MG4492','MG4493','MG4521','MG4551','MG4552','MG4554','MG4558','MG4562','MG4563','MG06109','MG06114','MG06117','MG06126'};
else
subjects = {'LN_C31','LN_C35','LN_C37','LN_C38','LN_C39','LN_C48','LN_C56','LN_C57','LN_C58','LN_C59','LN_C60','LN_C61','LN_C62'};...
end

for s = 1:13
   
    if TFmeth ==3
        %        eval(['load(''',subj{s},'_tfcoh_timecourse_BF.mat'')']);%TF and COH on con't data
        %        eval(['load(''',subj{s},'_tf_timecourse_BF.mat'')']);%TF  on con't data
        eval(['load(''',subjects{s},'_tf_fixed_ds_long_timecourse_BF_cluster1_2_CONJ_P_FWE.mat'')']);%TF  on con't data (match LFP/MEG)
    else
        %     eval(['load(''P',num2str(s,'%02d'),'_raw_source_timecourse.mat'')']);%epoched, raw, longer than what we need: -1 to 1s, for TF windowing (epoch is -500 to 800ms around cue); based on beta effect in source (combined P+C) based on control sensor effect
        %     eval(['load(''',subjects{s},'_cluster1_2_CONJ_P_FWE_raw_long_timecourse_BF.mat'')']);%epoched, raw, longer than what we need: -1 to 1.6s, for TF windowing and up to cue i+1
        eval(['load(''',subjects{s},'_notsss_cluster2_raw_long_timecourse_BF.mat'')']);%epoched, raw, longer than what we need: -1 to 1.6s, for TF windowing and up to cue i+1
        %     eval(['load(''',subjects{s},'_notsss_BA6_raw_long_timecourse_BF.mat'')']);%epoched, raw, longer than what we need: -1 to 1.6s, for TF windowing and up to cue i+1
        %     eval(['load(''',subjects{s},'_cluster1_CONJ_P_FWEBA6_raw_long_timecourse_BF.mat'')']);%epoched, raw, longer than what we need: -1 to 1.6s, for TF windowing and up to cue i+1
        %     eval(['load(''',subjects{s},'_uncorrpeaks_raw_timecourse_BF.mat'')']);%nch=5;based  on beta effect in source in patients based on control sensor effect
        %     eval(['load(''',subjects{s},'_fwethetapatientstats_raw_long_timecourse_BF.mat'')']);  %nch=4;based  on theta effect in source in patients based on patient sensor effect
        %     eval(['load(''',subjects{s},'_fwe_combined_thetapatientstats_raw_long_timecourse_BF.mat'')']);  %nch=7;based  on theta effect in source in combined group based on patient sensor effect
    end
    
    %rename cortical sources
    data.label{1}='BA7';
    data.label{2}='BA23';
    data.label{3}='BA123';
    data.label{4}='BA6';
    
    %%%%sort conditions based on cue sequence
    same = [];same_1=[];same_0=[];same_pen=[];
    diff = [];diff_1 = [];diff_0 = [];diff_pen=[];
    bad = [];    % already removed
    for t = 1:length(seq)
        
        %by accuracy
        temp = seq{t}.resp;if strcmp(temp,'Left_correct')==1 || strcmp(temp,'Right_correct')==1 ; acc(t)=1;else acc(t)=0;end
        
        index = seq{t}.ind; %index of cue within a trial
        if index > 1 && index < length (seq{t}.seq)  %ignore the first and the last trial
%         if index > 1 && index == length (seq{t}.seq) -1   % penultimate
%         if index > 1 && index == 2   % random cue number

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
     
    if TFmeth==1 || TFmeth==2
        
        %frequency decomp with phase
        %     fsample=1/(data.time(2)-data.time(1));%make sure raw data in FT format (dimord: 'rpt_chan_time')
        fsample=1/(data.time{1}(2)-data.time{1}(1));%make sure raw data in FT format (dimord: 'rpt_chan_time')
        cfg = [];
        cfg.output ='powandcsd';
        cfg.keeptapers='no';
        cfg.taper = 'dpss';
        cfg.method  ='mtmconvol';
        cfg.pad = 'nextpow2';

        if COH==1 

        cfg.keeptrials = 'no';
        cfg.channel      = 'all';
        cfg.channelcmb   = 'all';

        elseif COH==2 %only take relevant channel combos otherwise memory issue
            
        cfg.keeptrials = 'yes';
        cfg.channel      = 'all';
        cfg.channelcmb =  {
%             'STN_R01' 'BA7'
%             'STN_R12' 'BA7'
%             'STN_R23' 'BA7'
%             'STN_R01' 'BA23'
%             'STN_R12' 'BA23'
%             'STN_R23' 'BA23'
%             'STN_R01' 'BA123'
%             'STN_R12' 'BA123'
%             'STN_R23' 'BA123'
            'STN_R01' 'BA6'
            'STN_R12' 'BA6'
            'STN_R23' 'BA6'};
        end
   
        
        if TFmeth==1
            
            cfg.foi        = 1:45;
            cfg.tapsmofrq  = 0.2*cfg.foi; %smoothing
            cfg.toi        = -0.5:1/fsample:0.8; %time window
            cfg.t_ftimwin  =  5./cfg.foi;   % number of cycles per timewindow
            
        elseif TFmeth==2
            cfg.foi        = 1:45;%45 but 30 for memory help in COH 2
            res  = 2.5*ones(size(cfg.foi));
            cfg.tapsmofrq  = res; %smoothing
            cfg.toi        = -0.5:1/fsample:1.6; %time window
            cfg.t_ftimwin  =  0.4*ones(size(cfg.foi));
        end
                
        
        cfg.trials =same;
        freqS = ft_freqanalysis(cfg, data);
        
        cfg.trials =[];cfg.trials =diff;
        freqD = ft_freqanalysis(cfg, data);       
        
        group_freqS_raw{s}=freqS;
        group_freqD_raw{s}=freqD;
        
        if COH==1
        %ITPC
        cfg.output ='fourier';
        cfg.keeptrials = 'yes';
        cfg.keeptapers = 'yes';
        cfg.trials =[];cfg.trials =same;
        freqS_itpc = ft_freqanalysis(cfg, data);
        
        cfg.trials =[];cfg.trials =diff;
        freqD_itpc = ft_freqanalysis(cfg, data);
        
            % compute inter-trial phase coherence (itpc)
            F = freqS_itpc.fourierspctrm;   % copy the Fourier spectrum
            N = size(F,1);           % number of trials
            itc           = [];
            itc.label     = freqS_itpc.label;
            itc.freq      = freqS_itpc.freq;
            itc.time      = freqS_itpc.time;
            itc.dimord    = 'chan_freq_time';
            itc.itpc      = F./abs(F);         % divide by amplitude
            itc.itpc      = sum(itc.itpc,1);   % sum angles
            itc.itpc      = abs(itc.itpc)/N;   % take the absolute value and normalize
            itc.itpc      = squeeze(itc.itpc); % remove the first singleton dimension
        
            group_freqS_itpc{s}=itc;
            
            F = freqD_itpc.fourierspctrm;   % copy the Fourier spectrum
            N = size(F,1);           % number of trials
            itc           = [];
            itc.label     = freqD_itpc.label;
            itc.freq      = freqD_itpc.freq;
            itc.time      = freqD_itpc.time;
            itc.dimord    = 'chan_freq_time';
            itc.itpc      = F./abs(F);         % divide by amplitude
            itc.itpc      = sum(itc.itpc,1);   % sum angles
            itc.itpc      = abs(itc.itpc)/N;   % take the absolute value and normalize
            itc.itpc      = squeeze(itc.itpc); % remove the first singleton dimension
        
            group_freqD_itpc{s}=itc;
          
        end
        
        %connectivity
        cfg = [];
        
        if COH==1
            
        cfg.method  = 'coh';
        cfg.complex = 'abs';%'imag' doens't look good but requested at some point
        
        elseif COH==2
            
        cfg.method  = 'wpli_debiased';%need trialwise data (request:JoN reviewer)
        
        end   
        cohS = ft_connectivityanalysis(cfg, freqS);
        cohD = ft_connectivityanalysis(cfg, freqD);
        
        group_cohS{s}=cohS;
        group_cohD{s}=cohD;        

        %TF rescale
        cfg = [];
        cfg.baselinetype ='db';
        cfg.baseline = [-0.5 1.6];%before did whole epoch
        freqSlog = ft_freqbaseline(cfg, freqS);
        freqDlog = ft_freqbaseline(cfg, freqD);
        
        group_freqS{s}=freqSlog;
        group_freqD{s}=freqDlog;
        

        
    else
        
        
        datafreq_same=data;
        datafreq_same.powspctrm=squeeze(mean(data.powspctrm(same,:,:,:),1));
        datafreq_diff=data;
        datafreq_diff.powspctrm=squeeze(mean(data.powspctrm(diff,:,:,:),1));
        
        datafreq_same.dimord =  'chan_freq_time';
        datafreq_diff.dimord =  'chan_freq_time';
        %         datafreq_same.time = [-0.45:0.05:0.75];
        %         datafreq_diff.time = [-0.45:0.05:0.75];
        group_freqS{s}=datafreq_same;
        group_freqD{s}=datafreq_diff;
        %     group_cohS{s}=datacoh_same;
        %     group_cohD{s}=datacoh_diff;
        
    end
    
    clearvars -except s subjects group* seqinfofiles incl TFmeth TF subj COH P
    
end


%% group (specify which STN channels to use and combine double vertices)
% 

for s=1:13
    
    ref = {'STN_R01','STN_R12','STN_R23'};%{'STN_L01','STN_L12','STN_L23'};
%     if s==11;continue;end%for L chans need to remove this subject
    
    nch=1;%cortical sources
    for cc=1:size(ref,2)
        temp=find(strcmp(group_cohS{s}.labelcmb(:,1),ref(cc)));
        if isempty(temp)==0;ch(cc,:)=temp(1:nch);else;end%only get ref to cortex
%         if isempty(temp)==0;ch(cc,:)=temp([1,6]);else;end%only occ and ba6 
    end
    if s==8; ch(1,:)=[];end%for R chans missing R01 

%     if s==8;ch=[58:63];else ch=[66:68,70:75];end%only L & R STN coh
       
    %SAME
    group_cohS_red{s}=group_cohS{s};
    %DIFF
    group_cohD_red{s}=group_cohD{s};
    
    group_cohS_red{s}.dimord='chan_freq_time';    
    group_cohS_red{s}.label={'BA6'};%{'BA7','BA23','BA123','BA6'};%{'BA6'};%{'BA23_L','BA6_L','BA47_L','BA23_R','BA6_R','BA6_RD','BA47_R'};%{'BA23','BA6','BA47'};%{'BA23','BA1.2.3','BA4','BA6_M','BA6_L'};%{'BA7_L','BA7_R','BA39'};%{'BA23','BA1.2.3','BA4','BA6_M','BA6_L'};%{'OCC','BA6'};%{'STN_L'};%{'occ','temp','frontal'};{'BA7','BA40','BA4','BA9'};%{'occ','temp','frontal'};%{'BA7','BA40','BA4','BA9'};%{'BA8','OCC','BA40','BA6','IFG'};%{'occ','temp','frontal'};
    group_cohS_red{s}.time=group_cohS{s}.time;
    group_cohD_red{s}.dimord='chan_freq_time';
    group_cohD_red{s}.label={'BA6'};%{'BA7','BA23','BA123','BA6'};%{'BA6'};%{'BA23_L','BA6_L','BA47_L','BA23_R','BA6_R','BA6_RD','BA47_R'};%{'BA7','BA23','BA40','BA6'};%{'BA23','BA6','BA47'};%{'BA7','BA23','BA40','BA6'};%%{'BA7_L','BA7_R','BA39'};%{'BA23','BA1.2.3','BA4','BA6_M','BA6_L'};%{'OCC','BA6'};%{'STN_L'};%{'BA23','BA1.2.3','BA4','BA6_M','BA6_L'};%{'OCC','BA6'};%{'STN_L'};%{'BA23','BA1.2.3','BA4','BA6_M','BA6_L'};%{'OCC','BA6'};%{'BA23','BA1.2.3','BA4','BA6_M','BA6_L'};%{'occ','temp','frontal'};%{'BA23','BA40','BA1.2.3','BA4','BA4_2','BA6'};%{'occ','temp','frontal'};%{'BA23','BA40','BA1.2.3','BA4','BA4_2','BA6'};%{'occ','temp','frontal'};%{'BA7','BA40','BA4','BA9'};%{'occ','temp','frontal'};%{'BA7','BA40','BA4','BA9'};%{'BA8','OCC','BA40','BA6','IFG'};%{'occ','temp','frontal'};
    group_cohD_red{s}.time=group_cohS{s}.time;
    
%     if size(ref,2)>1
%       nch=1;%avg over BA6
        if COH==1
            siz=size(group_cohS_red{s}.cohspctrm);
        elseif COH==2
            siz=size(group_cohS_red{s}.wpli_debiasedspctrm);
        end
        group_cohS_red{s}.cohspctrm=nan(nch,siz(2),siz(3));
        group_cohD_red{s}.cohspctrm=nan(nch,siz(2),siz(3));
        group_cohS_red{s}.wpli_debiasedspctrm=nan(nch,siz(2),siz(3));
        group_cohD_red{s}.wpli_debiasedspctrm=nan(nch,siz(2),siz(3));

%         k=1;%nch=3;
%         for cc=1:nch %avg over refs and double-vertices
%             group_cohS_red{s}.cohspctrm(cc,:,:)=nanmean(group_cohS{s}.cohspctrm(ch(:,k:k+1),:,:));
%             group_cohD_red{s}.cohspctrm(cc,:,:)=nanmean(group_cohD{s}.cohspctrm(ch(:,k:k+1),:,:));k=k+2;
%         end
%       
        if COH==1
            for cc=1:nch %avg over refs only (no doubles)
                group_cohS_red{s}.cohspctrm(cc,:,:)=nanmean(group_cohS{s}.cohspctrm(ch(:,cc),:,:));
                group_cohD_red{s}.cohspctrm(cc,:,:)=nanmean(group_cohD{s}.cohspctrm(ch(:,cc),:,:));
            end
        elseif COH==2
            for cc=1:nch %avg over refs only (no doubles)
                group_cohS_red{s}.wpli_debiasedspctrm(cc,:,:)=nanmean(group_cohS{s}.wpli_debiasedspctrm(ch,:,:));
                group_cohD_red{s}.wpli_debiasedspctrm(cc,:,:)=nanmean(group_cohD{s}.wpli_debiasedspctrm(ch,:,:));
            end
        end
%         for cc=1:nch %avg over all refs (within STN)
%             group_cohS_red{s}.cohspctrm(cc,:,:)=nanmean(group_cohS{s}.cohspctrm(ch,:,:));
%             group_cohD_red{s}.cohspctrm(cc,:,:)=nanmean(group_cohD{s}.cohspctrm(ch,:,:));
%         end
%     else
%         group_cohS_red{s}.cohspctrm=group_cohS{s}.cohspctrm(ch,:,:);
%         group_cohD_red{s}.cohspctrm=group_cohD{s}.cohspctrm(ch,:,:);
%     end
end

clear ch temp*

if mean(strcmp(ref,'STN_R01'))==1
group_cohS_red{8}=[];
group_cohD_red{8}=[];
end

group_cohS_red=group_cohS_red(~cellfun('isempty',group_cohS_red));
group_cohD_red=group_cohD_red(~cellfun('isempty',group_cohD_red));


%% data for stats
group_cohSvsD=group_cohS_red;group_cohDum=group_cohD_red;
for s=1:size(group_cohSvsD,2)
    if COH==1
    group_cohSvsD{s}.cohspctrm = atanh(group_cohD_red{s}.cohspctrm) - atanh(group_cohS_red{s}.cohspctrm);%%DvsS
    group_cohDum{s}.cohspctrm(:)=0;
%     %or could just compare atans
%     group_cohSvsD{s}.cohspctrm = atanh(group_cohD_red{s}.cohspctrm);%D
%     group_cohDum{s}.cohspctrm = atanh(group_cohS_red{s}.cohspctrm);%S
    
    elseif COH==2
    group_cohSvsD{s}.wpli_debiasedspctrm = atanh(group_cohD_red{s}.wpli_debiasedspctrm) - atanh(group_cohS_red{s}.wpli_debiasedspctrm);%%DvsS
    group_cohDum{s}.wpli_debiasedspctrm(:)=0;
    end
    
end
%% stats
cfg=[];
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.method ='montecarlo';
cfg.correctm = 'cluster';
cfg.clusteralpha  = 0.05;
cfg.alpha       = 0.05;
cfg.tail        = 0;%0 default
if COH==1
cfg.parameter = 'cohspctrm';
elseif COH==2
cfg.parameter = 'wpli_debiasedspctrm';
end
cfg.numrandomization = 1000;
cfg.latency = [0 1.6];
cfg.frequency = [1 30];
cfg.neighbours       = [];
cfg.avgovertime='no';
cfg.avgoverfreq='no';
cfg.spmversion = 'spm12' ;
% cfg.channel='BA6';
% cfg.channel={'BA23_L','BA23_R','BA6_R','BA47_R'};%CONJ FWE areas only

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

for i=1:size(stat.mask,1)
    figure;set(gcf,'color','w');ylabel('frequency');xlabel('time');
    imagesc(squeeze(stat.prob(i,:,:)));
    set(gca,'YDir','normal','XTickLabel',round(stat.time(1:31:end),1),'XTick',[1:31:240])
end
 xlabel ('Time from cue onset (s)', 'FontSize',15);
 ylabel ('Frequency (Hz)', 'FontSize',15);
%%
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
%% 

cfg=[];
if COH==1
    cfg.parameter       = 'cohspctrm';
elseif COH==2
    cfg.parameter       = 'wpli_debiasedspctrm';
end
group_cohSvsD_avg=ft_freqgrandaverage(cfg,group_cohSvsD{:});
group_cohS_avg=ft_freqgrandaverage(cfg,group_cohS_red{:});
group_cohD_avg=ft_freqgrandaverage(cfg,group_cohD_red{:});

%plot COH TF with stat mask (rather than plotting stat as above)
cfg = [];
cfg.zlim=[-0.05 0.05];
cfg.maskparameter = 'mask';
cfg.maskstyle = 'outline';
if COH==1
    cfg.parameter       = 'cohspctrm';
    siz=size(group_cohSvsD_avg.cohspctrm);
elseif COH==2
    cfg.parameter       = 'wpli_debiasedspctrm';
    siz=size(group_cohSvsD_avg.wpli_debiasedspctrm);
end

tempmask=padarray(stat.mask,[0,0,find(group_cohSvsD_avg.time==stat.time(1))-1],'pre');%time pad (eg:stat doesnt have baseline)
tempmask=padarray(tempmask,[0,0,siz(3)-find(group_cohSvsD_avg.time==stat.time(end))],'post');%time pad 
tempmask=padarray(tempmask,[0,siz(2)-size(stat.mask,2),0],'post');%freq pad
tempmask=padarray(tempmask,[siz(1)-size(stat.mask,1),0,0],'post');%channel pad
group_cohSvsD_avg.mask=logical(tempmask>0);

% figure;set(gcf,'color','w');ylabel('frequency');xlabel('time');set(gca,'FontSize',15);
labels={'BA7','BA23','BA123','BA6'};%{'occ','temp','frontal'};%{'BA6'};%{'BA23_L','BA23_R','BA6_R','BA47_R'};%{'BA23_L','BA6_L','BA47_L','BA23_R','BA6_R','BA6_RD','BA47_R'};%{'BA23','BA6','BA47'};%{'BA23','BA1.2.3','BA4','BA6_M','BA6_L'};%{'BA7_L','BA7_R','BA39'};%{'OCC','BA6'};%{'STN_L'};%{'BA23','BA1.2.3','BA4','BA6'};%{'occ','temp','frontal'};%{'BA7','BA40','BA4','BA9'};%{'BA8','OCC','BA40','BA6','IFG'};%
for i=1:length(labels)
    cfg.channel=labels{i};
    cfg.xlim = [-0.5 1.6];
    figure;set(gcf,'color','w');ylabel('frequency');xlabel('time');set(gca,'FontSize',15);%subplot(1,2,i);
    ft_singleplotTFR(cfg, group_cohSvsD_avg);
    ylabel ('coherence frequency', 'FontSize',15);xlabel('time','FontSize',15);
    if length(ref)==1
        title([group_cohSvsD_avg.label(i),' <> ',ref],'FontSize',15);
    else
        title([labels{i},'<>right STN'],'FontSize',15);
%         title([labels{i},'<>left STN'],'FontSize',15);
    end
    set(gca,'Fontsize',15);
end

%% trialnumbers
for i=1:13;trl_num_s(i)=size(group_freqS_raw{1,i}.powspctrm,1);end
for i=1:13;trl_num_d(i)=size(group_freqD_raw{1,i}.powspctrm,1);end
mean(trl_num_d)
mean(trl_num_s)
std(trl_num_d)/sqrt(13)
std(trl_num_s)/sqrt(13)
[h p stats]=ttest(trl_num_d,trl_num_s)
%% plot effect over different conditions
n=13;

sigtimeIND_adj=sigtimeIND+75;%stats not run on precue timepoints so need to pad to include up to time zero;

sigtimeIND_adj_c1=sigtimeIND_adj(1:59);%alpha in right dorsal BA6
sigtimeIND_adj_c2=sigtimeIND_adj(59:end);%beta in right dorsal BA6
sigfreqIND_c1=sigfreqIND(1:4);
sigfreqIND_c2=sigfreqIND(5:end);


%%
sigfreqIND_adj=sigfreqIND_c2;
sigtimeIND_adj=sigtimeIND_adj_c2;

for s=1:size(group_cohS_red,2)        
    
        S_spectra(s,:)=nanmean(squeeze(nanmean(group_cohS_red{s}.cohspctrm(sigchIND,:,sigtimeIND_adj),1)),2);
        D_spectra(s,:)=nanmean(squeeze(nanmean(group_cohD_red{s}.cohspctrm(sigchIND,:,sigtimeIND_adj),1)),2);
        
        S_betacoh_timecourse(s,:)=nanmean(squeeze(nanmean(group_cohS_red{s}.cohspctrm(sigchIND,sigfreqIND_adj,:),1)),1);
        D_betacoh_timecourse(s,:)=nanmean(squeeze(nanmean(group_cohD_red{s}.cohspctrm(sigchIND,sigfreqIND_adj,:),1)),1);    
    
    clear temp*
end


%remove NaNs at edge of epoch
[r c]=find(isnan(S_betacoh_timecourse)==1);rm=unique(c);S_betacoh_timecourse(:,rm)=[];D_betacoh_timecourse(:,rm)=[];
group_cohS_red{1}.time(rm)=[];

figure;set(gcf,'color','w');
% subplot(1,2,1);
plot(nanmean(S_spectra),'k','linewidth',5);hold on;plot(nanmean(D_spectra),'b','linewidth',5);ylim([0 0.14])%xticks(1:5:75);
ylabel('Coherence');xlabel('Frequency');set(gca,'Fontsize',15);
% subplot(1,2,2);
% plotpatch(S_betacoh_timecourse,group_cohS_red{1}.time,'k');hold on;plotpatch(D_betacoh_timecourse,group_cohS_red{1}.time,'b');
% ylabel('coherence');xlabel('time');set(gca,'Fontsize',15);
% % suptitle('right STN <> right frontal');set(gca,'Fontsize',15);
% suptitle('left STN <> left BA4');set(gca,'Fontsize',15);

%  with confidence intervals
x = group_cohS_red{1}.time;
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
xlabel ('Time from cue onset (s)', 'FontSize',15);
ylabel ('BETA coherence Diff-Same (+95%CI)', 'FontSize',15);
xlim([-0.5 1.6])
set(gca,'Fontsize',15);
% xline(0,'k');xline(0.8,'k');yline(0,'k');
ylim([-0.05 0.05])

% % % %individual
% figure;set(gcf,'color','w');
% for s=1:13
%     subplot(4,4,s);
%     plot(S_spectra(s,:),'k','linewidth',3);hold on;plot(D_spectra(s,:),'b','linewidth',3);xticks(2:10:75);ylim([0 0.2])
%     ylabel('coherence');xlabel('frequency');
% end
% suptitle('right STN <> frontal');
% 
%% by contact
% R01_coh=yMean;
% R12_coh=yMean;
% R23_coh=yMean;

figure;set(gcf,'color','w');
plot(x, R01_coh,'k','LineWidth',2)
hold on
plot(x, R12_coh,'b','LineWidth',2)
plot(x, R23_coh,'c','LineWidth',2)
plot(x, mean(cat(1,R01_coh,R12_coh,R23_coh)),'r','LineWidth',2)
hold off
xlabel ('Time from cue onset (s)', 'FontSize',15);
ylabel ('BETA coherence Diff-Same (+95%CI)', 'FontSize',15);
xlim([-0.5 1.6])
set(gca,'Fontsize',15);


%% effect size
cfg=[];
cfg.parameter       = 'cohspctrm';
cfg.keepindividual ='yes';
group_cohS_red_avg=ft_freqgrandaverage(cfg,group_cohS_red{:});
group_cohD_red_avg=ft_freqgrandaverage(cfg,group_cohD_red{:});

nsub=13;
siz=size(group_cohD_red_avg.cohspctrm);
x1=squeeze(mean(mean(mean(group_cohS_red_avg.cohspctrm(:,sigchIND,sigfreqIND,sigtimeIND_adj),2),3),4));
x2=squeeze(mean(mean(mean(group_cohD_red_avg.cohspctrm(:,sigchIND,sigfreqIND,sigtimeIND_adj),2),3),4));
     
%or for within subj/paired samples
% cohensd = mean(x1-x2) ./ std(x1-x2);
cohensd = mean(abs(x1-x2)) ./ std(abs(x1-x2));
disp(cohensd)
%% check out plain TF plots

STN=0;%plot for STN or Cortical sites
nch=7;

for s=1:13
    
   SvsD_group{s}= group_freqS{s};
   
   if STN==1
       SvsD_group{s}.powspctrm = group_freqD{s}.powspctrm(7:end,:,:)-group_freqS{s}.powspctrm(7:end,:,:);
       SvsD_group{s}.label={'STN R01';'STN R12';'STN R23';'STN L01';'STN L12';'STN L23'};
   else
       SvsD_group{s}.powspctrm = group_freqD{s}.powspctrm(1:nch,:,:)-group_freqS{s}.powspctrm(1:nch,:,:);
%        temp = group_freqD{s}.powspctrm(1:nch,:,:)-group_freqS{s}.powspctrm(1:nch,:,:);
%        SvsD_group{s}.powspctrm(1,:,:) = nanmean(temp([1:2],:,:),1);
%        SvsD_group{s}.powspctrm(2,:,:) = nanmean(temp([3:4],:,:),1);
%        SvsD_group{s}.powspctrm(3,:,:) = nanmean(temp([5:6],:,:),1);%nan(1,size(SvsD_group{s}.powspctrm,2),size(SvsD_group{s}.powspctrm,3));%
%        SvsD_group{s}.powspctrm(4:end,:,:) =[];
       SvsD_group{s}.label={'BA23_L','BA6_L','BA47_L','BA23_R','BA6_R','BA6_RD','BA47_R'};%{'BA7','BA23','BA40','BA6'};%{'BA23','BA6','BA47'};%{'BA7','BA23','BA40','BA6'};%{'OCC';'TEMP';'BA6'};%{'lateral BA6';'dorsal BA6'};%
   end
   
end

% %add in NaNs for missing channels in 2 subjects
% SvsD_group{8}.powspctrm(6,:,:)=NaN;SvsD_group{8}.powspctrm=circshift(SvsD_group{8}.powspctrm,1);
% if STN==1;SvsD_group{8}.label={'STN_R01';'STN_R12';'STN_R23';'STN_L01';'STN_L12';'STN_L23'};else SvsD_group{8}.label={'O1';'O2';'T1';'T2';'F1';'F2'};end
% SvsD_group{11}.powspctrm(4:6,:,:)=NaN;
% if STN==1;SvsD_group{11}.label={'STN_R01';'STN_R12';'STN_R23';'STN_L01';'STN_L12';'STN_L23'};else SvsD_group{11}.label={'O1';'O2';'T1';'T2';'F1';'F2'};end
% 

% cfg=[];
% cfg.baseline = [0 0.8];
% cfg.baselinetype ='db';
% for i=1:13;SvsD_group{i} = ft_freqbaseline(cfg, SvsD_group{i});end

cfg=[];
cfg.keepindividual='no';
% SvsD_group_avg=ft_freqgrandaverage(cfg,SvsD_group{[1:7,9,10,12,13]});
SvsD_group_avg=ft_freqgrandaverage(cfg,SvsD_group{:});

%ALL CHANS
figure;set(gcf,'color','w');
chans=[1,4,5,7];%ids of the 4 FWE patient voxels from coonjunction
for ch=1:size(chans,2);%(SvsD_group_avg.label,1)
subplot(2,2,ch)    
cfg=[];
cfg.channel=SvsD_group_avg.label(chans(ch));
cfg.zlim=[-0.1 0.1];
cfg.xlim=[-0.5 1.6];
ft_singleplotTFR(cfg,SvsD_group_avg);
set(gca,'Fontsize',15);
end
suptitle('different vs same');
ylabel('frequency', 'FontSize',15);
xlabel('time from cue onset', 'FontSize',15);
%%
for s=1:13
       temp = group_freqD{s}.powspctrm(1:nch,:,:);
%        group_freqD{s}.powspctrm(1,:,:) = nanmean(temp([1:2],:,:),1);
%        group_freqD{s}.powspctrm(2,:,:) = nanmean(temp([3:4],:,:),1);
%        group_freqD{s}.powspctrm(3,:,:) = nanmean(temp([5:6],:,:),1);%nan(1,size(SvsD_group{s}.powspctrm,2),size(SvsD_group{s}.powspctrm,3));%
       group_freqD{s}.powspctrm(nch+1:end,:,:) =[];
       group_freqD{s}.label={'BA23_L','BA6_L','BA47_L','BA23_R','BA6_R','BA6_RD','BA47_R'};%{'BA7','BA23','BA40','BA6'};%{'BA23','BA6','BA47'};%{'OCC';'TEMP';'BA6'};%{'BA7','BA23','BA40','BA6'};%{'lateral BA6';'dorsal BA6'};%
end

for s=1:13
       temp = group_freqS{s}.powspctrm(1:nch,:,:);
%        group_freqS{s}.powspctrm(1,:,:) = nanmean(temp([1:2],:,:),1);
%        group_freqS{s}.powspctrm(2,:,:) = nanmean(temp([3:4],:,:),1);
%        group_freqS{s}.powspctrm(3,:,:) = nanmean(temp([5:6],:,:),1);%nan(1,size(SvsD_group{s}.powspctrm,2),size(SvsD_group{s}.powspctrm,3));%
       group_freqS{s}.powspctrm(nch+1:end,:,:) =[];
       group_freqS{s}.label={'BA23_L','BA6_L','BA47_L','BA23_R','BA6_R','BA6_RD','BA47_R'};%{'BA7','BA23','BA40','BA6'};%{'BA23','BA6','BA47'};%{'OCC';'TEMP';'BA6'};%{'BA7','BA23','BA40','BA6'};%{'lateral BA6';'dorsal BA6'};%
end

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
Nsub=13;
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

[stat] = ft_freqstatistics(cfg, group_freqS{:}, group_freqD{:});

for i=1:length(group_freqD{1}.label);figure;imagesc(squeeze(stat.prob(i,:,:)));end


%% between BA6 dorsal(SMA) and lateral(premotor) & BA23 (PCC)
%lateral BA6 52,-7,44; dorsal BA6 7 2 69
load('COH_TFmeth2long_cuei+1_notsss_cluster1&2_conj_patients.mat')

group_cohD_red=group_cohD;group_cohS_red=group_cohS;
for s=1:13
    if s==8
        group_cohD_red{s}.cohspctrm=group_cohD_red{s}.cohspctrm(39,:,:);
        group_cohS_red{s}.cohspctrm=group_cohS_red{s}.cohspctrm(39,:,:);
    elseif s==11
        group_cohD_red{s}.cohspctrm=group_cohD_red{s}.cohspctrm(31,:,:);
        group_cohS_red{s}.cohspctrm=group_cohS_red{s}.cohspctrm(31,:,:);
    else
        group_cohD_red{s}.cohspctrm=group_cohD_red{s}.cohspctrm(43,:,:);
        group_cohS_red{s}.cohspctrm=group_cohS_red{s}.cohspctrm(43,:,:);
    end
    group_cohS_red{s}.label={'BA6_LvsD'};
    group_cohD_red{s}.label={'BA6_LvsD'};
    group_cohS_red{s}.dimord='chan_freq_time';    
    group_cohD_red{s}.dimord='chan_freq_time'; 
end

group_cohD_red=group_cohD;group_cohS_red=group_cohS;
for s=1:13
    group_cohD_red{s}.cohspctrm=group_cohD_red{s}.cohspctrm(5,:,:);
    group_cohS_red{s}.cohspctrm=group_cohS_red{s}.cohspctrm(5,:,:);
    group_cohS_red{s}.label={'BA6RD_BA23L'};
    group_cohD_red{s}.label={'BA6RD_BA23L'};
end

group_cohD_red=group_cohD;group_cohS_red=group_cohS;
for s=1:13
    group_cohD_red{s}.cohspctrm=group_cohD_red{s}.cohspctrm(4,:,:);
    group_cohS_red{s}.cohspctrm=group_cohS_red{s}.cohspctrm(4,:,:);
    group_cohS_red{s}.label={'BA6RL_BA23L'};
    group_cohD_red{s}.label={'BA6RL_BA23L'};
    group_cohS_red{s}.dimord='chan_freq_time';    
    group_cohD_red{s}.dimord='chan_freq_time'; 
end

%% itpc

for s=1:13
    
    SvsD_group_itpc{s}= group_freqD_itpc{s};
    SvsD_group_itpc{s}.itpc = group_freqD_itpc{s}.itpc(1:nch,:,:)-group_freqS_itpc{s}.itpc(1:nch,:,:);
    SvsD_group_itpc{s}.label={'BA23_L','BA6_L','BA47_L','BA23_R','BA6_R','BA6_RD','BA47_R'};%
    
    cfg=[];
    cfg.keepindividual='no';
    cfg.parameter='itpc';
    SvsD_group_itpc_avg=ft_freqgrandaverage(cfg,SvsD_group_itpc{:});

end

%ALL CHANS
figure;set(gcf,'color','w');
for ch=1:size(SvsD_group_itpc_avg.label,1)
subplot(3,3,ch)    
cfg=[];
cfg.channel=SvsD_group_itpc_avg.label(ch);
cfg.zlim=[-0.1 0.1];
cfg.xlim=[-0.5 1.6];
cfg.parameter='itpc';
ft_singleplotTFR(cfg,SvsD_group_itpc_avg);
set(gca,'Fontsize',15);
end
suptitle('ITPC different vs same');
ylabel('frequency', 'FontSize',15);
xlabel('time from cue onset', 'FontSize',15);


for s=1:13
       group_freqD_itpc{s}.itpc(nch+1:end,:,:)=[];
       group_freqD_itpc{s}.label={'BA23_L','BA6_L','BA47_L','BA23_R','BA6_R','BA6_RD','BA47_R'};%{'BA7','BA23','BA40','BA6'};%{'BA23','BA6','BA47'};%{'OCC';'TEMP';'BA6'};%{'BA7','BA23','BA40','BA6'};%{'lateral BA6';'dorsal BA6'};%
end

for s=1:13
       group_freqS_itpc{s}.itpc(nch+1:end,:,:)=[];
       group_freqS_itpc{s}.label={'BA23_L','BA6_L','BA47_L','BA23_R','BA6_R','BA6_RD','BA47_R'};%{'BA7','BA23','BA40','BA6'};%{'BA23','BA6','BA47'};%{'OCC';'TEMP';'BA6'};%{'BA7','BA23','BA40','BA6'};%{'lateral BA6';'dorsal BA6'};%
end

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
cfg.latency = [0 1.6];
cfg.avgovertime = 'no';
cfg.neighbours    =[];
cfg.minnbchan =0;
cfg.spmversion = 'spm12' ;
Nsub=13;
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

cfg.parameter='itpc';
[stat_itpc] = ft_freqstatistics(cfg, group_freqS_itpc{:}, group_freqD_itpc{:});

for i=1:length(group_freqS_itpc{1}.label);figure;imagesc(squeeze(stat_itpc.prob(i,:,:)));end



%% granger
   
subjects = {'LN_C31','LN_C35','LN_C37','LN_C38','LN_C39','LN_C48','LN_C56','LN_C57','LN_C58','LN_C59','LN_C60','LN_C61','LN_C62'};...

for s = 1:13

     eval(['load(''',subjects{s},'_cluster1_2_CONJ_P_FWE_raw_long_timecourse_BF.mat'')']);%epoched, raw, longer than what we need: -1 to 1.6s, for TF windowing and up to cue i+1

    %%%%sort conditions based on cue sequence
    same = [];same_1=[];same_0=[];same_pen=[];
    diff = [];diff_1 = [];diff_0 = [];diff_pen=[];
    bad = [];    % already removed
    for t = 1:length(seq)
        
        %by accuracy
        temp = seq{t}.resp;if strcmp(temp,'Left_correct')==1 || strcmp(temp,'Right_correct')==1 ; acc(t)=1;else acc(t)=0;end
        
        index = seq{t}.ind; %index of cue within a trial
        if index > 1 && index < length (seq{t}.seq)  %ignore the first and the last trial
%         if index > 1 && index == length (seq{t}.seq) -1   % penultimate
%         if index > 1 && index == 2   % random cue number

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
    cfg=[];
    cfg.trials =same;
    dataS=ft_selectdata(cfg,data);
    cfg.trials =[];cfg.trials =diff;
    dataD=ft_selectdata(cfg,data);
    
    
    cfg         = [];
    cfg.order   = 5;
    cfg.toolbox = 'bsmart';    
    mdataS       = ft_mvaranalysis(cfg, dataS);
    mdataD       = ft_mvaranalysis(cfg, dataD);

    cfg        = [];
    cfg.method = 'mvar';
    mfreqS      = ft_freqanalysis(cfg, mdataS);
    mfreqD      = ft_freqanalysis(cfg, mdataD);
    
    cfg           = [];
    cfg.method    = 'granger';
    grangerS       = ft_connectivityanalysis(cfg, mfreqS);
    grangerD       = ft_connectivityanalysis(cfg, mfreqD);
    
    group_cohS_granger{s}=grangerS;
    group_cohD_granger{s}=grangerD; 
        
    
end   
   

cfg=[];
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.method ='montecarlo';
cfg.correctm = 'cluster';
cfg.clusteralpha  = 0.05;
cfg.alpha       = 0.05;
cfg.tail        = 0;%0 default
cfg.parameter = 'grangerspctrm';
cfg.numrandomization = 1000;
cfg.frequency = [1 30];
cfg.neighbours       = [];
cfg.avgoverfreq='no';
cfg.spmversion = 'spm12' ;

Nsub=length(group_cohS_granger);

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

[stat] = ft_freqstatistics(cfg, group_cohS_granger{:}, group_cohD_granger{:});%zscore diff vs dummy

for i=1:size(stat.mask,1);figure;imagesc(squeeze(stat.prob(i,:,:)));end

group_cohSvsD=group_cohS_granger;

for s=1:size(group_cohSvsD,2)
    
    group_cohSvsD{s}.grangerspctrm = group_cohD_granger{s}.grangerspctrm-group_cohS_granger{s}.grangerspctrm;%%DvsS
    
end

for s=1:13
    figure;
    cfg           = [];
    cfg.parameter = 'grangerspctrm';
%     cfg.zlim      = [0 1];
    ft_connectivityplot(cfg, group_cohSvsD{s});
    
end 
    
%%FT: The presence of transient activity such as evoked responses is a major violation of the assumptions in Granger analysis. Therefore we will focus on time intervals before stimulus presentation that is not contaminated by evoked responses.
    
    
    
 
        
