% consolidate behavioural data


clear
cd V:\

P=0;


if P==1
wanted={'LN_C31','LN_C35','LN_C37','LN_C38','LN_C39','LN_C43','LN_C46','LN_C48','LN_C56','LN_C57','LN_C58','LN_C59','LN_C60','LN_C61','LN_C62'};
else
wanted={'MG4492','MG4493','MG4521','MG4551','MG4552','MG4554','MG4558','MG4562','MG4563','MG06109','MG06114','MG06117','MG06126'};
end

for s=1:size(wanted,2)

    cd(fullfile(wanted{s}))
    
    if P==1 && s==13 || s==15
    listing=dir('*\Behavioural\2*time*.mat');%LN_C62 different spelling and only 1 file
    else
    listing=dir('*\Behavioral\2*time*.mat');
    end
    
    sessfix=1;
    for sess=1:size(listing,1)
     
    cd(listing(sess).folder)
    load(listing(sess).name);
    
    if isnan(choice)==1 | isnan(correct)==1 | size(correct,2) < 10
        continue
    else
    behav{s}.session(sessfix).acc=correct;
    behav{s}.session(sessfix).choice=choice;
    behav{s}.session(sessfix).RT=RT;
    behav{s}.session(sessfix).stimuli=stimuli;
    sessfix=sessfix+1;
    end
    
    end
            
    cd V:\

end

clearvars -except behav wanted

% some fixes

for s=1:size(behav,2)
    
    for i=1:size(behav{s}.session,2)
    
            for n=1:size(behav{s}.session(i).stimuli,2)
            behav{s}.session(i).stimnum(n)=size(behav{s}.session(i).stimuli{n},2);
            end
  
    end
end

% behav{15}.session(2)=[];% repeated data LN_C62

%%
cd D:\
save behavioural_data_controls behav wanted
%% analysis
clear
load('behavioural_data.mat');new_RT_cutoff=278;%207;%PATIENTS
% load('behavioural_data_controls.mat');new_RT_cutoff=305%234;%CONTROLS
% k=1;
P=1;
fig=0;
figRT=0;
figcorr=0;
sumfig=0;
summary_dump1=[];summary_dump2=[];summary_dump3=[];summary_dump4=[];summary_dump5=[];summary_dump6=[];summary_dump7=[];summary_dump8=[];summary_dump9=[];summary_dump10=[];summary_dumpRT_diff=[];summary_dumpRT_same=[];
summary_dump_bayes=[];

for s=[1:5,8:size(behav,2)]
    
    for i=1:size(behav{s}.session,2)

        for n=1:size(behav{s}.session(i).stimuli,2)
        
            A(n)=sum(behav{s}.session(i).stimuli{n}==1);
            B(n)=sum(behav{s}.session(i).stimuli{n}==2);
            
            C(n)=behav{s}.session(i).choice(n);
            
            %convert to choice vs alternative (from L/R)
            if C(n)==1
                CHOICE(n)=A(n)-B(n);
                ALT(n)=B(n);
                ALLSTIM(n)=A(n)+B(n);
            else
                CHOICE(n)=B(n)-A(n);
                ALT(n)=A(n);
                ALLSTIM(n)=A(n)+B(n);
            end
          
            
            if behav{s}.session(i).choice(n)==behav{s}.session(i).stimuli{n}(end)
                prob_choice_laststim(n)=1;
            else
                prob_choice_laststim(n)=0;
            end
            
            if behav{s}.session(i).choice(n)==behav{s}.session(i).stimuli{n}(1)
                choice_match1stim(n)=1;
            else
                choice_match1stim(n)=0;
            end
            if size(behav{s}.session(i).stimuli{n},2)>1
                if behav{s}.session(i).choice(n)==behav{s}.session(i).stimuli{n}(2)
                    choice_match2stim(n)=1;
                else
                    choice_match2stim(n)=0;
                end
            else
                choice_match2stim(n)=nan;
            end
            if size(behav{s}.session(i).stimuli{n},2)>2
                
                if behav{s}.session(i).choice(n)==behav{s}.session(i).stimuli{n}(3)
                    choice_match3stim(n)=1;
                else
                    choice_match3stim(n)=0;
                end
            else
                choice_match3stim(n)=nan;
            end
            if size(behav{s}.session(i).stimuli{n},2)>3
                
                if behav{s}.session(i).choice(n)==behav{s}.session(i).stimuli{n}(4)
                    choice_match4stim(n)=1;
                else
                    choice_match4stim(n)=0;
                end
            else
                choice_match4stim(n)=nan;
            end
            if size(behav{s}.session(i).stimuli{n},2)>4
                
                if behav{s}.session(i).choice(n)==behav{s}.session(i).stimuli{n}(5)
                    choice_match5stim(n)=1;
                else
                    choice_match5stim(n)=0;
                end
            else
                choice_match5stim(n)=nan;
            end
            
            
        end
   
        [bA, mA, nA] = unique(CHOICE);
%         [bB, mB, nB] = unique(ALT);
        [bB, mB, nB] = unique(ALLSTIM);
        counter = accumarray([nA(:), nB(:)],1);
        
        %make a matrix with occurances for scatter size
        counterT=[counter;bB];
        counterT=[counterT cat(2,bA,999)'];
      
        if fig
        figure;set(gcf,'color','w');
        for n=1:size(behav{s}.session(i).stimuli,2)
            
            if behav{s}.session(i).acc(n)==1; c=[0 1 0];else 
                c=[0 0 0];
            end
            
            temp=find(counterT(:,end)==CHOICE(n));
%             temp2=find(counterT(end,:)==ALT(n));
            temp2=find(counterT(end,:)==ALLSTIM(n));
            
%             subplot(2,14,k)
%             scatter(CHOICE(n),ALT(n),counter(temp,temp2)*10,c,'filled')
            scatter(ALLSTIM(n),CHOICE(n),counter(temp,temp2)*10,c,'filled')
            hold on
            
        end
        errorbar(mean(ALLSTIM),mean(CHOICE),std(CHOICE),std(CHOICE),std(ALLSTIM),std(ALLSTIM),'or','CapSize',0);
        
%         yl=ylim;xl=xlim;if  yl(2)>xl(2) ; x=yl(2); else x=xl(2); end
%         xlabel('Evidence for Choice', 'FontSize',15);ylabel('Alternative', 'FontSize',15);set(gca,'FontSize',15);xlim([1,x]);ylim([1,x]);
        xlabel('# Stimuli seen', 'FontSize',15);ylabel('Evidence for Chosen', 'FontSize',15);set(gca,'FontSize',15);%xlim([1,x]);ylim([1,x]);
        title([wanted{s},' session ',num2str(i)]);
        end
        

        if figRT
            figure;set(gcf,'color','w');
            if P==1 && s==1 && i==2
            histogram(behav{s}.session(i).RT - (behav{s}.session(i).stimnum -1)*1000,10); title('RT norm by stimnum', 'FontSize',15);
            else
            histogram(behav{s}.session(i).RT - (behav{s}.session(i).stimnum -1)*800,10); title('RT norm by stimnum', 'FontSize',15);
            end
        end
%         %get RTs based on elapsed time and decide if 'last' stimulus really
%         %last or t-1
        if P==1 && s==1 && i==2
            behav{s}.session(i).newRTs = behav{s}.session(i).RT - (behav{s}.session(i).stimnum -1)*1000;
        else
            behav{s}.session(i).newRTs = behav{s}.session(i).RT - (behav{s}.session(i).stimnum -1)*800;
        end
        
        %check if sometimes after 'last' stim not n-1
        RT_long=find(behav{s}.session(i).newRTs>800);
        if sum(RT_long)>0
            if P==1 && s==1 && i==2
            behav{s}.session(i).newRTs(RT_long)=behav{s}.session(i).newRTs(RT_long)-1000;
            else
            behav{s}.session(i).newRTs(RT_long)=behav{s}.session(i).newRTs(RT_long)-800;
            end
        end
        
        behav{s}.session(i).last(behav{s}.session(i).newRTs <= 0) = -1;
        behav{s}.session(i).last(behav{s}.session(i).newRTs > 0) = 1;
     
        
        
        
%       %match LFP analysis - same vs different
        for n=1:size(behav{s}.session(i).stimuli,2)
                      
            stimuli= behav{s}.session(i).stimuli{n};if length(stimuli)==1;continue;else;end
            choice= behav{s}.session(i).choice(n);
            
            [evidence_chosen, evidence_left, evidence_current, evidence_absolute, current_left, same, agree, evidence_total] = trial_predictors_stim_by_choice (stimuli, choice);
     
            end_same(n)=stimuli(end)==stimuli(end-1);
            end_evidence_for_chosen(n)=evidence_chosen(end);
            end_sum_same(n)=sum(same==1);%#of same pairs in sequence
            end_same_ratio(n)=sum(same==1)/behav{s}.session(i).stimnum(n);%#of same pairs in sequence compared to whole
            
            acc(n)=behav{s}.session(i).acc(n);
            RT(n)=behav{s}.session(i).newRTs(n);
            stimnum(n)=behav{s}.session(i).stimnum(n);
            
            bayes_given_stim=mean(stimuli);
            if bayes_given_stim == 1.5
                bayes_given_stim = datasample([1,2],1,'Replace',true);
            else
                bayes_given_stim =round(mean(stimuli));
            end
            if bayes_given_stim == choice
                acc_bayes(n)=acc(n);
            elseif bayes_given_stim ~= choice && acc(n)==0
                acc_bayes(n)=1;
            elseif bayes_given_stim ~= choice && acc(n)==1
                acc_bayes(n)=0;
            end
        end
        
        if figcorr
        figure;set(gcf,'color','w');scatter(end_sum_same,RT,'k');hold on;scatter(end_evidence_for_chosen,RT,'r');
        figure;set(gcf,'color','w');scatter(end_same_ratio,RT,'g');
        end
        
        summary(s,i,1) = nanmean(behav{s}.session(i).stimnum);
        summary(s,i,2) = nanstd(behav{s}.session(i).stimnum);
        summary(s,i,3) = nanmean(behav{s}.session(i).acc);     
        
        
        %%PUT INTO OUTPUT STRUCTURE
        %remove very fast RTs
        RT_excl=find(RT<=new_RT_cutoff);RT(RT_excl)=[];end_same_RT=end_same;end_same_RT(RT_excl)=[];
        summary(s,i,17)=size(RT_excl,2);
        
        size(RT_excl)
        
        summary(s,i,7)=nanmean(RT(end_same_RT==0));summary(s,i,8)=nanmean(RT(end_same_RT==1));%RT x whether end is on 'same'
        summary_dumpRT_diff = [summary_dumpRT_diff, RT(end_same_RT==0) ]; 
        summary_dumpRT_same = [summary_dumpRT_same, RT(end_same_RT==1) ]; 
       
        summary(s,i,4) = nanmean(RT);     
        summary(s,i,5) = nanstd(RT); 
        
        summary_dump1 = [summary_dump1, RT]; 
        stimnum2=stimnum;stimnum2(RT_excl)=[];
        summary_dump2 = [summary_dump2, stimnum2]; 
        summary_dump3 = [summary_dump3, end_same]; 
        summary_dump4 = [summary_dump4, prob_choice_laststim];
        summary_dump5 = [summary_dump5, choice_match1stim];
        summary_dump6 = [summary_dump6, choice_match2stim];
        summary_dump7 = [summary_dump7, choice_match3stim];
        summary_dump8 = [summary_dump8, choice_match4stim];
        summary_dump9 = [summary_dump9, choice_match5stim];
        summary_dump10 = [summary_dump10, ones(1,size(choice_match5stim,2))*s];
        summary_dump_bayes = [summary_dump_bayes, acc_bayes];
        
        summary(s,i,6) = ttest(behav{s}.session(i).acc,0.5);  
       
     
        
        summary(s,i,9)=nanmean(acc(end_same==0));summary(s,i,10)=nanmean(acc(end_same==1));%ACC x whether end is on 'same'
        summary(s,i,11)=nanmean(stimnum(end_same==0));summary(s,i,12)=nanmean(stimnum(end_same==1));%stimnum x whether end is on 'same'
        summary(s,i,13)=nanmean(end_same);
        
        summary(s,i,14)=nanmean(end_sum_same);
        summary(s,i,15)=nanmean(end_same_ratio);
        summary(s,i,16)=nanmean(end_evidence_for_chosen);
        
        summary(s,i,18)=size(behav{s}.session(i).stimuli,2);
        
        summary(s,i,19)=nanmean(acc_bayes);
        
        summary_RT_trlwise{s,i,1}=RT(end_same_RT==1);
        summary_RT_trlwise{s,i,2}=RT(end_same_RT==0);
        summary_RT_trlwise{s,i,3}=RT;
        summary_RT_trlwise{s,i,4}=RT_excl;
        
        stinum_excl=stimnum;stinum_excl(RT_excl)=[];
        summary_stimnum_trlwise{s,i,1}=stinum_excl;

        
        clearvars -except summary* s i behav *fig* wanted new_RT_cutoff P
        
    end
    
end
%remove without 2ndsession
summary(summary==0)=nan;

%get trialnums
summarytrlnum=nansum(cat(2,summary(:,1,18),summary(:,2,18)),2);

%avg over sessions
summary=squeeze(nanmean(summary(:,:,1:19),2));

%remove s6/7 as chance performance (summary column 6 is ttest of acc)
if P==1
summary([6,7],:)=[];
summarytrlnum([6,7])=[];
end

bins=30;figure;set(gcf,'color','w');hist(summary_dump1,bins);
[X,Y,B]=histcounts(summary_dump1,bins);[X2,Y2,B2]=histcounts(summary_dump2,bins);
for i=1:bins
    probchoice_RTbins(i)=mean(summary_dump4(B==i));
    probchoice_RTbins_sem(i)=std(summary_dump4(B==i))/sqrt(sum(B==i));
    probchoice_Sbins(i)=mean(summary_dump4(B2==i));
    probchoice_Sbins_sem(i)=std(summary_dump4(B2==i))/sqrt(sum(B2==i));
end
figure;set(gcf,'color','w');barwitherr(probchoice_RTbins_sem,probchoice_RTbins);
ylabel('probability choice == last stim', 'FontSize',15);xlabel('RT bins', 'FontSize',15);set(gca,'FontSize',15);
short_RTs=summary_dump1(B==1|B==2|B==3|B==4);new_RT_cutoff=max(short_RTs);

figure;set(gcf,'color','w');barwitherr(probchoice_Sbins_sem,probchoice_Sbins);
ylabel('probability choice == last stim', 'FontSize',15);xlabel('Stim# bins', 'FontSize',15);set(gca,'FontSize',15);

% figure;set(gcf,'color','w');histogram(summary_dump2,max(summary_dump2));summary_dump2_bin=mod(summary_dump2,2)==0;

if sumfig
%accuracy by session and average number of stimuli before choice
figure;set(gcf,'color','w');
BAR=bar(summary(:,3),'k');%BAR(1).FaceColor=[0 0 0];BAR(2).FaceColor=[0.5 0.5 0.5];
ylabel('accuracy', 'FontSize',15);xlabel('subjects', 'FontSize',15);set(gca,'FontSize',15);
hold on;plot(xlim,[0.5 0.5],'r');

figure;set(gcf,'color','w');
B2=barwitherr(summary(:,2),summary(:,1),'k');%B2(1).FaceColor=[0 0 0];B2(2).FaceColor=[0.5 0.5 0.5];
ylabel('number of stimuli before choice', 'FontSize',15);xlabel('subjects', 'FontSize',15);set(gca,'FontSize',15);

figure;set(gcf,'color','w');
B3=barwitherr(summary(:,5),summary(:,4),'k');%B3(1).FaceColor=[0 0 0];B3(2).FaceColor=[0.5 0.5 0.5];
ylabel('RT (ms)', 'FontSize',15);xlabel('subjects', 'FontSize',15);set(gca,'FontSize',15);


figure;set(gcf,'color','w');
B3=barwitherr(cat(2,std(nanmean(summary(:,8),2))/sqrt(s),std(nanmean(summary(:,7),2))/sqrt(s)),cat(2,mean(nanmean(summary(:,8),2)),mean(nanmean(summary(:,7),2))),'k');
ylabel('RT (ms) end same or different', 'FontSize',15);xlabel('last pair of stimuli', 'FontSize',15);xticklabels({'same','different'}); set(gca,'FontSize',15);

figure;set(gcf,'color','w');
B3=barwitherr(cat(2,std(nanmean(summary(:,10),2))/sqrt(s),std(nanmean(summary(:,9),2))/sqrt(s)),cat(2,mean(nanmean(summary(:,10),2)),mean(nanmean(summary(:,9),2))),'k');
% B3(1).FaceColor=[0 0 0];B3(2).FaceColor=[0.5 0.5 0.5];
ylabel('ACC end same or different', 'FontSize',15);xlabel('last pair of stimuli', 'FontSize',15);xticklabels({'same','different'}); set(gca,'FontSize',15);


figure;set(gcf,'color','w');
scatter(summary(:,14),summary(:,4),'k');% B3(1).FaceColor=[0 0 0];B3(2).FaceColor=[0.5 0.5 0.5];
xlabel('sum of ''same'' pairs in sequence', 'FontSize',15);ylabel('RT', 'FontSize',15);set(gca,'FontSize',15);set(gca,'FontSize',15);

figure;set(gcf,'color','w');
scatter(summary(:,15),summary(:,4),'k');
% B3(1).FaceColor=[0 0 0];B3(2).FaceColor=[0.5 0.5 0.5];
xlabel('ratio of ''same'' in whole sequence', 'FontSize',15);ylabel('RT', 'FontSize',15);set(gca,'FontSize',15);set(gca,'FontSize',15);

figure;set(gcf,'color','w');
scatter(summary(:,16),summary(:,4),'k');
% B3(1).FaceColor=[0 0 0];B3(2).FaceColor=[0.5 0.5 0.5];
xlabel('evidence for chosen at end', 'FontSize',15);ylabel('RT', 'FontSize',15);set(gca,'FontSize',15);

figure;set(gcf,'color','w');
scatter(summary(:,13),summary(:,4),'k');
% B3(1).FaceColor=[0 0 0];B3(2).FaceColor=[0.5 0.5 0.5];
xlabel('ratio of responses after ''same'' pair at end', 'FontSize',15);ylabel('RT', 'FontSize',15);set(gca,'FontSize',15);set(gca,'FontSize',15);

figure;set(gcf,'color','w');
B2=bar(summary(:,13),'k');%B2(1).FaceColor=[0 0 0];B2(2).FaceColor=[0.5 0.5 0.5];
ylabel('ratio of responses after ''same'' pair at end', 'FontSize',15);xlabel('subjects', 'FontSize',15);set(gca,'FontSize',15);
hold on;plot(xlim,[0.58 0.58],'r');
% figure;set(gcf,'color','w');
% B4=bar(summary(:,:,5));B4(1).FaceColor=[0 0 0];B4(2).FaceColor=[0.5 0.5 0.5];
% ylabel('IES score', 'FontSize',15);xlabel('subjects', 'FontSize',15);set(gca,'FontSize',15);
end

%% samples hist
load('summary_patients.mat')

figure;set(gcf,'color','w');
H=histogram(summary_dump2,25,'FaceAlpha',1,'FaceColor','r','FaceAlpha',0.4);
% H=histogram(summary_dump2/sum(summarytrlnum),25,'FaceAlpha',1,'FaceColor','r','FaceAlpha',0.6);
xlabel('number of stimuli before choice', 'FontSize',15);
ylabel('occurences', 'FontSize',15)

hold on

samples=size(summary_dump2,2);

clearvars -except H samples

load('summary_controls.mat')

% histogram(summary_dump2/sum(summarytrlnum),15,'FaceAlpha',1,'FaceAlpha',0.4,'FaceColor','b')%normaize by trials completed
histogram(datasample(summary_dump2,samples,'Replace',false),H.BinEdges,'FaceAlpha',0.4,'FaceColor','b')%samples same amount of cues randomly
% histogram(summary_dump2,H.BinEdges,'FaceAlpha',0.4,'FaceColor','b')
set(gca,'FontSize',15);


figure;set(gcf,'color','w');
p1=histogram(summary_dump1(summary_dump2<6));
hold on;
p2=histogram(summary_dump1(summary_dump2>6),'FaceAlpha',1,'FaceColor','r','FaceAlpha',0.3);
xlabel('RT', 'FontSize',15);
ylabel('occurences', 'FontSize',15)
lgd=legend ([p1(1) p2(1) ],'Box', 'off');
lgd.String={'less than 6 cues', 'more than 6 cues'};


%% correlations (stim#, acc, RT, forgetting)

[h p]=ttest(summary(:,13),0.5)

for i=3:4%acc&RT
    for ii=13:16%various same ratios
        [r p]=corr(cat(2,summary(:,ii),summary(:,i)));
        R(i,ii)=r(2,1);P(i,ii)=p(2,1);
    end
end

load('R_dc_choice_patients.mat')
[r p]=corrcoef(cat(2,mean(summary(:,:,1),2),mean(summary(:,:,3),2),mean(summary(:,:,4),2),mean(1-abs(LAMBDA),2)))


load('D:\Data\DBS-MEG\BETA_samediff_timecourses.mat')
BETA_samevs_diff=avsame-avdiff;BETA_samevs_diff=nanmean(nanmean(BETA_samevs_diff(:,:,13:20),3),2);%sig timepoint with ttest(only)
RT_samevsdiff=summary(:,8)-summary(:,7);


for s=1:9%[1:5,8:15]
tempS=cat(2,summary_stimnum_trlwise{s,1,1},summary_stimnum_trlwise{s,2,1});%stimnum\
tempRT=cat(2,summary_RT_trlwise{s,1,3},summary_RT_trlwise{s,2,3});%RTs 
temp=find(isnan(tempRT));tempRT(temp)=[];tempS(temp)=[];
[r p]=corrcoef(tempS,tempRT);
Rr(s)=r(2);Pp(s)=p(2);
end


%% LFP correl

[ avsame  avsame2 avsame3 avdiff avdiff2 avdiff3,sigA]=plot_same_same(1,0,2)%need to be in Github dir

RT_samevsdiff=summary(:,8)-summary(:,7);
%cue i
BETA_samevs_diff=avdiff-avsame;BETA_samevs_diff=nanmean(squeeze(nanmean(BETA_samevs_diff(:,:,13:20),2)),2);%sig timepoint with ttest(only)

[r p]=corrcoef(BETA_samevs_diff,RT_samevsdiff)

[r p]=corrcoef(BETA_samevs_diff,summary(:,4))%rt
[r p]=corrcoef(BETA_samevs_diff,summary(:,1))%sample
[r p]=corrcoef(BETA_samevs_diff,summary(:,3))%acc
[r p]=corrcoef(BETA_samevs_diff,summary(:,13))%end same

%cue i+1
BETA_samevs_diff2=avdiff2-avsame2;BETA_samevs_diff2=nanmean(nanmean(BETA_samevs_diff2(:,:,[10:16]),2),3);%sig timepoint with ttest(only)
% BETA_samevs_diff2=avdiff2-avsame2;BETA_samevs_diff2=nanmean(nanmean(BETA_samevs_diff2(:,:,[11:17,21:26]),2),3);%sig timepoint with ttest(only)
% BETA_samevs_diff2=avdiff2-avsame2;BETA_samevs_diff2=nanmean(nanmean(BETA_samevs_diff2(:,:,[20:26]),2),3);%sig timepoint with ttest(only)
% BETA_samevs_diff2=abs(avdiff2-avsame2);BETA_samevs_diff2=nanmean(nanmean(BETA_samevs_diff2(:,:,[11:17,21:26]),2),3);%sig timepoint with ttest(only)
% BETA_samevs_diff2=abs(avdiff2-avsame2);BETA_samevs_diff2=nanmean(nanmean(BETA_samevs_diff2(:,:,[10:17]),2),3)%sig timepoint with ttest(only)
BETA_samevs_diff3=avdiff2-avsame2;BETA_samevs_diff3=nanmean(nanmean(BETA_samevs_diff3(:,:,[20:26]),2),3);%sig timepoint with ttest(only)
BETA_samevs_diff3=BETA_samevs_diff3*-1;

BETA_samevs_diffcuei=(BETA_samevs_diff2+BETA_samevs_diff3)/2;

[r p]=corrcoef(BETA_samevs_diff2,RT_samevsdiff)

[r p]=corrcoef(BETA_samevs_diff2,summary(:,4))%rt
[r p]=corrcoef(BETA_samevs_diff2,summary(:,1))%sample
[r p]=corrcoef(BETA_samevs_diff2,summary(:,3))%acc
[r p]=corrcoef(BETA_samevs_diff2,summary(:,13))%end same

figure;scatter(BETA_samevs_diff2,summary(:,4))
figure;scatter(BETA_samevs_diff2,summary(:,1))
figure;scatter(BETA_samevs_diff2,summary(:,3))


 mdl = fitlm(BETA_samevs_diff2,summary(:,4))
 
 
[r p]=corrcoef(BETA_samevs_diffcuei,summary(:,4))%rt
[r p]=corrcoef(BETA_samevs_diffcuei,summary(:,1))%sample
[r p]=corrcoef(BETA_samevs_diffcuei,summary(:,3))%acc


%%%%%for norm models need plot_control_analysis
%with norm model etc (taking peak timepoint)
[r p]=corrcoef(Bsame(:,17),summary(:,4))%rt
[r p]=corrcoef(Bsame(:,17),summary(:,1))%sample
[r p]=corrcoef(Bsame(:,17),summary(:,3))%acc
[r p]=corrcoef(Bsame(:,17),summary(:,13))%end same
[r p]=corrcoef(BMwinner(:,20),summary(:,4))%rt
[r p]=corrcoef(BMwinner(:,20),summary(:,1))%sample
[r p]=corrcoef(BMwinner(:,20),summary(:,3))%acc
[r p]=corrcoef(BMwinner(:,20),summary(:,13))%end same

%with norm model etc (taking sig periods )
[r p]=corrcoef(mean(Bsame(:,13:18),2),summary(:,4))%rt
[r p]=corrcoef(mean(Bsame(:,13:18),2),summary(:,1))%sample
[r p]=corrcoef(mean(Bsame(:,13:18),2),summary(:,3))%acc
[r p]=corrcoef(mean(Bsame(:,13:18),2),summary(:,13))%end same
[r p]=corrcoef(mean(BMwinner(:,18:22),2),summary(:,4))%rt
[r p]=corrcoef(mean(BMwinner(:,18:22),2),summary(:,1))%sample
[r p]=corrcoef(mean(BMwinner(:,18:22),2),summary(:,3))%acc
[r p]=corrcoef(mean(BMwinner(:,18:22),2),summary(:,13))%end same

%% correlation of behaviour with norm model etc (fullGLM, as opposed to pairwise betas)
load('summary_patients.mat')
%run lfp_regression_allcontacts_simple
[r p]=corrcoef(mean(Bsame_fullglm(:,15:18),2),summary(:,4))%rt
[r p]=corrcoef(mean(Bsame_fullglm(:,15:18),2),summary(:,1))%sample
[r p]=corrcoef(mean(Bsame_fullglm(:,15:18),2),summary(:,3))%acc
[r p]=corrcoef(mean(Bsame_fullglm(:,15:18),2),summary(:,13))%end same
[r p]=corrcoef(mean(Bsame_fullglm(:,24:27),2),summary(:,4))%rt
[r p]=corrcoef(mean(Bsame_fullglm(:,24:27),2),summary(:,1))%sample
[r p]=corrcoef(mean(Bsame_fullglm(:,24:27),2),summary(:,3))%acc
[r p]=corrcoef(mean(Bsame_fullglm(:,24:27),2),summary(:,13))%end same
[r p]=corrcoef(Bsame_fullglm(:,17),summary(:,4))%rt
[r p]=corrcoef(Bsame_fullglm(:,17),summary(:,1))%sample
[r p]=corrcoef(Bsame_fullglm(:,17),summary(:,3))%acc
[r p]=corrcoef(Bsame_fullglm(:,17),summary(:,13))%end same
[r p]=corrcoef(BEV_fullglm(:,16),summary(:,4))%rt
[r p]=corrcoef(BEV_fullglm(:,16),summary(:,1))%sample
[r p]=corrcoef(BEV_fullglm(:,16),summary(:,3))%acc
[r p]=corrcoef(BEV_fullglm(:,16),summary(:,13))%end same
[r p]=corrcoef(BEV_fullglm(:,24),summary(:,4))%rt
[r p]=corrcoef(BEV_fullglm(:,24),summary(:,1))%sample
[r p]=corrcoef(BEV_fullglm(:,24),summary(:,3))%acc
[r p]=corrcoef(BEV_fullglm(:,24),summary(:,13))%end same
[r p]=corrcoef(mean(BEV_fullglm(:,22:25),2),summary(:,4))%rt
[r p]=corrcoef(mean(BEV_fullglm(:,22:25),2),summary(:,1))%sample
[r p]=corrcoef(mean(BEV_fullglm(:,22:25),2),summary(:,3))%acc
[r p]=corrcoef(mean(BEV_fullglm(:,22:25),2),summary(:,13))%end same
[r p]=corrcoef(BU_fullglm(:,13),summary(:,4))%rt
[r p]=corrcoef(BU_fullglm(:,13),summary(:,1))%sample
[r p]=corrcoef(BU_fullglm(:,13),summary(:,3))%acc
[r p]=corrcoef(BU_fullglm(:,13),summary(:,13))%end same
[r p]=corrcoef(BU_fullglm(:,22),summary(:,4))%rt
[r p]=corrcoef(BU_fullglm(:,22),summary(:,1))%sample
[r p]=corrcoef(BU_fullglm(:,22),summary(:,3))%acc
[r p]=corrcoef(BU_fullglm(:,22),summary(:,13))%end same
[r p]=corrcoef(mean(BU_fullglm(:,11:16),2),summary(:,4))%rt
[r p]=corrcoef(mean(BU_fullglm(:,11:16),2),summary(:,1))%sample
[r p]=corrcoef(mean(BU_fullglm(:,11:16),2),summary(:,3))%acc
[r p]=corrcoef(mean(BU_fullglm(:,11:16),2),summary(:,13))%end same
[r p]=corrcoef(mean(BU_fullglm(:,21:24),2),summary(:,4))%rt
[r p]=corrcoef(mean(BU_fullglm(:,21:24),2),summary(:,1))%sample
[r p]=corrcoef(mean(BU_fullglm(:,21:24),2),summary(:,3))%acc
[r p]=corrcoef(mean(BU_fullglm(:,21:24),2),summary(:,13))%end same
[r p]=corrcoef(BNORM_fullglm(:,21),summary(:,4))%rt
[r p]=corrcoef(BNORM_fullglm(:,21),summary(:,1))%sample
[r p]=corrcoef(BNORM_fullglm(:,21),summary(:,3))%acc
[r p]=corrcoef(BNORM_fullglm(:,21),summary(:,13))%end same

figure;set(gcf,'color','w');
scatter(nanmean(Bsame_fullglm(:,15:18),2),summary(:,4),30,'filled','k')%rt
title('Cue identity (Diff-Same) regressor at 200-350ms', 'FontSize',15);
xlabel ('BETA Power Regression coefficient', 'FontSize',15);
ylabel ('RT', 'FontSize',15);

figure;scatter(Bsame_fullglm(:,17),summary(:,1))%sample
figure;scatter(Bsame_fullglm(:,17),summary(:,3))%acc

%is this caused by an outlier?
isoutlier(Bsame_fullglm(:,17),'quartiles')%good for non-normally distributed data
isoutlier(mean(Bsame_fullglm(:,15:18),2),'quartiles')%good for non-normally distributed data


[r p]=corrcoef(Bsame_fullglm(isoutlier(Bsame_fullglm(:,17))==0,17),summary(isoutlier(Bsame_fullglm(:,17))==0,4))%rt
[r p]=corrcoef(nanmean(Bsame_fullglm(isoutlier(Bsame_fullglm(:,17))==0,15:18),2),summary(isoutlier(Bsame_fullglm(:,17))==0,4))%rt
%15:18;24:27

BETA_samevs_diff=avdiff-avsame;
BETA_samevs_diff_meanacrosschans=squeeze(nanmean(BETA_samevs_diff,2));
for s=1:13;plot(BETA_samevs_diff_meanacrosschans(s,:));hold on;end


%% with coherence
load('COH_TFmeth2long_cuei+1_notsss_cluster2_patients')

COH_samevs_diff=D_betacoh_timecourse-S_betacoh_timecourse;%for the MEGLFP_coherence_final to get these values

COH_samevs_diff=nanmean(COH_samevs_diff(:,sigtimeIND_adj),2);

[r p]=corrcoef(COH_samevs_diff,summary(:,4))%rt
[r p]=corrcoef(COH_samevs_diff,summary(:,1))%sample
[r p]=corrcoef(COH_samevs_diff,summary(:,3))%acc
[r p]=corrcoef(COH_samevs_diff,summary(:,13))%end same#

isoutlier(COH_samevs_diff,'quartiles')%good for non-normally distributed data

figure;set(gcf,'color','w');
scatter(COH_samevs_diff,summary(:,1),30,'filled','k')%rt
title('Coherence STN-dPM', 'FontSize',15);
xlabel ('Coherence Diff-Same', 'FontSize',15);
ylabel ('cues sampled', 'FontSize',15);

%% does RT congruency change depending on stim num (de Lage 2010 JoN)

for s=[1:5,8:15]
RTeffect_stim1(s)=nanmean(summary_dump1(summary_dump5==1 | summary_dump10==s))-nanmean(summary_dump1(summary_dump5==0 | summary_dump10==s));
RTeffect_stim2(s)=nanmean(summary_dump1(summary_dump6==1 | summary_dump10==s))-nanmean(summary_dump1(summary_dump6==0 | summary_dump10==s));
RTeffect_stim3(s)=nanmean(summary_dump1(summary_dump7==1 | summary_dump10==s))-nanmean(summary_dump1(summary_dump7==0 | summary_dump10==s));
RTeffect_stim4(s)=nanmean(summary_dump1(summary_dump8==1 | summary_dump10==s))-nanmean(summary_dump1(summary_dump8==0 | summary_dump10==s));
RTeffect_stim5(s)=nanmean(summary_dump1(summary_dump9==1 | summary_dump10==s))-nanmean(summary_dump1(summary_dump9==0 | summary_dump10==s));
end

figure;set(gcf,'color','w');
B3=barwitherr(cat(2,std(RTeffect_stim1([1:5,8:end]))/sqrt(13),std(RTeffect_stim2([1:5,8:end]))/sqrt(13),std(RTeffect_stim3([1:5,8:end]))/sqrt(13),std(RTeffect_stim4([1:5,8:end]))/sqrt(13),std(RTeffect_stim5([1:5,8:end]))/sqrt(13)),cat(2,mean(RTeffect_stim1([1:5,8:end])),mean(RTeffect_stim2([1:5,8:end])),mean(RTeffect_stim3([1:5,8:end])),mean(RTeffect_stim4([1:5,8:end])),mean(RTeffect_stim5([1:5,8:end]))),'k');
% B3(1).FaceColor=[0 0 0];B3(2).FaceColor=[0.5 0.5 0.5];
ylabel('RT(diff) congruency to chosen direction', 'FontSize',15);xlabel('cue number in sequence', 'FontSize',15);xticklabels({'1','2','3','4','5'}); set(gca,'FontSize',15);


%% regression as in Leimbach 2018

%get stimulus sequence (t)
%get evidence sequence (x)
%get action sequence (a)
%example:
%stimulus: 1 1 2 2 1 1 1 1 1 
%action: 0 0 0 0 0 0 0 0 1
%evidence for chosen: 1 2 1 0 1 2 3 4 5

THRESH=[];
for s=1:size(behav,2)
    
    for i=1:size(behav{s}.session,2)
        
        STIMSEQ=[];ACTIONSEQ=[];EVIDENCESEQ=[];
        
        for n=1:size(behav{s}.session(i).stimuli,2)
            
        stimuliseen=behav{s}.session(i).stimuli{n};
        
        choice=behav{s}.session(i).choice(n);

        action=ones(1,size(stimuliseen,2)); action(end)=2;
%         action=zeros(1,size(stimuliseen,2));action(end)=1;
        
        [evidence_chosen, evidence_left, evidence_current, evidence_absolute, current_left, same, agree, evidence_total] = trial_predictors_stim_by_choice (stimuliseen, choice);
        
        STIMSEQ= [STIMSEQ,1:1:size(stimuliseen,2)];
        ACTIONSEQ= [ACTIONSEQ,action];
        EVIDENCESEQ= [EVIDENCESEQ,evidence_chosen];
%         EVIDENCESEQ= [EVIDENCESEQ,evidence_current];
        
        clear choise stimuliseen action

        end
        
        %logistic regression
        B=mnrfit(cat(2,STIMSEQ',EVIDENCESEQ'),ACTIONSEQ');
        
        %slope of the decrease of decision threshold  taken as arctan(?t/?x). 
        THRESH(s,i)=atan(B(2)/B(3))*180/pi;
        
    end
end

THRESH(THRESH==0)=NaN;

p_THRESH=THRESH;
c_THRESH=THRESH;

figure;set(gcf,'color','w');
subplot(1,2,1)
BAR=barwitherr(nanstd(p_THRESH)/sqrt(size(p_THRESH,1)),nanmean(p_THRESH));BAR(1).FaceColor=[0 0 0];
ylabel('threshold slope', 'FontSize',15);xlabel('session', 'FontSize',15);set(gca,'FontSize',15);yy=ylim;
title('patients');
clear ylim
subplot(1,2,2)
BAR2=barwitherr(nanstd(c_THRESH)/sqrt(size(c_THRESH,1)),nanmean(c_THRESH));BAR2(1).FaceColor=[0.5 0.5 0.5];
ylabel('threshold slope', 'FontSize',15);xlabel('session', 'FontSize',15);set(gca,'FontSize',15);ylim([yy]);
title('controls');


%% calculate P(S) for sample sequences
load('behavioural_data.mat')
STN_norm_term_cuewise=struct('s',15,'i',2,'n',{});
for s=1:size(behav,2)
    
    for i=1:size(behav{s}.session,2)
        
        for n=1:size(behav{s}.session(i).stimuli,2)
            
            stimseen=behav{s}.session(i).stimuli{n};choice= behav{s}.session(i).choice(n);
            [evidence_chosen, evidence_left, evidence_current, evidence_absolute, current_left, same, agree, evidence_total] = trial_predictors_stim_by_choice (stimseen, choice);
            evidence_chosen_end(n) = evidence_chosen(end);
            
            
            pAL=0.5;pAR=0.5;%actions equally likely @ start            
            for j=1:size(stimseen,2)
                    if stimseen(j)==1
                        NORM(j)=pAL*0.7+pAR*0.3;
                        BAYES(j) = pAL*0.7/NORM(j);
                        pAL=BAYES(j);
                        pAR=1-pAL;
                    elseif stimseen(j)==2                        
                        NORM(j)=pAR*0.7+pAL*0.3;
                        BAYES(j) = pAR*0.7/NORM(j);
                        pAR=BAYES(j);
                        pAL=1-pAR;
                    end                    
            end
            
            STN_norm_term(s,i,n) = sum(NORM);
            STN_norm_term_cuewise{s,i,n} = NORM;
            BAYES_final(s,i,n)=BAYES(end);
            
            clearvars -except STN_norm_term STN_norm_term_cuewise s i n behav BAYES_*
        end  
    end
end


STN_norm_term([5,7,15],2,:)=nan;
BAYES_final([5,7,15],2,:)=nan;


%%

load('D:\TaskCode\unbiased_6.mat')
for j=1:1000
for i=1:size(stimseq,1)
    
    stimseq_sampled=stimseq(i,1:datasample(summary_dump2,1,'Replace',true));
    
    if abs(mean(stimseq_sampled) - direction(i)) <0.5
        acc(j,i)=1;
    else
        acc(j,i)=0;
    end

    
end
end