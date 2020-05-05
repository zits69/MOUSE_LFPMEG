clear
%% user input

P=0;%1 for PATIENTS, 0 for CONTROLS

% %should figured be displayed?

fig=0;%recreates individual plots like in Leimbach 2018, JoCN
sumfig=1;%group level plots

%% load data
if P==1
    load('behavioural_data.mat');new_RT_cutoff=207;%PATIENTS - these cutoffs were established by looking for very fast RTs and by relating the choice to the last stimulus presented (see: line 272)
else
    load('behavioural_data_controls.mat');new_RT_cutoff=234;%CONTROLS
end
%%

summary_dump1=[];summary_dump2=[];summary_dump3=[];summary_dump4=[];summary_dump5=[];summary_dump6=[];summary_dump7=[];summary_dump8=[];summary_dump9=[];summary_dump10=[];
summary_dump_bayes=[];

for s=1:size(behav,2)
    
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
            temp2=find(counterT(end,:)==ALLSTIM(n));
       
            scatter(ALLSTIM(n),CHOICE(n),counter(temp,temp2)*10,c,'filled')
            hold on
            
        end
        errorbar(mean(ALLSTIM),mean(CHOICE),std(CHOICE),std(CHOICE),std(ALLSTIM),std(ALLSTIM),'or','CapSize',0);
        xlabel('# Stimuli seen', 'FontSize',15);ylabel('Evidence for Chosen', 'FontSize',15);set(gca,'FontSize',15);%xlim([1,x]);ylim([1,x]);
        title([wanted{s},' session ',num2str(i)]);
        end
        

% % get RTs based on elapsed time and decide if 'last' stimulus really last or t-1
        if P==1 && s==1 && i==2
            behav{s}.session(i).newRTs = behav{s}.session(i).RT - (behav{s}.session(i).stimnum -1)*1000;
        else
            behav{s}.session(i).newRTs = behav{s}.session(i).RT - (behav{s}.session(i).stimnum -1)*800;
        end
        behav{s}.session(i).last(behav{s}.session(i).newRTs <= 0) = -1;
        behav{s}.session(i).last(behav{s}.session(i).newRTs > 0) = 1;
     
        
% % match LFP analysis - same vs different analysis
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
            
            bayes_given_stim=mean(stimuli); %what a bayesian agent would do given the stimuli sampled
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
        
    
% %output into summary and also dump everything  into trialwise summary   
        summary(s,i,1) = nanmean(behav{s}.session(i).stimnum);
        summary(s,i,2) = nanstd(behav{s}.session(i).stimnum);
        summary(s,i,3) = nanmean(behav{s}.session(i).acc);     
        summary(s,i,4) = nanmean(behav{s}.session(i).newRTs);     
        summary(s,i,5) = nanstd(behav{s}.session(i).newRTs); 
        summary(s,i,6) = ttest(behav{s}.session(i).acc,0.5);  
  
% % remove very fast RTs
        RT_excl=find(RT<=new_RT_cutoff);RT(RT_excl)=[];end_same_RT=end_same;end_same_RT(RT_excl)=[];
        summary(s,i,17)=size(RT_excl,2);
        
        summary(s,i,7)=nanmean(RT(end_same_RT==0));summary(s,i,8)=nanmean(RT(end_same_RT==1));%RT x whether end is on 'same'
        
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

        
        summary_dump1 = [summary_dump1, behav{s}.session(i).newRTs]; 
        summary_dump2 = [summary_dump2, behav{s}.session(i).stimnum]; 
        summary_dump3 = [summary_dump3, end_same]; 
        summary_dump4 = [summary_dump4, prob_choice_laststim];
        summary_dump5 = [summary_dump5, choice_match1stim];
        summary_dump6 = [summary_dump6, choice_match2stim];
        summary_dump7 = [summary_dump7, choice_match3stim];
        summary_dump8 = [summary_dump8, choice_match4stim];
        summary_dump9 = [summary_dump9, choice_match5stim];
        summary_dump10 = [summary_dump10, ones(1,size(choice_match5stim,2))*s];
        summary_dump_bayes = [summary_dump_bayes, acc_bayes];
        clearvars -except summary* s i behav *fig* wanted new_RT_cutoff P
        
    end
    
end
%remove those without 2ndsession (1 session is 100 trials) 
summary(summary==0)=nan;

%get trialnums
summarytrlnum=nansum(cat(2,summary(:,1,18),summary(:,2,18)),2);

%avg over sessions
summary=squeeze(nanmean(summary(:,:,1:19),2));

%remove s6/7 as at chance performance (summary column 6 is ttest of accuracy)
if P==1
summary([6,7],:)=[];
summarytrlnum([6,7])=[];
end

%% plots

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
