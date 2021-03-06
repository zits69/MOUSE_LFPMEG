%load in behavioural data and calculate decision variables for incoming
%stimuli trains; input 'P' == 1 for doing on patient data, and 'meth'
%determined the type of DV calculation to display (display only as
%summary)
function [LogL1, trialnum] = M1_parametersearch(P)

if P==1
    load('behavioural_data.mat');
else
    load('behavioural_data_controls.mat');
end

for s=1:size(behav,2)
    
    if P==1
        %combine 2 sessions into one for estimation
        if s==5 || s==7 || s==15
            combosession=struct();
            combosession.acc=behav{s}.session(1).acc;
            combosession.RT=behav{s}.session(1).RT;
            combosession.choice=behav{s}.session(1).choice;
            combosession.stimuli=behav{s}.session(1).stimuli;
            combosession.stimnum=behav{s}.session(1).stimnum;
        else
            combosession=struct();
            combosession.acc=[behav{s}.session(1).acc,behav{s}.session(2).acc];
            combosession.RT=[behav{s}.session(1).RT,behav{s}.session(2).RT];
            combosession.choice=[behav{s}.session(1).choice,behav{s}.session(2).choice];
            combosession.stimuli=[behav{s}.session(1).stimuli,behav{s}.session(2).stimuli];
            combosession.stimnum=[behav{s}.session(1).stimnum,behav{s}.session(2).stimnum];
        end
    elseif P~=1 && s==13
        combosession=struct();
        combosession.acc=[behav{s}.session(1).acc(51:100),behav{s}.session(2).acc,behav{s}.session(3).acc];
        combosession.RT=[behav{s}.session(1).RT(51:100),behav{s}.session(2).RT,behav{s}.session(3).RT];
        combosession.choice=[behav{s}.session(1).choice(51:100),behav{s}.session(2).choice,behav{s}.session(3).choice];
        combosession.stimuli=[behav{s}.session(1).stimuli{51:100},behav{s}.session(2).stimuli,behav{s}.session(3).stimuli];
        combosession.stimnum=[behav{s}.session(1).stimnum(51:100),behav{s}.session(2).stimnum,behav{s}.session(3).stimnum];
    else
        combosession=struct();
        combosession.acc=[behav{s}.session(1).acc,behav{s}.session(2).acc];
        combosession.RT=[behav{s}.session(1).RT,behav{s}.session(2).RT];
        combosession.choice=[behav{s}.session(1).choice,behav{s}.session(2).choice];
        combosession.stimuli=[behav{s}.session(1).stimuli,behav{s}.session(2).stimuli];
        combosession.stimnum=[behav{s}.session(1).stimnum,behav{s}.session(2).stimnum];
    end
    
    trialnum(s)=size(combosession.acc,2);
    
    %run DV models with L estimates
    for n=1:size(combosession.stimuli,2)
        
        stimuliseen = combosession.stimuli{n};
        choice(n) = combosession.choice(n);
        
        [evidence_chosen, evidence_left, evidence_current, evidence_absolute, current_left, same, agree, evidence_total] = trial_predictors_stim_by_choice (stimuliseen, choice(n));
        
        DV1(n)=0;DVsame(n)=0;DVsameend(n)=0;DVdiff(n)=0;DVdiffend(n)=0;
        
        for ss=1:size(current_left,2)
            
            DV1(n)= DV1(n) + current_left(ss)*-1;%M1
            
        end
        
      
        
    end
    
    %estimate logL within subj along with the lambda and weights
    Y1=pdf('Logistic',DV1,choice); LogL1(s) = -sum(log(Y1));
    
    
    fprintf('.');
    
end


end


