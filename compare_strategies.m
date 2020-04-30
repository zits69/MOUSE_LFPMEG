function [integrators, samers, combo]=compare_strategies (P)

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
    
    DV = [];
    decide = [];
    
    %run DV models with L estimates
    for n=1:size(combosession.stimuli,2)
        
        stimuliseen = combosession.stimuli{n};
        choice = combosession.choice(n);
        
        [evidence_chosen, evidence_left, evidence_current, evidence_absolute, current_left, same, agree, evidence_total] = trial_predictors_stim_by_choice (stimuliseen, choice);
        
        DVch=0;DVev=0;
        num_stim = size(current_left,2);
        
        for ss=1:num_stim
            
            %RAFAL
            %based on choice (could also use 'evidence chosen' variable, i
            %added this at som epoint to the function)
            if stimuliseen(ss) == choice
                DVch = DVch + 1;
            else
                DVch = DVch - 1;
            end

            %Z from previous scripts (eg M1)
            %based on incoming evidence
            DVev = DVev + current_left(ss)*-1;%M1 
%             DVev = DVev + evidence_current(ss);%M1 
           
            if ss>1
%                 DV = [DV, [DVev; DVch; same(ss)]];%FULL MODEL
                DV = [DV, [DVch; same(ss)]];% RB RESPONSE MODEL
                if ss == num_stim
                    decide = [decide, 2];
                else
                    decide = [decide, 1];
                end
            end
        end
               
    end
    
    %estimate logL within subj along with the lambda and weights
    [~,~,stats] = mnrfit(DV',decide');
    pev(s) = stats.p(2);
%     pch(s) = stats.p(3);
    psame(s) = stats.p(3);
    fprintf('.');
    
end
fprintf('\nPatients relying only on integration\n');
disp (find (pev<0.05 & psame >= 0.05));integrators=find (pev<0.05 & psame >= 0.05);
fprintf('\nPatients relying only on responding after same stimuli\n');
disp (find (pev>=0.05 & psame < 0.05));samers=find (pev>=0.05 & psame < 0.05);
fprintf('\nPatients relying on integration with higher chance of responding after same\n');
disp (find (pev<0.05 & psame < 0.05));combo=find (pev<0.05 & psame < 0.05);
fprintf('\nPatients without a significant criterion\n');
disp (find (pev>=0.05 & psame >= 0.05));

