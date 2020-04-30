
function [BONUS,  LogL3, trialnum] = M3_parametersearch(P)


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
        
        %initiate search for values
        params = fminsearch (@negcor, 1, [], combosession);
        
        BONUS(s) = params(1);
      
        %run DV models with L estimates
        for n=1:size(combosession.stimuli,2)
            
            stimuliseen = combosession.stimuli{n};
            choice(n) = combosession.choice(n);
            
            [evidence_chosen, evidence_left, evidence_current, evidence_absolute, current_left, same, agree, evidence_total] = trial_predictors_stim_by_choice (stimuliseen, choice(n));
                        
            %M3
            DV3(n)=0;
            for ss=1:size(current_left,2)
                if same(ss)==1
                    DV3(n) = DV3(n) + BONUS(s)*current_left(ss)*-1;
                else
                    DV3(n) = DV3(n) + 1*current_left(ss)*-1; %no weight
                end
            end
                    
        end
        
       %estimate logL within subj along with the lambda and weights
        Y3=pdf('Logistic',DV3,choice); LogL3(s) = -sum(log(Y3));

        fprintf('.');

    end
end

function [R3]= negcor (params, combosession)

    bonus = params(1);
    
    if bonus < 0 || bonus > 10 
        R3 = 1000000000;%make sure fminsearch is bounded to certain values
    else
       
                for n=1:size(combosession.stimuli,2)
                    
                    stimuliseen = combosession.stimuli{n};
                    choice(n) = combosession.choice(n);
                    
                    [evidence_chosen, evidence_left, evidence_current, evidence_absolute, current_left, same, agree, evidence_total] = trial_predictors_stim_by_choice (stimuliseen,  choice(n));
                    
                    DV3(n)=0;
                    for ss=1:size(current_left,2)
                        if same(ss)==1
                            DV3(n) = DV3(n) + bonus*current_left(ss)*-1;
                        else
                            DV3(n) = DV3(n) + 1*current_left(ss)*-1; %no weight
                        end
                    end

                end
                
         %remove choice == 3 (mistake)
         rm=find(choice==3);choice(rm)=[];DV3(rm)=[];
                                       
         temp = -corrcoef(cat(2,DV3',choice'));R3 = temp(1,2);   
%            temp = - mnrfit(DV2,choice);R2 = temp(2,1);   

    end
end
 