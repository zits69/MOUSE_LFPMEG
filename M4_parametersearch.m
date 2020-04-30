
function [LAMBDA_m4, BONUS_m4, LogL4, trialnum] = M4_parametersearch(P)

    
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
        params = fminsearch (@negcor, [0,1], [], combosession);
        
        LAMBDA_m4(s) = params(1);
        BONUS_m4(s) = params(2);
        
        %run DV models with L&W estimates
        for n=1:size(combosession.stimuli,2)
            
            stimuliseen = combosession.stimuli{n};
            choice(n) = combosession.choice(n);
            
            [evidence_chosen, evidence_left, evidence_current, evidence_absolute, current_left, same, agree, evidence_total] = trial_predictors_stim_by_choice (stimuliseen,  choice(n));
            
            DV4(n)=0;
            for ss=1:size(current_left,2)
            
                if same(ss)==1
                    DV4(n)= (1-LAMBDA_m4(s))*DV4(n) + BONUS_m4(s)*(current_left(ss)*-1);%M4
                else
                    DV4(n)= (1-LAMBDA_m4(s))*DV4(n) + 1*(current_left(ss)*-1);%M4
                end
                
            end
            
        end
        
       %estimate logL within subj along with the lambda and weights
        Y4=pdf('Logistic',DV4,choice); LogL4(s) = -sum(log(Y4));

        fprintf('.');

    end
end

function [R4]= negcor (params, combosession)

    lambda = params(1);
    bonus = params(2);
    
    if lambda > 1 || bonus < 0 || bonus > 10 
        R4 = 1000000000;%make sure fminsearch is bounded to certain values
    else
       
                for n=1:size(combosession.stimuli,2)
                    
                    stimuliseen = combosession.stimuli{n};
                    choice(n) = combosession.choice(n);
                    
                    [evidence_chosen, evidence_left, evidence_current, evidence_absolute, current_left, same, agree, evidence_total] = trial_predictors_stim_by_choice (stimuliseen, choice(n));
                                       
              
                    %M4
                    DV4(n)=0;                    
                    for ss=1:size(current_left,2)  
                        if same(ss)==1
                            DV4(n) = (1-lambda) * DV4(n) + bonus*current_left(ss)*-1;
                        else
                            DV4(n) = (1-lambda) * DV4(n) + 1*current_left(ss)*-1; %no weight
                        end
                    end
                    
                end
                
           %remove choice == 3 (mistake)
           rm=find(choice==3);choice(rm)=[];DV4(rm)=[];
                                       
           temp = -corrcoef(cat(2,DV4',choice'));R4 = temp(1,2);     
%            temp = - mnrfit(DV4,choice);R4 = temp(2,1);     
           

    end
end
 