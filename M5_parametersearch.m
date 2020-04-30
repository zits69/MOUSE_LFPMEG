

function [BONUSLFP,  ALPHA, LogL5, trialnum] = M5_parametersearch(P)

    [BONUS,  LogL3, trlnum3] = M3_parametersearch(1);

    %load trialwise data
    load('D:\Data\DBS-MEG\beta_trialwise_weightsM5.mat');
    %zscore across group
    beta_group=[];for i=[1:5,8:15];beta_group=[beta_group,beta_trlwise{i}'];beta_trl_subj(i)=size(beta_trlwise{i},1);end
    beta_group_z=zscore(beta_group);
    k=1;for i=[1:5,8:15];beta_group_z_subj{i}=beta_group_z(1,k:k+beta_trl_subj(i)-1);k=k+beta_trl_subj(i);disp(i),disp(k);end

    if P==1
        load('D:\behavioural_data.mat');
    else
        load('D:\behavioural_data_controls.mat');
    end

        for s=[1:5,8:15]

            %combine 2 sessions into one for estimation
            if s==1 || s== 10 || s==5 || s==7 || s==9 %either 1 sess only or 2nd was rejected %s9 trial 101 has 7 stim but only 6 got recorded for seq so only do 1st sess because stim don't line up
                combosession=struct();
                combosession.acc=behav{s}.session(1).acc;
                combosession.RT=behav{s}.session(1).RT;
                combosession.choice=behav{s}.session(1).choice;
                combosession.stimuli=behav{s}.session(1).stimuli;
                combosession.stimnum=behav{s}.session(1).stimnum;
            elseif s==6
                combosession=struct();
                combosession.acc=behav{s}.session(2).acc;
                combosession.RT=behav{s}.session(2).RT;
                combosession.choice=behav{s}.session(2).choice;
                combosession.stimuli=behav{s}.session(2).stimuli;
                combosession.stimnum=behav{s}.session(2).stimnum;
            elseif s==15
                combosession=struct();
                combosession.acc=behav{s}.session(1).acc(1:end-2);
                combosession.RT=behav{s}.session(1).RT(1:end-2);
                combosession.choice=behav{s}.session(1).choice(1:end-2);
                combosession.stimuli=behav{s}.session(1).stimuli(1:end-2);
                combosession.stimnum=behav{s}.session(1).stimnum(1:end-2);
            else
                combosession=struct();
                combosession.acc=[behav{s}.session(1).acc,behav{s}.session(2).acc];
                combosession.RT=[behav{s}.session(1).RT,behav{s}.session(2).RT];
                combosession.choice=[behav{s}.session(1).choice,behav{s}.session(2).choice];
                combosession.stimuli=[behav{s}.session(1).stimuli,behav{s}.session(2).stimuli];
                combosession.stimnum=[behav{s}.session(1).stimnum,behav{s}.session(2).stimnum];
            end
            trialnum(s)=size(combosession.acc,2);

            BONUSLFP = beta_group_z_subj{s};
            bonus=BONUS(s);

            disp(s)
            %initiate search for values
            params = fminsearch (@negcor, 1, [], combosession, BONUSLFP, bonus);

            ALPHA (s) = params(1);

            %run DV models with L estimates
            sb=1;
            for n=1:size(combosession.stimuli,2)

                stimuliseen = combosession.stimuli{n};
                choice(n) = combosession.choice(n);

                [evidence_chosen, evidence_left, evidence_current, evidence_absolute, current_left, same, agree, evidence_total] = trial_predictors_stim_by_choice (stimuliseen, choice(n));

                %M5A
                %             DV5(n)=0;
                %             for ss=1:size(current_left,2)
                %
                %                     DV5(n) = DV5(n) + (1+ALPHA(s)*BONUSLFP(sb))*current_left(ss)*-1;
                %
                %             sb=sb+1;%counter for LFP bonus which is cue-wise
                %             end

                %M5B
                DV5(n)=0;
                for ss=1:size(current_left,2)
                    if same(ss)==1
                        DV5(n) = DV5(n) + (bonus+ALPHA(s)*BONUSLFP(sb))*current_left(ss)*-1;
                    else
                        DV5(n) = DV5(n) + current_left(ss)*-1;
                    end
                    sb=sb+1;%counter for LFP bonus which is cue-wise
                end



            end

            %estimate logL within subj along with the lambda and weights
            Y5=pdf('Logistic',DV5,choice); LogL5(s) = -sum(log(Y5));

            fprintf('.');

        end
end

function [R5]= negcor (params, combosession, BONUSLFP,bonus)


        alpha = params(1);

        if alpha < 0 || alpha > 10
            R5 = 1000000000;%make sure fminsearch is bounded to certain values
        else
            sb=1;
            for n=1:size(combosession.stimuli,2)

                stimuliseen = combosession.stimuli{n};
                choice(n) = combosession.choice(n);

                [evidence_chosen, evidence_left, evidence_current, evidence_absolute, current_left, same, agree, evidence_total] = trial_predictors_stim_by_choice (stimuliseen,  choice(n));

                %M5A
                %                     DV5(n)=0;
                %                     for ss=1:size(current_left,2)
                %
                %                         DV5(n) = DV5(n) + (1+alpha*BONUSLFP(sb))*current_left(ss)*-1;
                %
                %                     sb=sb+1;%counter for LFP bonus which is cue-wise
                %                     end

                %M5B
                DV5(n)=0;
                for ss=1:size(current_left,2)

                    if same(ss)==1
                        DV5(n) = DV5(n) + (bonus+alpha*BONUSLFP(sb))*current_left(ss)*-1;
                    else
                        DV5(n) = DV5(n) + current_left(ss)*-1;
                    end

                    sb=sb+1;%counter for LFP bonus which is cue-wise
                end
            end

        
        %                 disp(sb)

        %remove choice == 3 (mistake)
        rm=find(choice==3);choice(rm)=[];DV5(rm)=[];

        temp = -corrcoef(cat(2,DV5',choice'));R5 = temp(1,2);
        %            temp = - mnrfit(DV2,choice);R2 = temp(2,1);
      end
end
