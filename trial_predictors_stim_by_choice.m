function [evidence_chosen, evidence_left, evidence_current, evidence_absolute, current_left, same, agree, evidence_total] = trial_predictors_stim_by_choice (stimuli, choice)

% function [evidence_left, evidence_current, evidence_absolute, current_left, same, agree, evidence_total] = trial_predictors_simple (stimuli)
%
% Function generates predictions on power of overall activity in the STN
% (measured by power of LFP or BOLD) according to the MSPRT model (version
% of Bogacz & Larsen, 2011)
% Input
%  stimuli - a sequence of stimuli: vector of elements: 1 (left) and 2
%            (right)
% Outputs:
%  evidence_left    = # stimuli for left - # stimuli for right
%  evidence_current = # stimuli for current before - #stimuli for other
%  evidence_absolute= absolute value of the evidence_left
%  current_left     = (current stimulus left) - (current stimulus right)
%  same - if the stimulus the same as before:
%        1 - if it is
%        -1 - if it is different
%        0 - if it is the first stimulus of the trial
%  agree - does the stimulus agree with the previous ones
%         1 - if the same as weight of previous
%         -1 - if different than weight of previous
%         0 - if previous stimuli had weight 0
%  evidence_total = # of stimuli so far (like the urgency)

% %evidence for chosen (added by EZP 2019)

nstim = length(stimuli);
evidence_chosen     = zeros(1,nstim);
evidence_left       = zeros(1,nstim);
evidence_current    = zeros(1,nstim);
evidence_absolute   = zeros(1,nstim);
current_left        = zeros(1,nstim);
same                = zeros(1,nstim);
agree               = zeros(1,nstim);
evidence_total      = zeros(1,nstim);

j=1;k=1;
for i = 1:nstim
    current_left(i) = sum(stimuli(i)==1) - sum(stimuli(i)==2);
    evidence_left(i) = sum(stimuli(1:i)==1) - sum(stimuli(1:i)==2);
    evidence_current(i) = sum(stimuli(1:i-1)==stimuli(i)) - sum(stimuli(1:i-1)~=stimuli(i));
    evidence_absolute_left(i) = abs(evidence_left(i));
    
    %evidence for chosen (added by EZP 2019)
    if choice==1 && stimuli(i)==1
        evidence_chosen(i)=j;
        j=j+1;
    elseif choice==1 && stimuli(i)==2
        if i==1
            evidence_chosen(i)=0;
        else
        evidence_chosen(i)=j-2;
        j=j-1;
        end
    elseif choice==2 && stimuli(i)==2
        evidence_chosen(i)=k;
        k=k+1;
    elseif choice==2 && stimuli(i)==1    
        if i==1
            evidence_chosen(i)=0;
        else
        evidence_chosen(i)=k-2;
        k=k-1;
        end
    end
    
    
    if i > 1
        if stimuli(i) == stimuli(i-1)
            same(i) = 1;
        else
            same(i) = -1;
        end
        asymmetry0 = sum(stimuli(1:i-1)==1) - sum(stimuli(1:i-1)==2);
        if abs(asymmetry0) < 0.0001
            agree(i) = 0;
        elseif (asymmetry0 > 0 && stimuli(i) == 1) || (asymmetry0 < 0 && stimuli(i) == 2)
            agree(i) = 1;
        else
            agree(i) = -1;
        end
    end
    evidence_total(i) = i;
end