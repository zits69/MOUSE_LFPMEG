%% HOW CAN WE BEST DESCRIBE BEHAVIOUR IN THE EXPANDED JUDGEMENT TASK?

%%DV = decision variable, X = stimulus (1 for chosen and -1 else), LAMBDA =
%%forgetting term (range:0-1), W = weight if 'same' stimulus as
%%before(range: 1 or alpha)

%%M1: DVt = DVt-1 + Xt
%%M2: DVt = (1-LAMBDA)DVt-1 + Xt
%%M3: DVt = DVt-1 + WtXt
%%M4: DVt = (1-LAMBDA)DVt-1 + WtXt

%% input
%if P ==1 then use patients, otherwise control data
P=1;
%% FMINSEACRH CALCULATIONS
[LogL1, trlnum1] = M1_parametersearch(P);
[LAMBDA,  LogL2, trlnum2] = M2_parametersearch(P);
[BONUS,  LogL3, trlnum3] = M3_parametersearch(P);
[LAMBDA_m4, BONUS_m4, LogL4, trlnum4] = M4_parametersearch(P);

% % % removeS43/46 (chance performance patients
if P==1
LogL1([6:7])=[];LogL2([6:7])=[];LogL3([6:7])=[];LogL4([6:7])=[];
trlnum1([6:7])=[];trlnum2([6:7])=[];trlnum3([6:7])=[];trlnum4([6:7])=[];
else
end
%% model comparison

models={LogL1,LogL2,LogL3,LogL4};
params=[1,2,2,3,1];
trlnums={trlnum1,trlnum2,trlnum3,trlnum4};

for i=1:4
    
logL=models{i};
numParam=params(i);
numObs=trlnums{i};

[aic(i,:),bic(i,:)] = aicbic(logL,numParam,numObs);

end

for s=1:size(bic,2)
winning_model(s) = find(bic(:,s)==min(bic(:,s),[],1));
end

[H,P,STATS]=chi2gof(winning_model);


