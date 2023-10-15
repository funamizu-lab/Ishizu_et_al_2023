function FigureS1c_NtrialInfluence

parent  ='G:/Ishizu_data';
outpath ='/Revise_ishizu/output/NtrialInfluence';

% hulistic parameter %
N=10;% N-th prior trails used for linear regression
figsaveTYPE='-dsvg';
%--------------------%

%%% collect the choice data %%%
behavePath='/Tokyo_ephys_ishizu/only_all_behaviors';
cd([parent,behavePath]);
list = dir('*_behave');

% mousename  = cell(length(list),1);
mouseTag   = cell(length(list),1);
Coef = cell(length(list),1);
for i=1:length(list)
    cd(list(i).name);
%     mousename{i,1} = list(i).name(1:3);
    Bpodlist = dir('Bpod*');
    mouseTag{i} = i*ones(length(Bpodlist),1);
     
    [Choice,~,L_Outcome,R_Outcome,Sound] = getChoiceData;
    [Coef{i},~,varname]= LenearRegression(Choice,L_Outcome,R_Outcome,Sound,N);    
    cd([parent,behavePath]);
end

mouseTagWhole = cell2mat(mouseTag);
CoefWhole = cell2mat(Coef');
allSessionNum = 0;
for i=1:length(Coef)
    allSessionNum = allSessionNum+size(Coef{i},2);
end

coeff_model = struct();
for s=1:allSessionNum
    coeff_model(s).session = s;
    coeff_model(s).mouse = mouseTagWhole(s);
    for v=1:length(varname)-1
        eval(['coeff_model(s).',varname{v},'=CoefWhole(v,s);']);
    end
end

%%% kentei %%%
model_result=getCoeffPval(coeff_model);

%%% save folder %%%
if(~exist([parent,outpath],'dir')), mkdir([parent,outpath]); end
cd([parent,outpath]);
save('savemat.mat','coeff_model','model_result','CoefWhole','N');

%% plot figure model1 %%
close all;
drawResults(coeff_model,model_result,N,'model',figsaveTYPE);

end

%-------------------------------------------------------------------------%
function drawResults(coeff_model,model_result,N,figname,figsaveTYPE)
%-------------------------------------------------------------------------%
%%% plot the each mouse data %%% 
f=fieldnames(coeff_model);
label=f(3:length(f));
eachlabelID_start=1:N:length(label);
eachlabelID_start(end)=[];
labelTYPE=label(eachlabelID_start);
for i=1:length(labelTYPE)
    labelTYPE{i}=labelTYPE{i}(1:end-1);
end

mouseval=zeros(length(coeff_model),1);
for i=1:length(coeff_model)
    mouseval(i)=coeff_model(i).mouse;
end

ticklabel=cell(1,N+1);
for i=1:N+1
    ticklabel{i}=num2str(i-1);
end

%%% plot the all mouse data %%%
estVal=cell2mat(model_result.estVal);
Pval  =cell2mat(model_result.Pval);
CIs   =cell2mat(model_result.CIs);
h=figure('Position',[20,200,2000,500]);
for j=1:length(labelTYPE)
    if(j==length(labelTYPE))
        k=1; tl=ticklabel;
    else
        k=0; tl=ticklabel(2:end);
    end
    id=eachlabelID_start(j):(eachlabelID_start(j)+N-1+k);
    Coef = estVal(id);
    P = Pval(id);
    Ps= floor(-log10(P));
    Ps(Ps==1)= P(Ps==1);
    plabel=[];
    for i=1:length(Ps)
        if(i==1)
            plabel=num2str(Ps(i));
        else
            plabel=[plabel,' / ',num2str(Ps(i))];
        end
    end
    ci= CIs(id,:);
    upper= ci(:,2)-Coef;
    lower= Coef-ci(:,1);
    subplot(1,length(labelTYPE),j);
    hold on;
    errorplot(1:length(Coef),Coef,upper,lower,'k',0.1,1.5);
    plot([0,N+1+k],[0,0],'--k','LineWidth',0.5);
    xlim([0,N+1+k]);
    xticks(1:N+k); xticklabels(tl);
    xlabel('N-th trail prior to the current choice');
    ylabel('Coefficient');
    title({labelTYPE{j},['-log10(Pval): ',plabel]});    
end
set(h,'PaperPositionMode','auto');
print(h,'-r0',['all mouse ',figname],figsaveTYPE);
end

%-------------------------------------------------------------------------%
function output=getCoeffPval(coeff_model)
%-------------------------------------------------------------------------%
mouseval=zeros(length(coeff_model),1);
sessionval=zeros(length(coeff_model),1);
for i=1:length(coeff_model)
    mouseval(i)=coeff_model(i).mouse;
    sessionval(i)=coeff_model(i).session;
end

f=fieldnames(coeff_model);
name=f(3:length(f));
estVal =cell(length(f)-2,1);
Pval   =cell(length(f)-2,1);
CIs =cell(length(f)-2,2);
for i=3:length(f)
    data=zeros(length(coeff_model),1);
    for k=1:length(coeff_model)
        eval(['data(k)=coeff_model(k).',f{i},';']);
    end
    tbl = table(data,mouseval,sessionval,'VariableNames',{'value','mouse','session'});
    lme1 = fitlme(tbl,'value ~ 1 + (1|mouse)+(1|session)+(mouse-1|session)'); %linear regression
    lme2 = fitlme(tbl,'value ~ 1 + (1|mouse)+(1|session)'); %linear regression
    if(lme1.ModelCriterion.BIC < lme2.ModelCriterion.BIC)
        lme=lme1;
    else
        lme=lme2;
    end
    estVal{i-2}= lme.Coefficients.Estimate;
    Pval{i-2}  = lme.Coefficients.pValue;
    CIs{i-2,1} = lme.Coefficients.Lower;
    CIs{i-2,2} = lme.Coefficients.Upper;
end

output.name  = name;
output.estVal= estVal;
output.Pval  = Pval;
output.CIs   = CIs;
end

%-------------------------------------------------------------------------%
function [Choice,Reward,L_Outcome,R_Outcome,Sound] = getChoiceData
%-------------------------------------------------------------------------%
list = dir('Bpod_*');

Choice=cell(length(list),1);
Reward=cell(length(list),1);
L_Outcome=cell(length(list),4);% left reward / left error / left small reward / left big reward
R_Outcome=cell(length(list),4);% right reward / right error / right small reward / right big reward
Sound=cell(length(list),1);
for i=1:length(list)
    load(list(i).name,'Outcome','Chosen_side','Correct_side','TrialBlock','BlockReward','EvidenceStrength','Tone_cloud');
    [RewardAmount,LeftReward,RightReward,smallLeftReward,smallRightReward,...
        bigLeftReward,bigRightReward,LeftError,RightError,binary_tone]=...
        getToneStructure(Outcome,Chosen_side,Correct_side,TrialBlock,BlockReward,EvidenceStrength,Tone_cloud);
    Choice_trial  = find(Outcome == 1 | Outcome == 2);
    useBlock_trial= find(TrialBlock>1);
    use_trial = intersect(Choice_trial,useBlock_trial);
    
    Choice{i}= Chosen_side(use_trial);
    Reward{i}= RewardAmount(use_trial);
    L_Outcome{i,1}= LeftReward(use_trial);
    L_Outcome{i,2}= LeftError(use_trial);
    L_Outcome{i,3}= smallLeftReward(use_trial);
    L_Outcome{i,4}= bigLeftReward(use_trial);
    R_Outcome{i,1}= RightReward(use_trial);
    R_Outcome{i,2}= RightError(use_trial);
    R_Outcome{i,3}= smallRightReward(use_trial);
    R_Outcome{i,4}= bigRightReward(use_trial);
    Sound{i}= binary_tone(use_trial);
end
end

%-------------------------------------------------------------------------%
function [RewardAmount,LeftReward,RightReward,smallLeftReward,smallRightReward,...
    bigLeftReward,bigRightReward,LeftError,RightError,binary_tone]=...
    getToneStructure(Outcome,Chosen_side,Correct_side,TrialBlock,BlockReward,EvidenceStrength,Tone_cloud)
%-------------------------------------------------------------------------%
% %Outcome
% outcome_EW     = 0; %early withdrawal
% outcome_IC     = 1; %incorrect choice
% outcome_reward = 2; %reward was dispensed (either automatically in early training, or after correct choice)
% outcome_NC     = 3; %no choice was made and time elapsed
% outcome_UN     = 4; %undefined or Free water:

%%% reward %%%
left  = find(Chosen_side == 0);
right = find(Chosen_side == 1);
error=find(Outcome==1);
correct=find(Outcome==2);
blocknum=unique(TrialBlock);
BlockReward2=zeros(length(blocknum),2);% reward size:[left/right]
for i=1:2
    BlockReward2(mod(blocknum,2)==1,i)=BlockReward(3,i);
    BlockReward2(mod(blocknum,2)==0,i)=BlockReward(2,i);
end
BlockReward2(1,:)=BlockReward(1,:);
LeftReward =zeros(length(TrialBlock),1); 
RightReward=zeros(length(TrialBlock),1); 
LeftReward(intersect(correct, left)) =BlockReward2(TrialBlock(intersect(correct,left)),1);
RightReward(intersect(correct,right))=BlockReward2(TrialBlock(intersect(correct,right)),2);
RewardAmount=LeftReward+RightReward;

LRamount = unique(LeftReward);
RRamount = unique(RightReward);
LRamount(LRamount==0)=[]; 
RRamount(RRamount==0)=[];
smallLeftReward =zeros(length(LeftReward),1);
smallRightReward=zeros(length(RightReward),1);
smallLeftReward(LeftReward==min(LRamount))=min(LRamount);
smallRightReward(RightReward==min(RRamount))=min(RRamount);
bigLeftReward =zeros(length(LeftReward),1);
bigRightReward=zeros(length(RightReward),1);
bigLeftReward(LeftReward==max(LRamount))=max(LRamount);
bigRightReward(RightReward==max(RRamount))=max(RRamount);

% scaling data to [0 1]%
% LeftReward = rescale(LeftReward);
% RightReward= rescale(RightReward);
% RewardAmount=rescale(RewardAmount);
    
LeftError =zeros(length(TrialBlock),1);
RightError=zeros(length(TrialBlock),1);
LeftError(intersect(error,left)) =1;
RightError(intersect(error,right)) =1;

%%% sound (evidence) %%%
temp_evi = unique(EvidenceStrength);
temp_evi_low  = 0.5 - temp_evi/2;
temp_evi_high = 0.5 + temp_evi/2;

%Put tone evidence in all trials;
trial_evidence = zeros(length(Outcome),1);
low  = find(Correct_side == 0);
high = find(Correct_side == 1);
for i = 1:length(temp_evi)
    temp = find(EvidenceStrength == temp_evi(i));
    temp_left  = intersect(temp,low);
    temp_right = intersect(temp,high);
    trial_evidence(temp_left)  = temp_evi_low(i);
    trial_evidence(temp_right) = temp_evi_high(i);
end

%Make the true tone cloud value
binary_tone=zeros(length(Tone_cloud),1);
for i = 1:length(Tone_cloud)
    temp_tone = Tone_cloud(i).matrix;
    %Get the data in all sound
    temp1 = find(temp_tone >= 9);
    binary_tone(i) = length(temp1) ./ length(temp_tone);
end

%Based on the correct trial, flip the tone cloud
if mean(binary_tone(Correct_side == 1)) < 0.5 %low for right correct
    % flipping tones %
    binary_tone = 1 - binary_tone;
else
end
end

%-------------------------------------------------------------------------%
function  [Coef,Pval,varname]= LenearRegression(Choice,L_Outcome,R_Outcome,Sound,N)
%-------------------------------------------------------------------------%
sessionNum = length(Choice);
Coef=zeros(5*N+1,sessionNum);
Pval=zeros(5*N+1,sessionNum);
for i=1:sessionNum
    dataC=Choice{i};
    dataLR=L_Outcome{i,1};% Left Reward
    dataLE=L_Outcome{i,2};% Left Error
    dataRR=R_Outcome{i,1};% Right Reward
    dataRE=R_Outcome{i,2};% Right Error
    dataS=Sound{i};
    varname=cell(1,5*N+2);
    
    % data X (predictor) %
    Xlr = zeros(length(dataC)-N,N);
    Xle = zeros(length(dataC)-N,N);
    Xrr = zeros(length(dataC)-N,N);
    Xre = zeros(length(dataC)-N,N);
    Xs1 = zeros(length(dataC)-N,N);
    for k=1:N
        Xlr(:,k) = dataLR(N+1-k:end-k);
        Xle(:,k) = dataLE(N+1-k:end-k);
        Xrr(:,k) = dataRR(N+1-k:end-k);
        Xre(:,k) = dataRE(N+1-k:end-k);
        Xs1(:,k) = dataS(N+1-k:end-k);
        varname{k}=['leftR',num2str(k)];
        varname{k+N}=['leftE',num2str(k)];
        varname{k+2*N}=['rightR',num2str(k)];
        varname{k+3*N}=['rightE',num2str(k)];
        varname{k+4*N+1}=['sound',num2str(k)];
    end
    varname{4*N+1}='sound0';
    Xs=[dataS(N+1:end),Xs1];
    X=[Xlr,Xle,Xrr,Xre,Xs];
    
    % data y (explain) %
    varname{end}='currentChoice';
    y = dataC(N+1:end);
    
    % linear regression %
    mdl = fitglm(X,y,'Distribution','Binomial','Link','logit','VarNames',varname);
    Coef(:,i) = mdl.Coefficients.Estimate(2:end);
    Pval(:,i) = mdl.Coefficients.pValue(2:end);
end
end