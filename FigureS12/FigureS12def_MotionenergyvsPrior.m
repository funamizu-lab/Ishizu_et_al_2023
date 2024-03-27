function FigureS12def_MotionenergyvsPrior

parent='G:/Ishizu_data';
workpath='/Revise_ishizu';
outpath='/Revise_ishizu/output/MotionvsPrior2';

figsaveTYPE='-dsvg';
%--------------------%

cd([parent,workpath]);

%%% collect the data %%%
% folder setting %
Dir1 = {'AC','FOF','mPFC'};
Dir2 = {'PPC','posPPC'};
mouse= {'i24','i34','i35','i43'};% mouse no.4,5,6,7
dataPath='/movie';

disp('data collecting...');
folderdata1=cell(length(Dir1),1);
folderlength1=zeros(length(Dir1),1);
for i=1:length(Dir1)
    folderdata1{i} = getadress(Dir1{i},mouse,parent,dataPath);
    folderlength1(i)=length(folderdata1{i});
end

% erase the overlapped data due to AC-PPC simultaneous recording % 
folderdata2=cell(length(Dir2),1);
folderlength2=zeros(length(Dir2),1);
for i=1:length(Dir2) 
    folderdata2{i}  = getadress(Dir2{i},mouse,parent,dataPath);    
    folderlength2(i)=length(folderdata2{i});
end

folderData=struct(); id=1;
for i=1:length(folderdata1)
    data=folderdata1{i};
    for k=1:length(data)
        mname= data(k).mouse;
        date = data(k).date;        
        folderData(id).mouse=mname;
        folderData(id).folder=data(k).folder;
        folderData(id).date =date;
        id=id+1;
    end
end
for i=1:length(folderdata2)
    data=folderdata2{i};
    for k=1:length(data)
        mname= data(k).mouse;
        date = data(k).date;        
        folderData(id).mouse=mname;
        folderData(id).folder=data(k).folder;
        folderData(id).date =date;
        id=id+1;
    end
end

%%% linear reggresion %%%
warning off;
Regression = getData(folderData);
warning on;

% %%% save folder %%%
if(~exist([parent,outpath],'dir')), mkdir([parent,outpath]); end
cd([parent,outpath]);
save('savemat.mat','Regression');

%% plot figure %%
close all;
Data = drawData(Regression,mouse,figsaveTYPE,folderData);
drawDataAllmouse(Data,'',figsaveTYPE);
drawDataAllmouseMEtrace(Data,figsaveTYPE);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function drawDataAllmouseMEtrace(Data,figsaveTYPE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LongCtrace_mousemean = cell2mat(Data.LongCtrace_mousemean);
LongEtrace_mousemean = cell2mat(Data.LongEtrace_mousemean);
LongCtrace_total = mean(LongCtrace_mousemean,1);
LongEtrace_total = mean(LongEtrace_mousemean,1);

%%% fig S12d %%%
n=size(LongCtrace_mousemean,2); 
h=figure('Position',[10 10 1500 700]);
subplot(1,3,1);hold on;
for i=1:size(LongCtrace_mousemean,1)
    plot(1:n,LongCtrace_mousemean(i,:),'Color',[1 0 0 0.5],'LineWidth',0.25);
    plot(1:n,LongEtrace_mousemean(i,:),'Color',[0 0 0 0.5],'LineWidth',0.25);
end
plot(1:n,LongCtrace_total,'r');
plot(1:n,LongEtrace_total,'k');
xticks(0:500:n)

subplot(1,3,2); hold on;
sdata = struct(); 
sdata.xtime_msec = (1:n)'-1500;
for i=1:size(LongCtrace_mousemean,1)
    plot(1:n,LongCtrace_mousemean(i,:),'Color',[1 0 0 0.5],'LineWidth',0.25);
    eval(['sdata.mouse',num2str(i+3),'_data = transpose(LongCtrace_mousemean(i,:));']);
end
plot(1:n,LongCtrace_total,'r');
xticks(0:500:n)

%%% source data %%%
sdata.average_data = LongCtrace_total';
T = struct2table(sdata);
writetable(T, 'source fig S12d middle.csv');

subplot(1,3,3); hold on;
sdata = struct();
for i=1:size(LongCtrace_mousemean,1)
    plot(1:n,LongEtrace_mousemean(i,:),'Color',[0 0 0 0.5],'LineWidth',0.25);
    eval(['sdata.mouse',num2str(i+3),'_data = transpose(LongEtrace_mousemean(i,:));']);
end
plot(1:n,LongEtrace_total,'k');
xticks(0:500:n)

%%% source data %%%
sdata.average_data = LongCtrace_total';
T = struct2table(sdata);
writetable(T, 'source fig S12d right.csv');

set(h,'PaperPositionMode','auto');
print(h,'-r0','all mouse ME trace',figsaveTYPE);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function drawDataAllmouse(Data,figname,figsaveTYPE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mouseTag =Data.mouseTag;
session=(1:length(mouseTag))';

% model : run~ sound+choice+prior %
run_prior =Data.run_prior;
run_sound =Data.run_sound;
run_choice=Data.run_choice;
model = {run_prior, run_sound, run_choice};

%%% Kentei: linear mixed model %%%
% each val=[model parameter, running period]
estVal =zeros(3,3); Pval =zeros(3,3); Lower =zeros(3,3); Upper =zeros(3,3);
for i=1:3
    data=model{i};
    for k=1:3 
        tbl = table(data(:,k),mouseTag,session,'VariableNames',{'value','mouse','session'});
        lme1 = fitlme(tbl,'value ~ 1 + (1|mouse)+(1|session)+(mouse-1|session)'); %linear regression
        lme2 = fitlme(tbl,'value ~ 1 + (1|mouse)+(1|session)'); %linear regression
        if(lme1.ModelCriterion.BIC < lme2.ModelCriterion.BIC)
            lme=lme1;
        else
            lme=lme2;
        end
        estVal(k,i)= lme.Coefficients.Estimate;
        Pval(k,i)  = lme.Coefficients.pValue;
        Lower(k,i) = lme.Coefficients.Lower;
        Upper(k,i) = lme.Coefficients.Upper;
    end
end

h=figure('Position',[20,200,2000,500]);% model1
titlelabel ={'prior period','sound period','choice period'};
% for i=1:3
i=1;
val= estVal(:,i);
p=Pval(:,i);
low=val - Lower(:,i);
upp=Upper(:,i)-val;

%%% Fig S12e %%%
hold on;
plot([0,4],[0,0],'--k','LineWidth',0.5);
errorplot(1:3,val,upp,low,'k',0.1,1.5);
xlim([0,4]);
xticks(1:3); xticklabels({'sound','choice','prior'});
xlabel(['p: ',num2str(p(1)),' / ',num2str(p(2)),' / ',num2str(p(3))]);
ylabel('Coefficient');
title(titlelabel{i});
% end
set(h,'PaperPositionMode','auto');
print(h,'-r0',['all mouse ',figname],figsaveTYPE);

%%% source data %%%
sdata = struct();
sdata.x={'sound';'choice';'prior'};
sdata.coefficient = val;
sdata.upper = Upper(:,i);
sdata.lower = Lower(:,i);
T = struct2table(sdata);
writetable(T, 'source fig S12e.csv');


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Data = drawData(Regression,mouse,figsaveTYPE,folderData)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prior =zeros(length(Regression),3);
prior_pval=zeros(length(Regression),3);
sound =zeros(length(Regression),3);
choice=zeros(length(Regression),3);
mouseTag  =zeros(length(Regression),1);
run=zeros(length(Regression),1);
all_longcorTrace = cell(length(Regression),1);
all_longerrTrace = cell(length(Regression),1);
for i=1:length(Regression)
    run(i)=Regression(i).MEChoicePrior.runflag;
    if(~isempty(Regression(i).Coef))
        data1=Regression(i).Coef;
        data2=Regression(i).Pval;
        mouseTag(i)= find(ismember(mouse,Regression(i).mouse));
        prior(i,:) = data1.prior';
        prior_pval(i,:) = data2.prior';
        sound(i,:) = data1.sound';
        choice(i,:)= data1.choice';
        all_longcorTrace{i}=Regression(i).MEChoicePrior.allME_longcor_trace;
        all_longerrTrace{i}=Regression(i).MEChoicePrior.allME_longerr_trace;
    end
end
LongCtrace = cell2mat(all_longcorTrace);
LongEtrace = cell2mat(all_longerrTrace);

%%% Fig S12 f %%%
h=figure('Position',[20,0,1000,1500]);
label={'tone','choice','prior'};
prior_sig = zeros(size(prior,1),1);
LongCtrace_mousemean = cell(length(mouse),1);
LongEtrace_mousemean = cell(length(mouse),1);

sdata = struct();% source data
sdata.x={'sound';'choice';'prior'};
for i=1:length(mouse)
    data_id = find(mouseTag==i);
    Pval_prior= prior_pval(mouseTag==i,:);
    Coef_prior= prior(mouseTag==i,:);
    n=size(Coef_prior,1);
    
    LongCtrace_mousemean{i} = mean(LongCtrace(data_id,:),1);
    LongEtrace_mousemean{i} = mean(LongEtrace(data_id,:),1);
    
    prior_sig(data_id(Pval_prior(:,3)<0.05))=1;
    Coef_prior_on =Coef_prior(Pval_prior(:,3)<0.05,:);
    Coef_prior_off=Coef_prior(Pval_prior(:,3)>=0.05,:);
    [estVal,Pval]=sigtest(Coef_prior);    
    
    subplot(4,1,i); hold on;       
    plot([0,4],[0,0],'--k','LineWidth',0.5);
    plot(Coef_prior_on','k','LineWidth',1.5);
    plot(Coef_prior_off','Color',[.5 .5 .5],'LineWidth',1);
    plot(estVal,'r','LineWidth',1.5);
    xlim([0,4]);
    xticks(1:3); xticklabels(label);
    ylabel('Coefficient');
    title({[mouse{i},' prior period / session ', num2str(n)],...
        ['prior sig:',num2str(length(find(Pval_prior(:,3)<0.05)))],...
        num2str(Pval)});
    
    eval(['sdata.mouse',num2str(i+3),'_average=transpose(estVal);']);
    for j=1:size(Coef_prior_on,1)
        eval(['sdata.mouse',num2str(i+3),'_sig',num2str(j),'=transpose(Coef_prior_on(j,:));']);
    end
    for j=1:size(Coef_prior_off,1)
        eval(['sdata.mouse',num2str(i+3),'_nonsig',num2str(j),'=transpose(Coef_prior_off(j,:));']);
    end
end
set(h,'PaperPositionMode','auto');
print(h,'-r0','each mouse model',figsaveTYPE);

%%% source data %%%
T = struct2table(sdata);
writetable(T, 'source fig S12f.csv');

prior_sig_session=prior_sig;
for i=1:length(folderData)
    folderData(i).prior_sig = prior_sig_session(i);
end
save('folderdata.mat','folderData');
    
Data.mouseTag =mouseTag;
Data.run_prior =prior;
Data.run_sound =sound;
Data.run_choice=choice;
Data.LongCtrace_mousemean=LongCtrace_mousemean;
Data.LongEtrace_mousemean=LongEtrace_mousemean;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [estVal,Pval]=sigtest(Data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

session=(1:size(Data,1))';

%%% Kentei: linear mixed model %%%
estVal =zeros(1,3); Pval =zeros(1,3);
for k=1:3
    tbl = table(Data(:,k),session,'VariableNames',{'value','session'});
    lme = fitlme(tbl,'value ~ 1 +(1|session)'); %linear regression
    estVal(k)= mean(Data(:,k));
    Pval(k)  = lme.Coefficients.pValue;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Data = getadress(Dir,mouse,parent,dataPath)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Data=struct(); id=1;
for m=1:length(mouse)
    cd([parent,dataPath]);
    if(exist(mouse{m},'dir'))
        cd(mouse{m});
        workpath=pwd;
        list = dir(['202*_',Dir]);
        for k=1:length(list)
            cd(workpath);
            cd(list(k).name);
            Data(id).region  =Dir;
            Data(id).mouse   =mouse{m};
            Data(id).dataname=list(k).name;
            Data(id).date    =extractBefore(list(k).name,'_');
            Data(id).folder  =pwd;
            id=id+1;
        end
    end
end
% 
% parfor i=1:length(Data)
%     cd(Data(i).folder);
% %     if(~exist('MotionEnergy.mat','file'))
%         getMotionEnergy
% %     end
% end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output = getData(folderData)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output = struct(); 
parfor i=1:length(folderData)
% for i=1:length(folderData)
    cd(folderData(i).folder);
    MEChoicePrior= getTaskStructure;
    output(i).mouse = folderData(i).mouse;
    output(i).MEChoicePrior = MEChoicePrior;
    [output(i).Coef,output(i).Pval] = LenearRegression(MEChoicePrior);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output = getTaskStructure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=4000;%ms / plot point (long-tone correct/error trial)

%% Get task parameter(Bpod)
temp = dir('Bpod*.mat');
Bpod_file=temp.name;
load(temp.name);

left  = find(Chosen_side == 0);
right = find(Chosen_side == 1);
stim_length = unique(StimDuration);
Long  = find(StimDuration == stim_length(2));
Short = find(StimDuration == stim_length(1));

Choice_trial = find(Outcome == 1 | Outcome == 2);
useBlock_trial= find(TrialBlock>1);
use_trial = intersect(Choice_trial,useBlock_trial);
use_trial(end) = [];

%%% Get chosen data
output.Choice = Chosen_side(use_trial);

%%% Get tone evidence in all trials form the true tone cloud value
binary_tone=zeros(length(Tone_cloud),1);
for i = 1:length(Tone_cloud)
    temp_tone = Tone_cloud(i).matrix;
    temp1 = find(temp_tone >= 9);  %Get the data in all sound
    binary_tone(i) = length(temp1) ./ length(temp_tone);
end

%Based on the correct trial, flip the tone cloud
if mean(binary_tone(Correct_side == 1)) < 0.5 %low for right correct
    binary_tone = 1 - binary_tone;    % flipping tones %
else
end
output.ToneES = binary_tone(use_trial);

LongCorrect= intersect(Long,find(Outcome==2));
LongError  = intersect(Long,find(Outcome==1));
longcor_usetrial = ismember(use_trial,LongCorrect);
longerr_usetrial = ismember(use_trial,LongError);

%% Get prior value form RL model
temp = dir('RL_20220818*'); 
load(temp.name,'para_max','N_trial');
Prior = Dual_RL_model_block1_20220314_para_determined(Bpod_file,para_max(3,:),N_trial); 
relativePrior=Prior(:,1)./(Prior(:,1)+Prior(:,2));
output.Prior=relativePrior(1:length(use_trial),:);

%% Get wheel speed(NI daq) 
temp = dir('task_frame_tokyo_ephys_20220210*');
load(temp.name,'frame_spout','frame_sound','frame_choice','frame_end','ave_velocity');

%frame_sound_off
frame_sound_off = frame_sound;
frame_sound_off(Long) = frame_sound_off(Long) + 1000;% Add 1000 ms
frame_sound_off(Short)= frame_sound_off(Short)+ round(stim_length(1)*10)*100; % Add short stim 200 or 400 ms

frame_sound = frame_sound(use_trial);
frame_sound_off = frame_sound_off(use_trial);

%frame_choice_select
frame_choice_select = nan(length(frame_choice),1);
frame_choice_select(left) = frame_choice(left,1);
frame_choice_select(right)= frame_choice(right,2);

frame_choice_select =frame_choice_select(use_trial);
frame_spout=frame_spout(use_trial,:);

% frame end
frame_end = frame_end(use_trial);
    
if(sum(ave_velocity)==0)
    runflag=0;    
else
    runflag=1;
    %normalize ave_velocity%
    ave_velocity = rescale(ave_velocity);
    %Get wheel speed value
    speed_prior = nan(length(use_trial),1);
    speed_sound = nan(length(use_trial),1);
    speed_choice= nan(length(use_trial),1);
    speedtime_longcor_trial  = cell(length(find(longcor_usetrial)),1); j=1;
    speedtime_longerr_trial  = cell(length(find(longerr_usetrial)),1); k=1;
    for i=1:length(use_trial)
        time_trial = frame_spout(i,1):frame_end(i);
        time_prior = frame_spout(i,1):frame_sound(i);
        time_sound = frame_sound(i):frame_sound_off(i);
        time_choice= frame_sound_off(i):frame_choice_select(i);
        
        speed_prior(i) = mean(ave_velocity(time_prior));
        speed_sound(i) = mean(ave_velocity(time_sound));
        speed_choice(i)= mean(ave_velocity(time_choice));
        tmp=ave_velocity(time_trial);
        if(longcor_usetrial(i))
            speedtime_longcor_trial{j} = tmp(1:n);
            j=j+1;
        end
        if(longerr_usetrial(i))
            speedtime_longerr_trial{k} = tmp(1:n);
            k=k+1;
        end
    end
    
    %Use only during use_trial
    output.speed_prior = speed_prior;
    output.speed_sound = speed_sound;
    output.speed_choice= speed_choice;
end
output.runflag=runflag;


%% Get motion energy (DLC data) 
temp = dir('Motion*');
load(temp.name,'video_t','trialID','spout','mouth_MEraw','nose_MEraw','tongue_MEraw');
allME = mouth_MEraw+nose_MEraw+tongue_MEraw;
allME = rescale(allME);

% sychronize data between taskframe & DLC 
FirstUseTrial = use_trial(1);
FisrtSpoutAwayTime_taskframe = frame_spout(1,1);

tmp=find(spout>0.5);
tmp2=find(diff(tmp)>10);% the timing that spout go across 0.5
if(tmp(1)==1)% spout away
    spoutAwayID = tmp(tmp2); 
else
    spoutAwayID = tmp(tmp2+1); 
end
tmp=find(trialID==FirstUseTrial-1);
FisrtSpoutAwayTime_video = video_t(tmp(ismember(tmp,spoutAwayID)));
video_t = video_t-FisrtSpoutAwayTime_video;

% ajust task frame time
frame_sound = frame_sound - FisrtSpoutAwayTime_taskframe;
frame_sound_off = frame_sound_off - FisrtSpoutAwayTime_taskframe;
frame_choice_select =frame_choice_select - FisrtSpoutAwayTime_taskframe;
frame_spout = frame_spout - FisrtSpoutAwayTime_taskframe;
frame_end = frame_end - FisrtSpoutAwayTime_taskframe;

%Get motion value
ME_prior  = nan(length(use_trial),1);
ME_sound  = nan(length(use_trial),1);
ME_choice = nan(length(use_trial),1);
ME_longcor_trial = cell(length(find(longcor_usetrial)),1);
ME_longerr_trial = cell(length(find(longerr_usetrial)),1);
j=1; k=1;
for i=1:length(use_trial)
    time_trial = frame_spout(i,1):frame_end(i);
    time_prior = frame_spout(i,1):frame_sound(i);
    time_sound = frame_sound(i):frame_sound_off(i);
    time_choice= frame_sound_off(i):frame_choice_select(i);
    
    Frame_trial = ismember(video_t,time_trial);
    Frame_prior = ismember(video_t,time_prior);
    Frame_sound = ismember(video_t,time_sound);
    Frame_choice= ismember(video_t,time_choice);   

    vTime = video_t(Frame_trial);
    
    ME_prior(i) = mean(allME(Frame_prior));
    ME_sound(i) = mean(allME(Frame_sound));
    ME_choice(i)= mean(allME(Frame_choice));
    if(longcor_usetrial(i))
        t=vTime -vTime(1);
        ME_longcor_trial{j} = spCorr(n,t,allME(Frame_trial));
        j=j+1;
    end
    if(longerr_usetrial(i))
        t=vTime -vTime(1);
        ME_longerr_trial{k} = spCorr(n,t,allME(Frame_trial));
        k=k+1;
    end
end
output.ME_prior = ME_prior;
output.ME_sound = ME_sound;
output.ME_choice= ME_choice;

%% Each
allME_longcor= cell2mat(ME_longcor_trial);
allME_longerr= cell2mat(ME_longerr_trial);

output.allME_longcor_trace= mean(allME_longcor,1);
output.allME_longerr_trace= mean(allME_longerr,1);

h=figure('Position',[10 10 1000 1000],'Visible','off');
subplot(1,2,1); hold on;
plot(1:n,allME_longcor(1,:),'k');
title('single correct trial')
xticks(0:500:n)

subplot(1,2,2); hold on;
p1=errorplot(1:n,mean(allME_longcor,1),std(allME_longcor,1),std(allME_longcor,1),'r',.5,1);
p2=errorplot(1:n,mean(allME_longerr,1),std(allME_longerr,1),std(allME_longerr,1),'k',.5,1);
legend([p1,p2],{'cor','err'});
title('trials')
xticks(0:500:n)

set(h,'PaperPositionMode','auto');
print(h,'-r0','move','-dpng');
print(h,'-r0','move','-dsvg');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  output= spCorr(n,t,v)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F = griddedInterpolant(t,v);
t_correct = 0:1:t(end);
v_correct = F(t_correct);
output= v_correct(1:n);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [Coef,Pval]= LenearRegression(MEChoicePrior)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Choice = MEChoicePrior.Choice;
ToneES = MEChoicePrior.ToneES;
Prior  = MEChoicePrior.Prior;
ME_prior = MEChoicePrior.ME_prior;
ME_sound = MEChoicePrior.ME_sound;
ME_choice= MEChoicePrior.ME_choice;

% linear regression %
[Coef.prior, Pval.prior] = LRmodel(ME_prior,ToneES,Choice,Prior);
[Coef.sound, Pval.sound] = LRmodel(ME_sound,ToneES,Choice,Prior);
[Coef.choice,Pval.choice]= LRmodel(ME_choice,ToneES,Choice,Prior);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [Coef,Pval]= LRmodel(Run,ToneES,Choice,Prior)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (regression mode2) Run ~ ToneES + Choice + Prior

% data X (predictor) %
X = [ToneES,Choice,Prior];

% data y (explain) %
y = Run;

tmp=[X,y];
test = tmp;

% table data
tbl= table(test(:,1),test(:,2),test(:,3),test(:,4),...
    'VariableNames',{'sound','choice','prior','run'});

% linear regression %
mdl = fitlm(tbl,'run ~ sound + choice + prior ');
Coef = mdl.Coefficients.Estimate(2:end);
Pval = mdl.Coefficients.pValue(2:end);
end