function FigureS12gh_WheelSpeedvsPrior
% relavant to the REVIER#1 major 1-3 comment
%
% The influence of task-irrelevant movement such as running speed to
% cortical activity
% The correlation with prior value and biased choice in subsequent trials 

parent='G:/Ishizu_data';
workpath='/Revise_ishizu';
outpath='/Revise_ishizu/output/WheelSpeedvsPrior2';

% hulistic parameter %
figsaveTYPE='-dpng';
%--------------------%

cd([parent,workpath]);

%%% collect the data %%%
% folder setting %
Dir1 = {'auditory','fof','mpfc'};
mouse= {'a04','a08','i20','i24','i34','i35','i43','i46'};
dataPath='/Tokyo_ephys_ishizu';

disp('data collecting...');
folderdata1=cell(length(Dir1),1);
folderlength1=zeros(length(Dir1),1);
for i=1:length(Dir1)
    folderdata1{i} = getadress(Dir1{i},mouse,parent,dataPath);
    folderlength1(i)=length(folderdata1{i});
end

folderData=struct(); id=1;
for i=1:length(folderdata1)
    data=folderdata1{i};
    for k=1:length(data)
        folderData(id).mouse =data(k).mouse;
        folderData(id).folder=data(k).folder;
        folderData(id).date=cell2mat(extract(data(k).date,digitsPattern)');
        id=id+1;
    end
end

%%% linear reggresion %%%
warning off;
% Regression = getData(folderData);
warning on;

% %%% save folder %%%
if(~exist([parent,outpath],'dir')), mkdir([parent,outpath]); end
cd([parent,outpath]);
% save('savemat.mat','Regression');

%% plot figure %%
close all;
Data = drawData(Regression,mouse,figsaveTYPE,folderData);
drawDataAllmouse(Data,'',figsaveTYPE);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function drawDataAllmouse(Data,figname,figsaveTYPE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mouseTag =Data.mouseTag;
session=(1:length(mouseTag))';

% model: run~ sound+choice+prior %
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

%%% Fig S12g %%%
h=figure('Position',[20,200,2000,500]);
titlelabel ={'prior period','sound period','choice period'};
% for i=1:3
i=1;
val= estVal(:,i);
p=Pval(:,i);
low =val - Lower(:,i);
upp =Upper(:,i)-val;

%     subplot(1,3,i); hold on;
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
writetable(T, 'source fig S12g.csv');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Data = drawData(Regression,mouse,figsaveTYPE,folderData)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run_prior =zeros(length(Regression),3);
run_prior_pval =zeros(length(Regression),3);
run_sound =zeros(length(Regression),3);
run_choice=zeros(length(Regression),3);
mouseTag   =zeros(length(Regression),1);
run=zeros(length(Regression),1);
for i=1:length(Regression)
    run(i)=Regression(i).WheelChoicePrior.runflag;
    if(~isempty(Regression(i).Coef))
        data1=Regression(i).Coef;
        data2=Regression(i).Pval;
        mouseTag(i)   = find(ismember(mouse,Regression(i).mouse));
        run_prior(i,:) = data1.run_prior';
        run_prior_pval(i,:) = data2.run_prior';
        run_sound(i,:) = data1.run_sound';
        run_choice(i,:)= data1.run_choice';
    end
end

%%% except for no run session %%%
mouseTag(run==0,:)=[];
run_prior(run==0,:)=[];
run_prior_pval(run==0,:)=[];
run_sound(run==0,:)=[];
run_choice(run==0,:)=[];

%%% fig S12h %%%
h=figure('Position',[20,0,1000,1500]);
label2={'tone','choice','prior'};
prior_sig = zeros(size(run_prior,1),1);

sdata = struct();% source data
sdata.x={'sound';'choice';'prior'};
for i=1:length(mouse)
    data_id = find(mouseTag==i);
    Pval_prior= run_prior_pval(mouseTag==i,:);
    Coef_prior= run_prior(mouseTag==i,:);
    Coef_sound= run_sound(mouseTag==i,:);
    Coef_choice= run_choice(mouseTag==i,:);
    n=size(Coef_prior,1);
    
    prior_sig(data_id(Pval_prior(:,3)<0.05))=1;
    Coef_prior_on =Coef_prior(Pval_prior(:,3)<0.05,:);
    Coef_prior_off=Coef_prior(Pval_prior(:,3)>=0.05,:);
    [estVal,Pval]=sigtest(Coef_prior);            
    subplot(8,3,1+3*(i-1)); hold on;
    plot([0,4],[0,0],'--k','LineWidth',0.5);
    plot(Coef_prior_on','k','LineWidth',1.5);
    plot(Coef_prior_off','Color',[.5 .5 .5],'LineWidth',1);
    plot(estVal,'r','LineWidth',1.5);
    xlim([0,4]);
    xticks(1:3); xticklabels(label2);
    ylabel('Coefficient');
    title({[mouse{i},' prior period / session ', num2str(n)],...
        ['prior sig:',num2str(length(find(Pval_prior(:,3)<0.05)))],num2str(Pval)});
    
    % source data %
    eval(['sdata.mouse',num2str(i),'_average=transpose(estVal);']);
    for j=1:size(Coef_prior_on,1)
        eval(['sdata.mouse',num2str(i),'_sig',num2str(j),'=transpose(Coef_prior_on(j,:));']);
    end
    for j=1:size(Coef_prior_off,1)
        eval(['sdata.mouse',num2str(i),'_nonsig',num2str(j),'=transpose(Coef_prior_off(j,:));']);
    end
    % ----------- %
    
    [estVal,Pval]=sigtest(Coef_sound);
    subplot(8,3,2+3*(i-1)); hold on;
    plot([0,4],[0,0],'--k','LineWidth',0.5);
    plot(Coef_sound','Color',[.5 .5 .5],'LineWidth',1);
    plot(estVal,'r','LineWidth',1.5);
    xlim([0,4]);
    xticks(1:3); xticklabels(label2);
    ylabel('Coefficient');
    title({[mouse{i},' sound period'],num2str(Pval)});
    
    [estVal,Pval]=sigtest(Coef_choice);
    subplot(8,3,3+3*(i-1)); hold on;
    plot([0,4],[0,0],'--k','LineWidth',0.5)
    plot(Coef_choice','Color',[.5 .5 .5],'LineWidth',1);
    plot(estVal,'r','LineWidth',1.5);
    xlim([0,4]);
    xticks(1:3); xticklabels(label2);
    ylabel('Coefficient');
    title({[mouse{i},' choice period'],num2str(Pval)});
end
set(h,'PaperPositionMode','auto');
print(h,'-r0','each mouse model',figsaveTYPE);

%%% source data %%%
T = struct2table(sdata);
writetable(T, 'source fig S12h.csv');


prior_sig_session=nan(length(folderData),1);
prior_sig_session(run==1)=prior_sig;
for i=1:length(folderData)
    folderData(i).prior_sig = prior_sig_session(i);
end
save('folderdata.mat','folderData');
    
Data.mouseTag =mouseTag;
Data.run_prior =run_prior;
Data.run_sound =run_sound;
Data.run_choice=run_choice;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [estVal,Pval]=sigtest(Data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

session=(1:size(Data,1))';

%%% Kentei: linear mixed model %%%
% each val=[model parameter, running period]
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
    cd([parent,dataPath,'/',Dir]);
    if(exist(mouse{m},'dir'))
        cd(mouse{m});
        workpath=pwd;
        list = dir('202*');
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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Data = getadressPPC(Dir,AudF,mouse,parent,dataPath)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Data=struct(); id=1;
for m=1:length(mouse)
    date_lap={};s=1;
    for i=1:length(AudF)
    if(strcmp(AudF(i).mouse,mouse{m}))
        date_lap{s}=AudF(i).date;
        s=s+1;
    end
    end
    cd([parent,dataPath,'/',Dir]);
    if(exist(mouse{m},'dir'))
        cd(mouse{m});
        workpath=pwd;
        list = dir('202*');
        for k=1:length(list)
            date = extractBefore(list(k).name,'_');
            out  = 0;
            if(~isempty(date_lap))
            for t=1:length(date_lap)
                if(strcmp(date_lap{t},date))
                    out=1;
                end
            end
            end
            if(out==0)
                cd(workpath);
                cd(list(k).name);
                Data(id).region  =Dir;
                Data(id).mouse   =mouse{m};
                Data(id).dataname=list(k).name;
                Data(id).date    =date;
                Data(id).folder  =pwd;
                id=id+1;
            end
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output = getData(folderData)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output = struct(); 
parfor i=1:length(folderData)
% for i=1:length(folderData)
    cd(folderData(i).folder);
    f=dir('*_task');
    cd(f.name);
    WheelChoicePrior = getTaskStructure;
    output(i).mouse = folderData(i).mouse;
    output(i).WheelChoicePrior = WheelChoicePrior;
    [output(i).Coef,output(i).Pval] = LenearRegression(WheelChoicePrior,0);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output = getTaskStructure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get task parameter(Bpod)
temp = dir('Bpod*');
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

%% Get prior value form RL model
temp = dir('RL_20220818*'); 
load(temp.name,'para_max','N_trial');
Prior = Dual_RL_model_block1_20220314_para_determined(Bpod_file,para_max(3,:),N_trial); 
relativePrior=Prior(:,1)./(Prior(:,1)+Prior(:,2));
output.Prior=relativePrior(1:length(use_trial),:);

%% Get wheel speed(NI daq) 
temp = dir('task_frame_tokyo_ephys_20220210*');
load(temp.name,'frame_spout','frame_sound','frame_choice','ave_velocity');

if(sum(ave_velocity)==0)
    runflag=0;
else
    runflag=1;
    
    %normalize ave_velocity%
    ave_velocity = rescale(ave_velocity);
    
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
    
    %Get wheel speed value
    speed_before_sound = nan(length(use_trial),1);
    speed_for_sound    = nan(length(use_trial),1);
    speed_before_choice= nan(length(use_trial),1);
    for i=1:length(use_trial)
        time_before_sound = frame_spout(i,1):frame_sound(i);
        time_for_sound = frame_sound(i):frame_sound_off(i);
        time_after_sound_before_choice = frame_sound_off(i):frame_choice_select(i);
        
        speed_before_sound(i) =mean(ave_velocity(time_before_sound));
        speed_for_sound(i)    =mean(ave_velocity(time_for_sound));
        speed_before_choice(i)=mean(ave_velocity(time_after_sound_before_choice));
    end
    
    %Use only during use_trial
    output.speed_before_sound = speed_before_sound;
    output.speed_for_sound = speed_for_sound;
    output.speed_before_choice = speed_before_choice;
end
output.runflag=runflag;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [Coef,Pval]= LenearRegression(WheelChoicePrior,transflag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

flag=WheelChoicePrior.runflag;
if(flag==1)
    Choice = WheelChoicePrior.Choice;
    ToneES = WheelChoicePrior.ToneES;
    Prior  = WheelChoicePrior.Prior;
    Run_prior = WheelChoicePrior.speed_before_sound;
    Run_sound = WheelChoicePrior.speed_for_sound;
    Run_choice= WheelChoicePrior.speed_before_choice;
        
    % linear regression %
    [Coef.run_prior, Pval.run_prior] = LRmodel(Run_prior,ToneES,Choice,Prior,transflag);
    [Coef.run_sound, Pval.run_sound] = LRmodel(Run_sound,ToneES,Choice,Prior,transflag);
    [Coef.run_choice,Pval.run_choice]= LRmodel(Run_choice,ToneES,Choice,Prior,transflag);
else
    Coef=[];
    Pval=[];
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [Coef,Pval]= LRmodel(Run,ToneES,Choice,Prior,transflag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (regression mode2) Run ~ ToneES + Choice + Prior

% data X (predictor) %
X = [ToneES,Choice,Prior];

% data y (explain) %
y = Run;

test =[X,y];

% table data
tbl= table(test(:,1),test(:,2),test(:,3),test(:,4),...
    'VariableNames',{'sound','choice','prior','run'});

% linear regression %
mdl = fitlm(tbl,'run ~ sound + choice + prior ');
Coef = mdl.Coefficients.Estimate(2:end);
Pval = mdl.Coefficients.pValue(2:end);
end