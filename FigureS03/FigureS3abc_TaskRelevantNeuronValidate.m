function FigureS3abc_TaskRelevantNeuronValidate
parent='G:/Ishizu_data';
outpath='/Revise_ishizu/output/TaskRelevantNeuronValidate';

% hulistic parameter %
%Differentially analyze between long sound and short sound
TW.time_window_base = 200; %ms
TW.time_window = 100; %ms
TW.time_long_pre  = 1500;
TW.time_long_pre2 = 500;
TW.time_long_post2= 2500;
TW.time_short_pre = 1500;
TW.time_short_pre2= 500;
TW.time_short_post2= 2500;
% 
p_threshold = 1.0e-10;  %10e-10
figsaveTYPE='-dpng';
%--------------------%

%%% collect the data of firing rate %%%
% folder setting %
region = {'auditory','fof','mpfc'};
mouse= {'a04','a08','i20','i24','i34','i35','i43','i46'};
dataPath='/Tokyo_ephys_ishizu';

% save folder %%%
if(~exist([parent,outpath],'dir')), mkdir([parent,outpath]); end

disp('data collecting...');
SpikeData_AC  = getData(p_threshold,region{1},mouse,TW,parent,dataPath);
cd([parent,outpath]);
save('savematAC.mat','SpikeData_AC');
save('TW.mat','TW');

SpikeData_FOF = getData(p_threshold,region{2},mouse,TW,parent,dataPath);
cd([parent,outpath]);
save('savematFOF.mat','SpikeData_FOF');

SpikeData_MPFC= getData(p_threshold,region{3},mouse,TW,parent,dataPath);
cd([parent,outpath]);
save('savematMPFC.mat','SpikeData_MPFC');

%% plot figure %%
cd([parent,outpath]);
outputAC  = getDistribution(SpikeData_AC,mouse);
outputFOF = getDistribution(SpikeData_FOF,mouse);
outputMPFC= getDistribution(SpikeData_MPFC,mouse);

%%% FigS3a %%%
[~,ac_taskPdist , ~,ac_notaskPdist]  =drawDistribution(outputAC,  p_threshold,'AC',figsaveTYPE);
[~,fof_taskPdist, ~,fof_notaskPdist] =drawDistribution(outputFOF, p_threshold,'FOF',figsaveTYPE);
[~,mpfc_taskPdist,~,mpfc_notaskPdist]=drawDistribution(outputMPFC,p_threshold,'mPFC',figsaveTYPE);

% ranksum(ac_taskPdist,ac_notaskPdist)
% ranksum(fof_taskPdist,fof_notaskPdist)
% ranksum(mpfc_taskPdist,mpfc_notaskPdist)

%% compare P distributuion across regions

x = [ac_notaskPdist; ac_taskPdist; fof_notaskPdist; fof_taskPdist; mpfc_notaskPdist; mpfc_taskPdist];


g1 = repmat({'ACnotask'}, length(ac_notaskPdist),1);
g2 = repmat({'ACtask'},   length(ac_taskPdist),1);
g3 = repmat({'FOFnotask'},length(fof_notaskPdist),1);
g4 = repmat({'FOFtask'},  length(fof_taskPdist),1);
g5 = repmat({'mPFCnotask'},length(mpfc_notaskPdist),1);
g6 = repmat({'mPFCtask'}, length(mpfc_taskPdist),1);
g = [g1; g2; g3; g4; g5; g6];

y1 = ones(length(ac_notaskPdist),1);
y2 = 2*ones(length(ac_taskPdist),1);
y3 = 3*ones(length(fof_notaskPdist),1);
y4 = 4*ones(length(fof_taskPdist),1);
y5 = 5*ones(length(mpfc_notaskPdist),1);
y6 = 6*ones(length(mpfc_taskPdist),1);
y = [y1; y2; y3; y4; y5; y6];

%%% FigS3b %%%
tag = {'ACnotask','ACtask ','FOFnotask','FOFtask','mPFCnotask','mPFCtask'};
drawCompareDistribution(x, g, y, tag, 'rawdist_',p_threshold,figsaveTYPE);

%% scatter plot(wilcoxon vs distribution)

%%% FigS3c %%%
th_dist=2;
drawScatter(outputAC,'AC task scatter',figsaveTYPE,th_dist);
drawScatter(outputFOF,'FOF task scatter',figsaveTYPE,th_dist);
drawScatter(outputMPFC,'MPFC task scatter',figsaveTYPE,th_dist);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function drawScatter(output,name,figsaveTYPE,th_dist)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_pdist= -log10(output.taskPdist);
t_pval = -log10(output.taskP);
t_mouse=  output.task_mouse;
t_session=output.task_session;

tbl = table(t_pdist,t_pval,t_mouse,t_session,'VariableNames',{'value','pval','mouse','session'});
lme = fitlme(tbl,'value ~ 1 + pval+(1|mouse)+(1|session)'); %linear regression
estVal= lme.Coefficients.Estimate;
Pval  = lme.Coefficients.pValue;
% Rcorr=corr(t_pdist,t_pval,'Type','Spearman');

%%% Fig S3c %%%
h=figure('Position',[100 100 1500 500]);
hold on
scatter(t_pval,t_pdist,5,'filled','b')
plot(0:80,estVal(2)*(0:80)+estVal(1),'r');
ylabel('P value (-log10) from firing rate distribution');
xlabel('P value (-log10) from wilcoxon criteria');
xlim([0 80])
% title({['a: (est)',num2str(estVal(2)),' / (p)',num2str(Pval(2))],...
%     ['b: (est)',num2str(estVal(1)),' / (p)',num2str(Pval(1))],...
%     ['R: ',num2str(Rcorr)]});
axis square

set(h,'PaperPositionMode','auto');
print(h,'-r0',name,figsaveTYPE);

% source data %
sdata = struct();
sdata.x_scatter =t_pval;
sdata.y_scatter =t_pdist;
T = struct2table(sdata);
writetable(T, ['source fig S3c ',name(1:4),'.csv']);

sdata = struct();
sdata.x_lmfit =(0:80)';
sdata.y_lmfit =(estVal(2)*(0:80)+estVal(1))';
T = struct2table(sdata);
writetable(T, ['source fig S3c regression',name(1:4),'.csv']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function drawCompareDistribution(x,g,y,tag,name,p_threshold,figsaveTYPE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% kruskal-wallis test %%%
[~,~,stats] = kruskalwallis(x,g,'off');
c = multcompare(stats,'Display','off'); pvalue = c(:,6);
id1 = find(c(:,6) < 0.05); id2 = setdiff(1:size(c,1),id1);
if(~isempty(id1))
    label1 = cell(1,length(id1));
    for i=1:length(id1),    tmp = id1(i);
        if(i==length(id1)), label1{i} = strcat(tag{c(tmp,1)},'-',tag{c(tmp,2)});
        else,               label1{i} = strcat(tag{c(tmp,1)},'-',tag{c(tmp,2)},' / ');
        end
    end
end
if(~isempty(id2))
    label2 = cell(1,length(id2));
    for i=1:length(id2),    tmp = id2(i);
        if(i==length(id2)), label2{i} = strcat(tag{c(tmp,1)},'-',tag{c(tmp,2)});
        else,               label2{i} = strcat(tag{c(tmp,1)},'-',tag{c(tmp,2)},' / ');
        end
    end
end
if(isempty(id1)),    XLabel = {['significant: ','NaN']};
elseif(isempty(id2)),XLabel = {['significant: ','ALL'],num2str(pvalue')};
else, XLabel = {['significant: ',cell2mat(label1)], num2str(pvalue(id1)'),['non-significant: ',cell2mat(label2)], num2str(pvalue(id2)')};
end

%%% Fig S3b %%%
h=figure('Position',[100 100 1500 500]);
swarmchart(y,x,2,'filled');
ylabel('P value (-log10) from firing rate distribution');
title(XLabel);
xticks(1:length(tag));
xticklabels(tag);
pval=-log10(p_threshold);

% source data %
sdata = struct();
sdata.neuralactivity =x;
sdata.labelid =y;
T = struct2table(sdata);
writetable(T, 'source fig S3b.csv');

sdata = struct();
sdata.label =tag;
T = struct2table(sdata);
writetable(T, 'source fig S3b label.csv');

set(h,'PaperPositionMode','auto');
print(h,'-r0',[name,'CpmpareAll_pThreshold_e-',num2str(pval)],figsaveTYPE);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [allPdist,taskPdist,taskP,notaskPdist]=...
    drawDistribution(output,p_threshold,figname,figsaveTYPE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pval=-log10(p_threshold);
Nbin=50;

allPdist =-log10(output.allPdist(:));
taskPdist=-log10(output.taskPdist(:));
taskP = -log10(output.taskP(:));
notaskPdist=-log10(output.notaskPdist(:));

%%% Fig S3a %%%
h=figure('Position',[100 100 500 500]);
h1=histogram(-log10(output.allP(:)),0:5:80);
% x=h1.BinEdges+h1.BinEdges/2;
x=h1.BinEdges;
x(end)=[];
counts1=h1.Values;
h2=histogram(-log10(output.taskP(:)),h1.BinEdges);
counts=h2.Values;
x2=[x(1:2),10,x(3:end)];
counts2=[counts(1:2),0,counts(3:end)];
a1=area(x,counts1,'FaceColor','w'); hold on;
a2=area(x2,counts2,'FaceColor','r');
legend([a1,a2],{'all neuron','task relevant'});
xlim([-1 80])
xlabel('p value (-log10)');
ylabel('# of neuron');
title('wilcoxon criteria');

% source data %
sdata = struct();
sdata.x_allNeuron =x';
sdata.y_allNeuron =counts1';
T = struct2table(sdata);
writetable(T, ['source fig S3a all neuron ',figname,'.csv']);

sdata = struct();
sdata.x_taskNeuron =x2';
sdata.y_taskNeuron =counts2';
T = struct2table(sdata);
writetable(T, ['source fig S3a task neuron ',figname,'.csv']);

set(h,'PaperPositionMode','auto');
print(h,'-r0',[figname,'pThreshold_e-',num2str(pval)],figsaveTYPE);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output = getDistribution(SpikeData,mouse)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

base_frate = cell(length(SpikeData),1);
task_frate = cell(length(SpikeData),1);
notask_frate = cell(length(SpikeData),1);
all_pval = cell(length(SpikeData),1);
task_pval= cell(length(SpikeData),1);
notask_pval= cell(length(SpikeData),1);
all_pdist= cell(length(SpikeData),1);
task_pdist= cell(length(SpikeData),1);
notask_pdist= cell(length(SpikeData),1);
task_mouse= cell(length(SpikeData),1);
notask_mouse= cell(length(SpikeData),1);
task_session= cell(length(SpikeData),1);
notask_session= cell(length(SpikeData),1);
for i=1:length(SpikeData)
    data =SpikeData(i).data;
    
    tmp_sig = nan(length(data),1);
    tmp_baseFR = nan(length(data),1);
    tmp_taskFR = nan(length(data),1);
    tmp_allPval  = nan(length(data),1);
    tmp_allPdist = nan(length(data),1);
    for j=1:length(data)
        tmp_sig(j) = data(j).task_sig;
        tmp_baseFR(j)  = data(j).baseActivity;
        tmp_taskFR(j)  = data(j).taskActivity;
        tmp_allPval(j) = data(j).min_pval;
        tmp_allPdist(j)= data(j).min_pdist;
    end
    base_frate{i} = tmp_baseFR;
    task_frate{i} = tmp_taskFR(tmp_sig==1);
    notask_frate{i} = tmp_taskFR(tmp_sig==0);
    all_pval{i}   = tmp_allPval;
    task_pval{i}  = tmp_allPval(tmp_sig==1);
    notask_pval{i}= tmp_allPval(tmp_sig==0);
    all_pdist{i}  = tmp_allPdist;
    task_pdist{i} = tmp_allPdist(tmp_sig==1);
    notask_pdist{i}=tmp_allPdist(tmp_sig==0);
    
    data_session=i;
    data_mouse= find(ismember(mouse,SpikeData(i).mouse));
    task_mouse{i}= data_mouse*ones(length(find(tmp_sig==1)),1);
    task_session{i}= data_session*ones(length(find(tmp_sig==1)),1);
    notask_mouse{i}= data_mouse*ones(length(find(tmp_sig==0)),1);
    notask_session{i}= data_session*ones(length(find(tmp_sig==0)),1);
end
output.baseSpike = cell2mat(base_frate);
output.allP = cell2mat(all_pval);
output.allPdist = cell2mat(all_pdist);
output.task_mouse = cell2mat(task_mouse);
output.task_session=cell2mat(task_session);
output.notask_mouse = cell2mat(notask_mouse);
output.notask_session=cell2mat(notask_session);

% x=cellfun(@isempty,task_frate);
% task_frate(x==1)=[];
output.taskSpike = cell2mat(task_frate);

x=cellfun(@isempty,notask_frate);
notask_frate(x==1)=[];
output.notaskSpike = cell2mat(notask_frate);

% x=cellfun(@isempty,task_pval);
% task_pval(x==1)=[];
output.taskP= cell2mat(task_pval);

x=cellfun(@isempty,notask_pval);
notask_pval(x==1)=[];
output.notaskPval= cell2mat(notask_pval);

% x=cellfun(@isempty,task_pdist);
% task_pdist(x==1)=[];
output.taskPdist= cell2mat(task_pdist);

x=cellfun(@isempty,notask_pdist);
notask_pdist(x==1)=[];
output.notaskPdist= cell2mat(notask_pdist);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SpikeData = getData(p_threshold,region,mouse,TW,parent,dataPath)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(region);
SpikeData=struct(); id=1;
for m=1:length(mouse)
    cd([parent,dataPath,'/',region]);
    if(exist(mouse{m},'dir'))
        cd(mouse{m});
        workpath=pwd;
        list = dir('202*');
        for k=1:length(list)
            cd(workpath);
            cd(list(k).name);
            SpikeData(id).mouse   =mouse{m};
            SpikeData(id).dataname=list(k).name;
            SpikeData(id).folder  =pwd;
            id=id+1;
        end
    end
end

for i=1:length(SpikeData)
% parfor i=1:length(SpikeData)
    cd(SpikeData(i).folder);
    f=dir('*_task');
    cd(f.name);
    frame = getTaskStructure(TW);
    SpikeData(i).data = getSpikePdist(p_threshold,frame,TW);
end
fprintf('\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output = getTaskStructure(TW)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Start of task
%Sound on
%Sound off
%Before choice
%After choice (0sec)
%After choice (1sec)
%After choice (2sec)

%%% Get task parameter(Bpod) %%%
temp = dir('Bpod*');
load(temp.name);

Choice_trial = find(Outcome == 1 | Outcome == 2);
useBlock_trial= find(TrialBlock>1);
use_trial = sort(intersect(Choice_trial,useBlock_trial));
% use_trial(end) = [];

left  = find(Chosen_side == 0);
right = find(Chosen_side == 1);
stim_length = unique(StimDuration);
Long  = find(StimDuration == stim_length(2));
Short = find(StimDuration == stim_length(1));

use_long  = intersect(use_trial,Long);
use_short = intersect(use_trial,Short);

%%% Get frame parameter(NI daq) %%%
temp = dir('task_frame_tokyo_ephys_20220210*');
load(temp.name,'frame_spout','frame_sound','frame_choice');

time_window = TW.time_window;
time_window_base = TW.time_window_base;
time_long_pre  = TW.time_long_pre;
time_long_pre2 = TW.time_long_pre2;
time_long_post2=TW.time_long_post2;
time_short_pre  =TW.time_short_pre;
time_short_pre2 =TW.time_short_pre2;
time_short_post2=TW.time_short_post2;

%sound_off
time_long_post= time_long_pre+1000;
time_short_post=time_short_pre+round(stim_length(1)*10)*100;% Add short stim 200 or 400 ms

%frame_choice_select
frame_choice_select = nan(length(frame_choice),1);
frame_choice_select(left) = frame_choice(left,1);
frame_choice_select(right)= frame_choice(right,2);
frame_choice_select_long = frame_choice_select(use_long,:);
frame_choice_select_short= frame_choice_select(use_short,:);

%%% long / short anlysis window %%%
time_long  = round((time_long_pre + time_long_post) ./ time_window);
time_long2 = round((time_long_pre2+ time_long_post2)./ time_window);
time_short  = round((time_short_pre + time_short_post) ./ time_window);
time_short2 = round((time_short_pre2+ time_short_post2)./ time_window);

frame_sound_long = frame_sound(use_long,:);
frame_sound_short= frame_sound(use_short,:);

frame_spout_long = frame_spout(use_long,:);
base_for_spout_on_long  = frame_spout_long(:,1)-time_window_base; %Before moving spout
base_for_spout_off_long = frame_spout_long(:,3)-time_window_base; %Before moving spout
frame_spout_short = frame_spout(use_short,:);
base_for_spout_on_short  = frame_spout_short(:,1)-time_window_base; %Before moving spout
base_for_spout_off_short = frame_spout_short(:,3)-time_window_base; %Before moving spout

frame_sound_use_on_long = zeros(length(frame_sound_long),length(time_long));
for i = 1:time_long
    frame_sound_use_on_long(:,i) = frame_sound_long - time_long_pre + (i-1)*time_window;
end
frame_sound_use_on_short= zeros(length(frame_sound_short),length(time_short));
for i = 1:time_short
    frame_sound_use_on_short(:,i) = frame_sound_short - time_short_pre + (i-1)*time_window;
end

frame_sound_use2_on_long = zeros(length(frame_choice_select_long),length(time_long));
for i = 1:time_long2
    frame_sound_use2_on_long(:,i) = frame_choice_select_long - time_long_pre2 + (i-1)*time_window;
end
frame_sound_use2_on_short= zeros(length(frame_choice_select_short),length(time_short));
for i = 1:time_short2
    frame_sound_use2_on_short(:,i) = frame_choice_select_short - time_short_pre2 + (i-1)*time_window;
end

%Use only during use_trial
output.use_long  = use_long;  output.use_short  = use_short;
output.time_long = time_long; output.time_long2 = time_long2;
output.time_short= time_short;output.time_short2= time_short2;

output.base_for_spout_on_long  = base_for_spout_on_long;
output.base_for_spout_off_long = base_for_spout_off_long;
output.base_for_spout_on_short = base_for_spout_on_short;
output.base_for_spout_off_short= base_for_spout_off_short;

output.frame_sound_long = frame_sound_long;
output.frame_sound_short= frame_sound_short;
output.frame_sound_use_on_long =frame_sound_use_on_long;
output.frame_sound_use2_on_long=frame_sound_use2_on_long;
output.frame_sound_use_on_short =frame_sound_use_on_short;
output.frame_sound_use2_on_short=frame_sound_use2_on_short;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function neuron = getSpikePdist(p_threshold,frame,TW)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

temp = dir('sig2_task_neurons_2022*');
load(temp.name,'p_task_long','p_task2_long','p_task_short','p_task2_short');

%%% get spike depth data %%%
temp = dir('depth_spike_20220517*');
if length(temp) == 1
    load(temp.name,'spike_depth','def_depth');
    depth_neuron = find(spike_depth <= def_depth);
else
    depth_neuron = 1:size(p_task_long,1); %Use all the neurons
end

%%% get significant %%%
p_task_long  = p_task_long(depth_neuron,:);
p_task2_long = p_task2_long(depth_neuron,:);
p_task_short = p_task_short(depth_neuron,:);
p_task2_short= p_task2_short(depth_neuron,:);
p = [p_task_long,p_task2_long,p_task_short,p_task2_short];
sig_neuron = find(sum(p < p_threshold,2)>0); % sig & depth
Pmin_task_long  = min(p_task_long,[],2);
Pmin_task2_long = min(p_task2_long,[],2);
Pmin_task_short = min(p_task_short,[],2);
Pmin_task2_short= min(p_task2_short,[],2);
Pmin = [Pmin_task_long, Pmin_task2_long, Pmin_task_short,Pmin_task2_short];
[~,pickup]=min(Pmin,[],2);

%%% spike data %%%
spike_dir = dir('spike_ch*');
cd(spike_dir.name);

neuron = struct();
for file_count = 1:length(depth_neuron)
    if(ismember(file_count,sig_neuron))
        neuron(file_count).task_sig = 1;
    else
        neuron(file_count).task_sig = 0;
    end
    temp_file = sprintf('task_spike_stripe20210520_%d',depth_neuron(file_count));
    data = load(temp_file); %spike_mark
    spike_mark = data.spike_mark;
    
    %%% get frame data %%%
    flag=pickup(file_count);
    switch flag
        case 1 % long tone
            use_trial=frame.use_long;
            base_for_spout = frame.base_for_spout_on_long;
            frame_sound_use_on =frame.frame_sound_use_on_long;
            p_task = p_task_long(file_count,:);
        case 2 % long choice
            use_trial=frame.use_long;
            base_for_spout= frame.base_for_spout_off_long;
            frame_sound_use_on=frame.frame_sound_use2_on_long;
            p_task= p_task2_long(file_count,:);
        case 3 % short tone
            use_trial=frame.use_short;
            base_for_spout = frame.base_for_spout_on_short;
            frame_sound_use_on =frame.frame_sound_use_on_short;
            p_task = p_task_short(file_count,:);
        case 4 % short choice
            use_trial=frame.use_short;
            base_for_spout = frame.base_for_spout_on_short;
            frame_sound_use_on =frame.frame_sound_use2_on_short;
            p_task = p_task2_short(file_count,:);
    end
    [neuron(file_count).baseActivity, neuron(file_count).taskActivity,...
        neuron(file_count).min_pval, neuron(file_count).min_pdist] = ...
        calcData(spike_mark,TW,use_trial,base_for_spout,frame_sound_use_on,p_task);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [baseActivity, taskActivity, min_pval, min_pdist] = ...
    calcData(spike_mark,TW,use_trial,base_for_spout,frame_sound_use_on,p_task)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[min_pval,id] = min(p_task);

spike_base = nan(length(use_trial),1);
spike_frame= nan(size(frame_sound_use_on,1),1);
for i = 1:length(use_trial) %trial
    base_frame = base_for_spout(i) : base_for_spout(i) + TW.time_window_base-1;
    temp_frame = frame_sound_use_on(i,id) : frame_sound_use_on(i,id) + TW.time_window-1;
    spike_base(i) = mean(spike_mark(base_frame));
    spike_frame(i)= mean(spike_mark(temp_frame));
end

%%% pvalue from spike distribution
spike=[mean(spike_frame);spike_base];
baseActivity = mean(spike_base);
taskActivity = mean(spike_frame);
[~,I]=sort(spike);
min_pdist =(length(spike)-find(I==1)+1)/length(spike);
end