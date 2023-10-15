function FigureS3d_NeuralSequenceTest_shuffleTimeCourse

parent='G:/Ishizu_data';
dataPath='/Revise_ishizu/output/SpikeTrace_longtone';
outPath='/Revise_ishizu/output/NeuralSequenceTest_shuffleTimeCourse';

% hulistic parameter %
shuffNum=100;
p_threshold = 1.0e-10; 
figsaveTYPE='-dsvg';
%--------------------%

%% collect the data of firing rate %%%
% folder setting %
region = {'auditory','fof','mpfc'};
folders= {'auc_ishizu','fof_ishizu','mpfc_ishizu'};
% save folder setting %
savefolder =[parent,outPath];
if(~exist(savefolder,'dir')), mkdir(savefolder); end

disp('data collecting...');
getData(p_threshold,region{1},folders{1},shuffNum,parent,dataPath,savefolder,'SpikeData_AC.mat');
getData(p_threshold,region{2},folders{2},shuffNum,parent,dataPath,savefolder,'SpikeData_FOF.mat');
getData(p_threshold,region{3},folders{3},shuffNum,parent,dataPath,savefolder,'SpikeData_MPFC.mat');

%% plot figure %%
CompareSequence('SpikeData_AC.mat',shuffNum,savefolder,'AC example','AC recuruit', figsaveTYPE);
CompareSequence('SpikeData_FOF.mat',shuffNum,savefolder,'FOF example','FOF recuruit', figsaveTYPE);
CompareSequence('SpikeData_MPFC.mat',shuffNum,savefolder,'MPFC example','MPFC recuruit', figsaveTYPE);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Rcoeff=CompareSequence(savemat, shuffNum, savefolder, figname1, figname2, figsaveTYPE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(savefolder);
load(savemat,'grand_average_data','rand_average_data','rand_average_remain_data','sound_pre');

soundPeriod = sound_pre+(-100:999);

%%% the reqruitment of all long-sound respond neuron in sound period %%% 
g_data=zscore(grand_average_data(:,soundPeriod),0,2);
[~,maxOrder] =max(g_data,[],2);
[~,sortOrder]=sort(maxOrder);
g_data = g_data(sortOrder,:);
[~,x_sound]=max(g_data,[],2);
y_sound=1:size(g_data,1);

Rcoeff= zeros(shuffNum,1);
x_sound_test1 = cell(shuffNum,2);
x_sound_test2 = cell(shuffNum,2);
for i=1:shuffNum
    
    test1=rand_average_data{i};
    test_data1=zscore(test1(:,soundPeriod),0,2);
    test_data1=test_data1(sortOrder,:);
    [~,x_max_test1]=max(test_data1(y_sound,:),[],2);
    x_sound_test1{i,2} =x_max_test1';
    
    test2=rand_average_remain_data{i};
    test_data2=zscore(test2(:,soundPeriod),0,2);
    test_data2=test_data2(sortOrder,:);
    [~,x_max_test2]=max(test_data2(y_sound,:),[],2);
    x_sound_test2{i,2} =x_max_test2';
    
    Rcoeff(i)=corr(x_max_test1,x_max_test2,'Type','Spearman');
        
    if(i==1)
        h=figure('Position',[100 100 500 500]);
        
        ax=subplot(1,1,1); hold on
        imagesc(g_data);
        plot(x_sound,y_sound,'r');
        xlim([0,size(g_data,2)]);
        ylim([1,size(g_data,1)]);
        yticks([1,size(g_data,1)]);
        xticks([0,100,size(g_data,2)]);
        xticklabels({'-100','0','1000'});
        xlabel('time from sound onset [ms]');
        ylabel('neuron id');
        title('grand average');
        ax.YDir='reverse';
        
        set(h,'PaperPositionMode','auto');
        print(h,'-r0',figname1,figsaveTYPE);
    end
end

testdata2=cell2mat(x_sound_test1(:,2));

h=figure('Position',[100,100,500,500]);hold on;
errorbar(y_sound,mean(testdata2),std(testdata2),'CapSize',0,'Color',[.5 .5 .5]);
scatter(y_sound,mean(testdata2),5,'k','filled');
plot(y_sound,x_sound,'r','LineWidth',1.5);
ylim([0 1100]);
yticks([0 100 1100]);
yticklabels({'-100','0','1000'});
xlim([1 size(testdata2,2)]);
xticks([1 size(testdata2,2)]);
title({['corrcoef/ average: ', num2str(round(mean(Rcoeff),3))],...
       ['          std:     ', num2str(round(std(Rcoeff),3))]});
xlabel('neuron id');
ylabel('time from sound onset [ms]');
view(90,90);

set(h,'PaperPositionMode','auto');
print(h,'-r0',figname2,figsaveTYPE);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function getData(p_threshold,region,folders,shuffNum,parent,dataPath,savefolder,savename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(region);
cd([parent,dataPath]);
datalist = dir([region,'*']);

%%% get target neuron id %%%
sig_neuron = getPval_long(p_threshold,folders);
cd([parent,dataPath]);
sigNeuronIDEachSession = cell(length(datalist),1);
for i=1:length(datalist)
    example = matfile(datalist(i).name);
    NeuronNum=size(example,'neuron',2);
    sigNeuronIDEachSession{i}=sig_neuron(sig_neuron-NeuronNum<=0);
    sig_neuron=sig_neuron-NeuronNum;
    sig_neuron(sig_neuron<=0)=[];
end

%%% spike trace data %%%
grand_average=cell(length(datalist),1); 
rand_aveData=cell(length(datalist),shuffNum);
rand_aveData_remain=cell(length(datalist),shuffNum);

load(datalist(1).name,'frame');
sound_pre=frame.sound_pre;
% for i=1:length(datalist)
parfor i=1:length(datalist)
    if(~isempty(sigNeuronIDEachSession{i}))
        Data= getSequence(datalist(i).name,shuffNum,sigNeuronIDEachSession{i});
        grand_average{i}= Data.grand_average_data;
        tmp_shuff=Data.randomshuff_average;
        tmp_shuff_remain=Data.randomshuff_average_remain;
        for j=1:shuffNum
            rand_aveData{i,j} = tmp_shuff{j};
            rand_aveData_remain{i,j} = tmp_shuff_remain{j};
        end
    end
end
grand_ave_data=cell2mat(grand_average);
[~,maxOrder] =max(grand_ave_data,[],2);
[~,sortOrder]=sort(maxOrder);
grand_average_data = grand_ave_data(sortOrder,:);

rand_average_data=cell(shuffNum,1);
rand_average_remain_data=cell(shuffNum,1);
for i=1:shuffNum
    tmpmat1 = cell2mat(rand_aveData(:,i));
    rand_average_data{i}=tmpmat1(sortOrder,:);
    tmpmat2 = cell2mat(rand_aveData_remain(:,i));
    rand_average_remain_data{i}=tmpmat2(sortOrder,:);
end

cd(savefolder);
save(savename,'grand_average_data','rand_average_data','rand_average_remain_data','sound_pre','-v7.3');
% save(savename,'grand_average_data','rand_average_data','-v7.3');

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Data = getSequence(dataname,shuffNum,sigid)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(dataname,'neuron');
% load(dataname,'neuron','sound_pre');

% get trial length %
for i=1:length(neuron)
    if(~isempty(neuron(i).spike_trace))
    trialNum= size(neuron(i).spike_trace,1);% use trial length
    end
end

shuffTrial = cell(shuffNum,1);% sharing the shuffling trials in the same session
for k=1:shuffNum
    shuffTrial{k}=randperm(trialNum,round(trialNum/2));
end

grand_average_data = cell(length(sigid),1);
randomshuff_average= cell(length(sigid),shuffNum);
randomshuff_average_remain= cell(length(sigid),shuffNum);
for i=1:length(sigid)
    spike_trace= neuron(sigid(i)).spike_trace;
    norm_spike = shuffleTimeCourse(spike_trace);
    grand_average_data{i} = mean(norm_spike,1);
    for n=1:shuffNum
        tmp=setdiff(1:trialNum,shuffTrial{n});
        randomshuff_average{i,n} = mean(norm_spike(shuffTrial{n},:),1);
        randomshuff_average_remain{i,n} = mean(norm_spike(tmp,:),1);
    end
end
randshuff = cell2mat(randomshuff_average);
randshuff_remain = cell2mat(randomshuff_average_remain);
neuronNum = length(grand_average_data);
time = size(grand_average_data{1},2);
Data.grand_average_data = cell2mat(grand_average_data);
Data.randomshuff_average= mat2cell(randshuff,neuronNum,time*ones(1,shuffNum));
Data.randomshuff_average_remain= mat2cell(randshuff_remain,neuronNum,time*ones(1,shuffNum));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function shuff_spike_trace = shuffleTimeCourse(spike_trace)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trial= size(spike_trace,1);
time = size(spike_trace,2);

shuff_spike_trace=zeros(trial,time);
for i=1:trial    
    shuff_spike_trace(i,:)=spike_trace(i,randperm(time));
end

end