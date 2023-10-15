function FigureS10a_drawEachNeuronTrace_TCchoice
parent='G:/Ishizu_data';
dataPath='/Revise_ishizu/output/SpikeTrace_longtone';
outPath='/Revise_ishizu/output/drawEachNeuronTrace_TCchoice';

% hulistic parameter %
p_threshold = 1.0e-10;
figsaveTYPE='-dpng';
%--------------------%


%% collect the data of firing rate %%%
% folder setting %
region = {'auditory','fof','mpfc'};
folders= {'auc_ishizu','fof_ishizu','mpfc_ishizu'};
% save folder setting %
savefolder =[parent,outPath];
if(~exist(savefolder,'dir')), mkdir(savefolder); end
cd(savefolder);

disp('data collecting...');
% getData(p_threshold,region{1},folders{1},parent,dataPath,savefolder,'SpikeData_AC.mat');

%% plot figure %%
% drawData('SpikeData_AC.mat',savefolder,'AC example', figsaveTYPE);

%% plot figure %%
drawData_single(141,'SpikeData_AC.mat',savefolder,'AC', figsaveTYPE);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function drawData_single(i, savemat, savefolder, figname, figsaveTYPE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(savefolder);
load(savemat,'sound_pre','neuron_label','evi_correct_ALL','evi_error_ALL',...
    'evi_correct_step_ALL','evi_correct_stepSE_ALL',...
    'evi_error_step_ALL','evi_error_stepSE_ALL','evi_step_pval_ALL','nanflag_all','flipflag_all');

% sound_pre=1500; %msec
neuron_label=neuron_label(nanflag_all==0);
flip_label=flipflag_all(nanflag_all==0);
evi_correct_ALL= evi_correct_ALL(nanflag_all==0,:);
evi_error_ALL  = evi_error_ALL(nanflag_all==0,:);
evi_correct_step_ALL  =evi_correct_step_ALL(nanflag_all==0,:);
evi_correct_stepSE_ALL=evi_correct_stepSE_ALL(nanflag_all==0,:);
evi_error_step_ALL  =evi_error_step_ALL(nanflag_all==0,:);
evi_error_stepSE_ALL=evi_error_stepSE_ALL(nanflag_all==0,:);
evi_step_pval_ALL   =evi_step_pval_ALL(nanflag_all==0,:);

t=25;
x=[0,0.25,0.45,0.55,0.75,1];
col=jet(6);
h=figure('Position',[100 100 1500 500]);
subplot(1,3,1); hold on;
for j=1:6
    if(flip_label(i)==0)
    plot(evi_correct_ALL{i,j}*1000,'Color',col(j,:));
    else        
    plot(evi_correct_ALL{i,6-j+1}*1000,'Color',col(j,:));
    end
end
title('correct')
xticks([sound_pre-500,sound_pre,sound_pre+500,sound_pre+1000,sound_pre+1500,sound_pre+2000]);
xticklabels({'','0','','1','','2'});
xlim([sound_pre-500,sound_pre+2000])
xlabel('Time from sound onset [s]');
ylabel('Spike (Hz)');

subplot(1,3,2); hold on;
for j=1:6
    plot(evi_error_ALL{i,j}*1000,'Color',col(j,:));
end
title('error')
xticks([sound_pre-500,sound_pre,sound_pre+500,sound_pre+1000,sound_pre+1500,sound_pre+2000]);
xticklabels({'','0','','1','','2'});
xlim([sound_pre-500,sound_pre+2000])
xlabel('Time from sound onset [s]');
ylabel('Spike (Hz)');

subplot(1,3,3); hold on;
corre=zeros(1,6);  c_se=zeros(1,6);
error=zeros(1,6);  e_se=zeros(1,6);
pval =zeros(1,6);
for j=1:6
    if(flip_label(i)==0)
    corre(j)= evi_correct_step_ALL{i,j}(t);
    error(j)= evi_error_step_ALL{i,j}(t);
    c_se(j) = evi_correct_stepSE_ALL{i,j}(t);
    e_se(j) = evi_error_stepSE_ALL{i,j}(t);
    pval(j) = evi_step_pval_ALL{i,j}(t);
    else        
    corre(j)= evi_correct_step_ALL{i,6-j+1}(t);
    error(j)= evi_error_step_ALL{i,6-j+1}(t);
    c_se(j) = evi_correct_stepSE_ALL{i,6-j+1}(t);
    e_se(j) = evi_error_stepSE_ALL{i,6-j+1}(t);
    pval(j) = evi_step_pval_ALL{i,6-j+1}(t);
    end
end
errorplot(x,corre,c_se,c_se,'r',0.5,1);
errorplot(x,error,e_se,e_se,'b',0.5,1);
xticks(0:0.2:1)
xlim([-0.1, 1.1])
title({num2str(t),neuron_label{i},['pval: ',num2str(pval)]})

set(h,'PaperPositionMode','auto');
print(h,'-r0',figname,figsaveTYPE);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function drawData(savemat, savefolder, figname, figsaveTYPE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(savefolder);
load(savemat,'sound_pre','evi_correct_ALL','evi_error_ALL',...
    'evi_correct_step_ALL','evi_correct_stepSE_ALL',...
    'evi_error_step_ALL','evi_error_stepSE_ALL','evi_step_pval_ALL','nanflag_all','flipflag_all');

% sound_pre=1500;

NeuronNum = length(find(nanflag_all==0));
flip_label=flipflag_all(nanflag_all==0);
evi_correct_ALL= evi_correct_ALL(nanflag_all==0,:);
evi_error_ALL  = evi_error_ALL(nanflag_all==0,:);
evi_correct_step_ALL  =evi_correct_step_ALL(nanflag_all==0,:);
evi_correct_stepSE_ALL=evi_correct_stepSE_ALL(nanflag_all==0,:);
evi_error_step_ALL  =evi_error_step_ALL(nanflag_all==0,:);
evi_error_stepSE_ALL=evi_error_stepSE_ALL(nanflag_all==0,:);
evi_step_pval_ALL   =evi_step_pval_ALL(nanflag_all==0,:);

t=25;
x=[0,0.25,0.45,0.55,0.75,1];

tmp=extractBetween(savemat,'_','.');
mkdir(tmp{1});
cd(tmp{1});
col=jet(6);
for i=1:NeuronNum
    %     h=figure('visible','off');
    h=figure('Position',[100 100 1500 500],'visible','off');
    subplot(1,3,1); hold on;
    for j=1:6
    if(flip_label(i)==0)
        plot(evi_correct_ALL{i,j}*1000,'Color',col(j,:));
    else        
        plot(evi_correct_ALL{i,6-j+1}*1000,'Color',col(j,:));
    end
    end
    title('correct')
    xticks([sound_pre-500,sound_pre,sound_pre+500,sound_pre+1000,sound_pre+1500,sound_pre+2000]);
    xticklabels({'','0','','1','','2'});
    xlim([sound_pre-500,sound_pre+2000])
    xlabel('Time from sound onset [s]');
    ylabel('Spike (Hz)');
    
    subplot(1,3,2); hold on;
    for j=1:6
    if(flip_label(i)==0)
        plot(evi_error_ALL{i,j}*1000,'Color',col(j,:));
    else
        plot(evi_error_ALL{i,6-j+1}*1000,'Color',col(j,:));
    end
    end
    title('error')
    xticks([sound_pre-500,sound_pre,sound_pre+500,sound_pre+1000,sound_pre+1500,sound_pre+2000]);
    xticklabels({'','0','','1','','2'});
    xlim([sound_pre-500,sound_pre+2000])
    xlabel('Time from sound onset [s]');
    ylabel('Spike (Hz)');
    
    subplot(1,3,3); hold on;
    corre=zeros(1,6);  c_se=zeros(1,6);
    error=zeros(1,6);  e_se=zeros(1,6);
    pval=zeros(1,6);
    for j=1:6
    if(flip_label(i)==0)
        corre(j)= evi_correct_step_ALL{i,j}(t);
        error(j)= evi_error_step_ALL{i,j}(t);
        c_se(j) = evi_correct_stepSE_ALL{i,j}(t);
        e_se(j) = evi_error_stepSE_ALL{i,j}(t);
        pval(j) = evi_step_pval_ALL{i,j}(t);
    else        
        corre(j)= evi_correct_step_ALL{i,6-j+1}(t);
        error(j)= evi_error_step_ALL{i,6-j+1}(t);
        c_se(j) = evi_correct_stepSE_ALL{i,6-j+1}(t);
        e_se(j) = evi_error_stepSE_ALL{i,6-j+1}(t);
        pval(j) = evi_step_pval_ALL{i,6-j+1}(t);
    end
    end
    errorplot(x,corre,c_se,c_se,'r',0.5,1);
    errorplot(x,error,e_se,e_se,'b',0.5,1);
    xticks(0:0.2:1)
    xlim([-0.1, 1.1])
    title({num2str(t),['pval: ',num2str(pval)]})
    
    set(h,'PaperPositionMode','auto');
    print(h,'-r0',[figname,'_',num2str(i)],figsaveTYPE);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function getData(p_threshold,region,folders,parent,dataPath,savefolder,savename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(region);

%%% get target neuron id %%%
sig_neuron = getPval_long(p_threshold,folders);
neuron_label = cell(length(sig_neuron),1);
cd([parent,dataPath]);
datalist = dir([region,'*']);
sigNeuronIDEachSession = cell(length(datalist),1);
for i=1:length(datalist)
    example = matfile(datalist(i).name);
    NeuronNum=size(example,'neuron',2);
    sigNeuronIDEachSession{i}=sig_neuron(sig_neuron-NeuronNum<=0);
    sig_neuron=sig_neuron-NeuronNum;
    sig_neuron(sig_neuron<=0)=[];
end

%%% spike trace data %%%
evi_correct= cell(length(datalist),1);
evi_error  = cell(length(datalist),1);
evi_correct_step_ave= cell(length(datalist),1);
evi_error_step_ave  = cell(length(datalist),1);
evi_correct_step_se= cell(length(datalist),1);
evi_error_step_se  = cell(length(datalist),1);
evi_step_pval  = cell(length(datalist),1);
nanflag= cell(length(datalist),1);
flipflag= cell(length(datalist),1);
load(datalist(1).name,'frame');
sound_pre=frame.sound_pre;
id=1;
for i=1:length(datalist)
    if(~isempty(sigNeuronIDEachSession{i}))
        for j=1:length(sigNeuronIDEachSession{i})
            neuron_label{id}=[extractBefore(datalist(i).name,'.'),'_n',num2str(sigNeuronIDEachSession{i}(j))];
            id=id+1;
        end
        [Data,Time,step]= getSequence(datalist(i).name,sigNeuronIDEachSession{i});
        
        evi_correct{i}= Data.evi_correct;
        evi_error{i}  = Data.evi_error;
        evi_correct_step_ave{i}= Data.evi_correct_step;
        evi_error_step_ave{i}  = Data.evi_error_step;
        evi_correct_step_se{i} = Data.evi_correct_stepSE;
        evi_error_step_se{i}   = Data.evi_error_stepSE;
        evi_step_pval{i}       = Data.evi_step_pval;
        nanflag{i} = Data.nanflag;
        flipflag{i} = Data.flipflag;
    end
end
evi_correct_all= cell2mat(evi_correct);
evi_error_all  = cell2mat(evi_error);
evi_correct_step_all= cell2mat(evi_correct_step_ave);
evi_error_step_all  = cell2mat(evi_error_step_ave);
evi_correct_stepSE_all= cell2mat(evi_correct_step_se);
evi_error_stepSE_all  = cell2mat(evi_error_step_se);
evi_step_pval_all     = cell2mat(evi_step_pval);
n=size(evi_correct_all,1);
evi_correct_ALL= mat2cell(evi_correct_all,ones(1,n),Time*ones(1,6));
evi_error_ALL  = mat2cell(evi_error_all,ones(1,n),Time*ones(1,6));
evi_correct_step_ALL= mat2cell(evi_correct_step_all,ones(1,n),step*ones(1,6));
evi_error_step_ALL  = mat2cell(evi_error_step_all,ones(1,n),step*ones(1,6));
evi_correct_stepSE_ALL= mat2cell(evi_correct_stepSE_all,ones(1,n),step*ones(1,6));
evi_error_stepSE_ALL  = mat2cell(evi_error_stepSE_all,ones(1,n),step*ones(1,6));
evi_step_pval_ALL     = mat2cell(evi_step_pval_all,ones(1,n),step*ones(1,6));
nanflag_all = cell2mat(nanflag);
flipflag_all = cell2mat(flipflag);

cd(savefolder);
save(savename,'neuron_label','sound_pre','evi_correct_ALL','evi_error_ALL','evi_correct_step_ALL','evi_correct_stepSE_ALL',...
    'evi_error_step_ALL','evi_error_stepSE_ALL','evi_step_pval_ALL','nanflag_all','flipflag_all','-v7.3');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Data,Time,step] = getSequence(dataname,sigid)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(dataname,'neuron','frame');

time_window=100;%ms
tmp= neuron(1).spike_trace;
Time = size(tmp,2);
step = Time/time_window;

sound_pre= frame.sound_pre;
checkPeriod=sound_pre:sound_pre + time_window;% msec
use_trial= frame.use_long;
outcome  = ismember(use_trial,frame.correct_trial_long);% 1: correct / 0:error
ES_trial = frame.binary_tone_long;

lowcorrect = intersect(find(sum(ES_trial==0,2)==1),find(outcome==1));
highcorrect= intersect(find(sum(ES_trial==1,2)==1),find(outcome==1));

ES = unique(ES_trial);
evi_correct_trial = cell(6,1);
evi_error_trial = cell(6,1);
for j=1:6
    evi_correct_trial{j}=intersect(find(ES_trial==ES(j)),find(outcome==1));
    evi_error_trial{j} = intersect(find(ES_trial==ES(j)),find(outcome==0));
end

evi_correct= cell(length(sigid),6);
evi_error  = cell(length(sigid),6);
evi_correct_step_ave= cell(length(sigid),6);
evi_error_step_ave  = cell(length(sigid),6);
evi_correct_step_se= cell(length(sigid),6);
evi_error_step_se  = cell(length(sigid),6);
evi_step_pval  = cell(length(sigid),6);
nanflag=zeros(length(sigid),1);
flipflag=zeros(length(sigid),1);
for i=1:length(sigid)
    %     spike_trace = zscore(neuron(sigid(i)).spike_trace,0,2);
    spike_trace = neuron(sigid(i)).spike_trace;
    round_spike = zeros(size(spike_trace,1),step);
    for k=1:step
        round_spike(:,k)=mean(spike_trace(:,(1:time_window)+(k-1)*time_window),2);
    end
    %%% check prefer tone side %%%
    high_spike = spike_trace(highcorrect,checkPeriod);
    low_spike  = spike_trace(lowcorrect,checkPeriod);
    if(mean(high_spike(:)) > mean(low_spike(:)))%high prefer
        evi_correct_trial2= evi_correct_trial;
        evi_error_trial2  = evi_error_trial;
    else%low prefer -> flip trial evidence
        flipflag(i)=1;
        evi_correct_trial2= flip(evi_correct_trial);
        evi_error_trial2  = flip(evi_error_trial);
    end
    
    for j=1:6
        if(length(evi_correct_trial2{j})>2)
            evi_correct{i,j}= mean(spike_trace(evi_correct_trial2{j},:),1);
            evi_correct_step_ave{i,j}= mean(round_spike(evi_correct_trial2{j},:),1);
            evi_correct_step_se{i,j} = std(round_spike(evi_correct_trial2{j},:),0,1)/sqrt(length(evi_correct_trial2{j}));
        else
            evi_correct{i,j}= nan(1,Time);
            evi_correct_step_ave{i,j}= nan(1,step);
            evi_correct_step_se{i,j}= nan(1,step);
            nanflag(i)=1;
        end
        if(length(evi_error_trial2{j})>2)
            evi_error{i,j}  = mean(spike_trace(evi_error_trial2{j},:),1);
            evi_error_step_ave{i,j}= mean(round_spike(evi_error_trial2{j},:),1);
            evi_error_step_se{i,j} = std(round_spike(evi_error_trial2{j},:),0,1)/sqrt(length(evi_error_trial2{j}));
        else
            evi_error{i,j}= nan(1,Time);
            evi_error_step_ave{i,j}= nan(1,step);
            evi_error_step_se{i,j}= nan(1,step);
            nanflag(i)=1;
        end
        try
            c= round_spike(evi_correct_trial2{j},:);
            e= round_spike(evi_error_trial2{j},:);
            p=zeros(1,step);
            for k=1:step
                p(k)=ranksum(c(:,k),e(:,k));
            end
            evi_step_pval{i,j}=p;
        catch
            evi_step_pval{i,j}=nan(1,step);
        end
    end
end
Data.evi_correct = cell2mat(evi_correct);
Data.evi_error = cell2mat(evi_error);
Data.evi_correct_step = cell2mat(evi_correct_step_ave);
Data.evi_error_step = cell2mat(evi_error_step_ave);
Data.evi_correct_stepSE = cell2mat(evi_correct_step_se);
Data.evi_error_stepSE = cell2mat(evi_error_step_se);
Data.evi_step_pval = cell2mat(evi_step_pval);
Data.nanflag = nanflag;
Data.flipflag= flipflag;
end
