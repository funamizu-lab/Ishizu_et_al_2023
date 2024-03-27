
%{
----------------------------------------------------------------------------
Determine the time window for task relevant neurons
%Start of task
%Sound on
%Sound off
%Before choice
%After choice (0sec)
%After choice (1sec)
%After choice (2sec)
%p = 0.001
%Each epoch, predict the prior, sensory and choice (integration)
----------------------------------------------------------------------------
%}
function Figure2ef_process_20220818_prior_process_value_all2

close all
analysis_folder = {
    'G:\Ishizu_data\Tokyo_ephys_ishizu\only_all_behaviors\a04_behave\ML_CumGauss_Fit\RL_20220818'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\only_all_behaviors\a08_behave\ML_CumGauss_Fit\RL_20220818'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\only_all_behaviors\i20_behave\ML_CumGauss_Fit\RL_20220818'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\only_all_behaviors\i24_behave\ML_CumGauss_Fit\RL_20220818'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\only_all_behaviors\i34_behave\ML_CumGauss_Fit\RL_20220818'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\only_all_behaviors\i35_behave\ML_CumGauss_Fit\RL_20220818'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\only_all_behaviors\i43_behave\ML_CumGauss_Fit\RL_20220818'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\only_all_behaviors\i46_behave\ML_CumGauss_Fit\RL_20220818'
    };

use_para = 3;

longRs = [];
longLs = [];
longR = [];
longL = [];
shortR = [];
shortL = [];
para = [];
prior_longR = [];
prior_longL = [];
prior_shortR = [];
prior_shortL = [];
for i = 1:length(analysis_folder)
    [longRs1,longLs1,longR1,longL1,shortR1,shortL1,para1,...
        prior_longR1,prior_longL1,prior_shortR1,prior_shortL1] = ...
        get_process_20220818_prior_process_value2(analysis_folder{i}, use_para);
    
    longRs = [longRs; longRs1]; %long sound trials @ short duration posterior
    longLs = [longLs; longLs1];
    longR = [longR; longR1];  %long sound trials @ long duration posterior
    longL = [longL; longL1];
    shortR = [shortR; shortR1]; %short sound trials @ short duration posterior
    shortL = [shortL; shortL1];
    para = [para; para1];
    prior_longR = [prior_longR; prior_longR1];
    prior_longL = [prior_longL; prior_longL1];
    prior_shortR = [prior_shortR; prior_shortR1];
    prior_shortL = [prior_shortL; prior_shortL1];
    
    [Ndata(i,1),~] = size(longRs1);
    sub_longRs(i,:) = mean(longRs1);
    sub_longLs(i,:) = mean(longLs1);
    sub_longR(i,:) = mean(longR1);
    sub_longL(i,:) = mean(longL1);
    sub_shortR(i,:) = mean(shortR1);
    sub_shortL(i,:) = mean(shortL1);
    sub_para(i,:) = mean(para1);
    sub_prior_longR(i,:) = mean(prior_longR1);
    sub_prior_longL(i,:) = mean(prior_longL1);
    sub_prior_shortR(i,:) = mean(prior_shortR1);
    sub_prior_shortL(i,:) = mean(prior_shortL1);
end
all_subject = [];
for i = 1:length(analysis_folder)
    temp = ones(Ndata(i),1) * i;
    all_subject = [all_subject; temp];
end

%% fig2e
cd('G:\upload_code\Figure2\Fig2ef');
temp_x = [0 0.25 0.45 0.55 0.75 1];
figure
sdata = struct();% source data 
subplot(1,3,1); hold on
plot_mean_se_bar(temp_x,prior_longR,[1 0 0],2)
plot_mean_se_bar(temp_x,prior_longL,[0 0 1],2)
set(gca,'xlim',[-0.1 1.1],'ylim',[0 1])
sdata.x=temp_x';
sdata.beforeSound_R_mean= transpose(mean(prior_longR));
sdata.beforeSound_R_se  = transpose(std(prior_longR)./ sqrt(size(prior_longR,1)));
sdata.beforeSound_L_mean= transpose(mean(prior_longL));
sdata.beforeSound_L_se  = transpose(std(prior_longL)./ sqrt(size(prior_longL,1)));

subplot(1,3,2); hold on
% plot_mean_se_bar(temp_x,shortR,[1 0 0],i)
% plot_mean_se_bar(temp_x,shortL,[0 0 1],i)
% set(gca,'xlim',[-0.1 1.1],'ylim',[0 1])
% sdata.short_R_mean= transpose(mean(shortR));
% sdata.short_R_se  = transpose(std(shortR)./ sqrt(size(shortR,1)));
% sdata.short_L_mean= transpose(mean(shortL));
% sdata.short_L_se  = transpose(std(shortL)./ sqrt(size(shortL,1)));
plot_mean_se_bar(temp_x,longRs,[1 0 0],i)
plot_mean_se_bar(temp_x,longLs,[0 0 1],i)
set(gca,'xlim',[-0.1 1.1],'ylim',[0 1])
sdata.short_R_mean= transpose(mean(longRs));
sdata.short_R_se  = transpose(std(longRs)./ sqrt(size(longRs,1)));
sdata.short_L_mean= transpose(mean(longLs));
sdata.short_L_se  = transpose(std(longLs)./ sqrt(size(longLs,1)));

subplot(1,3,3); hold on
plot_mean_se_bar(temp_x,longR,[1 0 0],2)
plot_mean_se_bar(temp_x,longL,[0 0 1],2)
set(gca,'xlim',[-0.1 1.1],'ylim',[0 1])
sdata.long_R_mean= transpose(mean(longR));
sdata.long_R_se  = transpose(std(longR)./ sqrt(size(longR,1)));
sdata.long_L_mean= transpose(mean(longL));
sdata.long_L_se  = transpose(std(longL)./ sqrt(size(longL,1)));
T = struct2table(sdata);
writetable(T, 'source fig2e.csv');

%% fig 2f
%About choice: Difference between left and right correct trials
longS = get_choice_prob(longRs,longLs);
long  = get_choice_prob(longR, longL);
short = get_choice_prob(shortR,shortL);
p_long= get_choice_prob(prior_longR, prior_longL);
% lme_choice_sabun1 = fitlme_analysis_20210520_0(p_long-longS,all_subject);
% lme_choice_sabun2 = fitlme_analysis_20210520_0(long-longS,all_subject);
% lme_choice_sabun3 = fitlme_analysis_20210520_0(short-long,all_subject);


%About prior: Difference between red and blue blocks
sabun_prior_long = get_prior_prob(prior_longR, prior_longL);
sabun_longs= get_prior_prob(longRs,longLs);
sabun_long = get_prior_prob(longR, longL);
sabun_short= get_prior_prob(shortR,shortL);
% lme_prior_sabun1 = fitlme_analysis_20210520_0(sabun_prior_long-sabun_longs,all_subject);
% lme_prior_sabun2 = fitlme_analysis_20210520_0(sabun_prior_long-sabun_long,all_subject);
% lme_prior_sabun3 = fitlme_analysis_20210520_0(sabun_longs-sabun_long,all_subject);
% lme_prior_sabun4 = fitlme_analysis_20210520_0(sabun_short-sabun_long,all_subject);

figure
sdata1 = struct();% source data 
sdata2 = struct();% source data 
subplot(1,2,1);hold on
boxplot([p_long,longS,long])
plot([p_long,longS,long]')
sdata1.before =p_long;
sdata1.short =longS;
sdata1.long  =long; 
T = struct2table(sdata1);
writetable(T, 'source fig2f choice.csv');

subplot(1,2,2);hold on
boxplot([sabun_prior_long,sabun_longs,sabun_long])
plot([sabun_prior_long,sabun_longs,sabun_long]')
sdata2.before =sabun_prior_long;
sdata2.short =sabun_longs;
sdata2.long  =sabun_long; 
T = struct2table(sdata2);
writetable(T, 'source fig2f prior.csv');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function longS = get_choice_prob(longRs, longLs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
longS = (longRs + longLs) ./ 2;
longS = mean(longS(:,[4:6]),2) - mean(longS(:,[1:3]),2);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sabun_prior_long = get_prior_prob(prior_longR, prior_longL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sabun_prior_long = prior_longR - prior_longL;
sabun_prior_long = mean(sabun_prior_long,2);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [longRs,longLs,longR,longL,shortR,shortL,para,...
    prior_longR,prior_longL,prior_shortR,prior_shortL] = ...
    get_process_20220818_prior_process_value2(folders, use_para)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(folders);
filename1 = dir('*.mat');
    
for filecount = 1:length(filename1)
    temp_folder = filename1(filecount).folder;
    temp_filename = filename1(filecount).name;
    RL_file{filecount} = [temp_folder,'\',temp_filename];
end

cd ../../ %Get bpod file
filename2 = dir('Bpod*.mat');
if length(filename1)~= length(filename2)
    hoge
else
    for filecount = 1:length(filename2)
        temp_folder = filename2(filecount).folder;
        temp_filename = filename2(filecount).name;
        Bpod_file{filecount} = [temp_folder,'\',temp_filename];
    end
end

longRs = [];
longLs = [];
longR = [];
longL = [];
shortR = [];
shortL = [];
para = [];
prior_longR = [];
prior_longL = [];
prior_shortR = [];
prior_shortL = [];

for i = 1:length(filename1)    
    [temp_longRs, temp_longLs, temp_longR,temp_longL,temp_shortR,temp_shortL,temp_para,...
        temp_prior_longR,temp_prior_longL,temp_prior_shortR,temp_prior_shortL] = ...
        Task_kaiseki_tokyo1_20220818_prior_process_value2(RL_file{i},Bpod_file{i}, use_para);
    close all

    longRs = [longRs; temp_longRs];
    longLs = [longLs; temp_longLs];
    longR = [longR; temp_longR];
    longL = [longL; temp_longL];
    shortR = [shortR; temp_shortR];
    shortL = [shortL; temp_shortL];
    para = [para; temp_para];
    
    prior_longR = [prior_longR; temp_prior_longR];
    prior_longL = [prior_longL; temp_prior_longL];
    prior_shortR = [prior_shortR; temp_prior_shortR];
    prior_shortL = [prior_shortL; temp_prior_shortL];
end
delete(gcp('nocreate'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_mean_se_bar(temp_x,prior_longR,use_color,std_base)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[size_y,size_x] = size(prior_longR);
mean_prior = mean(prior_longR);
std_prior = std(prior_longR);
se_prior = std_prior ./ sqrt(size_y);

plot(temp_x,mean_prior,'color',use_color)
hold on
if std_base == 1
    for i = 1:size_x
        plot([temp_x(i) temp_x(i)],[mean_prior(i)+std_prior(i), mean_prior(i)-std_prior(i)],'color',use_color)
    end
else
    for i = 1:size_x
        plot([temp_x(i) temp_x(i)],[mean_prior(i)+se_prior(i), mean_prior(i)-se_prior(i)],'color',use_color)
    end
end