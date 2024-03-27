%{
----------------------------------------------------------------------------
Analyzing behavioral data
At least for the correct rate
----------------------------------------------------------------------------
%}

function B = Figure1e_Dual_behave_analysis_block1_220818_lick_short

analysis_folder = {
    'G:\Ishizu_data\IntermediateFiles\only_all_behaviors\a04_behave\'
    'G:\Ishizu_data\IntermediateFiles\only_all_behaviors\a08_behave\'
    'G:\Ishizu_data\IntermediateFiles\only_all_behaviors\i20_behave\'
    'G:\Ishizu_data\IntermediateFiles\only_all_behaviors\i24_behave\'
    'G:\Ishizu_data\IntermediateFiles\only_all_behaviors\i34_behave\'
    'G:\Ishizu_data\IntermediateFiles\only_all_behaviors\i35_behave\'
    'G:\Ishizu_data\IntermediateFiles\only_all_behaviors\i43_behave\'
    'G:\Ishizu_data\IntermediateFiles\only_all_behaviors\i46_behave\'
    };

%Long: 1 - 100
%Short: 1 - 40
time_bin = [10;10;10;10;20;10;20;10];
x_bin = [0 3];

B = [];
B_prior = [];
B_min = [];
B_prior_min = [];
all_subject = [];
for i = 1:length(analysis_folder)
    temp_folder = analysis_folder{i};
    cd(temp_folder);
    filename1 = dir('Bpod*.mat');
    
    clear temp_file
    for j = 1:length(filename1)
        temp_file{j} = filename1(j).name;
    end
    
    clear file_B file_B_prior file_B_min file_B_prior_min
    for filecount = 1 : length(temp_file)
        
        fpath = temp_file(filecount);
        fpath = cell2mat(fpath);
        if i == 8 && filecount == 1 %This is for i46, Jun01
            [temp_B_max, temp_B_prior_max, stats_max, temp_B_min, temp_B_prior_min, stats_min] = ...
                behave_analysis_block1_201028_lick_each(fpath, 20);
        else
            [temp_B_max, temp_B_prior_max, stats_max, temp_B_min, temp_B_prior_min, stats_min] = ...
                behave_analysis_block1_201028_lick_each(fpath, time_bin(i));
        end
        
        close all
        file_B(filecount,:) = temp_B_max';
        file_B_prior(filecount,:) = temp_B_prior_max';
        file_B_min(filecount,:) = temp_B_min';
        file_B_prior_min(filecount,:) = temp_B_prior_min';
        %B(filecount,:) = -B(filecount,:);
    end
    B = [B; file_B];
    B_prior = [B_prior; file_B_prior];
    B_min = [B_min; file_B_min];
    B_prior_min = [B_prior_min; file_B_prior_min];
    
    %Get the subject number
    temp_subject = ones(length(temp_file),1) * i;
    all_subject = [all_subject; temp_subject];
end

plot_B_prior2(B_prior_min,x_bin,all_subject);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_B_prior2(B_prior,x_bin,all_subject)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot about B_prior
[size_session,size_B] = size(B_prior);

%%% fig 1e short %%%
cd('G:\upload_code\Figure1\Fig1e');
figure
[mean_trace,~,se_trace] = plot_mean_se_moto(B_prior(:,2:size_B-1),[0 0 0],2); %only sound
set(gca,'xlim',x_bin)
sdata.mean=mean_trace';
sdata.se =se_trace';
T = struct2table(sdata);
writetable(T, 'source fig1e short.csv');

% p = kruskalwallis(B_prior(:,[2:size_B-1]),[],'off');
% for i = 1:max(all_subject)
%     temp_B = B_prior(all_subject == i,:);
%     B_subject(i,:) = mean(temp_B);
% end
% 
% for i = 2:size_B-1
%     p_each_sound(i-1) = signrank(B_prior(:,i));
%     lme_each_sound(i-1).matrix = fitlme_analysis_20210520_0(B_prior(:,i),all_subject);
% 
%     p_each_subject(i-1) = signrank(B_subject(:,i));
% end
% p_prior = signrank(B_prior(:,size_B));
% lme_prior = fitlme_analysis_20210520_0(B_prior(:,size_B),all_subject);
% %p_each_prior
% p_prior_subject = signrank(B_subject(:,size_B));

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [B_max, B_prior_max, stats_max, B_min, B_prior_min, stats_min] = ...
    behave_analysis_block1_201028_lick_each(filename1, time_bin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%[filename1, pathname1]=uigetfile('*.mat','Block_mat');
load(filename1)

[minD_trial,maxD_trial,Choice_trial,tone_evidence,trial_evidence,use_trial2,use_trial3,use_trial_all,...
    low,high,correct,error,flip_tone,number_use_trial,...
    binary_tone,right_trial_all,number_trial_all,right_trial,number_trial] ...
 = Dual_get_basic_task_structure_20210204(filename1);

%Based on the binary tone decide the pseudo tone evidence
% temp_tone(1).matrix = find(binary_tone == 0);
% temp_tone(2).matrix = find(binary_tone > 0 & binary_tone <= 0.35);
% temp_tone(3).matrix = find(binary_tone > 0.35 & binary_tone < 0.5);
% temp_tone(4).matrix = find(binary_tone > 0.5 & binary_tone < 0.65);
% temp_tone(5).matrix = find(binary_tone >= 0.65 & binary_tone < 1);
% temp_tone(6).matrix = find(binary_tone == 1);

%Make block bias
block_bias = zeros(length(Choice_trial),1);
%Correct Choice Evidence
if BlockReward(2,1) < BlockReward(2,2) % Right -> Left
    block_bias(use_trial2) = 1;
    block_bias(use_trial3) = -1;
else % Left -> Right
    block_bias(use_trial2) = -1;
    block_bias(use_trial3) = 1;
end

%Pick the trials to use the regression
min_trial = intersect(use_trial_all, minD_trial);

[B_min, B_prior_min, stats_min] = ...
    get_time_based_choice_regression(min_trial, Tone_cloud, Correct_side, Chosen_side, trial_evidence, block_bias, time_bin);
B_max = nan;
B_prior_max = nan;
stats_max = nan;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [B, B_prior, stats] = ...
    get_time_based_choice_regression(max_trial, Tone_cloud, Correct_side, Chosen_side, trial_evidence, block_bias, time_bin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tone_cloud = Tone_cloud(max_trial);
Correct_trial = Correct_side(max_trial);
Chosen_trial = Chosen_side(max_trial);
trial_evidence = trial_evidence(max_trial);
block_bias = block_bias(max_trial);

clear temp_tone
temp_tone(1).matrix = find(trial_evidence == 0);
temp_tone(2).matrix = find(trial_evidence > 0 & trial_evidence <= 0.35);
temp_tone(3).matrix = find(trial_evidence > 0.35 & trial_evidence < 0.5);
temp_tone(4).matrix = find(trial_evidence > 0.5 & trial_evidence < 0.65);
temp_tone(5).matrix = find(trial_evidence >= 0.65 & trial_evidence < 1);
temp_tone(6).matrix = find(trial_evidence == 1);
mid_trial  = [temp_tone(2).matrix; temp_tone(5).matrix];
dif_trial  = [temp_tone(3).matrix; temp_tone(4).matrix];

%update the tone cloud
for i = 1:length(Tone_cloud)
    tone_cloud(i,:) = Tone_cloud(i).matrix;
end
clear Tone_cloud

%Use time bin to make regression value
%binarize the tone cloud
%low  1-6
%high 13-18
[size_y,size_x] = size(tone_cloud);
temp_size = ceil(size_x/time_bin);
for i = 1:temp_size
    temp_min = (i-1)*time_bin+1;
    temp_max = i*time_bin;
    temp_max = min(temp_max, size_x);
    bin_use(i).matrix = temp_min : temp_max;
end
length_bin = length(bin_use);

for i = 1:size_y
    temp_tone = tone_cloud(i,:);
    for j = 1:length_bin
        temp_temp_tone = temp_tone(bin_use(j).matrix);
        temp0 = find(temp_temp_tone <= 8);
        temp1 = find(temp_temp_tone >= 9);
        binary_tone(i,j) = (length(temp1)-length(temp0)) ./ length(bin_use(j).matrix);
    end
    
    %Get the data in all sound
    temp1 = find(temp_tone >= 9);
    binary_all_tone(i,1) = length(temp1) ./ length(temp_tone);
end

%Based on the correct trial, flip the tone cloud
if mean(binary_all_tone(Correct_trial == 1)) < 0.5 %low for right correct
    disp('flip tones')
    binary_tone = -binary_tone;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%start regression
%Use difficult stimuli
Evi_diff = sort([mid_trial; dif_trial]);

binary_tone = binary_tone(Evi_diff,:);
Chosen_trial = Chosen_trial(Evi_diff);
block_bias = block_bias(Evi_diff); %-1 or 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Normalize each parameter
[size_y,size_x] = size(binary_tone);
temp_binary_tone = reshape(binary_tone,size_y*size_x,1);
mean_tone = mean(temp_binary_tone);
std_tone = std(temp_binary_tone);
for j = 1:length_bin
    binary_tone(:,j) = (binary_tone(:,j) - mean_tone) ./ std_tone;
end
block_bias = (block_bias - mean(block_bias)) ./ std(block_bias);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Pick only difficult trials
B = glmfit(binary_tone,Chosen_trial,'binomial','link','logit');
[B_prior,~,stats] = glmfit([binary_tone,block_bias],Chosen_trial,'binomial','link','logit');
return
