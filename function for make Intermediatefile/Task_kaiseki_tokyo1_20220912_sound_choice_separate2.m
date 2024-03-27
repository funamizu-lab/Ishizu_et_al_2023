
%{
----------------------------------------------------------------------------
First_take number of frames in each tif files
Analyzing imaging data simply
At least for the correct rate
----------------------------------------------------------------------------
%}
function Task_kaiseki_tokyo1_20220912_sound_choice_separate2(pathname)

switch nargin
    case 0
        pathname = pwd;
    case 1
        disp('OK to analyze')
    otherwise
        hoge
end
cd(pathname)
%spike_dir = 'spike_ch1';
spike_dir = dir('spike_ch*');
if length(spike_dir) ~= 1
    hoge
end
spike_dir = spike_dir.name

%new_p_thre = 100;
new_p_thre = 10;
%new_p_thre = 5;
%new_p_thre = 0;

% [filename1, pathname1,findex]=uigetfile('*.*','frame file');
% filename1 = [pathname1,filename1];
% load(filename1)
temp = dir('task_frame*');
if length(temp) ~= 1
    temp
    hoge
end
load(temp.name);
%frame_start
%frame_sound
%frame_end

temp = dir('Bpod*');
if length(temp) ~= 1
    temp
    hoge
end
load(temp.name);
Bpod_file = temp.name;

[minD_trial,maxD_trial,Choice_trial,tone_evidence,trial_evidence,block2,block3,use_trial,...
    low,high,correct,error,flip_tone,number_use_trial,...
    binary_tone,right_trial_all,number_trial_all,right_trial,number_trial] ...
 = Dual_get_basic_task_structure_20210204(Bpod_file);

temp = dir('RL_20220818*'); %RL_20220118 frame
if length(temp) ~= 1
    temp
    hoge
end
load(temp.name);
%'ave_likeli' ,'BIC','log_likeli','para_max','N_trial'
%Get parameter of RL model
[prior,posterior,ave_likeli,likelihood,Long_posterior,Short_posterior,...
    prior_value,long_value,short_value,para] = ...
    Dual_RL_model_block1_20220314_para_determined(Bpod_file,para_max(3,:),N_trial);

prior = prior(:,2) ./ sum(prior(:,1)+prior(:,2));
posterior = posterior(:,2) ./ sum(posterior(:,1)+posterior(:,2));
Q = prior_value;
if length(Q) ~= length(prior)
    [length(Q), length(prior)]
    hoge
end
prior_value = prior_value(:,2) ./ (prior_value(:,1)+prior_value(:,2)); %relative value
binary_prior = zeros(length(prior_value),1);
temp = find(prior_value >= 0.5);
binary_prior(temp) = 1;

%Setup the prior and posterior
if length(binary_prior) ~= length(use_trial)
    hoge
end
% if length(binary_posterior) ~= length(use_trial)
%     hoge
% end
all_binary_prior = nan(length(Choice_trial),1);
all_binary_prior(use_trial) = binary_prior;
% temp = dir('sig_task_neurons_2021*');
% if length(temp) ~= 1
%     temp
%     hoge
% end
% load(temp.name);

%Get task parameter
Choice_trial = find(Outcome == 1 | Outcome == 2);
low  = find(Correct_side == 0);
high = find(Correct_side == 1);
stim_length = unique(StimDuration);
Long  = find(StimDuration == stim_length(2));
Short = find(StimDuration == stim_length(1));
left = find(Chosen_side == 0);
right = find(Chosen_side == 1);
correct_error = Correct_side == Chosen_side;
correct = find(correct_error == 1);
error = find(correct_error == 0);

%Make the choice timing frames
frame_choice_select = nan(length(frame_choice),1);
frame_choice_select(left) = frame_choice(left,1);
frame_choice_select(right) = frame_choice(right,2);

%Adjust the trials to use
correct = intersect(correct, use_trial);
use_long = intersect(Long, use_trial);
use_short = intersect(Short, use_trial);
if BlockReward(2,1) > BlockReward(2,2) %Left -> Right
    Long_R = intersect(Long, block3);
    Short_R = intersect(Short, block3);
    Long_L = intersect(Long, block2);
    Short_L = intersect(Short, block2);
else %Right -> Left
    Long_R = intersect(Long, block2);
    Short_R = intersect(Short, block2);
    Long_L = intersect(Long, block3);
    Short_L = intersect(Short, block3);
end

%Get tone evidence. Put tone evidence in all trials;
for i = 1:length(tone_evidence)
    temp = find(trial_evidence == tone_evidence(i));
    temp = intersect(temp, correct);
    evi_trial(i).matrix = temp; %Only correct trials
end

evidence_group = [];
evidence_group_long = [];
evidence_group_short = [];
for j = 1:length(tone_evidence)
    temp = ones(length(evi_trial(j).matrix),1) * j;
    evidence_group = [evidence_group; temp];
    
    temp_long  = intersect(evi_trial(j).matrix, Long);
    temp = ones(length(temp_long),1) * j;
    evidence_group_long = [evidence_group_long; temp];
    
    temp_short = intersect(evi_trial(j).matrix, Short);
    temp = ones(length(temp_short),1) * j;
    evidence_group_short = [evidence_group_short; temp];
end

time_window_base = 200;
time_window = 100; %ms
time_long_pre  = 1500; %4sec
time_long_post = 2500; %4sec
time_long_pre2  = 500; %4sec
time_long_post2 = 2500; %4sec

time_short_pre = 1500; %3.2sec
time_short_post = 1700; %3.2sec
time_short_pre2 = 500; %3.2sec
time_short_post2 = 2500; %3.2sec

cd(spike_dir);

tif_name = dir('task_spike_stripe*.mat'); %get all the tif files
length_tif = length(tif_name);
max_tif = length_tif;

[long1_correct,long1_error,long1_sound,long1_choice,long1_prior,...
 long2_correct,long2_error,long2_sound,long2_choice,long2_prior, ...
 long1_correct_p,long1_error_p,long1_sound_p,long1_choice_p,long1_prior_p,...
 long2_correct_p,long2_error_p,long2_sound_p,long2_choice_p,long2_prior_p] = ...
    get_task_activ_separate(use_long,frame_sound,frame_choice_select,frame_spout,max_tif,...
    time_long_pre,time_long_post,time_long_pre2,time_long_post2,time_window,time_window_base, ...
    Correct_side, Chosen_side, correct_error, all_binary_prior);

[short1_correct,short1_error,short1_sound,short1_choice,short1_prior,...
 short2_correct,short2_error,short2_sound,short2_choice,short2_prior, ...
 short1_correct_p,short1_error_p,short1_sound_p,short1_choice_p,short1_prior_p,...
 short2_correct_p,short2_error_p,short2_sound_p,short2_choice_p,short2_prior_p] = ...
    get_task_activ_separate(use_short,frame_sound,frame_choice_select,frame_spout,max_tif,...
    time_short_pre,time_short_post,time_short_pre2,time_short_post2,time_window,time_window_base, ...
    Correct_side, Chosen_side, correct_error, all_binary_prior);

cd(pathname)

save Tokyo2_20220912_sound_choice_separate ...
    long1_correct long1_error long1_sound long1_choice long1_prior ...
    long2_correct long2_error long2_sound long2_choice long2_prior ...
    long1_correct_p long1_error_p long1_sound_p long1_choice_p long1_prior_p ...
    long2_correct_p long2_error_p long2_sound_p long2_choice_p long2_prior_p ...
    short1_correct short1_error short1_sound short1_choice short1_prior ...
    short2_correct short2_error short2_sound short2_choice short2_prior ...
    short1_correct_p short1_error_p short1_sound_p short1_choice_p short1_prior_p ...
    short2_correct_p short2_error_p short2_sound_p short2_choice_p short2_prior_p
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stim1_index_correct,stim1_index_error,stim1_index_sound,stim1_index_choice,stim1_index_prior,...
          stim2_index_correct,stim2_index_error,stim2_index_sound,stim2_index_choice,stim2_index_prior, ...
          stim1_p_correct,stim1_p_error,stim1_p_sound,stim1_p_choice,stim1_p_prior,...
          stim2_p_correct,stim2_p_error,stim2_p_sound,stim2_p_choice,stim2_p_prior] = ...
    get_task_activ_separate(use_long,frame_sound,frame_choice_select,frame_spout,max_tif,...
    time_long_pre,time_long_post,time_long_pre2,time_long_post2,time_window,time_window_base, ...
    Correct_side, Chosen_side, correct_error, binary_prior)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Correct_side = Correct_side(use_long);
Chosen_side = Chosen_side(use_long);
correct_error = correct_error(use_long);
binary_prior = binary_prior(use_long);

low  = find(Correct_side == 0);
high = find(Correct_side == 1);
left = find(Chosen_side == 0);
right = find(Chosen_side == 1);
correct = find(correct_error == 1);
error = find(correct_error == 0);
prior1 = find(binary_prior == 1);
prior0 = find(binary_prior == 0);

trial_L_correct = intersect(low,correct);
trial_L_error = intersect(low,error);
trial_H_correct = intersect(high,correct);
trial_H_error = intersect(high,error);

time_long  = round((time_long_pre + time_long_post) ./ time_window);
time_long2 = round((time_long_pre2 + time_long_post2) ./ time_window);

frame_sound = frame_sound(use_long,:);
frame_choice_select = frame_choice_select(use_long,:);

frame_spout = frame_spout(use_long,:);
base_for_spout_on  = frame_spout(:,1)-time_window_base; %Before moving spout
base_for_spout_off = frame_spout(:,3)-time_window_base; %Before moving spout

for i = 1:time_long
    frame_sound_use_on(:,i) = frame_sound - time_long_pre + (i-1)*time_window;
    %frame_sound_use_off(:,i) = frame_sound - time_long_pre + i*time_window - 1;
end
for i = 1:time_long2
    frame_sound_use2_on(:,i) = frame_choice_select - time_long_pre2 + (i-1)*time_window;
    %frame_sound_use2_off(:,i) = frame_choice_select - time_long_pre2 + i*time_window - 1;
end

parfor file_count = 1:max_tif
    %temp_file = tif_name(file_count).name;
    %temp_file
    temp_file = sprintf('task_spike_stripe20210520_%d',file_count);
    temp_file
    
    %clear data
    data = load(temp_file); %spike_mark
    
    spike_mark = data.spike_mark;
    spike_multi(file_count) = data.single_multi;
    
    spike_frame1 = nan(length(frame_sound),time_long);
    spike_frame2 = nan(length(frame_sound),time_long2);
    for i = 1:length(use_long) %trial
        for j = 1:time_long
            if isnan(frame_sound_use_on(i,j)) == 0
                temp_frame = [frame_sound_use_on(i,j) : frame_sound_use_on(i,j)+time_window-1];
                spike_frame1(i,j) = mean(spike_mark(temp_frame));
            end
        end
        for j = 1:time_long2
            if isnan(frame_sound_use2_on(i,j)) == 0
                temp_frame = [frame_sound_use2_on(i,j) : frame_sound_use2_on(i,j)+time_window-1];
                spike_frame2(i,j) = mean(spike_mark(temp_frame));
            end
        end
    end
    for j = 1:time_long
        temp_H = mean(spike_frame1(trial_H_correct,j));
        temp_L = mean(spike_frame1(trial_L_correct,j));
        stim1_index_correct(file_count,j) = (temp_H - temp_L) ./ (temp_H + temp_L);
        stim1_p_correct(file_count,j) = ranksum(spike_frame1(trial_H_correct,j),spike_frame1(trial_L_correct,j));
        temp_H = mean(spike_frame1(trial_H_error,j));
        temp_L = mean(spike_frame1(trial_L_error,j));
        stim1_index_error(file_count,j) = (temp_H - temp_L) ./ (temp_H + temp_L);
        if ~isempty(trial_H_error) && ~isempty(trial_L_error)
            stim1_p_error(file_count,j) = ranksum(spike_frame1(trial_H_error,j),spike_frame1(trial_L_error,j));
        else
            stim1_p_error(file_count,j) = nan;
        end
        temp_H = mean(spike_frame1(high,j));
        temp_L = mean(spike_frame1(low,j));
        stim1_index_sound(file_count,j) = (temp_H - temp_L) ./ (temp_H + temp_L);
        stim1_p_sound(file_count,j) = ranksum(spike_frame1(high,j),spike_frame1(low,j));
        temp_H = mean(spike_frame1(right,j));
        temp_L = mean(spike_frame1(left,j));
        stim1_index_choice(file_count,j) = (temp_H - temp_L) ./ (temp_H + temp_L);
        stim1_p_choice(file_count,j) = ranksum(spike_frame1(right,j),spike_frame1(left,j));
        temp_H = mean(spike_frame1(prior1,j));
        temp_L = mean(spike_frame1(prior0,j));
        stim1_index_prior(file_count,j) = (temp_H - temp_L) ./ (temp_H + temp_L);
        if ~isempty(prior1) && ~isempty(prior0)
            stim1_p_prior(file_count,j) = ranksum(spike_frame1(prior1,j),spike_frame1(prior0,j));
        else
            stim1_p_prior(file_count,j) = nan;
        end
    end
    for j = 1:time_long2
        temp_H = mean(spike_frame2(trial_H_correct,j));
        temp_L = mean(spike_frame2(trial_L_correct,j));
        stim2_index_correct(file_count,j) = (temp_H - temp_L) ./ (temp_H + temp_L);
        stim2_p_correct(file_count,j) = ranksum(spike_frame2(trial_H_correct,j),spike_frame2(trial_L_correct,j));
        temp_H = mean(spike_frame2(trial_H_error,j));
        temp_L = mean(spike_frame2(trial_L_error,j));
        stim2_index_error(file_count,j) = (temp_H - temp_L) ./ (temp_H + temp_L);
        if ~isempty(trial_H_error) && ~isempty(trial_L_error)
            stim2_p_error(file_count,j) = ranksum(spike_frame2(trial_H_error,j),spike_frame2(trial_L_error,j));
        else
            stim2_p_error(file_count,j) = nan;
        end
        temp_H = mean(spike_frame2(high,j));
        temp_L = mean(spike_frame2(low,j));
        stim2_index_sound(file_count,j) = (temp_H - temp_L) ./ (temp_H + temp_L);
        stim2_p_sound(file_count,j) = ranksum(spike_frame2(high,j),spike_frame2(low,j));
        temp_H = mean(spike_frame2(right,j));
        temp_L = mean(spike_frame2(left,j));
        stim2_index_choice(file_count,j) = (temp_H - temp_L) ./ (temp_H + temp_L);
        stim2_p_choice(file_count,j) = ranksum(spike_frame2(right,j),spike_frame2(left,j));
        temp_H = mean(spike_frame2(prior1,j));
        temp_L = mean(spike_frame2(prior0,j));
        stim2_index_prior(file_count,j) = (temp_H - temp_L) ./ (temp_H + temp_L);
        if ~isempty(prior1) && ~isempty(prior0)
            stim2_p_prior(file_count,j) = ranksum(spike_frame2(prior1,j),spike_frame2(prior0,j));
        else
            stim2_p_prior(file_count,j) = nan;
        end
    end
end

return



