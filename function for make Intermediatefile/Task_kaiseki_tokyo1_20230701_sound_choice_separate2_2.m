
%{
----------------------------------------------------------------------------
First_take number of frames in each tif files
Analyzing imaging data simply
At least for the correct rate
----------------------------------------------------------------------------
%}
function Task_kaiseki_tokyo1_20230701_sound_choice_separate2_2(pathname)

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get parameter of RL model
[prior,posterior,ave_likeli,likelihood,Long_posterior,Short_posterior,...
    prior_value,long_value,short_value,para,conf_Q] = ...
    Dual_RL_model_block1_20220703_para_determined(Bpod_file,para_max(3,:),N_trial);
prior = prior_value(:,2) ./ (prior_value(:,1) + prior_value(:,2));
posterior = posterior(:,2) ./ (posterior(:,1)+posterior(:,2));
conf_Q = conf_Q(:,2) ./ (conf_Q(:,1)+conf_Q(:,2)); %relative
binary_prior = zeros(length(prior),1);
binary_posterior = zeros(length(posterior),1);
binary_conf = zeros(length(conf_Q),1);
temp = find(prior > 0.5);
binary_prior(temp) = 1;
temp = find(posterior > 0.5);
binary_posterior(temp) = 1;
temp = find(conf_Q > 0.5);
binary_conf(temp) = 1;
%prior posterior conf_Q binary_prior binary_posterior binary_conf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Setup the prior and posterior
if length(binary_prior) ~= length(use_trial)
    hoge
end
% all_binary_prior = nan(length(Choice_trial),1);
% all_binary_prior(use_trial) = binary_prior;

stim_length = unique(StimDuration);
%Focus on the use_trial
StimDuration = StimDuration(use_trial);
Chosen_side = Chosen_side(use_trial);
Correct_side = Correct_side(use_trial);
frame_choice = frame_choice(use_trial,:);

%Get task parameter
low  = find(Correct_side == 0);
high = find(Correct_side == 1);
Long  = find(StimDuration == stim_length(2));
Short = find(StimDuration == stim_length(1));
left = find(Chosen_side == 0);
right = find(Chosen_side == 1);
correct_error = Correct_side == Chosen_side;
correct = find(correct_error == 1);
error = find(correct_error == 0);

%Make the choice timing frames
frame_sound = frame_sound(use_trial);
frame_spout = frame_spout(use_trial,:);

frame_choice_select = nan(length(frame_choice),1);
frame_choice_select(left) = frame_choice(left,1);
frame_choice_select(right) = frame_choice(right,2);

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

%Just focus on the 3 phases, before, initial and end of sound
%save the raw data with and the Q-values
long_use_timing = [15, 17, 25]; %before, first and end of sound
short_use_timing = [15, 17]; %before, first and end of sound

tif_name = dir('task_spike_stripe*.mat'); %get all the tif files
length_tif = length(tif_name);
max_tif = length_tif;

[long_spike, long_prior, long_posterior, long_conf_Q, long_binary_prior, long_binary_posterior, long_binary_conf,long_rotation,wheel_good(1)] = ...
    get_task_activ_separate_20230701(long_use_timing, Long,frame_sound,frame_choice_select,frame_spout,max_tif,...
    time_long_pre,time_long_post,time_long_pre2,time_long_post2,time_window,time_window_base, ...
    prior, posterior, conf_Q, binary_prior, binary_posterior, binary_conf, ave_velocity);

[short_spike, short_prior, short_posterior, short_conf_Q, short_binary_prior, short_binary_posterior, short_binary_conf,short_rotation,wheel_good(2)] = ...
    get_task_activ_separate_20230701(short_use_timing, Short,frame_sound,frame_choice_select,frame_spout,max_tif,...
    time_short_pre,time_short_post,time_short_pre2,time_short_post2,time_window,time_window_base, ...
    prior, posterior, conf_Q, binary_prior, binary_posterior, binary_conf, ave_velocity);

cd(pathname)

save Tokyo2_20230701_prior_conf_activity ...
    long_spike long_prior long_posterior long_conf_Q long_binary_prior long_binary_posterior long_binary_conf ...
    short_spike short_prior short_posterior short_conf_Q short_binary_prior short_binary_posterior short_binary_conf ...
    long_rotation short_rotation wheel_good

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spike_frame1, prior, posterior, conf_Q, binary_prior, binary_posterior, binary_conf,rotation1,wheel_good] = ...
    get_task_activ_separate_20230701(long_use_timing, use_long,frame_sound,frame_choice_select,frame_spout,max_tif,...
    time_long_pre,time_long_post,time_long_pre2,time_long_post2,time_window,time_window_base, ...
    prior, posterior, conf_Q, binary_prior, binary_posterior, binary_conf, ave_velocity)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_long  = round((time_long_pre + time_long_post) ./ time_window);
time_long2 = round((time_long_pre2 + time_long_post2) ./ time_window);

frame_sound = frame_sound(use_long,:);
frame_choice_select = frame_choice_select(use_long,:);

frame_spout = frame_spout(use_long,:);
% base_for_spout_on  = frame_spout(:,1)-time_window_base; %Before moving spout
% base_for_spout_off = frame_spout(:,3)-time_window_base; %Before moving spout

prior = prior(use_long);
posterior = posterior(use_long);
conf_Q = conf_Q(use_long);
binary_prior = binary_prior(use_long);
binary_posterior = binary_posterior(use_long);
binary_conf = binary_conf(use_long);

for i = 1:time_long
    frame_sound_use_on(:,i) = frame_sound - time_long_pre + (i-1)*time_window;
    %frame_sound_use_off(:,i) = frame_sound - time_long_pre + i*time_window - 1;
end
for i = 1:time_long2
    frame_sound_use2_on(:,i) = frame_choice_select - time_long_pre2 + (i-1)*time_window;
    %frame_sound_use2_off(:,i) = frame_choice_select - time_long_pre2 + i*time_window - 1;
end

ave_velocity = double(ave_velocity);

rotation1 = zeros(length(use_long),length(long_use_timing));
wheel_good = 1;
%if length(ave_velocity) ~= 1
if length(ave_velocity) > frame_sound(end)
    for i = 1:length(use_long) %trial
        for j = 1:length(long_use_timing) %time
            use_j = long_use_timing(j);
            
            if isnan(frame_sound_use_on(i,use_j)) == 0
                temp_frame = [frame_sound_use_on(i,use_j) : frame_sound_use_on(i,use_j)+time_window-1];
                rotation1(i,j) = mean(ave_velocity(temp_frame));
            end
        end
    end
else
    wheel_good = 0;
end

for i = 1:length(long_use_timing)
    spike_frame1(i).matrix = nan(length(frame_sound),max_tif);
    %spike_frame2 = nan(length(frame_sound),time_long2);
end

    
for file_count = 1:max_tif
    %temp_file = tif_name(file_count).name;
    %temp_file
    temp_file = sprintf('task_spike_stripe20210520_%d',file_count);
    temp_file
    
    %clear data
    data = load(temp_file); %spike_mark
    
    spike_mark = data.spike_mark;
    spike_multi(file_count) = data.single_multi;
    
    %Normalize the spike_mark
    mean_spike = mean(spike_mark);
    std_spike = std(spike_mark);
    if std_spike ~= 0
        spike_mark = (spike_mark - mean_spike) ./ std_spike; %normalized the activity
    else
        spike_mark = zeros(1,length(spike_mark));
    end
    
    temp_spike_frame1 = nan(length(frame_sound),length(long_use_timing));
    
    for i = 1:length(use_long) %trial
        for j = 1:length(long_use_timing) %trial
            use_j = long_use_timing(j);
            
            if isnan(frame_sound_use_on(i,use_j)) == 0
                temp_frame = [frame_sound_use_on(i,use_j) : frame_sound_use_on(i,use_j)+time_window-1];
                temp_spike_frame1(i,j) = mean(spike_mark(temp_frame));
            end
        end
%         for j = 1:time_long2
%             if isnan(frame_sound_use2_on(i,j)) == 0
%                 temp_frame = [frame_sound_use2_on(i,j) : frame_sound_use2_on(i,j)+time_window-1];
%                 spike_frame2(i,j) = mean(spike_mark(temp_frame));
%             end
%         end
    end
    for j = 1:length(long_use_timing)
        spike_frame1(j).matrix(:,file_count) = temp_spike_frame1(:,j);
    end
end

return





