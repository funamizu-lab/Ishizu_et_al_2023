
%{
----------------------------------------------------------------------------
First_take number of frames in each tif files
Analyzing imaging data simply
At least for the correct rate
----------------------------------------------------------------------------
%}
function Task_kaiseki_tokyo_decode_single_neuron_20240110_save_auto(pathname)

%ROC in previous task parameters

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
prior = prior_value(:,2) ./ (prior_value(:,1) + prior_value(:,2));
posterior = posterior(:,2) ./ sum(posterior(:,1)+posterior(:,2));
binary_prior = zeros(length(prior),1);
binary_posterior = zeros(length(posterior),1);
temp = find(prior > 0.5);
binary_prior(temp) = 1;
temp = find(posterior > 0.5);
binary_posterior(temp) = 1;

%Setup the prior and posterior
if length(binary_prior) ~= length(use_trial)
    hoge
end
all_binary_prior = nan(length(Choice_trial),1);
all_binary_prior(use_trial) = binary_prior;


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
%roc_long = Long; %previous trial
%roc_short = Short; %previous trial
%roc_long = intersect(roc_long, use_trial);
%roc_short = intersect(roc_short, use_trial);
add_trial_L = 0;
for i = 1:length(use_long)
    temp = use_long(i);
    if temp == use_trial(1)
        roc_long(i,1) = use_trial(1)-1;
        add_trial_L = 1;
    else
        temp = find(use_trial < temp);
        roc_long(i,1) = use_trial(temp(end));
    end
end
add_trial_S = 0;
for i = 1:length(use_short)
    temp = use_short(i);
    if temp == use_trial(1)
        roc_short(i,1) = use_trial(1)-1;
        add_trial_S = 1;
    else
        temp = find(use_trial < temp);
        roc_short(i,1) = use_trial(temp(end));
    end
end
test = intersect(roc_long,use_trial);
if length(test)+add_trial_L ~= length(roc_long)
    [length(test), length(roc_long)]
    hoge
end
test = intersect(roc_short,use_trial);
if length(test)+add_trial_S ~= length(roc_short)
    [length(test), length(roc_short)]
    hoge
end
% [use_long-roc_long]
% hoge

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

[long_sound_ROC1,long_choice_ROC1,long_prior_ROC1, ...
 long_sound_ROC2,long_choice_ROC2,long_prior_ROC2] = ...
    ROC_task_activ_separate_20240110(use_long,frame_sound,frame_choice_select,frame_spout,max_tif,...
    time_long_pre,time_long_post,time_long_pre2,time_long_post2,time_window,time_window_base, ...
    Correct_side, Chosen_side, correct_error, all_binary_prior, roc_long);

[short_sound_ROC1,short_choice_ROC1,short_prior_ROC1, ...
 short_sound_ROC2,short_choice_ROC2,short_prior_ROC2] = ...
    ROC_task_activ_separate_20240110(use_short,frame_sound,frame_choice_select,frame_spout,max_tif,...
    time_short_pre,time_short_post,time_short_pre2,time_short_post2,time_window,time_window_base, ...
    Correct_side, Chosen_side, correct_error, all_binary_prior, roc_short);

cd(pathname)

save ROC_20240110_long_short ...
    long_sound_ROC1 long_choice_ROC1 long_prior_ROC1 ...
    long_sound_ROC2 long_choice_ROC2 long_prior_ROC2 ...
    short_sound_ROC1 short_choice_ROC1 short_prior_ROC1 ...
    short_sound_ROC2 short_choice_ROC2 short_prior_ROC2
    
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sound_ROC1,choice_ROC1,prior_ROC1, ...
          sound_ROC2,choice_ROC2,prior_ROC2] = ...
    ROC_task_activ_separate_20240110(use_long,frame_sound,frame_choice_select,frame_spout,max_tif,...
    time_long_pre,time_long_post,time_long_pre2,time_long_post2,time_window,time_window_base, ...
    Correct_side, Chosen_side, correct_error, binary_prior, roc_trial)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Correct_side = Correct_side(roc_trial);
Chosen_side = Chosen_side(roc_trial);
correct_error = correct_error(roc_trial);
binary_prior = binary_prior(roc_trial);

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

sound_ROC1 = nan(max_tif,time_long);
choice_ROC1 = nan(max_tif,time_long);
prior_ROC1 = nan(max_tif,time_long);
sound_ROC2 = nan(max_tif,time_long2);
choice_ROC2 = nan(max_tif,time_long2);
prior_ROC2 = nan(max_tif,time_long2);
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
        sound_ROC1(file_count,j) = make_auROC_20210620(spike_frame1(:,j), low, high);
        choice_ROC1(file_count,j) = make_auROC_20210620(spike_frame1(:,j), left, right);
        prior_ROC1(file_count,j) = make_auROC_20210620(spike_frame1(:,j), prior1, prior0);
    end
    for j = 1:time_long2
        sound_ROC2(file_count,j) = make_auROC_20210620(spike_frame2(:,j), low, high);
        choice_ROC2(file_count,j) = make_auROC_20210620(spike_frame2(:,j), left, right);
        prior_ROC2(file_count,j) = make_auROC_20210620(spike_frame2(:,j), prior1, prior0);
    end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mean_prob = make_auROC_20210620(activity, low, high) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

low_activity = sort(activity(low));
high_activity = sort(activity(high));
low_activity = round(low_activity,3);
high_activity = round(high_activity,3);

if isempty(low) || isempty(high)
    mean_prob = nan;
else
    for i = 1:length(low_activity)
        temp = find(high_activity > low_activity(i));
        temp05 = find(high_activity == low_activity(i));
        prob(i) = (length(temp)+length(temp05)/2) ./ length(high_activity);
    end
    mean_prob = mean(prob);
    if mean_prob < 0.5
        mean_prob = 1 - mean_prob;
    end
end

max_activity = max(activity);
min_activity = min(activity);
if mean_prob > 0.9
    hist_low = histcounts(low_activity,[min_activity:0.001:max_activity]);
    hist_high = histcounts(high_activity,[min_activity:0.001:max_activity]);
    %[min_activity,max_activity]
% %     hist_low
% %     hist_high
% %     prob
%     figure
%     plot(hist_low,'b')
%     hold on
%     plot(hist_high,'r')
%     %hoge
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [spike_count,p,spike_trace] = get_sound_response3_count_only(spike_mark, spike_filter, frame_sound, pre_frame, post_frame, pre_frame2, post_frame2)
function [spike_count,p] = get_sound_response3_count_only(spike_mark, frame_sound, pre_frame, post_frame)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spike_count = nan(length(frame_sound),2);
%spike_trace = nan(length(frame_sound),pre_frame2+post_frame2);

for i = 1:length(frame_sound)
    
    if isnan(frame_sound(i)) == 0
        temp_pre  = [frame_sound(i)-pre_frame : frame_sound(i)-1];
        temp_post = [frame_sound(i) : frame_sound(i)+post_frame-1];
        %temp_all = [frame_sound(i)-pre_frame2 : frame_sound(i)+post_frame2-1];
    
        temp_pre = spike_mark(temp_pre);
        temp_post = spike_mark(temp_post);
        spike_count(i,:) = [sum(temp_pre), sum(temp_post)];
    else
        spike_count(i,:) = [nan, nan];
    end
    %spike_trace(i,:) = spike_filter(temp_all);
end

p = signrank(spike_count(:,1),spike_count(:,2));

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [frame_behave] = thre_detection_frame(behave_time, blue_scan)    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    temp = find(isnan(behave_time) == 1);
    if isempty(temp),
        for i = 1:length(behave_time),
            temp_time = behave_time(i);
            %frame_behave(i) = find(blue_scan >= temp_time, 1);
            temp_behave = find(blue_scan >= temp_time, 1);
            %Rotary encoder does not stop by task!!
            if length(temp_behave) == 0,
                temp_time
                blue_scan(length(blue_scan))
            end
            frame_behave(i) = temp_behave;
        end
    else
        frame_behave = nan;
    end

    return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [length_scan, scan_time, first_scan] = thre_detection(temp_scan, t, thre_scan)    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %temp_scan = ch(i,:);
    temp = find(temp_scan > thre_scan); %above thre
    temp_plus = [-1, temp]; %Add the first value
    temp1 = temp_plus(1:length(temp_plus)-1);
    temp2 = temp_plus(2:length(temp_plus));
    temp_sabun = temp2 - temp1;
    temp_sabun1 = find(temp_sabun > 1);
    temp_sabun2 = find(temp_sabun == 1);
    temp_sabun2 = temp_sabun2 - 1;
    temp_use = intersect(temp_sabun1, temp_sabun2);
    
    if ~isempty(temp_use),
        temp_use = temp(temp_use);
    %     length_scan(i) = length(temp_use);
    %     scan_time(i).matrix = t(temp_use);
    %     first_scan(i) = t(temp_use(1));
        length_scan = length(temp_use);
        scan_time = t(temp_use);
        first_scan = t(temp_use(1));
    else
        length_scan = nan;
        scan_time = nan;
        first_scan = nan;
    end
    
    return
