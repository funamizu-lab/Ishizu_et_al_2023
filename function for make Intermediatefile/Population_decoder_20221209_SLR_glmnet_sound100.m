
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
function Population_decoder_20221209_SLR_glmnet_sound100(pathname)

switch nargin
    case 0
        pathname   = pwd;
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

temp = dir('RL_20220818*'); %RL_20220118 frame
if length(temp) ~= 1
    temp
    hoge
end
load(temp.name);
%'ave_likeli' ,'BIC','log_likeli','para_max','N_trial'

temp = dir('task_frame*'); %Task frame
if length(temp) ~= 1
    temp
    hoge
end
load(temp.name);

temp = dir('Bpod*'); %Bpod
if length(temp) ~= 1
    temp
    hoge
end
load(temp.name);
Bpod_file = temp.name;

cd(spike_dir)
tif_name = dir('task_spike*.mat'); %get all the tif files
length_tif = length(tif_name);
cd(pathname)

temp = dir('depth_spike_20220517*');
if length(temp) == 1
    load(temp.name);
    %spike_depth def_depth length_neuron
    if length_tif ~= length_neuron
        [length_tif length_neuron]
        hoge
    end
    depth_neuron = find(spike_depth <= def_depth);
elseif length(temp) == 0
    depth_neuron = [1:length_tif]; %Use all the neurons
else
    temp
    hoge
end

%Get parameter of RL model
% [prior,posterior,ave_likeli,likelihood,Q,Choice_prob,para] = ...
%     Dual_RL_model_block1_20220118_para_determined(Bpod_file,para_max(2,:),N_trial);
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

[minD_trial,maxD_trial,Choice_trial,tone_evidence,trial_evidence,use_trial2,use_trial3,use_trial_all,...
    low,high,correct,error,flip_tone,number_use_trial,...
    binary_tone,right_trial_all,number_trial_all,right_trial,number_trial] ...
 = Dual_get_basic_task_structure_20210204_2(Bpod_file);

%Chech the number of trials 
if length(all_trial_time) ~= length(Outcome)
    [length(all_trial_time), length(Outcome)]
    hoge
end

%Get task parameter
%Choice_trial = find(Outcome == 1 | Outcome == 2);
%low  = find(Correct_side == 0);
%high = find(Correct_side == 1);
Reward = Correct_side == Chosen_side;
Reward = double(Reward);
left  = find(Chosen_side == 0);
right = find(Chosen_side == 1);
stim_length = unique(StimDuration);
long  = find(StimDuration == stim_length(2));
short = find(StimDuration == stim_length(1));

% %Get tone evidence
% temp_evi = unique(EvidenceStrength);
% temp_evi_low  = 0.5 - temp_evi/2;
% temp_evi_high = 0.5 + temp_evi/2;
% temp_evi_all = [temp_evi_low', temp_evi_high'];
% tone_evidence = sort(temp_evi_all);

% %Put tone evidence in all trials;
% trial_evidence = zeros(length(Outcome),1);
% for i = 1:length(temp_evi),
%     temp = find(EvidenceStrength == temp_evi(i));
%     temp_left  = intersect(temp,low);
%     temp_right = intersect(temp,high);
%     trial_evidence(temp_left)  = temp_evi_low(i);
%     trial_evidence(temp_right) = temp_evi_high(i);
% end
for i = 1:length(tone_evidence)
    temp = find(trial_evidence == tone_evidence(i));
    evi_trial(i).matrix = temp;
end

% %TrialBlock
% for i = 1:max(TrialBlock)
%     block(i).matrix = find(TrialBlock == i);
% end
% block2 = [];
% block3 = [];
% for i = 2:max(TrialBlock)
%     if rem(i,2) == 0
%         block2 = [block2; block(i).matrix];
%     else
%         block3 = [block3; block(i).matrix];
%     end
% end
% block2 = sort(block2);
% block3 = sort(block3);
% %block2 = sort([block(2).matrix; block(4).matrix]);
% %block3 = sort([block(3).matrix; block(5).matrix]);
% block2 = intersect(block2, Choice_trial);
% block3 = intersect(block3, Choice_trial);
% use_trial = sort([block2;block3]);

block2 = use_trial2;
block3 = use_trial3;
use_trial = use_trial_all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Setup the prior and posterior
if length(binary_prior) ~= length(use_trial)
    hoge
end
if length(binary_posterior) ~= length(use_trial)
    hoge
end
StimDuration23 = StimDuration(use_trial);
long23_prior  = find(StimDuration23 == stim_length(2));
short23_prior = find(StimDuration23 == stim_length(1));

binary_prior_long23 = binary_prior(long23_prior);
binary_prior_short23 = binary_prior(short23_prior);
binary_posterior_long23 = binary_posterior(long23_prior);
binary_posterior_short23 = binary_posterior(short23_prior);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% low23  = intersect(low,use_trial);
% high23 = intersect(high,use_trial);
% left23  = intersect(left,use_trial);
% right23 = intersect(right,use_trial);
% use_frame = frame_sound(use_trial);
short23 = intersect(short,use_trial);
long23 = intersect(long,use_trial);

cd(spike_dir)

tif_name = dir('task_spike*.mat'); %get all the tif files
length_tif = length(tif_name);
%max_tif = length_tif;
max_tif = length(depth_neuron);

%Short trials
LowHigh_s = Correct_side(short23);
LeftRight_s = Chosen_side(short23);
SoundEvi_s = trial_evidence(short23);
Reward_s = Reward(short23);
block_s = TrialBlock(short23);
s_low  = find(LowHigh_s == 0);
s_high = find(LowHigh_s == 1);
s_left  = find(LeftRight_s == 0);
s_right = find(LeftRight_s == 1);
s_low_left = intersect(s_low,s_left);
s_low_right = intersect(s_low,s_right);
s_high_left = intersect(s_high,s_left);
s_high_right = intersect(s_high,s_right);
%This is for choice decoding
for i = 1:length(tone_evidence)
    temp = find(SoundEvi_s == tone_evidence(i));
    %evi_trial(i).matrix = temp;
    s_evi_left(i).matrix = intersect(temp,s_left);
    s_evi_right(i).matrix = intersect(temp,s_right);
end

%Long trials
LowHigh_l = Correct_side(long23);
LeftRight_l = Chosen_side(long23);
SoundEvi_l = trial_evidence(long23);
Reward_l = Reward(long23);
block_l = TrialBlock(long23);
l_low  = find(LowHigh_l == 0);
l_high = find(LowHigh_l == 1);
l_left  = find(LeftRight_l == 0);
l_right = find(LeftRight_l == 1);
l_low_left = intersect(l_low,l_left);
l_low_right = intersect(l_low,l_right);
l_high_left = intersect(l_high,l_left);
l_high_right = intersect(l_high,l_right);
%This is for choice decoding
for i = 1:length(tone_evidence)
    temp = find(SoundEvi_l == tone_evidence(i));
    %evi_trial(i).matrix = temp;
    l_evi_left(i).matrix = intersect(temp,l_left);
    l_evi_right(i).matrix = intersect(temp,l_right);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make time for short sound and long sound trials
time_window = 200;
pre_sound  = 2000;
post_sound = 3000;
% time_long_pre2  = 500;
% time_long_post2 = 2500;
window_length = round((post_sound + pre_sound) ./ time_window);
%all_sound_short = nan(length(short23),max_tif);
for j = 1:window_length
    all_sound_short(j).matrix = nan(length(short23),max_tif);
    all_sound_long(j).matrix = nan(length(long23),max_tif);
end
    
for file_count = 1:max_tif
    [file_count, max_tif]
    %temp_file = tif_name(file_count).name;
    %temp_file
    %temp_file = sprintf('task_spike_stripe20210520_%d',sig_neuron(file_count));
    %temp_file = sprintf('task_spike_stripe20210520_%d', file_count);
    temp_file = sprintf('task_spike_stripe20210520_%d',depth_neuron(file_count));
    temp_file
    
    %clear data
    data = load(temp_file); %spike_mark
    spike_mark = data.spike_mark;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Short sound trials
    sound_frame_short = nan(length(short23),window_length);
    for i = 1:length(short23)
        for j = 1:window_length
            temp_frame = frame_sound(short23(i));
            temp1 = temp_frame + time_window*(j-1) - pre_sound;
            temp2 = temp_frame + time_window*(j) - pre_sound - 1;
            %temp_frame = [temp_frame : temp_frame + time_window-1];
            temp_frame = spike_mark([temp1:temp2]);
            sound_frame_short(i,j) = mean(temp_frame);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Long sound trials
    sound_frame_long = nan(length(long23),window_length);
    for i = 1:length(long23)
        for j = 1:window_length
            temp_frame = frame_sound(long23(i));
            temp1 = temp_frame + time_window*(j-1) - pre_sound;
            temp2 = temp_frame + time_window*(j) - pre_sound - 1;
            %temp_frame = [temp_frame : temp_frame + time_window-1];
            temp_frame = spike_mark([temp1:temp2]);
            sound_frame_long(i,j) = mean(temp_frame);
        end
    end    
    
    for j = 1:window_length
        all_sound_short(j).matrix(:,file_count) = sound_frame_short(:,j);
        all_sound_long(j).matrix(:,file_count) = sound_frame_long(:,j);
    end
end
cd(pathname)

%each sound make SLR
CV_para = 10;
CV_repeat = 1000;
use_trial = 180;

use_neuron = 100;
selected_window = [10 11 12 15 17];

%each sound make SLR
if length(LowHigh_s) >= use_trial && max_tif >= use_neuron
%if max_tif >= use_neuron
    %Check whether the binary_priors are OK
    test0 = find(LowHigh_s(1:use_trial) == 0);
    test1 = find(LowHigh_s(1:use_trial) == 1);
    if length(test0) < 2 ||  length(test1) < 2
        short_correct_rate = [];
        short_ave_likeli = [];
        short_dp = [];
    else
        %for j = 1:window_length
        for j = 1:length(selected_window)
            %[j,window_length]
            %[short_correct_rate(:,j),short_ave_likeli(:,j)] = SparseLogisticRegression_CV_20220404_glmnet_sound2(all_sound_short(selected_window(j)).matrix,LowHigh_s,LeftRight_s,SoundEvi_s,Reward_s,block_s,CV_para,CV_repeat,use_neuron,use_trial);
            [short_correct_rate(:,j),short_ave_likeli(:,j),short_dp(:,j)] = SparseLogisticRegression_CV_20221201_glmnet_sound2(all_sound_short(selected_window(j)).matrix,LowHigh_s,LeftRight_s,SoundEvi_s,Reward_s,block_s,CV_para,CV_repeat,use_neuron,use_trial);
        end
    end
else
    short_correct_rate = [];
    short_ave_likeli = [];
    short_dp = [];
end
if length(LowHigh_l) >= use_trial && max_tif >= use_neuron
%if max_tif >= use_neuron
    test0 = find(LowHigh_l(1:use_trial) == 0);
    test1 = find(LowHigh_l(1:use_trial) == 1);
    if length(test0) < 2 ||  length(test1) < 2
        long_correct_rate = [];
        long_ave_likeli = [];
        long_dp = [];
    else
        %for j = 1:window_length
        for j = 1:length(selected_window)
            %[j,window_length]
            %[long_correct_rate(:,j),long_ave_likeli(:,j)] = SparseLogisticRegression_CV_20220404_glmnet_sound2(all_sound_long(selected_window(j)).matrix,LowHigh_l,LeftRight_l,SoundEvi_l,Reward_l,block_l,CV_para,CV_repeat,use_neuron,use_trial);
            [long_correct_rate(:,j),long_ave_likeli(:,j),long_dp(:,j)] = SparseLogisticRegression_CV_20221201_glmnet_sound2(all_sound_long(selected_window(j)).matrix,LowHigh_l,LeftRight_l,SoundEvi_l,Reward_l,block_l,CV_para,CV_repeat,use_neuron,use_trial);
        end
    end
else
    long_correct_rate = [];
    long_ave_likeli = [];   
    long_dp = [];
end
length_trial = [length(LowHigh_s),length(LowHigh_l)];

save SLR100_20221209_glmnet_sound short_correct_rate short_ave_likeli long_correct_rate long_ave_likeli use_neuron use_trial length_trial selected_window long_dp short_dp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
function [ROC_sound,p_sound,temp_shuffle_ROC] = get_ROC_shuffle(sound_frame, low, high, shuffle_time)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ROC_sound = make_auROC_20210620(sound_frame, low, high); %left choice

temp_shuffle_ROC  = [];
parfor i = 1:shuffle_time
%for i = 1:shuffle_time
    shuffle_frame = sound_frame(randperm(length(sound_frame)));
    temp_shuffle_ROC(1,i) = make_auROC_20210620(shuffle_frame, low, high);
end
temp_sound = find(temp_shuffle_ROC > ROC_sound); %how many random beyound the ROC
p_sound = length(temp_sound) ./ shuffle_time;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sound_choice_index = get_sound_choice_index(spike_sound,trial_low,trial_high,trial_left,trial_right)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    spike_low = mean(spike_sound(trial_low,2));
    spike_high = mean(spike_sound(trial_high,2));
    spike_left = mean(spike_sound(trial_left,2));
    spike_right = mean(spike_sound(trial_right,2));
    sound_choice_index(1) = (spike_low - spike_high) ./ (spike_low + spike_high);
    sound_choice_index(2) = (spike_left - spike_right) ./ (spike_left + spike_right);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spike_count,p,spike_trace] = get_sound_response3(spike_mark, spike_filter, frame_sound, pre_frame, post_frame, pre_frame2, post_frame2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spike_count = nan(length(frame_sound),2);
spike_trace = nan(length(frame_sound),pre_frame2+post_frame2);

for i = 1:length(frame_sound)
    
    if isnan(frame_sound(i)) == 0
        temp_pre  = [frame_sound(i)-pre_frame : frame_sound(i)-1];
        temp_post = [frame_sound(i) : frame_sound(i)+post_frame-1];
        temp_all = [frame_sound(i)-pre_frame2 : frame_sound(i)+post_frame2-1];

        temp_pre = spike_mark(temp_pre);
        temp_post = spike_mark(temp_post);
        spike_count(i,:) = [sum(temp_pre), sum(temp_post)];
    
        spike_trace(i,:) = spike_filter(temp_all);
    else
        spike_count(i,:) = nan(1,2);
        spike_trace(i,:) = nan(1,post_frame2+pre_frame2);
    end        
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
