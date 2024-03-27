
%{
----------------------------------------------------------------------------
First_take number of frames in each tif files
Analyzing imaging data simply
At least for the correct rate
----------------------------------------------------------------------------
%}
function Task_kaiseki_tokyo1_20220221_make_ave_trace(pathname)

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

% [filename2, pathname2,findex]=uigetfile('*.*','bpod file');
% filename2 = [pathname2,filename2];
% load(filename2)
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

%Get task parameter
Choice_trial = find(Outcome == 1 | Outcome == 2);
Reward = Correct_side == Chosen_side;
low  = find(Correct_side == 0);
high = find(Correct_side == 1);
correct = find(Reward == 1);
error = find(Reward == 0);
stim_length = unique(StimDuration);
Long  = find(StimDuration == stim_length(2));
Short = find(StimDuration == stim_length(1));
Long = intersect(Long, use_trial);
Short = intersect(Short, use_trial);

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

%This is for choice decoding
for i = 1:length(tone_evidence)
    temp = find(trial_evidence == tone_evidence(i));
    evi_stim(i).matrix = intersect(temp,use_trial);
    s_evi_stim(i).matrix = intersect(temp, Short);
    l_evi_stim(i).matrix = intersect(temp, Long);
    
    %Only with correct trials
    s_evi_correct(i).matrix = intersect(s_evi_stim(i).matrix, correct);
    l_evi_correct(i).matrix = intersect(l_evi_stim(i).matrix, correct);
end
%Block dependent correct and error check
%Block_R
LongR_detail(1).matrix = get_trial_block_type(low, correct, Long_R);
LongR_detail(2).matrix = get_trial_block_type(low, error, Long_R);
LongR_detail(3).matrix = get_trial_block_type(high, correct, Long_R);
LongR_detail(4).matrix = get_trial_block_type(high, error, Long_R);
LongL_detail(1).matrix = get_trial_block_type(low, correct, Long_L);
LongL_detail(2).matrix = get_trial_block_type(low, error, Long_L);
LongL_detail(3).matrix = get_trial_block_type(high, correct, Long_L);
LongL_detail(4).matrix = get_trial_block_type(high, error, Long_L);

ShortR_detail(1).matrix = get_trial_block_type(low, correct, Short_R);
ShortR_detail(2).matrix = get_trial_block_type(low, error, Short_R);
ShortR_detail(3).matrix = get_trial_block_type(high, correct, Short_R);
ShortR_detail(4).matrix = get_trial_block_type(high, error, Short_R);
ShortL_detail(1).matrix = get_trial_block_type(low, correct, Short_L);
ShortL_detail(2).matrix = get_trial_block_type(low, error, Short_L);
ShortL_detail(3).matrix = get_trial_block_type(high, correct, Short_L);
ShortL_detail(4).matrix = get_trial_block_type(high, error, Short_L);

%Get the sound frame
% pre_frame = 500; %1sec
% post_frame = 500; %1sec
pre_frame = 200; %1sec
post_frame = 200; %1sec
pre_frame2 = 2000; %1sec
post_frame2 = 4000; %1sec

cd(spike_dir);

tif_name = dir('task_spike_stripe*.mat'); %get all the tif files
length_tif = length(tif_name);
max_tif = length_tif;

use_first_trial = use_trial(1);
use_end_trial = use_trial(end);
frame_start = frame_sound(use_first_trial) - pre_frame2;
frame_end = frame_sound(use_end_trial) + post_frame2;
nan_trace = nan(1,pre_frame2+post_frame2);

for i = 1:max_tif
%for i = 3:3
%    [temp(i), neuron_base(temp(i))]
    [i,max_tif]
    temp_file = sprintf('task_spike_stripe20210520_%d',i);
    
    clear data
    data = load(temp_file); %spike_mark
    spike_mark = data.spike_mark;
    
    [~,~,spike_trace] = get_sound_response3(data.spike_mark, data.spike_filter, frame_sound, pre_frame, post_frame, pre_frame2, post_frame2);
    
    spike_mark = spike_mark(frame_start:frame_end);
    mean_spike(i,1) = mean(spike_mark);
    std_spike(i,1) = std(spike_mark);
    if std_spike(i) ~= 0
        norm_spike_trace = (spike_trace - mean_spike(i)) ./ std_spike(i);
    else
        norm_spike_trace = spike_trace;
    end

    %Use Normalized activity
    for j = 1:length(tone_evidence)
        evi_spike(j).matrix(i,:) = make_ave_trace_trial(norm_spike_trace, evi_stim(j).matrix, nan_trace);
        s_evi_spike(j).matrix(i,:) = make_ave_trace_trial(norm_spike_trace, s_evi_stim(j).matrix, nan_trace);
        l_evi_spike(j).matrix(i,:) = make_ave_trace_trial(norm_spike_trace, l_evi_stim(j).matrix, nan_trace);
        
        %Only with correct trials
        s_evi_correct_spike(j).matrix(i,:) = make_ave_trace_trial(norm_spike_trace, s_evi_correct(j).matrix, nan_trace);
        l_evi_correct_spike(j).matrix(i,:) = make_ave_trace_trial(norm_spike_trace, l_evi_correct(j).matrix, nan_trace);
    end
    for j = 1:4
        LongR_spike(j).matrix(i,:) = make_ave_trace_trial(norm_spike_trace, LongR_detail(j).matrix, nan_trace);
        LongL_spike(j).matrix(i,:) = make_ave_trace_trial(norm_spike_trace, LongL_detail(j).matrix, nan_trace);
        ShortR_spike(j).matrix(i,:) = make_ave_trace_trial(norm_spike_trace, ShortR_detail(j).matrix, nan_trace);
        ShortL_spike(j).matrix(i,:) = make_ave_trace_trial(norm_spike_trace, ShortL_detail(j).matrix, nan_trace);
    end
end

cd(pathname)

save Tokyo1_20220221_make_ave_trace ...
    evi_spike s_evi_spike l_evi_spike s_evi_correct_spike l_evi_correct_spike ...
    LongR_spike LongL_spike ShortR_spike ShortL_spike


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
function ave_trace = make_ave_trace_trial(norm_spike_trace, use_trial, nan_trace)    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

if length(use_trial) ~= 0
    ave_trace = mean(norm_spike_trace(use_trial,:));
else
    ave_trace = nan_trace;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial = get_trial_block_type(low, correct, Long_R)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trial = intersect(low, correct);
trial = intersect(trial, Long_R);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spike_count,p,spike_trace] = get_sound_response3(spike_mark, spike_filter, frame_sound, pre_frame, post_frame, pre_frame2, post_frame2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spike_count = nan(length(frame_sound),2);
spike_trace = nan(length(frame_sound),pre_frame2+post_frame2);

for i = 1:length(frame_sound)
    temp_pre  = [frame_sound(i)-pre_frame : frame_sound(i)-1];
    temp_post = [frame_sound(i) : frame_sound(i)+post_frame-1];
    temp_all = [frame_sound(i)-pre_frame2 : frame_sound(i)+post_frame2-1];
    
    temp_pre = spike_mark(temp_pre);
    temp_post = spike_mark(temp_post);
    spike_count(i,:) = [sum(temp_pre), sum(temp_post)];
    
    spike_trace(i,:) = spike_filter(temp_all);
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
