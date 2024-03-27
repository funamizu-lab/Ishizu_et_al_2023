
%{
----------------------------------------------------------------------------
Determine the time window for task relevant neurons
%Before start of task
%Start of task
%Sound on
%Sound off
%Before choice
%After choice (0sec)
%After choice (1sec)
%After choice (2sec)
%p = 0.001

----------------------------------------------------------------------------
%}
function Task_kaiseki_tokyo1_20220201_2(pathname)

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
temp = dir('task_frame_tokyo_ephys_20220210*');
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
low  = find(Correct_side == 0);
high = find(Correct_side == 1);
left  = find(Chosen_side == 0);
right = find(Chosen_side == 1);
stim_length = unique(StimDuration);
Long  = find(StimDuration == stim_length(2));
Short = find(StimDuration == stim_length(1));

%Get tone evidence
temp_evi = unique(EvidenceStrength);
temp_evi_low  = 0.5 - temp_evi/2;
temp_evi_high = 0.5 + temp_evi/2;
temp_evi_all = [temp_evi_low', temp_evi_high'];
tone_evidence = sort(temp_evi_all);

%Put tone evidence in all trials;
trial_evidence = zeros(length(Outcome),1);
for i = 1:length(temp_evi),
    temp = find(EvidenceStrength == temp_evi(i));
    temp_left  = intersect(temp,low);
    temp_right = intersect(temp,high);
    trial_evidence(temp_left)  = temp_evi_low(i);
    trial_evidence(temp_right) = temp_evi_high(i);
end
for i = 1:length(tone_evidence)
    temp = find(trial_evidence == tone_evidence(i));
    evi_trial(i).matrix = temp;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Use to determine the threshold
%Start of task
%Sound on
%Sound off
%Before choice
%After choice (0sec)
%After choice (1sec)
%After choice (2sec)
%p = 0.001

%frame_sound_off
frame_sound_off = frame_sound;
frame_sound_off(Long) = frame_sound_off(Long) + 1000; %Add 1000 ms
frame_sound_off(Short) = frame_sound_off(Short) + 200; %Add 200 ms 
%frame_choice_select
frame_choice_select = nan(length(frame_choice),1);
frame_choice_select(left) = frame_choice(left,1);
frame_choice_select(right) = frame_choice(right,2);
temp = frame_choice_select(Choice_trial);
if max(isnan(temp)) == 1
    frame_choice_select(Choice_trial)
    test_nan = isnan(temp);
    %sum(isnan(temp))
    test_nan = find(test_nan == 1);
    Choice_trial(test_nan)
    
    frame_choice(test_nan,:)
    Chosen_side(test_nan,:)
    Outcome(test_nan)
    trial_lick_left(test_nan).matrix
    trial_lick_right(test_nan).matrix
    hoge
end

% %CHANGE frame_sound_off
% frame_sound_off = frame_choice_select - 500;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%frame_spout: move on off and move on off
time_window = 200; %ms
base_for_before_sound = frame_spout(:,1)-time_window-time_window; %1sec before moving spout
base_for_sound = frame_sound-time_window; %1sec before moving spout
%base_for_after_sound = frame_spout(:,3)-time_window; %1sec before moving spout
base_for_after_sound = frame_sound_off-time_window; %1sec before moving spout

frame_thre(:,1) = frame_spout(:,1)-time_window; %Before moving spout
frame_thre(:,2) = frame_spout(:,1); %move start of spout
frame_thre(:,3) = frame_sound; %Sound_on
frame_thre(:,4) = frame_sound_off-time_window; %Sound_end
frame_thre(:,5) = frame_sound_off; %Sound_end
frame_thre(:,6) = frame_spout(:,3); %Spout move
frame_thre(:,7) = frame_choice_select-time_window; %Before choice
frame_thre(:,8) = frame_choice_select; %After choice 0 sec
frame_thre(:,9) = frame_choice_select + 1000; %After choice 1 sec
frame_thre(:,10) = frame_choice_select + 2000; %After choice 2 sec
frame_thre(:,11) = frame_sound - time_window; %Before sound

%Use only during use_trial
base_for_before_sound = base_for_before_sound(use_trial,:);
base_for_sound = base_for_sound(use_trial,:);
base_for_after_sound = base_for_after_sound(use_trial,:);
frame_thre = frame_thre(use_trial,:);

[length_frame_thre, test_frame_thre] = size(frame_thre)

% for i = 1:7
%     temp = frame_thre(:,i+1) - frame_thre(:,i);
%     temp = find(temp < 0);
%     [i,length(temp)]
%     pause(0.5)
% %     if length(temp) ~= 0
% %     end
% end
p_thre = 3;

cd(spike_dir);

tif_name = dir('task_spike_stripe*.mat'); %get all the tif files
length_tif = length(tif_name);
max_tif = length_tif;

spike_multi = nan(max_tif,1);
p_task = nan(max_tif,test_frame_thre);
p_task_binary = zeros(max_tif,test_frame_thre);

for file_count = 1:max_tif
    %temp_file = tif_name(file_count).name;
    %temp_file
    temp_file = sprintf('task_spike_stripe20210520_%d',file_count);
    temp_file
    
    clear data
    data = load(temp_file); %spike_mark
    
    spike_mark = data.spike_mark;
    spike_multi(file_count) = data.single_multi;
    
    spike_base1 = nan(length(frame_thre),1);
    spike_base2 = nan(length(frame_thre),1);
    spike_base3 = nan(length(frame_thre),1);
    spike_frame = nan(length(frame_thre),test_frame_thre);
    for i = 1:length_frame_thre %trial
        base_frame1 = [base_for_before_sound(i) : base_for_before_sound(i)+time_window-1];
        base_frame2 = [base_for_sound(i) : base_for_sound(i)+time_window-1];
        base_frame3 = [base_for_after_sound(i) : base_for_after_sound(i)+time_window-1];
        spike_base1(i) = mean(spike_mark(base_frame1));
        spike_base2(i) = mean(spike_mark(base_frame2));
        spike_base3(i) = mean(spike_mark(base_frame3));
        for j = 1:test_frame_thre
            if isnan(frame_thre(i,j)) == 0
                temp_frame = [frame_thre(i,j) : frame_thre(i,j)+time_window-1];
                spike_frame(i,j) = mean(spike_mark(temp_frame));
            end
        end
    end
    %Kentei
    p_task(file_count,1) = signrank(spike_frame(:,1),spike_base1,'tail','right'); %Before sound
    p_task(file_count,2) = signrank(spike_frame(:,2),spike_base1,'tail','right'); %Before sound
    p_task(file_count,3) = signrank(spike_frame(:,3),spike_base2,'tail','right'); %Sound
    p_task(file_count,4) = signrank(spike_frame(:,4),spike_base2,'tail','right'); %Sound
    p_task(file_count,5) = signrank(spike_frame(:,5),spike_base2,'tail','right'); %Sound
    for j = 6:test_frame_thre-1
        temp_p1 = signrank(spike_frame(:,j),spike_base2);
        temp_p2 = signrank(spike_frame(:,j),spike_base3);
        %p_task(file_count,j) = min(temp_p1, temp_p2);
        p_task(file_count,j) = temp_p2;
    end
    p_task(file_count,11) = signrank(spike_frame(:,11),spike_base1,'tail','right'); %Before sound
    
    neuron_base(file_count).matrix = [spike_base1,spike_base2,spike_base3];
    neuron_frame(file_count).matrix = spike_frame;
end
log_p = -log10(p_task);
log_p = find(log_p > p_thre);
p_task_binary(log_p) = 1;

[length(spike_multi), sum(abs(spike_multi-2))]
sum(p_task_binary)

cd(pathname)

save sig_task_neurons_20220201 p_thre p_task_binary p_task spike_multi neuron_base neuron_frame

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
