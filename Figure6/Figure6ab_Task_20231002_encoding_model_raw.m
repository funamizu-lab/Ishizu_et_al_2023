
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
function Figure6ab_Task_20231002_encoding_model_raw

close all
pathname='G:\Ishizu_data\IntermediateFiles\auditory\a04\2021-02-18_a04_1AC2PPC_Right\recording1_task';
cd(pathname)

Kernel_folder = 'KernelEncoding_20231002';
Kernel_file = 'KernelEncoding_20231002_kernel.mat';

temp = dir('task_frame*'); %task file
if length(temp) ~= 1
    hoge
end
load(temp.name);
temp = dir('Bpod*'); %Bpod file
if length(temp) ~= 1
    hoge
end
load(temp.name);
Bpod_file = temp.name;

[~,~,~,~,~,~,~,use_trial_all] = Dual_get_basic_task_structure_20210204(Bpod_file);

temp = dir('sig_task*');
if length(temp) ~= 1
    hoge
end
load(temp.name);

Correct_side = Correct_side(use_trial_all);
Chosen_side  = Chosen_side(use_trial_all);
% StimDuration = StimDuration(use_trial_all);
%Get task parameter
low  = find(Correct_side == 0);
high = find(Correct_side == 1);
left  = find(Chosen_side == 0);
right = find(Chosen_side == 1);
%Get correct and error trials
low_left = intersect(low,left);
low_right = intersect(low,right);
high_left = intersect(high,left);
high_right = intersect(high,right);

%Make Sig neuron
p_thre = 100;
use_p_sound = p_task(:,[3,4]);
use_p_sound = min(use_p_sound,[],2);
use_p_sound = -log10(use_p_sound);
use_p_sound = find(use_p_sound > p_thre);
if length(use_p_sound) > 15
    sig_neuron = use_p_sound(1:15);
else
    sig_neuron = use_p_sound;
end

cd(Kernel_folder);
load(Kernel_file); %kernel_all, kernel_size_y, moto_frame5

filter_sound_on = moto_frame_sound; %6 tone evidence with long and short
filter_sound_on_frame = find(filter_sound_on == 1);
% if length(filter_sound_on_frame) ~= length(use_trial_all)
%     [length(filter_sound_on_frame), length(use_trial_all)]
% end

%Start analyze the regression parameters
pre_sound = 4; %0.5sec
post_sound = 20; %2sec

trial_start = 18;
trial_end = 21;

plot_time = filter_sound_on_frame(trial_start)-pre_sound : filter_sound_on_frame(trial_end) + post_sound;
plot_sound_time = filter_sound_on_frame(trial_start:trial_end) - filter_sound_on_frame(trial_start)+pre_sound + 1;

cd(pathname) %back to the task folder
%% Fig 6ab
file_count = 8;
% for file_count = 1:length(sig_neuron)
cd(Kernel_folder)
temp_file = sprintf('KernelEncoding_20231002_glmnet_%d.mat',sig_neuron(file_count));
data = load(temp_file); %spike_mark

model_spike = data.Ridge_neuron.y_all;
filter_spike = data.spike_filter_new; %normalized spikes

%Need to change to moto spike numbers
min_filter_spike = min(filter_spike);
model_spike = exp(model_spike);
model_spike = model_spike - 1 + min_filter_spike;

%Filtered raw spikes
model_trace_sound = nan(length(filter_sound_on_frame), pre_sound+post_sound+1);
filter_trace_sound = nan(length(filter_sound_on_frame), pre_sound+post_sound+1);
for i = 1:length(filter_sound_on_frame)
    temp_pre  = filter_sound_on_frame(i) - pre_sound;
    temp_post = filter_sound_on_frame(i) + post_sound;
    
    model_trace_sound(i,:) = model_spike(temp_pre:temp_post);
    filter_trace_sound(i,:) = filter_spike(temp_pre:temp_post);
end

% Fig 6a right
%Plot the Result of kernel modeling
cd('G:\upload_code\Figure6\Fig6a');
figure; hold on
sdata = struct();% source data 
plot(filter_spike(plot_time),'k') %Test and filter should be tested!!
plot(model_spike(plot_time),'b') %Test and filter should be tested!!
for i = 1:length(plot_sound_time)
    plot([plot_sound_time(i),plot_sound_time(i)],[0 1],'r')
end
sdata.activity = filter_spike(plot_time);
sdata.model = model_spike(plot_time);
T = struct2table(sdata);
writetable(T, 'source fig6a right.csv');

% Fig 5b
cd('G:\upload_code\Figure6\Fig6b');
figure
sdata = struct();% source data 
subplot(1,2,1); hold on
plot(mean(model_trace_sound(low_left,:)),'color',[0 0 1])
plot(mean(model_trace_sound(low_right,:)),':','color',[0 0 1])
plot(mean(model_trace_sound(high_right,:)),'color',[1 0 0])
plot(mean(model_trace_sound(high_left,:)),':','color',[1 0 0])
sdata.model_toneLeft_correct  = mean(model_trace_sound(low_left,:))';
sdata.model_toneLeft_error= mean(model_trace_sound(low_right,:))';
sdata.model_toneRight_correct  = mean(model_trace_sound(high_right,:))';
sdata.model_toneRight_error= mean(model_trace_sound(high_left,:))';

subplot(1,2,2)
hold on
plot(mean(filter_trace_sound(low_left,:)),'color',[0 0 1])
plot(mean(filter_trace_sound(low_right,:)),':','color',[0 0 1])
plot(mean(filter_trace_sound(high_right,:)),'color',[1 0 0])
plot(mean(filter_trace_sound(high_left,:)),':','color',[1 0 0])
sdata.smooth_toneLeft_correct  = mean(filter_trace_sound(low_left,:))';
sdata.smooth_toneLeft_error= mean(filter_trace_sound(low_right,:))';
sdata.smooth_toneRight_correct  = mean(filter_trace_sound(high_right,:))';
sdata.smooth_toneRight_error= mean(filter_trace_sound(high_left,:))';

T = struct2table(sdata);
writetable(T, 'source fig6b left.csv');