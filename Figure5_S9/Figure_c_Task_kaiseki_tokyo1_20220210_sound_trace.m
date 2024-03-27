
%{
----------------------------------------------------------------------------
First_take number of frames in each tif files
Analyzing imaging data simply
At least for the correct rate
----------------------------------------------------------------------------
%}
function Figure_c_Task_kaiseki_tokyo1_20220210_sound_trace

close all
flag=2;% auc:1(Fig5) / fof:2(FigS9)

if(flag==1)
    Fig ='5';
    id =338;% spike id
    directry = 'G:\Ishizu_data\Tokyo_ephys_ishizu\auditory\i24\2021-10-06_i24_AC_left_OK\recording1_task';
elseif(flag==2)
    Fig ='S9';% auc:5 / fof:S9
    id =526;% spike id
    directry = 'G:\Ishizu_data\Tokyo_ephys_ishizu\fof\i24\2021-10-01_i24_FOF_left_OK\recording1_task';
end

cd(directry)

temp = dir('task_frame*');
if length(temp) ~= 1
    hoge
end
load(temp.name);
%frame_start / frame_sound / frame_end

temp = dir('Bpod*');
if length(temp) ~= 1
    hoge
end
load(temp.name);

spike_dir = dir('spike_ch*');
if length(spike_dir) ~= 1
    hoge
end
spike_dir = spike_dir.name;

%Get task parameter
Choice_trial = find(Outcome == 1 | Outcome == 2);
low  = find(Correct_side == 0);
high = find(Correct_side == 1);
stim_length = unique(StimDuration);
Long  = find(StimDuration == stim_length(2));
left = find(Chosen_side == 0);
right = find(Chosen_side == 1);

trial_LL = intersect(low,Long);
trial_HL = intersect(high,Long);

trial_LL = intersect(trial_LL,Choice_trial);
trial_HL = intersect(trial_HL,Choice_trial);

trial_LL_correct = intersect(trial_LL,left);
trial_LL_error = intersect(trial_LL,right);
trial_HL_correct = intersect(trial_HL,right);
trial_HL_error = intersect(trial_HL,left);

%Get the sound frame
pre_frame = 200; %1sec
post_frame = 200; %1sec
pre_frame2 = 500; %1sec
post_frame2 = 2000; %1sec

cd(spike_dir);
temp_file = sprintf('task_spike_stripe20210520_%d',id);
data = load(temp_file); %spike_mark
[~,~,spike_trace] = get_sound_response3(data.spike_mark, data.spike_filter, frame_sound, pre_frame, post_frame, pre_frame2, post_frame2);

%% for Fig d
cd(directry);
% significance 
temp = dir('sig2_task_neurons_2022*');
if length(temp) ~= 1
    hoge
end
load(temp.name,'p_task_long');
p_neuron=p_task_long(id,:);
 
% Activity
temp = dir('Tokyo2_20220912_sound_choice_separate*');
if length(temp) ~= 1
    hoge
end
load(temp.name,'long1_correct','long1_error');
correct_neuron =long1_correct(id,:);
error_neuron =long1_error(id,:);

%ROC
temp = dir('ROC_20230701_long_short*');
if length(temp) ~= 1
    hoge
end
load(temp.name,'long_sound_ROC1','long_choice_ROC1');
sROC_neuron =long_sound_ROC1(id,:);
cROC_neuron=long_choice_ROC1(id,:);

%% Fig c
cd('G:\upload_code\Figure5_S9\c');
sdata = struct();% source data 
figure; hold on
plot(mean(spike_trace(trial_LL_correct,:)),'color',[0 0 1]) %Blue
plot(mean(spike_trace(trial_LL_error,:)),  'color',[34 139 34]./255) %Green
plot(mean(spike_trace(trial_HL_correct,:)),'color',[255 0 0]./255) %Red
plot(mean(spike_trace(trial_HL_error,:)),  'color',[255 143 34]./255) %Orange

sdata.RightRight_mean= mean(spike_trace(trial_HL_correct,:))'*1e3;
sdata.RightLeft_mean = mean(spike_trace(trial_HL_error,:))'*1e3;
sdata.LeftLeft_mean  = mean(spike_trace(trial_LL_correct,:))'*1e3;
sdata.LeftRight_mean = mean(spike_trace(trial_LL_error,:))'*1e3;
T = struct2table(sdata);
writetable(T, ['source fig',Fig,'c.csv']);

%%% save data
save(['fig',Fig,'cNeuron.mat'],'p_neuron','correct_neuron','error_neuron',...
    'sROC_neuron','cROC_neuron');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spike_count,p,spike_trace] = get_sound_response3(spike_mark, spike_filter, frame_sound, pre_frame, post_frame, pre_frame2, post_frame2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spike_count = nan(length(frame_sound),2);
spike_trace = nan(length(frame_sound),pre_frame2+post_frame2);

for i = 1:length(frame_sound)
    temp_pre = frame_sound(i)-pre_frame : frame_sound(i)-1;
    temp_post= frame_sound(i) : frame_sound(i)+post_frame-1;
    temp_all = frame_sound(i)-pre_frame2 : frame_sound(i)+post_frame2-1;
    
    temp_pre = spike_mark(temp_pre);
    temp_post = spike_mark(temp_post);
    spike_count(i,:) = [sum(temp_pre), sum(temp_post)];
    
    spike_trace(i,:) = spike_filter(temp_all);
end

p = signrank(spike_count(:,1),spike_count(:,2));
return
