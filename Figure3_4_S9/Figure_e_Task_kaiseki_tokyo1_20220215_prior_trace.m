
%{
----------------------------------------------------------------------------
First_take number of frames in each tif files
Analyzing imaging data simply
At least for the correct rate
----------------------------------------------------------------------------
%}
function Figure_e_Task_kaiseki_tokyo1_20220215_prior_trace

close all
id =260;% spike id
Fig ='S8';% mpfc:3 / auc:4 / fof:S8

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
Bpod_file = temp.name;

[~,~,~,~,~,block2,block3,use_trial]= Dual_get_basic_task_structure_20210204(Bpod_file);

spike_dir = dir('spike_ch*');
if length(spike_dir) ~= 1
    hoge
end
spike_dir = spike_dir.name;

%Check block change
TrialBlock_use = TrialBlock(use_trial);
temp = TrialBlock_use(2:length(TrialBlock_use)) - TrialBlock_use(1:length(TrialBlock_use)-1);
Blockchange_use = find(temp ~= 0) + 0.5; %start of new block

%Get task parameter
stim_length = unique(StimDuration);
Long  = find(StimDuration == stim_length(2));
left = find(Chosen_side == 0);
right= find(Chosen_side == 1);

%Adjust the trials to use
%Left right for each block
Long = intersect(Long, use_trial);
left = intersect(left, use_trial);
right = intersect(right, use_trial);
if BlockReward(2,1) > BlockReward(2,2) %Left -> Right
    left_R = intersect(left, block3);
    right_R = intersect(right, block3);
    left_L = intersect(left, block2);
    right_L = intersect(right, block2);
else %Right -> Left
    left_R = intersect(left, block2);
    right_R = intersect(right, block2);
    left_L = intersect(left, block3);
    right_L = intersect(right, block3);
end
left_R_long = intersect(left_R, Long);
left_L_long = intersect(left_L, Long);
right_R_long = intersect(right_R, Long);
right_L_long = intersect(right_L, Long);

%Get the sound frame
pre_frame = 200; %1sec
post_frame = 200; %1sec
pre_frame2 = 2000; %1sec
post_frame2 = 1000; %1sec

cd(spike_dir);
temp_file = sprintf('task_spike_stripe20210520_%d',id);
data = load(temp_file); %spike_mark
[spike_count,~,spike_trace] = get_sound_response3(data.spike_mark, data.spike_filter, frame_sound, pre_frame, post_frame, pre_frame2, post_frame2);

%Use spike pre to make the spike_count per trial
spike_count = spike_count(use_trial,1);

%% Fig e top
cd('G:\upload_code\Figure3_4_S8\e');
sdata = struct();% source data 
figure; hold on
plot(mean(spike_trace(left_R_long,:)),'color',[0 0 1]) %Blue
plot(mean(spike_trace(left_L_long,:)),  'color',[34 139 34]./255) %Green
plot(mean(spike_trace(right_R_long,:)),'color',[255 0 0]./255) %Red
plot(mean(spike_trace(right_L_long,:)),  'color',[255 143 34]./255) %Orange

% Block-Chosen
sdata.RightRight_mean= mean(spike_trace(right_R_long,:))'*1e3;
sdata.RightLeft_mean = mean(spike_trace(left_R_long,:))'*1e3;
sdata.LeftLeft_mean  = mean(spike_trace(right_L_long,:))'*1e3;
sdata.LeftRight_mean = mean(spike_trace(left_L_long,:))'*1e3;
T = struct2table(sdata);
writetable(T, ['source fig',Fig,'e top.csv']);

%% Fig e bottom
sdata = struct();% source data 
figure; hold on
for j = 1:length(Blockchange_use)
    line([Blockchange_use(j),Blockchange_use(j)], [0 1],'color',[0 0 0],'LineWidth',0.5);
end
plot(spike_count,'k.')
set(gca,'xlim',[0 length(spike_count)+1])
sdata.SpikeRate = spike_count*5; 
T = struct2table(sdata);
writetable(T, ['source fig',Fig,'e bottom.csv']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spike_count,p,spike_trace] = get_sound_response3(spike_mark, spike_filter, frame_sound, pre_frame, post_frame, pre_frame2, post_frame2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spike_count = nan(length(frame_sound),2);
spike_trace = nan(length(frame_sound),pre_frame2+post_frame2);

for i = 1:length(frame_sound)
    temp_pre  = frame_sound(i)-pre_frame : frame_sound(i)-1;
    temp_post = frame_sound(i) : frame_sound(i)+post_frame-1;
    temp_all  = frame_sound(i)-pre_frame2 : frame_sound(i)+post_frame2-1;
    
    temp_pre = spike_mark(temp_pre);
    temp_post = spike_mark(temp_post);
    spike_count(i,:) = [sum(temp_pre), sum(temp_post)];
    
    spike_trace(i,:) = spike_filter(temp_all);
end

p = signrank(spike_count(:,1),spike_count(:,2));

return