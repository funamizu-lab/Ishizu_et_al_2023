
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
function Figure_b_process_20220516_norm_ave_trace_depth_control(folders)

close all
analysis_dir = eval(folders);

norm_long_all = [];
norm_short_all = [];
long_01all = [];
short_01all = [];
long_time_all = [];
short_time_all = [];
for i = 1:length(analysis_dir)
    [i,length(analysis_dir)]
    
    [norm_spike_long,norm_spike_short,spike_long01,spike_short01,...
        long_max_time,short_max_time,neuron_number(i,:)] = ...
        Task_kaiseki_tokyo1_20220516_make_ave_depth_control(analysis_dir{i});

    norm_long_all = [norm_long_all; norm_spike_long];
    norm_short_all = [norm_short_all; norm_spike_short];
    long_01all = [long_01all; spike_long01];
    short_01all = [short_01all; spike_short01];
    long_time_all = [long_time_all; long_max_time];
    short_time_all = [short_time_all; short_max_time];
end
delete(gcp('nocreate'))

[time_long,sort_long] = sort(long_time_all);

%% fig b heatmap
cd('G:\upload_code\Figure5_S9\b');
figure
sdata = struct();% source data 
imagesc(long_01all(sort_long,:))
sdata.time = long_01all(sort_long,:);
T = struct2table(sdata);
% writetable(T, 'source fig5b heatmap.csv');
writetable(T, 'source figS8b heatmap.csv');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_basic_tone_choice_index(Sound_prefer_index, use_neuron)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
subplot(1,3,1) %Short sound
hold on
plot(Sound_prefer_index(use_neuron,4), Sound_prefer_index(use_neuron,1),'b.')
plot([0 0],[-1 1],'k--')
plot([-1 1],[0 0],'k--')
set(gca,'xlim',[-1 1],'ylim',[-1 1])

subplot(1,3,2) %Long sound
hold on
plot(Sound_prefer_index(use_neuron,5), Sound_prefer_index(use_neuron,2),'b.')
plot([0 0],[-1 1],'k--')
plot([-1 1],[0 0],'k--')
set(gca,'xlim',[-1 1],'ylim',[-1 1])
subplot(1,3,3) %Long sound end
hold on
plot(Sound_prefer_index(use_neuron,6), Sound_prefer_index(use_neuron,3),'b.')
plot([0 0],[-1 1],'k--')
plot([-1 1],[0 0],'k--')
set(gca,'xlim',[-1 1],'ylim',[-1 1])

signrank(Sound_prefer_index(:,4))
signrank(Sound_prefer_index(:,5))
signrank(Sound_prefer_index(:,6))

figure %Onset and Offset of sound @ correct trials
hold on
plot(Sound_prefer_index(use_neuron,3),Sound_prefer_index(use_neuron,2),'b.')
plot([-1 1],[-1 1],'k')
plot([0 0],[-1 1],'k--')
plot([-1 1],[0 0],'k--')
set(gca,'xlim',[-1 1],'ylim',[-1 1])
signrank(Sound_prefer_index(use_neuron,2), Sound_prefer_index(use_neuron,3))
signrank(abs(Sound_prefer_index(use_neuron,2)), abs(Sound_prefer_index(use_neuron,3)))
[median(Sound_prefer_index(use_neuron,2)), median(Sound_prefer_index(use_neuron,3))]

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_sigR_sig_L(sig_all,sig_R,sig_L,tuning0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sig_all1=sig_all(1).matrix(tuning0,:);
sig_all2=sig_all(2).matrix(tuning0,:);
sig_all3=sig_all(3).matrix(tuning0,:);
sig_R1=sig_R(1).matrix(tuning0,:);
sig_R2=sig_R(2).matrix(tuning0,:);
sig_R3=sig_R(3).matrix(tuning0,:);
sig_L1=sig_L(1).matrix(tuning0,:);
sig_L2=sig_L(2).matrix(tuning0,:);
sig_L3=sig_L(3).matrix(tuning0,:);

figure
subplot(1,3,1)
hold on
plot(nanmean(sig_all1),'k')
plot(nanmean(sig_R1),'r')
plot(nanmean(sig_L1),'b')

subplot(1,3,2)
hold on
plot(nanmean(sig_all2),'k')
plot(nanmean(sig_R2),'r')
plot(nanmean(sig_L2),'b')

subplot(1,3,3)
hold on
plot(nanmean(sig_all3),'k')
plot(nanmean(sig_R3),'r')
plot(nanmean(sig_L3),'b')

clear p
for i = 1:6
    p(1,i) = signrank(sig_R1(:,i),sig_L1(:,i));
    p(2,i) = signrank(sig_R2(:,i),sig_L2(:,i));
    p(3,i) = signrank(sig_R3(:,i),sig_L3(:,i));
end
p