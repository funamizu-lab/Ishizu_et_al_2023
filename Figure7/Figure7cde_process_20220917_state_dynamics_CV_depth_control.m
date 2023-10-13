
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
function Figure7cde_process_20220917_state_dynamics_CV_depth_control(folders)

close all
analysis_dir = eval(folders);
for i = 1:4
    prior_mode_all(i).matrix = [];
    sound_mode_all(i).matrix = [];
end

count = 0;
for i = 1:length(analysis_dir)
    [i,length(analysis_dir)]
    
    [mean_prior_mode,mean_sound_mode,temp_cosine,temp_prior_EV,temp_stim_EV,temp_base_EV] = ...
        State_dynamics_20220917_QR_plot_depth_control(analysis_dir{i});
    
    if length(temp_cosine) ~= 0
        count = count + 1;
        
        prior_EV(count,1) = temp_prior_EV;
        stim_EV(count,1) = temp_stim_EV;
        base_EV(count,1) = temp_base_EV;
        for j = 1:4
            prior_mode_all(j).matrix(count,:) = mean_prior_mode(j,:);
            sound_mode_all(j).matrix(count,:) = mean_sound_mode(j,:);
        end
    end
end
% signrank(prior_EV,base_EV)
% signrank(stim_EV,base_EV)
% signrank(prior_EV, stim_EV)

%% Fig 7d 
cd('G:\upload_code\Figure7\Fig7cde');
name = extractBefore(folders,'_');
sdata = struct();% source data 
temp_prior = (rand(length(prior_EV),1)-0.5)*0.2 + 1;
temp_stim = (rand(length(stim_EV),1)-0.5)*0.2 + 2;

figure; hold on
boxplot([prior_EV, stim_EV])
plot(temp_prior,prior_EV,'k.')
plot(temp_stim,stim_EV,'k.')
plot([0.5 2.5],[mean(base_EV), mean(base_EV)])
ylim([0 0.1])
sdata.block_mode=prior_EV; 
sdata.soundchoice_mode=stim_EV;
 
T = struct2table(sdata);
writetable(T, ['source fig 7d_',name,'.csv']);

%% fig 7e right
time_window2 = 2300:2499;

prior11 = prior_mode_all(1).matrix;
prior10 = prior_mode_all(2).matrix;
prior01 = prior_mode_all(3).matrix;
prior00 = prior_mode_all(4).matrix;

sound11 = sound_mode_all(1).matrix;
sound10 = sound_mode_all(2).matrix;
sound01 = sound_mode_all(3).matrix;
sound00 = sound_mode_all(4).matrix;

% top
get_difference_priormode(prior11, prior10, prior01, prior00, time_window2,name);

% bottom
get_difference_soundmode(sound11, sound10, sound01, sound00, time_window2,name);

%% fig 7e left
figure; hold on% top
sdata = struct();% source data 
[mean_trace,~,se_trace] = plot_mean_se_moto(prior11,[1 0 0],2);
sdata.block_Right_Right_mean=mean_trace';
sdata.block_Right_Right_se  =se_trace';

[mean_trace,~,se_trace] = plot_mean_se_moto(prior01,[0 102 0]./255,2);
sdata.block_Left_Right_mean =mean_trace'; 
sdata.block_Left_Right_se   =se_trace'; 

[mean_trace,~,se_trace] = plot_mean_se_moto(prior10,[255 122 0]./255,2);
sdata.block_Right_Left_mean =mean_trace';
sdata.block_Right_Left_se   =se_trace';

[mean_trace,~,se_trace] = plot_mean_se_moto(prior00,[0 0 1],2);
sdata.block_Left_Left_mean  =mean_trace';
sdata.block_Left_Left_se    =se_trace';
set(gca,'xlim',[0 3001])


figure;hold on% bottom
[mean_trace,~,se_trace] = plot_mean_se_moto(sound11,[1 0 0],2);
sdata.sound_Right_Right_mean=mean_trace';
sdata.sound_Right_Right_se  =se_trace';

[mean_trace,~,se_trace] = plot_mean_se_moto(sound01,[0 102 0]./255,2);
sdata.sound_Left_Right_mean =mean_trace'; 
sdata.sound_Left_Right_se   =se_trace'; 

[mean_trace,~,se_trace] = plot_mean_se_moto(sound10,[255 122 0]./255,2);
sdata.sound_Right_Left_mean =mean_trace';
sdata.sound_Right_Left_se   =se_trace';

[mean_trace,~,se_trace] = plot_mean_se_moto(sound00,[0 0 1],2);
sdata.sound_Left_Left_mean  =mean_trace';
sdata.sound_Left_Left_se    =se_trace';
set(gca,'xlim',[0 3001])

T = struct2table(sdata);
writetable(T, ['source fig 7e left_',name,'.csv']);

%% fig 7c
if(ismember({'auc'},name))
figure; hold on
sdata = struct();% source data 
plot(mean(prior11),mean(sound11),'color',[1 0 0])
sdata.Right_Right_block=mean(prior11)';
sdata.Right_Right_sound=mean(sound11)';

plot(mean(prior01),mean(sound01),'color',[0 102 0]./255);
sdata.Left_Right_block =mean(prior01)'; 
sdata.Left_Right_sound =mean(sound01)'; 

plot(mean(prior10),mean(sound10),'color',[255 122 0]./255)
sdata.Right_Left_block =mean(prior10)';
sdata.Right_Left_sound =mean(sound10)';

plot(mean(prior00),mean(sound00),'color',[0 0 1])
sdata.Left_Left_block  =mean(prior00)';
sdata.Left_Left_sound  =mean(sound00)';

% plot(mean(prior11(:,1)),mean(sound11(:,1)),'*','color',[1 0 0])
% plot(mean(prior01(:,1)),mean(sound01(:,1)),'*','color',[0 102 0]./255)
% plot(mean(prior10(:,1)),mean(sound10(:,1)),'*','color',[255 122 0]./255)
% plot(mean(prior00(:,1)),mean(sound00(:,1)),'*','color',[0 0 1])

T = struct2table(sdata);
writetable(T, ['source fig7c ',name,'.csv']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function get_difference_priormode(prior11, prior10, prior01, prior00, time_window2,name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ave_prior11 = mean(prior11(:,time_window2),2);
ave_prior10 = mean(prior10(:,time_window2),2);
ave_prior01 = mean(prior01(:,time_window2),2);
ave_prior00 = mean(prior00(:,time_window2),2);

figure; 
sdata = struct();% source data 
subplot(1,2,1); hold on
boxplot([ave_prior11,ave_prior10])
plot([ave_prior11,ave_prior10]')

subplot(1,2,2); hold on
boxplot([ave_prior01,ave_prior00])
plot([ave_prior01,ave_prior00]')

sdata.block_Right_Right=ave_prior11; 
sdata.block_Right_Left =ave_prior10;
sdata.block_Left_Right =ave_prior01; 
sdata.block_Left_Left  =ave_prior00;
 
T = struct2table(sdata);
writetable(T, ['source fig 7e right top_',name,'.csv']);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function get_difference_soundmode(prior11, prior10, prior01, prior00, time_window2,name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ave_prior11 = mean(prior11(:,time_window2),2);
ave_prior10 = mean(prior10(:,time_window2),2);
ave_prior01 = mean(prior01(:,time_window2),2);
ave_prior00 = mean(prior00(:,time_window2),2);

figure; 
sdata = struct();% source data 
subplot(1,2,1); hold on
boxplot([ave_prior11,ave_prior01])
plot([ave_prior11,ave_prior01]')

subplot(1,2,2); hold on
boxplot([ave_prior10,ave_prior00])
plot([ave_prior10,ave_prior00]')

sdata.block_Right_Right=ave_prior11; 
sdata.block_Right_Left =ave_prior10;
sdata.block_Left_Right =ave_prior01; 
sdata.block_Left_Left  =ave_prior00;
 
T = struct2table(sdata);
writetable(T, ['source fig 7e right bottom_',name,'.csv']);
return
