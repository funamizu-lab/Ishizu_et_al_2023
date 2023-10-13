
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
function Figure6e_process_20220520_0912_SLR_limit_neuron_depth_control(folders)
close all
[analysis_dir,short_length] = eval(folders);

%selected_window = [10 11 12 15 17];
SLR_sound = 'SLR_20220520_glmnet_sound_depth.mat';
SLR_choice= 'SLR_20220520_glmnet_choice_depth.mat';
SLR_prior = 'SLR_20220520_0912_glmnet_priorV_para2_depth.mat';
SLR_value = 'SLR_20220520_0912_glmnet_value_para2_depth.mat';

thre_neuron = 20;

for i = 1:length(short_length)
        use_frame(i,:) = [1,2,4];
end
% sound_on_frame = 11;
% sound_long_off_frame = 15;
% prior_frame = 10;
% use_frame = [prior_frame, sound_on_frame, sound_long_off_frame];
[~,nan_data(:,4)] = ...
    get_average_correct_rate2(analysis_dir, SLR_value, use_frame, 0, thre_neuron);

[correct_l_sound,nan_data(:,1)] = ...
    get_average_correct_rate2(analysis_dir, SLR_sound, use_frame, 0, thre_neuron);

[~,nan_data(:,2)] = ...
    get_average_correct_rate2(analysis_dir, SLR_choice, use_frame, 0, thre_neuron);

[~,nan_data(:,3)] = ...
    get_average_correct_rate2(analysis_dir, SLR_prior, use_frame, 0, thre_neuron);

%in each use frame change the value for long trials
correct_l_sound = get_use_frame_correct_rate(correct_l_sound,use_frame);

%Check nan, sound and choice decode should not have NAN!!
temp = max(nan_data);
if temp(1) == 1 || temp(2) == 1 || temp(4) == 1%Sound or choice
    hoge
end

%% fig6e left
cd('G:\upload_code\Figure6\Fig6e');
sdata = struct();% source data 
get_correct_plot_20220520(correct_l_sound);

sdata.initial_sound=correct_l_sound(:,2);
sdata.late_sound =correct_l_sound(:,3);
T = struct2table(sdata);
writetable(T, 'source fig6e left.csv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function new_correct_sound= ...
    get_use_frame_correct_rate(correct_l_sound,use_frame)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[size_mouse,~] = size(correct_l_sound);

for i = 1:size_mouse
    for j = 1:3
        new_correct_sound(i,j) = correct_l_sound(i,use_frame(i,j));
    end
end
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function get_correct_plot_20220520(correct_l_sound)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
use_frame = 1:3;

for i = 1:length(use_frame)
    l_sound(:,i) = correct_l_sound(:,use_frame(i));
end

figure; hold on
plot(l_sound(:,2),l_sound(:,3),'b.')
plot([0.3 1],[0.3 1],'k')
set(gca,'xlim',[0.3 1],'ylim',[0.3 1])
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [correct_l_all, nan_data] = ...
    get_average_correct_rate2(analysis_dir, SLR_name, use_frame, psycho_sign, thre_neuron)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

count = 0;
count_neuron = 0;
for i = 1:length(analysis_dir)
    [i,length(analysis_dir)]

    cd(analysis_dir{i});
    load(SLR_name);
    
    %Check the number of neurons
    temp = dir('sig2_task_neurons_2022*');
    if length(temp) ~= 1
        hoge
    end
    load(temp.name);

    temp = dir('depth_spike_20220517*');
    if length(temp) == 1
        load(temp.name);
        %spike_depth def_depth length_neuron
        if size(p_task_long,1) ~= length_neuron
            hoge
        end
        depth_neuron = find(spike_depth <= def_depth);
    elseif length(temp) == 0
        depth_neuron = 1:size(p_task_long,1); %Use all the neurons
    else
        hoge
    end
    data = length(depth_neuron);
    
    if data >= thre_neuron %Ready to analyze for all the definition
        count_neuron = count_neuron + 1;
    if length(lasso_sound_l) ~= 1 %not nan && neuron > thre
        count = count + 1;
    
        [~,correct_l] = ...
        Population_decoder_20220518_SLR_process(analysis_dir{i},SLR_name, use_frame(i,:), psycho_sign);

        correct_l_all(count,:) = correct_l;
        
        nan_data(count_neuron,1) = 0;
    else
        nan_data(count_neuron,1) = 1;
    end
    end
end

% %More easy way is to analyze the correct rate with last-tone neurons and 
% %first-tone neurons. and see how much the correct rate deteriorate with the
% %neuron group!!

return

