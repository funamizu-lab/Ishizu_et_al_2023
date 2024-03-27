
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
function FigureS13ac_process_SLR_limit_neuron_depth_control_dprime(folders)

close all
[analysis_dir,short_length] = eval(folders);

%selected_window = [10 11 12 15 17];
SLR_name =  'SLR_20220520_glmnet_sound_depth.mat';
SLR_choice= 'SLR_20220520_glmnet_choice_depth.mat';
SLR_prior = 'SLR_20220520_0912_glmnet_priorV_para2_depth.mat';
SLR_value = 'SLR_20220520_0912_glmnet_value_para2_depth.mat';

thre_neuron = 20;

for i = 1:length(short_length)
    use_frame(i,:) = [1,2,4];
end

[nan_data(:,4),dp_l_value, dp_ls_value] = ...
    get_average_correct_rate3(analysis_dir, SLR_value, use_frame, 0, thre_neuron, 3); %Value

[nan_data(:,1),dp_l_sound] = ...
    get_average_correct_rate3(analysis_dir, SLR_name, use_frame, 0, thre_neuron, 0); %Sound

[nan_data(:,2),dp_l_choice] = ...
    get_average_correct_rate3(analysis_dir, SLR_choice, use_frame, 0, thre_neuron, 1); %Choice

[nan_data(:,3),dp_l_prior] = ...
    get_average_correct_rate3(analysis_dir, SLR_prior, use_frame, 0, thre_neuron, 2); %Prior

%in each use frame change the value for long trials
[dp_l_value,dp_l_sound,dp_l_choice,dp_l_prior] = ...
    get_use_frame_correct_rate(dp_l_prior,dp_ls_value,dp_l_value,dp_l_sound,dp_l_choice,use_frame);


%Check nan, sound and choice decode should not have NAN!!
temp = max(nan_data);
if temp(1) == 1 || temp(2) == 1 || temp(4) == 1%Sound or choice
    hoge
end
max_nan_data = max(nan_data,[],2);
max_nan_data = find(max_nan_data == 1); %Only from prior

%%% FigS12a %%%
cd('G:\upload_code\FigureS13\FigS13a');
name = extractBefore(folders,'_');
get_correct_prior_include3(dp_l_sound,dp_l_choice,dp_l_prior,dp_l_value,...
    use_frame, max_nan_data, [-0.5 3],name);

%%% FigS12c %%% 
cd('G:\upload_code\FigureS13\FigS13c');
get_correct_plot_20221206(dp_l_sound,[-1 4],name);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [new_correct_value,new_correct_sound,new_correct_choice,new_correct_prior] = ...
    get_use_frame_correct_rate(correct_l_prior,correct_ls_value,correct_l_value,correct_l_sound,correct_l_choice,use_frame)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[size_mouse,~] = size(correct_l_prior);

for i = 1:size_mouse
    new_correct_value(i,1) = correct_l_prior(i,use_frame(i,1));
    new_correct_value(i,2) = correct_ls_value(i,use_frame(i,2));
    new_correct_value(i,3) = correct_l_value(i,use_frame(i,3));
    for j = 1:3
        new_correct_sound(i,j) = correct_l_sound(i,use_frame(i,j));
        new_correct_choice(i,j) = correct_l_choice(i,use_frame(i,j));
        new_correct_prior(i,j) = correct_l_prior(i,use_frame(i,j));
    end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function get_correct_prior_include3(l_sound,l_choice,l_prior,l_value,...
    use_frame, max_nan_data, plot_scale, name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,size_frame] = size(use_frame);

%Assuming that the nan happens only at the prior decoding!!
l_sound(max_nan_data,:) = [];
l_choice(max_nan_data,:) = [];

[size_y,~] = size(l_prior);
temp_box = [zeros(size_y,1); ones(size_y,1); ones(size_y,1) * 2];

%%% fig S13a %%%
figure
tag={'before sound','initial sound','late sound'};
for i = 1:size_frame
    sdata = struct();% source data
    temp = [l_sound(:,i); l_choice(:,i); l_prior(:,i)];
    subplot(1,size_frame,i); hold on
    boxplot(temp, temp_box)
    plot([l_sound(:,i), l_choice(:,i), l_prior(:,i)]')
    set(gca,'ylim',plot_scale)
%     p_SCP_compare(i,1) = signrank(l_sound(:,i), l_choice(:,i));
%     p_SCP_compare(i,2) = signrank(l_sound(:,i), l_prior(:,i));
%     p_SCP_compare(i,3) = signrank(l_prior(:,i), l_choice(:,i));
%     
%     p_SCP_compare(i,4) = signrank(l_sound(:,i), l_value(:,i));
%     p_SCP_compare(i,5) = signrank(l_choice(:,i), l_value(:,i));
%     p_SCP_compare(i,6) = signrank(l_prior(:,i), l_value(:,i));

    %%% source data %%%
    sdata.sound = l_sound(:,i);
    sdata.choice= l_choice(:,i);
    sdata.prior = l_prior(:,i);
    T = struct2table(sdata);
    writetable(T, ['source figS13a ',name,' ',tag{i},'.csv']);
 
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function get_correct_plot_20221206(correct_l_sound, plot_scale, name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
use_frame = 1:3;

for i = 1:length(use_frame)
    l_sound(:,i) = correct_l_sound(:,use_frame(i));
end

figure; hold on
plot(l_sound(:,2),l_sound(:,3),'b.')
plot(plot_scale,plot_scale,'k')
set(gca,'xlim',plot_scale,'ylim',plot_scale)

%%% source data %%%
sdata.x = l_sound(:,2);
sdata.y = l_sound(:,3);
T = struct2table(sdata);
writetable(T, ['source figS13c ',name,'.csv']);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nan_data, dp_l_all, dp_ls_all] = ...
    get_average_correct_rate3(analysis_dir, SLR_name, use_frame, psycho_sign, thre_neuron, sign_category)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sign_category
%0: sound
%1: choice
%2: prior
%3: value

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

        [~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~, ...
        dp_l,~,dp_ls] = ...
        Population_decoder_20221206_SLR_process_dprime(analysis_dir{i},SLR_name, use_frame(i,:), psycho_sign, sign_category);

        dp_l_all(count,:) = dp_l;
        dp_ls_all(count,:) = dp_ls;
        
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