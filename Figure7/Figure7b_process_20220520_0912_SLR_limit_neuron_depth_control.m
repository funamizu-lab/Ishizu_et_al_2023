
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
function Figure7b_process_20220520_0912_SLR_limit_neuron_depth_control(folders)

close all
[analysis_dir,short_length] = eval(folders);

%selected_window = [10 11 12 15 17];
SLR_name = 'SLR_20220520_glmnet_sound_depth.mat';
SLR_choice = 'SLR_20220520_glmnet_choice_depth.mat';
SLR_prior = 'SLR_20220520_0912_glmnet_priorV_para2_depth.mat';
SLR_value = 'SLR_20220520_0912_glmnet_value_para2_depth.mat';

thre_neuron = 20;

for i = 1:length(short_length)
%     if short_length(i) == 0
        use_frame(i,:) = [1,2,4];
%     else
%         use_frame(i,:) = [1,2,4];
%     end
end
% sound_on_frame = 11;
% sound_long_off_frame = 15;
% prior_frame = 10;
% use_frame = [prior_frame, sound_on_frame, sound_long_off_frame];
[~,correct_l_value,~,~,~,~,~,~,~,~,~,correct_ls_value,~,~,nan_data(:,4)] = ...
    get_average_correct_rate2(analysis_dir, SLR_value, use_frame, 0, thre_neuron);

[~,correct_l_sound,~,~,~,~,~,~,~,~,~,~,~,~,nan_data(:,1)] = ...
    get_average_correct_rate2(analysis_dir, SLR_name, use_frame, 0, thre_neuron);

[~,correct_l_choice,~,~,~,~,~,~,~,~,~,~,~,~,nan_data(:,2)] = ...
    get_average_correct_rate2(analysis_dir, SLR_choice, use_frame, 0, thre_neuron);

[~,correct_l_prior,~,~,~,~,~,~,~,~,~,~,~,~,nan_data(:,3)] = ...
    get_average_correct_rate2(analysis_dir, SLR_prior, use_frame, 0, thre_neuron);

%in each use frame change the value for long trials
[correct_l_value,correct_l_sound,correct_l_choice,correct_l_prior] = ...
    get_use_frame_correct_rate(correct_l_prior,correct_ls_value,correct_l_value,correct_l_sound,correct_l_choice,use_frame);


%Check nan, sound and choice decode should not have NAN!!
temp = max(nan_data);
if temp(1) == 1 || temp(2) == 1 || temp(4) == 1%Sound or choice
    hoge
end
max_nan_data = max(nan_data,[],2);
max_nan_data = find(max_nan_data == 1); %Only from prior

%% Fig6b
cd('G:\upload_code\Figure7\Fig7b');
sdata = struct();% source data 
name = extractBefore(folders,'_');

get_correct_prior_include2(correct_l_sound,correct_l_choice,correct_l_prior,correct_l_value,...
    use_frame, max_nan_data)
sdata.leftpanel_sound =correct_l_sound(:,1);
sdata.leftpanel_choice=correct_l_choice(:,1);
sdata.leftpanel_prior =correct_l_prior(:,1);
sdata.middlepanel_sound =correct_l_sound(:,2);
sdata.middlepanel_choice=correct_l_choice(:,2);
sdata.middlepanel_prior =correct_l_prior(:,2);
sdata.rightpanel_sound =correct_l_sound(:,3);
sdata.rightpanel_choice=correct_l_choice(:,3);
sdata.rightpanel_prior =correct_l_prior(:,3);
T = struct2table(sdata);
writetable(T, ['source fig7b ',name,'.csv']);

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
function get_correct_prior_include2(l_sound,l_choice,l_prior,l_value,...
    use_frame, max_nan_data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,size_frame] = size(use_frame);

%Assuming that the nan happens only at the prior decoding!!
l_sound(max_nan_data,:) = [];
l_choice(max_nan_data,:) = [];
% l_value(max_nan_data,:) = [];

[size_y,~] = size(l_prior);
% temp_box = [zeros(size_y,1); ones(size_y,1); ones(size_y,1) * 2; ones(size_y,1) * 3];
temp_box = [zeros(size_y,1); ones(size_y,1); ones(size_y,1) * 2]; 
figure; 
for i = 1:size_frame
%     temp = [l_sound(:,i); l_choice(:,i); l_prior(:,i); l_value(:,i)];
    temp = [l_sound(:,i); l_choice(:,i); l_prior(:,i)];
    subplot(1,size_frame,i);hold on
    boxplot(temp, temp_box)
    plot([l_sound(:,i), l_choice(:,i), l_prior(:,i)]')
    set(gca,'ylim',[0.45 0.9])    
    
%     p_SCP_compare(i,1) = signrank(l_sound(:,i), l_choice(:,i));
%     p_SCP_compare(i,2) = signrank(l_sound(:,i), l_prior(:,i));
%     p_SCP_compare(i,3) = signrank(l_prior(:,i), l_choice(:,i));
%     
%     p_SCP_compare(i,4) = signrank(l_sound(:,i), l_value(:,i));
%     p_SCP_compare(i,5) = signrank(l_choice(:,i), l_value(:,i));
%     p_SCP_compare(i,6) = signrank(l_prior(:,i), l_value(:,i));
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [correct_s_all, correct_l_all, likeli_s_all, likeli_l_all, opt_s_all, opt_l_all, count, ...
    correct_each_s_all, likeli_each_s_all, correct_each_l_all, likeli_each_l_all, ...
    correct_ls_all, likeli_ls_all, opt_ls_all, nan_data] = ...
    get_average_correct_rate2(analysis_dir, SLR_name, use_frame, psycho_sign, thre_neuron)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,size_frame] = size(use_frame);

%spike_ch1 = 'spike_ch1';
%parpool
count = 0;
count_neuron = 0;
for i = 1:length(analysis_dir)
    [i,length(analysis_dir)]

    cd(analysis_dir{i})
    load(SLR_name)
    
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
    
        [correct_s, correct_l, opt_s, low_count_s, high_count_s, ...
        opt_l, low_count_l, high_count_l, x_edge_plot, evi_x, ...
        low_easy_s, low_mid_s, low_dif_s, high_easy_s, high_mid_s, high_dif_s, ...
        low_easy_l, low_mid_l, low_dif_l, high_easy_l, high_mid_l, high_dif_l, ...
        likelihood_s, likelihood_l, window_length, ...
        correct_each_s, likeli_each_s, correct_each_l, likeli_each_l, ...
        correct_ls, opt_ls, low_count_ls, high_count_ls, ...
        low_easy_ls, low_mid_ls, low_dif_ls, high_easy_ls, high_mid_ls, high_dif_ls, ...
        likelihood_ls] = ...
        Population_decoder_20220518_SLR_process(analysis_dir{i},SLR_name, use_frame(i,:), psycho_sign);

        correct_s_all(count,:) = correct_s;
        correct_l_all(count,:) = correct_l;
        likeli_s_all(count,:) = likelihood_s;
        likeli_l_all(count,:) = likelihood_l;
        
        correct_ls_all(count,:) = correct_ls;
        likeli_ls_all(count,:) = likelihood_ls;
        
        %for j = 1:window_length
        for j = 1:5
            low_s_all(j).matrix(count,:) = low_count_s(j,:);
            high_s_all(j).matrix(count,:) = high_count_s(j,:);
        
            low_easy_s_all(j).matrix(count,:) = low_easy_s(j,:);
            low_mid_s_all(j).matrix(count,:) = low_mid_s(j,:);
            low_dif_s_all(j).matrix(count,:) = low_dif_s(j,:);
            high_easy_s_all(j).matrix(count,:) = high_easy_s(j,:);
            high_mid_s_all(j).matrix(count,:) = high_mid_s(j,:);
            high_dif_s_all(j).matrix(count,:) = high_dif_s(j,:);

            low_l_all(j).matrix(count,:) = low_count_l(j,:);
            high_l_all(j).matrix(count,:) = high_count_l(j,:);
        
            low_easy_l_all(j).matrix(count,:) = low_easy_l(j,:);
            low_mid_l_all(j).matrix(count,:) = low_mid_l(j,:);
            low_dif_l_all(j).matrix(count,:) = low_dif_l(j,:);
            high_easy_l_all(j).matrix(count,:) = high_easy_l(j,:);
            high_mid_l_all(j).matrix(count,:) = high_mid_l(j,:);
            high_dif_l_all(j).matrix(count,:) = high_dif_l(j,:);
            
            correct_each_s_all(j).matrix(count,:) = correct_each_s(j,:);
            likeli_each_s_all(j).matrix(count,:) = likeli_each_s(j,:);
            correct_each_l_all(j).matrix(count,:) = correct_each_l(j,:);
            likeli_each_l_all(j).matrix(count,:) = likeli_each_l(j,:);
        end
        
        for j = 1:size_frame
            opt_s_all(j).matrix(count,:) = opt_s(j,:);
            opt_l_all(j).matrix(count,:) = opt_l(j,:);
            opt_ls_all(j).matrix(count,:) = opt_ls(j,:);
        end
        
        nan_data(count_neuron,1) = 0;
    else
        nan_data(count_neuron,1) = 1;
    end
    end
end

% 
% %More easy way is to analyze the correct rate with last-tone neurons and 
% %first-tone neurons. and see how much the correct rate deteriorate with the
% %neuron group!!

return

