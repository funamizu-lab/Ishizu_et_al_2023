
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
function FigureS5c_simple_process_test4_depth_control(folders)

close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Determine the sig_neuron at certain time window
sig_time_window = 0; %Before + During sound neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

analysis_dir = eval(folders);

use_frame = 6:40;
for i = 1:40
    sig_trace_prior_all(i).matrix = [];
    sig_trace_nonprior_all(i).matrix = [];
    sig_anova_all(i).matrix = [];
    for j = 1:6
        evi_trace_all(i,j).matrix = [];
        evi_prefer_all(i,j).matrix = [];
        evi_nonprefer_all(i,j).matrix = [];
    end
end
for i = 1:length(analysis_dir)
    [i,length(analysis_dir)]
   
    [~,~,~,~,~,~,~,~,~,~,~,evi_trace, evi_prefer, evi_nonprefer, sig_trace_prior, sig_trace_nonprior, ...
        ~,~,~,~,~,~,~,~,sig_anova_long] = ...
        Task_kaiseki_tokyo1_20230711_sound_choice_process4A_depth(analysis_dir{i},sig_time_window);
    

    for j = 1:40
        sig_anova_all(j).matrix = [sig_anova_all(j).matrix; sig_anova_long(j).matrix];     
        sig_trace_prior_all(j).matrix = [sig_trace_prior_all(j).matrix; sig_trace_prior(j).matrix]; %mean activity of block
        sig_trace_nonprior_all(j).matrix = [sig_trace_nonprior_all(j).matrix; sig_trace_nonprior(j).matrix];
        
        for k = 1:6
            evi_trace_all(j,k).matrix = [evi_trace_all(j,k).matrix; evi_trace(j,k).matrix];
            evi_prefer_all(j,k).matrix = [evi_prefer_all(j,k).matrix; evi_prefer(j,k).matrix];  %Correct trials with evidence prefer block
            evi_nonprefer_all(j,k).matrix = [evi_nonprefer_all(j,k).matrix; evi_nonprefer(j,k).matrix]; %Correct trials with evidence nonprefer block
        end
    end
    
end

for i = 1:length(use_frame)
    [temp_p,temp_n,temp_anova] = plot_prior_block_trace2(use_frame(i),evi_prefer_all,evi_nonprefer_all,sig_anova_all);

    temp_p_anova = zeros(length(temp_anova(:,1)),1);
    temp_p_anova(temp_anova(:,1) < 0.01) = 1;
    p_anova_trace(:,i) = temp_p_anova;
    
    temp_prior = mean(temp_p,2) - mean(temp_n,2);
    sabun_prior_trace(:,i) = temp_prior;
end
%make sort averaging
sort_sabun = sabun_prior_trace;
ave_sabun = mean(sort_sabun,2);
[~,ave_sabun] = sort(ave_sabun,'descend');
sort_sabun = sort_sabun(ave_sabun,:);
anova_sabun = p_anova_trace(ave_sabun,:);

sig_window = find(anova_sabun == 1);
sig_sort = zeros(size(sort_sabun));
sig_sort(sig_window) = sort_sabun(sig_window);
sig_sort(sig_sort > 0) = 1;
sig_sort(sig_sort < 0) = -1;

%%% FigureS5c left panels %%%
figure
imagesc(sort_sabun)
% make_red_blue_colormap_20230710(sort_sabun,color_step_size)

figure
imagesc(sig_sort)
% make_red_blue_colormap_20230710(sig_sort,color_step_size)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [trace,length_neuron] = get_non_nan_trace_new(evi_trace_all,use_frame,use_neuron)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get the nan frame
length_neuron = length(use_neuron);

for k = 1:6
    moto_data = evi_trace_all(use_frame,k).matrix;
    trace(k).matrix = moto_data(use_neuron,:);
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [prefer_activ,nonprefer_activ,sig_anova] = plot_prior_block_trace2(prior_frame,evi_prefer_all,evi_nonprefer_all,sig_anova)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get the nan frame
for k = 1:6
    moto_data = evi_prefer_all(prior_frame,k).matrix;
    nan_check1(:,k) = isnan(moto_data); %detect_non_nan
    moto_data = evi_nonprefer_all(prior_frame,k).matrix;
    nan_check2(:,k) = isnan(moto_data); %detect_non_nan
end
nan_check = [nan_check1,nan_check2];
nan_check = max(nan_check,[],2);
use_neuron = find(nan_check == 0);

for k = 1:6
    moto_data = evi_prefer_all(prior_frame,k).matrix;
    moto_data = moto_data(use_neuron);
    prefer_activ(:,k) = moto_data;
        
    moto_data = evi_nonprefer_all(prior_frame,k).matrix;
    moto_data = moto_data(use_neuron);
    nonprefer_activ(:,k) = moto_data;
end
sig_anova = sig_anova(prior_frame).matrix;
sig_anova = sig_anova(use_neuron,:);
return