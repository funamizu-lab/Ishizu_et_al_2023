
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
function FigureS14e_process_SLR_long_detail_limit_neuron2(folders)

close all
[analysis_dir,short_length] = eval(folders);

SLR_choice = 'SLR_20220520_glmnet_choice_depth.mat';

thre_neuron = 20;

for i = 1:length(short_length)
    use_frame(i,:) = [2,4];
end

%%% figS13e right 2 panels %%%
name = extractBefore(folders,'_');
make_psychometric_function(analysis_dir,SLR_choice, use_frame, thre_neuron, name);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function make_psychometric_function(analysis_dir,SLR_name, use_frame, thre_neuron, name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[opt_l_R_all, opt_l_L_all] = ...
    get_average_correct_rate(analysis_dir, SLR_name, use_frame, 1, thre_neuron);

cd('G:\upload_code\FigureS14\FigS14e');
[~,size_frame] = size(use_frame);
sdata = struct();% source data
evi_x = 0:0.01:1;
sdata.x=evi_x';
tag={'initial sound','late sound'};
for j = 1:size_frame
    temp_R = opt_l_R_all(j).matrix;
    temp_L = opt_l_L_all(j).matrix;
    
    figure;hold on
    [mean_trace_R,~,se_trace_R] =plot_mean_se_moto_x_axis(temp_R,evi_x,[1 0 0],2);
    [mean_trace_L,~,se_trace_L] =plot_mean_se_moto_x_axis(temp_L,evi_x,[0 0 1],2);
    set(gca,'xlim',[-0.1 1.1],'ylim',[0 1])
    
    %%% source data%%%
    sdata.Right=mean_trace_R';
    sdata.Right_se=se_trace_R';
    sdata.Left=mean_trace_L';
    sdata.Left_se=se_trace_L';
    T = struct2table(sdata);
    writetable(T, ['source figS14e right ',name,' ',tag{j},'.csv']);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [opt_l_R_all, opt_l_L_all] = ...
    get_average_correct_rate(analysis_dir, SLR_name, use_frame, psycho_sign, thre_neuron)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,size_frame] = size(use_frame);

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
        if length(lasso_sound_l) ~= 1 %not nan
            count = count + 1;
            [~,opt_l_R,opt_l_L] = ...
                Population_decoder_20220121_SLR_process2_long_detail(analysis_dir{i}, SLR_name, use_frame(i,:), psycho_sign);
            
            for j = 1:size_frame
                opt_l_R_all(j).matrix(count,:) = opt_l_R(j,:);
                opt_l_L_all(j).matrix(count,:) = opt_l_L(j,:);
            end
        end
    end
end

return