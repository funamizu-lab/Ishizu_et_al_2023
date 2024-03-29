
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
function Figure7e_process_20220520_0912_SLR_long_detail_limit_neuron2(folders)

close all
[analysis_dir,short_length] = eval(folders);

SLR_sound = 'SLR_20220520_glmnet_sound_depth.mat';
% SLR_choice= 'SLR_20220520_glmnet_choice_depth.mat';

thre_neuron = 20;

for i = 1:length(short_length)
    use_frame(i,:) =[2,4];
end

%% Fig 6e right 2 panels
[sabun,mean_traceR,mean_traceL,se_traceR,se_traceL,trace_x]=...
    make_psychometric_function(analysis_dir,SLR_sound, use_frame, thre_neuron);

cd('G:\upload_code\Figure7\Fig7e');
name = extractBefore(folders,'_');
sdata = struct();% source data 
sdata.sabun = sabun;
T = struct2table(sdata);
writetable(T, ['source fig7e right ',name,'.csv']);

if(ismember({'auc'},name))
    sdata = struct();% source data
    sdata.x = trace_x';
    sdata.rightblock_mean= mean_traceR';
    sdata.rightblock_se  = se_traceR';
    sdata.leftblock_mean= mean_traceL';
    sdata.leftblock_se  = se_traceL';
    T = struct2table(sdata);
    writetable(T, ['source fig7e middle ',name,'.csv']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sabun,mean_traceR,mean_traceL,se_traceR,se_traceL,trace_x]=...
    make_psychometric_function(analysis_dir,SLR_name, use_frame, thre_neuron)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[opt_l_R_all,opt_l_L_all]=get_average_correct_rate(analysis_dir,SLR_name,use_frame,1,thre_neuron);

% [~,size_frame] = size(use_frame);

evi_x = 0:0.01:1;
%Opt_evi
% for j = 1:size_frame
j=2;
temp_R = opt_l_R_all(j).matrix;
temp_L = opt_l_L_all(j).matrix;

figure; hold on
[mean_traceR,~,se_traceR,trace_x] = plot_mean_se_moto_x_axis(temp_R,evi_x,[1 0 0],2);
[mean_traceL,~,se_traceL] = plot_mean_se_moto_x_axis(temp_L,evi_x,[0 0 1],2);
set(gca,'xlim',[-0.1 1.1],'ylim',[0 1])

%Plot difference with behavior
sabun = mean(temp_R,2) - mean(temp_L,2);
% end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [opt_l_R_all, opt_l_L_all] = ...
    get_average_correct_rate(analysis_dir,SLR_name,use_frame,psycho_sign,thre_neuron)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,size_frame] = size(use_frame);

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
    if length(lasso_sound_l) ~= 1 %not nan
        count = count + 1;%        
        [~,opt_l_R, opt_l_L] = ...
            Population_decoder_20220121_SLR_process2_long_detail(analysis_dir{i},SLR_name,use_frame(i,:),psycho_sign);

        for j = 1:size_frame
%             opt_l_all(j).matrix(count,:) = opt_l(j,:);
            opt_l_R_all(j).matrix(count,:) = opt_l_R(j,:);
            opt_l_L_all(j).matrix(count,:) = opt_l_L(j,:);
            
        end
    end
    end
end

return