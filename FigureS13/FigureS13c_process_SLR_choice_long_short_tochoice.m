
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
function  FigureS13c_process_SLR_choice_long_short_tochoice(folders)

close all
[analysis_dir,short_length] = eval(folders);

SLR_choice = 'SLR_20230920_glmnet_choice_depth.mat';
thre_neuron = 20;
for i = 1:length(short_length)
    use_frame(i,:) = [1,2,4]; %Before sound, init sound, end sound
    
    frame_long(i,:)  = 1:6;
    frame_short(i,:) = 1:6;
end

[correct_each_s_choice, correct_each_l_choice, use_session] = ...
    get_average_correct_rate3(analysis_dir, SLR_choice, use_frame, 0, thre_neuron);


frame_long = frame_long(use_session,:);
frame_short = frame_short(use_session,:);

size_mouse = size(frame_long,1);
for i = 1:size_mouse
    for j = 1:6
        %correct rate for each tone difficulties
        each_l_choice(j).matrix(i,:) = correct_each_l_choice(frame_long(i,j)).matrix(i,:);
        each_s_choice(j).matrix(i,:) = correct_each_s_choice(frame_short(i,j)).matrix(i,:);
        
    end
end

%Flip the each correct rate
for i = 1:6
    each_l_choice(i).matrix = fliplr(each_l_choice(i).matrix);
    each_s_choice(i).matrix = fliplr(each_s_choice(i).matrix);
end

%%% FigS13c right 3 panels %%%
figure
subplot(1,3,1);hold on
plot_mean_se_moto(each_s_choice(2).matrix,[0 0 1],2) 
plot_mean_se_moto(each_l_choice(2).matrix,[1 0 0],2)
set(gca,'xlim',[0.5 3.5])
set(gca,'ylim',[0.65 1])

subplot(1,3,2);hold on
plot_mean_se_moto(each_s_choice(3).matrix,[0 0 1],2) 
plot_mean_se_moto(each_l_choice(3).matrix,[1 0 0],2)
set(gca,'xlim',[0.5 3.5])
set(gca,'ylim',[0.65 1])

subplot(1,3,3);hold on
plot_mean_se_moto(each_s_choice(4).matrix,[0 0 1],2) 
plot_mean_se_moto(each_l_choice(4).matrix,[1 0 0],2)
set(gca,'xlim',[0.5 3.5])
set(gca,'ylim',[0.65 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [correct_each_s_all, correct_each_l_all, use_session] = ...
    get_average_correct_rate3(analysis_dir, SLR_name, use_frame, psycho_sign, thre_neuron)
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
            
            %Use frame is for psychometric function only
            [~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,correct_each_s,~,correct_each_l] = ...
                Population_decoder_20220518_SLR_process(analysis_dir{i},SLR_name, use_frame(i,:), psycho_sign);
            
            for j = 1:8
                
                correct_each_s_all(j).matrix(count,:) = correct_each_s(j,:);
                correct_each_l_all(j).matrix(count,:) = correct_each_l(j,:);
            end
            
            
            use_session(count) = i;
        else
        end
    end
end
return

