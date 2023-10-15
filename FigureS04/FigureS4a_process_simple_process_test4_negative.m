
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
function FigureS4a_process_simple_process_test4_negative(folders)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Determine the sig_neuron at certain time window
sig_time_window = 0; %Before + During sound neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
analysis_dir = eval(folders);

use_frame = 6:40;

for i = 1:40
    for j = 1:6
        evi_trace_all(i,j).matrix = [];
        evi_prefer_all(i,j).matrix = [];
        evi_nonprefer_all(i,j).matrix = [];
    end
end
for i = 1:length(analysis_dir)
    [i,length(analysis_dir)]
   
    [~,~,~,~,~,~,~,~,~,evi_trace,evi_prefer,evi_nonprefer] = ...
        Task_kaiseki_tokyo1_20230711_sound_choice_process4A_nega(analysis_dir{i},sig_time_window);
    
    for j = 1:40
        for k = 1:6
            evi_trace_all(j,k).matrix = [evi_trace_all(j,k).matrix; evi_trace(j,k).matrix];
            evi_prefer_all(j,k).matrix = [evi_prefer_all(j,k).matrix; evi_prefer(j,k).matrix];  %Correct trials with evidence prefer block
            evi_nonprefer_all(j,k).matrix = [evi_nonprefer_all(j,k).matrix; evi_nonprefer(j,k).matrix]; %Correct trials with evidence nonprefer block
        end
    end
end

for j = 1:length(use_frame)
    temp_frame = use_frame(j);
    %Get the nan frame
    for k = 1:6
        moto_data = evi_prefer_all(temp_frame,k).matrix;
        nan_check_prefer(:,k) = isnan(moto_data); %detect_non_nan
        moto_data = evi_nonprefer_all(temp_frame,k).matrix;
        nan_check_nonprefer(:,k) = isnan(moto_data); %detect_non_nan
    end
    %nan_check = max(nan_check,[],2);
    nan_check_prefer = max(nan_check_prefer,[],2);
    nan_check_nonprefer = max(nan_check_nonprefer,[],2);
    nan_check = max([nan_check_prefer,nan_check_nonprefer],[],2);
    use_neuron = find(nan_check == 0);
    
    [evi_trace,evi_trace_neuron(j)] = get_non_nan_trace_new(evi_trace_all,temp_frame,use_neuron);
    for k = 1:6
        plot_evi_trace(k).matrix(j).matrix = evi_trace(k).matrix;
    end
end

%% Fig S4a
use_color = jet(6);
figure; hold on
for k = 1:6
    plot_mean_se_moto_x_axis_matrix(plot_evi_trace(k).matrix,1:length(use_frame),use_color(k,:),2);
end
set(gca,'xlim',[0 length(use_frame)])
set(gca,'xtick',0:2:length(use_frame))




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
function [prefer_activ,nonprefer_activ,length_neuron,p_anovan,sig_anova] = plot_prior_block_trace2(prior_frame,evi_prefer_all,evi_nonprefer_all,temp_x,sig_anova)
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

%for anavan
ano_activ = [];
ano_sound = [];
ano_block = [];
for k = 1:6
    moto_data = evi_prefer_all(prior_frame,k).matrix;
    moto_data = moto_data(use_neuron);
    median_prefer(k) = median(moto_data);
    std_prefer(k) = median(abs(moto_data-median_prefer(k)));
    se_prefer(k) = 1.4826 * std_prefer(k) ./ (sqrt(length(use_neuron)));

    prefer_activ(:,k) = moto_data;
    ano_activ = [ano_activ;moto_data];
    sound0 = ones(length(moto_data),1).*k;
    block0 = zeros(length(moto_data),1);
        
    moto_data = evi_nonprefer_all(prior_frame,k).matrix;
    moto_data = moto_data(use_neuron);
    median_nonprefer(k) = median(moto_data);
    std_nonprefer(k) = median(abs(moto_data-median_nonprefer(k)));
    se_nonprefer(k) = 1.4826 * std_nonprefer(k) ./ (sqrt(length(use_neuron)));
    
    nonprefer_activ(:,k) = moto_data;
    ano_activ = [ano_activ;moto_data];
    sound1 = ones(length(moto_data),1).*k;
    block1 = ones(length(moto_data),1);

    ano_sound = [ano_sound;sound0;sound1];
    ano_block = [ano_block;block0;block1];
end
length_neuron = length(use_neuron);
p_anovan = anovan(ano_activ,{ano_sound,ano_block},'display','off');
    
plot_tuning_curve(temp_x, median_prefer, se_prefer, [1 0 0])
hold on
plot_tuning_curve(temp_x, median_nonprefer, se_nonprefer, [0 0 1])
set(gca,'xlim',[-0.1 1.1])

sig_anova = sig_anova(prior_frame).matrix;
sig_anova = sig_anova(use_neuron,:);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_tuning_curve(temp_x, median_evi_trace, se_evi_trace, use_color)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(temp_x,median_evi_trace,'color',use_color) %sound start
hold on
for i = 1:6
    plot([temp_x(i),temp_x(i)],[median_evi_trace(i)+se_evi_trace(i),median_evi_trace(i)-se_evi_trace(i)],'color',use_color)
    hold on
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [median_trace,std_trace,se_trace] = matrix2medians(moto_data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Extract nan trial
moto_data = moto_data(isnan(moto_data) == 0);

median_trace = median(moto_data);
std_trace  = median(abs(moto_data-median_trace));
se_trace = 1.4826 * std_trace ./ (sqrt(length(moto_data)));

return
