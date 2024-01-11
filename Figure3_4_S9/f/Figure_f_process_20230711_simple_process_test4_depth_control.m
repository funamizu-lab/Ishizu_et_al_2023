
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
function Figure_f_process_20230711_simple_process_test4_depth_control(folders)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Determine the sig_neuron at certain time window
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sig_time_window = 15; %Before sound neurons 
%sig_time_window = 25; %Only during sound neurons
sig_time_window = 0; %Before + During sound neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

analysis_dir = eval(folders);

use_frame = 6:40;

sig_start_all = [];
choice_prior_all = [];
check_overlap = [];
for i = 1:40
    sig_correct_all(i).matrix = [];
    sig_error_all(i).matrix = [];
    sig_sound_all(i).matrix = [];
    sig_choice_all(i).matrix = [];
    sig_prior_all(i).matrix = [];
    sig_sin_all(i).matrix = [];
    sig_choice1_all(i).matrix = [];
    sig_choice2_all(i).matrix = [];
    sig_prior1_all(i).matrix = [];
    sig_prior2_all(i).matrix = [];
    sig_trace_prior_all(i).matrix = [];
    sig_trace_nonprior_all(i).matrix = [];
    sig_anova_all(i).matrix = [];
    for j = 1:6
        evi_trace_all(i,j).matrix = [];
        evi_prefer_all(i,j).matrix = [];
        evi_nonprefer_all(i,j).matrix = [];
    end
end
for i = 1:32
    short_correct_all(i).matrix = [];
    s_trace_prior_all(i).matrix = [];
    s_trace_nonprior_all(i).matrix = [];
    s_sig_anova_all(i).matrix = [];
    for j = 1:6
        s_evi_trace_all(i,j).matrix = [];
        s_evi_prefer_all(i,j).matrix = [];
        s_evi_nonprefer_all(i,j).matrix = [];
    end
end
for i = 1:length(analysis_dir)
    [i,length(analysis_dir)]
   
    [prop_sig(i,:), number_sig(i,:), sig_correct, sig_error, sig_sound, sig_choice, sig_prior, sig_sin, ...
        choice_prior_matrix, sig_start, length_neuron(i), ...
        evi_trace, evi_prefer, evi_nonprefer, sig_trace_prior, sig_trace_nonprior, ...
        sig_correct_short, sig_sin_session(i,:), ...
        s_sig_trace_prior,s_sig_trace_nonprior,...
        s_evi_trace,s_evi_prefer,s_evi_nonprefer, check_neuron, sig_anova_long, sig_anova_short] = ...
        Task_kaiseki_tokyo1_20230711_sound_choice_process4A_depth(analysis_dir{i},sig_time_window);
    
    choice_prior_all = [choice_prior_all; choice_prior_matrix];
    sig_start_all = [sig_start_all; sig_start];
    check_overlap = [check_overlap; check_neuron];
    for j = 1:40
        sig_correct_all(j).matrix = [sig_correct_all(j).matrix; sig_correct(j).matrix];
        sig_error_all(j).matrix = [sig_error_all(j).matrix; sig_error(j).matrix];
        sig_sound_all(j).matrix = [sig_sound_all(j).matrix; sig_sound(j).matrix];
        sig_choice_all(j).matrix = [sig_choice_all(j).matrix; sig_choice(j).matrix];
        sig_prior_all(j).matrix = [sig_prior_all(j).matrix; sig_prior(j).matrix];
        sig_sin_all(j).matrix = [sig_sin_all(j).matrix; sig_sin(j).matrix];
        sig_anova_all(j).matrix = [sig_anova_all(j).matrix; sig_anova_long(j).matrix];     
        sig_trace_prior_all(j).matrix = [sig_trace_prior_all(j).matrix; sig_trace_prior(j).matrix]; %mean activity of block
        sig_trace_nonprior_all(j).matrix = [sig_trace_nonprior_all(j).matrix; sig_trace_nonprior(j).matrix];
        
        for k = 1:6
            evi_trace_all(j,k).matrix = [evi_trace_all(j,k).matrix; evi_trace(j,k).matrix];
            evi_prefer_all(j,k).matrix = [evi_prefer_all(j,k).matrix; evi_prefer(j,k).matrix];  %Correct trials with evidence prefer block
            evi_nonprefer_all(j,k).matrix = [evi_nonprefer_all(j,k).matrix; evi_nonprefer(j,k).matrix]; %Correct trials with evidence nonprefer block
        end
    end
    
    %for short sound
    for j = 1:32
        short_correct_all(j).matrix = [short_correct_all(j).matrix; sig_correct_short(j).matrix];
        s_trace_prior_all(j).matrix = [s_trace_prior_all(j).matrix; s_sig_trace_prior(j).matrix];
        s_trace_nonprior_all(j).matrix = [s_trace_nonprior_all(j).matrix; s_sig_trace_nonprior(j).matrix];
        s_sig_anova_all(j).matrix = [s_sig_anova_all(j).matrix; sig_anova_short(j).matrix];
        for k = 1:6
            s_evi_trace_all(j,k).matrix = [s_evi_trace_all(j,k).matrix; s_evi_trace(j,k).matrix];
            s_evi_prefer_all(j,k).matrix = [s_evi_prefer_all(j,k).matrix; s_evi_prefer(j,k).matrix];
            s_evi_nonprefer_all(j,k).matrix = [s_evi_nonprefer_all(j,k).matrix; s_evi_nonprefer(j,k).matrix];
        end
    end
end

for j = 1:length(use_frame)
    temp_frame = use_frame(j);
    [median_sig_prior(j),std_sig_prior(j),se_sig_prior(j)] = ...
        matrix2medians(sig_trace_prior_all(temp_frame).matrix);
    [median_sig_nonprior(j),std_sig_nonprior(j),se_sig_nonprior(j)] = ...
        matrix2medians(sig_trace_nonprior_all(temp_frame).matrix);
    
    %Get the nan frame
    for k = 1:6
        moto_data = evi_prefer_all(temp_frame,k).matrix;
        nan_check_prefer(:,k) = isnan(moto_data); %detect_non_nan
        moto_data = evi_nonprefer_all(temp_frame,k).matrix;
        nan_check_nonprefer(:,k) = isnan(moto_data); %detect_non_nan
    end
    nan_check_prefer = max(nan_check_prefer,[],2);
    nan_check_nonprefer = max(nan_check_nonprefer,[],2);
    nan_check = max([nan_check_prefer,nan_check_nonprefer],[],2);
    use_neuron = find(nan_check == 0);
    
    plot_sig_prior(j).matrix = sig_trace_prior_all(temp_frame).matrix(use_neuron,:); %Mean activity for each block all trials
    plot_sig_nonprior(j).matrix = sig_trace_nonprior_all(temp_frame).matrix(use_neuron,:);

    plot_sabun(j).matrix = plot_sig_prior(j).matrix - plot_sig_nonprior(j).matrix;
    [size_sabun(j,1),size_sabun(j,2)] = size(plot_sabun(j).matrix);
    
    [evi_trace,evi_trace_neuron(j)] = get_non_nan_trace_new(evi_trace_all,temp_frame,use_neuron);
    [evi_prefer_trace,evi_prefer_trace_neuron(j)] = get_non_nan_trace_new(evi_prefer_all,temp_frame,use_neuron);
    [evi_nonprefer_trace,evi_nonprefer_trace_neuron(j)] = get_non_nan_trace_new(evi_nonprefer_all,temp_frame,use_neuron);
    for k = 1:6
        plot_evi_trace(k).matrix(j).matrix = evi_trace(k).matrix;
        plot_evi_prefer(k).matrix(j).matrix = evi_prefer_trace(k).matrix;
        plot_evi_nonprefer(k).matrix(j).matrix = evi_nonprefer_trace(k).matrix;
    end
end

temp_x = [0 0.25 0.45 0.55 0.75 1];

%Make sabun trace based on the correct trials with all evidence 
for i = 1:length(use_frame)
    [temp_p,temp_n,temp_prefer_length(i),temp_p_anovan(i,:),temp_anova]=...
        plot_prior_block_trace2(use_frame(i),evi_prefer_all,evi_nonprefer_all,temp_x,sig_anova_all);

    temp_p_anova = zeros(length(temp_anova(:,1)),1);
    temp_p_anova(temp_anova(:,1) < 0.01) = 1;
    p_anova_trace(:,i) = temp_p_anova;
    
    temp_prior = mean(temp_p,2) - mean(temp_n,2);
    sabun_prior_trace(:,i) = temp_prior;
    p_prior_trace(:,i) = mean(temp_p,2);
    n_prior_trace(:,i) = mean(temp_n,2);
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
clear sig_prop_neuron
for i = 1:length(use_frame)
    temp1 = find(sig_sort(:,i) == 1);
    temp2 = find(sig_sort(:,i) == -1);
    sig_prop_neuron(1,i) = length(temp1) ./ length(sig_sort(:,i));
    sig_prop_neuron(2,i) = length(temp2) ./ length(sig_sort(:,i));
end

%% fig f left
cd('G:\upload_code\Figure3_4_S8\f');
region_num={'3','4','S8'};
name = cell2mat(region_num(ismember({'mpfc','auc','fof'},extractBefore(folders,'_'))));
use_color = jet(6);
figure; hold on
sdata = struct();% source data 
tag={'prefer100','prefer75','prefer55','prefer45','prefer25','prefer0'};
for k = 1:6
    [~,mean_trace,~,se_trace] = ...
        plot_mean_se_moto_x_axis_matrix(plot_evi_trace(k).matrix,1:length(use_frame),use_color(k,:),2);    
    eval(['sdata.',tag{end-k+1},'_mean=transpose(mean_trace)']);
    eval(['sdata.',tag{end-k+1},'_se=transpose(se_trace)']);
end
set(gca,'xlim',[0 length(use_frame)])
set(gca,'xtick',0:2:length(use_frame))
T = struct2table(sdata);
writetable(T, ['source fig',name,'f left.csv']);

%% fig f reight top
figure; hold on
sdata = struct();% source data 
[mean_trace1,~,se_trace1]=plot_mean_se_moto_x_axis(p_prior_trace,1:length(use_frame),[1 0 0],2);
[mean_trace2,~,se_trace2]=plot_mean_se_moto_x_axis(n_prior_trace,1:length(use_frame),[0 0 1],2);
sdata.prefer_mean= mean_trace1';
sdata.prefer_se  = se_trace1'; 
sdata.nonprefer_mean= mean_trace2';
sdata.nonprefer_se  = se_trace2'; 
T = struct2table(sdata);
writetable(T, ['source fig',name,'f right top.csv']);

%% fig f reight bottom
figure; hold on;
sdata = struct();% source data 
plot(sig_prop_neuron(1,:),'r')
plot(sig_prop_neuron(2,:),'b')
sdata.prefer= sig_prop_neuron(1,:)';
sdata.nonprefer  = sig_prop_neuron(2,:)'; 
T = struct2table(sdata);
writetable(T, ['source fig',name,'f right bottom.csv']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [trace,length_neuron] = get_non_nan_trace_new(evi_trace_all,use_frame,use_neuron)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

sig_anova = sig_anova(prior_frame).matrix;
sig_anova = sig_anova(use_neuron,:);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_tuning_curve(temp_x, median_evi_trace, se_evi_trace, use_color)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hold on
plot(temp_x,median_evi_trace,'color',use_color) %sound start
for i = 1:6
    plot([temp_x(i),temp_x(i)],[median_evi_trace(i)+se_evi_trace(i),median_evi_trace(i)-se_evi_trace(i)],'color',use_color)
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
