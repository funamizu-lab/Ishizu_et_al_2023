
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
function FigureS8a_process_20230704_simple_process_test4_depth_control(folders)

close all
analysis_dir = eval(folders);

dir_filename = 'Tokyo2_Correct_20230701_prior_conf_activity*';

hist_count = 6;
long_use_timing = [15, 17, 25]; %before, first and end of sound

for i = 1:length(long_use_timing)
    L_prior(i).matrix = [];
    L_posterior(i).matrix = [];
    L_conf(i).matrix = [];
    L_r(i).matrix = [];
    L_partial_r(i).matrix = [];
end

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
    for j = 1:6
        evi_trace_all(i,j).matrix = [];
        evi_prefer_all(i,j).matrix = [];
        evi_nonprefer_all(i,j).matrix = [];
    end
end

for i = 1:length(analysis_dir)
    [i,length(analysis_dir)]
    
    [L_time_prior,L_time_posterior,L_time_conf,L_length_trial,~,L_r_prior,L_r_partial_prior] = ...
        Task_kaiseki_tokyo1_20230704_sound_choice_process4_depth(analysis_dir{i}, hist_count, dir_filename);
    
    for j = 1:length(long_use_timing)
        L_prior(j).matrix = [L_prior(j).matrix; L_time_prior(j).matrix];
        L_posterior(j).matrix = [L_posterior(j).matrix; L_time_posterior(j).matrix];
        L_conf(j).matrix = [L_conf(j).matrix; L_time_conf(j).matrix];
        L_r(j).matrix = [L_r(j).matrix; L_r_prior(j).matrix];
        L_partial_r(j).matrix = [L_partial_r(j).matrix; L_r_partial_prior(j).matrix];
    end
    
    L_prior_trial(i,:) = L_length_trial(:,1)';
    L_posterior_trial(i,:) = L_length_trial(:,2)';
    L_conf_trial(i,:) = L_length_trial(:,3)';
    
    %Add previous neurons
    [prop_sig(i,:), number_sig(i,:), sig_correct,~,~,~,~,~,~,~,length_neuron(i), ...
        evi_trace, evi_prefer, evi_nonprefer,~,~,~, sig_sin_session(i,:),s_sig_trace_prior] = ...
        Task_kaiseki_tokyo1_20220516_sound_choice_process4_depth(analysis_dir{i});
    for j = 1:40
        for k = 1:6
            evi_trace_all(j,k).matrix = [evi_trace_all(j,k).matrix; evi_trace(j,k).matrix];
            evi_prefer_all(j,k).matrix = [evi_prefer_all(j,k).matrix; evi_prefer(j,k).matrix];
            evi_nonprefer_all(j,k).matrix = [evi_nonprefer_all(j,k).matrix; evi_nonprefer(j,k).matrix];
        end
    end
end
% 
% %Get the significant neurons
temp_x = [0 0.25 0.45 0.55 0.75 1];
prior_frame = 15;
short_frame = 17;
long_frame = 25;
use_neuron(1).matrix = plot_prior_block_trace_20230704(prior_frame,evi_prefer_all,evi_nonprefer_all,temp_x);
use_neuron(2).matrix = plot_prior_block_trace_20230704(short_frame,evi_prefer_all,evi_nonprefer_all,temp_x);
use_neuron(3).matrix = plot_prior_block_trace_20230704(long_frame,evi_prefer_all,evi_nonprefer_all,temp_x);


%Adjust the number of neuron
for i = 1:3
    L_prior(i).matrix = L_prior(i).matrix(use_neuron(i).matrix,:);
    L_posterior(i).matrix = L_posterior(i).matrix(use_neuron(i).matrix,:);
    L_conf(i).matrix = L_conf(i).matrix(use_neuron(i).matrix,:);
    L_r(i).matrix = L_r(i).matrix(use_neuron(i).matrix,:);
    L_partial_r(i).matrix = L_partial_r(i).matrix(use_neuron(i).matrix,:);
end

%% figure i
cd('G:\upload_code\FigureS8\i');
region_num={'3','4','S8'};
name = cell2mat(region_num(ismember({'mpfc','auc','fof'},extractBefore(folders,'_'))));

size_bin = size(L_conf(1).matrix,2);
figure
sdata = struct();% source data
sdata.x = temp_x';% source data
tag={'left','middle','right'};
for i = 1:3
    subplot(1,3,i)
    hold on
    [median_trace,~,se_trace]= plot_median_se_each_nan(L_prior(i).matrix,hist_count,[0 0 0],2);
    eval(['sdata.',tag{i},'_priorval_median=transpose(median_trace)']);
    eval(['sdata.',tag{i},'_priorval_se=transpose(se_trace)']);
    [median_trace,~,se_trace]= plot_median_se_each_nan(L_conf(i).matrix,hist_count,[0 154 205]./255,2);
    eval(['sdata.',tag{i},'_actionval_median=transpose(median_trace)']);
    eval(['sdata.',tag{i},'_actionval_se=transpose(se_trace)']);
    set(gca,'xlim',[0 size_bin+1])
end
T = struct2table(sdata);
writetable(T, ['source fig',name,'i.csv']);
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mean_trace,std_trace,se_trace]= plot_median_se_each_nan(prior,hist_count,trace_color,std_se)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

length_trace = length(prior(:,1));
for i = 1:hist_count
    temp = prior(:,i);
    test1 = find(isnan(temp) == 1); %nan
    temp(test1) = [];
    mean_trace(i) = median(temp);
    temp_mean_trace = repmat(mean_trace(i),length_trace-length(test1),1);
    std_trace(i)  = median(abs(temp-temp_mean_trace),1);
    se_trace(i)   = 1.4826 .* std_trace(i) ./ (sqrt(length_trace));
    length_nan(i) = length(test1);
end

length_nan = [length_nan, length_trace];

std_plus  = mean_trace + std_trace;
std_minus = mean_trace - std_trace;
se_plus  = mean_trace + se_trace;
se_minus = mean_trace - se_trace;

if std_se == 0 %mean only
    plot(mean_trace,'color',trace_color,'LineWidth',1)
    box off
    
elseif std_se == 1 %std
    plot(mean_trace,'color',trace_color,'LineWidth',1)
    % box off
    for i = 1:hist_count
        hold on
        plot([i,i],[std_plus(i),std_minus(i)],'color',trace_color)
    end
    
elseif std_se == 2 %se
    plot(mean_trace,'color',trace_color,'LineWidth',1)
    for i = 1:hist_count
        hold on
        plot([i,i],[se_plus(i),se_minus(i)],'color',trace_color)
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function use_neuron = plot_prior_block_trace_20230704(prior_frame,evi_prefer_all,evi_nonprefer_all,temp_x)
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
return
