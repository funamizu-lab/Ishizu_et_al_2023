
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
function Figure_i_process_20240102_simple_process_test4_depth_control(folders)

savefolder = 'G:\upload_code\Figure5_S9\i';
analysis_dir = eval(folders);

dir_filename = 'Tokyo2_Correct_20230701_prior_conf_activity*';

hist_count = 6;
long_use_timing = [15, 17, 25]; %before, first and end of sound

for i = 1:length(long_use_timing)
    L_prior(i).matrix = [];
    L_add(i).matrix = [];
end

for i = 1:40
    for j = 1:6
        evi_prefer_all(i,j).matrix = [];
        evi_nonprefer_all(i,j).matrix = [];
    end
end
for i = 1:length(analysis_dir)
    [i,length(analysis_dir)]
    
    [L_time_prior,~,~,~,~,~,~,L_time_add] = ...
        Task_kaiseki_tokyo1_20240102_sound_choice_process4_depth(analysis_dir{i}, hist_count, dir_filename);

    for j = 1:length(long_use_timing)
        L_prior(j).matrix = [L_prior(j).matrix; L_time_prior(j).matrix];
        L_add(j).matrix = [L_add(j).matrix; L_time_add(j).matrix];
    end
        
    %Add previous neurons
    [~,~,~,~,~,~,~,~,~,~,~,~, evi_prefer, evi_nonprefer] = ...
        Task_kaiseki_tokyo1_20220516_sound_choice_process4_depth(analysis_dir{i});
    for j = 1:40
        for k = 1:6
            evi_prefer_all(j,k).matrix = [evi_prefer_all(j,k).matrix; evi_prefer(j,k).matrix];
            evi_nonprefer_all(j,k).matrix = [evi_nonprefer_all(j,k).matrix; evi_nonprefer(j,k).matrix];
        end
    end
end
delete(gcp('nocreate'))

% Get the significant neurons
prior_frame = 15;
short_frame = 17;
long_frame = 25;
use_neuron(1).matrix = plot_prior_block_trace_20230704(prior_frame,evi_prefer_all,evi_nonprefer_all);
use_neuron(2).matrix = plot_prior_block_trace_20230704(short_frame,evi_prefer_all,evi_nonprefer_all);
use_neuron(3).matrix = plot_prior_block_trace_20230704(long_frame,evi_prefer_all,evi_nonprefer_all);

%Adjust the number of neuron
for i = 1:3
    L_prior(i).matrix = L_prior(i).matrix(use_neuron(i).matrix,:);
    L_add(i).matrix = L_add(i).matrix(use_neuron(i).matrix,:);
end

%%% fig i %%%
label = {'left','mid','right'};
cd(savefolder);
sdata = plot_Q_traces_all(L_prior, L_add, hist_count,label);
T = struct2table(sdata);
writetable(T, [folders(1:3), ' source fig i.csv']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function use_neuron = plot_prior_block_trace_20230704(prior_frame,evi_prefer_all,evi_nonprefer_all)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sdata=plot_Q_traces_all(L_prior, L_add, hist_count,label)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = (1:hist_count)/(hist_count+1);

sdata=struct();
sdata.x=x';
figure
for i = 1:3
    subplot(1,3,i)   
    hold on
    [med_prior,se_prior]=plot_median_se_each_nan(x,L_prior(i).matrix,hist_count,[0 0 0]);
    [med_add,  se_add]=plot_median_se_each_nan(x,L_add(i).matrix,hist_count,[0 154 205]./255);
    set(gca,'xlim',[0 1])
    eval(['sdata.prior_median_',label{i},'=transpose(med_prior);']);
    eval(['sdata.prior_se_',label{i},'=transpose(se_prior);']);
    eval(['sdata.additive_median_',label{i},'=transpose(med_add);']);
    eval(['sdata.additive_se_',label{i},'=transpose(se_add);']);
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [median_trace,se_trace] = plot_median_se_each_nan(x,data,hist_count,trace_color)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

length_trace = length(data(:,1));
for i = 1:hist_count
    temp = data(:,i);
    test1 = find(isnan(temp) == 1); %nan
    temp(test1) = [];
    median_trace(i) = median(temp);
    temp_median_trace = repmat(median_trace(i),length_trace-length(test1),1);
    std_trace(i)  = median(abs(temp-temp_median_trace),1);
    se_trace(i)   = 1.4826 .* std_trace(i) ./ (sqrt(length_trace));
end

se_plus  = median_trace + se_trace;
se_minus = median_trace - se_trace;

plot(x,median_trace,'color',trace_color,'LineWidth',1)
for i = 1:hist_count
    plot([x(i),x(i)],[se_plus(i),se_minus(i)],'color',trace_color)
end
