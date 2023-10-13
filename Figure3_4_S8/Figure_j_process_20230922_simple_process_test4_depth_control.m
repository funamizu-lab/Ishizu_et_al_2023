
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
function Figure_j_process_20230922_simple_process_test4_depth_control(folders)

analysis_dir = eval(folders);

% sig_start_all = [];
% choice_prior_all = [];
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
    
    long_time_RMSE(i).matrix = [];
    for j = 1:6
        evi_trace_all(i,j).matrix = [];
        evi_prefer_all(i,j).matrix = [];
        evi_nonprefer_all(i,j).matrix = [];
    end
end
for i = 1:length(analysis_dir)
    [i,length(analysis_dir)]

    [prop_sig(i,:), number_sig(i,:), sig_correct, sig_error, sig_sound, sig_choice, sig_prior, sig_sin, ...
        choice_prior_matrix, sig_start, length_neuron(i), ...
        evi_trace, evi_prefer, evi_nonprefer, sig_trace_prior, sig_trace_nonprior, ...
        ~, sig_sin_session(i,:),s_sig_trace_prior,~,~,~,~,temp_RMSE] = ...
        Task_kaiseki_tokyo1_20230920_sound_choice_process4_depth(analysis_dir{i});    
    
%     choice_prior_all = [choice_prior_all; choice_prior_matrix];
%     sig_start_all = [sig_start_all; sig_start];
    for j = 1:40
        sig_correct_all(j).matrix = [sig_correct_all(j).matrix; sig_correct(j).matrix];
        sig_error_all(j).matrix = [sig_error_all(j).matrix; sig_error(j).matrix];
        sig_sound_all(j).matrix = [sig_sound_all(j).matrix; sig_sound(j).matrix];
        sig_choice_all(j).matrix = [sig_choice_all(j).matrix; sig_choice(j).matrix];
        sig_prior_all(j).matrix = [sig_prior_all(j).matrix; sig_prior(j).matrix];
        sig_sin_all(j).matrix = [sig_sin_all(j).matrix; sig_sin(j).matrix];
        
        sig_trace_prior_all(j).matrix = [sig_trace_prior_all(j).matrix; sig_trace_prior(j).matrix];
        sig_trace_nonprior_all(j).matrix = [sig_trace_nonprior_all(j).matrix; sig_trace_nonprior(j).matrix];
        
        long_time_RMSE(j).matrix = [long_time_RMSE(j).matrix, temp_RMSE(j).matrix];
        
        for k = 1:6
            evi_trace_all(j,k).matrix = [evi_trace_all(j,k).matrix; evi_trace(j,k).matrix];
            evi_prefer_all(j,k).matrix = [evi_prefer_all(j,k).matrix; evi_prefer(j,k).matrix];
            evi_nonprefer_all(j,k).matrix = [evi_nonprefer_all(j,k).matrix; evi_nonprefer(j,k).matrix];
        end
    end
end

temp_x = [0 0.25 0.45 0.55 0.75 1];

%Sound on: 16
%Sound end: 25
prior_frame = 15;
short_frame = 17;
long_frame  = 25;
use_neuron1 = plot_prior_block_trace2(prior_frame,evi_prefer_all,evi_nonprefer_all,temp_x);
use_neuron2 = plot_prior_block_trace2(short_frame,evi_prefer_all,evi_nonprefer_all,temp_x);
use_neuron3 = plot_prior_block_trace2(long_frame,evi_prefer_all,evi_nonprefer_all,temp_x);

%make boxplot for each RMSE
%choice, prior, conf_Q, binary_tone, correct_side
temp1 = long_time_RMSE(prior_frame).matrix';
temp2 = long_time_RMSE(short_frame).matrix';
temp3 = long_time_RMSE(long_frame).matrix';

%Update the number of neurons with use_neuron
temp1 = temp1(use_neuron1,:);
temp2 = temp2(use_neuron2,:);
temp3 = temp3(use_neuron3,:);

%16 and 17 is empty...
%Make the sabun plot
[~, p(1,:), min_number(1), min_value(1), ...
    temp1_WO_conf, p_WO_conf(1,:), min_WO(1), min_WO_value(1)] = analyze_regression_each_time(temp1);
[~, p(2,:), min_number(2), min_value(2), ...
    temp2_WO_conf, p_WO_conf(2,:), min_WO(2), min_WO_value(2)] = analyze_regression_each_time(temp2);
[~, p(3,:), min_number(3), min_value(3), ...
    temp3_WO_conf, p_WO_conf(3,:), min_WO(3), min_WO_value(3)] = analyze_regression_each_time(temp3);

%% Fig i
cd('G:\upload_code\Figure3_4_S8\j');
sdata = struct();% source data 

temp_pick = [1,2,3,6];%Pick only choice, prior sound, prior+sound
figure
subplot(1,3,1)
boxplot(temp1_WO_conf(:,temp_pick),'symbol','')
rec1 =temp1_WO_conf(:,temp_pick);
group1= ones(size(rec1,1),1);

subplot(1,3,2)
boxplot(temp2_WO_conf(:,temp_pick),'symbol','')
rec2 =temp3_WO_conf(:,temp_pick);
group2= ones(size(rec2,1),1)*2;

subplot(1,3,3)
boxplot(temp3_WO_conf(:,temp_pick),'symbol','')
rec3 =temp3_WO_conf(:,temp_pick);
group3= ones(size(rec3,1),1)*3;

rec_data=[rec1;rec2;rec3];
group=[group1;group2;group3];
sdata.fig = group;
sdata.choice = rec_data(:,1);
sdata.priorval =rec_data(:,2);
sdata.sound =rec_data(:,3);
sdata.sound_prior =rec_data(:,4);

name = extractBefore(folders,'_');
T = struct2table(sdata);
writetable(T, ['source figj ',name,'.csv']);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [temp1, p, min1, min1_value, temp1_WO_conf, p_WO_conf, min1_WO, min1_WO_value] = analyze_regression_each_time(temp1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp = temp1(:,1);
temp = repmat(temp,1,17);
temp1 = temp1(:,[2:end]);
temp1 = temp1 - temp;
temp1(:,[15,16]) = [];

temp_compare(1,:) = [1,2];
temp_compare(2,:) = [1,3];
temp_compare(3,:) = [1,4];
temp_compare(4,:) = [2,3];
temp_compare(5,:) = [2,4];
temp_compare(6,:) = [3,4];
%compare between RMSE
for i = 1:6
    p(i) = signrank(temp1(:,temp_compare(i,1)), temp1(:,temp_compare(i,2)));
end

%Find the minimum median
temp = mean(temp1);
min1 = find(temp == min(temp)) + 1;
min1_value = temp(min1-1);

%remove the conf_Q parts and just focus on the choice prior sound
temp = [2,3,5, 6,8,10, 13] - 1;
temp1_WO_conf = temp1(:,temp);

temp = mean(temp1_WO_conf);
min1_WO = find(temp == min(temp));
min1_WO_value = temp(min1_WO);

%temp_pick = [1,2,3,6];
clear temp_compare
temp_compare(1,:) = [1,2];
temp_compare(2,:) = [1,3];
temp_compare(3,:) = [1,6];
temp_compare(4,:) = [2,3];
temp_compare(5,:) = [2,6];
temp_compare(6,:) = [3,6];

%compare between RMSE
for i = 1:6
    p_WO_conf(i) = signrank(temp1_WO_conf(:,temp_compare(i,1)), temp1(:,temp_compare(i,2)));
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [trace,length_neuron] = get_non_nan_trace(evi_trace_all,use_frame)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get the nan frame
for k = 1:6
    moto_data = evi_trace_all(use_frame,k).matrix;
    nan_check(:,k) = isnan(moto_data); %detect_non_nan
end
nan_check = max(nan_check,[],2);
use_neuron = find(nan_check == 0);
length_neuron = length(use_neuron);

for k = 1:6
    moto_data = evi_trace_all(use_frame,k).matrix;
    trace(k).matrix = moto_data(use_neuron,:);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function use_neuron = plot_prior_block_trace2(prior_frame,evi_prefer_all,evi_nonprefer_all,temp_x)
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [median_trace,std_trace,se_trace] = matrix2medians(moto_data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Extract nan trial
moto_data = moto_data(isnan(moto_data) == 0);

median_trace = median(moto_data);
std_trace  = median(abs(moto_data-median_trace));
se_trace = 1.4826 * std_trace ./ (sqrt(length(moto_data)));

return
