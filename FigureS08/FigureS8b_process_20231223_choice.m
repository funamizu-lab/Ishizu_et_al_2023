
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
function FigureS8b_process_20231223_choice(folders)

analysis_dir = eval(folders);
use_para = [2,3,6,4]; %for choice

for i = 1:40    
    long_time_RMSE(i).matrix = [];
    for j = 1:6
        evi_prefer_all(i,j).matrix = [];
        evi_nonprefer_all(i,j).matrix = [];
    end
end
for i = 1:length(analysis_dir)
    [i,length(analysis_dir)]
   
    [~,~,~,~,~,~,~,~,~,~,~,~,evi_prefer,evi_nonprefer,~,~,~,~,~,~,~,~,~,temp_RMSE] = ...
        Task_kaiseki_tokyo1_20230920_sound_choice_process4_depth(analysis_dir{i});    
    
    for j = 1:40
        long_time_RMSE(j).matrix = [long_time_RMSE(j).matrix, temp_RMSE(j).matrix];        
        for k = 1:6
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
long_frame = 25;
use_neuron1 = plot_prior_block_trace2(prior_frame,evi_prefer_all,evi_nonprefer_all);
use_neuron2 = plot_prior_block_trace2(short_frame,evi_prefer_all,evi_nonprefer_all);
use_neuron3 = plot_prior_block_trace2(long_frame,evi_prefer_all,evi_nonprefer_all);
 
%
%make boxplot for each RMSE
%choice, prior, conf_Q, binary_tone, correct_side
temp1 = long_time_RMSE(prior_frame).matrix';
temp2 = long_time_RMSE(short_frame).matrix';
temp3 = long_time_RMSE(long_frame).matrix';

%Update the number of neurons with use_neuron
temp1 = temp1(use_neuron1,:);
temp2 = temp2(use_neuron2,:);
temp3 = temp3(use_neuron3,:);

temp1_WO_conf = analyze_regression_each_time(temp1, use_para);
temp2_WO_conf = analyze_regression_each_time(temp2, use_para);
temp3_WO_conf = analyze_regression_each_time(temp3, use_para);

figure
ax(1) = subplot(1,3,1);
boxplot(ax(1), temp1_WO_conf(:,1:3), 'symbol', '')
ax(2) = subplot(1,3,2);
boxplot(ax(2), temp2_WO_conf(:,1:3), 'symbol', '')
ax(3) = subplot(1,3,3);
boxplot(ax(3), temp3_WO_conf, 'symbol', '')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function temp1_WO_conf = analyze_regression_each_time(temp1, use_para)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

temp = temp1(:,1);
temp = repmat(temp,1,17);
temp1 = temp1(:,2:end);
temp1 = temp1 - temp;
temp1(:,[15,16]) = [];

%remove the conf_Q parts and just focus on the choice prior sound
temp = use_para - 1;
temp1_WO_conf = temp1(:,temp);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function use_neuron = plot_prior_block_trace2(prior_frame,evi_prefer_all,evi_nonprefer_all)
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