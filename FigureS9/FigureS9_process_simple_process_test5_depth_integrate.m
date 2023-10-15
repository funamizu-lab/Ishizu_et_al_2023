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
function FigureS9_process_simple_process_test5_depth_integrate

close all
use_frame = 15:25; % 15:25 (-0.1-1.0 s) / 16:25 (0-1.0 s) in long-sound trials 

[mpfc_prior,mpfc_choice]= value_process_20220516_simple_process_test5_depth_control('mpfc_ishizu');
[fof_prior, fof_choice] = value_process_20220516_simple_process_test5_depth_control('fof_ishizu');

figure
subplot(2,2,1)
plot_median_se_moto_x_axis_matrix(mpfc_prior, 1:length(use_frame), [0 0 0],2)
subplot(2,2,2)
plot_median_se_moto_x_axis_matrix(mpfc_choice,1:length(use_frame), [0 0 0],2)
subplot(2,2,3)
plot_median_se_moto_x_axis_matrix(fof_prior, 1:length(use_frame), [0 0 0],2)
subplot(2,2,4)
plot_median_se_moto_x_axis_matrix(fof_choice,1:length(use_frame), [0 0 0],2)

length_data = length(mpfc_prior);

disp('prior all')
BIC_computation_mpfc_fof(mpfc_prior, fof_prior, [1,length_data]);
disp('prior without first')
BIC_computation_mpfc_fof(mpfc_prior, fof_prior, [2,length_data]);

disp('choice all')
BIC_computation_mpfc_fof(mpfc_choice, fof_choice, [1,length_data]);
disp('choice without first')
BIC_computation_mpfc_fof(mpfc_choice, fof_choice, [2,length_data]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function BIC_computation_mpfc_fof(mpfc_prior, fof_prior, use_data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[y_mpfc_prior, x_mpfc_prior] = reshape_neuron_data(mpfc_prior,use_data);
[y_fof_prior, x_fof_prior] = reshape_neuron_data(fof_prior,use_data);

moto_b(:,1) = glmfit(x_mpfc_prior,y_mpfc_prior,'normal','link','identity');
moto_b(:,2) = glmfit(x_fof_prior,y_fof_prior,'normal','link','identity');
moto_b(:,3) = glmfit([x_mpfc_prior; x_fof_prior],[y_mpfc_prior; y_fof_prior],'normal','link','identity');

opt = optimset('Display','off');

length_data = length(y_mpfc_prior)+length(y_fof_prior);

%Use two parameters for the computation.
for i = 1:3
    [b2(:,i),FCAL_all2(i)] = fminsearch(@RMSE_computation2,moto_b(:,i),opt,...
        x_mpfc_prior, y_mpfc_prior, x_fof_prior, y_fof_prior);
    
    [b3b(:,i),FCAL_all3b(i)] = fminsearch(@RMSE_computation3_b,[moto_b(1,i),moto_b(1,i),moto_b(2,i)],opt,...
        x_mpfc_prior, y_mpfc_prior, x_fof_prior, y_fof_prior);

    [b3a(:,i),FCAL_all3a(i)] = fminsearch(@RMSE_computation3_a,[moto_b(1,i),moto_b(2,i),moto_b(2,i)],opt,...
        x_mpfc_prior, y_mpfc_prior, x_fof_prior, y_fof_prior);

    [b4(:,i),FCAL_all4(i)] = fminsearch(@RMSE_computation4,[moto_b(1,i),moto_b(1,i),moto_b(2,i),moto_b(2,i)],opt,...
        x_mpfc_prior, y_mpfc_prior, x_fof_prior, y_fof_prior);
end
use_FCAL = find(FCAL_all2 == min(FCAL_all2),1);
FCAL_all2 = FCAL_all2(use_FCAL);
b2 = b2(:,use_FCAL)';

use_FCAL = find(FCAL_all3b == min(FCAL_all3b),1);
FCAL_all3b = FCAL_all3b(use_FCAL);
b3b = b3b(:,use_FCAL)';

use_FCAL = find(FCAL_all3a == min(FCAL_all3a),1);
FCAL_all3a = FCAL_all3a(use_FCAL);
b3a = b3a(:,use_FCAL)';

use_FCAL = find(FCAL_all4 == min(FCAL_all4),1);
FCAL_all4 = FCAL_all4(use_FCAL);
b4 = b4(:,use_FCAL)';

BIC(1) = length_data * log(FCAL_all2/length_data) + 2 * log(length_data);
BIC(2) = length_data * log(FCAL_all3b/length_data)+ 3 * log(length_data);
BIC(3) = length_data * log(FCAL_all3a/length_data)+ 3 * log(length_data);
BIC(4) = length_data * log(FCAL_all4/length_data) + 4 * log(length_data);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sabun_integ = RMSE_computation2(para, x_mpfc_prior, y_mpfc_prior, x_fof_prior, y_fof_prior)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
est_mpfc = para(2) * x_mpfc_prior + para(1);
est_fof = para(2) * x_fof_prior + para(1);

sabun_mpfc = y_mpfc_prior - est_mpfc;
sabun_fof = y_fof_prior - est_fof;

sabun_mpfc = sabun_mpfc .* sabun_mpfc;
sabun_fof = sabun_fof .* sabun_fof;

sabun_integ = (sum(sabun_mpfc) + sum(sabun_fof));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sabun_integ = RMSE_computation3_b(para, x_mpfc_prior, y_mpfc_prior, x_fof_prior, y_fof_prior)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
est_mpfc = para(3) * x_mpfc_prior + para(1);
est_fof = para(3) * x_fof_prior + para(2);

sabun_mpfc = y_mpfc_prior - est_mpfc;
sabun_fof = y_fof_prior - est_fof;

sabun_mpfc = sabun_mpfc .* sabun_mpfc;
sabun_fof = sabun_fof .* sabun_fof;

sabun_integ = (sum(sabun_mpfc) + sum(sabun_fof));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sabun_integ = RMSE_computation3_a(para, x_mpfc_prior, y_mpfc_prior, x_fof_prior, y_fof_prior)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
est_mpfc = para(2) * x_mpfc_prior + para(1);
est_fof = para(3) * x_fof_prior + para(1);

sabun_mpfc = y_mpfc_prior - est_mpfc;
sabun_fof = y_fof_prior - est_fof;

sabun_mpfc = sabun_mpfc .* sabun_mpfc;
sabun_fof = sabun_fof .* sabun_fof;

sabun_integ = (sum(sabun_mpfc) + sum(sabun_fof));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sabun_integ = RMSE_computation4(para, x_mpfc_prior, y_mpfc_prior, x_fof_prior, y_fof_prior)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
est_mpfc = para(3) * x_mpfc_prior + para(1);
est_fof = para(4) * x_fof_prior + para(2);

sabun_mpfc = y_mpfc_prior - est_mpfc;
sabun_fof = y_fof_prior - est_fof;

sabun_mpfc = sabun_mpfc .* sabun_mpfc;
sabun_fof = sabun_fof .* sabun_fof;

sabun_integ = (sum(sabun_mpfc) + sum(sabun_fof));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [all_y, all_x] = reshape_neuron_data(mpfc_prior, use_data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%reshape the data
all_y = [];
all_x = [];
for i = use_data(1):use_data(2)
    temp_y = mpfc_prior(i).matrix;
    temp_x = ones(length(temp_y),1) * i;
    
    all_y = [all_y; temp_y];
    all_x = [all_x; temp_x];
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [neuron_prior, neuron_choice] = value_process_20220516_simple_process_test5_depth_control(folders)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

analysis_dir = eval(folders);

for i = 1:40
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
   
    [~,~,~,~,~,~,~,~,~,~,~,evi_trace,evi_prefer,evi_nonprefer,sig_trace_prior,sig_trace_nonprior] = ...
        Task_kaiseki_tokyo1_20220516_sound_choice_process4_depth(analysis_dir{i});
    
    for j = 1:40
        
        sig_trace_prior_all(j).matrix = [sig_trace_prior_all(j).matrix; sig_trace_prior(j).matrix];
        sig_trace_nonprior_all(j).matrix = [sig_trace_nonprior_all(j).matrix; sig_trace_nonprior(j).matrix];
        
        for k = 1:6
            evi_trace_all(j,k).matrix = [evi_trace_all(j,k).matrix; evi_trace(j,k).matrix];
            evi_prefer_all(j,k).matrix = [evi_prefer_all(j,k).matrix; evi_prefer(j,k).matrix];
            evi_nonprefer_all(j,k).matrix = [evi_nonprefer_all(j,k).matrix; evi_nonprefer(j,k).matrix];
        end
    end    
end

%Sound on: 16
%Sound end: 25
use_frame = 15:25;
for i = 1:length(use_frame)
    [prefer_activ(i).matrix,nonprefer_activ(i).matrix] = plot_prior_block_trace(use_frame(i),evi_prefer_all,evi_nonprefer_all);
    
    temp_p = prefer_activ(i).matrix;
    temp_n = nonprefer_activ(i).matrix;
    
    neuron_prior(i).matrix = mean(temp_p,2) - mean(temp_n,2);
    
    temp = (temp_p + temp_n) ./ 2;
    neuron_choice(i).matrix = mean(temp(:,4:6),2)-mean(temp(:,1:3),2);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [prefer_activ,nonprefer_activ] = plot_prior_block_trace(prior_frame,evi_prefer_all,evi_nonprefer_all)
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
end
