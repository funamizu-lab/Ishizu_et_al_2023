
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
function process_20230923_simple_process_test5_depth_integrate_CV2
%number of training data are now same in all the conditions

use_frame = 15:25;

[mpfc_prior, mpfc_choice] = value_process_20220516_simple_process_test5_depth_control('mpfc_ishizu');
[fof_prior, fof_choice] = value_process_20220516_simple_process_test5_depth_control('fof_ishizu');

close all
figure
subplot(2,2,1)
plot_median_se_moto_x_axis_matrix(mpfc_prior, [1:length(use_frame)], [0 0 0],2)
subplot(2,2,2)
plot_median_se_moto_x_axis_matrix(mpfc_choice, [1:length(use_frame)], [0 0 0],2)
subplot(2,2,3)
plot_median_se_moto_x_axis_matrix(fof_prior, [1:length(use_frame)], [0 0 0],2)
subplot(2,2,4)
plot_median_se_moto_x_axis_matrix(fof_choice, [1:length(use_frame)], [0 0 0],2)

length_data = length(mpfc_prior);

for i = 1:length_data
    p(i) = ranksum(mpfc_prior(i).matrix,fof_prior(i).matrix);
end

CV_para = 10;

disp('prior all')
CV_computation_mpfc_fof(mpfc_prior, fof_prior, [1,length_data], CV_para);
disp('prior without first')
CV_computation_mpfc_fof(mpfc_prior, fof_prior, [2,length_data], CV_para);

disp('choice all')
CV_computation_mpfc_fof(mpfc_choice, fof_choice, [1,length_data], CV_para);
disp('choice without first')
CV_computation_mpfc_fof(mpfc_choice, fof_choice, [2,length_data], CV_para);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CV_computation_mpfc_fof(mpfc_prior, fof_prior, use_data, CV_para)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[y_mpfc_prior, x_mpfc_prior] = reshape_neuron_data(mpfc_prior,use_data);
[y_fof_prior, x_fof_prior] = reshape_neuron_data(fof_prior,use_data);

% length_data = length(y_mpfc_prior)+length(y_fof_prior);
length_mpfc = length(y_mpfc_prior);
length_fof = length(y_fof_prior);

CV_repeat = 1000;
for k = 1:CV_repeat

clear vali_mpfc vali_fof train_mpfc train_fof
temp_mpfc = randperm(length_mpfc);
temp_mpfc = rem(temp_mpfc,CV_para) + 1;
temp_fof = randperm(length_fof);
temp_fof = rem(temp_fof,CV_para) + 1;

for i = 1:CV_para
    vali_mpfc(i).matrix = find(temp_mpfc == i);
    vali_fof(i).matrix = find(temp_fof == i);
    
    train_mpfc(i).matrix = setdiff([1:length_mpfc], vali_mpfc(i).matrix);
    train_fof(i).matrix = setdiff([1:length_fof], vali_fof(i).matrix);
end

for j = 1:CV_para
    v_x_mpfc = x_mpfc_prior(vali_mpfc(j).matrix);
    v_y_mpfc = y_mpfc_prior(vali_mpfc(j).matrix);
    v_x_fof = x_fof_prior(vali_fof(j).matrix);
    v_y_fof = y_fof_prior(vali_fof(j).matrix);

    t_x_mpfc = x_mpfc_prior(train_mpfc(j).matrix);
    t_y_mpfc = y_mpfc_prior(train_mpfc(j).matrix);
    t_x_fof = x_fof_prior(train_fof(j).matrix);
    t_y_fof = y_fof_prior(train_fof(j).matrix);

    %Use two parameters for the computation.
    clear b2 b3b b3a b4 
    clear FCAL_all2 FCAL_all3b FCAL_all3a FCAL_all4
    
    min_train = min(length(t_x_mpfc), length(t_x_fof));
    integ_t_x = [t_x_mpfc; t_x_fof];
    integ_t_y = [t_y_mpfc; t_y_fof];

    use_integ = randperm(length(integ_t_x));
    use_integ = use_integ([1:min_train]);

    use_mpfc = randperm(length(t_x_mpfc));
    use_mpfc = use_mpfc([1:min_train]);

    use_fof = randperm(length(t_x_fof));
    use_fof = use_fof([1:min_train]);

    %2 parameters
    b2 = robustfit(integ_t_x(use_integ),integ_t_y(use_integ));
    b4_mpfc = robustfit(t_x_mpfc(use_mpfc),t_y_mpfc(use_mpfc));
    b4_fof = robustfit(t_x_fof(use_fof),t_y_fof(use_fof));
    b4 = [b4_mpfc(1); b4_fof(1); b4_mpfc(2); b4_fof(2)];
    %Validata the data
    sabun_integ(j,1) = RMSE_computation2(b2, v_x_mpfc, v_y_mpfc, v_x_fof, v_y_fof);
    sabun_integ(j,2) = RMSE_computation4(b4, v_x_mpfc, v_y_mpfc, v_x_fof, v_y_fof);
end

sabun_all(k,:) = sum(sabun_integ);
end

mean(sabun_all)

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

return

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

return

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

return

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

return

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
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [neuron_prior, neuron_choice] = value_process_20220516_simple_process_test5_depth_control(folders)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

analysis_dir = eval(folders);

for i = 1:40
    for j = 1:6
        evi_prefer_all(i,j).matrix = [];
        evi_nonprefer_all(i,j).matrix = [];
    end
end
for i = 1:length(analysis_dir)
    [i,length(analysis_dir)]
   
    [~,~,~,~,~,~,~,~,~,~,~,~, evi_prefer, evi_nonprefer] = ...
        Task_kaiseki_tokyo1_20220516_sound_choice_process4_depth(analysis_dir{i});
    
    for j = 1:40      
        for k = 1:6
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

%for anavan
for k = 1:6
    moto_data = evi_prefer_all(prior_frame,k).matrix;
    moto_data = moto_data(use_neuron);
    prefer_activ(:,k) = moto_data;

    moto_data = evi_nonprefer_all(prior_frame,k).matrix;
    moto_data = moto_data(use_neuron);
    nonprefer_activ(:,k) = moto_data;
end
return