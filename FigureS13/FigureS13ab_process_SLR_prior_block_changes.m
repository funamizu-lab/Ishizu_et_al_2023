
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
function FigureS13ab_process_SLR_prior_block_changes

close all
auc_b = test_process_20221216_SLR_prior_block_changes('auc_ishizu');
fof_b = test_process_20221216_SLR_prior_block_changes('fof_ishizu');
mpfc_b= test_process_20221216_SLR_prior_block_changes('mpfc_ishizu');


%%% FigureS13b left %%%
auc_b(auc_b(:,2) < 0.004,:) = [];
fof_b(fof_b(:,2) < 0.004,:) = [];
mpfc_b(mpfc_b(:,2) < 0.004,:) = [];

figure; hold on;
plot(auc_b(:,2),auc_b(:,1),'b.')
plot(mpfc_b(:,2),mpfc_b(:,1),'g.')
plot(fof_b(:,2),fof_b(:,1),'r.')

%%% FigureS13b right %%%
auc_sabun_b = auc_b(:,1) ./ auc_b(:,2);
mpfc_sabun_b= mpfc_b(:,1) ./ mpfc_b(:,2);
fof_sabun_b = fof_b(:,1) ./ fof_b(:,2);
% p_sabun_b(1) = ranksum(auc_sabun_b, mpfc_sabun_b);
% p_sabun_b(2) = ranksum(fof_sabun_b, mpfc_sabun_b);
% p_sabun_b(3) = ranksum(auc_sabun_b, fof_sabun_b);

data = [auc_sabun_b;fof_sabun_b;mpfc_sabun_b];
group= [repmat({'auc'},length(auc_sabun_b),1);...
    repmat({'m2'},length(fof_sabun_b),1);...
    repmat({'mpfc'},length(mpfc_sabun_b),1)];

figure
hold on
boxplot(data,group)
temp_x = ones(length(auc_sabun_b),1) + 0.2 * (0.5 - rand(length(auc_sabun_b),1));
plot(temp_x,auc_sabun_b,'b.')

temp_x = ones(length(fof_sabun_b),1)*2 + 0.2 * (0.5 - rand(length(fof_sabun_b),1));
plot(temp_x,fof_sabun_b,'r.')

temp_x = ones(length(mpfc_sabun_b),1)*3 + 0.2 * (0.5 - rand(length(mpfc_sabun_b),1));
plot(temp_x,mpfc_sabun_b,'y.')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function b = test_process_20221216_SLR_prior_block_changes(folders)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

analysis_dir = eval(folders);

N_session = 0;
for i = 1:length(analysis_dir)
    [i,length(analysis_dir)]
    [block_change_prior, RL_prior] = ...
        Population_decoder_20220520_1216_SLR_one_session_prior(analysis_dir{i});
    
    if length(RL_prior) ~= 0
        N_session = N_session + 1; 
        ave_block_change_prior(N_session,:) = mean(block_change_prior);
        ave_RL_prior(N_session,:) = mean(RL_prior);
        
        
        %Get the regression parameter for each session here
        b(N_session,1) = get_slope(block_change_prior(:,21:end));
        b(N_session,2) = get_slope(RL_prior(:,21:end));
    end        
end
norm_parameter(ave_block_change_prior, ave_RL_prior, 0);
% norm_parameter(ave_block_change_prior, ave_RL_prior, 4);
% norm_parameter(ave_block_change_prior, ave_RL_prior, 9);

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [prior_change,RL_change,session_not_use] = norm_parameter(ave_block_change_prior, ave_RL_prior, last_trial)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[size_y,size_x] = size(ave_block_change_prior);
for i = 1:size_x-last_trial
    use_trial = i:i+last_trial;
    
    temp = ave_block_change_prior(:,use_trial);
    norm_prior(:,i) = mean(temp,2);

    temp = ave_RL_prior(:,use_trial);
    norm_RL(:,i) = mean(temp,2);
end

%When did the value below 0.5
for i = 1:size_y
    temp = norm_prior(i,21-last_trial:end);
    temp = find(temp > 0.5,1);
    if length(temp) ~= 0
        prior_change(i,1) = temp;
    else
        prior_change(i,1) = 40;
    end
    temp = norm_RL(i,21-last_trial:end);
    temp = find(temp > 0.5,1);
    if length(temp) ~= 0
        RL_change(i,1) = temp;
    else
        RL_change(i,1) = 40;
    end
end

session_not_use = [];

%%% FigS13a %%%
figure; hold on;
plot_mean_se_moto(norm_prior,[0 0 1],2)
plot_mean_se_moto(norm_RL,[1 0 0],2)
set(gca,'ylim',[0.2 0.8])

return

%Get the slope value for prior and RL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function b = get_slope(temp_y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[size_y,size_x] = size(temp_y);
temp_time = 1:size_x;
temp_time = repmat(temp_time,size_y,1);
temp_y = reshape(temp_y,size_y*size_x,1);
temp_time = reshape(temp_time,size_y*size_x,1);
temp_const = ones(length(temp_time),1);
b = regress(temp_y,[temp_time,temp_const]);
b = b(1);

return