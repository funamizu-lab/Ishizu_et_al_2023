
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
function Figure7a3_Population_decoder_SLR_one_session_prior

close all
pathname='G:\Ishizu_data\IntermediateFiles\mpfc\a04\2021-03-04_a04_mPFC_Right_OK\recording1_task';
cd(pathname)

SLR_prior = 'SLR_20220520_0912_glmnet_priorV_para2_depth.mat';

use_frame = 1;

%Get number of trials
load(SLR_prior);
lasso_sound_l = lasso_sound_l(use_frame).matrix;
%b_CV CVerr_CV y_test y_binary correct_error likelihood lambda x_opt c_opt
lasso_sound_s = lasso_sound_s(use_frame).matrix;

temp = dir('task_frame*'); %Task frame
if length(temp) ~= 1
    hoge
end
load(temp.name);

temp = dir('Bpod*'); %Bpod
if length(temp) ~= 1
    hoge
end
load(temp.name);
Bpod_file = temp.name;

[~,~,~,~,~,~,~,use_trial_all] = Dual_get_basic_task_structure_20210204(Bpod_file);

%Get task parameter
stim_length = unique(StimDuration);
long  = find(StimDuration == stim_length(2));
use_trial = use_trial_all;

%Check block change
TrialBlock_use = TrialBlock(use_trial);
temp = TrialBlock_use(2:length(TrialBlock_use)) - TrialBlock_use(1:length(TrialBlock_use)-1);
Blockchange_use = find(temp ~= 0) + 0.5; %start of new block

long23 = intersect(long,use_trial);
use_long_short = zeros(length(use_trial),1);
for i = 1:length(use_trial)
    temp = find(long23 == use_trial(i));
    if ~isempty(temp)
        use_long_short(i) = 1; %long
    end
end

prior_SLR = nan(length(use_trial),1);
prior_SLR(use_long_short == 1) = lasso_sound_l.y_test;
prior_SLR(use_long_short == 0) = lasso_sound_s.y_test;

%% Fig6a bottom
cd('G:\upload_code\Figure7\Fig7a');
figure;hold on
sdata = struct();% source data 
for i = 1:length(Blockchange_use)
    line([Blockchange_use(i),Blockchange_use(i)], [0 1],'color',[0 0 0],'LineWidth',0.5);
end
plot([0 length(prior_SLR)],[0.5 0.5],'r')
plot(prior_SLR,'k.')
sdata.trials = (1:length(prior_SLR))';
sdata.p_rightprior=prior_SLR;
T = struct2table(sdata);
writetable(T, 'source fig 7a bottm.csv');


