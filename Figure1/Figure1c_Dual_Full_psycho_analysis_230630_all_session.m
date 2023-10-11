%{
----------------------------------------------------------------------------
Analyzing behavioral data
At least for the correct rate
----------------------------------------------------------------------------
%}

function Figure1c_Dual_Full_psycho_analysis_230630_all_session

analysis_folder = {
    'G:\Ishizu_data\Tokyo_ephys_ishizu\only_all_behaviors\a04_behave\'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\only_all_behaviors\a08_behave\'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\only_all_behaviors\i20_behave\'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\only_all_behaviors\i24_behave\'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\only_all_behaviors\i34_behave\'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\only_all_behaviors\i35_behave\'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\only_all_behaviors\i43_behave\'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\only_all_behaviors\i46_behave\'
    };

clear ave_likeli BIC log_likeli para_max
all_min_L.chosen = [];
all_min_R.chosen = [];
all_max_L.chosen = [];
all_max_R.chosen = [];

all_min_L.Evi = [];
all_min_R.Evi = [];
all_max_L.Evi = [];
all_max_R.Evi = [];

min_L_all = [];
min_R_all = [];
max_L_all = [];
max_R_all = [];
r_min_L_all = [];
r_min_R_all = [];
r_max_L_all = [];
r_max_R_all = [];

evi_x = 0:0.01:1;
sound_evi = [0 0.25 0.45 0.55 0.75 1];

for i = 1:length(analysis_folder)
    temp_folder = analysis_folder{i};
    cd(temp_folder);
    filename1 = dir('Bpod*.mat');
    %     length_session(i) = length(filename1);
    
    clear temp_file
    for j = 1:length(filename1)
        temp_file{j} = filename1(j).name;
    end
    
    [min_L,min_R,max_L,max_R,opt_min_L,opt_min_R,opt_max_L,opt_max_R, ...
        r_min_L, r_min_R, r_max_L, r_max_R]...
        = get_full_psychometric_all_session(temp_file, temp_folder);
    
    %Make psychometric function in each session
    [mouse_min_L(i,:),~,~,conf_min_L] = make_full_gauss_psychometric(min_L.chosen, min_L.Evi);
    [mouse_min_R(i,:),~,~,conf_min_R] = make_full_gauss_psychometric(min_R.chosen, min_R.Evi);
    [mouse_max_L(i,:),~,~,conf_max_L] = make_full_gauss_psychometric(max_L.chosen, max_L.Evi);
    [mouse_max_R(i,:),~,~,conf_max_R] = make_full_gauss_psychometric(max_R.chosen, max_R.Evi);
    
    figure
    subplot(2,2,1);  hold on
    plot(evi_x, opt_min_L', 'b')
    plot(evi_x, opt_min_R', 'r')
    set(gca,'xlim',[-0.1 1.1], 'ylim',[0 1])
    subplot(2,2,2);  hold on
    plot(evi_x, opt_max_L', 'b')
    plot(evi_x, opt_max_R', 'r')
    set(gca,'xlim',[-0.1 1.1], 'ylim',[0 1])
%     subplot(2,2,3)
%     plot_psychometric_session(evi_x,mouse_min_L(i,:),mouse_min_R(i,:),sound_evi,conf_min_L,conf_min_R,r_min_L,r_min_R);
%     subplot(2,2,4)
%     plot_psychometric_session(evi_x,mouse_max_L(i,:),mouse_max_R(i,:),sound_evi,conf_max_L,conf_max_R,r_max_L,r_max_R);
    
    %Integrate all trials
    min_L_all = [min_L_all; mouse_min_L];
    min_R_all = [min_R_all; mouse_min_R];
    max_L_all = [max_L_all; mouse_max_L];
    max_R_all = [max_R_all; mouse_max_R];
    
    all_min_L = integrate_min_trials(all_min_L, min_L);
    all_min_R = integrate_min_trials(all_min_R, min_R);
    all_max_L = integrate_min_trials(all_max_L, max_L);
    all_max_R = integrate_min_trials(all_max_R, max_R);
    
    r_min_L_all = [r_min_L_all; r_min_L];
    r_min_R_all = [r_min_R_all; r_min_R];
    r_max_L_all = [r_max_L_all; r_max_L];
    r_max_R_all = [r_max_R_all; r_max_R];
end

min_L = make_full_gauss_psychometric(all_min_L.chosen, all_min_L.Evi);
min_R = make_full_gauss_psychometric(all_min_R.chosen, all_min_R.Evi);
max_L = make_full_gauss_psychometric(all_max_L.chosen, all_max_L.Evi);
max_R = make_full_gauss_psychometric(all_max_R.chosen, all_max_R.Evi);

figure
subplot(2,2,1); hold on
% plot(evi_x, min_L_all', 'b')
% plot(evi_x, min_R_all', 'r')
plot(evi_x, mouse_min_L', 'b')
plot(evi_x, mouse_min_R', 'r')
set(gca,'xlim',[-0.1 1.1], 'ylim',[0 1])

subplot(2,2,2); hold on
% plot(evi_x, max_L_all', 'b')
% plot(evi_x, max_R_all', 'r')
plot(evi_x, mouse_max_L', 'b')
plot(evi_x, mouse_max_R', 'r')
set(gca,'xlim',[-0.1 1.1], 'ylim',[0 1])

subplot(2,2,3)
plot_psychometric_session(evi_x,min_L,min_R,sound_evi,r_min_L_all,r_min_R_all)

subplot(2,2,4)
plot_psychometric_session(evi_x,max_L,max_R,sound_evi,r_max_L_all,r_max_R_all)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_psychometric_session(evi_x,min_L,min_R,sound_evi, r_L, r_R)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mean_r_L = mean(r_L);
mean_r_R = mean(r_R);
std_r_L = std(r_L);
std_r_R = std(r_R);

hold on
plot(evi_x, min_L, 'b')
plot(evi_x, min_R, 'r')
for j = 1:6
    temp_x = [sound_evi(j), sound_evi(j)];
    temp_y_L = [mean_r_L(j)-std_r_L(j),mean_r_L(j)+std_r_L(j)];
    temp_y_R = [mean_r_R(j)-std_r_R(j),mean_r_R(j)+std_r_R(j)];
    
    plot(temp_x,temp_y_L,'b')
    hold on
    plot(temp_x,temp_y_R,'r')
    hold on
end

set(gca,'xlim',[-0.1 1.1], 'ylim',[0 1])

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [min_L,min_R,max_L,max_R,opt_min_L,opt_min_R,opt_max_L,opt_max_R, ...
    r_min_L, r_min_R, r_max_L, r_max_R]...
    = get_full_psychometric_all_session(filename1, temp_path)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

min_L.chosen = [];
min_R.chosen = [];
max_L.chosen = [];
max_R.chosen = [];

min_L.Evi = [];
min_R.Evi = [];
max_L.Evi = [];
max_R.Evi = [];

for filecount = 1 : length(filename1)
    [filecount, length(filename1)]
    temp_filename = filename1{filecount};
    fpath = fullfile(temp_path, temp_filename);
    
    [temp_min_L,temp_min_R,temp_max_L,temp_max_R] = integrate_behave_data_230630(fpath);
    
    min_L = integrate_min_trials(min_L, temp_min_L);
    min_R = integrate_min_trials(min_R, temp_min_R);
    max_L = integrate_min_trials(max_L, temp_max_L);
    max_R = integrate_min_trials(max_R, temp_max_R);
    
    %Make psychometric function in each session
    [opt_min_L(filecount,:),~,r_min_L(filecount,:)] = make_full_gauss_psychometric(temp_min_L.chosen, temp_min_L.Evi);
    [opt_min_R(filecount,:),~,r_min_R(filecount,:)] = make_full_gauss_psychometric(temp_min_R.chosen, temp_min_R.Evi);
    [opt_max_L(filecount,:),~,r_max_L(filecount,:)] = make_full_gauss_psychometric(temp_max_L.chosen, temp_max_L.Evi);
    [opt_max_R(filecount,:),~,r_max_R(filecount,:)] = make_full_gauss_psychometric(temp_max_R.chosen, temp_max_R.Evi);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [opt,X,p_R,conf_R] = make_full_gauss_psychometric(Chosen_side, trial_evidence)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opt = optimset('Display','off');

evi_x = 0:0.01:1;
para = [0.5 0.01 0 0
    0.5 0.1  0 0
    0.5 0.2  0 0
    0.5 0.5  0 0
    0.5 10    0 0
    0.5 0.01 0.1 0.1
    0.5 0.1  0.1 0.1
    0.5 0.2  0.1 0.1
    0.5 0.5  0.1 0.1
    0.5 10    0.1 0.1];

[right_trial, number_trial] = get_right_choice_trials(Chosen_side,trial_evidence);
prob_right = right_trial ./ number_trial;
lapse = [prob_right(1), 1-prob_right(6)]; %limit for lapse

for i = 1:6
    [p_R(i,1),conf_R(i,:)] = binofit(right_trial(i),number_trial(i));
end

parfor i = 1:10
    [X(i,:),FCAL(i)] = fminsearch(@Opt_psychometric_Gauss,para(i,:),opt,Chosen_side, trial_evidence, evi_x, lapse);
end
min_L = find(FCAL == min(FCAL),1);
X = X(min_L,:);

[~,opt,X] = Opt_psychometric_Gauss_max(X, Chosen_side, trial_evidence, evi_x, lapse);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [right_trial, number_trial] = get_right_choice_trials(use_choice,use_tone)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

use_tone = round(use_tone,3);

tone_evidence = unique(use_tone);

for i = 1:length(tone_evidence)
    temp_trial = find(use_tone == tone_evidence(i));
    temp_choice = use_choice(temp_trial);
    number_trial(i) = length(temp_trial);
    right_trial(i) = sum(temp_choice);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function log_likeli = Opt_psychometric_Gauss(para, chosen_side, tone_evi, evi_x, lapse_limit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[log_likeli,~] = Opt_psychometric_Gauss_max(para, chosen_side, tone_evi, evi_x, lapse_limit);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [log_likeli,neurometric,para] = Opt_psychometric_Gauss_max(para, chosen_side, tone_evi, evi_x, lapse_limit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%para(1): threthold: should be around 0.5
%para(2): standard deviation of Gaussian: from 0.01 to 1
%para(3): lambda1
%para(4): lambda2

%4 parameters
%Yn = lambda1 + (1-lambda1-lambda2) .* norminv(x,std)

if para(1) < 0 %para(2) shoud be positive
    para(1) = 0;
elseif para(1) > 1
    para(1) = 1;
end
if para(2) < 0 %para(2) shoud be positive
    para(2) = eps;
end
if para(3) < 0
    para(3) = 0;
elseif para(3) > lapse_limit(1)
    para(3) = lapse_limit(1);
end
if para(4) < 0
    para(4) = 0;
elseif para(4) > lapse_limit(2)
    para(4) = lapse_limit(2);
end

N_trial = size(tone_evi,1);
likelihood = zeros(1,N_trial);

%Trial by trial, get the likelihood
clear temp_p temp_exp
temp_right = normcdf(1, tone_evi, para(2)) - normcdf(para(1), tone_evi, para(2)); %Right choice probablity
temp_left  = normcdf(para(1), tone_evi, para(2)) - normcdf(0, tone_evi, para(2)); %Left choice probablity
temp_gauss = temp_right ./ (temp_right + temp_left); %normalized with truncated part
temp_p = para(3) + (1-para(3)-para(4)) .* temp_gauss;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp1 = find(chosen_side == 0);
temp2 = find(chosen_side == 1);
likelihood(temp1) = 1-temp_p(temp1); %left choice
likelihood(temp2) = temp_p(temp2);

%likelihood keisan
log_likeli = sum(log(likelihood));
log_likeli = -log_likeli;

%get the tuning function with evi_x
temp_right = normcdf(1, evi_x, para(2)) - normcdf(para(1), evi_x, para(2)); %Right choice probablity
temp_left  = normcdf(para(1), evi_x, para(2)) - normcdf(0, evi_x, para(2)); %Left choice probablity
temp_gauss = temp_right ./ (temp_right + temp_left); %normalized with truncated part
neurometric = para(3) + (1-para(3)-para(4)) .* temp_gauss;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function min_L = integrate_min_trials(min_L, temp_min_L)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

min_L.chosen = [min_L.chosen; temp_min_L.chosen];
min_L.Evi = [min_L.Evi; temp_min_L.Evi];

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [use_min_L, use_min_R, use_max_L, use_max_R] = integrate_behave_data_230630(filename1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(filename1)

[minD_trial,maxD_trial,~,~,~,use_trial2,use_trial3,~,~,~,~,~,~,~,binary_tone,~,~,~,~] ...
    = Dual_get_basic_task_structure_20210204(filename1);

if BlockProb(2) ~= BlockProb(3) %Block change task
    if BlockProb(2) > BlockProb(3) % Right -> Left
        block_R = use_trial2;
        block_L = use_trial3;
        
    else % Left -> Right
        block_L = use_trial2;
        block_R = use_trial3;
        
    end
else %Reward change task
    if BlockReward(2,1) < BlockReward(2,2) % Right -> Left
        block_R = use_trial2;
        block_L = use_trial3;
    else % Left -> Right
        block_L = use_trial2;
        block_R = use_trial3;
    end
end

%min stimulus
min_L = intersect(block_L, minD_trial);
min_R = intersect(block_R, minD_trial);
max_L = intersect(block_L, maxD_trial);
max_R = intersect(block_R, maxD_trial);

use_min_L.chosen = Chosen_side(min_L);
use_min_R.chosen = Chosen_side(min_R);
use_max_L.chosen = Chosen_side(max_L);
use_max_R.chosen = Chosen_side(max_R);

use_min_L.Evi = binary_tone(min_L);
use_min_R.Evi = binary_tone(min_R);
use_max_L.Evi = binary_tone(max_L);
use_max_R.Evi = binary_tone(max_R);
return


