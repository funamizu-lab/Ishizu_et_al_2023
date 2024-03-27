%{
----------------------------------------------------------------------------
%Stimulus bias changes by block
%Just model analysis
%Based on Rao 2010 Front Comp Neurosci
%Value updating + Prior updating
----------------------------------------------------------------------------
%}

function [ave_likeli,BIC_all,log_likeli,para_max,N_trial] = Dual_RL_model_block1_20220314_para_search(filename1)

switch nargin
    case 0
        [filename1, pathname1]=uigetfile('*.mat','Block_mat');
        filename1 = [pathname1, filename1];
        load(filename1)
    case 1
        load(filename1)
    otherwise
        hoge
end

% %Outcome
% outcome_EW     = 0; %early withdrawal
% outcome_IC     = 1; %incorrect choice
% outcome_reward = 2; %reward was dispensed (either automatically in early training, or after correct choice)
% outcome_NC     = 3; %no choice was made and time elapsed
% outcome_UN     = 4; %undefined or Free water:

[minD_trial,maxD_trial,Choice_trial,tone_evidence,trial_evidence,use_trial2,use_trial3,use_trial_all,...
low,high,correct,error,flip_tone,number_use_trial,...
binary_tone,right_trial_all,number_trial_all,right_trial,number_trial] ...
    = Dual_get_basic_task_structure_20210204(filename1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make trial type vector
Tone_type = nan(length(Correct_side),1);
Tone_type(maxD_trial) = 1;
Tone_type(minD_trial) = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Take out chose success trials:
N_trial = length(use_trial_all);
Correct_trial = Correct_side(use_trial_all);
Chosen_trial = Chosen_side(use_trial_all);
binary_tone = binary_tone(use_trial_all);
Reward_LCR = Reward_LCR(use_trial_all,[1,3]); %choose chosen trials
Tone_type = Tone_type(use_trial_all);

clear init_Q
init_Q = BlockReward(1,:); %2.4 2.4
    
%reward amount for each trial

reward_trial = zeros(N_trial,1);
for i = 1:N_trial
    %Reward trials
    if Correct_trial(i) == Chosen_trial(i),
        reward_trial(i) = Reward_LCR(i,Chosen_trial(i)+1);
    end
end

%Parameter
bin_x = [0:0.001:1];
stim_LR   = [0,0.25,0.45,0.55,0.75,1]; %3stim
stim_prob_LR   = [0.5,0.25,0.25,0.25,0.25,0.5]; %3stim
stim_prob_kitei = [0.25,0.125,0.125,0.125,0.125,0.25]; %3stim
%Sense_std = [0.01:0.05:1];
%reward = [3,1]; %left-right
%stim_prob = [0.5, 0.5];

% %Parameter
% para(1) = 0.1; % forgetting value for long
% para(2) = 0.1; % forgetting value for short
% para(3) = 0.2;  % std of long sound
% para(4) = 0.2;  % std of short sound
% para(5) = 0;  % bias for long
% para(6) = 0;  % bias for short
% para(7) = 0.3683;  % beta for long
% para(8) = 0.3683;  % beta for short
% %para(5) = 0.1;  %inverse temperature for prior

use_para(1).matrix = [1 3 5 7];
use_para(2).matrix = [1 3 4 5 7];
% use_para(3).matrix = [1 3 4 5 6 7];
% use_para(4).matrix = [1 3 4 5 7 8];
% use_para(5).matrix = [1 3 4 5 6 7 8];

para_temp = [0   0   0.2 0.2 0 0 0.5 0.5;
             0   0   0.5 0.5 0 0 1   1;
             0   0   1   1   0 0 0.5 0.5;
             0.1 0.1 0.2 0.2 0 0 1   1;
             0.1 0.1 0.5 0.5 0 0 0.5 0.5;
             0.1 0.1 1   1   0 0 1   1;
             0.3 0.3 0.2 0.2 0 0 0.5 0.5;
             0.3 0.3 0.5 0.5 0 0 1   1;
             0.6 0.6 0.2 0.2 0 0 0.5 0.5;
             0.6 0.6 0.5 0.5 0 0 1   1;
             ];

%Get maximum prediction
opt = optimset('Display','off');
clear ave_likeli BIC_all log_likeli
for j = 1:length(use_para)
    j
    temp_use_para = use_para(j).matrix;
    parfor i = 1:10
        i
        [X1(i,:),FCAL1(i,:),EXITFLAG,OUTPUT] = ...
            fminsearch(@Thre_update_201108_all,para_temp(i,:), opt, reward_trial, Chosen_trial, binary_tone, Tone_type, bin_x, stim_LR, stim_prob_LR, init_Q, temp_use_para);
    end
    min_FCAL = find(FCAL1 == min(FCAL1),1);
    temp_para_max = X1(min_FCAL,:);
    para_temp(j,:) = temp_para_max;
    
    [ave_likeli(1,j),likelihood,Q,Choice_prob,para_max(j,:)] = ...
        Thre_update_201108_max_all(temp_para_max, reward_trial, Chosen_trial, binary_tone, Tone_type, bin_x, stim_LR, stim_prob_LR, init_Q, temp_use_para);
    
    sumlog_likeli = sum(log(likelihood));
    log_likeli(1,j) = sumlog_likeli;
    BIC_all(1,j) = -2 * sum(log(likelihood)) + length(temp_use_para) * log(N_trial);
        %log_likelihood(filecount).matrix = log(likelihood);
end

ave_likeli
BIC_all
log_likeli
para_max
N_trial

figure
plot(Q(:,1),'b')
hold on
plot(Q(:,2),'r')
set(gca,'xlim',[1 N_trial])
set(gca,'xtick',[0:40:N_trial])
set(gca,'ylim',[0 4])

%hoge

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ave_likeli,likelihood,Q,Choice_prob,para] = ...
    Thre_update_201108_max_all(para, reward, action, Stim_freq, Tone_type, bin_x, stim_LR, stim_prob_LR, init_Q, use_para)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% use_para(1).matrix = [1 3 5];
% use_para(2).matrix = [1 3 4 5];
% use_para(3).matrix = [1 2 3 5];
% use_para(4).matrix = [1 2 3 4 5];
% %Parameter
% para(1) = 0.1; % forgetting value for long
% para(2) = 0.1; % forgetting value for short
% para(3) = 0.2;  % std of long sound
% para(4) = 0.2;  % std of short sound
% para(5) = 0;  % bias for long
% para(6) = 0;  % bias for short
% para(7) = 0.3683;  % beta for long
% para(8) = 0.3683;  % beta for short

%para shusei
if para(1) < 0
    para(1) = 0;
% elseif para(1) > 1
%     para(1) = 1;
end
if para(2) < 0
    para(2) = 0;
elseif para(2) > 1
    para(2) = 1;
end
if para(3) < 0
    para(3) = 0;
end
if para(4) < 0
    para(4) = 0;
end
if para(7) < 0
    para(7) = 0;
end
if para(8) < 0
    para(8) = 0;
end

temp = find(use_para == 2, 1);
if isempty(temp)
    para(2) = para(1);
end
temp = find(use_para == 4, 1);
if isempty(temp)
    para(4) = para(3);
end
temp = find(use_para == 5, 1);
if isempty(temp)
    para(5) = 0;
end
temp = find(use_para == 6, 1);
if isempty(temp)
    para(6) = para(5);
end

temp = find(use_para == 7, 1);
if isempty(temp)
    softmax_check = 0;
else
    softmax_check = 1;
end
temp = find(use_para == 8, 1);
if isempty(temp)
    para(8) = para(7);
end

%Plot the gaussian based class for left and right
%Write the figures for bayes computation
N_trial = length(reward);
%long stimulus
for i = 1:length(stim_LR)
    temp_x = normpdf(bin_x,stim_LR(i),para(3));
    temp_x = temp_x ./ sum(temp_x);
    LR_stim(i,:) = stim_prob_LR(i) .* temp_x;
    norm_pdf_long(i,:) = temp_x;
end
%Combination of left and right
L_long = sum(LR_stim([1:3],:));
R_long = sum(LR_stim([4:6],:));
%short stimulus
for i = 1:length(stim_LR)
    temp_x = normpdf(bin_x,stim_LR(i),para(4));
    temp_x = temp_x ./ sum(temp_x);
    LR_stim(i,:) = stim_prob_LR(i) .* temp_x;
    norm_pdf_short(i,:) = temp_x;
end
%Combination of left and right
L_short = sum(LR_stim([1:3],:));
R_short = sum(LR_stim([4:6],:));


clear Prior Q_left Q_right Choice_prob
%Belief_basis for left choice reward probability
Q = zeros(N_trial,2);
%Q(1,:)  = init_Q;
Q(1,:)  = init_Q .* 0.5; %Prior with 0.5 in both side
likelihood = zeros(N_trial,1);
Choice_prob = zeros(N_trial,2);
action = action + 1; %0,1 -> 1,2

for i = 1:N_trial
    clear Posterior_LR Decision_LR
    %Tone_type: 1(long) or 2(short)
    temp_alpha = para(1); %1
    temp_forget = para(2); %1
    temp_std = para(Tone_type(i) + 2); %3 or 4
    temp_bias = para(Tone_type(i) + 4); %5 or 6
    temp_soft = para(Tone_type(i) + 6); %7 or 8
    
    if softmax_check,
        choice_prior = softmax(Q(i,:),temp_soft);
    else
        choice_prior = Q(i,:)./sum(Q(i,:));
    end
    decision_x = choice_prior(1);
    
    left_choice  = normcdf(decision_x + temp_bias,Stim_freq(i),temp_std) - normcdf(0,Stim_freq(i),temp_std);
    right_choice = normcdf(1,Stim_freq(i),temp_std) - normcdf(decision_x + temp_bias,Stim_freq(i),temp_std);
    left_choice = left_choice / (left_choice + right_choice);
    right_choice = 1 - left_choice;

    if left_choice > 1
        left_choice = 1;
        right_choice = 0;
    elseif right_choice > 1
        left_choice = 0;
        right_choice = 1;
    end

    %Confidence
    left_conf  = normcdf(0.5,Stim_freq(i),temp_std) - normcdf(0,Stim_freq(i),temp_std);
    right_conf = normcdf(1,Stim_freq(i),temp_std) - normcdf(0.5,Stim_freq(i),temp_std);
    left_conf = left_conf / (left_conf + right_conf);
    right_conf = 1 - left_conf;
    
    Choice_prob(i,:) = [left_choice,right_choice];
    %Get probability from action and reward;
    likelihood(i) = Choice_prob(i,action(i));
    
    %calculate the Q_stim
    clear temp_Q_stim
%     temp_Q_stim(1) = Q(i,1) .* left_choice;
%     temp_Q_stim(2) = Q(i,2) .* right_choice;
%    temp_Q_stim(1) = 2*Q(i,1) .* left_choice; %update Q with stim prob
%    temp_Q_stim(2) = 2*Q(i,2) .* right_choice;  %update Q with stim prob
    temp_Q_stim(1) = 2*Q(i,1) .* left_conf; %update Q with stim prob
    temp_Q_stim(2) = 2*Q(i,2) .* right_conf;  %update Q with stim prob
    
%     temp = find(use_para == 2, 1);
%     if isempty(temp)
%         temp_Q_stim = Q(i,:);
%     end
    
    %Forgetting Q learning
    %Q(i+1,:) = Qlearning_ori(Q(i,:), action(i), reward(i), temp_alpha);
    Q(i+1,:) = Qlearning_confidence2(Q(i,:), temp_Q_stim, action(i), reward(i), temp_alpha, temp_forget);
end

%likelihood keisan
log_likeli = sum(log(likelihood));
log_likeli = log_likeli / N_trial;
ave_likeli = exp(log_likeli);
ave_likeli = -ave_likeli;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function post_prob = Qlearning_confidence2(Q, Q_stim, choice, reward, alpha, para_forget)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Control alpha
temp_Q = Q(:,choice) + alpha .* (reward - Q_stim(:,choice)); 

if choice == 1
    post_prob = [temp_Q, (1-para_forget) * Q(:,2)];
    %post_prob = [temp_Q, Q(:,2)];
else
    post_prob = [(1-para_forget) * Q(:,1), temp_Q];
    %post_prob = [Q(:,1), temp_Q];
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ave_likeli_mean = Thre_update_201108_all(para, reward_trial, Chosen_trial, binary_tone, Tone_type, bin_x, stim_LR, stim_prob_LR, init_Q, use_para)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ave_likeli_mean,~] = ...
    Thre_update_201108_max_all(para, reward_trial, Chosen_trial, binary_tone, Tone_type, bin_x, stim_LR, stim_prob_LR, init_Q, use_para);
ave_likeli_mean
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function post_prob = Qlearning_ori(Q, choice, reward, alpha)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Control alpha
temp_Q = Q(:,choice) + alpha .* (reward - Q(:,choice)); 

if choice == 1
    post_prob = [temp_Q, (1-alpha) * Q(:,2)];
else
    post_prob = [(1-alpha) * Q(:,1), temp_Q];
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stim_prob = softmax(temp_stim,beta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear stime_prob
stim_prob(1) = exp(beta * temp_stim(1)) / (exp(beta * temp_stim(1)) + exp(beta * temp_stim(2)));
stim_prob(2) = exp(beta * temp_stim(2)) / (exp(beta * temp_stim(1)) + exp(beta * temp_stim(2)));

return

