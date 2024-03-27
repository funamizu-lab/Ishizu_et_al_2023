%{
----------------------------------------------------------------------------
%Stimulus bias changes by block
%Just model analysis
%Based on Rao 2010 Front Comp Neurosci
%Value updating + Prior updating
----------------------------------------------------------------------------
%}

function [prior,posterior,ave_likeli,likelihood,Long_posterior,Short_posterior,...
    prior_value,long_value,short_value,para] = ...
    Dual_RL_model_block1_20220314_para_determined(filename1,para,ori_trial)

switch nargin
    case 0
        hoge
    case 1
        hoge
    case 2
        hoge
    case 3
%         [filename1, pathname1]=uigetfile('*.mat','Block_mat');
%         filename1 = [pathname1, filename1];
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
if ori_trial ~= N_trial
    hoge
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

[prior,posterior,ave_likeli,likelihood,Long_posterior,Short_posterior,...
    prior_value,long_value,short_value,para] = ...
    Thre_update_210310_max_value(para, reward_trial, Chosen_trial, binary_tone, Tone_type, bin_x, stim_LR, stim_prob_LR, init_Q);

% figure
% plot(Q(:,1),'b')
% hold on
% plot(Q(:,2),'r')
% set(gca,'xlim',[1 N_trial])
% set(gca,'xtick',[0:40:N_trial])
% set(gca,'ylim',[0 4])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [prior,posterior,ave_likeli,likelihood,Long_posterior,Short_posterior,...
    prior_value,long_value,short_value,para] = ...
    Thre_update_210310_max_value(para, reward, action, Stim_freq, Tone_type, bin_x, stim_LR, stim_prob_LR, init_Q)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

softmax_check = 1;

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
prior_value = zeros(N_trial,2);
long_value = zeros(N_trial,2);
short_value = zeros(N_trial,2);

prior = zeros(N_trial,2);
posterior = zeros(N_trial,2);
Long_posterior = zeros(N_trial,2);
Short_posterior = zeros(N_trial,2);
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
    
    %Long tone
    Long_posterior(i,:) = get_left_right_prob(1,decision_x,Stim_freq(i),para);
    %Short tone
    Short_posterior(i,:) = get_left_right_prob(2,decision_x,Stim_freq(i),para);
    
    %Value for prior, long, short
    prior_value(i,:) = Q(i,:);
    long_value(i,:)  = get_left_right_value(1,Q(i,:),Stim_freq(i),para);
    short_value(i,:) = get_left_right_value(2,Q(i,:),Stim_freq(i),para);
    
    Choice_prob(i,:) = [left_choice,right_choice];
    %Get probability from action and reward;
    likelihood(i) = Choice_prob(i,action(i));
    
    prior(i,:) = choice_prior;
    posterior(i,:) = Choice_prob(i,:);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Q_stim = get_left_right_value(Tone_type,Q,Stim_freq,para)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    temp_std = para(Tone_type + 2); %3 or 4
    
    temp_left  = normcdf(0.5,Stim_freq,temp_std) - normcdf(0,Stim_freq,temp_std);
    temp_right = normcdf(1,Stim_freq,temp_std) - normcdf(0.5,Stim_freq,temp_std);
    temp_left = temp_left / (temp_left + temp_right);
    temp_right = 1 - temp_left;
    if temp_left > 1
        temp_left = 1;
        temp_right = 0;
    elseif temp_right > 1
        temp_left = 0;
        temp_right = 1;
    end
    
    Q_stim(1) = 2*Q(1) .* temp_left; %update Q with stim prob
    Q_stim(2) = 2*Q(2) .* temp_right;  %update Q with stim prob
    return
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function temp_choice = get_left_right_prob(Tone_type,decision_x,Stim_freq,para)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    temp_std = para(Tone_type + 2); %3 or 4
    temp_bias = para(Tone_type + 4); %5 or 6
    
    temp_left  = normcdf(decision_x + temp_bias,Stim_freq,temp_std) - normcdf(0,Stim_freq,temp_std);
    temp_right = normcdf(1,Stim_freq,temp_std) - normcdf(decision_x + temp_bias,Stim_freq,temp_std);
    temp_left = temp_left / (temp_left + temp_right);
    temp_right = 1 - temp_left;
    if temp_left > 1
        temp_left = 1;
        temp_right = 0;
    elseif temp_right > 1
        temp_left = 0;
        temp_right = 1;
    end
    temp_choice = [temp_left temp_right];
    return
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ave_likeli_mean = Thre_update_201108_all(para, reward_trial, Chosen_trial, binary_tone, Tone_type, bin_x, stim_LR, stim_prob_LR, init_Q, use_para)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ave_likeli_mean,~] = ...
    Thre_update_210310_max_value(para, reward_trial, Chosen_trial, binary_tone, Tone_type, bin_x, stim_LR, stim_prob_LR, init_Q, use_para);
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

