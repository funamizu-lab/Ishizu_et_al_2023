%{
----------------------------------------------------------------------------
%Stimulus bias changes by block
%Just model analysis
%Based on Rao 2010 Front Comp Neurosci
%Value updating + Prior updating
----------------------------------------------------------------------------
%}

function Figure2b_Dual_RL_model_block1_220818_plot_values2(filename1)

% data = a04_2021_0227
% ~\IntermediateFiles\only_all_behaviors\a04_behave\Bpod_mat_210202_Feb27-2021_Session1.mat
% ~\IntermediateFiles\only_all_behaviors\a04_behave\ML_CumGauss_Fit\RL_20220818\RL_20220818_confidence2_Bpod_mat_210202_Feb27-2021_Session1.mat
switch nargin
    case 0
        [filename1, pathname1]=uigetfile('*.mat','Block_mat');
        filename1 = [pathname1, filename1];
        load(filename1);
        [filename2, pathname2]=uigetfile('*.mat','RL_mat');
        filename2 = [pathname2, filename2];
        data = load(filename2);
    case 1
        hoge
    otherwise
        hoge
end

[minD_trial,maxD_trial,Choice_trial,tone_evidence,~,~,~,~,~,~,~,~,~,~,binary_tone] ...
 = Dual_get_basic_task_structure_20210204(filename1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make trial type vector
Tone_type = nan(length(Correct_side),1);
Tone_type(maxD_trial) = 1;
Tone_type(minD_trial) = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Take out chose success trials:

%Use all trials for analysis
N_trial = length(Choice_trial);
Correct_trial = Correct_side(Choice_trial);
Chosen_trial = Chosen_side(Choice_trial);
binary_tone = binary_tone(Choice_trial);
Reward_LCR = Reward_LCR(Choice_trial,[1,3]); %choose chosen trials
Tone_type = Tone_type(Choice_trial);
TrialBlock = TrialBlock(Choice_trial);
StimDuration = StimDuration(Choice_trial);
minD_trial = find(StimDuration == min(StimDuration));
maxD_trial = find(StimDuration == max(StimDuration));

clear init_Q
init_Q = BlockReward(1,:); %2.4 2.4
    
%reward amount for each trial

reward_trial = zeros(N_trial,1);
for i = 1:N_trial
    %Reward trials
    if Correct_trial(i) == Chosen_trial(i)
        reward_trial(i) = Reward_LCR(i,Chosen_trial(i)+1);
    end
end

%Parameter
% bin_x = [0:0.001:1];
% stim_LR   = [0,0.25,0.45,0.55,0.75,1]; %3stim
% stim_prob_LR   = [0.5,0.25,0.25,0.25,0.25,0.5]; %3stim
% stim_prob_kitei = [0.25,0.125,0.125,0.125,0.125,0.25]; %3stim
%Sense_std = [0.01:0.05:1];
%reward = [3,1]; %left-right
%stim_prob = [0.5, 0.5];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
%Based on the binary tone decide the pseudo tone evidence
temp_tone(1).matrix = find(binary_tone == 0);
temp_tone(2).matrix = find(binary_tone > 0 & binary_tone <= 0.35);
temp_tone(3).matrix = find(binary_tone > 0.35 & binary_tone < 0.5);
temp_tone(4).matrix = find(binary_tone > 0.5 & binary_tone < 0.65);
temp_tone(5).matrix = find(binary_tone >= 0.65 & binary_tone < 1);
temp_tone(6).matrix = find(binary_tone == 1);

new_tone_evi = nan(length(binary_tone),1);
for i = 1:6
    new_tone_evi(temp_tone(i).matrix) = tone_evidence(i);
end
%Check nan
temp = unique(isnan(new_tone_evi));
temp = find(temp == 1);
if ~isempty(temp)
    [binary_tone, new_tone_evi]
    disp('nan detected')
    hoge
end
trial_evidence = new_tone_evi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get block structure
use_trial2 = [];
use_trial3 = [];
for i = 2:max(TrialBlock)
    temp = find(TrialBlock == i);
    if rem(i,2) == 0
        use_trial2 = [use_trial2; temp];
    else
        use_trial3 = [use_trial3; temp];
    end
end
if BlockReward(2,1) < BlockReward(2,2) % Right -> Left
    block_R = use_trial2;
    block_L = use_trial3;
else % Left -> Right
    block_L = use_trial2;
    block_R = use_trial3;
end
min_R = intersect(minD_trial, block_R);
min_L = intersect(minD_trial, block_L);
max_R = intersect(maxD_trial, block_R);
max_L = intersect(maxD_trial, block_L);


use_parameter = 3; 
temp_use_para = [1 3 4 5 6 7];
para_max = data.para_max(use_parameter,:);

[~,~,Q,~,~]=Thre_update_201108_max(para_max,reward_trial,Chosen_trial,binary_tone,Tone_type,init_Q,temp_use_para);

%% fig 2b left bottom
figure;
cd('G:\upload_code\Figure2\Fig2b');
sdata1 = struct();% source data 
hold on
plot(Q(:,1),'b')
plot(Q(:,2),'r')
set(gca,'xlim',[1 N_trial])
set(gca,'xtick',0:40:N_trial)
set(gca,'ylim',[0 4])
sdata1.Toneforleft = Q(:,1);
sdata1.Toneforright = Q(:,2);
T = struct2table(sdata1);
writetable(T, 'source fig2b left bottom.csv');

%%
disp('start simulation')
%Simulate choice and make psychometric function
parfor i = 1:100
    [i,100]
    
    [sim_action(i,:),Q] = ...
        Thre_update_201108_simulation(para_max, Reward_LCR, binary_tone, Tone_type, init_Q, temp_use_para);
    
    [opt_L_min(i,:),opt_R_min(i,:),opt_L_max(i,:),opt_R_max(i,:),right_minL(i,:),...
        right_minR(i,:),right_maxL(i,:),right_maxR(i,:),right_minL05] = ...
        simulated_choice_full_psycho2(sim_action(i,:),trial_evidence,tone_evidence, ...
        binary_tone, min_L,min_R,max_L,max_R);

    [left_prob2(i,:), right_prob2(i,:),...
     left_prob_long(i,:), right_prob_long(i,:), ...
     left_prob_short(i,:),right_prob_short(i,:)] = ...
     gauss_stim_choice(Correct_trial, sim_action(i,:)', StimDuration);
end
delete(gcp('nocreate'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot choice behavior
%Find block change
temp = TrialBlock(2:length(TrialBlock)) - TrialBlock(1:length(TrialBlock)-1);
Blockchange = find(temp ~= 0) + 0.5; %start of new block
left_reward = find(Correct_trial' == 0);
right_reward = find(Correct_trial' == 1);

%% fig 2b left center
y_limit = [0 1];
figure; hold on
sdata1 = struct();% source data 
sdata2 = struct();% source data 
for i = 1:length(Blockchange)
    line([Blockchange(i),Blockchange(i)],y_limit,'color',[0 0 0])
end
[mean_trace1,std_trace1,~,x1] = plot_mean_se_moto_x_axis(left_prob2,left_reward,[0 0 1],1);
[mean_trace2,std_trace2,~,x2] = plot_mean_se_moto_x_axis(right_prob2,right_reward,[1 0 0],1);
set(gca,'xlim',[1,length(Choice_trial)],'ylim',y_limit)
set(gca,'xtick',0:40:length(Choice_trial))
sdata1.ToneforLeft_x = x1';
sdata1.ToneforLeft_y = mean_trace1';
sdata1.ToneforLeft_sd= std_trace1';
sdata2.ToneforRight_x = x2';
sdata2.ToneforRight_y = mean_trace2';
sdata2.ToneforRight_sd= std_trace2';
T = struct2table(sdata1);
writetable(T, 'source fig2b left center1.csv');
T = struct2table(sdata2);
writetable(T, 'source fig2b left center2.csv');


%% Plot psychometric
evi_x = 0:0.01:1;
figure
sdata1 = struct();% source data 
sdata2 = struct();% source data 
subplot(1,2,1); hold on;%Max trials
[mean_trace1,std_trace1,~,x1] =plot_mean_se_moto_x_axis(opt_L_max,evi_x,[0 0 1],1);
[mean_trace2,std_trace2] =plot_mean_se_moto_x_axis(opt_R_max,evi_x,[1 0 0],1);
set(gca,'xlim',[-0.1 1.1],'ylim',[0 1])
sdata1.x=x1';
sdata1.leftblock_y=mean_trace1';
sdata1.leftblock_sd=std_trace1';
sdata1.rightblock_y=mean_trace2';
sdata1.rightblock_sd=std_trace2';
T = struct2table(sdata1);
writetable(T, 'source fig2b right long.csv');


subplot(1,2,2); hold on%Max trials
[mean_trace1,std_trace1,~,x1] =plot_mean_se_moto_x_axis(opt_L_min,evi_x,[0 0 1],1);
[mean_trace2,std_trace2] =plot_mean_se_moto_x_axis(opt_R_min,evi_x,[1 0 0],1);
set(gca,'xlim',[-0.1 1.1],'ylim',[0 1])
sdata2.x=x1';
sdata2.leftblock_y=mean_trace1';
sdata2.leftblock_sd=std_trace1';
sdata2.rightblock_y=mean_trace2';
sdata2.rightblock_sd=std_trace2';
T = struct2table(sdata2);
writetable(T, 'source fig2b right short.csv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [left_prob2, right_prob2,...
          left_prob_long, right_prob_long, ...
          left_prob_short,right_prob_short] = gauss_stim_choice(Correct_trial, Chosen_trial, StimDuration)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Chosen_trial = Chosen_trial - 1;

StimLength = unique(StimDuration);
StimCategory = ones(1,length(StimDuration));
StimCategory(StimDuration == max(StimLength)) = 2;
short_trial = find(StimCategory == 1);
long_trial = find(StimCategory == 2);
%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot the correct rate of short and long stimulus

%Make plot for the choice probability
%Get the closest 9 trials to calculate the probability of left or right
%choice probability
left_reward = find(Correct_trial == 0);
right_reward = find(Correct_trial == 1);
left_reward_short = intersect(left_reward, short_trial);
right_reward_short = intersect(right_reward, short_trial);
left_reward_long = intersect(left_reward, long_trial);
right_reward_long = intersect(right_reward, long_trial);
std_norm2 = 10;
%Chosen_trial
%Gaussian is required

left_prob2 = get_gauss_correct_rate(left_reward,Correct_trial,Chosen_trial,std_norm2);
right_prob2 = get_gauss_correct_rate(right_reward,Correct_trial,Chosen_trial,std_norm2);

left_prob_long = get_gauss_correct_rate(left_reward_long,Correct_trial,Chosen_trial,std_norm2);
right_prob_long = get_gauss_correct_rate(right_reward_long,Correct_trial,Chosen_trial,std_norm2);
left_prob_short = get_gauss_correct_rate(left_reward_short,Correct_trial,Chosen_trial,std_norm2);
right_prob_short = get_gauss_correct_rate(right_reward_short,Correct_trial,Chosen_trial,std_norm2);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function left_prob = get_gauss_correct_rate(left_reward,Correct_trial,Chosen_trial,std_norm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Correct_trial = Correct_trial(left_reward);
Chosen_trial = Chosen_trial(left_reward);
%Get correct trial
correct_side = Chosen_trial == Correct_trial;
correct_side = double(correct_side);

for i = 1:length(left_reward)
    temp_x = left_reward - left_reward(i);
    filter_left = normpdf(temp_x,0,std_norm);
    filter_left = filter_left ./ sum(filter_left);
    left_prob(i) = sum(correct_side .* filter_left);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [opt_L_min,opt_R_min,opt_L_max,opt_R_max,...
    right_minL,right_minR,right_maxL,right_maxR, ...
    right_minL05,right_minR05,right_maxL05,right_maxR05] = ...
    simulated_choice_full_psycho2(simu_choice,trial_evidence, tone_evidence, ...
    binary_tone, min_L,min_R,max_L,max_R)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Make artificial choice based on the Choice_prob
simu_choice = simu_choice - 1;

[right_minL, right_minL05] = get_right_prob(trial_evidence,simu_choice,tone_evidence,min_L);
[right_minR, right_minR05] = get_right_prob(trial_evidence,simu_choice,tone_evidence,min_R);
[right_maxL, right_maxL05] = get_right_prob(trial_evidence,simu_choice,tone_evidence,max_L);
[right_maxR, right_maxR05] = get_right_prob(trial_evidence,simu_choice,tone_evidence,max_R);
%%%%%%%%%%%%%%%%%

[opt_L_min,opt_R_min] = get_Gauss_standard_fit2(min_L, min_R, ...
    simu_choice, binary_tone, right_minL, right_minR);
[opt_L_max,opt_R_max] = get_Gauss_standard_fit2(max_L, max_R, ...
    simu_choice, binary_tone, right_maxL, right_maxR);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [min_right,min_right05] = get_right_prob(trial_evidence,simu_choice,tone_evidence,minD_trial)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Min trials
min_evi = trial_evidence(minD_trial);
min_choice = simu_choice(minD_trial);
for i = 1:length(tone_evidence)
    temp_choice = min_choice(min_evi == tone_evidence(i));    
    [min_right(i),min_right05(i,:)] = binofit(sum(temp_choice),length(temp_choice));
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ave_likeli,likelihood,Q,Choice_prob,para] = ...
    Thre_update_201108_max(para, reward, action, Stim_freq, Tone_type, init_Q, use_para)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

clear Prior Q_left Q_right Choice_prob
%Belief_basis for left choice reward probability
Q = zeros(N_trial,2);
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
    
    if softmax_check
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
    temp_Q_stim(1) = 2*Q(i,1) .* left_conf; %update Q with stim prob
    temp_Q_stim(2) = 2*Q(i,2) .* right_conf;  %update Q with stim prob
    
    %Forgetting Q learning
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
else
    post_prob = [(1-para_forget) * Q(:,1), temp_Q];
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sim_action,Q,Choice_prob,para] = ...
    Thre_update_201108_simulation(para, reward_LR, Stim_freq, Tone_type, init_Q, use_para)
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
N_trial = length(Stim_freq);
binary_tone = ones(N_trial,1);
binary_tone(Stim_freq > 0.5) = 2;

clear Prior Q_left Q_right Choice_prob
%Belief_basis for left choice reward probability
Q = zeros(N_trial,2);
Q(1,:)  = init_Q .* 0.5; %Prior with 0.5 in both side
Choice_prob = zeros(N_trial,2);

for i = 1:N_trial
    clear Posterior_LR Decision_LR
    %Tone_type: 1(long) or 2(short)
    temp_alpha = para(1); %1
    temp_forget = para(2); %1
    temp_std = para(Tone_type(i) + 2); %3 or 4
    temp_bias = para(Tone_type(i) + 4); %5 or 6
    temp_soft = para(Tone_type(i) + 6); %7 or 8
    
    if softmax_check
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
    if left_choice > rand
        sim_action(i) = 1; %left
    else
        sim_action(i) = 2; %right
    end
    
    %Confidence
    left_conf  = normcdf(0.5,Stim_freq(i),temp_std) - normcdf(0,Stim_freq(i),temp_std);
    right_conf = normcdf(1,Stim_freq(i),temp_std) - normcdf(0.5,Stim_freq(i),temp_std);
    left_conf = left_conf / (left_conf + right_conf);
    right_conf = 1 - left_conf;
 
    Choice_prob(i,:) = [left_choice,right_choice];
    %Based on the Choice prob, decide sim_action
    %Get probability from action and reward;
    %likelihood(i) = Choice_prob(i,action(i));
    
    %calculate the Q_stim
    clear temp_Q_stim
    temp_Q_stim(1) = 2*Q(i,1) .* left_conf; %update Q with stim prob
    temp_Q_stim(2) = 2*Q(i,2) .* right_conf;  %update Q with stim prob

    %Forgetting Q learning
    if sim_action(i) == binary_tone(i)
        temp_reward = reward_LR(i,sim_action(i));
    else
        temp_reward = 0;
    end
    Q(i+1,:) = Qlearning_confidence2(Q(i,:), temp_Q_stim, sim_action(i), temp_reward, temp_alpha, temp_forget);
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stim_prob = softmax(temp_stim,beta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear stime_prob
stim_prob(1) = exp(beta * temp_stim(1)) / (exp(beta * temp_stim(1)) + exp(beta * temp_stim(2)));
stim_prob(2) = exp(beta * temp_stim(2)) / (exp(beta * temp_stim(1)) + exp(beta * temp_stim(2)));

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [opt_L,opt_R,X_L,X_R] = get_Gauss_standard_fit2(block_L, block_R, ...
    Chosen_side, trial_evidence, right_prob_L, right_prob_R)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

evi_x = 0:0.01:1;

lapse_L = [right_prob_L(1), 1-right_prob_L(6)]; %limit for lapse
lapse_R = [right_prob_R(1), 1-right_prob_R(6)]; %limit for lapse

opt = optimset('Display','off');
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

for i = 1:10
    [X_L(i,:),FCAL_L(i)] = fminsearch(@Opt_psychometric_Gauss,para(i,:),opt,Chosen_side(block_L), trial_evidence(block_L), evi_x, lapse_L);
    [X_R(i,:),FCAL_R(i)] = fminsearch(@Opt_psychometric_Gauss,para(i,:),opt,Chosen_side(block_R), trial_evidence(block_R), evi_x, lapse_R);
    %data fitting is not least square but with maximum likelihood
end

min_L = find(FCAL_L == min(FCAL_L),1);
min_R = find(FCAL_R == min(FCAL_R),1);
X_L = X_L(min_L,:);
X_R = X_R(min_R,:);

[log_likelli,opt_L,X_L] = Opt_psychometric_Gauss_max(X_L, Chosen_side(block_L), trial_evidence(block_L), evi_x, lapse_L);
[log_likelli,opt_R,X_R] = Opt_psychometric_Gauss_max(X_R, Chosen_side(block_R), trial_evidence(block_R), evi_x, lapse_R);

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

if length(temp1)+length(temp2) ~= N_trial
    [length(temp1),length(temp2),N_trial]
    hoge
end

%likelihood keisan
log_likeli = sum(log(likelihood));
log_likeli = -log_likeli;

%get the tuning function with evi_x
temp_right = normcdf(1, evi_x, para(2)) - normcdf(para(1), evi_x, para(2)); %Right choice probablity
temp_left  = normcdf(para(1), evi_x, para(2)) - normcdf(0, evi_x, para(2)); %Left choice probablity
temp_gauss = temp_right ./ (temp_right + temp_left); %normalized with truncated part
neurometric = para(3) + (1-para(3)-para(4)) .* temp_gauss;

return


