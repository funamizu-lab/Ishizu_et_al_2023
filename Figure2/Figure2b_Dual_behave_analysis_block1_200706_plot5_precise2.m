%{
----------------------------------------------------------------------------
Analyzing behavioral data
At least for the correct rate
----------------------------------------------------------------------------
%}

function Figure2b_Dual_behave_analysis_block1_200706_plot5_precise2

close all
% data = a04_2021_0227
% ~\IntermediateFiles\only_all_behaviors\a04_behave\Bpod_mat_210202_Feb27-2021_Session1.mat
[filename1, pathname1]=uigetfile('*.mat','Block_mat');
filename1 = [pathname1, filename1];
load(filename1);

[~,~,~,~,~,~,~,~,~,~,~,~,~,~,binary_tone,~,~,~,~] ...
    = Dual_get_basic_task_structure_20210204(filename1);

%Plot trial series
plot_stim_choice(Outcome, Correct_side, Chosen_side, TrialBlock)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_stim_choice(Outcome, Correct_side, Chosen_side, TrialBlock)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Choice_trial = find(Outcome == 1 | Outcome == 2);
Correct_trial = Correct_side(Choice_trial);
Chosen_trial = Chosen_side(Choice_trial);
TrialBlock = TrialBlock(Choice_trial);
%Find block change
temp = TrialBlock(2:length(TrialBlock)) - TrialBlock(1:length(TrialBlock)-1);
Blockchange = find(temp ~= 0) + 0.5; %start of new block


%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot the correct rate of short and long stimulus

%Make plot for the choice probability
%Get the closest 9 trials to calculate the probability of left or right
%choice probability
left_reward = find(Correct_trial == 0);
right_reward = find(Correct_trial == 1);
std_norm2 = 10;

%Chosen_trial
%Gaussian is required
left_prob2 = get_gauss_correct_rate(left_reward,Correct_trial,Chosen_trial,std_norm2);
right_prob2= get_gauss_correct_rate(right_reward,Correct_trial,Chosen_trial,std_norm2);

y_limit = [0 1];

%%% fig 2b left top
figure;
hold on
cd('G:\upload_code\Figure2\Fig2b');
sdata1 = struct();% source data 
sdata2 = struct();% source data 
%Plot the block change
for i = 1:length(Blockchange)
    line([Blockchange(i),Blockchange(i)],y_limit,'color',[0 0 0])
end
plot(left_reward,left_prob2,'b')
plot(right_reward,right_prob2,'r')
set(gca,'xlim',[1,length(Choice_trial)],'ylim',y_limit)
set(gca,'xtick',0:40:length(Choice_trial))
sdata1.Toneforleft_x = left_reward;
sdata1.Toneforleft_y = left_prob2';
sdata2.Toneforleft_x = right_reward;
sdata2.Toneforleft_y = right_prob2';

T = struct2table(sdata1);
writetable(T, 'source fig2b left top1.csv');
T = struct2table(sdata2);
writetable(T, 'source fig2b left top2.csv');
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