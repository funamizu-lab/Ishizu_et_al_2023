%{
----------------------------------------------------------------------------
Analyzing behavioral data
At least for the correct rate
----------------------------------------------------------------------------
%}

function Figure1b_Dual_behave_analysis_block1_200706_plot5_precise2

close all
% data = a04_2021_0227
[filename1, pathname1]=uigetfile('*.mat','Block_mat');
filename1 = [pathname1, filename1];
cd(pathname1);
load(filename1);

[~,~,~,~,~,~,~,~,~,~,~,~,~,~,binary_tone,~,~,~,~] ...
    = Dual_get_basic_task_structure_20210204(filename1);

%Plot trial series
cd('G:\upload_code\Figure1\Fig1b');
plot_stim_choice(Outcome, Correct_side, Chosen_side, binary_tone, StimDuration, TrialBlock)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_stim_choice(Outcome, Correct_side, Chosen_side, trial_evidence, StimDuration, TrialBlock)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Choice_trial = find(Outcome == 1 | Outcome == 2);
Correct_trial = Correct_side(Choice_trial);
Chosen_trial = Chosen_side(Choice_trial);
trial_evidence = trial_evidence(Choice_trial);
StimDuration = StimDuration(Choice_trial);
TrialBlock = TrialBlock(Choice_trial);
%Find block change
temp = TrialBlock(2:length(TrialBlock)) - TrialBlock(1:length(TrialBlock)-1);
Blockchange = find(temp ~= 0) + 0.5; %start of new block

StimLength = unique(StimDuration);
StimCategory = ones(1,length(StimDuration));
StimCategory(StimDuration == max(StimLength)) = 2;
short_trial = find(StimCategory == 1);
long_trial = find(StimCategory == 2);

length(Choice_trial)

trial_color2 = [0.8,0,0.8; 0,0.8,0];

%%% Fig1b top %%%
figure; hold on
sdata1 = struct;% source data
%Make block lines
for i = 1:length(Blockchange)
    line([Blockchange(i),Blockchange(i)], [-0.6 0.6],'color',[0 0 0],'LineWidth',0.5);
end
for j = 1:length(Correct_trial) %0 or 1
    sdata1(j).trial=j;
    sdata1(j).y=trial_evidence(j)-0.5;
    sdata1(j).stimcategory=StimCategory(j);
    
    temp_x = [j, j];
    temp_y = [0, trial_evidence(j)-0.5];
    line(temp_x, temp_y,'color',trial_color2(StimCategory(j),:),'LineWidth',1);
end
set(gca,'xlim',[1,length(Choice_trial)],'ylim',[-0.6,0.6])
set(gca,'xtick',0:40:length(Choice_trial))
T = struct2table(sdata1);
writetable(T, 'source fig1b top.csv');

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
right_prob2= get_gauss_correct_rate(right_reward,Correct_trial,Chosen_trial,std_norm2);

left_prob_long = get_gauss_correct_rate(left_reward_long,Correct_trial,Chosen_trial,std_norm2);
right_prob_long= get_gauss_correct_rate(right_reward_long,Correct_trial,Chosen_trial,std_norm2);
left_prob_short = get_gauss_correct_rate(left_reward_short,Correct_trial,Chosen_trial,std_norm2);
right_prob_short= get_gauss_correct_rate(right_reward_short,Correct_trial,Chosen_trial,std_norm2);

y_limit = [0 1];


%%% Fig1b center bottom %%%
figure
sdata2 = struct;% source data
sdata3 = struct;% source data
subplot(2,1,1);
hold on
for i = 1:length(Blockchange)
    line([Blockchange(i),Blockchange(i)],y_limit,'color',[0 0 0])
end
plot(left_reward_short,left_prob_short,'b')
plot(right_reward_short,right_prob_short,'r')
set(gca,'xlim',[1,length(Choice_trial)],'ylim',y_limit)
set(gca,'xtick',0:40:length(Choice_trial))
sdata2.ToneforLeft_x = left_reward_short;
sdata2.ToneforLeft_y = left_prob_short';
sdata3.ToneforRight_x = right_reward_short;
sdata3.ToneforRight_y = right_prob_short';
T = struct2table(sdata2);
writetable(T, 'source fig1b center1.csv');
T = struct2table(sdata3);
writetable(T, 'source fig1b center2.csv');

sdata4 = struct;% source data
sdata5 = struct;% source data
subplot(2,1,2)
hold on
for i = 1:length(Blockchange)
    line([Blockchange(i),Blockchange(i)],y_limit,'color',[0 0 0])
end
plot(left_reward_long,left_prob_long,'b')
plot(right_reward_long,right_prob_long,'r')
set(gca,'xlim',[1,length(Choice_trial)],'ylim',y_limit)
set(gca,'xtick',0:40:length(Choice_trial))
sdata4.ToneforLeft_x = left_reward_long;
sdata4.ToneforLeft_y = left_prob_long';
sdata5.ToneforRight_x = right_reward_long;
sdata5.ToneforRight_y = right_prob_long';

T = struct2table(sdata4);
writetable(T, 'source fig1b bottom1.csv');
T = struct2table(sdata5);
writetable(T, 'source fig1b bottom2.csv');
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