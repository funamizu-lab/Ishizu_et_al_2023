%{
----------------------------------------------------------------------------
Analyzing behavioral data
At least for the correct rate
----------------------------------------------------------------------------
%}

function Dual_ML_CumGauss_Fit_20220114_lapse

[filename1, pathname1]=uigetfile('*.mat','Block_mat,','MultiSelect','on');

%CV_repeat = 10;
CV_repeat = 1;
%CV_para = 10;
CV_para = 1; %80% 20%

clear ave_likeli BIC log_likeli para_max
for filecount = 1 : length(filename1)
        
    clear temp_filename temp_pass fpath
    temp_filename = filename1(filecount); 
    temp_filename = cell2mat(temp_filename);
    temp_path = pathname1;
    fpath = fullfile(temp_path, temp_filename);
        
    [ave_likeli(filecount,:),BIC(filecount,:),log_likeli(filecount,:),para_max(filecount).matrix] = ...
        behave_analysis_block1_20210114_CumGauss_lapse(fpath,CV_para,CV_repeat);
    close all
    filename{filecount} = temp_filename;
end
delete(gcp('nocreate'))
save ML_CumGauss_Fit_20220114_lapse ave_likeli BIC log_likeli para_max filename

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ave_ave_likeli,ave_BIC_all,ave_log_likeli,para_max] = behave_analysis_block1_20210114_CumGauss_lapse(filename1,CV_para,CV_repeat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[filename1, pathname1]=uigetfile('*.mat','Block_mat');
load(filename1)

[minD_trial,maxD_trial,Choice_trial,tone_evidence,trial_evidence,use_trial2,use_trial3,use_trial_all,...
low,high,correct,error,flip_tone,number_use_trial,...
binary_tone,right_trial_all,number_trial_all,right_trial,number_trial] ...
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
block_LR = [block_L; block_R];
block_LR = sort(block_LR);

Correct_R = Correct_side(block_R);
Chosen_R = Chosen_side(block_R);
Evi_R = binary_tone(block_R);

%Flip the block_L to adjust to prefer
Correct_L = Correct_side(block_L);
Chosen_L = Chosen_side(block_L);
Evi_L = binary_tone(block_L);
% Correct_L = double(~Correct_L);
% Chosen_L = double(~Chosen_L);
% Evi_L = 1-Evi_L;

Correct_side = [Correct_R; Correct_L];
use_choice23 = [Chosen_R; Chosen_L];
use_evidence23 = [Evi_R; Evi_L];

bias2 = ones(length(block_R),1) * -1;
bias3 = ones(length(block_L),1);
bias23 = [bias2;bias3];
sense23 = bias23;

StimDuration = ones(length(binary_tone),1);
StimDuration(minD_trial) = -1;
StimDuration = [StimDuration(block_R); StimDuration(block_L)];
sense_Stim = StimDuration;

temp_thre  = ones(length(block_R)+length(block_L),1);
temp_sense  = ones(length(block_R)+length(block_L),1);
temp_lapse = ones(length(block_R)+length(block_L),1);
bias_lapse = bias23;
temp_x = [temp_thre, temp_sense, bias23, sense23, temp_lapse, temp_lapse];

use_para = [1 1 0 0 0 0; %1
            1 1 1 0 0 0; %2
            1 1 1 1 0 0; %3
            1 1 1 0 1 0; %4
            1 1 1 0 0 1; %5
            1 1 1 1 1 0; %6
            1 1 1 1 0 1; %7
            1 1 1 1 1 1; %8
            ];         
        
[size_y,size_x] = size(use_para);
previous_parameter(1).matrix = [];
previous_parameter(2).matrix = 1;
previous_parameter(3).matrix = 2;
previous_parameter(4).matrix = 2;
previous_parameter(5).matrix = 2;
previous_parameter(6).matrix = 3;
previous_parameter(7).matrix = 3;
previous_parameter(8).matrix = [6,7];
% previous_parameter(6).matrix = 2;
% previous_parameter(7).matrix = 2;
% previous_parameter(8).matrix = [6,7];
% previous_parameter(9).matrix = [3,6];
% previous_parameter(10).matrix = [3,7];
% previous_parameter(11).matrix = [8,9,10];
% previous_parameter(12).matrix = [4,6];
% previous_parameter(13).matrix = [4,7];
% previous_parameter(14).matrix = [8,12,13];
% previous_parameter(15).matrix = [5,9,12];
% previous_parameter(16).matrix = [5,10,13];
% previous_parameter(17).matrix = [11,14,15,16];

clear ave_ave_likeli ave_BIC_all ave_log_likeli para_max
for i = 1:size_y
    [i,size_y]
    %para search using previous parameters
    temp_para = previous_parameter(i).matrix;
    if isempty(temp_para)
        temp_para = zeros(1,6);
    else
        temp_likeli = find(ave_ave_likeli(temp_para) == max(ave_ave_likeli(temp_para)),1);
        temp_likeli = temp_para(temp_likeli);
        temp_para = para_max(temp_likeli,:);
    end
    non_zero_para = length(find(use_para(i,:) == 1));

    clear ave_likeli likelihood BIC_all log_likeilihood
    for j = 1:CV_repeat
        [ave_likeli(j,1),likelihood] = Cross_vali_20201105(temp_x,use_choice23,use_evidence23,use_para(i,:),temp_para,CV_para);
        %[ave_likeli(j,1),likelihood(:,j),~] = Cross_vali_20201027(temp_x,use_choice23,use_para(i,:),temp_para,CV_para);
        BIC_all(j,1) = -2 * sum(log(likelihood)) + non_zero_para * log(length(use_choice23));
        log_likeilihood(j,1) = sum(log(likelihood));
    end
    %Get parameter with all data
    temp_para_max = Train_para_CumGauss_20201105(temp_x,use_choice23,use_evidence23,use_para(i,:),temp_para);
    [~,~,train_para] = Opt_psychometric_max_tokyo_20201106_lapse(temp_para_max, temp_x, use_choice23,use_evidence23,use_para(i,:));
    %Update the use_para with real parameter
    para_max(i,:) = train_para;
    ave_ave_likeli(1,i) = mean(ave_likeli,1);
    ave_BIC_all(1,i) = mean(BIC_all,1);
    ave_log_likeli(1,i) = mean(log_likeilihood,1);
end


return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ave_likeli,likelihood] = Cross_vali_20201105(temp_x,use_choice23,evidence23,use_para,init_para,CV_para)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Determine parameter using training data
length_trial = length(use_choice23);
number_trial = [1:length_trial];

if length_trial > CV_para
    temp_group = rem(number_trial,CV_para) + 1;
    temp_group = temp_group(randperm(length_trial));
else
    temp_group = number_trial;
end
%temp_group

likelihood = nan(length_trial,1);
for i = 1:CV_para
%     temp_test = find(temp_group == i);
%     temp_train = setdiff(number_trial, temp_test);
%     
%     train_x = temp_x(temp_train,:);
%     train_choice = use_choice23(temp_train,:);
%     test_x = temp_x(temp_test,:);
%     test_choice = use_choice23(temp_test,:);
    
    %Train parameter
    para_max = Train_para_CumGauss_20201105(temp_x,use_choice23,evidence23,use_para,init_para);
%     para_max = Train_para_CumGauss(train_x,train_choice,use_para,init_para);
    %Get maximum prediction with test data
    [~,likelihood] = Opt_psychometric_max_tokyo_20201106_lapse(para_max, temp_x, use_choice23,evidence23,use_para);
%     [~,temp_likelihood] = Opt_psychometric_max_tokyo(para_max, test_x, test_choice, use_para);
%     likelihood(temp_test) = temp_likelihood;
%     %Update the use_para with real parameter
%     temp_para = find(use_para ~= 0);
%     use_para(temp_para) = para_max;
end

%likelihood keisan
log_likeli = sum(log(likelihood));
log_likeli = log_likeli / length_trial;
ave_likeli = exp(log_likeli);
%ave_likeli = -ave_likeli;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function para_max = Train_para_CumGauss_20201105(temp_x,use_choice23,evidence23,use_para,init_para)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Log(likelihood) = Tn*logYn + (1-Tn)*log(1-Yn)
%Yn = 1/(1+exp(-bx))
%Tn is the answer
%Above equations are from Bishopè„, pp. 205
para_search = 50;

number_para = find(use_para == 1);
init_para = init_para(number_para);
init_para0 = find(init_para == 0);

opt = optimset('Display','off');
%for i = 1 : 3,
parfor i = 1 : para_search
    if i < para_search/2
        rand_para = init_para;
        rand_para(init_para0) = 0.1 * rand(1,length(init_para0));
    else
        %rand_para = 0.1 * rand(1,length(number_para));
        rand_para = 0.5 + 0.2 * (rand(1,length(number_para))-0.5);
    end
    
    [X_all(i,:),FCAL_all(i),EXITFLAG,OUTPUT] = fminsearch(@Opt_psychometric_20201105,rand_para,opt,...
        temp_x, use_choice23, evidence23, use_para);
    
%     rand_para = rand(1,length(init_para0));
%     temp_para = init_para;
%     temp_para(init_para0) = rand_para;
    %Detect non-real vale
    if isreal(FCAL_all(i)) == 0,
        FCAL_all(i) = 0;
    end
end
temp_FCAL = find(FCAL_all == min(FCAL_all),1);
para_max = X_all(temp_FCAL,:);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ave_likeli,likelihood,temp_para] = Opt_psychometric_max_tokyo_20201106_lapse(para, temp_x, use_choice23, evidence23, use_para)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

para_non_zero = find(use_para ~= 0);
temp_para = zeros(1,length(use_para));
temp_para(para_non_zero) = para;
if temp_para(5) < 0
    temp_para(5) = 0;
end
if temp_para(6) < 0
    temp_para(6) = 0;
end
if temp_para(5) + temp_para(6) > 1
    temp_para(5) = temp_para(5) ./ (temp_para(5) + temp_para(6));
    temp_para(6) = temp_para(6) ./ (temp_para(5) + temp_para(6));
end

%temp_x = [temp_thre, temp_sense, bias23, sense23, temp_lapse, temp_lapse];
%Average for para 1 3 5
%Std para 2 4 6

[N_trial,size_x] = size(temp_x);

likelihood = zeros(1,N_trial);

para_data = repmat(temp_para,N_trial,1);

mean_para = para_data(:,[1,3]) .* temp_x(:,[1,3]);
std_para  = para_data(:,[2,4]) .* temp_x(:,[2,4]);
mean_para = sum(mean_para,2);
std_para = sum(std_para,2);
std_para = abs(std_para);

temp_right = normcdf(1,evidence23,std_para) - normcdf(mean_para,evidence23,std_para);
temp_left = normcdf(mean_para,evidence23,std_para) - normcdf(0,evidence23,std_para);
temp_exp = temp_right ./ (temp_right + temp_left);

%Without lapse
if use_para(5) ~= 0 && use_para(6) ~= 0
    temp_p = temp_para(5) + (1-temp_para(5)-temp_para(6)) .* temp_exp;
elseif use_para(5) ~= 0
    temp_p = temp_para(5) + (1-2.*temp_para(5)) .* temp_exp;
elseif use_para(6) ~= 0 %non-prefer lapse
    block_R = find(temp_x(:,3) == -1);
    block_L = find(temp_x(:,3) == 1);
    temp_p = temp_exp;
    temp_p(block_R) = temp_para(6) + (1-temp_para(6)) .* temp_exp(block_R);
    temp_p(block_L) = (1-temp_para(6)) .* temp_exp(block_L);
else
    temp_p = temp_exp;
end

temp1 = find(use_choice23 == 0);
temp2 = find(use_choice23 == 1);
likelihood(temp1) = 1-temp_p(temp1);
likelihood(temp2) = temp_p(temp2);

if length(temp1)+length(temp2) ~= N_trial
    [length(temp1),length(temp2),N_trial]
    hoge
end

%likelihood keisan
log_likeli = sum(log(likelihood));
log_likeli = log_likeli / N_trial;
ave_likeli = exp(log_likeli);
ave_likeli = -ave_likeli;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ave_likeli_mean = Opt_psychometric_20201105(para, temp_x, use_choice23, evidence23, use_para)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ave_likeli_mean,~] = Opt_psychometric_max_tokyo_20201106_lapse(para, temp_x, use_choice23, evidence23, use_para);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_stim_choice(Outcome, Correct_side, Chosen_side, EvidenceStrength, use_evidence, trials)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Choice_trial = find(Outcome == 1 | Outcome == 2);
Choice_trial = intersect(Choice_trial, trials); %trials selection
Correct_trial = Correct_side(Choice_trial);
Chosen_trial = Chosen_side(Choice_trial);
Evidence_trial = EvidenceStrength(Choice_trial);
%stim_strength
for i = 1:length(use_evidence)-1,
    temp = find(Evidence_trial > use_evidence(i) & Evidence_trial < use_evidence(i+1));
    Evidence_stim(temp) = i;
end

length(Choice_trial)

trial_color2 = [0,0,1;1,0,0];
figure
for i = 1:2,
    temp_trial = find(Correct_trial == i-1); %0 Left 1 Right
    temp_stim = Evidence_stim(temp_trial);
    temp_choice = Chosen_trial(temp_trial);
    
    for j = 1:length(temp_trial),
        temp_x = [temp_trial(j), temp_trial(j)];
        temp_y = [1,0]-temp_choice(j);
        temp_y = temp_y * temp_stim(j); %Bar hight shows the difficulty
        line(temp_x, temp_y,'color',trial_color2(i,:),'LineWidth',1);
        hold on
    end
end
set(gca,'xlim',[1,length(Choice_trial)],'ylim',[-4,4])
%set(gca,'xtick',[1,60,100,200,210,300,360,400,500])
set(gca,'xtick',[0:40:length(Choice_trial)])

% trial_color = [0, 176/255, 240/255; 51/255, 51/255, 1; 0, 0, 102/255; ...
%                1, 102/255, 0;       192/255, 0, 0;     102/255, 0, 51/255];
% figure
% for i = 1:2,
%     temp_trial = find(Correct_trial == i-1); %0 Left 1 Right
%     temp_stim = Evidence_stim(temp_trial);
%     temp_choice = Chosen_trial(temp_trial);
%     
%     for j = 1:length(temp_trial),
%         temp_x = [temp_trial(j), temp_trial(j)];
%         temp_y = [1,0]-temp_choice(j);
%         temp_color = 3 * (i-1) + temp_stim(j);
%         line(temp_x, temp_y,'color',trial_color(temp_color,:),'LineWidth',1);
%         hold on
%     end
% end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_psycho_curve2(p_right_all, conf_all, p_fit_all, tone_evidence, x_evi_plot, plot_color)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(tone_evidence', p_right_all','.','color',plot_color) %Mean
for i = 1:length(tone_evidence),
    hold on
    plot([tone_evidence(i),tone_evidence(i)], [conf_all(i,1),conf_all(i,2)],'-','color',plot_color)
end
hold on
plot(x_evi_plot,p_fit_all,'-','LineWidth',1,'color',plot_color)

set(gca,'xlim',[-0.1 1.1])
set(gca,'ylim',[0 1])
%set(gca,'xtick',[0:0.05:1])
set(gca,'xtick',[0:0.1:1])
set(gca,'ytick',[0:0.1:1])
xlabel('[high tones - low tones]/s')
ylabel('Fraction rightward')

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_psycho_curve(p_right_all, conf_all, p_fit_all, tone_evidence, plot_color)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(tone_evidence', p_right_all','.','color',plot_color) %Mean
for i = 1:length(tone_evidence),
    hold on
    plot([tone_evidence(i),tone_evidence(i)], [conf_all(i,1),conf_all(i,2)],'-','color',plot_color)
end
hold on
plot(tone_evidence',p_fit_all,'-','LineWidth',1,'color',plot_color)

set(gca,'xlim',[-0.1 1.1])
set(gca,'ylim',[0 1])
%set(gca,'xtick',[0:0.05:1])
set(gca,'xtick',[0:0.1:1])
set(gca,'ytick',[0:0.1:1])
xlabel('[high tones - low tones]/s')
ylabel('Fraction rightward')

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Select use trials and make figures
function [choice_stim, right_number_trial, number_trial] = ...
    psycho_plot2(Outcome, Correct_side, Chosen_side, EvidenceStrength, use_evidence, trials)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Select choice
%Center lick error should not count for the block trial!!
Choice_trial = find(Outcome == 1 | Outcome == 2);
Choice_trial = intersect(Choice_trial, trials); %trials selection

Correct_trial = Correct_side(Choice_trial);
Chosen_trial = Chosen_side(Choice_trial);
Evidence_trial = EvidenceStrength(Choice_trial);

length(Choice_trial)
clear right_prob
clear choice_stim
choice_stim = [];
for i = 1:2,
    temp = find(Correct_trial == i-1);
    temp_select = Chosen_trial(temp);
    temp_evidence = Evidence_trial(temp);
    
    for j = 1:length(use_evidence)-1,
        temp_trial = find(temp_evidence > use_evidence(j) & temp_evidence < use_evidence(j+1));
        temp_trial = temp_select(temp_trial);
        temp_correct = find(temp_trial == 1);
        right_prob(i,j) = length(temp_correct) / length(temp_trial);
        right_trial(i,j) = length(temp_correct);
        number_trial(i,j) = length(temp_trial);
        
        %Keep the choice stim record
        temp_stim = (length(use_evidence)-1) * (i-1) + j;
        temp_choice_stim = [temp_trial, ones(length(temp_trial),1) * temp_stim];
        choice_stim = [choice_stim; temp_choice_stim];
    end
end

plot_right_prob = [right_prob(1,3),right_prob(1,2),right_prob(1,1),right_prob(2,1),right_prob(2,2),right_prob(2,3)]
%Left reward to Right reward
right_number_trial = [right_trial(1,3),right_trial(1,2),right_trial(1,1),right_trial(2,1),right_trial(2,2),right_trial(2,3)]
number_trial = [number_trial(1,3),number_trial(1,2),number_trial(1,1),number_trial(2,1),number_trial(2,2),number_trial(2,3)]
% figure
% plot(plot_left_prob)
% set(gca,'xlim',[0 7],'ylim',[0 1])

return