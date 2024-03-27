%{
----------------------------------------------------------------------------
Analyzing behavioral data
At least for the correct rate
----------------------------------------------------------------------------
%}
function Figure1f_Dual_behave_analysis_block1_220818_short_long

analysis_folder = {
    'G:\Ishizu_data\IntermediateFiles\only_all_behaviors\a04_behave\'
    'G:\Ishizu_data\IntermediateFiles\only_all_behaviors\a08_behave\'
    'G:\Ishizu_data\IntermediateFiles\only_all_behaviors\i20_behave\'
    'G:\Ishizu_data\IntermediateFiles\only_all_behaviors\i24_behave\'
    'G:\Ishizu_data\IntermediateFiles\only_all_behaviors\i34_behave\'
    'G:\Ishizu_data\IntermediateFiles\only_all_behaviors\i35_behave\'
    'G:\Ishizu_data\IntermediateFiles\only_all_behaviors\i43_behave\'
    'G:\Ishizu_data\IntermediateFiles\only_all_behaviors\i46_behave\'
    };


right_minD2 = [];
right_maxD2 = [];
left_minD2 = [];
left_maxD2 = [];
easy_right = [];
easy_left = [];
dif_right = [];
dif_left = [];
for i = 1:length(analysis_folder)
    temp_folder = analysis_folder{i};
    cd(temp_folder);
    filename1 = dir('Bpod*.mat');
    
    clear prob_right_minD2 prob_right_maxD2
    clear prob_left_minD2 prob_left_maxD2
    clear prob_easy_right prob_easy_left
    clear prob_dif_right prob_dif_left
    Ndata(i,1) = length(filename1);
    
    for filecount = 1 : length(filename1)
        temp_filename = filename1(filecount).name;
        temp_file{filecount} = temp_filename;
        fpath = temp_filename;
    
        [prob_right_minD2(filecount,:),prob_right_maxD2(filecount,:),...
         prob_left_minD2(filecount,:),prob_left_maxD2(filecount,:), ...
         prob_easy_right(filecount,:),prob_easy_left(filecount,:),...
         prob_dif_right(filecount,:),prob_dif_left(filecount,:)] = ...
            get_Dual_behave_analysis_block1_220228_short_long(fpath);
    end
    
    right_minD2 = [right_minD2; prob_right_minD2];
    right_maxD2 = [right_maxD2; prob_right_maxD2];
    left_minD2 = [left_minD2; prob_left_minD2];
    left_maxD2 = [left_maxD2; prob_left_maxD2];
    
    easy_right = [easy_right; prob_easy_right];
    easy_left = [easy_left; prob_easy_left];
    dif_right = [dif_right; prob_dif_right];
    dif_left = [dif_left; prob_dif_left];
end

all_subject = [];
for i = 1:length(analysis_folder)
    temp = ones(Ndata(i),1) * i;
    all_subject = [all_subject; temp];
end

nan_R_minD2 = mean(right_minD2,2);
nan_R_minD2 = find(isnan(nan_R_minD2) == 1);
nan_R_maxD2 = mean(right_maxD2,2);
nan_R_maxD2 = find(isnan(nan_R_maxD2) == 1);
nan_L_minD2 = mean(left_minD2,2);
nan_L_minD2 = find(isnan(nan_L_minD2) == 1);
nan_L_maxD2 = mean(left_maxD2,2);
nan_L_maxD2 = find(isnan(nan_L_maxD2) == 1);

%nan session
nan_sessionR = union(nan_R_minD2,nan_R_maxD2);
nan_sessionL = union(nan_L_minD2,nan_L_maxD2);
nan_session = union(nan_sessionR,nan_sessionL);

%kentei
disp('right short long')
fitlme_analysis_20220310(right_minD2,right_maxD2,all_subject,nan_session);
disp('left short long')
fitlme_analysis_20220310(left_minD2,left_maxD2,all_subject,nan_session);

right_minD2(nan_session,:) = [];
left_minD2(nan_session,:) = [];
right_maxD2(nan_session,:) = [];
left_maxD2(nan_session,:) = [];

temp_x = [0 0.25 0.45 0.55 0.75 1];
%%% fig 1f right %%%
cd('G:\upload_code\Figure1\Fig1f');
figure; hold on
sdata1 = struct();
[rs_mean,~,rs_se]=plot_mean_se_moto_x_axis(right_minD2,temp_x,[1 0 1],2);
[ls_mean,~,ls_se]=plot_mean_se_moto_x_axis(left_minD2,temp_x,[1 0 1],2);
[rl_mean,~,rl_se]=plot_mean_se_moto_x_axis(right_maxD2,temp_x,[108 173 119]./255,2);
[ll_mean,~,ll_se]=plot_mean_se_moto_x_axis(left_maxD2,temp_x,[108 173 119]./255,2);
set(gca,'xlim',[-0.1 1.1],'ylim',[0 1])
sdata1.x=temp_x';
sdata1.RightShort_mean=rs_mean';
sdata1.RightLong_mean =rl_mean';
sdata1.LeftShort_mean =ls_mean';
sdata1.LeftLong_mean  =ll_mean';
sdata1.RightShort_se=rs_se';
sdata1.RightLong_se =rl_se';
sdata1.LeftShort_se =ls_se';
sdata1.LeftLong_se  =ll_se';
T = struct2table(sdata1);
writetable(T, 'source fig1f right.csv');


nan_R_minD2 = mean(dif_right,2);
nan_R_minD2 = find(isnan(nan_R_minD2) == 1);
nan_R_maxD2 = mean(easy_right,2);
nan_R_maxD2 = find(isnan(nan_R_maxD2) == 1);
nan_L_minD2 = mean(dif_left,2);
nan_L_minD2 = find(isnan(nan_L_minD2) == 1);
nan_L_maxD2 = mean(easy_left,2);
nan_L_maxD2 = find(isnan(nan_L_maxD2) == 1);

%nan session
nan_sessionR = union(nan_R_minD2,nan_R_maxD2);
nan_sessionL = union(nan_L_minD2,nan_L_maxD2);
nan_session = union(nan_sessionR,nan_sessionL);

%Kentei Dif and easy
disp('right easy dif')
fitlme_analysis_20220310(dif_right,easy_right,all_subject,nan_session);
disp('left easy dif')
fitlme_analysis_20220310(dif_left,easy_left,all_subject,nan_session);

dif_right(nan_session,:) = [];
dif_left(nan_session,:) = [];
easy_right(nan_session,:) = [];
easy_left(nan_session,:) = [];

%%% fig 1f left %%%
figure; hold on
sdata2 = struct();
[rd_mean,~,rd_se]=plot_mean_se_moto_x_axis(dif_right,temp_x,[1 0 0],2);
[ld_mean,~,ld_se]=plot_mean_se_moto_x_axis(dif_left,temp_x,[0 0 1],2);
[re_mean,~,re_se]=plot_mean_se_moto_x_axis(easy_right,temp_x,[1 0 1],2);
[le_mean,~,le_se]=plot_mean_se_moto_x_axis(easy_left,temp_x,[0 0.5 1],2);
set(gca,'xlim',[-0.1 1.1],'ylim',[0 1])
sdata2.x=temp_x';
sdata2.RightDiff_mean=rd_mean';
sdata2.RightEasy_mean=re_mean';
sdata2.LeftDiff_mean =ld_mean';
sdata2.LeftEasy_mean =le_mean';
sdata2.RightDiff_se=rd_se';
sdata2.RightEasy_se=re_se';
sdata2.LeftDiff_se =ld_se';
sdata2.LeftEasy_se =le_se';
T = struct2table(sdata2);
writetable(T, 'source fig1f left.csv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fitlme_analysis_20220310(right_minD2,right_maxD2,all_subject,nan_session)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

right_minD2(nan_session,:) = [];
right_maxD2(nan_session,:) = [];
all_subject(nan_session,:) = [];

temp_minD2 = mean(right_minD2,2);
temp_maxD2 = mean(right_maxD2,2);
lme_right = fitlme_analysis_20210520_0(temp_minD2-temp_maxD2,all_subject);
lme_right(1).lme
lme_right(2).lme

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [prob_right_minD2,prob_right_maxD2,prob_left_minD2,prob_left_maxD2,...
    prob_easy_right,prob_easy_left,prob_dif_right,prob_dif_left] = ...
    get_Dual_behave_analysis_block1_220228_short_long(filename1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(filename1)

[minD_trial,maxD_trial,~,tone_evidence,trial_evidence,~,~,use_trial_all,low,high,correct] ...
 = Dual_get_basic_task_structure_20210204(filename1);

%min stimulus
minD_trial = intersect(minD_trial, use_trial_all);
maxD_trial = intersect(maxD_trial, use_trial_all);

right_correct = intersect(high,correct);
right_correct = intersect(right_correct, use_trial_all);
left_correct = intersect(low,correct);
left_correct = intersect(left_correct,use_trial_all);

right_minD = intersect(right_correct,minD_trial);
right_maxD = intersect(right_correct,maxD_trial);
left_minD = intersect(left_correct,minD_trial);
left_maxD = intersect(left_correct,maxD_trial);

prob_right_minD2 = get_right_trial_prob(right_minD,use_trial_all,Chosen_side,trial_evidence,tone_evidence);
prob_right_maxD2 = get_right_trial_prob(right_maxD,use_trial_all,Chosen_side,trial_evidence,tone_evidence);
prob_left_minD2 = get_right_trial_prob(left_minD,use_trial_all,Chosen_side,trial_evidence,tone_evidence);
prob_left_maxD2 = get_right_trial_prob(left_maxD,use_trial_all,Chosen_side,trial_evidence,tone_evidence);

if length(tone_evidence) ~= 6
    unique(tone_evidence)
    hoge
end
for i = 1:6
    evi_trial(i).matrix = find(trial_evidence == tone_evidence(i));
end
easy_trial = sort([evi_trial(1).matrix; evi_trial(6).matrix]);
dif_trial = sort([evi_trial(3).matrix; evi_trial(4).matrix]);

easy_right_correct = intersect(easy_trial,right_correct);
easy_left_correct = intersect(easy_trial,left_correct);
dif_right_correct = intersect(dif_trial,right_correct);
dif_left_correct = intersect(dif_trial,left_correct);

prob_easy_right = get_right_trial_prob(easy_right_correct,use_trial_all,Chosen_side,trial_evidence,tone_evidence);
prob_easy_left = get_right_trial_prob(easy_left_correct,use_trial_all,Chosen_side,trial_evidence,tone_evidence);
prob_dif_right = get_right_trial_prob(dif_right_correct,use_trial_all,Chosen_side,trial_evidence,tone_evidence);
prob_dif_left = get_right_trial_prob(dif_left_correct,use_trial_all,Chosen_side,trial_evidence,tone_evidence);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function prob_right_minD2 = get_right_trial_prob(right_minD,use_trial_all,Chosen_side,trial_evidence,tone_evidence)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

right_minD2 = right_minD(right_minD+1 <= max(use_trial_all)) + 1;

[temp_right, temp_trial] = get_right_choice_trials(right_minD2,Chosen_side,trial_evidence,tone_evidence);
prob_right_minD2 = temp_right ./ temp_trial;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [right_trial, number_trial] = get_right_choice_trials(use_trials,Chosen_side,trial_evidence,tone_evidence)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

use_choice = Chosen_side(use_trials);
use_tone = trial_evidence(use_trials);

for i = 1:length(tone_evidence)
    temp_trial = find(use_tone == tone_evidence(i));
    temp_choice = use_choice(temp_trial);
    number_trial(i) = length(temp_trial);
    right_trial(i) = sum(temp_choice);
end

return
