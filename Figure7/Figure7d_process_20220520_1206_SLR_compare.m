
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
function Figure7d_process_20220520_1206_SLR_compare

close all
[auc_correct_s, auc_correct_c, auc_correct_p] = process_20221206_SLR('auc_folders');
[mpfc_correct_s,mpfc_correct_c,mpfc_correct_p]= process_20221206_SLR('mpfc_folders');
[fof_correct_s, fof_correct_c, fof_correct_p] = process_20221206_SLR('fof_folders');


prior_frame = 1;
sound_on_frame = 2;
sound_long_off_frame = 4;

%% Fig 6d
cd('G:\upload_code\Figure7\Fig7d');
correct_rate_y_values = [0.4 0.9];
sdata = struct();% source data 
%Prior
[auc,fof,mpfc,group] = get_p_correct_20221206(auc_correct_p, fof_correct_p, mpfc_correct_p, prior_frame, correct_rate_y_values);
sdata.prior_group=group;
sdata.prior_data =[auc;fof;mpfc];

%Sound
[auc,fof,mpfc,group] = get_p_correct_20221206(auc_correct_s, fof_correct_s, mpfc_correct_s, sound_on_frame, correct_rate_y_values);
sdata.sond_group=group;
sdata.sound_data =[auc;fof;mpfc];

%Choice
[auc,fof,mpfc,group] = get_p_correct_20221206(auc_correct_c, fof_correct_c, mpfc_correct_c, sound_long_off_frame, correct_rate_y_values);
sdata.choice_group=group;
sdata.choice_data =[auc;fof;mpfc];

T = struct2table(sdata);
writetable(T, 'source fig7d.csv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [auc,fof,mpfc,group] = get_p_correct_20221206(auc_correct_s, fof_correct_s, mpfc_correct_s, use_frame, y_axis_values)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% p(1) = ranksum(auc_correct_s(:,use_frame), fof_correct_s(:,use_frame));
% p(2) = ranksum(auc_correct_s(:,use_frame), mpfc_correct_s(:,use_frame));
% p(3) = ranksum(mpfc_correct_s(:,use_frame),fof_correct_s(:,use_frame));
% 
% median_correct = [median(auc_correct_s(:,use_frame)),median(fof_correct_s(:,use_frame)),median(mpfc_correct_s(:,use_frame))];

all_correct = [auc_correct_s(:,use_frame);fof_correct_s(:,use_frame);mpfc_correct_s(:,use_frame)];
temp_x = [ones(length(auc_correct_s(:,use_frame)),1);
          ones(length(fof_correct_s(:,use_frame)),1)*2;
          ones(length(mpfc_correct_s(:,use_frame)),1)*3];

figure
hold on
boxplot(all_correct,temp_x)
plot( (rand(length(auc_correct_s(:,use_frame)), 1)-0.5)*0.2 + 1, auc_correct_s(:,use_frame),'.')
plot( (rand(length(fof_correct_s(:,use_frame)), 1)-0.5)*0.2 + 2, fof_correct_s(:,use_frame),'.')
plot( (rand(length(mpfc_correct_s(:,use_frame)),1)-0.5)*0.2 + 3, mpfc_correct_s(:,use_frame),'.')
set(gca,'ylim',y_axis_values)

auc = auc_correct_s(:,use_frame);
fof = fof_correct_s(:,use_frame);
mpfc= mpfc_correct_s(:,use_frame);
group=temp_x;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [long_correct_s, long_correct_c, long_correct_p] = process_20221206_SLR(folders)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

analysis_dir = eval(folders);

% SLR_name = 'SLR50_20221209_glmnet_sound.mat';
% SLR_choice= 'SLR50_20221209_glmnet_choice.mat';
% SLR_prior = 'SLR50_20221209_glmnet_prior.mat';
SLR_choice = 'SLR100_20221209_glmnet_choice.mat';
SLR_prior = 'SLR100_20221209_glmnet_prior.mat';
SLR_name = 'SLR100_20221209_glmnet_sound.mat';

long_correct_s = get_average_correct_rate_20221206(analysis_dir, SLR_name);
long_correct_c = get_average_correct_rate_20221206(analysis_dir, SLR_choice);
long_correct_p = get_average_correct_rate_20221206(analysis_dir, SLR_prior);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function long_correct = get_average_correct_rate_20221206(analysis_dir, SLR_name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%SLR_20220228_glmnet 
%short_correct_rate short_ave_likeli 
%long_correct_rate long_ave_likeli use_neuron use_trial

count_l = 0;
count_s = 0;
for i = 1:length(analysis_dir)
    [i,length(analysis_dir)]

    cd(analysis_dir{i});
    load(SLR_name);
    
    if length(long_correct_rate) ~= 0
        count_l = count_l + 1;
        long_correct(count_l,:)= mean(long_correct_rate);
%         long_likeli(count_l,:) = mean(long_ave_likeli);
%         new_long_dp(count_l,:) = mean(long_dp);
%         long_trial(count_l,1) = length_trial(2);
        
        [size_CV,~] = size(long_correct_rate);

        if size_CV~=1000
            cd(analysis_dir{i})
            hoge
        end
    end
end    

return