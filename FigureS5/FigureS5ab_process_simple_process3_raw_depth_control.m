
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
function FigureS5ab_process_simple_process3_raw_depth_control(folders)

close all
analysis_dir = eval(folders);

sig_correct25_all17 = [];
for i = 1:40
    sig_correct_all(i).matrix = [];
    sig_error_all(i).matrix = [];
    sig_choice_all(i).matrix = [];
    sig_prior_all(i).matrix = [];
    short_correct_all(i).matrix = [];
    ROC_choice_all(i).matrix = [];
    ROC_prior_all(i).matrix = [];
    for j = 1:6
        evi_prefer_all(i,j).matrix = [];
        evi_nonprefer_all(i,j).matrix = [];
    end
end
for i = 1:length(analysis_dir)
    [i,length(analysis_dir)]
    
    %Activity
    [~,~,sig_correct, sig_error,~, sig_choice, sig_prior,~,~,~,~,~,...
        evi_prefer, evi_nonprefer,~,~,sig_correct_short,sig_correct25_from17] = ...
        Task_kaiseki_tokyo1_20221122_sound_choice_process3_raw_depth(analysis_dir{i});
    
    sig_correct25_all17 = [sig_correct25_all17; sig_correct25_from17];
    
    %ROC
    [~, ~, ~, ROC_choice, ROC_prior] = ...
        Task_kaiseki_tokyo1_20230701_sound_choice_ROC(analysis_dir{i});
    
    for j = 1:40
        sig_correct_all(j).matrix = [sig_correct_all(j).matrix; sig_correct(j).matrix];
        sig_error_all(j).matrix = [sig_error_all(j).matrix; sig_error(j).matrix];
        sig_choice_all(j).matrix = [sig_choice_all(j).matrix; sig_choice(j).matrix];
        sig_prior_all(j).matrix = [sig_prior_all(j).matrix; sig_prior(j).matrix];
        short_correct_all(j).matrix = [short_correct_all(j).matrix; sig_correct_short(j).matrix];
        
        %ROC
        ROC_choice_all(j).matrix = [ROC_choice_all(j).matrix; ROC_choice(j).matrix];
        ROC_prior_all(j).matrix = [ROC_prior_all(j).matrix; ROC_prior(j).matrix];
        for k = 1:6
            evi_prefer_all(j,k).matrix = [evi_prefer_all(j,k).matrix; evi_prefer(j,k).matrix];
            evi_nonprefer_all(j,k).matrix = [evi_nonprefer_all(j,k).matrix; evi_nonprefer(j,k).matrix];
        end
    end
end

%Based on the evi_prefer and evi_non_prefer, determine which
%neurons to use
for j = 1:25
    clear nan_check1 nan_check2
    for k = 1:6
        moto_data = evi_prefer_all(j,k).matrix;
        nan_check1(:,k) = isnan(moto_data); %detect_non_nan
        moto_data = evi_nonprefer_all(j,k).matrix;
        nan_check2(:,k) = isnan(moto_data); %detect_non_nan
    end
    nan_check = [nan_check1,nan_check2];
    nan_check = max(nan_check,[],2);
    use_neuron = find(nan_check == 0);
    temp_check(j) = length(find(nan_check == 1));
    
    sig_prior_all(j).matrix = sig_prior_all(j).matrix(use_neuron);
    sig_choice_all(j).matrix = sig_choice_all(j).matrix(use_neuron);
    sig_correct_all(j).matrix = sig_correct_all(j).matrix(use_neuron);
    sig_error_all(j).matrix = sig_error_all(j).matrix(use_neuron);
    short_correct_all(j).matrix = short_correct_all(j).matrix(use_neuron);
    
    ROC_choice_all(j).matrix = ROC_choice_all(j).matrix(use_neuron);
    ROC_prior_all(j).matrix = ROC_prior_all(j).matrix(use_neuron);
    
    if j == 17
        sig_correct25_all17 = sig_correct25_all17(use_neuron);
    end
end    

%Sound on: 16
%Sound end: 25

%%% Fig S5a %%%
figure
subplot(1,3,1); hold on
plot_correct_error(short_correct_all,sig_correct_all,17);

subplot(1,3,2); hold on
plot(short_correct_all(17).matrix,sig_correct25_all17,'b.')
plot([-1 1],[0 0],'k')
plot([0 0],[-1 1],'k')
set(gca,'xlim',[-1 1],'ylim',[-1 1])

%Compare between the 17 and 25
sabun_sound_end = abs(short_correct_all(17).matrix-sig_correct_all(17).matrix);
sabun_sound_initial = abs(short_correct_all(17).matrix-sig_correct25_all17);

subplot(1,3,3); hold on
plot(sabun_sound_initial,sabun_sound_end,'k.');
plot([0 1],[0 1],'k')

% disp('short and long: correlation difference and ABS significance')
% corr(short_correct_all(17).matrix, sig_correct_all(17).matrix)
% corr(short_correct_all(17).matrix, sig_correct25_all17)
% signrank(sabun_sound_end,sabun_sound_initial)

%%% Fig S5b left %Before sound, compare the prior and choice index %%%
figure
abs_plot_prior_choice(sig_prior_all,sig_choice_all,15);

%Compare between the initial and end correct tone index
% correct17 = abs(sig_correct_all(17).matrix);
% correct25 = abs(sig_correct_all(25).matrix);
% correct17_25 = abs(sig_correct25_all17);
% ranksum(correct17, correct25)
% signrank(correct17, correct17_25)

%%% Fig S5b right %%%
figure; hold on
plot_ROC(ROC_prior_all,ROC_choice_all,15);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_ROC(sig_sound_all,sig_choice_all,use_frame)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(sig_sound_all(use_frame).matrix,sig_choice_all(use_frame).matrix,'.','color',[0.4 0.4 0.4])
plot([0.48 1],[0.48 1],'k')
set(gca,'xlim',[0.48 1],'ylim',[0.48 1])

% signrank(sig_sound_all(use_frame).matrix,sig_choice_all(use_frame).matrix)
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function abs_plot_prior_choice(sig_error_all,sig_correct_all,use_frame)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp1 = abs(sig_error_all(use_frame).matrix);
temp2 = abs(sig_correct_all(use_frame).matrix);

plot(temp1,temp2,'b.')
hold on
plot([-0.1 1],[-0.1 1],'k')
set(gca,'xlim',[-0.1 1],'ylim',[-0.1 1])
signrank(temp1,temp2)
length(temp1)
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_correct_error(sig_error_all,sig_correct_all,use_frame)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(sig_error_all(use_frame).matrix,sig_correct_all(use_frame).matrix,'b.')
plot([-1 1],[0 0],'k')
plot([0 0],[-1 1],'k')
set(gca,'xlim',[-1 1],'ylim',[-1 1])

temp_x = length(sig_error_all(use_frame).matrix);

[b,bint,r,rint,stats] = regress(sig_correct_all(use_frame).matrix,[sig_error_all(use_frame).matrix,ones(temp_x,1)]);
stats(3)
return