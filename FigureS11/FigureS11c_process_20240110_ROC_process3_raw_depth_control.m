
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
function FigureS11c_process_20240110_ROC_process3_raw_depth_control(folders)

analysis_dir = eval(folders);
analysis_dir

Sound_prefer_index = [];
sig_start_all = [];
p_all = [];
sig_correct25_all17 = [];
choice_prior_all = [];
for i = 1:40
    sig_correct_all(i).matrix = [];
    sig_error_all(i).matrix = [];
    sig_sound_all(i).matrix = [];
    sig_choice_all(i).matrix = [];
    sig_prior_all(i).matrix = [];
    sig_sin_all(i).matrix = [];
    sig_choice1_all(i).matrix = [];
    sig_choice2_all(i).matrix = [];
    sig_prior1_all(i).matrix = [];
    sig_prior2_all(i).matrix = [];
    sig_trace_prior_all(i).matrix = [];
    sig_trace_nonprior_all(i).matrix = [];
    short_correct_all(i).matrix = [];
    
    ROC_sound_all(i).matrix = [];
    ROC_choice_all(i).matrix = [];
    ROC_prior_all(i).matrix = [];

    ROC_sound_all_pre(i).matrix = [];
    ROC_choice_all_pre(i).matrix = [];
    ROC_prior_all_pre(i).matrix = [];
for j = 1:6
        evi_trace_all(i,j).matrix = [];
        evi_prefer_all(i,j).matrix = [];
        evi_nonprefer_all(i,j).matrix = [];
    end
end
for i = 1:length(analysis_dir)
    [i,length(analysis_dir)]
   
    %Activity
    [prop_sig(i,:), number_sig(i,:), sig_correct, sig_error, sig_sound, sig_choice, sig_prior, sig_sin, ...
        choice_prior_matrix, sig_start, length_neuron(i,1), ...
        evi_trace, evi_prefer, evi_nonprefer, sig_trace_prior, sig_trace_nonprior, ...
        sig_correct_short, sig_correct25_from17] = ...
        Task_kaiseki_tokyo1_20221122_sound_choice_process3_raw_depth(analysis_dir{i});
    
    %ROC
    [~, ~, ROC_sound, ROC_choice, ROC_prior, ...
        ROC_short_sound, ROC_short_choice, ROC_short_prior, ~, ...
        ROC_sound25_from17, ROC_choice25_from17] = ...
        Task_kaiseki_tokyo1_20230701_sound_choice_ROC(analysis_dir{i});

    [~, ~, ROC_sound_pre, ROC_choice_pre, ROC_prior_pre, ...
        ~, ~, ~, ~, ...
        ~, ~] = ...
        Task_kaiseki_tokyo1_20240110_sound_choice_ROC(analysis_dir{i});

    choice_prior_all = [choice_prior_all; choice_prior_matrix];
    sig_start_all = [sig_start_all; sig_start];
    sig_correct25_all17 = [sig_correct25_all17; sig_correct25_from17];
    
    for j = 1:40
        sig_correct_all(j).matrix = [sig_correct_all(j).matrix; sig_correct(j).matrix];
        sig_error_all(j).matrix = [sig_error_all(j).matrix; sig_error(j).matrix];
        sig_sound_all(j).matrix = [sig_sound_all(j).matrix; sig_sound(j).matrix];
        sig_choice_all(j).matrix = [sig_choice_all(j).matrix; sig_choice(j).matrix];
        sig_prior_all(j).matrix = [sig_prior_all(j).matrix; sig_prior(j).matrix];
        sig_sin_all(j).matrix = [sig_sin_all(j).matrix; sig_sin(j).matrix];
        
        sig_trace_prior_all(j).matrix = [sig_trace_prior_all(j).matrix; sig_trace_prior(j).matrix];
        sig_trace_nonprior_all(j).matrix = [sig_trace_nonprior_all(j).matrix; sig_trace_nonprior(j).matrix];
        short_correct_all(j).matrix = [short_correct_all(j).matrix; sig_correct_short(j).matrix];
        
        %ROC
        ROC_sound_all(j).matrix = [ROC_sound_all(j).matrix; ROC_sound(j).matrix];
        ROC_choice_all(j).matrix = [ROC_choice_all(j).matrix; ROC_choice(j).matrix];
        ROC_prior_all(j).matrix = [ROC_prior_all(j).matrix; ROC_prior(j).matrix];
        
        ROC_sound_all_pre(j).matrix = [ROC_sound_all_pre(j).matrix; ROC_sound_pre(j).matrix];
        ROC_choice_all_pre(j).matrix = [ROC_choice_all_pre(j).matrix; ROC_choice_pre(j).matrix];
        ROC_prior_all_pre(j).matrix = [ROC_prior_all_pre(j).matrix; ROC_prior_pre(j).matrix];

        for k = 1:6
            evi_trace_all(j,k).matrix = [evi_trace_all(j,k).matrix; evi_trace(j,k).matrix];
            evi_prefer_all(j,k).matrix = [evi_prefer_all(j,k).matrix; evi_prefer(j,k).matrix];
            evi_nonprefer_all(j,k).matrix = [evi_nonprefer_all(j,k).matrix; evi_nonprefer(j,k).matrix];
        end
    end
end
delete(gcp('nocreate'))

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
    
    ROC_sound_all(j).matrix = ROC_sound_all(j).matrix(use_neuron);
    ROC_choice_all(j).matrix = ROC_choice_all(j).matrix(use_neuron);
    ROC_prior_all(j).matrix = ROC_prior_all(j).matrix(use_neuron);

    ROC_sound_all_pre(j).matrix = ROC_sound_all_pre(j).matrix(use_neuron);
    ROC_choice_all_pre(j).matrix = ROC_choice_all_pre(j).matrix(use_neuron);
    ROC_prior_all_pre(j).matrix = ROC_prior_all_pre(j).matrix(use_neuron);

    if j == 17
        sig_correct25_all17 = sig_correct25_all17(use_neuron);
    end
end    
temp_check

% use_color = jet(6);
% temp_x = [0 0.25 0.45 0.55 0.75 1];
% 
% [~,sort_sig] = sort(sig_start_all);
% 
% %figure
% %imagesc(choice_prior_all(sort_sig,:))
% %Sound on: 16
% %Sound end: 25
% 
% temp_prop_sig = mean(prop_sig,2);
% temp_prop_sig = find(isnan(temp_prop_sig) == 0);
% length(temp_prop_sig)
% 
% %Check the number of neurons
% if length(sig_correct_all(17).matrix) ~= length(ROC_choice_all(17).matrix)
%     hoge
% end
% if length(sig_correct_all(25).matrix) ~= length(ROC_choice_all(25).matrix)
%     hoge
% end
% 
% disp('compare correct and error')
% [length(sig_correct_all(17).matrix), length(sig_correct_all(25).matrix)]
% figure
% subplot(2,2,1)
% plot_correct_error(sig_error_all,sig_correct_all,17);
% subplot(2,2,2)
% plot_correct_error(sig_error_all,sig_correct_all,25);
% subplot(2,2,3)
% plot_ROC(ROC_sound_all,ROC_choice_all,17);
% subplot(2,2,4)
% plot_ROC(ROC_sound_all,ROC_choice_all,25);
% 
% disp('compare short and long')
% figure
% subplot(1,3,1)
% plot_correct_error(short_correct_all,sig_correct_all,17);
% %plot_correct_error(sig_error_all,sig_correct_all,29);
% % signrank(sig_error_all(17).matrix)
% % signrank(sig_error_all(25).matrix)
% % signrank(short_correct_all(17).matrix, sig_correct_all(17).matrix)
% [length(sig_correct_all(17).matrix), length(sig_correct_all(25).matrix)]
% 
% subplot(1,3,2)
% plot(short_correct_all(17).matrix,sig_correct25_all17,'b.')
% hold on
% plot([-1 1],[0 0],'k')
% hold on
% %plot([0 0],[-0.1 1],'k')
% plot([0 0],[-1 1],'k')
% set(gca,'xlim',[-1 1],'ylim',[-1 1])
% 
% %Compare between the 17 and 25
% sabun_sound_end = abs(short_correct_all(17).matrix-sig_correct_all(17).matrix);
% sabun_sound_initial = abs(short_correct_all(17).matrix-sig_correct25_all17);
% 
% subplot(1,3,3)
% %boxplot([sabun_sound_end,sabun_sound_initial]);
% plot(sabun_sound_initial,sabun_sound_end,'k.');
% hold on
% plot([0 1],[0 1],'k')
% 
% disp('short and long: correlation difference and ABS significance')
% corr(short_correct_all(17).matrix, sig_correct_all(17).matrix)
% corr(short_correct_all(17).matrix, sig_correct25_all17)
% signrank(sabun_sound_end,sabun_sound_initial)

% %Compare between the initial and end correct tone index
% correct17 = abs(sig_correct_all(17).matrix);
% correct25 = abs(sig_correct_all(25).matrix);
% correct17_25 = abs(sig_correct25_all17);
% 
% ranksum(correct17, correct25)
% signrank(correct17, correct17_25)
% [median(correct17),median(correct25),median(correct17_25)]

disp('compare choice and prior before sound')
figure
%Before sound, compare the prior and choice index
subplot(1,2,1)
abs_plot_prior_choice(sig_prior_all,sig_choice_all,15);
subplot(1,2,2)
plot_ROC(ROC_prior_all,ROC_choice_all,15);

%%% Fig S11c %%%
cd('G:\upload_code\FigureS11');
figure
sdata = struct(); % source data
subplot(1,2,1)
[xdata,ydata] = plot_ROC(ROC_prior_all,ROC_sound_all_pre,15);
sdata.x = xdata;
sdata.y = ydata;
%%% source data %%%
T = struct2table(sdata);
writetable(T, 'source fig S11c.csv');

subplot(1,2,2)
plot_ROC(ROC_prior_all,ROC_choice_all_pre,15);

hoge

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
function plot_prior_choice(sig_error_all,sig_correct_all,use_frame)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(sig_error_all(use_frame).matrix,sig_correct_all(use_frame).matrix,'b.')
hold on
plot([-1 1],[0 0],'k')
hold on
plot([0 0],[-1 1],'k')
hold on
plot([-1 1],[-1 1],'k')
set(gca,'xlim',[-1 1],'ylim',[-1 1])

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xdata,ydata] = plot_ROC(sig_sound_all,sig_choice_all,use_frame)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xdata = sig_sound_all(use_frame).matrix;
ydata = sig_choice_all(use_frame).matrix;
plot(xdata,ydata,'.','color',[0.4 0.4 0.4])
hold on
%plot([0.5 1],[0.5 1],'k')
plot([0.48 1],[0.48 1],'k')
% hold on
% plot([0 0],[-1 1],'k')
set(gca,'xlim',[0.48 1],'ylim',[0.48 1])

temp_x = length(sig_sound_all(use_frame).matrix);

[median(sig_sound_all(use_frame).matrix), median(sig_choice_all(use_frame).matrix)]
signrank(sig_sound_all(use_frame).matrix,sig_choice_all(use_frame).matrix)
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_correct_error(sig_error_all,sig_correct_all,use_frame)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(sig_error_all(use_frame).matrix,sig_correct_all(use_frame).matrix,'b.')
hold on
plot([-1 1],[0 0],'k')
hold on
%plot([0 0],[-0.1 1],'k')
plot([0 0],[-1 1],'k')
%set(gca,'xlim',[-1 1],'ylim',[-0.1 1])
set(gca,'xlim',[-1 1],'ylim',[-1 1])

temp_x = length(sig_error_all(use_frame).matrix);

[b,bint,r,rint,stats] = regress(sig_correct_all(use_frame).matrix,[sig_error_all(use_frame).matrix,ones(temp_x,1)]);
b
stats(3)
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_basic_tone_choice_index(Sound_prefer_index, use_neuron)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
subplot(1,3,1) %Short sound
plot(Sound_prefer_index(use_neuron,4), Sound_prefer_index(use_neuron,1),'b.')
hold on
plot([0 0],[-1 1],'k--')
hold on
plot([-1 1],[0 0],'k--')
set(gca,'xlim',[-1 1],'ylim',[-1 1])
%signrank(Sound_prefer_index(:,4), Sound_prefer_index(:,1))
subplot(1,3,2) %Long sound
plot(Sound_prefer_index(use_neuron,5), Sound_prefer_index(use_neuron,2),'b.')
hold on
plot([0 0],[-1 1],'k--')
hold on
plot([-1 1],[0 0],'k--')
set(gca,'xlim',[-1 1],'ylim',[-1 1])
%signrank(Sound_prefer_index(:,5), Sound_prefer_index(:,2))
subplot(1,3,3) %Long sound end
plot(Sound_prefer_index(use_neuron,6), Sound_prefer_index(use_neuron,3),'b.')
hold on
plot([0 0],[-1 1],'k--')
hold on
plot([-1 1],[0 0],'k--')
set(gca,'xlim',[-1 1],'ylim',[-1 1])
%signrank(Sound_prefer_index(:,6), Sound_prefer_index(:,3))

signrank(Sound_prefer_index(use_neuron,4))
signrank(Sound_prefer_index(use_neuron,5))
signrank(Sound_prefer_index(use_neuron,6))

figure %Onset and Offset of sound @ correct trials
plot(Sound_prefer_index(use_neuron,3),Sound_prefer_index(use_neuron,2),'b.')
hold on
plot([-1 1],[-1 1],'k')
hold on
plot([0 0],[-1 1],'k--')
hold on
plot([-1 1],[0 0],'k--')
set(gca,'xlim',[-1 1],'ylim',[-1 1])
signrank(Sound_prefer_index(use_neuron,2), Sound_prefer_index(use_neuron,3))
%signrank(abs(Sound_prefer_index(use_neuron,2)), abs(Sound_prefer_index(use_neuron,3)))
[median(Sound_prefer_index(use_neuron,2)), median(Sound_prefer_index(use_neuron,3))]

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_sigR_sig_L(sig_all,sig_R,sig_L,tuning0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sig_all1=sig_all(1).matrix(tuning0,:);
sig_all2=sig_all(2).matrix(tuning0,:);
sig_all3=sig_all(3).matrix(tuning0,:);
sig_R1=sig_R(1).matrix(tuning0,:);
sig_R2=sig_R(2).matrix(tuning0,:);
sig_R3=sig_R(3).matrix(tuning0,:);
sig_L1=sig_L(1).matrix(tuning0,:);
sig_L2=sig_L(2).matrix(tuning0,:);
sig_L3=sig_L(3).matrix(tuning0,:);

figure
subplot(1,3,1)
plot(nanmean(sig_all1),'k')
hold on
plot(nanmean(sig_R1),'r')
hold on
plot(nanmean(sig_L1),'b')
subplot(1,3,2)
plot(nanmean(sig_all2),'k')
hold on
plot(nanmean(sig_R2),'r')
hold on
plot(nanmean(sig_L2),'b')
subplot(1,3,3)
plot(nanmean(sig_all3),'k')
hold on
plot(nanmean(sig_R3),'r')
hold on
plot(nanmean(sig_L3),'b')

% figure
% for i = 1:6
%     max_sig = max(max([sig_R1(:,i),sig_L1(:,i)]));
%     subplot(2,3,i)
%     plot(sig_R1(:,i),sig_L1(:,i),'b.')
%     hold on
%     plot([0 max_sig],[0 max_sig],'k')
% end
% figure
% for i = 1:6
%     max_sig = max(max([sig_R2(:,i),sig_L2(:,i)]));
%     subplot(2,3,i)
%     plot(sig_R2(:,i),sig_L2(:,i),'b.')
%     hold on
%     plot([0 max_sig],[0 max_sig],'k')
% end
% figure
% for i = 1:6
%     max_sig = max(max([sig_R3(:,i),sig_L3(:,i)]));
%     subplot(2,3,i)
%     plot(sig_R3(:,i),sig_L3(:,i),'b.')
%     hold on
%     plot([0 max_sig],[0 max_sig],'k')
% end

clear p
for i = 1:6
    p(1,i) = signrank(sig_R1(:,i),sig_L1(:,i));
    p(2,i) = signrank(sig_R2(:,i),sig_L2(:,i));
    p(3,i) = signrank(sig_R3(:,i),sig_L3(:,i));
end
p