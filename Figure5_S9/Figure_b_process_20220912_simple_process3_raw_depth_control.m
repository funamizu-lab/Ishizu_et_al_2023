
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
function Figure_b_process_20220912_simple_process3_raw_depth_control(folders)

close all
analysis_dir = eval(folders);

sig_start_all = [];
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
    for j = 1:6
        evi_trace_all(i,j).matrix = [];
        evi_prefer_all(i,j).matrix = [];
        evi_nonprefer_all(i,j).matrix = [];
    end
end
for i = 1:length(analysis_dir)
   
    [prop_sig(i,:), number_sig(i,:), sig_correct, sig_error, sig_sound, sig_choice, sig_prior, sig_sin, ...
        choice_prior_matrix, sig_start, length_neuron(i,1), ...
        evi_trace, evi_prefer, evi_nonprefer, sig_trace_prior, sig_trace_nonprior,sig_correct_short] = ...
        Task_kaiseki_tokyo1_20220912_sound_choice_process3_raw_depth(analysis_dir{i});
    
    choice_prior_all = [choice_prior_all; choice_prior_matrix];
    sig_start_all = [sig_start_all; sig_start];
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
        
        for k = 1:6
            evi_trace_all(j,k).matrix = [evi_trace_all(j,k).matrix; evi_trace(j,k).matrix];
            evi_prefer_all(j,k).matrix = [evi_prefer_all(j,k).matrix; evi_prefer(j,k).matrix];
            evi_nonprefer_all(j,k).matrix = [evi_nonprefer_all(j,k).matrix; evi_nonprefer(j,k).matrix];
        end
    end
end
delete(gcp('nocreate'))

for j = 1:40
    mean_sig_prior(j) = nanmean(sig_trace_prior_all(j).matrix);
    mean_sig_nonprior(j) = nanmean(sig_trace_nonprior_all(j).matrix);
    for k = 1:6
        mean_evi_trace(k,j) = nanmean(evi_trace_all(j,k).matrix);
        mean_evi_prefer(k,j) = nanmean(evi_prefer_all(j,k).matrix);
        mean_evi_nonprefer(k,j) = nanmean(evi_nonprefer_all(j,k).matrix);
    end
end

%Sound on: 16
%Sound end: 25

%% fig b bottom panels
cd('G:\upload_code\Figure5_S9\b');
temp_prop_sig = mean(prop_sig,2);
temp_prop_sig = find(isnan(temp_prop_sig) == 0);
figure
sdata = struct();% source data 
subplot(2,1,1)
[mean_trace,~,se_trace]=plot_mean_se_moto(prop_sig(temp_prop_sig,:),[0 0 0],2);
set(gca,'ylim',[0 1])
set(gca,'xtick',0:2:40)
sdata.taskrelevant_mean = mean_trace';
sdata.taskrelevant_se = se_trace';

subplot(2,1,2)
[mean_trace,~,se_trace]=plot_mean_se_moto(number_sig(temp_prop_sig,:),[0 0 0],2);
set(gca,'ylim',[0 0.14],'ytick',0:0.02:0.14)
set(gca,'xtick',0:2:40)
sdata.deltaactivate_mean = mean_trace';
sdata.deltaactivate_se = se_trace';
T = struct2table(sdata);
writetable(T, 'source fig4b bottom.csv');
% writetable(T, 'source figS9b bottom.csv');


