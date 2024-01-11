
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
function Figure_g_process_20220516_simple_process_test4_depth_control(folders)

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
    for j = 1:6
        evi_trace_all(i,j).matrix = [];
        evi_prefer_all(i,j).matrix = [];
        evi_nonprefer_all(i,j).matrix = [];
    end
end
for i = 1:32
    short_correct_all(i).matrix = [];
    s_trace_prior_all(i).matrix = [];
    s_trace_nonprior_all(i).matrix = [];
    for j = 1:6
        s_evi_trace_all(i,j).matrix = [];
        s_evi_prefer_all(i,j).matrix = [];
        s_evi_nonprefer_all(i,j).matrix = [];
    end
end
for i = 1:length(analysis_dir)
    [i,length(analysis_dir)]
   
    [prop_sig(i,:), number_sig(i,:), sig_correct, sig_error, sig_sound, sig_choice, sig_prior, sig_sin, ...
        choice_prior_matrix, sig_start, length_neuron(i), ...
        evi_trace, evi_prefer, evi_nonprefer, sig_trace_prior, sig_trace_nonprior, ...
        sig_correct_short, sig_sin_session(i,:), ...
        s_sig_trace_prior,s_sig_trace_nonprior,...
        s_evi_trace,s_evi_prefer,s_evi_nonprefer] = ...
        Task_kaiseki_tokyo1_20220516_sound_choice_process4_depth(analysis_dir{i});
    
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
        
        for k = 1:6
            evi_trace_all(j,k).matrix = [evi_trace_all(j,k).matrix; evi_trace(j,k).matrix];
            evi_prefer_all(j,k).matrix = [evi_prefer_all(j,k).matrix; evi_prefer(j,k).matrix];
            evi_nonprefer_all(j,k).matrix = [evi_nonprefer_all(j,k).matrix; evi_nonprefer(j,k).matrix];
        end
    end
    
    %for short sound
    for j = 1:32
        short_correct_all(j).matrix = [short_correct_all(j).matrix; sig_correct_short(j).matrix];
        s_trace_prior_all(j).matrix = [s_trace_prior_all(j).matrix; s_sig_trace_prior(j).matrix];
        s_trace_nonprior_all(j).matrix = [s_trace_nonprior_all(j).matrix; s_sig_trace_nonprior(j).matrix];
        for k = 1:6
            s_evi_trace_all(j,k).matrix = [s_evi_trace_all(j,k).matrix; s_evi_trace(j,k).matrix];
            s_evi_prefer_all(j,k).matrix = [s_evi_prefer_all(j,k).matrix; s_evi_prefer(j,k).matrix];
            s_evi_nonprefer_all(j,k).matrix = [s_evi_nonprefer_all(j,k).matrix; s_evi_nonprefer(j,k).matrix];
        end
    end
end



%% fig g
cd('G:\upload_code\Figure3_4_S8\g');
region_num={'3','4','S8'};
name = cell2mat(region_num(ismember({'mpfc','auc','fof'},extractBefore(folders,'_'))));

%Sound on: 16
%Sound end: 25
prior_frame = 15;
short_frame = 17;
long_frame = 25;

figure
sdata = struct();% source data 
temp_x = [0 0.25 0.45 0.55 0.75 1];
sdata.x=temp_x';
subplot(1,3,1)
[median_prefer,se_prefer,median_nonprefer,se_nonprefer]=...
    plot_prior_block_trace(prior_frame,evi_prefer_all,evi_nonprefer_all,temp_x);
sdata.leftpanel_prefer_median=median_prefer';
sdata.leftpanel_prefer_se =se_prefer';
sdata.leftpanel_nonprefer_median=median_nonprefer';
sdata.leftpanel_nonprefer_se =se_nonprefer';

subplot(1,3,2)
[median_prefer,se_prefer,median_nonprefer,se_nonprefer]=...
    plot_prior_block_trace(short_frame,evi_prefer_all,evi_nonprefer_all,temp_x);
sdata.middlepanel_prefer_median=median_prefer';
sdata.middlepanel_prefer_se =se_prefer';
sdata.middlepanel_nonprefer_median=median_nonprefer';
sdata.middlepanel_nonprefer_se =se_nonprefer';

subplot(1,3,3)
[median_prefer,se_prefer,median_nonprefer,se_nonprefer]=...
    plot_prior_block_trace(long_frame,evi_prefer_all,evi_nonprefer_all,temp_x);
sdata.rightpanel_prefer_median=median_prefer';
sdata.rightpanel_prefer_se =se_prefer';
sdata.rightpanel_nonprefer_median=median_nonprefer';
sdata.rightpanel_nonprefer_se =se_nonprefer';
T = struct2table(sdata);
writetable(T, ['source fig',name,'g.csv']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [median_prefer,se_prefer,median_nonprefer,se_nonprefer]=...
    plot_prior_block_trace(prior_frame,evi_prefer_all,evi_nonprefer_all,temp_x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get the nan frame
for k = 1:6
    moto_data = evi_prefer_all(prior_frame,k).matrix;
    nan_check1(:,k) = isnan(moto_data); %detect_non_nan
    moto_data = evi_nonprefer_all(prior_frame,k).matrix;
    nan_check2(:,k) = isnan(moto_data); %detect_non_nan
end
nan_check = [nan_check1,nan_check2];
nan_check = max(nan_check,[],2);
use_neuron = find(nan_check == 0);

%for anavan
ano_activ = [];
ano_sound = [];
ano_block = [];
for k = 1:6
    moto_data = evi_prefer_all(prior_frame,k).matrix;
    moto_data = moto_data(use_neuron);
    median_prefer(k) = median(moto_data);
    std_prefer(k) = median(abs(moto_data-median_prefer(k)));
    se_prefer(k) = 1.4826 * std_prefer(k) ./ (sqrt(length(use_neuron)));

    prefer_activ(:,k) = moto_data;
    ano_activ = [ano_activ;moto_data];
    sound0 = ones(length(moto_data),1).*k;
    block0 = zeros(length(moto_data),1);
        
    moto_data2 = evi_nonprefer_all(prior_frame,k).matrix;
    moto_data2 = moto_data2(use_neuron);
    median_nonprefer(k) = median(moto_data2);
    std_nonprefer(k) = median(abs(moto_data2-median_nonprefer(k)));
    se_nonprefer(k) = 1.4826 * std_nonprefer(k) ./ (sqrt(length(use_neuron)));
    
    nonprefer_activ(:,k) = moto_data2;
    ano_activ = [ano_activ;moto_data2];
    sound1 = ones(length(moto_data2),1).*k;
    block1 = ones(length(moto_data2),1);

    ano_sound = [ano_sound;sound0;sound1];
    ano_block = [ano_block;block0;block1];
    
    p_each_tone(k) = signrank(moto_data,moto_data2);
end
length_neuron = length(use_neuron);
p_anovan = anovan(ano_activ,{ano_sound,ano_block},'display','off');
    
hold on
plot_tuning_curve(temp_x, median_prefer, se_prefer, [1 0 0])
plot_tuning_curve(temp_x, median_nonprefer, se_nonprefer, [0 0 1])
set(gca,'xlim',[-0.1 1.1])

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_tuning_curve(temp_x, median_evi_trace, se_evi_trace, use_color)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(temp_x,median_evi_trace,'color',use_color) %sound start
for i = 1:6
    plot([temp_x(i),temp_x(i)],[median_evi_trace(i)+se_evi_trace(i),median_evi_trace(i)-se_evi_trace(i)],'color',use_color)
end
return
