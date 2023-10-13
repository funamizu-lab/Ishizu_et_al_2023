
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
function Figure_d_process_20230701_simple_process3_raw_depth_control(folders)

close all
analysis_dir = eval(folders);

sig_start_all = [];
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
    for j = 1:6
        evi_trace_all(i,j).matrix = [];
        evi_prefer_all(i,j).matrix = [];
        evi_nonprefer_all(i,j).matrix = [];
    end
end
for i = 1:length(analysis_dir)
    [i,length(analysis_dir)]
   
    %Activity
    [~, ~, sig_correct, sig_error, sig_sound, sig_choice, sig_prior, sig_sin, ...
        choice_prior_matrix, sig_start,~, ...
        evi_trace, evi_prefer, evi_nonprefer, sig_trace_prior, sig_trace_nonprior, ...
        sig_correct_short, sig_correct25_from17] = ...
        Task_kaiseki_tokyo1_20221122_sound_choice_process3_raw_depth(analysis_dir{i});
    
    %ROC
    [~, ~, ROC_sound, ROC_choice, ROC_prior] = ...
        Task_kaiseki_tokyo1_20230701_sound_choice_ROC(analysis_dir{i});

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
        
        for k = 1:6
            evi_trace_all(j,k).matrix = [evi_trace_all(j,k).matrix; evi_trace(j,k).matrix];
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
    
    sig_prior_all(j).matrix = sig_prior_all(j).matrix(use_neuron);
    sig_choice_all(j).matrix = sig_choice_all(j).matrix(use_neuron);
    sig_correct_all(j).matrix = sig_correct_all(j).matrix(use_neuron);
    sig_error_all(j).matrix = sig_error_all(j).matrix(use_neuron);
    short_correct_all(j).matrix = short_correct_all(j).matrix(use_neuron);
    
    ROC_sound_all(j).matrix = ROC_sound_all(j).matrix(use_neuron);
    ROC_choice_all(j).matrix = ROC_choice_all(j).matrix(use_neuron);
    ROC_prior_all(j).matrix = ROC_prior_all(j).matrix(use_neuron);
    
    if j == 17
        sig_correct25_all17 = sig_correct25_all17(use_neuron);
    end
end    

%Check the number of neurons
if length(sig_correct_all(17).matrix) ~= length(ROC_choice_all(17).matrix)
    hoge
end
if length(sig_correct_all(25).matrix) ~= length(ROC_choice_all(25).matrix)
    hoge
end

%% Fig d
name = extractBefore(folders,'_');
region_num={'3','4','S8'};
tag = cell2mat(region_num(ismember({'mpfc','auc','fof'},name)));
cd('G:\upload_code\Figure3_4_S8\c');
tmp = dir(['fig',tag,'*.mat']);
load(tmp.name);

cd('G:\upload_code\Figure3_4_S8\d');
figure
subplot(2,2,1)
sdata = struct();% source data 
plot_correct_error(sig_error_all,sig_correct_all,17);
recval = [sig_error_all(17).matrix,sig_correct_all(17).matrix];
if(p_neuron(17)<1e-10)
    repval = [error_neuron(17),correct_neuron(17)];
else
    repval = [NaN, NaN];
end
recval(ismember(recval,repval,'rows'),:)=[];
recval = [repval;recval];
sdata.toneindex_error  = recval(:,1);
sdata.toneindex_correct= recval(:,2);
T = struct2table(sdata);
writetable(T, ['source fig',tag,'d left top.csv']);

subplot(2,2,2)
sdata = struct();% source data 
plot_correct_error(sig_error_all,sig_correct_all,25);
recval = [sig_error_all(25).matrix,sig_correct_all(25).matrix];
if(p_neuron(25)<1e-10)
    repval = [error_neuron(25),correct_neuron(25)];
else
    repval = [NaN, NaN];
end
recval(ismember(recval,repval,'rows'),:)=[];
recval = [repval;recval];
sdata.toneindex_error  = recval(:,1);
sdata.toneindex_correct= recval(:,2);
T = struct2table(sdata);
writetable(T, ['source fig',tag,'d right top.csv']);

subplot(2,2,3)
sdata = struct();% source data 
plot_ROC(ROC_sound_all,ROC_choice_all,17);
recval = [ROC_sound_all(17).matrix,ROC_choice_all(17).matrix];
if(p_neuron(17)<1e-10)
    repval = [sROC_neuron(17),cROC_neuron(17)];
else
    repval = [NaN, NaN];
end
recval(ismember(recval,repval,'rows'),:)=[];
recval = [repval;recval];
sdata.ROC_sound = recval(:,1);
sdata.ROC_choice= recval(:,2);
T = struct2table(sdata);
writetable(T, ['source fig',tag,'d left bottom.csv']);

subplot(2,2,4)
sdata = struct();% source data 
plot_ROC(ROC_sound_all,ROC_choice_all,25);
recval = [ROC_sound_all(25).matrix,ROC_choice_all(25).matrix];
if(p_neuron(25)<1e-10)
    repval = [sROC_neuron(25),cROC_neuron(25)];
else
    repval = [NaN, NaN];
end
recval(ismember(recval,repval,'rows'),:)=[];
recval = [repval;recval];
sdata.ROC_sound  = recval(:,1);
sdata.ROC_choice= recval(:,2);
T = struct2table(sdata);
writetable(T, ['source fig',tag,'d right bottom.csv']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_ROC(sig_sound_all,sig_choice_all,use_frame)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
plot(sig_sound_all(use_frame).matrix,sig_choice_all(use_frame).matrix,'.','color',[0.4 0.4 0.4])
plot([0.48 1],[0.48 1],'k')
set(gca,'xlim',[0.48 1],'ylim',[0.48 1])

% signrank(sig_sound_all(use_frame).matrix,sig_choice_all(use_frame).matrix)
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_correct_error(sig_error_all,sig_correct_all,use_frame)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
plot(sig_error_all(use_frame).matrix,sig_correct_all(use_frame).matrix,'b.')
plot([-1 1],[0 0],'k')
plot([0 0],[-1 1],'k')
set(gca,'xlim',[-1 1],'ylim',[-1 1])

% temp_x = length(sig_error_all(use_frame).matrix);
% [b,bint,r,rint,stats] = regress(sig_correct_all(use_frame).matrix,[sig_error_all(use_frame).matrix,ones(temp_x,1)]);
% stats(3)
return
