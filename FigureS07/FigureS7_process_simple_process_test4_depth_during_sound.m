
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
function FigureS7_process_simple_process_test4_depth_during_sound(folders)

close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Determine the sig_neuron at certain time window
sig_time_window = 25; %Only during sound neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

analysis_dir = eval(folders);

check_overlap = [];
for i = 1:40
    for j = 1:6
        evi_trace_all(i,j).matrix = [];
        evi_prefer_all(i,j).matrix = [];
        evi_nonprefer_all(i,j).matrix = [];
    end
end
for i = 1:length(analysis_dir)
    [i,length(analysis_dir)]
   
    [~,~,~,~,~,~,~,~,~,~,~,evi_trace,evi_prefer,evi_nonprefer,~,~,~,~,~,~,~,~,~,check_neuron] = ...
        Task_kaiseki_tokyo1_20230707_sound_choice_process4_depth(analysis_dir{i},sig_time_window);
    
    check_overlap = [check_overlap; check_neuron];
    for j = 1:40    
        for k = 1:6
            evi_trace_all(j,k).matrix = [evi_trace_all(j,k).matrix; evi_trace(j,k).matrix];
            evi_prefer_all(j,k).matrix = [evi_prefer_all(j,k).matrix; evi_prefer(j,k).matrix];
            evi_nonprefer_all(j,k).matrix = [evi_nonprefer_all(j,k).matrix; evi_nonprefer(j,k).matrix];
        end
    end
end

use_frame = 6:40;
for j = 1:length(use_frame)
    temp_frame = use_frame(j);
    
    %Get the nan frame
    for k = 1:6
        moto_data = evi_prefer_all(temp_frame,k).matrix;
        nan_check_prefer(:,k) = isnan(moto_data); %detect_non_nan
        moto_data = evi_nonprefer_all(temp_frame,k).matrix;
        nan_check_nonprefer(:,k) = isnan(moto_data); %detect_non_nan
    end
    nan_check_prefer = max(nan_check_prefer,[],2);
    nan_check_nonprefer = max(nan_check_nonprefer,[],2);
    nan_check = max([nan_check_prefer,nan_check_nonprefer],[],2);
    use_neuron = find(nan_check == 0);
    
    evi_trace = get_non_nan_trace_new(evi_trace_all,temp_frame,use_neuron);
    for k = 1:6
        plot_evi_trace(k).matrix(j).matrix = evi_trace(k).matrix;
    end
end

temp_x = [0 0.25 0.45 0.55 0.75 1];

%%% FifS7 right 3 panels : neurometric plots %%%
cd('G:\upload_code\FigureS7');
name = extractBefore(folders,'_');
%Sound on: 16
%Sound end: 25
prior_frame = 15;
short_frame = 17;
long_frame = 25;
figure
subplot(1,3,1)
sdata = plot_prior_block_trace(prior_frame,evi_prefer_all,evi_nonprefer_all,temp_x);
%%% source data %%%
T = struct2table(sdata);
writetable(T, ['source fig ',name,' S7 rightbottompannel left.csv']);

subplot(1,3,2)
sdata = plot_prior_block_trace(short_frame,evi_prefer_all,evi_nonprefer_all,temp_x);
%%% source data %%%
T = struct2table(sdata);
writetable(T, ['source fig ',name,' S7 rightbottompannel mid.csv']);

subplot(1,3,3)
sdata = plot_prior_block_trace(long_frame,evi_prefer_all,evi_nonprefer_all,temp_x);
%%% source data %%%
T = struct2table(sdata);
writetable(T, ['source fig ',name,' S7 rightbottompannel right.csv']);


%%% FigS7 middle panel : 6tone %%%
use_color = jet(6);
figure; hold on

sdata = struct();% source data 
tag={'prefer100','prefer75','prefer55','prefer45','prefer25','prefer0'};
for k = 1:6
    [~,mean_trace,~,se_trace] = plot_mean_se_moto_x_axis_matrix(plot_evi_trace(k).matrix,[1:length(use_frame)],use_color(k,:),2);
   if(k==1)
        sdata.x_time=transpose((1:length(mean_trace))/10)-1;
    end
    eval(['sdata.',tag{end-k+1},'_mean=transpose(mean_trace)']);
    eval(['sdata.',tag{end-k+1},'_se=transpose(se_trace)']);
end
set(gca,'xlim',[0 length(use_frame)])
set(gca,'xtick',0:2:length(use_frame))

%%% source data %%%
T = struct2table(sdata);
writetable(T, ['source fig ',name,' S7 midlle bottom.csv']);


%%% Fig S7 left : pie chart %%%
other_neuron = find(check_overlap == 0);
pre_neuron = find(check_overlap == 1);
post_only = find(check_overlap == 2);
both = find(check_overlap == 3);
figure
pie([length(pre_neuron),length(both),length(post_only)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [trace,length_neuron] = get_non_nan_trace_new(evi_trace_all,use_frame,use_neuron)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get the nan frame
length_neuron = length(use_neuron);

for k = 1:6
    moto_data = evi_trace_all(use_frame,k).matrix;
    trace(k).matrix = moto_data(use_neuron,:);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sdata, prefer_activ,nonprefer_activ,length_neuron,p_anovan] = plot_prior_block_trace(prior_frame,evi_prefer_all,evi_nonprefer_all,temp_x)
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
        
    moto_data = evi_nonprefer_all(prior_frame,k).matrix;
    moto_data = moto_data(use_neuron);
    median_nonprefer(k) = median(moto_data);
    std_nonprefer(k) = median(abs(moto_data-median_nonprefer(k)));
    se_nonprefer(k) = 1.4826 * std_nonprefer(k) ./ (sqrt(length(use_neuron)));
    
    nonprefer_activ(:,k) = moto_data;
    ano_activ = [ano_activ;moto_data];
    sound1 = ones(length(moto_data),1).*k;
    block1 = ones(length(moto_data),1);

    ano_sound = [ano_sound;sound0;sound1];
    ano_block = [ano_block;block0;block1];
end
length_neuron = length(use_neuron);
p_anovan = anovan(ano_activ,{ano_sound,ano_block},'display','off');
  
hold on
plot_tuning_curve(temp_x, median_prefer, se_prefer, [1 0 0])
plot_tuning_curve(temp_x, median_nonprefer, se_nonprefer, [0 0 1])
set(gca,'xlim',[-0.1 1.1])

sdata = struct();% source data 
sdata.x=temp_x';
sdata.prefer=median_prefer';
sdata.prefer_se=se_prefer';
sdata.nonprefer=median_nonprefer';
sdata.nonprefer_se=se_nonprefer';

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_tuning_curve(temp_x, median_evi_trace, se_evi_trace, use_color)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(temp_x,median_evi_trace,'color',use_color) %sound start
for i = 1:6
    plot([temp_x(i),temp_x(i)],[median_evi_trace(i)+se_evi_trace(i),median_evi_trace(i)-se_evi_trace(i)],'color',use_color)
end
return