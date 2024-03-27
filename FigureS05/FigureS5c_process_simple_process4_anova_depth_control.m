
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
function FigureS5c_process_simple_process4_anova_depth_control(folders)

close all
analysis_dir = eval(folders);

for i = 1:40
    sig_anova_all(i).matrix = [];
    for j = 1:6
        evi_trace_all(i,j).matrix = [];
        evi_prefer_all(i,j).matrix = [];
        evi_nonprefer_all(i,j).matrix = [];
    end
end
for i = 1:length(analysis_dir)
    [i,length(analysis_dir)]
   
    [~,~,~,~,~,~,~,~,~,~,~,evi_trace,evi_prefer,evi_nonprefer,~,~, ...
        ~,~,~,~,~,~,~, sig_anova_long] = ...
        Task_kaiseki_tokyo1_20220516_sound_choice_process4A_depth(analysis_dir{i});
    for j = 1:40
        sig_anova_all(j).matrix = [sig_anova_all(j).matrix; sig_anova_long(j).matrix];        
        for k = 1:6
            evi_trace_all(j,k).matrix = [evi_trace_all(j,k).matrix; evi_trace(j,k).matrix];
            evi_prefer_all(j,k).matrix = [evi_prefer_all(j,k).matrix; evi_prefer(j,k).matrix];
            evi_nonprefer_all(j,k).matrix = [evi_nonprefer_all(j,k).matrix; evi_nonprefer(j,k).matrix];
        end
    end
end

temp_x = [0 0.25 0.45 0.55 0.75 1];

%Sound on: 16
%Sound end: 25
prior_frame = 15;
short_frame = 17;
long_frame = 25;
[prefer_activ(1).matrix,nonprefer_activ(1).matrix,~,~,sig_anova(1).matrix] = plot_prior_block_trace2(prior_frame,evi_prefer_all,evi_nonprefer_all,temp_x,sig_anova_all);
[prefer_activ(2).matrix,nonprefer_activ(2).matrix,~,~,sig_anova(2).matrix] = plot_prior_block_trace2(short_frame,evi_prefer_all,evi_nonprefer_all,temp_x,sig_anova_all);
[prefer_activ(3).matrix,nonprefer_activ(3).matrix,~,~,sig_anova(3).matrix] = plot_prior_block_trace2(long_frame,evi_prefer_all,evi_nonprefer_all,temp_x,sig_anova_all);

neuron_prior_all = [];
neuron_number_all = [];

temp_x = [-inf, -0.095:0.01:0.095, inf];
temp_x_tick = -0.1:0.01:0.1;
size(temp_x)
size(temp_x_tick)
tag={'presound','soundstart','soundend'};
sdata1 = struct();% source data 
sdata2 = struct();% source data 
for i = 1:3
    temp_p = prefer_activ(i).matrix;
    temp_n = nonprefer_activ(i).matrix;
    temp_anova = sig_anova(i).matrix;
    temp_anova_prior  = find(temp_anova(:,1) < 0.01);
    
    temp_prior = mean(temp_p,2) - mean(temp_n,2);
    neuron_prior(i).matrix = temp_prior;
    
    neuron_prior_all = [neuron_prior_all; neuron_prior(i).matrix];
    neuron_number_all = [neuron_number_all; ones(length(neuron_prior(i).matrix),1) * i];
    
    sig_prior = histcounts(temp_prior(temp_anova_prior),temp_x);
    all_prior = histcounts(temp_prior,temp_x);
    
    %%% figure S5c %%%
    figure;    hold on
    plot(temp_x_tick,all_prior,'k')   
    eval(['sdata1.',tag{i},'_x=transpose(temp_x_tick)']);
    eval(['sdata1.',tag{i},'_y=transpose(all_prior)']);
    
    x = [temp_x_tick(1), temp_x_tick,temp_x_tick(end)];
    y = [0,sig_prior,0];
    fill(x,y,[0.5 0.5 0.5],'EdgeColor',[0 0 0])
    eval(['sdata2.',tag{i},'_x=transpose(x)']);
    eval(['sdata2.',tag{i},'_y=transpose(y)']);
    
    set(gca,'xlim',[-0.105,0.105])
        
    sig_neuron_number(i,1) = length(temp_anova_prior) ./ length(temp_prior);
end

%%% source data %%%
cd('G:\upload_code\FigureS5\FigS5c');
name = extractBefore(folders,'_');
T = struct2table(sdata1);
writetable(T, ['source fig ',name,' S5c rightpannel black.csv']);
T = struct2table(sdata2);
writetable(T, ['source fig ',name,' S5c rightpannel grey.csv']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [trace,length_neuron] = get_non_nan_trace(evi_trace_all,use_frame)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get the nan frame
for k = 1:6
    moto_data = evi_trace_all(use_frame,k).matrix;
    nan_check(:,k) = isnan(moto_data); %detect_non_nan
end
nan_check = max(nan_check,[],2);
use_neuron = find(nan_check == 0);
length_neuron = length(use_neuron);

for k = 1:6
    moto_data = evi_trace_all(use_frame,k).matrix;
    trace(k).matrix = moto_data(use_neuron,:);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [prefer_activ,nonprefer_activ,length_neuron,p_anovan,sig_anova] = plot_prior_block_trace2(prior_frame,evi_prefer_all,evi_nonprefer_all,temp_x, sig_anova)
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
    
sig_anova = sig_anova(prior_frame).matrix;
sig_anova = sig_anova(use_neuron,:);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [median_trace,std_trace,se_trace] = matrix2medians(moto_data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Extract nan trial
moto_data = moto_data(isnan(moto_data) == 0);

median_trace = median(moto_data);
std_trace  = median(abs(moto_data-median_trace));
se_trace = 1.4826 * std_trace ./ (sqrt(length(moto_data)));

return