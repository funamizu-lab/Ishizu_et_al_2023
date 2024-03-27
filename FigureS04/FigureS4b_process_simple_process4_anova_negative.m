
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
function FigureS4b_process_simple_process4_anova_negative(folders)

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
   
    [~,~,~,~,~,~,~,~,~,evi_trace,evi_prefer,evi_nonprefer,~,~,~,~,~,~,~,~,~,sig_anova_long] = ...
        Task_kaiseki_tokyo1_20220516_sound_choice_process4A_nega(analysis_dir{i});
    
    for j = 1:40
        sig_anova_all(j).matrix = [sig_anova_all(j).matrix; sig_anova_long(j).matrix];
        for k = 1:6
            evi_trace_all(j,k).matrix = [evi_trace_all(j,k).matrix; evi_trace(j,k).matrix];
            evi_prefer_all(j,k).matrix = [evi_prefer_all(j,k).matrix; evi_prefer(j,k).matrix];
            evi_nonprefer_all(j,k).matrix = [evi_nonprefer_all(j,k).matrix; evi_nonprefer(j,k).matrix];
        end
    end
end


%Sound on: 16
%Sound end: 25
use_frame = 15:25;
disp_frame= [15,17,25];
save_frame= 16:25;
temp_x = [0 0.25 0.45 0.55 0.75 1];
for i=1:length(use_frame)
[prefer_activ(i).matrix,nonprefer_activ(i).matrix,~,~,sig_anova(i).matrix]=...
    plot_prior_block_trace2(use_frame(i),evi_prefer_all,evi_nonprefer_all,temp_x,sig_anova_all);
end
disp_id = find(ismember(use_frame,disp_frame));

neuron_prior_all = [];
neuron_choice_all = [];
neuron_number_all = [];

temp_x = [-inf, -0.095:0.01:0.095, inf];
temp_x_tick = -0.1:0.01:0.1;
for i = 1:length(use_frame)
    temp_p = prefer_activ(i).matrix;
    temp_n = nonprefer_activ(i).matrix;
    temp_anova = sig_anova(i).matrix;
    temp_anova_prior  = find(temp_anova(:,1) < 0.01);
    
    temp_prior = mean(temp_p,2) - mean(temp_n,2);
    neuron_prior(i).matrix = temp_prior;
    
    neuron_prior_all = [neuron_prior_all; neuron_prior(i).matrix];
    neuron_number_all = [neuron_number_all; ones(length(neuron_prior(i).matrix),1) * i];
    
    sig_prior2(i,:) = histcounts(temp_prior(temp_anova_prior),temp_x);
    all_prior(i,:) = histcounts(temp_prior,temp_x);    
    sig_neuron_number(i,1) = length(temp_anova_prior) ./ length(temp_prior);
end


max_prior = max(all_prior(:)); max_prior = ceil(max_prior/10)*10;
tag={'presound','soundstart','soundend'};
sdata1 = struct();% source data 
sdata2 = struct();% source data 
for i = 1:length(disp_id)
    t=disp_id(i);
    figure; hold on
    plot(temp_x_tick,all_prior(t,:),'k')
    eval(['sdata1.',tag{i},'_x=transpose(temp_x_tick)']);
    eval(['sdata1.',tag{i},'_y=transpose(all_prior(t,:))']);
    
    x=[temp_x_tick(1), temp_x_tick,temp_x_tick(end)];
    y=[0,sig_prior2(t,:),0];
    fill(x,y,[0.5 0.5 0.5],'EdgeColor',[0 0 0])
    eval(['sdata2.',tag{i},'_x=transpose(x)']);
    eval(['sdata2.',tag{i},'_y=transpose(y)']);
    
    set(gca,'xlim',[-0.105,0.105])
    ylim([0 max_prior]);
    title({['prior: ',num2str(use_frame(t))],[num2str(sum(all_prior(t,:))),' neuron'],[num2str(sig_neuron_number(i,1)), '*100%']});
end


%%% source data %%%
cd('G:\upload_code\FigureS4\FigS4b');
name = extractBefore(folders,'_');
T = struct2table(sdata1);
writetable(T, ['source fig ',name,' S4b black.csv']);
T = struct2table(sdata2);
writetable(T, ['source fig ',name,' S4b grey.csv']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [prefer_activ,nonprefer_activ,length_neuron,p_anovan,sig_anova] = plot_prior_block_trace2(use_frame,evi_prefer_all,evi_nonprefer_all,temp_x, sig_anova)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get the nan frame
for k = 1:6
    moto_data = evi_prefer_all(use_frame,k).matrix;
    nan_check1(:,k) = isnan(moto_data); %detect_non_nan
    moto_data = evi_nonprefer_all(use_frame,k).matrix;
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
    moto_data = evi_prefer_all(use_frame,k).matrix;
    moto_data = moto_data(use_neuron);
    median_prefer(k) = median(moto_data);
    std_prefer(k) = median(abs(moto_data-median_prefer(k)));
    se_prefer(k) = 1.4826 * std_prefer(k) ./ (sqrt(length(use_neuron)));

    prefer_activ(:,k) = moto_data;
    ano_activ = [ano_activ;moto_data];
    sound0 = ones(length(moto_data),1).*k;
    block0 = zeros(length(moto_data),1);
        
    moto_data = evi_nonprefer_all(use_frame,k).matrix;
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
sig_anova = sig_anova(use_frame).matrix;
sig_anova = sig_anova(use_neuron,:);

return