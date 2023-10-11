
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
function Figure6a_Population_decoder_20220520_SLR_process_one_session(pathname)

close all
switch nargin
    case 0
        pathname = pwd;
    case 1
    otherwise
        hoge
end
cd(pathname)

%selected_window = [10 11 12 15 17];
SLR_name = 'SLR_20220520_glmnet_sound_depth.mat';
SLR_choice = 'SLR_20220520_glmnet_choice_depth.mat';

thre_neuron = 20;
use_frame = [2, 4];

[~,~,~,prob_sound_all] = get_average_correct_rate_one_session(SLR_name, use_frame, 1, thre_neuron);

[~,~,~,prob_choice_all,~,~,~,~,~,evi_trial_all,Block_trial_all] = ...
    get_average_correct_rate_one_session(SLR_choice, use_frame, 1, thre_neuron);

cd('G:\upload_code\Figure6\Fig6a');
%% for sound decording / Fig6a top
sdata = struct();% source data 
[sdata.tonegroup,sdata.block,sdata.p_right]=...
    get_prob_SLR_one_session(prob_sound_all(:,2),evi_trial_all,Block_trial_all);

T = struct2table(sdata);
writetable(T, 'source fig 6a top.csv');

%% for choice decording / Fig6a middle
sdata = struct();% source data 
[sdata.tonegroup,sdata.block,sdata.p_right]=...
    get_prob_SLR_one_session(prob_choice_all(:,2),evi_trial_all,Block_trial_all);

T = struct2table(sdata);
writetable(T, 'source fig 6a middle.csv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tone,block,prob]=get_prob_SLR_one_session(prob_sound,evi_trial,Block_trial)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[evi_trial, Block_trial]
figure; hold on
use_color = [0 0 1; 1 0 0];
boxplot(prob_sound,evi_trial)

tone_evi=[0,0.25,0.45,0.55,0.75,1];
% block_label={'L','R'};

tone =nan(length(prob_sound),1);
block=nan(length(prob_sound),1);
prob =nan(length(prob_sound),1);
for i = 1:6
    temp_sound = find(evi_trial == i);
    
    for j = 0:1
        temp_block = find(Block_trial == j);
        temp = intersect(temp_sound, temp_block);
        temp_x = (rand(length(temp),1)-0.5) * 0.2 + i;
        plot(temp_x, prob_sound(temp,1),'.','color',use_color(j+1,:));
        
        tone(temp) =tone_evi(i);
        block(temp)=j;
        prob(temp) =prob_sound(temp);
    end
end
set(gca,'ylim',[0 1])

return
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [opt_l, opt_l_R, opt_l_L, ...
    prob_SLR_all, prob_SLR_R, prob_SLR_L, ...
    conf_all, conf_R, conf_L, evi_trial, trial_all] = ...
    get_average_correct_rate_one_session(SLR_name, use_frame, psycho_sign, thre_neuron)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(SLR_name);
analysis_dir = pwd;
%Check the number of neurons
temp = dir('sig_task_neurons_*'); %Bpod
if length(temp) ~= 1
    hoge
end
data = load(temp.name);
data = length(data.neuron_base);

if data < thre_neuron %Ready to analyze for all the definition
    hoge
end
if length(lasso_sound_l) ~= 1 %not nan
    [opt_l, opt_l_R, opt_l_L, detail_block_l] = ...
        Population_decoder_20220121_SLR_process2_long_detail(analysis_dir, SLR_name, use_frame, psycho_sign);
    
    %use likeilihood
    for i = 1:length(use_frame)
        Long_R = detail_block_l(i).matrix.Long_R;
        Long_L = detail_block_l(i).matrix.Long_L;
        Long_all = union(Long_L, Long_R);
        
        likelihood = detail_block_l(i).matrix.likelihood;
        evi_sound = detail_block_l(i).matrix.evi_trial;
        y_binary = detail_block_l(i).matrix.y_binary;
        
        if length(evi_sound) ~= 6
            hoge
        end
        evi_trial = nan(length(likelihood),1);
        for j = 1:6
            evi_trial(evi_sound(j).matrix) = j;
        end
        %nan check
        temp = isnan(evi_trial);
        if max(temp) == 1
            hoge
        end
        
        [prob_SLR_all(:,i), conf_all(i).matrix] = get_basic_SLR_results(likelihood(Long_all), evi_trial(Long_all), y_binary(Long_all));
        [prob_SLR_R(:,i), conf_R(i).matrix] = get_basic_SLR_results(likelihood(Long_R), evi_trial(Long_R), y_binary(Long_R));
        [prob_SLR_L(:,i), conf_L(i).matrix] = get_basic_SLR_results(likelihood(Long_L), evi_trial(Long_L), y_binary(Long_L));
        
        trial_all = nan(length(Long_all),1);
        for j = 1:length(Long_all)
            temp = find(Long_R == Long_all(j));
            if ~isempty(temp) %Right block
                trial_all(j) = 1;
            else %Left block
                trial_all(j) = 0;
            end
        end
    end
    evi_trial = evi_trial(Long_all);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [prob_SLR, conf_R] = get_basic_SLR_results(likelihood, evi_trial, y_test)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prob_SLR = nan(length(likelihood),1);
for i = 1:6
    temp = find(evi_trial == i);
    if i <= 3
        prob_SLR(temp) = 1-likelihood(temp);
    else
        prob_SLR(temp) = likelihood(temp);
    end
    temp_y = y_test(temp);
    
    %Get confidence interval
    [p_R(i,1),conf_R(i,:)] = binofit(sum(temp_y),length(temp));
end

return

