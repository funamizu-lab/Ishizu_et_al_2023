
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

function Figure5cd_process_20231002_encoding_model_compare_depth_control

close all
[auc_sound,auc_choice,auc_prior, auc_prob, auc_neuron] = process_20231002_encoding_model2_depth('auc_ishizu');
[fof_sound,fof_choice,fof_prior, fof_prob, fof_neuron] = process_20231002_encoding_model2_depth('fof_ishizu');
[mpfc_sound,mpfc_choice,mpfc_prior, mpfc_prob, mpfc_neuron] = process_20231002_encoding_model2_depth('mpfc_ishizu');

%% fig 5d
cd('G:\upload_code\Figure5\Fig5cd');
sdata = struct();% source data 
 [hist_x,hist_auc,hist_fof,hist_mpfc] = ranksum_areas(auc_sound,fof_sound,mpfc_sound);
sdata.hist_x= hist_x';
sdata.sound_ac=hist_auc'; 
sdata.sound_m2=hist_fof';
sdata.sound_mpfc=hist_mpfc'; 

[~,hist_auc,hist_fof,hist_mpfc] = ranksum_areas(auc_choice,fof_choice,mpfc_choice);
sdata.choice_ac=hist_auc'; 
sdata.choice_m2=hist_fof';
sdata.choice_mpfc=hist_mpfc'; 

[~,hist_auc,hist_fof,hist_mpfc] = ranksum_areas(auc_prior,fof_prior,mpfc_prior);
sdata.priorval_ac=hist_auc'; 
sdata.priorval_m2=hist_fof';
sdata.priorval_mpfc=hist_mpfc'; 
 
T = struct2table(sdata);
writetable(T, 'source fig 5d.csv');

%% fig 5c
figure
sdata = struct();% source data 
bar([auc_prob;fof_prob;mpfc_prob]')
sdata.sound_ac  =auc_prob(1); 
sdata.sound_m2  =fof_prob(1);
sdata.sound_mpfc=mpfc_prob(1); 
sdata.choice_ac  =auc_prob(2); 
sdata.choice_m2  =fof_prob(2);
sdata.choice_mpfc=mpfc_prob(2); 
sdata.priorval_ac=auc_prob(3); 
sdata.priorval_m2=fof_prob(3);
sdata.priorval_mpfc=mpfc_prob(3); 
T = struct2table(sdata);
writetable(T, 'source fig 5d.csv');

% Significant test with kai 2 test
% sound, choice, prior 
p_sound = kai2test_3regions(auc_neuron,fof_neuron,mpfc_neuron,1);
p_choice= kai2test_3regions(auc_neuron,fof_neuron,mpfc_neuron,2);
p_prior = kai2test_3regions(auc_neuron,fof_neuron,mpfc_neuron,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = kai2test_3regions(auc_neuron,fof_neuron,mpfc_neuron,use_number)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
use_auc =  [auc_neuron(use_number), auc_neuron(4)- auc_neuron(use_number)];
use_fof =  [fof_neuron(use_number), fof_neuron(4)- fof_neuron(use_number)];
use_mpfc = [mpfc_neuron(use_number),mpfc_neuron(4)-mpfc_neuron(use_number)];
p(1) = kai2_test(use_auc,use_fof);
p(2) = kai2_test(use_auc,use_mpfc);
p(3) = kai2_test(use_mpfc,use_fof);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hist_x,hist_auc,hist_fof,hist_mpfc] = ranksum_areas(auc_sound,fof_sound,mpfc_sound)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hist_x = [-inf,-600:10:-400,inf];
hist_auc = histcounts(auc_sound,hist_x);
hist_fof = histcounts(fof_sound,hist_x);
hist_mpfc = histcounts(mpfc_sound,hist_x);
hist_auc =hist_auc./sum(hist_auc);
hist_fof = hist_fof./sum(hist_fof);
hist_mpfc = hist_mpfc./sum(hist_mpfc);

figure;hold on
plot(hist_auc./sum(hist_auc),'r')
plot(hist_fof./sum(hist_fof),'g')
plot(hist_mpfc./sum(hist_mpfc),'b')
set(gca,'xlim',[0 length(hist_x)])
set(gca,'xtick',1:2:length(hist_x))

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [BIC_sound,BIC_choice,BIC_prior,prob_sig,neuron_sig] = process_20231002_encoding_model2_depth(folders)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

analysis_dir = eval(folders);
BIC_all = [];
LL_all = [];
y_sound_all = [];
y_choice_all = [];
y_prior_all = [];
for i = 1:30
    BIC_each(i).matrix = [];
    LL_each(i).matrix = [];
    
    y_sound_each(i).matrix = [];
    y_choice_each(i).matrix = [];
    y_prior_each(i).matrix = [];
end

for i = 1:length(analysis_dir)
    [i,length(analysis_dir)]
   
    [temp_BIC_each,temp_LL_each,temp_BIC_all,temp_LL_all, ...
        parameter_size,parameter_sound,parameter_choice,parameter_prior,...
        temp_sound_each,temp_choice_each,temp_prior_each, ...
        temp_sound_all,temp_choice_all,temp_prior_all] = ...
        Task_kaiseki_tokyo1_20231002_encoding_model_depth_control(analysis_dir{i});
    
    BIC_all = [BIC_all; temp_BIC_all];
    LL_all = [LL_all; temp_LL_all];
    y_sound_all = [y_sound_all; temp_sound_all];
    y_choice_all = [y_choice_all; temp_choice_all];
    y_prior_all = [y_prior_all; temp_prior_all];
    for j = 1:30
        BIC_each(j).matrix = [BIC_each(j).matrix; temp_BIC_each(j).matrix];
        LL_each(j).matrix = [LL_each(j).matrix; temp_LL_each(j).matrix];
        y_sound_each(j).matrix = [y_sound_each(j).matrix; temp_sound_each(j).matrix];
        y_choice_each(j).matrix = [y_choice_each(j).matrix; temp_choice_each(j).matrix];
        y_prior_each(j).matrix = [y_prior_each(j).matrix; temp_prior_each(j).matrix];        
    end
end
delete(gcp('nocreate'))

[BIC_sound,BIC_choice,BIC_prior,prob_sig,neuron_sig] = get_BIC_loglikeli_plot(BIC_all, LL_all,parameter_size,parameter_sound,parameter_choice,parameter_prior);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [BIC_sound,BIC_choice,BIC_prior,prob_sig,neuron_sig] = get_BIC_loglikeli_plot(BIC, part_LL,parameter_size,parameter_sound,parameter_choice,parameter_prior)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nan check
nanBIC = mean(BIC,2);
nanBIC = find(isnan(nanBIC) == 0);
BIC = BIC(nanBIC,:);
part_LL = part_LL(nanBIC,:);

BIC_sound = BIC(:,2) - BIC(:,1);
BIC_choice = BIC(:,3) - BIC(:,1);
BIC_prior = BIC(:,4) - BIC(:,1);

LL_sound = part_LL(:,1) - part_LL(:,2);
LL_choice = part_LL(:,1) - part_LL(:,3);
LL_prior = part_LL(:,1) - part_LL(:,4);

p_LL_sound = 1-chi2cdf(2*LL_sound, parameter_size-parameter_sound);
p_LL_sound = -log10(p_LL_sound);
p_LL_choice = 1-chi2cdf(2*LL_choice, parameter_size-parameter_choice);
p_LL_choice = -log10(p_LL_choice);
p_LL_prior = 1-chi2cdf(2*LL_prior, parameter_size-parameter_prior);
p_LL_prior = -log10(p_LL_prior);

sig_sound = find(p_LL_sound > 2);
sig_choice = find(p_LL_choice > 2);
sig_prior = find(p_LL_prior > 2);

temp_sig = [length(sig_sound),length(sig_choice),length(sig_prior)];
neuron_sig = [temp_sig,length(p_LL_sound)];
prob_sig = temp_sig ./ length(p_LL_sound);

return
