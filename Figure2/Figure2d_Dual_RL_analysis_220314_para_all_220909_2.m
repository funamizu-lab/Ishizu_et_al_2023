%{
----------------------------------------------------------------------------
Analyzing behavioral data
At least for the correct rate
----------------------------------------------------------------------------
%}

function Figure2d_Dual_RL_analysis_220314_para_all_220909_2
close all
analysis_folder = {
    'G:\Ishizu_data\IntermediateFiles\only_all_behaviors\a04_behave\ML_CumGauss_Fit\RL_20220818'
    'G:\Ishizu_data\IntermediateFiles\only_all_behaviors\a08_behave\ML_CumGauss_Fit\RL_20220818'
    'G:\Ishizu_data\IntermediateFiles\only_all_behaviors\i20_behave\ML_CumGauss_Fit\RL_20220818'
    'G:\Ishizu_data\IntermediateFiles\only_all_behaviors\i24_behave\ML_CumGauss_Fit\RL_20220818'
    'G:\Ishizu_data\IntermediateFiles\only_all_behaviors\i34_behave\ML_CumGauss_Fit\RL_20220818'
    'G:\Ishizu_data\IntermediateFiles\only_all_behaviors\i35_behave\ML_CumGauss_Fit\RL_20220818'
    'G:\Ishizu_data\IntermediateFiles\only_all_behaviors\i43_behave\ML_CumGauss_Fit\RL_20220818'
    'G:\Ishizu_data\IntermediateFiles\only_all_behaviors\i46_behave\ML_CumGauss_Fit\RL_20220818'
    };

AIC_para = [4, 5, 6, 6, 7];
use_para = 3;
[~,~,~,~,para_max_all,N_para,all_subject] = ...
    get_likelihood_value(analysis_folder,use_para,AIC_para);

%% figure 2d
%Standard deviation of Gaussian
compare_parameters_RL([4,3],para_max_all,N_para,all_subject); %short long
%Bias for short and long
% compare_parameters_RL([6,5],para_max_all,N_para,all_subject); %short long

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function compare_parameters_RL(use_para,para_max_all,N_para,all_subject)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lme_para_sabun = fitlme_analysis_20210520_0(para_max_all(:,use_para(1))-para_max_all(:,use_para(2)),all_subject);
lme_para_sabun(1).lme
lme_para_sabun(2).lme

signrank(N_para(:,use_para(1)),N_para(:,use_para(2)))

cd('G:\upload_code\Figure2\Fig2d');
figure; hold on
sdata = struct();% source data 
plot(para_max_all(:,use_para)')
boxplot(para_max_all(:,use_para))
tmp=para_max_all(:,use_para);
sdata.short = tmp(:,1);
sdata.long = tmp(:,2);
T = struct2table(sdata);
writetable(T, 'source fig2d.csv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [BIC_all,ave_likeli_all,log_likeli_all,...
    sub_para,para_max_all,N_para,all_subject,N_loglikeli,N_BIC,N_AIC] = ...
    get_likelihood_value(analysis_folder,use_para,AIC_para)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

count = 0;
para_max_all = [];
for i = 1:length(analysis_folder)
    temp_folder = analysis_folder{i};
    cd(temp_folder);
    filename1 = dir('*.mat');
    Ndata(i,1) = length(filename1);
    
    clear sub_para sub_likeli sub_BIC sub_AIC
    for filecount = 1 : length(filename1)
        temp_filename = filename1(filecount).name;
        temp_file{filecount} = temp_filename;
        fpath = temp_filename;
    
        load(fpath)
    %'ave_likeli' ,'BIC','log_likeli','para_max','N_trial'
        %Calculate AIC
        count = count + 1;
        AIC = -2 * log_likeli + 2 * AIC_para;
        AIC_all(count,:) = AIC;
        BIC_all(count,:) = BIC;
        ave_likeli_all(count,:) = ave_likeli;
        log_likeli_all(count,:) = log_likeli;
        sub_para(filecount,:) = para_max(use_para,:);
        sub_likeli(filecount,:) = log_likeli;
        sub_BIC(filecount,:) = BIC;
        sub_AIC(filecount,:) = AIC;
    end
    para_max_all = [para_max_all; sub_para];
    N_para(i,:) = mean(sub_para);
    N_loglikeli(i,:) = mean(sub_likeli);
    N_BIC(i,:) = mean(sub_BIC);
    N_AIC(i,:) = mean(sub_AIC);
end

all_subject = [];
for i = 1:length(analysis_folder)
    temp = ones(Ndata(i),1) * i;
    all_subject = [all_subject; temp];
end

return