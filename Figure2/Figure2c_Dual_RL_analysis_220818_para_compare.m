%{
----------------------------------------------------------------------------
Analyzing behavioral data
At least for the correct rate
----------------------------------------------------------------------------
%}

function Figure2c_Dual_RL_analysis_220818_para_compare
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

use_para = 3;
count = 0;
BIC_all = [];
Psy_BIC_all = [];
for i = 1:length(analysis_folder)
    temp_folder = analysis_folder{i};
    cd(temp_folder);
    filename1 = dir('*.mat');
    Ndata(i,1) = length(filename1);
    
    clear sub_BIC
    for filecount = 1 : length(filename1)
        temp_filename = filename1(filecount).name;
        temp_file{filecount} = temp_filename;
        fpath = temp_filename;
    
        load(fpath); %'ave_likeli' ,'BIC','log_likeli','para_max','N_trial'

        count = count + 1;
        sub_BIC(filecount,:) = BIC;
        ave_likeli_all(count,:) = ave_likeli;
        log_likeli_all(count,:) = log_likeli;
        para_max_all(count,:) = para_max(use_para,:);
    end
    BIC_all = [BIC_all; sub_BIC];
    N_BIC(i,:) = mean(sub_BIC);
    
    %Get the BIC of psychometric function
    cd ../
    data = load('ML_CumGauss_Fit_20220114.mat');
    temp_BIC = data.BIC;
    Psy_BIC_all = [Psy_BIC_all; temp_BIC];
    N_Psy_BIC(i,:) = mean(temp_BIC);
end

all_subject = [];
for i = 1:length(analysis_folder)
    temp = ones(Ndata(i),1) * i;
    all_subject = [all_subject; temp];
end

use_RL = 3;
use_psycho = 5;

mean(N_BIC)
mean(N_Psy_BIC)
%Compare with session base
BIC_all = BIC_all(:,use_RL);
Psy_BIC_all = Psy_BIC_all(:,use_psycho);

lme_opt_sabun = fitlme_analysis_20210520_0(BIC_all-Psy_BIC_all,all_subject);
lme_opt_sabun(1).lme
lme_opt_sabun(2).lme

%Compare with mouse base
N_BIC = N_BIC(:,use_RL);
N_Psy_BIC = N_Psy_BIC(:,use_psycho);
signrank(N_BIC,N_Psy_BIC)

%% fig 2c
cd('G:\upload_code\Figure2\Fig2c');
sdata = struct();% source data 
temp = Psy_BIC_all - BIC_all;
temp = sort(temp);
figure
plot(temp,1:length(temp))
set(gca,'ylim',[0 length(temp)+1])
sdata.x = temp;
sdata.y = (1:length(temp))';
T = struct2table(sdata);
writetable(T, 'source fig2c.csv');

