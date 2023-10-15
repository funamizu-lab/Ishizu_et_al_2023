%{
----------------------------------------------------------------------------
Analyzing behavioral data
At least for the correct rate
----------------------------------------------------------------------------
%}

function FigureS1g_Dual_RL_analysis_220314_para_compare_220912
close all
analysis_folder = {
    'G:\Ishizu_data\Tokyo_ephys_ishizu\only_all_behaviors\a04_behave\ML_CumGauss_Fit\RL_20220818'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\only_all_behaviors\a08_behave\ML_CumGauss_Fit\RL_20220818'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\only_all_behaviors\i20_behave\ML_CumGauss_Fit\RL_20220818'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\only_all_behaviors\i24_behave\ML_CumGauss_Fit\RL_20220818'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\only_all_behaviors\i34_behave\ML_CumGauss_Fit\RL_20220818'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\only_all_behaviors\i35_behave\ML_CumGauss_Fit\RL_20220818'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\only_all_behaviors\i43_behave\ML_CumGauss_Fit\RL_20220818'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\only_all_behaviors\i46_behave\ML_CumGauss_Fit\RL_20220818'
    };
analysis_folder2 = {
    'G:\Ishizu_data\Tokyo_ephys_ishizu\only_all_behaviors\a04_behave\ML_CumGauss_Fit\RL_20220118'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\only_all_behaviors\a08_behave\ML_CumGauss_Fit\RL_20220118'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\only_all_behaviors\i20_behave\ML_CumGauss_Fit\RL_20220118'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\only_all_behaviors\i24_behave\ML_CumGauss_Fit\RL_20220118'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\only_all_behaviors\i34_behave\ML_CumGauss_Fit\RL_20220118'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\only_all_behaviors\i35_behave\ML_CumGauss_Fit\RL_20220118'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\only_all_behaviors\i43_behave\ML_CumGauss_Fit\RL_20220118'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\only_all_behaviors\i46_behave\ML_CumGauss_Fit\RL_20220118'
    };

use_para=3;
log_likeli_all = get_likelihood_value(analysis_folder);
log_likeli_all2= get_likelihood_value(analysis_folder2);

log_likeli_all = [log_likeli_all(:,use_para),log_likeli_all2(:,use_para)];
% signrank(log_likeli_all(:,1),log_likeli_all(:,2))

test = log_likeli_all(:,1)-log_likeli_all(:,2);

%% FigS1 g
figure
hold on
boxplot(test)
plot((rand(length(test),1)-0.5) * 0.1 + 1, test, 'k.')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function log_likeli_all =  get_likelihood_value(analysis_folder)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

count = 0;
for i = 1:length(analysis_folder)
    temp_folder = analysis_folder{i};
    cd(temp_folder);
    filename1 = dir('*.mat');
    Ndata(i,1) = length(filename1);
    
    clear sub_para sub_likeli sub_BIC
    for filecount = 1 : length(filename1)
        temp_filename = filename1(filecount).name;
        fpath = temp_filename;
    
        load(fpath)

        count = count + 1;
        log_likeli_all(count,:) = log_likeli;
    end
end

all_subject = [];
for i = 1:length(analysis_folder)
    temp = ones(Ndata(i),1) * i;
    all_subject = [all_subject; temp];
end

return