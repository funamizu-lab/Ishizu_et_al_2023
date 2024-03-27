%{
----------------------------------------------------------------------------
Analyzing behavioral data
At least for the correct rate
----------------------------------------------------------------------------
%}

function Figure1d_Dual_Full_psycho_analysis_220818_all

analysis_folder = {
    'G:\Ishizu_data\IntermediateFiles\only_all_behaviors\a04_behave\Full_psycho_20220114'
    'G:\Ishizu_data\IntermediateFiles\only_all_behaviors\a08_behave\Full_psycho_20220114'
    'G:\Ishizu_data\IntermediateFiles\only_all_behaviors\i20_behave\Full_psycho_20220114'
    'G:\Ishizu_data\IntermediateFiles\only_all_behaviors\i24_behave\Full_psycho_20220114'
    'G:\Ishizu_data\IntermediateFiles\only_all_behaviors\i34_behave\Full_psycho_20220114'
    'G:\Ishizu_data\IntermediateFiles\only_all_behaviors\i35_behave\Full_psycho_20220114'
    'G:\Ishizu_data\IntermediateFiles\only_all_behaviors\i43_behave\Full_psycho_20220114'
    'G:\Ishizu_data\IntermediateFiles\only_all_behaviors\i46_behave\Full_psycho_20220114'
    };


count = 0;
opt_sabun = [];
for i = 1:length(analysis_folder)
    temp_folder = analysis_folder{i};
    cd(temp_folder);
    filename1 = dir('*.mat');
    Ndata(i,1) = length(filename1);
    
    clear sub_sabun
    for filecount = 1 : length(filename1)
        temp_filename = filename1(filecount).name;
        temp_file{filecount} = temp_filename;
        fpath = temp_filename;
    
        load(fpath)
    %'ave_likeli' ,'BIC','log_likeli','para_max','N_trial'
    
        count = count + 1;
        mean_opt_L_min = mean(full_psycho.opt_L_min);
        mean_opt_R_min = mean(full_psycho.opt_R_min);
        mean_opt_L_max = mean(full_psycho.opt_L_max);
        mean_opt_R_max = mean(full_psycho.opt_R_max);
        sub_sabun(filecount,:) = [mean_opt_R_min-mean_opt_L_min,mean_opt_R_max-mean_opt_L_max];
    end
    N_opt_sabun(i,:) = mean(sub_sabun);
    opt_sabun = [opt_sabun; sub_sabun];
end

all_subject = [];
for i = 1:length(analysis_folder)
    temp = ones(Ndata(i),1) * i;
    all_subject = [all_subject; temp];
end

lme_opt_sabun = fitlme_analysis_20210520_0(opt_sabun(:,1)-opt_sabun(:,2),all_subject);
lme_opt_sabun(1).lme
lme_opt_sabun(2).lme

signrank(N_opt_sabun(:,1),N_opt_sabun(:,2))

%%% Fig 1d %%%
cd('G:\upload_code\Figure1\Fig1d');
figure
sdata = struct;% source data
hold on
boxplot(opt_sabun)
plot(opt_sabun')
sdata.short=opt_sabun(:,1);
sdata.long =opt_sabun(:,2);
T = struct2table(sdata);
writetable(T, 'source fig1d top.csv');
