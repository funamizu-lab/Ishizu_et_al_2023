%{
----------------------------------------------------------------------------
Analyzing behavioral data
At least for the correct rate
----------------------------------------------------------------------------
%}

function FigureS1d_Dual_ML_analysis_220818_para_all

analysis_folder = {
    'G:\Ishizu_data\Tokyo_ephys_ishizu\only_all_behaviors\a04_behave\ML_CumGauss_Fit'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\only_all_behaviors\a08_behave\ML_CumGauss_Fit'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\only_all_behaviors\i20_behave\ML_CumGauss_Fit'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\only_all_behaviors\i24_behave\ML_CumGauss_Fit'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\only_all_behaviors\i34_behave\ML_CumGauss_Fit'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\only_all_behaviors\i35_behave\ML_CumGauss_Fit'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\only_all_behaviors\i43_behave\ML_CumGauss_Fit'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\only_all_behaviors\i46_behave\ML_CumGauss_Fit'
    };

likeli1 = 'ML_CumGauss_Fit_20220114.mat';
likeli2 = 'ML_CumGauss_Fit_20220114_lapse.mat';

log_likeli_all = [];
BIC_all = [];
for i = 1:length(analysis_folder)
    cd(analysis_folder{i})
    data1 = load(likeli1);
    data2 = load(likeli2);
    
    data1_LL = data1.log_likeli;
    data2_LL = data2.log_likeli;
    data1_BIC = data1.BIC;
    data2_BIC = data2.BIC;
    
    use_data = [
        data1_LL(:,1),... %unbiased
        data1_LL(:,2),... %block threshold
        data1_LL(:,5),... %block thre + duration sensitivity
        data1_LL(:,4),... %block thre + duration thre
        data1_LL(:,3),... %block thre + block sensitibity
        data2_LL(:,4),... %block thre + lapse
        data1_LL(:,6),... %block thre + duration sensitivity + duration thre
        data1_LL(:,7),... %test1
        data1_LL(:,8),... %test2
        ];
    use_BIC = [
        data1_BIC(:,1),... %unbiased
        data1_BIC(:,2),... %block threshold
        data1_BIC(:,5),... %block thre + duration sensitivity
        data1_BIC(:,4),... %block thre + duration thre
        data1_BIC(:,3),... %block thre + block sensitibity
        data2_BIC(:,4),... %block thre + lapse
        data1_BIC(:,6),... %block thre + duration sensitivity + duration thre
        data1_BIC(:,7),... %test1
        data1_BIC(:,8),... %test2
        ];

    log_likeli_session(i,:) = mean(use_data);
    log_likeli_all = [log_likeli_all; use_data];
    BIC_session(i,:) = mean(use_BIC);
    BIC_all = [BIC_all; use_BIC];
end

%Plot the difference of log likelihood
plot_log_likeli(:,1) = log_likeli_all(:,2) - log_likeli_all(:,1);
plot_log_likeli(:,2) = log_likeli_all(:,3) - log_likeli_all(:,2);
plot_log_likeli(:,3) = log_likeli_all(:,4) - log_likeli_all(:,2);
plot_log_likeli(:,4) = log_likeli_all(:,5) - log_likeli_all(:,2);
plot_log_likeli(:,5) = log_likeli_all(:,6) - log_likeli_all(:,2);
plot_log_likeli(:,6) = log_likeli_all(:,7) - log_likeli_all(:,3);

%%% Fig S1d
cd('G:\upload_code\FigureS1\FigS1d');
temp_x = (rand(length(plot_log_likeli(:,1)),1) - 0.5) * 0.1 + 1;
figure
sdata = struct();% source data
subplot(1,3,1)
boxplot(plot_log_likeli(:,1))
hold on
plot(temp_x,plot_log_likeli(:,1),'k.')
sdata.data = plot_log_likeli(:,1);
T = struct2table(sdata);
writetable(T, 'source figS1d leftpannel.csv');

sdata = struct();% source data
subplot(1,3,2)
plot(plot_log_likeli(:,2:5)')
hold on
boxplot(plot_log_likeli(:,2:5))
sdata.data1 = plot_log_likeli(:,2);
sdata.data2 = plot_log_likeli(:,3);
sdata.data3 = plot_log_likeli(:,4);
sdata.data4 = plot_log_likeli(:,5);
T = struct2table(sdata);
writetable(T, 'source figS1d middlepannel.csv');

sdata = struct();% source data
subplot(1,3,3)
boxplot(plot_log_likeli(:,6))
hold on
plot(temp_x,plot_log_likeli(:,6),'k.')
sdata.data = plot_log_likeli(:,6);
T = struct2table(sdata);
writetable(T, 'source figS1d rightpannel.csv');



% %Based on the mean log likelihood, calculate the likelihood ratio test
% [~,p] = lratiotest(mean_N_loglikeli(2),mean_N_loglikeli(1),1)
% [~,p] = lratiotest(mean_N_loglikeli(3),mean_N_loglikeli(2),1)
% [~,p] = lratiotest(mean_N_loglikeli(4),mean_N_loglikeli(2),1)
% [~,p] = lratiotest(mean_N_loglikeli(5),mean_N_loglikeli(2),1)
% [~,p] = lratiotest(mean_N_loglikeli(6),mean_N_loglikeli(2),1)
% [~,p] = lratiotest(mean_N_loglikeli(7),mean_N_loglikeli(3),1)
