% a04
auc_a04 = {
    'G:\Ishizu_data\Tokyo_ephys_ishizu\auditory\a04\2021-02-13_a04_1PPC2AC_Left_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\auditory\a04\2021-02-18_a04_1AC2PPC_Right\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\auditory\a04\2021-02-19_a04_1AC2PPC_Right\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\auditory\a04\2021-02-20_a04_AC_Right_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\auditory\a04\2021-02-21_a04_AC_Right_OK\recording1_task'
};
fof_a04 = {
    'G:\Ishizu_data\Tokyo_ephys_ishizu\fof\a04\2021-02-27_a04_FOF_Right_chunkerror_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\fof\a04\2021-02-28_a04_FOF_Right_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\fof\a04\2021-03-01_a04_FOF_Right_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\fof\a04\2021-03-10_a04_FOF_Left_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\fof\a04\2021-03-11_a04_FOF_Left_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\fof\a04\2021-03-12_a04_FOF_Left_OK\recording1_task'
};
mpfc_a04 = {
    'G:\Ishizu_data\Tokyo_ephys_ishizu\mpfc\a04\2021-03-04_a04_mPFC_Right_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\mpfc\a04\2021-03-05_a04_mPFC_Right_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\mpfc\a04\2021-03-06_a04_mPFC_Right_OK\recording2_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\mpfc\a04\2021-03-07_a04_mPFC_Right_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\mpfc\a04\2021-03-15_a04_mPFC_Left_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\mpfc\a04\2021-03-16_a04_mPFC_Left_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\mpfc\a04\2021-03-17_a04_mPFC_Left_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\mpfc\a04\2021-03-18_a04_mPFC_Left_OK\recording1_task'
};

[auc_sound,auc_choice,auc_prior, auc_prob, auc_neuron] = process_20231002_encoding_model2_depth(auc_a04);
[fof_sound,fof_choice,fof_prior, fof_prob, fof_neuron] = process_20231002_encoding_model2_depth(fof_a04);
[mpfc_sound,mpfc_choice,mpfc_prior, mpfc_prob, mpfc_neuron] = process_20231002_encoding_model2_depth(mpfc_a04);
close all

prob_1 = [auc_prob;fof_prob;mpfc_prob];

% a08
auc_a08 = {
    'G:\Ishizu_data\Tokyo_ephys_ishizu\auditory\a08\2021-02-09_a08_1PPC2AC_Left_OK\recording2_task'
};
fof_a08 = {
    'G:\Ishizu_data\Tokyo_ephys_ishizu\fof\a08\2021-02-12_a08_FOF_Right_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\fof\a08\2021-02-13_a08_FOF_Right_OK\recording2_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\fof\a08\2021-02-19_a08_FOF_Left_OK\recording2_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\fof\a08\2021-02-20_a08_FOF_Left_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\fof\a08\2021-02-23_a08_FOF_Left_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\fof\a08\2021-02-26_a08_FOF_Left_OK\recording3_task'
};
mpfc_a08 = {
    'G:\Ishizu_data\Tokyo_ephys_ishizu\mpfc\a08\2021-02-16_a08_mPFC_Right_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\mpfc\a08\2021-02-17_a08_mPFC_Right_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\mpfc\a08\2021-02-21_a08_mPFC_Left_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\mpfc\a08\2021-02-22_a08_mPFC_Left_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\mpfc\a08\2021-02-24_a08_mPFC_Left_OK\recording1_task'
};

[auc_sound,auc_choice,auc_prior, auc_prob, auc_neuron] = process_20231002_encoding_model2_depth(auc_a08);
[fof_sound,fof_choice,fof_prior, fof_prob, fof_neuron] = process_20231002_encoding_model2_depth(fof_a08);
[mpfc_sound,mpfc_choice,mpfc_prior, mpfc_prob, mpfc_neuron] = process_20231002_encoding_model2_depth(mpfc_a08);
close all

prob_2 = [auc_prob;fof_prob;mpfc_prob];

% i20
auc_i20 = {};
fof_i20 = {
    'G:\Ishizu_data\Tokyo_ephys_ishizu\fof\i20\2021-04-29_i20_FOF_Right_chunkerror_OK\recording4_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\fof\i20\2021-04-30_i20_FOF_Right_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\fof\i20\2021-05-13_i20_FOF_Left_OK\recording2_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\fof\i20\2021-05-14_i20_FOF_Left_chunkerror_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\fof\i20\2021-05-15_i20_FOF_Left_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\fof\i20\2021-05-16_i20_FOF_Left_OK\recording2_task'
};
mpfc_i20 = {
    'G:\Ishizu_data\Tokyo_ephys_ishizu\mpfc\i20\2021-05-08_i20_mPFC_Right_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\mpfc\i20\2021-05-09_i20_mPFC_Right_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\mpfc\i20\2021-05-10_i20_mPFC_Right_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\mpfc\i20\2021-05-18_i20_mPFC_Left_OK\recording2_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\mpfc\i20\2021-05-19_i20_mPFC_Left_OK\recording3_task'
};

[fof_sound,fof_choice,fof_prior, fof_prob, fof_neuron] = process_20231002_encoding_model2_depth(fof_i20);
[mpfc_sound,mpfc_choice,mpfc_prior, mpfc_prob, mpfc_neuron] = process_20231002_encoding_model2_depth(mpfc_i20);
close all

prob_3 = [NaN    NaN    NaN;fof_prob;mpfc_prob];

% i24
auc_i24 = {
    'G:\Ishizu_data\Tokyo_ephys_ishizu\auditory\i24\2021-10-05_i24_AC_Left_chunk_error_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\auditory\i24\2021-10-06_i24_AC_left_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\auditory\i24\2021-10-07_i24_AC_left_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\auditory\i24\2021-10-08_i24_AC_left_OK\recording2_task'
};
fof_i24 = {
    'G:\Ishizu_data\Tokyo_ephys_ishizu\fof\i24\2021-08-03_i24_FOF_Right_OK\recording2_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\fof\i24\2021-08-04_i24_FOF_Right_OK\recording2_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\fof\i24\2021-08-07_i24_FOF_Right_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\fof\i24\2021-09-29_i24_FOF_left_OK\recording2_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\fof\i24\2021-10-01_i24_FOF_left_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\fof\i24\2021-10-02_i24_FOF_left_OK\recording1_task'
};
mpfc_i24 = {
    'G:\Ishizu_data\Tokyo_ephys_ishizu\mpfc\i24\2021-08-05_i24_mPFC_Right_OK\recording2_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\mpfc\i24\2021-08-06_i24_mPFC_Right_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\mpfc\i24\2021-09-30_i24_mPFC_left_OK\recording1_task'
};

[auc_sound,auc_choice,auc_prior, auc_prob, auc_neuron] = process_20231002_encoding_model2_depth(auc_i24);
[fof_sound,fof_choice,fof_prior, fof_prob, fof_neuron] = process_20231002_encoding_model2_depth(fof_i24);
[mpfc_sound,mpfc_choice,mpfc_prior, mpfc_prob, mpfc_neuron] = process_20231002_encoding_model2_depth(mpfc_i24);
close all

prob_4 = [auc_prob;fof_prob;mpfc_prob];

% i34
auc_i34 = {
    'G:\Ishizu_data\Tokyo_ephys_ishizu\auditory\i34\2021-11-17_i34_AC_left_OK\recording2_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\auditory\i34\2021-11-18_i34_AC_left_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\auditory\i34\2021-11-19_i34_AC_left_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\auditory\i34\2021-12-08_i34_AC_Right_OK\recording2_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\auditory\i34\2021-12-09_i34_AC_Right_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\auditory\i34\2021-12-10_i34_AC_Right_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\auditory\i34\2021-12-11_i34_AC_Right_OK\recording1_task'
};
fof_i34 = {
    'G:\Ishizu_data\Tokyo_ephys_ishizu\fof\i34\2021-11-10_i34_FOF_left_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\fof\i34\2021-12-01_i34_FOF_Right_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\fof\i34\2021-12-02_i34_FOF_Right_OK\recording2_task'
};
mpfc_i34 = {
    'G:\Ishizu_data\Tokyo_ephys_ishizu\mpfc\i34\2021-11-14_i34_mPFC_Left_OK\recording2_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\mpfc\i34\2021-12-03_i34_mPFC_Right_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\mpfc\i34\2021-12-05_i34_mPFC_Right_OK\recording1_task'
};

[auc_sound,auc_choice,auc_prior, auc_prob, auc_neuron] = process_20231002_encoding_model2_depth(auc_i34);
[fof_sound,fof_choice,fof_prior, fof_prob, fof_neuron] = process_20231002_encoding_model2_depth(fof_i34);
[mpfc_sound,mpfc_choice,mpfc_prior, mpfc_prob, mpfc_neuron] = process_20231002_encoding_model2_depth(mpfc_i34);
close all

prob_5 = [auc_prob;fof_prob;mpfc_prob];

% i35
auc_i35 = {
    'G:\Ishizu_data\Tokyo_ephys_ishizu\auditory\i35\2021-11-10_i35_AC_Left_OK\recording2_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\auditory\i35\2021-11-11_i35_AC_left_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\auditory\i35\2021-11-12_i35_AC_left_OK\recording3_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\auditory\i35\2021-11-13_i35_AC_Left_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\auditory\i35\2021-12-01_i35_AC_Right_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\auditory\i35\2021-12-02_i35_AC_Right_OK\recording3_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\auditory\i35\2021-12-03_i35_AC_Right_OK\recording1_task'
};
fof_i35 = {
    'G:\Ishizu_data\Tokyo_ephys_ishizu\fof\i35\2021-11-03_i35_FOF_Left_OK\recording2_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\fof\i35\2021-11-04_i35_FOF_Left_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\fof\i35\2021-11-26_i35_FOF_Right_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\fof\i35\2021-11-27_i35_FOF_Right_OK\recording1_task'
};
mpfc_i35 = {
    'G:\Ishizu_data\Tokyo_ephys_ishizu\mpfc\i35\2021-11-05_i35_mPFC_Left_OK\recording2_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\mpfc\i35\2021-11-06_i35_mPFC_Left_OK\recording2_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\mpfc\i35\2021-11-24_i35_mPFC_Right_OK\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\mpfc\i35\2021-11-25_i35_mPFC_Right_OK\recording1_task'
};

[auc_sound,auc_choice,auc_prior, auc_prob, auc_neuron] = process_20231002_encoding_model2_depth(auc_i35);
[fof_sound,fof_choice,fof_prior, fof_prob, fof_neuron] = process_20231002_encoding_model2_depth(fof_i35);
[mpfc_sound,mpfc_choice,mpfc_prior, mpfc_prob, mpfc_neuron] = process_20231002_encoding_model2_depth(mpfc_i35);
close all

prob_6 = [auc_prob;fof_prob;mpfc_prob];

% i43
auc_i43 = {
    'G:\Ishizu_data\Tokyo_ephys_ishizu\auditory\i43\2022-06-11_i43_1AC2PPC_Right\recording2_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\auditory\i43\2022-06-15_i43_AC_Left\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\auditory\i43\2022-06-16_i43_1AC2PPC_Left\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\auditory\i43\2022-06-17_i43_1AC2PPC_Left\recording1_task'
};
fof_i43 = {
    'G:\Ishizu_data\Tokyo_ephys_ishizu\fof\i43\2022-06-22_i43_FOF_Left\recording2_task'
};
mpfc_i43 = {
    'G:\Ishizu_data\Tokyo_ephys_ishizu\mpfc\i43\2022-06-03_i43_mPFC_Right\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\mpfc\i43\2022-06-04_i43_mPFC_Right\recording4_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\mpfc\i43\2022-06-24_i43_mPFC_Left\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\mpfc\i43\2022-06-27_i43_mPFC_Left\recording1_task'
};

[auc_sound,auc_choice,auc_prior, auc_prob, auc_neuron] = process_20231002_encoding_model2_depth(auc_i43);
[fof_sound,fof_choice,fof_prior, fof_prob, fof_neuron] = process_20231002_encoding_model2_depth(fof_i43);
[mpfc_sound,mpfc_choice,mpfc_prior, mpfc_prob, mpfc_neuron] = process_20231002_encoding_model2_depth(mpfc_i43);
close all

prob_7 = [auc_prob;fof_prob;mpfc_prob];

% i46
auc_i46 = {
    'G:\Ishizu_data\Tokyo_ephys_ishizu\auditory\i46\2022-06-08_i46_AC_Left\recording1_task'
};
fof_i46 = {
    'G:\Ishizu_data\Tokyo_ephys_ishizu\fof\i46\2022-06-01_i46_FOF_Left\recording1_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\fof\i46\2022-06-02_i46_FOF_Left\recording2_task'
};
mpfc_i46 = {
    'G:\Ishizu_data\Tokyo_ephys_ishizu\mpfc\i46\2022-06-03_i46_mPFC_Left\recording2_task'
    'G:\Ishizu_data\Tokyo_ephys_ishizu\mpfc\i46\2022-06-04_i46_mPFC_Left\recording1_task'
};

[auc_sound,auc_choice,auc_prior, auc_prob, auc_neuron] = process_20231002_encoding_model2_depth(auc_i46);
[fof_sound,fof_choice,fof_prior, fof_prob, fof_neuron] = process_20231002_encoding_model2_depth(fof_i46);
[mpfc_sound,mpfc_choice,mpfc_prior, mpfc_prob, mpfc_neuron] = process_20231002_encoding_model2_depth(mpfc_i46);
close all

prob_8 = [auc_prob;fof_prob;mpfc_prob];


prob_all(:,:,1) = prob_1;
prob_all(:,:,2) = prob_2;
prob_all(:,:,3) = prob_3;
prob_all(:,:,4) = prob_4;
prob_all(:,:,5) = prob_5;
prob_all(:,:,6) = prob_6;
prob_all(:,:,7) = prob_7;
prob_all(:,:,8) = prob_8;

[num_areas, num_groups, num_samples] = size(prob_all);
prob_sum = sum(prob_all,3);

for i = 1:num_samples
    auc_prob(i,:) = prob_all(1,:,i);
    fof_prob(i,:) = prob_all(2,:,i);
    mpfc_prob(i,:) = prob_all(3,:,i);
end
mean_auc_prob = mean(auc_prob,1,"omitnan");
mean_fof_prob = mean(fof_prob,1);
mean_mpfc_prob = mean(mpfc_prob,1);

%%% FigS12a %%%
cd('G:\upload_code\FigureS12\FigS12a');
figure;
sdata=struct();
sdata.label={'sound';'choice';'prior'};
sdata.barvalue_ac =mean_auc_prob';
sdata.barvalue_m2 =mean_fof_prob';
sdata.barvalue_mpfc =mean_mpfc_prob';
bar([mean_auc_prob;mean_fof_prob;mean_mpfc_prob]')
ylim([0 0.5])
hold on;
% plot each data
for group = 1:num_groups
    x = [group-0.225, group, group+0.225];
    for sample = 1:num_samples
        plot(x,[auc_prob(sample,group), fof_prob(sample, group), mpfc_prob(sample, group)], 'k.-', 'MarkerSize',10);
    end
end
for sample = 1:num_samples
    eval(['sdata.mouse',num2str(sample),'_ac=transpose(auc_prob(sample,:));']);
    eval(['sdata.mouse',num2str(sample),'_m2=transpose(fof_prob(sample,:));']);
    eval(['sdata.mouse',num2str(sample),'_mpfc=transpose(mpfc_prob(sample,:));']);
end
hold off;

p_sound = signrank_3regions(auc_prob,fof_prob,mpfc_prob,1);
p_choice = signrank_3regions(auc_prob,fof_prob,mpfc_prob,2);
p_prior = signrank_3regions(auc_prob,fof_prob,mpfc_prob,3);

p_sound
p_choice
p_prior

%%% source data %%%
T = struct2table(sdata);
writetable(T, 'source fig S12a.csv');

%%%%%%%%%%%%%%%%%%%%%%%%
function p = signrank_3regions(auc_prob,fof_prob,mpfc_prob,use_number)
%%%%%%%%%%%%%%%%%%%%%%%%
a = auc_prob(:,use_number);
b = fof_prob(:,use_number);
c = mpfc_prob(:,use_number);
p(1) = signrank(a,b);
p(2) = signrank(a,c);
p(3) = signrank(b,c);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = kai2test_3regions(auc_neuron,fof_neuron,mpfc_neuron,use_number)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
use_auc =  [auc_neuron(use_number), auc_neuron(4)- auc_neuron(use_number)];
use_fof =  [fof_neuron(use_number), fof_neuron(4)- fof_neuron(use_number)];
use_mpfc = [mpfc_neuron(use_number),mpfc_neuron(4)-mpfc_neuron(use_number)];
p(1) = kai2_test(use_auc,use_fof);
p(2) = kai2_test(use_auc,use_mpfc);
p(3) = kai2_test(use_mpfc,use_fof);
%p
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p,median_sound] = ranksum_areas(auc_sound,fof_sound,mpfc_sound)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p(1)=ranksum(auc_sound,fof_sound);
p(2)=ranksum(auc_sound,mpfc_sound);
p(3)=ranksum(fof_sound,mpfc_sound);

length_auc = length(auc_sound);
length_fof = length(fof_sound);
length_mpfc = length(mpfc_sound);
length_auc = [1:length_auc] ./ length_auc;
length_fof = [1:length_fof] ./ length_fof;
length_mpfc = [1:length_mpfc] ./ length_mpfc;

hist_x = [-inf,-600:10:-400,inf];
hist_auc = histcounts(auc_sound,hist_x);
hist_fof = histcounts(fof_sound,hist_x);
hist_mpfc = histcounts(mpfc_sound,hist_x);

figure
subplot(2,1,1)
plot(sort(auc_sound),length_auc,'r')
hold on
plot(sort(fof_sound),length_fof,'g')
hold on
plot(sort(mpfc_sound),length_mpfc,'b')
subplot(2,1,2)
plot(hist_auc./sum(hist_auc),'r')
hold on
plot(hist_fof./sum(hist_fof),'g')
hold on
plot(hist_mpfc./sum(hist_mpfc),'b')
set(gca,'xlim',[0 length(hist_x)])
set(gca,'xtick',[1:2:length(hist_x)])
% boxplot([auc_sound;fof_sound;mpfc_sound],[ones(length(auc_sound),1);ones(length(fof_sound),1)*2;ones(length(mpfc_sound),1)*3])

median_sound = [median(auc_sound),median(fof_sound),median(mpfc_sound)];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [BIC_sound,BIC_choice,BIC_prior,prob_sig,neuron_sig] = process_20231002_encoding_model2_depth(folders)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% analysis_dir = eval(folders);
analysis_dir = folders;
analysis_dir

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
get_regress_plot(BIC_all, y_sound_all,y_choice_all,y_prior_all)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function get_regress_plot(BIC, y_sound,y_choice,y_prior)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nan check
nanBIC = mean(BIC,2);
nanBIC = find(isnan(nanBIC) == 0);
BIC = BIC(nanBIC,:);
y_sound = y_sound(nanBIC,:);
y_choice = y_choice(nanBIC,:);
y_prior = y_prior(nanBIC,:);

%regression contains both the left and right side, need to split into half.
y_sound = get_value_regressor(y_sound);
y_choice = get_value_regressor(y_choice);
y_prior = get_value_regressor(y_prior);

y_sound = abs(y_sound);
y_choice = abs(y_choice);
y_prior = abs(y_prior);

figure
plot(mean(y_sound),'r')
hold on
plot(mean(y_choice),'g')
hold on
plot(mean(y_prior),'b')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function use_y = get_value_regressor(y_sound)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[length_BIC,length_y] = size(y_sound);
for i = 1:length_BIC
    temp_y = y_sound(i,:);
    temp1 = temp_y(1:length_y/2);
    temp2 = temp_y(length_y/2 + 1:length_y);
    
    if max(temp1) >= max(temp2)
        use_y(i,:) = temp1;
    else
        use_y(i,:) = temp2;
    end
end
end

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

temp_x = [1:length(LL_sound)];

figure
subplot(2,2,1)
plot(sort(LL_sound),temp_x,'r')
hold on
plot(sort(LL_choice),temp_x,'g')
hold on
plot(sort(LL_prior),temp_x,'b')
subplot(2,2,2)
plot(sort(BIC_sound),temp_x,'r')
hold on
plot(sort(BIC_choice),temp_x,'g')
hold on
plot(sort(BIC_prior),temp_x,'b')
subplot(2,2,3)
plot(sort(LL_sound),'r')
hold on
plot(sort(LL_choice),'g')
hold on
plot(sort(LL_prior),'b')
subplot(2,2,4)
plot(sort(BIC_sound),'r')
hold on
plot(sort(BIC_choice),'g')
hold on
plot(sort(BIC_prior),'b')

p_LL_sound = 1-chi2cdf(2*LL_sound, parameter_size-parameter_sound);
p_LL_sound = -log10(p_LL_sound);
p_LL_choice = 1-chi2cdf(2*LL_choice, parameter_size-parameter_choice);
p_LL_choice = -log10(p_LL_choice);
p_LL_prior = 1-chi2cdf(2*LL_prior, parameter_size-parameter_prior);
p_LL_prior = -log10(p_LL_prior);

figure
plot(sort(p_LL_sound),'r')
hold on
plot(sort(p_LL_choice),'g')
hold on
plot(sort(p_LL_prior),'b')

sig_sound = find(p_LL_sound > 2);
sig_choice = find(p_LL_choice > 2);
sig_prior = find(p_LL_prior > 2);

sig_sc = intersect(sig_sound,sig_choice);
sig_cp = intersect(sig_prior,sig_choice);
sig_sp = intersect(sig_prior,sig_sound);
sig_scp = intersect(sig_sc,sig_prior);

length_scp = length(sig_scp);
length_sc = length(sig_sc) - length(sig_scp);
length_sp = length(sig_sp) - length(sig_scp);
length_cp = length(sig_cp) - length(sig_scp);
length_s = length(sig_sound)-length_sc-length_sp-length_scp;
length_c = length(sig_choice)-length_sc-length_cp-length_scp;
length_p = length(sig_prior)-length_sp-length_cp-length_scp;

[length_s, length_c, length_p]
[length_sc, length_sp, length_cp, length_scp]

temp_sig = [length(sig_sound),length(sig_choice),length(sig_prior)];
neuron_sig = [temp_sig,length(p_LL_sound)];
prob_sig = temp_sig ./ length(p_LL_sound)

end
