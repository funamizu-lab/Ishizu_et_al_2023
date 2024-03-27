
%{
----------------------------------------------------------------------------
First_take number of frames in each tif files
Analyzing imaging data simply
At least for the correct rate
----------------------------------------------------------------------------
%}
function Task_kaiseki_tokyo_encode_20230915_GLM_each_window(pathname)

switch nargin
    case 0
        pathname = pwd;
    case 1
        disp('OK to analyze')
    otherwise
        hoge
end
cd(pathname)
%spike_dir = 'spike_ch1';
spike_dir = dir('spike_ch*');
if length(spike_dir) ~= 1
    hoge
end
spike_dir = spike_dir.name

%new_p_thre = 100;
new_p_thre = 10;
%new_p_thre = 5;
%new_p_thre = 0;

% [filename1, pathname1,findex]=uigetfile('*.*','frame file');
% filename1 = [pathname1,filename1];
% load(filename1)
temp = dir('task_frame*');
if length(temp) ~= 1
    temp
    hoge
end
load(temp.name);
%frame_start
%frame_sound
%frame_end

temp = dir('Bpod*');
if length(temp) ~= 1
    temp
    hoge
end
load(temp.name);
Bpod_file = temp.name;

[minD_trial,maxD_trial,Choice_trial,tone_evidence,trial_evidence,block2,block3,use_trial,...
    low,high,correct,error,flip_tone,number_use_trial,...
    binary_tone,right_trial_all,number_trial_all,right_trial,number_trial] ...
 = Dual_get_basic_task_structure_20210204(Bpod_file);

temp = dir('RL_20220818*'); %RL_20220118 frame
if length(temp) ~= 1
    temp
    hoge
end
load(temp.name);
%'ave_likeli' ,'BIC','log_likeli','para_max','N_trial'
%Get parameter of RL model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get parameter of RL model
[prior,posterior,ave_likeli,likelihood,Long_posterior,Short_posterior,...
    prior_value,long_value,short_value,para,conf_Q] = ...
    Dual_RL_model_block1_20220703_para_determined(Bpod_file,para_max(3,:),N_trial);
prior = prior_value(:,2) ./ (prior_value(:,1) + prior_value(:,2));
posterior = posterior(:,2) ./ (posterior(:,1)+posterior(:,2));
conf_Q = conf_Q(:,2) ./ (conf_Q(:,1)+conf_Q(:,2)); %relative
binary_prior = zeros(length(prior),1);
binary_posterior = zeros(length(posterior),1);
binary_conf = zeros(length(conf_Q),1);
temp = find(prior > 0.5);
binary_prior(temp) = 1;
temp = find(posterior > 0.5);
binary_posterior(temp) = 1;
temp = find(conf_Q > 0.5);
binary_conf(temp) = 1;
%prior posterior conf_Q binary_prior binary_posterior binary_conf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Setup the prior and posterior
if length(binary_prior) ~= length(use_trial)
    hoge
end
% all_binary_prior = nan(length(Choice_trial),1);
% all_binary_prior(use_trial) = binary_prior;

[length(binary_tone), length(StimDuration)]; %same number of trials

stim_length = unique(StimDuration);
%Focus on the use_trial
StimDuration = StimDuration(use_trial);
Chosen_side = Chosen_side(use_trial);
Correct_side = Correct_side(use_trial);
frame_choice = frame_choice(use_trial,:);
binary_tone = binary_tone(use_trial);
block = TrialBlock(use_trial);

%Get task parameter
low  = find(Correct_side == 0);
high = find(Correct_side == 1);
Long  = find(StimDuration == stim_length(2));
Short = find(StimDuration == stim_length(1));
left = find(Chosen_side == 0);
right = find(Chosen_side == 1);
correct_error = Correct_side == Chosen_side;
correct = find(correct_error == 1);
error = find(correct_error == 0);

Long_correct = intersect(Long, correct);
Short_correct = intersect(Short, correct);

%Make the choice timing frames
frame_sound = frame_sound(use_trial);
frame_spout = frame_spout(use_trial,:);

frame_choice_select = nan(length(frame_choice),1);
frame_choice_select(left) = frame_choice(left,1);
frame_choice_select(right) = frame_choice(right,2);

time_window_base = 200;
time_window = 100; %ms
time_long_pre  = 1500; %4sec
time_long_post = 2500; %4sec
time_long_pre2  = 500; %4sec
time_long_post2 = 2500; %4sec

time_short_pre = 1500; %3.2sec
time_short_post = 1700; %3.2sec
time_short_pre2 = 500; %3.2sec
time_short_post2 = 2500; %3.2sec

cd(spike_dir);

%Just focus on the 3 phases, before, initial and end of sound
%save the raw data with and the Q-values
long_use_timing = [15, 17, 25]; %before, first and end of sound
short_use_timing = [15, 17]; %before, first and end of sound

tif_name = dir('task_spike_stripe*.mat'); %get all the tif files
length_tif = length(tif_name);
max_tif = length_tif;

%Use only correct trial (Figure 3G)
%check whether the activity is correlated to choice, prior, action-value,
%sensory evidence, sound category (5 parameters)

[long_RMSE,wheel_good(1)] = ...
    regress_each_time_20230915(long_use_timing, Long,frame_sound,frame_choice_select,frame_spout,max_tif,...
    time_long_pre,time_long_post,time_long_pre2,time_long_post2,time_window,time_window_base, ...
    prior, posterior, conf_Q, binary_prior, binary_posterior, binary_conf, ...
    Chosen_side, binary_tone, ave_velocity, block, Correct_side);

[short_RMSE,wheel_good(2)] = ...
    regress_each_time_20230915(short_use_timing, Short,frame_sound,frame_choice_select,frame_spout,max_tif,...
    time_short_pre,time_short_post,time_short_pre2,time_short_post2,time_window,time_window_base, ...
    prior, posterior, conf_Q, binary_prior, binary_posterior, binary_conf, ...
    Chosen_side, binary_tone, ave_velocity, block, Correct_side);

[long_RMSE_correct,wheel_good(3)] = ...
    regress_each_time_20230915(long_use_timing, Long_correct,frame_sound,frame_choice_select,frame_spout,max_tif,...
    time_long_pre,time_long_post,time_long_pre2,time_long_post2,time_window,time_window_base, ...
    prior, posterior, conf_Q, binary_prior, binary_posterior, binary_conf, ...
    Chosen_side, binary_tone, ave_velocity, block, Correct_side);

[short_RMSE_correct,wheel_good(4)] = ...
    regress_each_time_20230915(short_use_timing, Short_correct,frame_sound,frame_choice_select,frame_spout,max_tif,...
    time_short_pre,time_short_post,time_short_pre2,time_short_post2,time_window,time_window_base, ...
    prior, posterior, conf_Q, binary_prior, binary_posterior, binary_conf, ...
    Chosen_side, binary_tone, ave_velocity, block, Correct_side);

cd(pathname)

save GLM_20230915_each_window ...
    long_RMSE short_RMSE long_RMSE_correct short_RMSE_correct wheel_good

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [all_neuron_error,wheel_good] = ...
    regress_each_time_20230915(long_use_timing, use_long,frame_sound,frame_choice_select,frame_spout,max_tif,...
    time_long_pre,time_long_post,time_long_pre2,time_long_post2,time_window,time_window_base, ...
    prior, posterior, conf_Q, binary_prior, binary_posterior, binary_conf, ...
    choice, binary_tone, ave_velocity, block, Correct_side)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_long  = round((time_long_pre + time_long_post) ./ time_window);
time_long2 = round((time_long_pre2 + time_long_post2) ./ time_window);

frame_sound = frame_sound(use_long,:);
frame_choice_select = frame_choice_select(use_long,:);

frame_spout = frame_spout(use_long,:);
% base_for_spout_on  = frame_spout(:,1)-time_window_base; %Before moving spout
% base_for_spout_off = frame_spout(:,3)-time_window_base; %Before moving spout

%Choice, binary_tone, rotation_speed, prior, conf_Q
prior = prior(use_long);
posterior = posterior(use_long);
conf_Q = conf_Q(use_long);
binary_prior = binary_prior(use_long);
binary_posterior = binary_posterior(use_long);
binary_conf = binary_conf(use_long);
choice = choice(use_long);
binary_tone = binary_tone(use_long);
block = block(use_long);
Correct_side = Correct_side(use_long);

%rotation_speed
x_const = [ones(length(use_long),1), binary_tone, choice, prior, conf_Q];
x_const = double(x_const);

for i = 1:time_long
    frame_sound_use_on(:,i) = frame_sound - time_long_pre + (i-1)*time_window;
    %frame_sound_use_off(:,i) = frame_sound - time_long_pre + i*time_window - 1;
end
for i = 1:time_long2
    frame_sound_use2_on(:,i) = frame_choice_select - time_long_pre2 + (i-1)*time_window;
    %frame_sound_use2_off(:,i) = frame_choice_select - time_long_pre2 + i*time_window - 1;
end

ave_velocity = double(ave_velocity);

rotation1 = zeros(length(use_long),time_long);
wheel_good = 1;
%if length(ave_velocity) ~= 1
if length(ave_velocity) > frame_sound(end)
    for i = 1:length(use_long) %trial
        for j = 1:time_long %time
            if isnan(frame_sound_use_on(i,j)) == 0
                temp_frame = [frame_sound_use_on(i,j) : frame_sound_use_on(i,j)+time_window-1];
                rotation1(i,j) = mean(ave_velocity(temp_frame));
            end
        end
    end
else
    wheel_good = 0;
end

CV_repeat = 10;
CV_para = 10; %10 fold
size_para = 6;

parfor file_count = 1:max_tif
%for file_count = 1:max_tif
    [file_count, max_tif]
    %temp_file = tif_name(file_count).name;
    %temp_file
    temp_file = sprintf('task_spike_stripe20210520_%d',file_count);
    temp_file
    
    %clear data
    data = load(temp_file); %spike_mark
    
    spike_mark = data.spike_mark;
    spike_multi(file_count) = data.single_multi;
    
    %Normalize the spike_mark
    mean_spike = mean(spike_mark);
    std_spike = std(spike_mark);
    if std_spike ~= 0
        spike_mark = (spike_mark - mean_spike) ./ std_spike; %normalized the activity
    else
        spike_mark = zeros(1,length(spike_mark));
    end
    
    temp_spike_frame1 = nan(length(use_long),time_long);
    
    for i = 1:length(use_long) %trial
        for j = 1:time_long %time
            if isnan(frame_sound_use_on(i,j)) == 0
                temp_frame = [frame_sound_use_on(i,j) : frame_sound_use_on(i,j)+time_window-1];
                temp_spike_frame1(i,j) = mean(spike_mark(temp_frame));
            end
        end
%         for j = 1:time_long2
%             if isnan(frame_sound_use2_on(i,j)) == 0
%                 temp_frame = [frame_sound_use2_on(i,j) : frame_sound_use2_on(i,j)+time_window-1];
%                 spike_frame2(i,j) = mean(spike_mark(temp_frame));
%             end
%         end
    end
    %regression analysis here
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %start of cross validation based on trial
    error1 = [];
    error2 = [];
    error3 = [];
    error4 = [];
    error5 = [];
    %repeat this CV_repeat time
    for n = 1:CV_repeat
        temp_group = make_CV_group(use_long, block,CV_para);
    
        all_error = [];
        for j = 1:time_long %time
            spike_data = temp_spike_frame1(:,j);
        
            if std(spike_data) ~= 0
                min_spike = min(spike_data);
                spike_data = spike_data - min_spike + 1; %Minimum becomes 1
                spike_data = log(spike_data); %Make log to fit with standard regression
            
                %Start CV
                %check whether the activity is correlated to choice, prior, action-value,
                %sensory evidence, sound category (5 parameters)
                %choice
                all_error(1,j) = make_CV_easy_regress(spike_data, temp_group, CV_para, choice);
                %prior
                all_error(2,j) = make_CV_easy_regress(spike_data, temp_group, CV_para, prior);
                %conf_Q
                all_error(3,j) = make_CV_easy_regress(spike_data, temp_group, CV_para, conf_Q);
                %sensory evidence
                all_error(4,j) = make_CV_easy_regress(spike_data, temp_group, CV_para, binary_tone);
                %sound category
                all_error(5,j) = make_CV_easy_regress(spike_data, temp_group, CV_para, Correct_side);
            else
                all_error(:,j) = nan(5,1);
            end
        end
        error1(n,:) = all_error(1,:); 
        error2(n,:) = all_error(2,:); 
        error3(n,:) = all_error(3,:); 
        error4(n,:) = all_error(4,:); 
        error5(n,:) = all_error(5,:); 
    end
    integ_error = [nanmean(error1); nanmean(error2); nanmean(error3); nanmean(error4); nanmean(error5)];
    all_neuron_error(file_count).matrix = integ_error;
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function all_error = make_CV_easy_regress(Activ, temp_group, CV_para, Sound)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

all_error = nan(length(Sound), 1);

temp_trial = [1:length(temp_group)];
for i = 1:CV_para
    clear ave_likeli correct_rate log_likeli
    vali_trial = find(temp_group == i);
    train_trial = setdiff(temp_trial,vali_trial);
   
    vali_activ = Activ(vali_trial,:);
    temp_activ = Activ(train_trial,:);

    vali_sound = Sound(vali_trial,:);
    temp_sound = Sound(train_trial);

    %vali_one = ones(length(vali_trial),1);
    temp_one = ones(length(train_trial),1);
    %regress
    b = regress(temp_activ, [temp_sound, temp_one]);
    
    temp = vali_sound * b(1) + b(2);
    %RSS
    temp = (temp - vali_activ) .* (temp - vali_activ);
    
    for j = 1:length(vali_trial)
        all_error(vali_trial(j)) = temp(j);
    end
end
all_error = sqrt(sum(all_error) ./ length(Sound));

return




%temp_group
all_likelihood = zeros(length_trial,1);
temp_trial = [1:length_trial];

clear b_block
for i = 1:CV_para

    
    rand0 = randperm(length(temp0));
    rand0 = rand0(1:min_temp);
    temp0 = temp0(rand0);
    
    rand1 = randperm(length(temp1));
    rand1 = rand1(1:min_temp);
    temp1 = temp1(rand1);
    use_temp = [temp0;temp1];
    use_temp = sort(use_temp);
    data_length_CV(i,:) = [length(temp_sound), length(use_temp)];
    
%     CVerr = cvglmnet(temp_activ,temp_sound,'binomial',[],'class',CV_para,[],[],[],[]);
    %unique(temp_sound(use_temp))
    size(Activ)
    size(temp_activ(use_temp,:))
    size(temp_sound(use_temp))
    CVerr = cvglmnet(temp_activ(use_temp,:),temp_sound(use_temp),'binomial',[],'class',CV_para,[],[],[],[]);
    temp_b = cvglmnetCoef(CVerr, 'lambda_min');
    CVpred = cvglmnetPredict(CVerr, vali_activ, 'lambda_min');
    CVerr_CV(i).matrix = CVerr;
    
    x = temp_b(2:length(temp_b));
    c = temp_b(1);
    

%     %test with train data
%     CVpred_train = cvglmnetPredict(CVerr, temp_activ, 0);
%     figure
%     plot(temp_sound, CVpred_train, 'b.')
    
    new_b(:,i) = temp_b;
end
return




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
function temp_group = make_CV_group(use_long, block,CV_para)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%Make CV block
length_trial = length(use_long);
temp_group = zeros(length_trial,1);
for i = min(block):max(block)
    temp_block = find(block == i);
    temp_length_trial = length(temp_block);
    
    clear block_group
    temp_trial = [1:temp_length_trial];
    if temp_length_trial > CV_para,
        block_group = rem(temp_trial,CV_para) + 1;
        block_group = block_group(randperm(temp_length_trial));
    else
        block_group = [1:temp_length_trial];
    end
    temp_group(temp_block) = block_group;
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function make_logGLM(spike_data,CV_para)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Makes the minimum value of spike_data as 1:
min_spike = min(spike_data);
spike_data = spike_data - min_spike + 1; %Minimum becomes 1
spike_data = log(spike_data); %Make log to fit with poisson regression
kernel_all = double(kernel_all);

%Get the best lambda by LOO cross validation and use the best lambda for
%the analysis
%[y_test,SumError,lambda,x_opt,CVerr_CV] = CV_RidgeKenelRegression(spike_data,kernel_all,lambda,CV_para,options);

%For all trials
CVerr_all = cvglmnet(kernel_all',spike_data,'gaussian',options,'deviance',CV_para,[],[],[],[]); %now the spike is LOG!!





