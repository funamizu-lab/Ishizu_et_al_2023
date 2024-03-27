
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
function State_dynamics_20220917_QR_CV_depth_control(pathname)

switch nargin
    case 0
        pathname   = pwd;
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

% temp = dir('RL_20220818*'); %RL_20220118 frame
% if length(temp) ~= 1
%     temp
%     hoge
% end
% load(temp.name);
% %'ave_likeli' ,'BIC','log_likeli','para_max','N_trial'

temp = dir('task_frame*'); %Task frame
if length(temp) ~= 1
    temp
    hoge
end
load(temp.name);

temp = dir('Bpod*'); %Bpod
if length(temp) ~= 1
    temp
    hoge
end
load(temp.name);
Bpod_file = temp.name;

cd(spike_dir)
tif_name = dir('task_spike*.mat'); %get all the tif files
length_tif = length(tif_name);
cd(pathname)

temp = dir('depth_spike_20220517*');
if length(temp) == 1
    load(temp.name);
    %spike_depth def_depth length_neuron
    if length_tif ~= length_neuron
        [length_tif length_neuron]
        hoge
    end
    depth_neuron = find(spike_depth <= def_depth);
elseif length(temp) == 0
    depth_neuron = [1:length_tif]; %Use all the neurons
else
    temp
    hoge
end

% %Get parameter of RL model
% % [prior,posterior,ave_likeli,likelihood,Q,Choice_prob,para] = ...
% %     Dual_RL_model_block1_20220118_para_determined(Bpod_file,para_max(2,:),N_trial);
% [prior,posterior,ave_likeli,likelihood,Long_posterior,Short_posterior,...
%     prior_value,long_value,short_value,para] = ...
%     Dual_RL_model_block1_20220314_para_determined(Bpod_file,para_max(3,:),N_trial);
% prior = prior_value(:,2) ./ (prior_value(:,1) + prior_value(:,2));
% posterior = posterior(:,2) ./ sum(posterior(:,1)+posterior(:,2));
% binary_prior = zeros(length(prior),1);
% binary_posterior = zeros(length(posterior),1);
% temp = find(prior > 0.5);
% binary_prior(temp) = 1;
% temp = find(posterior > 0.5);
% binary_posterior(temp) = 1;

[minD_trial,maxD_trial,Choice_trial,tone_evidence,trial_evidence,use_trial2,use_trial3,use_trial_all,...
    low,high,correct,error,flip_tone,number_use_trial,...
    binary_tone,right_trial_all,number_trial_all,right_trial,number_trial] ...
 = Dual_get_basic_task_structure_20210204(Bpod_file);

%Chech the number of trials 
if length(all_trial_time) ~= length(Outcome)
    [length(all_trial_time), length(Outcome)]
    hoge
end

%Get task parameter
%Choice_trial = find(Outcome == 1 | Outcome == 2);
%low  = find(Correct_side == 0);
%high = find(Correct_side == 1);
Reward = Correct_side == Chosen_side;
Reward = double(Reward);
correct = find(Reward == 1);
error = find(Reward == 0);
left  = find(Chosen_side == 0);
right = find(Chosen_side == 1);
stim_length = unique(StimDuration);
Long  = find(StimDuration == stim_length(2));
Short = find(StimDuration == stim_length(1));

frame_choice2 = nan(length(frame_choice),1);
frame_choice2(left) = frame_choice(left,1);
frame_choice2(right) = frame_choice(right,2);

%frame_sound_off
frame_sound_off = frame_sound;
frame_sound_off(Long) = frame_sound_off(Long) + 1000; %Add 1000 ms
frame_sound_off(Short) = frame_sound_off(Short) + 200; %Add 200 ms 

for i = 1:length(tone_evidence)
    temp = find(trial_evidence == tone_evidence(i));
    evi_trial(i).matrix = temp;
end

%Analyze the block LR
if BlockProb(2) ~= BlockProb(3) %Block change task
    if BlockProb(2) > BlockProb(3) % Right -> Left
        block_R = use_trial2;
        block_L = use_trial3;
        
    else % Left -> Right
        block_L = use_trial2;
        block_R = use_trial3;
        
    end
else %Reward change task
    if BlockReward(2,1) < BlockReward(2,2) % Right -> Left
        block_R = use_trial2;
        block_L = use_trial3;
    else % Left -> Right
        block_L = use_trial2;
        block_R = use_trial3;
    end
end   

use_block_sequence = zeros(length(use_trial_all),1);
for i = 1:length(block_R)
    temp = block_R(i);
    temp = find(use_trial_all == temp);
    if length(temp) ~= 0
        use_block_sequence(temp) = 1;
    end
end

% block2 = use_trial2;
% block3 = use_trial3;
use_trial = use_trial_all;
binary_prior = use_block_sequence;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Setup the prior and posterior
if length(binary_prior) ~= length(use_trial)
    hoge
end
% if length(binary_posterior) ~= length(use_trial)
%     hoge
% end
StimDuration23 = StimDuration(use_trial);
long23_prior  = find(StimDuration23 == stim_length(2));
short23_prior = find(StimDuration23 == stim_length(1));

% binary_prior_long23 = binary_prior(long23_prior);
% binary_prior_short23 = binary_prior(short23_prior);
% binary_posterior_long23 = binary_posterior(long23_prior);
% binary_posterior_short23 = binary_posterior(short23_prior);
% binary_prior0 = find(binary_prior == 0);
% binary_prior1 = find(binary_prior == 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% low23  = intersect(low,use_trial);
% high23 = intersect(high,use_trial);
% left23  = intersect(left,use_trial);
% right23 = intersect(right,use_trial);
% use_frame = frame_sound(use_trial);
short23 = intersect(Short,use_trial);
long23 = intersect(Long,use_trial);
long23_correct = intersect(long23,correct);

cd(spike_dir)

tif_name = dir('task_spike*.mat'); %get all the tif files
length_tif = length(tif_name);
%max_tif = length_tif;
max_tif = length(depth_neuron);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make time for short sound and long sound trials
time_window_prior = 200;
time_window_sound = 200;
time_window_choice_pre = 200;
time_window_choice_post = 500;
pre_sound  = 2000;
post_sound = 4000;
    
for i = 1:length(use_trial)
    trial_activity(i).matrix = nan(max_tif,post_sound+pre_sound);
end
for file_count = 1:max_tif
    neuron_norm_activ(file_count).matrix = nan(length(use_trial),3);
end

for file_count = 1:max_tif
    [file_count, max_tif]
    %temp_file = tif_name(file_count).name;
    %temp_file
    %temp_file = sprintf('task_spike_stripe20210520_%d',sig_neuron(file_count));
    %temp_file = sprintf('task_spike_stripe20210520_%d', file_count);
    temp_file = sprintf('task_spike_stripe20210520_%d',depth_neuron(file_count));
    temp_file
    
    %clear data
    data = load(temp_file); %spike_mark
    spike_mark = data.spike_mark;
    if std(spike_mark) ~= 0
        norm_spike_mark = (spike_mark - mean(spike_mark)) ./ std(spike_mark); %Normalized activity
    else
        norm_spike_mark = spike_mark;
    end
    ave_norm_spike_mark = norm_spike_mark(1:length(norm_spike_mark)-100+1);
    for i = 2:100
        ave_norm_spike_mark = ave_norm_spike_mark + norm_spike_mark(i:length(norm_spike_mark)-100+i);
    end
    ave_norm_spike_mark = ave_norm_spike_mark ./ 100;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %All trials for QR decomposition
    norm_activ = nan(length(use_trial),4);
    %norm_time_activ = nan(length(use_trial),post_sound+pre_sound);
    for i = 1:length(use_trial)
        temp_frame1 = frame_sound(use_trial(i));
        temp_frame2 = frame_choice2(use_trial(i));
        temp_frame3 = frame_sound_off(use_trial(i));
        temp_prior = norm_spike_mark(temp_frame1-time_window_prior : temp_frame1-1);
        temp_sound = norm_spike_mark(temp_frame1 : temp_frame1+time_window_sound-1);
        temp_sound_off = norm_spike_mark(temp_frame3-time_window_sound : temp_frame3-1);
        temp_choice = norm_spike_mark(temp_frame2-time_window_choice_pre : temp_frame2+time_window_choice_post-1);
        
        norm_activ(i,:) = [mean(temp_prior),mean(temp_sound),mean(temp_sound_off),mean(temp_choice)];
        %trial_activity(i).matrix(file_count,:) = norm_spike_mark(temp_frame1-pre_sound : temp_frame1+post_sound-1);
        trial_activity(i).matrix(file_count,:) = ave_norm_spike_mark(temp_frame1-pre_sound : temp_frame1+post_sound-1);
    end
    neuron_norm_activ(file_count).matrix = norm_activ;
%     prior_sabun  = mean(norm_activ(binary_prior1,1)) - mean(norm_activ(binary_prior0,1)); %Get Prior sabun
%     sound_sabun  = mean(norm_activ(sound23_1,2))     - mean(norm_activ(sound23_0,2)); %Get Sound sabun
%     choice_sabun = mean(norm_activ(choice23_1,3))    - mean(norm_activ(choice23_0,3)); %Get Choice sabun
%     QR_activity(file_count,:) = [prior_sabun,sound_sabun,choice_sabun];
end
cd(pathname)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make Cross validation to make Q and R and 
%Make CV block
%Use only trial with long and correct
block = TrialBlock(long23_correct);
length_trial = length(long23_correct);
temp_group = zeros(length_trial,1);

evi_long23_correct = trial_evidence(long23_correct);
clear evi_trial
for i = 1:length(tone_evidence)
    temp = find(evi_long23_correct == tone_evidence(i));
    evi_trial(i).matrix = temp;
end

binary_prior_long23_correct = nan(length(long23_correct),1);
use_long23_correct = nan(length(long23_correct),1);
for i = 1:length(long23_correct)
    temp = find(use_trial == long23_correct(i));
    binary_prior_long23_correct(i) = binary_prior(temp);
    use_long23_correct(i) = temp;
end
trial_activity = trial_activity(use_long23_correct);
binary_prior0 = find(binary_prior_long23_correct == 0);
binary_prior1 = find(binary_prior_long23_correct == 1);

%Get sound and choice during use_trial
temp_sound = Correct_side(long23_correct);
sound23_0 = find(temp_sound == 0);
sound23_1 = find(temp_sound == 1);
% temp_choice = Chosen_side(long23_correct); %Sound and choice is equal
% choice23_0 = find(temp_choice == 0);
% choice23_1 = find(temp_choice == 1);

CV_repeat = 100; %Repeat CV for 100 times
CV_para = 10;
temp_trial = [1:length_trial];
if max_tif >= 2

parfor CV = 1:CV_repeat
%for CV = 1:CV_repeat
    [CV, CV_repeat]
    %Make CV trials
    temp_group = make_CV_block(block,CV_para);
    
    prior_mode = nan(length_trial,post_sound+pre_sound);
    sound_mode = nan(length_trial,post_sound+pre_sound);
    %choice_mode = nan(length_trial,post_sound+pre_sound);
    for i = 1:CV_para
        %clear ave_likeli correct_rate log_likeli
        ave_likeli = [];
        correct_rate = [];
        log_likeli = [];
        vali_trial = find(temp_group == i);
        train_trial = setdiff(temp_trial,vali_trial);
    
        CV_prior1 = intersect(binary_prior1, train_trial);
        CV_prior0 = intersect(binary_prior0, train_trial);
        CV_sound1 = intersect(sound23_1, train_trial);
        CV_sound0 = intersect(sound23_0, train_trial);
        % CV_choice1 = intersect(choice23_1, train_trial);
        % CV_choice0 = intersect(choice23_0, train_trial);
    
        QR_activity = nan(max_tif,2);
        for file_count = 1:max_tif
            norm_activ = neuron_norm_activ(file_count).matrix;
            norm_activ = norm_activ(use_long23_correct,:);
            prior_sabun  = mean(norm_activ(CV_prior1,1)) - mean(norm_activ(CV_prior0,1)); %Get Prior sabun
            sound_sabun  = mean(norm_activ(CV_sound1,3)) - mean(norm_activ(CV_sound0,3)); %Get Sound sabun
            %choice_sabun = mean(norm_activ(CV_choice1,4)) - mean(norm_activ(CV_choice0,4)); %Get Choice sabun
            %QR_activity(file_count,:) = [prior_sabun,sound_sabun,choice_sabun];
            QR_activity(file_count,:) = [prior_sabun,sound_sabun];
        end
    
        [Q,R] = qr(QR_activity);
        %Q_CV(i).matrix = Q;
    
        for j = 1:length(vali_trial)
            prior_mode(vali_trial(j),:) = Q(:,1)' * trial_activity(vali_trial(j)).matrix;
            sound_mode(vali_trial(j),:) = Q(:,2)' * trial_activity(vali_trial(j)).matrix;
            %choice_mode(vali_trial(j),:) = Q(:,3)' * trial_activity(vali_trial(j)).matrix;
        end
    end
    prior_CV(CV).matrix = prior_mode;
    sound_CV(CV).matrix = sound_mode;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get the trial_based prior_mode_CV
prior_mode_CV = nan(length_trial,post_sound+pre_sound);
sound_mode_CV = nan(length_trial,post_sound+pre_sound);
for j = 1:length_trial
    prior_trial = nan(CV_repeat,post_sound+pre_sound);
    sound_trial = nan(CV_repeat,post_sound+pre_sound);
    for CV = 1:CV_repeat
        prior_trial(CV,:) = prior_CV(CV).matrix(j,:);
        sound_trial(CV,:) = sound_CV(CV).matrix(j,:);
    end
    prior_mode_CV(j,:) = mean(prior_trial);
    sound_mode_CV(j,:) = mean(sound_trial);
end

%QR for the training data
QR_activity = nan(max_tif,2);
for file_count = 1:max_tif
    norm_activ = neuron_norm_activ(file_count).matrix;
    norm_activ = norm_activ(use_long23_correct,:);
    prior_sabun  = mean(norm_activ(binary_prior1,1)) - mean(norm_activ(binary_prior0,1)); %Get Prior sabun
    sound_sabun  = mean(norm_activ(sound23_1,3)) - mean(norm_activ(sound23_0,3)); %Get Sound sabun
    %choice_sabun = mean(norm_activ(CV_choice1,4)) - mean(norm_activ(CV_choice0,4)); %Get Choice sabun
    %QR_activity(file_count,:) = [prior_sabun,sound_sabun,choice_sabun];
    QR_activity(file_count,:) = [prior_sabun,sound_sabun];
end
    
[Q,R] = qr(QR_activity);
size(Q)
%Prior mode: Q(:,1)
%Stimulus mode: Q(:,2)
%Choice mode: Q(:,3)
        
prior_mode = nan(length_trial,post_sound+pre_sound);
sound_mode = nan(length_trial,post_sound+pre_sound);
% choice_mode = nan(length_trial,post_sound+pre_sound);
for i = 1:length_trial
    prior_mode(i,:) = Q(:,1)' * trial_activity(i).matrix;
    sound_mode(i,:) = Q(:,2)' * trial_activity(i).matrix;
%     choice_mode(i,:) = Q(:,3)' * trial_activity(i).matrix;
end

else
    max_tif
    prior_mode = [];
    sound_mode = [];
    prior_mode_CV = [];
    sound_mode_CV = [];
    Q = [];
    Q_CV = [];
    R = [];
    QR_activity = [];
end
% save State_dynamics_20220519_QR_CV_depth prior_mode sound_mode Q Q_CV R CV_para QR_activity ...
%     binary_prior1 binary_prior0 sound23_1 sound23_0 evi_trial
save State_dynamics_20220917_QR_CV_depth prior_mode sound_mode prior_mode_CV sound_mode_CV ...
    Q R CV_para CV_repeat QR_activity ...
    binary_prior1 binary_prior0 sound23_1 sound23_0 evi_trial

% figure
% plot(mean(prior_mode(binary_prior1,:)),'r')
% hold on
% plot(mean(prior_mode(binary_prior0,:)),'b')
% figure
% plot(mean(sound_mode(sound23_1,:)),'r')
% hold on
% plot(mean(sound_mode(sound23_0,:)),'b')
% figure
% plot(mean(choice_mode(choice23_1,:)),'r')
% hold on
% plot(mean(choice_mode(choice23_0,:)),'b')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function temp_group = make_CV_block(block,CV_para)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
function [ROC_sound,p_sound,temp_shuffle_ROC] = get_ROC_shuffle(sound_frame, low, high, shuffle_time)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ROC_sound = make_auROC_20210620(sound_frame, low, high); %left choice

temp_shuffle_ROC  = [];
parfor i = 1:shuffle_time
%for i = 1:shuffle_time
    shuffle_frame = sound_frame(randperm(length(sound_frame)));
    temp_shuffle_ROC(1,i) = make_auROC_20210620(shuffle_frame, low, high);
end
temp_sound = find(temp_shuffle_ROC > ROC_sound); %how many random beyound the ROC
p_sound = length(temp_sound) ./ shuffle_time;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mean_prob = make_auROC_20210620(activity, low, high) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

low_activity = sort(activity(low));
high_activity = sort(activity(high));
low_activity = round(low_activity,3);
high_activity = round(high_activity,3);

if isempty(low) || isempty(high)
    mean_prob = nan;
else
    for i = 1:length(low_activity)
        temp = find(high_activity > low_activity(i));
        temp05 = find(high_activity == low_activity(i));
        prob(i) = (length(temp)+length(temp05)/2) ./ length(high_activity);
    end
    mean_prob = mean(prob);
    if mean_prob < 0.5
        mean_prob = 1 - mean_prob;
    end
end

max_activity = max(activity);
min_activity = min(activity);
if mean_prob > 0.9
    hist_low = histcounts(low_activity,[min_activity:0.001:max_activity]);
    hist_high = histcounts(high_activity,[min_activity:0.001:max_activity]);
    %[min_activity,max_activity]
% %     hist_low
% %     hist_high
% %     prob
%     figure
%     plot(hist_low,'b')
%     hold on
%     plot(hist_high,'r')
%     %hoge
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sound_choice_index = get_sound_choice_index(spike_sound,trial_low,trial_high,trial_left,trial_right)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    spike_low = mean(spike_sound(trial_low,2));
    spike_high = mean(spike_sound(trial_high,2));
    spike_left = mean(spike_sound(trial_left,2));
    spike_right = mean(spike_sound(trial_right,2));
    sound_choice_index(1) = (spike_low - spike_high) ./ (spike_low + spike_high);
    sound_choice_index(2) = (spike_left - spike_right) ./ (spike_left + spike_right);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spike_count,p,spike_trace] = get_sound_response3(spike_mark, spike_filter, frame_sound, pre_frame, post_frame, pre_frame2, post_frame2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spike_count = nan(length(frame_sound),2);
spike_trace = nan(length(frame_sound),pre_frame2+post_frame2);

for i = 1:length(frame_sound)
    
    if isnan(frame_sound(i)) == 0
        temp_pre  = [frame_sound(i)-pre_frame : frame_sound(i)-1];
        temp_post = [frame_sound(i) : frame_sound(i)+post_frame-1];
        temp_all = [frame_sound(i)-pre_frame2 : frame_sound(i)+post_frame2-1];

        temp_pre = spike_mark(temp_pre);
        temp_post = spike_mark(temp_post);
        spike_count(i,:) = [sum(temp_pre), sum(temp_post)];
    
        spike_trace(i,:) = spike_filter(temp_all);
    else
        spike_count(i,:) = nan(1,2);
        spike_trace(i,:) = nan(1,post_frame2+pre_frame2);
    end        
end

p = signrank(spike_count(:,1),spike_count(:,2));

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [frame_behave] = thre_detection_frame(behave_time, blue_scan)    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    temp = find(isnan(behave_time) == 1);
    if isempty(temp),
        for i = 1:length(behave_time),
            temp_time = behave_time(i);
            %frame_behave(i) = find(blue_scan >= temp_time, 1);
            temp_behave = find(blue_scan >= temp_time, 1);
            %Rotary encoder does not stop by task!!
            if length(temp_behave) == 0,
                temp_time
                blue_scan(length(blue_scan))
            end
            frame_behave(i) = temp_behave;
        end
    else
        frame_behave = nan;
    end

    return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [length_scan, scan_time, first_scan] = thre_detection(temp_scan, t, thre_scan)    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %temp_scan = ch(i,:);
    temp = find(temp_scan > thre_scan); %above thre
    temp_plus = [-1, temp]; %Add the first value
    temp1 = temp_plus(1:length(temp_plus)-1);
    temp2 = temp_plus(2:length(temp_plus));
    temp_sabun = temp2 - temp1;
    temp_sabun1 = find(temp_sabun > 1);
    temp_sabun2 = find(temp_sabun == 1);
    temp_sabun2 = temp_sabun2 - 1;
    temp_use = intersect(temp_sabun1, temp_sabun2);
    
    if ~isempty(temp_use),
        temp_use = temp(temp_use);
    %     length_scan(i) = length(temp_use);
    %     scan_time(i).matrix = t(temp_use);
    %     first_scan(i) = t(temp_use(1));
        length_scan = length(temp_use);
        scan_time = t(temp_use);
        first_scan = t(temp_use(1));
    else
        length_scan = nan;
        scan_time = nan;
        first_scan = nan;
    end
    
    return
