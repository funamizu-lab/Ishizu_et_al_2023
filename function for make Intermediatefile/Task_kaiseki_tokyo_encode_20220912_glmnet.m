
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
%Compare sound, choice prior during long sound (-100 ms to 1500 ms for long)
%Compare sound, choice prior during short sound (-100 ms to 700 ms for short)
%p = 0.001
%Each epoch, predict the prior, sensory and choice (integration)
----------------------------------------------------------------------------
%}
function Task_kaiseki_tokyo_encode_20220912_glmnet(pathname)
%p = genpath('/home/funamizu/Tokyo_ephys/Tokyo_ephys_program')
%addpath(p);

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

Time_stretch = 100; %100 ms
duration_long_sound  = 1000;
duration_short_sound = 200;

% %Open each spike data
% [pathname1] = uigetdir
% temp_cd = ['cd ',pathname1];
% eval(temp_cd); %Go to record folder

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

%Get parameter of RL model
%Use number 3 which seems to achive the best BIC
temp = dir('RL_20220818*'); %RL_20220118 frame
if length(temp) ~= 1
    temp
    hoge
end
load(temp.name);
%'ave_likeli' ,'BIC','log_likeli','para_max','N_trial'

% [prior,posterior,ave_likeli,likelihood,Q,Choice_prob,para] = ...
%     Dual_RL_model_block1_20220118_para_determined(Bpod_file,para_max(2,:),N_trial);
[prior,posterior,ave_likeli,likelihood,Long_posterior,Short_posterior,...
    prior_value,long_value,short_value,para] = ...
    Dual_RL_model_block1_20220314_para_determined(Bpod_file,para_max(3,:),N_trial);

prior = prior(:,2) ./ sum(prior(:,1)+prior(:,2));
posterior = posterior(:,2) ./ sum(posterior(:,1)+posterior(:,2));
Q = prior_value;
if length(Q) ~= length(prior)
    [length(Q), length(prior)]
    hoge
end
prior_value = prior_value(:,2) ./ (prior_value(:,1)+prior_value(:,2)); %relative value
binary_prior = zeros(length(prior_value),1);
temp = find(prior_value >= 0.5);
binary_prior(temp) = 1;
% binary_posterior = zeros(length(posterior),1);
% temp = find(prior > 0.5);
% binary_prior(temp) = 1;
% temp = find(posterior > 0.5);
% binary_posterior(temp) = 1;

[minD_trial,maxD_trial,Choice_trial,tone_evidence,trial_evidence,block2,block3,use_trial,...
    low,high,correct,error,flip_tone,number_use_trial,...
    binary_tone,right_trial_all,number_trial_all,right_trial,number_trial] ...
 = Dual_get_basic_task_structure_20210204_2(Bpod_file);

all_binary_prior = nan(length(Correct_side),1);
for i = 1:length(use_trial)
    all_binary_prior(use_trial(i)) = binary_prior(i);
end
prior0 = find(all_binary_prior == 0);
prior1 = find(all_binary_prior == 1);

%Chech the number of trials 
if length(all_trial_time) ~= length(Outcome)
    [length(all_trial_time), length(Outcome)]
    hoge
end
%Setup the prior and posterior
if length(prior) ~= length(use_trial)
    hoge
end
if length(posterior) ~= length(use_trial)
    hoge
end
%Make Q_choice
temp_Chosen = Chosen_side(use_trial); %0or1
for i = 1:length(temp_Chosen)
    Q_choice(i,1) = Q(i,temp_Chosen(i)+1);
end

%Get task parameter
%Choice_trial = find(Outcome == 1 | Outcome == 2);
low  = find(Correct_side == 0);
high = find(Correct_side == 1);
left  = find(Chosen_side == 0);
right = find(Chosen_side == 1);
CorrectError = Correct_side == Chosen_side;
%correct = find(CorrectError == 1);
%error = find(CorrectError == 0);
leftC = intersect(left, correct);
leftE = intersect(left, error);
rightC = intersect(right, correct);
rightE = intersect(right, error);
leftC = intersect(leftC, Choice_trial);
leftE = intersect(leftE, Choice_trial);
rightC = intersect(rightC, Choice_trial);
rightE = intersect(rightE, Choice_trial);

stim_length = unique(StimDuration);
long  = find(StimDuration == stim_length(2));
short = find(StimDuration == stim_length(1));
lowL = intersect(low,long);
lowS = intersect(low,short);
highL = intersect(high,long);
highS = intersect(high,short);
leftL = intersect(left, long);
leftS = intersect(left, short);
rightL = intersect(right, long);
rightS = intersect(right, short);

prior0L = intersect(prior0, long);
prior0S = intersect(prior0, short);
prior1L = intersect(prior1, long);
prior1S = intersect(prior1, short);

for i = 1:length(tone_evidence)
    temp = find(trial_evidence == tone_evidence(i));
    evi_trial(i).matrix = temp;
end

%Make kernels
%6 Stimuli
%Stimulus category
%Choice
%Reward
%Choice x Reward
%Prior
%Previous stimulus category
%Previous choice
%Previous Reward
%Previous choice x reward
%ave_velocity

%Bin is made at every 5 ms
frame5_time = ceil(length(ave_velocity) / Time_stretch);
check_time = ceil(frame_end(end) / Time_stretch);
[frame5_time, check_time]
if check_time > frame5_time
    disp('change frame5_time')
    frame5_time = check_time + 10; %1sec
end

frame_sound5 = ceil(frame_sound ./ Time_stretch);
% frame_sound_off5 = frame_sound5;
% frame_sound_off5(long) = frame_sound_off5(long) + duration_long_sound/Time_stretch;
% frame_sound_off5(short) = frame_sound_off5(short) + duration_short_sound/Time_stretch;
frame_choice5 = ceil(frame_choice ./ Time_stretch);
temp_frame_choice = nan(length(frame_choice5),1);
temp_frame_choice(left) = frame_choice5(left,1);
temp_frame_choice(right) = frame_choice5(right,2);
frame_choice5 = temp_frame_choice;
clear temp_frame_choice

frame_start5 = ceil(frame_start ./ Time_stretch);
frame_spout5 = ceil(frame_spout ./ Time_stretch); %on off on off

%Check the time between frame_choice to spout move
temp = frame_spout5([2:length(frame_spout5)],1) - frame_choice5([1:length(frame_spout5)-1]);
[min(temp), max(temp)]

%Make task state kernels
frame5_state = false(6,frame5_time); 
for i = 1:length(frame_start)
    state_time(i,1).matrix = [frame_spout5(i,1) : frame_sound5(i)-1]; %Start to sound start
    state_time(i,2).matrix = [frame_sound5(i) : frame_spout5(i,3)-1]; %Sound start to sound move
%    state_time(i,2).matrix = [frame_sound5(i) : frame_sound_off5(i)-1]; %Sound start to off
%    state_time(i,3).matrix = [frame_sound_off5(i) : frame_spout5(i,3)-1]; %Sound off to spout move
    state_time(i,3).matrix = [frame_spout5(i,3) : frame_choice5(i)-1]; %Spout move to choice
    state_time(i,4).matrix = [frame_choice5(i) : frame_choice5(i)+500/Time_stretch]; %Choice to Choice + 0.5 sec
    if i ~= length(frame_start)
        state_time(i,5).matrix = [frame_choice5(i)+1+500/Time_stretch : frame_spout5(i+1,1)-1]; %Choice + 0.5sec to ITI
    else
        state_time(i,5).matrix = [frame_choice5(i)+1+500/Time_stretch : frame_choice5(i)+1000/Time_stretch]; %+0.5sec
    end
    
    for j = 1:5
        if max(isnan(state_time(i,j).matrix)) ~= 1
            frame5_state(j,state_time(i,j).matrix) = 1; 
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sound kernels:
%frame5_low(frame_sound5(low)) = 1;
%frame5_high(frame_sound5(high)) = 1;
% %Sound off
% frame5_off = false(1,frame5_time);
% frame5_off(frame_sound_off5) = 1;
% % figure
% % plot(frame5_low,'b')
% % hold on
% % plot(frame5_high,'r')

%Sound for tone clouds
frame5_evi = false(length(tone_evidence),frame5_time);
frame5_evi_long = false(length(tone_evidence),frame5_time);
frame5_evi_short = false(length(tone_evidence),frame5_time);
frame5_evi_off = false(length(tone_evidence),frame5_time);
for i = 1:length(tone_evidence)
    frame5_evi(i,frame_sound5(evi_trial(i).matrix)) = 1;
    %long only
    temp_trial = intersect(evi_trial(i).matrix, long);
    frame5_evi_long(i,frame_sound5(temp_trial)) = 1;
    frame5_evi_off(i,frame_sound5(temp_trial) + duration_long_sound/Time_stretch) = 1;
    %short only
    temp_trial = intersect(evi_trial(i).matrix, short);
    frame5_evi_short(i,frame_sound5(temp_trial)) = 1;
    frame5_evi_off(i,frame_sound5(temp_trial) + duration_short_sound/Time_stretch) = 1;
end
% use_color = jet(length(tone_evidence));
% figure
% for i = 1:length(tone_evidence)
%     plot(frame5_evi(i,:),'color',use_color(i,:));
%     hold on
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sound duration kernel, Make temporal different windows
%Long stimulus: -0.2 sec to 1.5 sec
time1 = round(-200 / Time_stretch);
time2 = round(1500 / Time_stretch);
time_count = [time1:1:time2];
frame5_sound_evi_long = logical([]);
for i = 1:length(tone_evidence)
    %frame_sound_on_evi(i).matrix = make_time_kernel(frame5_evi(i,:), time_count);
    temp = make_time_kernel(frame5_evi_long(i,:), time_count);
    frame5_sound_evi_long = [frame5_sound_evi_long; temp];
end

%Short stimulus: -0.2 sec to 0.7 sec
time1 = round(-200 / Time_stretch);
time2 = round(700 / Time_stretch);
time_count = [time1:1:time2];
frame5_sound_evi_short = logical([]);
for i = 1:length(tone_evidence)
    %frame_sound_on_evi(i).matrix = make_time_kernel(frame5_evi(i,:), time_count);
    temp = make_time_kernel(frame5_evi_short(i,:), time_count);
    frame5_sound_evi_short = [frame5_sound_evi_short; temp];
end

[frame5_sound_category_long,moto_frame_long] = get_frame_choice_kernel(lowL, highL, frame_sound5, frame5_time, Time_stretch, -200, 1500);
[frame5_sound_category_short,moto_frame_short] = get_frame_choice_kernel(lowS, highS, frame_sound5, frame5_time, Time_stretch, -200, 700);
moto_frame_sound = moto_frame_long + moto_frame_short;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Choice and Reward kernels
%-1.0 sec to 1.0 sec
frame5_leftC = false(1,frame5_time);
frame5_leftE = false(1,frame5_time);
frame5_rightC = false(1,frame5_time);
frame5_rightE = false(1,frame5_time);
frame5_leftC(frame_choice5(leftC)) = 1;
frame5_leftE(frame_choice5(leftE)) = 1;
frame5_rightC(frame_choice5(rightC)) = 1;
frame5_rightE(frame_choice5(rightE)) = 1;

time1 = round(-100 / Time_stretch);
time2 = round(1000 / Time_stretch);
time_count = [time1:1:time2];
frame5_leftC_all = make_time_kernel(frame5_leftC, time_count);
frame5_leftE_all = make_time_kernel(frame5_leftE, time_count);
frame5_rightC_all = make_time_kernel(frame5_rightC, time_count);
frame5_rightE_all = make_time_kernel(frame5_rightE, time_count);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make choice kernels during sound
frame5_choice_long = get_frame_choice_kernel(leftL, rightL, frame_sound5, frame5_time, Time_stretch, -200, 1500);
frame5_choice_short = get_frame_choice_kernel(leftS, rightS, frame_sound5, frame5_time, Time_stretch, -200, 700);

%Make prior kernels during sound
frame5_prior_long = get_frame_choice_kernel(prior0L, prior1L, frame_sound5, frame5_time, Time_stretch, -200, 1500);
frame5_prior_short = get_frame_choice_kernel(prior0S, prior1S, frame_sound5, frame5_time, Time_stretch, -200, 700);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Make the task state dependent kernels
% frame5_state_leftC = get_state_behave_kernel(leftC, state_time, frame5_time);
% frame5_state_leftE = get_state_behave_kernel(leftE, state_time, frame5_time);
% frame5_state_rightC = get_state_behave_kernel(rightC, state_time, frame5_time);
% frame5_state_rightE = get_state_behave_kernel(rightE, state_time, frame5_time);
% %Make the task state dependent kernels based on previous trial
% frame5_state_leftC_next = get_state_behave_kernel_next_trial(leftC, state_time, frame5_time);
% frame5_state_leftE_next = get_state_behave_kernel_next_trial(leftE, state_time, frame5_time);
% frame5_state_rightC_next = get_state_behave_kernel_next_trial(rightC, state_time, frame5_time);
% frame5_state_rightE_next = get_state_behave_kernel_next_trial(rightE, state_time, frame5_time);
% % figure
% % imagesc(frame5_state_leftE_next)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Prior state kernels
% frame5_state_prior = get_state_behave_kernel_01(block2, block3, state_time, frame5_time); %6 different state with prior
% random_kernel = true(1,frame5_time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Kernel from RL model: Prior, Posterior, Q_choice
%state_time is for all trials!!
%frame5_prior_RL = get_state_behave_kernel_RL_model(prior, use_trial, state_time, frame5_time);
%frame5_posterior_RL = get_state_behave_kernel_RL_model(posterior, use_trial, state_time, frame5_time);
frame5_Q_choice_RL = get_state_behave_kernel_RL_model(Q_choice, use_trial, state_time, frame5_time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Average speed of mice
mean_speed = mean(ave_velocity);
speed_new = zeros(1,length(ave_velocity(1:length(ave_velocity)-Time_stretch+1)));
for i = 1:Time_stretch
    speed_new = speed_new + ave_velocity(i:length(ave_velocity)-Time_stretch+i);
end
speed_new = speed_new(1:Time_stretch:length(speed_new));
if length(speed_new) >= frame5_time
    frame5_speed = speed_new(1:frame5_time);
else
    frame5_speed = ones(1,frame5_time) .* mean_speed;
    frame5_speed(1:length(speed_new)) = speed_new;
end

%Kernel all
%frame5_state frame5_state_prior
%frame_sound_on_evi frame_sound_on_evi_long frame_sound_off_evi
%frame5_leftC_all frame5_leftE_all frame5_rightC_all frame5_rightE_all
%frame5_state_leftC frame5_state_leftE frame5_state_rightC frame5_state_rightE
%frame5_state_leftC_next frame5_state_leftE_next frame5_state_rightC_next frame5_state_rightE_next

% kernel_all = [frame5_state; random_kernel];

%frame state is included in the frame_state_leftC...
% kernel_all = [frame5_sound_evi_long; frame5_sound_evi_short; ...
%               frame5_leftC_all; frame5_leftE_all; frame5_rightC_all; frame5_rightE_all; ...
%               frame5_state_leftC; frame5_state_leftE; frame5_state_rightC; frame5_state_rightE; ...
%               frame5_state_leftC_next; frame5_state_leftE_next; frame5_state_rightC_next; frame5_state_rightE_next; ...
%               frame5_speed; frame5_prior_RL; frame5_Q_choice_RL];
%frame state is included in the frame_state_leftC...
kernel_all = [frame5_sound_category_long; frame5_sound_category_short; ...
              frame5_choice_long; frame5_choice_short; ...
              frame5_prior_long; frame5_prior_short; ...
              frame5_leftC_all; frame5_leftE_all; frame5_rightC_all; frame5_rightE_all; ...
              frame5_speed];
%islogical(kernel_all)

kernel_all = double(kernel_all);

kernel_size_y(1)  = size(frame5_sound_category_long,1);
kernel_size_y(2)  = size(frame5_sound_category_short,1);
kernel_size_y(3)  = size(frame5_choice_long,1);
kernel_size_y(4)  = size(frame5_choice_short,1);
kernel_size_y(5)  = size(frame5_prior_long,1);
kernel_size_y(6)  = size(frame5_prior_short,1);
kernel_size_y(7)  = size(frame5_leftC_all,1);
kernel_size_y(8)  = size(frame5_leftE_all,1);
kernel_size_y(9)  = size(frame5_rightC_all,1);
kernel_size_y(10)  = size(frame5_rightE_all,1);
kernel_size_y(11) = size(frame5_speed,1);
% kernel_size_y(7)  = size(frame5_state_leftC,1);
% kernel_size_y(8)  = size(frame5_state_leftE,1);
% kernel_size_y(9) = size(frame5_state_rightC,1);
% kernel_size_y(10) = size(frame5_state_rightE,1);
% kernel_size_y(11) = size(frame5_state_leftC_next,1);
% kernel_size_y(12) = size(frame5_state_leftE_next,1);
% kernel_size_y(13) = size(frame5_state_rightC_next,1);
% kernel_size_y(14) = size(frame5_state_rightE_next,1);
%kernel_size_y(16) = size(frame5_prior_RL,1);
%kernel_size_y(17) = size(frame5_Q_choice_RL,1);

kernel_size_y
sum(kernel_size_y)
size(kernel_all)
%hoge

%moto_frame5 = [frame5_evi_long; frame5_evi_short; frame5_leftC; frame5_leftE; frame5_rightC; frame5_rightE]; %6timings

%Decide the time window to analyze the data                           
%Use during the block 2 to 5
%use_trial
start_frame = state_time(use_trial(1),1).matrix; %use_trial(1) @ spout_start
start_frame = start_frame(1) - round(1000/Time_stretch); %1sec before trial start
end_frame = state_time(use_trial(end),5).matrix; %End time
end_frame = end_frame(end) + round(1000/Time_stretch); %1sec after trial start

kernel_all = kernel_all(:,[start_frame:end_frame]);
moto_frame_sound = moto_frame_sound([start_frame:end_frame]);
%moto_frame5 = moto_frame5(:,[start_frame:end_frame]);
size(kernel_all)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make vector for spike_density_function
gauss_std = 25 / Time_stretch; %25ms std

width = 1000 / Time_stretch; %1s to analze the gaussian filter
for i = 1:frame5_time
    use_time = [i : i + width];
    temp0 = find(use_time > 0 & use_time <= frame5_time);
    use_time = use_time(temp0);
    
    gauss_pdf = normpdf(use_time,i,gauss_std);
    gauss_pdf = gauss_pdf ./ sum(gauss_pdf);
    SDF_filter(i).time = use_time;
    SDF_filter(i).gauss_pdf = gauss_pdf;
end

CV_para = 10; %10 fold

%Open each spike data
% [pathname1] = uigetdir
% temp_cd = ['cd ',pathname1];
% eval(temp_cd);
cd(spike_dir)
mkdir('KernelEncoding_20220912');

%save KernelEncoding_20220320/KernelEncoding_20220320_kernel kernel_all kernel_size_y moto_frame5 CV_para
save KernelEncoding_20220912/KernelEncoding_20220912_kernel kernel_all kernel_size_y CV_para moto_frame_sound

tif_name = dir('task_spike*.mat'); %get all the tif files
length_tif = length(tif_name);
max_tif = length_tif
parfor file_count = 1:max_tif
%for file_count = 1:max_tif
    %temp_file = tif_name(file_count).name;
    %temp_file
    temp_file = sprintf('task_spike_stripe20210520_%d',file_count);
    [file_count, max_tif]
    temp_file
    
    data = load(temp_file); %spike_mark
    spike_count = data.spike_mark;
    
    %Make half Gaussin with the time window of 5ms.
    spike_new = zeros(1,length(spike_count(1:length(spike_count)-Time_stretch+1)));
    for i = 1:Time_stretch
        spike_new = spike_new + spike_count(i:length(spike_count)-Time_stretch+i);
    end
    spike_new = spike_new(1:Time_stretch:length(spike_new));
    spike_new(frame5_time) = 0; %Make sure that the length is enough.
    
    %Based on the spike_mark, make spike density function
    spike_filter_new = zeros(1,frame5_time); %ms
    for j = 1:frame5_time
        spike_filter_new(j) = sum(spike_new(SDF_filter(j).time) .* SDF_filter(j).gauss_pdf);
    end
    
    spike_filter_new = spike_filter_new(start_frame:end_frame)';
    %lasso_left = SparseLogisticRegression_CV_20210620(neuron_block(j).matrix(left,:), decode_sound(left), CV_para);
    
    %normalized the spike_filter_new
    mean_spike = mean(spike_filter_new);
    std_spike = std(spike_filter_new);
    if std_spike ~= 0
        spike_filter_new = (spike_filter_new - mean_spike) ./ std_spike;
    else
        spike_filter_new = zeros(size(spike_filter_new)); %zero
    end
    
    temp_test = unique(spike_filter_new);
    unique_spike = length(temp_test);
    %if length(temp_test) == 1
    if unique_spike <= CV_para*2
        temp_test
        Ridge_neuron = [];
    else
    %%%%%%%%%%%%%%%
    %Because the spike rate is filtered, it is not already poisson
    %distribution
    %Instead, just use the log(spike) and deal as Gaussian function
    %Ridge_neuron = EncodingKernel_CV_Ridge_20211216(spike_filter_new, kernel_all, CV_para);
    Ridge_neuron = EncodingKernel_Ridge_20220215_glmnet(spike_filter_new, kernel_all, CV_para);
    %Ridge_neuron = EncodingKernel_CV_Ridge_20211215_glmnet_poisson(spike_filter_new, kernel_all, CV_para);
    end    
    save_file = sprintf('KernelEncoding_20220912/KernelEncoding_20220912_glmnet_%d',file_count);
    
    save_parfor(save_file,Ridge_neuron,spike_filter_new,CV_para,Time_stretch,kernel_size_y,unique_spike);
end
delete(gcp('nocreate'))

cd(pathname)
%hoge

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [frame5_choice_long,moto_frame5] = get_frame_choice_kernel(leftL, rightL, frame_sound5, frame5_time, Time_stretch, pre_sound, post_sound)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frame5_left_long = false(1,frame5_time);
frame5_left_long(frame_sound5(leftL)) = 1;
frame5_right_long = false(1,frame5_time);
frame5_right_long(frame_sound5(rightL)) = 1;

time1 = round(pre_sound / Time_stretch);
time2 = round(post_sound / Time_stretch);
time_count = [time1:1:time2];
temp1 = make_time_kernel(frame5_left_long, time_count);
temp2 = make_time_kernel(frame5_right_long, time_count);
frame5_choice_long = [temp1; temp2];

moto_frame5 = frame5_right_long + frame5_left_long;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function save_parfor(save_file,Ridge_neuron,spike_filter_new,CV_para,Time_stretch,kernel_size_y,unique_spike)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save(save_file,'Ridge_neuron','spike_filter_new','CV_para','Time_stretch','kernel_size_y','unique_spike') 
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function frame5_state_leftC = get_state_behave_kernel_next_trial(leftC, state_time, frame5_time)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frame5_state_leftC = false(6,frame5_time);
[size_y,size_x] = size(state_time);
for i = 1:length(leftC)
    if leftC(i)+1 <= size_y
        for j = 1:6
            temp = state_time(leftC(i)+1,j).matrix;
            if max(isnan(temp)) ~= 1
                frame5_state_leftC(j,temp) = 1;
            end
        end
    end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function frame5_state = get_state_behave_kernel_01(left, right, state_time, frame5_time)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frame5_state = zeros(6,frame5_time);
for i = 1:length(left)
    for j = 1:6
        temp = state_time(left(i),j).matrix;
        frame5_state(j,temp) = 1;
    end
end
for i = 1:length(right)
    for j = 1:6
        temp = state_time(right(i),j).matrix;
        frame5_state(j,temp) = -1;
    end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function frame5_state = get_state_behave_kernel_RL_model(prior, use_trial, state_time, frame5_time)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(prior) ~= length(use_trial)
    [length(prior),length(use_trial),length(state_time)]
    hoge
end
frame5_state = ones(5,frame5_time);
frame5_state = frame5_state .* mean(prior);

for i = 1:length(prior)
    for j = 1:5
        temp = state_time(use_trial(i),j).matrix;
        frame5_state(j,temp) = prior(i);
    end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function frame5_state_leftC = get_state_behave_kernel(leftC, state_time, frame5_time)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frame5_state_leftC = false(6,frame5_time);
for i = 1:length(leftC)
    for j = 1:6
        temp = state_time(leftC(i),j).matrix;
        frame5_state_leftC(j,temp) = 1;
    end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function new_frame = make_time_kernel(use_frame, time_count)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
new_frame = false(length(time_count),length(use_frame));
temp_frame = find(use_frame == 1);
for i = 1:length(time_count)
    temp_time = time_count(i);
    temp = temp_frame + temp_time;
    temptemp = find(temp >= 1 & temp <= length(use_frame));
    temp = temp(temptemp);
    
    new_frame(i,temp) = 1;
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
