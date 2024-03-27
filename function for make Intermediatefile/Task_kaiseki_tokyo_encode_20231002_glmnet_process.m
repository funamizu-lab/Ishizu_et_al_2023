
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
function Task_kaiseki_tokyo_encode_20231002_glmnet_process(pathname)
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

folder = '\KernelEncoding_20231002';
folder_sound = '\KernelEncoding_20231002_WO_sound';
folder_choice = '\KernelEncoding_20231002_WO_choice';
folder_prior = '\KernelEncoding_20231002_WO_prior';

folder = [pathname,folder];
folder_sound = [folder,folder_sound];
folder_choice = [folder,folder_choice];
folder_prior = [folder,folder_prior];

% %Open each spike data
% [pathname1] = uigetdir
% temp_cd = ['cd ',pathname1];
% eval(temp_cd); %Go to record folder

%Get each encoding file
cd(folder)
[tif_name,max_tif,parameter_size] = get_kernel_data;
cd(folder_sound)
[tif_sound,max_sound,parameter_sound] = get_kernel_data;
cd(folder_choice)
[tif_choice,max_choice,parameter_choice] = get_kernel_data;
cd(folder_prior)
[tif_prior,max_prior,parameter_prior] = get_kernel_data;

if max_tif ~= max_sound
    hoge
end
if max_tif ~= max_choice
    hoge
end
if max_tif ~= max_prior
    hoge
end

% kernel_size_y(1)  = size(frame5_sound_category_long,1);
% kernel_size_y(2)  = size(frame5_sound_category_short,1);
% kernel_size_y(3)  = size(frame5_choice_long,1);
% kernel_size_y(4)  = size(frame5_choice_short,1);
% kernel_size_y(5)  = size(frame5_prior_long,1);
% kernel_size_y(6)  = size(frame5_prior_short,1);
% kernel_size_y(7)  = size(frame5_leftC_all,1);
% kernel_size_y(8)  = size(frame5_leftE_all,1);
% kernel_size_y(9)  = size(frame5_rightC_all,1);
% kernel_size_y(10)  = size(frame5_rightE_all,1);
% kernel_size_y(11) = size(frame5_speed,1);

cd(pathname)

%parfor file_count = 1:max_tif
for file_count = 1:max_tif
    [file_count, max_tif]

    temp_file = tif_name(file_count).name;
    temp_file = [folder,'\',temp_file];
    data = load(temp_file);
    
    temp_file = tif_sound(file_count).name;
    temp_file = [folder_sound,'\',temp_file];
    data_sound = load(temp_file);

    temp_file = tif_choice(file_count).name;
    temp_file = [folder_choice,'\',temp_file];
    data_choice = load(temp_file);

    temp_file = tif_prior(file_count).name;
    temp_file = [folder_prior,'\',temp_file];
    data_prior = load(temp_file);

    [BIC(file_count,1),part_LL(file_count,1)] = get_BIC(data, parameter_size);
    [BIC(file_count,2),part_LL(file_count,2)] = get_BIC(data_sound, parameter_sound);
    [BIC(file_count,3),part_LL(file_count,3)] = get_BIC(data_choice, parameter_choice);
    [BIC(file_count,4),part_LL(file_count,4)] = get_BIC(data_prior, parameter_prior);
    
    %Get parameter from the all encoding model
    
    kernel_size_y = data.kernel_size_y;
    cum_y = cumsum(kernel_size_y);
    if length(data.Ridge_neuron) ~= 0
        p_REE(file_count,1) = signrank(data.Ridge_neuron.SumError, data_sound.Ridge_neuron.SumError);
        p_REE(file_count,2) = signrank(data.Ridge_neuron.SumError, data_choice.Ridge_neuron.SumError);
        p_REE(file_count,3) = signrank(data.Ridge_neuron.SumError, data_prior.Ridge_neuron.SumError);
        %median_error(file_count,:) = [median(data.Ridge_neuron.SumError), median(data_sound.Ridge_neuron.SumError)];
         median_error(file_count,:) = [median(data.Ridge_neuron.SumError), median(data_sound.Ridge_neuron.SumError), ...
                                       median(data_choice.Ridge_neuron.SumError), median(data_prior.Ridge_neuron.SumError)];
                                   
        x_all = data.Ridge_neuron.x_all;
        x_all(1) = []; %remove the constant term
        %Sound long parameter (-0.2sec to 1.5 sec)
        y_sound_l(file_count,:) = x_all([1:cum_y(1)]);
        y_sound_s(file_count,:) = x_all([cum_y(1)+1:cum_y(2)]);
        y_choice_l(file_count,:) = x_all([cum_y(2)+1:cum_y(3)]);
        y_choice_s(file_count,:) = x_all([cum_y(3)+1:cum_y(4)]);
        y_prior_l(file_count,:) = x_all([cum_y(4)+1:cum_y(5)]);
        y_prior_s(file_count,:) = x_all([cum_y(5)+1:cum_y(6)]);
    else
        % p_REE(file_count,:) = nan(1,1);
        % median_error(file_count,:) = nan(1,2);
        p_REE(file_count,:) = nan(1,3);
        median_error(file_count,:) = nan(1,4);
        
        y_sound_l(file_count,:) = nan(1,cum_y(1));
        y_sound_s(file_count,:) = nan(1,cum_y(2)-cum_y(1));
        y_choice_l(file_count,:) = nan(1,cum_y(3)-cum_y(2));
        y_choice_s(file_count,:) = nan(1,cum_y(4)-cum_y(3));
        y_prior_l(file_count,:) = nan(1,cum_y(5)-cum_y(4));
        y_prior_s(file_count,:) = nan(1,cum_y(6)-cum_y(5));
    end
end
delete(gcp('nocreate'))

%Save with all neurons
cd(pathname)
save glmnet_20231002_process BIC part_LL parameter_size parameter_sound parameter_choice parameter_prior ...
    y_sound_l y_sound_s y_choice_l y_choice_s y_prior_l y_prior_s
% 
% %nan check
% nanBIC = mean(BIC,2);
% nanBIC = find(isnan(nanBIC) == 0);
% BIC = BIC(nanBIC,:);
% part_LL = part_LL(nanBIC,:);
% 
% BIC_sound = BIC(:,2) - BIC(:,1);
% BIC_choice = BIC(:,3) - BIC(:,1);
% BIC_prior = BIC(:,4) - BIC(:,1);
% 
% LL_sound = part_LL(:,1) - part_LL(:,2);
% LL_choice = part_LL(:,1) - part_LL(:,3);
% LL_prior = part_LL(:,1) - part_LL(:,4);
% 
% figure
% subplot(1,2,1)
% plot(sort(LL_sound),'r')
% hold on
% plot(sort(LL_choice),'g')
% hold on
% plot(sort(LL_prior),'b')
% subplot(1,2,2)
% plot(sort(BIC_sound),'r')
% hold on
% plot(sort(BIC_choice),'g')
% hold on
% plot(sort(BIC_prior),'b')
% 
% LL_sound = 1-chi2cdf(2*LL_sound, parameter_size-parameter_sound);
% LL_sound = -log10(LL_sound);
% LL_choice = 1-chi2cdf(2*LL_choice, parameter_size-parameter_choice);
% LL_choice = -log10(LL_choice);
% LL_prior = 1-chi2cdf(2*LL_prior, parameter_size-parameter_prior);
% LL_prior = -log10(LL_prior);
% 
% figure
% plot(sort(LL_sound),'r')
% hold on
% plot(sort(LL_choice),'g')
% hold on
% plot(sort(LL_prior),'b')
% 
% % hoge
% % 
% % 
% % 
% % figure
% % plot(LL_sound,BIC_sound,'b.')
% % 
% % BIC_sound = sort(BIC_sound);
% % BIC_choice = sort(BIC_choice);
% % BIC_prior = sort(BIC_prior);
% % 
% % temp_x = [1:length(nanBIC)] ./ length(nanBIC);
% % 
% % figure
% % plot(BIC_sound, temp_x, 'r')
% % hold on
% % plot(BIC_choice, temp_x, 'g')
% % hold on
% % plot(BIC_prior, temp_x, 'b')
% % hold on
% % 
% % hoge

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tif_name,max_tif,parameter_size] = get_kernel_data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tif_name = dir('KernelEncoding_20231002_glmnet_*.mat'); %get all the tif files
max_tif = length(tif_name);
temp = dir('KernelEncoding_20231002_kernel*'); %get all the tif files
if length(temp) ~= 1
    temp
    hoge
end
data = load(temp.name);
kernel_size_y = data.kernel_size_y;
parameter_size = sum(kernel_size_y) + 1; %including constant

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [BIC,part_log_likeli] = get_BIC(data, parameter_size)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if length(data.Ridge_neuron) ~= 0
    SumError = data.Ridge_neuron.SumError;
    length_data = length(SumError);
    SumError = sum(SumError);
    BIC = length_data * log(SumError/length_data) + parameter_size * log(length_data);
    
    part_log_likeli = -0.5 * length_data * log(SumError/length_data);
 else
     BIC = nan;
     part_log_likeli = nan;
 end
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
frame5_state = ones(6,frame5_time);
frame5_state = frame5_state .* mean(prior);

for i = 1:length(prior)
    for j = 1:6
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
