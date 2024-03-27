
%{
----------------------------------------------------------------------------
norm_Task_kaiseki_tokyo1_20220316_sound_choice_separate2_detail
----------------------------------------------------------------------------
%}
function norm_Task_kaiseki_tokyo1_20220316_sound_choice_ishizu_TCchoice(folders)



analysis_dir = eval(folders);
for i = 1:length(analysis_dir)
    [i,length(analysis_dir)]
    pathname=analysis_dir{i};
    cd(pathname)
    %spike_dir = 'spike_ch1';
    spike_dir = dir('spike_ch*');
    if length(spike_dir) ~= 1
        hoge
    end
    spike_dir = spike_dir.name;
    
    new_p_thre = 10;
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
        = Dual_get_basic_task_structure_20210204_2(Bpod_file);
    
    %Instead of model prediction, use the block trials
    if BlockReward(2,1) > BlockReward(2,2) %Left -> Right
        BlockL = block2;
        BlockR = block3;
    else
        BlockL = block3;
        BlockR = block2;
    end
    [~,use_R] = intersect(use_trial,BlockR);
    [~,use_L] = intersect(use_trial,BlockL);
    all_binary_block = nan(length(Choice_trial),1);
    all_binary_block(use_trial(use_R)) = 1;
    all_binary_block(use_trial(use_L)) = 0;
    
    %Get task parameter
    Choice_trial = find(Outcome == 1 | Outcome == 2);
    low  = find(Correct_side == 0);
    high = find(Correct_side == 1);
    stim_length = unique(StimDuration);
    Long  = find(StimDuration == stim_length(2));
    Short = find(StimDuration == stim_length(1));
    left = find(Chosen_side == 0);
    right = find(Chosen_side == 1);
    correct_error = Correct_side == Chosen_side;
    correct = find(correct_error == 1);
    error = find(correct_error == 0);
    
    %Make the choice timing frames
    frame_choice_select = nan(length(frame_choice),1);
    frame_choice_select(left) = frame_choice(left,1);
    frame_choice_select(right) = frame_choice(right,2);
    
    %Adjust the trials to use
    correct = intersect(correct, use_trial);
    use_long = intersect(Long, use_trial);
    use_short = intersect(Short, use_trial);
    if BlockReward(2,1) > BlockReward(2,2) %Left -> Right
        Long_R = intersect(Long, block3);
        Short_R = intersect(Short, block3);
        Long_L = intersect(Long, block2);
        Short_L = intersect(Short, block2);
    else %Right -> Left
        Long_R = intersect(Long, block2);
        Short_R = intersect(Short, block2);
        Long_L = intersect(Long, block3);
        Short_L = intersect(Short, block3);
    end
    
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
    
    tif_name = dir('task_spike_stripe*.mat'); %get all the tif files
    length_tif = length(tif_name);
    max_tif = length_tif;
    %Long trial
    
    [long1_correct,long1_error,long1_stim,long1_choice,long1_prior, ...
        long1_evi,long1_evi_correct,long1_evi_error,long1_evi_prior1,long1_evi_prior0,...
        long1_evi_correct_prior1, long1_evi_correct_prior0,...
        long2_correct,long2_error,long2_long,long2_choice,long2_prior, ...
        long2_evi,long2_evi_correct,long2_evi_error,long2_evi_prior1,long2_evi_prior0,...
        long2_evi_correct_prior1, long2_evi_correct_prior0] = ...
        get_task_activ_separate_detail(use_long,frame_sound,frame_choice_select,frame_spout,max_tif,...
        time_long_pre,time_long_post,time_long_pre2,time_long_post2,time_window,time_window_base, ...
        Correct_side, Chosen_side, correct_error, all_binary_block, trial_evidence);
    
    [short1_correct,short1_error,short1_stim,short1_choice,short1_prior, ...
        short1_evi,short1_evi_correct,short1_evi_error,short1_evi_prior1,short1_evi_prior0,...
        short1_evi_correct_prior1, short1_evi_correct_prior0,...
        short2_correct,short2_error,short2_short,short2_choice,short2_prior, ...
        short2_evi,short2_evi_correct,short2_evi_error,short2_evi_prior1,short2_evi_prior0,...
        short2_evi_correct_prior1, short2_evi_correct_prior0] = ...
        get_task_activ_separate_detail(use_short,frame_sound,frame_choice_select,frame_spout,max_tif,...
        time_short_pre,time_short_post,time_short_pre2,time_short_post2,time_window,time_window_base, ...
        Correct_side, Chosen_side, correct_error, all_binary_block, trial_evidence);
    
    cd(pathname)
    
    save norm3_Tokyo2_detail_20220316_sound_choice_separate ...
        long1_correct long1_error long1_stim long1_choice long1_prior ...
        long1_evi long1_evi_correct long1_evi_error long1_evi_prior1 long1_evi_prior0 ...
        long1_evi_correct_prior1 long1_evi_correct_prior0 ...
        long2_correct long2_error long2_long long2_choice long2_prior ...
        long2_evi long2_evi_correct long2_evi_error long2_evi_prior1 long2_evi_prior0 ...
        long2_evi_correct_prior1 long2_evi_correct_prior0 ...
        short1_correct short1_error short1_stim short1_choice short1_prior ...
        short1_evi short1_evi_correct short1_evi_error short1_evi_prior1 short1_evi_prior0 ...
        short1_evi_correct_prior1 short1_evi_correct_prior0 ...
        short2_correct short2_error short2_short short2_choice short2_prior ...
        short2_evi short2_evi_correct short2_evi_error short2_evi_prior1 short2_evi_prior0 ...
        short2_evi_correct_prior1 short2_evi_correct_prior0;
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stim1_correct,stim1_error,stim1_stim,stim1_choice,stim1_prior, ...
    stim1_evi,stim1_evi_correct,stim1_evi_error,stim1_evi_prior1,stim1_evi_prior0,...
    stim1_evi_correct_prior1, stim1_evi_correct_prior0,...
    stim2_correct,stim2_error,stim2_stim,stim2_choice,stim2_prior, ...
    stim2_evi,stim2_evi_correct,stim2_evi_error,stim2_evi_prior1,stim2_evi_prior0,...
    stim2_evi_correct_prior1, stim2_evi_correct_prior0] = ...
    get_task_activ_separate_detail(use_long,frame_sound,frame_choice_select,frame_spout,max_tif,...
    time_long_pre,time_long_post,time_long_pre2,time_long_post2,time_window,time_window_base, ...
    Correct_side, Chosen_side, correct_error, binary_prior, trial_evidence)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Correct_side = Correct_side(use_long);
Chosen_side = Chosen_side(use_long);
correct_error = correct_error(use_long);
binary_prior = binary_prior(use_long);
trial_evidence = trial_evidence(use_long);

tone_evidence = unique(trial_evidence);
if length(tone_evidence) ~= 6
    tone_evidence
end

low  = find(Correct_side == 0);
high = find(Correct_side == 1);
left = find(Chosen_side == 0);
right = find(Chosen_side == 1);
correct = find(correct_error == 1);
error = find(correct_error == 0);
prior1 = find(binary_prior == 1);
prior0 = find(binary_prior == 0);

trial_L_correct = intersect(low,correct);
trial_L_error = intersect(low,error);
trial_H_correct = intersect(high,correct);
trial_H_error = intersect(high,error);

for i = 1:length(tone_evidence)
    temp = find(trial_evidence == tone_evidence(i));
    evi_stim(i).matrix = temp;
    %Only with correct trials
    evi_correct(i).matrix = intersect(temp, correct);
    evi_error(i).matrix = intersect(temp, error);
    evi_prior1(i).matrix = intersect(temp, prior1);
    evi_prior0(i).matrix = intersect(temp, prior0);
    
    evi_correct_prior1(i).matrix = intersect(evi_correct(i).matrix, prior1);
    evi_correct_prior0(i).matrix = intersect(evi_correct(i).matrix, prior0);
end

time_long  = round((time_long_pre + time_long_post) ./ time_window);
time_long2 = round((time_long_pre2 + time_long_post2) ./ time_window);

frame_sound = frame_sound(use_long,:);
frame_choice_select = frame_choice_select(use_long,:);

frame_spout = frame_spout(use_long,:);
base_for_spout_on  = frame_spout(:,1)-time_window_base; %Before moving spout
base_for_spout_off = frame_spout(:,3)-time_window_base; %Before moving spout

for i = 1:time_long
    frame_sound_use_on(:,i) = frame_sound - time_long_pre + (i-1)*time_window;
    %frame_sound_use_off(:,i) = frame_sound - time_long_pre + i*time_window - 1;
end
for i = 1:time_long2
    frame_sound_use2_on(:,i) = frame_choice_select - time_long_pre2 + (i-1)*time_window;
    %frame_sound_use2_off(:,i) = frame_choice_select - time_long_pre2 + i*time_window - 1;
end

[stim1_correct,stim1_error,stim1_stim,stim1_choice,stim1_prior, ...
    stim1_evi,stim1_evi_correct,stim1_evi_error,stim1_evi_prior1,stim1_evi_prior0,...
    stim1_evi_correct_prior1, stim1_evi_correct_prior0] = ...
    get_basic_activity(max_tif, time_long, frame_sound_use_on, tone_evidence, frame_sound,use_long,time_window,...
    trial_H_correct,trial_L_correct,trial_H_error,trial_L_error,high,low,right,left,prior1,prior0, ...
    evi_stim,evi_correct,evi_error,evi_prior1,evi_prior0,evi_correct_prior1,evi_correct_prior0);

[stim2_correct,stim2_error,stim2_stim,stim2_choice,stim2_prior, ...
    stim2_evi,stim2_evi_correct,stim2_evi_error,stim2_evi_prior1,stim2_evi_prior0,...
    stim2_evi_correct_prior1, stim2_evi_correct_prior0] = ...
    get_basic_activity(max_tif, time_long2, frame_sound_use2_on, tone_evidence, frame_sound,use_long,time_window,...
    trial_H_correct,trial_L_correct,trial_H_error,trial_L_error,high,low,right,left,prior1,prior0, ...
    evi_stim,evi_correct,evi_error,evi_prior1,evi_prior0,evi_correct_prior1,evi_correct_prior0);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stim_correct,stim_error,stim1,choice1,stim1_prior, ...
    stim1_evi,stim1_correct,stim1_error,stim1_evi_prior1,stim1_evi_prior0,stim1_correct_prior1, stim1_correct_prior0] = ...
    get_basic_activity(max_tif, time_long, frame_sound_use_on, tone_evidence, frame_sound, use_long,time_window,...
    trial_H_correct,trial_L_correct,trial_H_error,trial_L_error,high,low,right,left,prior1,prior0, ...
    evi_stim,evi_correct,evi_error,evi_prior1,evi_prior0,evi_correct_prior1,evi_correct_prior0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stim1_stim_correct1 = nan(max_tif,time_long);
stim1_stim_correct2 = nan(max_tif,time_long);
stim1_stim_error1 = nan(max_tif,time_long);
stim1_stim_error2 = nan(max_tif,time_long);
stim1_stim1 = nan(max_tif,time_long);
stim1_stim2 = nan(max_tif,time_long);
stim1_choice1 = nan(max_tif,time_long);
stim1_choice2 = nan(max_tif,time_long);
stim1_prior1 = nan(max_tif,time_long);
stim1_prior2 = nan(max_tif,time_long);

stim1_evi = nan(max_tif,length(tone_evidence),time_long);
stim1_correct = nan(max_tif,length(tone_evidence),time_long);
stim1_error = nan(max_tif,length(tone_evidence),time_long);
stim1_evi_prior1 = nan(max_tif,length(tone_evidence),time_long);
stim1_evi_prior0 = nan(max_tif,length(tone_evidence),time_long);
stim1_correct_prior1 = nan(max_tif,length(tone_evidence),time_long);
stim1_correct_prior0 = nan(max_tif,length(tone_evidence),time_long);

parfor file_count = 1:max_tif
    %temp_file = tif_name(file_count).name;
    %temp_file
    temp_file = sprintf('task_spike_stripe20210520_%d',file_count);
    temp_file
    
    %clear data
    data = load(temp_file); %spike_mark
    
    spike_mark = data.spike_mark;
    spike_multi(file_count) = data.single_multi;
    spike_frame1 = nan(length(frame_sound),time_long);
    
    %Normalize the spike_mark
    mean_spike = mean(spike_mark);
    std_spike = std(spike_mark);
    if std_spike ~= 0
        spike_mark = (spike_mark - mean_spike) ./ std_spike; %normalized the activity
    else
        spike_mark = zeros(1,length(spike_mark));
    end
    
    for i = 1:length(use_long) %trial
        for j = 1:time_long
            if isnan(frame_sound_use_on(i,j)) == 0
                temp_frame = [frame_sound_use_on(i,j) : frame_sound_use_on(i,j)+time_window-1];
                spike_frame1(i,j) = mean(spike_mark(temp_frame));
            end
        end
    end
    for j = 1:time_long
        stim1_stim_correct1(file_count,j) = mean(spike_frame1(trial_H_correct,j));
        stim1_stim_correct2(file_count,j) = mean(spike_frame1(trial_L_correct,j));
        stim1_stim_error1(file_count,j) = mean(spike_frame1(trial_H_error,j));
        stim1_stim_error2(file_count,j) = mean(spike_frame1(trial_L_error,j));
        stim1_stim1(file_count,j) = mean(spike_frame1(high,j));
        stim1_stim2(file_count,j) = mean(spike_frame1(low,j));
        stim1_choice1(file_count,j) = mean(spike_frame1(right,j));
        stim1_choice2(file_count,j) = mean(spike_frame1(left,j));
        stim1_prior1(file_count,j) = mean(spike_frame1(prior1,j));
        stim1_prior2(file_count,j) = mean(spike_frame1(prior0,j));
        
        for k = 1:6
            stim1_evi(file_count,k,j) = mean(spike_frame1(evi_stim(k).matrix,j));
            stim1_correct(file_count,k,j) = mean(spike_frame1(evi_correct(k).matrix,j));
            stim1_error(file_count,k,j) = mean(spike_frame1(evi_error(k).matrix,j));
            stim1_evi_prior1(file_count,k,j) = mean(spike_frame1(evi_prior1(k).matrix,j));
            stim1_evi_prior0(file_count,k,j) = mean(spike_frame1(evi_prior0(k).matrix,j));
            stim1_correct_prior1(file_count,k,j) = mean(spike_frame1(evi_correct_prior1(k).matrix,j));
            stim1_correct_prior0(file_count,k,j) = mean(spike_frame1(evi_correct_prior0(k).matrix,j));
        end
    end
end
stim_correct(1).matrix = stim1_stim_correct1;
stim_correct(2).matrix = stim1_stim_correct2;
stim_error(1).matrix = stim1_stim_error1;
stim_error(2).matrix = stim1_stim_error2;
stim1(1).matrix = stim1_stim1;
stim1(2).matrix = stim1_stim2;
choice1(1).matrix = stim1_choice1;
choice1(2).matrix = stim1_choice2;
stim1_prior(1).matrix = stim1_prior1;
stim1_prior(2).matrix = stim1_prior2;

end

