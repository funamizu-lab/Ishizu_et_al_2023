
%{
----------------------------------------------------------------------------
norm_Task_kaiseki_tokyo1_20220316_sound_choice_separate2_detail
----------------------------------------------------------------------------
%}
function norm_Task_kaiseki_tokyo1_20220316_sound_presound_ishizu(folders)



analysis_dir = eval(folders);
for i = 1:length(analysis_dir)
    [i,length(analysis_dir)];
    pathname=analysis_dir{i};
    cd(pathname)
    %spike_dir = 'spike_ch1';
    spike_dir = dir('spike_ch*');
    if length(spike_dir) ~= 1
        hoge
    end
    spike_dir = spike_dir.name;
    
%     new_p_thre = 10;
    temp = dir('task_frame*');
    if length(temp) ~= 1
        hoge
    end
    load(temp.name);
    %frame_start
    %frame_sound
    %frame_end
    
    temp = dir('Bpod*');
    if length(temp) ~= 1
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
    
    %Get task parameter
    Choice_trial = find(Outcome == 1 | Outcome == 2);
    low  = find(Correct_side == 0);
    high = find(Correct_side == 1);
    stim_length = unique(StimDuration);
    Long  = find(StimDuration == stim_length(2));
    Short = find(StimDuration == stim_length(1));
    correct_error = Correct_side == Chosen_side;
    correct = find(correct_error == 1);
    
    
    %Adjust the trials to use
    correct = intersect(correct, use_trial);
    use_long = intersect(Long, use_trial);
    use_short = intersect(Short, use_trial);
%     if BlockReward(2,1) > BlockReward(2,2) %Left -> Right
%         Long_R = intersect(Long, block3);
%         Short_R = intersect(Short, block3);
%         Long_L = intersect(Long, block2);
%         Short_L = intersect(Short, block2);
%     else %Right -> Left
%         Long_R = intersect(Long, block2);
%         Short_R = intersect(Short, block2);
%         Long_L = intersect(Long, block3);
%         Short_L = intersect(Short, block3);
%     end
    
    time_window_base = 200;
    time_window = 100; %ms
    time_long_pre  = 1500; %4sec
    time_long_post = 2500; %4sec
    
    time_short_pre = 1500; %3.2sec
    time_short_post = 1700; %3.2sec
    
    cd(spike_dir);
    
    tif_name = dir('task_spike_stripe*.mat'); %get all the tif files
    length_tif = length(tif_name);
    max_tif = length_tif;
    %Long trial
    
    [long1_evi,long1_evi_correct,long1_evi_pre,long1_evi_correct_pre] = ...
        get_task_activ_separate_detail(use_long,frame_sound,max_tif,...
        time_long_pre,time_long_post,time_window, ...
        Correct_side, correct_error,trial_evidence);
    
% [short1_correct,short1_error,short1_stim,short1_choice,short1_prior, ...
%         short1_evi,short1_evi_correct,short1_evi_error,short1_evi_prior1,short1_evi_prior0,...
%         short1_evi_correct_prior1, short1_evi_correct_prior0,...
%         short2_correct,short2_error,short2_short,short2_choice,short2_prior, ...
%         short2_evi,short2_evi_correct,short2_evi_error,short2_evi_prior1,short2_evi_prior0,...
%         short2_evi_correct_prior1, short2_evi_correct_prior0] = ...
%         get_task_activ_separate_detail(use_short,frame_sound,frame_choice_select,frame_spout,max_tif,...
%         time_short_pre,time_short_post,time_short_pre2,time_short_post2,time_window,time_window_base, ...
%         Correct_side, Chosen_side, correct_error, all_binary_block, trial_evidence);
    
    cd(pathname)
    
    save norm_detail_20220316_sound_presound ...
        long1_evi long1_evi_correct long1_evi_pre long1_evi_correct_pre;
%     save norm_detail_20220316_sound_presound ...
%         long1_correct long1_error long1_stim long1_choice long1_prior ...
%         long1_evi long1_evi_correct long1_evi_error long1_evi_prior1 long1_evi_prior0 ...
%         long1_evi_correct_prior1 long1_evi_correct_prior0 ...
%         long2_correct long2_error long2_long long2_choice long2_prior ...
%         long2_evi long2_evi_correct long2_evi_error long2_evi_prior1 long2_evi_prior0 ...
%         long2_evi_correct_prior1 long2_evi_correct_prior0 ...
%         short1_correct short1_error short1_stim short1_choice short1_prior ...
%         short1_evi short1_evi_correct short1_evi_error short1_evi_prior1 short1_evi_prior0 ...
%         short1_evi_correct_prior1 short1_evi_correct_prior0 ...
%         short2_correct short2_error short2_short short2_choice short2_prior ...
%         short2_evi short2_evi_correct short2_evi_error short2_evi_prior1 short2_evi_prior0 ...
%         short2_evi_correct_prior1 short2_evi_correct_prior0;
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stim1_evi,stim1_correct,stim1_evi_pre,stim1_correct_pre] = ...
    get_task_activ_separate_detail(use_long,frame_sound,max_tif,...
    time_long_pre,time_long_post,time_window, ...
    Correct_side, correct_error, trial_evidence)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Correct_side = Correct_side(use_long);
% Chosen_side = Chosen_side(use_long);
correct_error = correct_error(use_long);
trial_evidence = trial_evidence(use_long);

tone_evidence = unique(trial_evidence);
if length(tone_evidence) ~= 6
    tone_evidence
end

low  = find(Correct_side == 0);
high = find(Correct_side == 1);
% left = find(Chosen_side == 0);
% right = find(Chosen_side == 1);
correct = find(correct_error == 1);

trial_L_correct = intersect(low,correct);
trial_H_correct = intersect(high,correct);
for i = 1:length(tone_evidence)
    temp = find(trial_evidence == tone_evidence(i));
    evi_stim(i).matrix = temp;
    evi_correct(i).matrix = intersect(temp, correct);
    
    % previous trial evidence %
    temp2= temp-1;
    if(temp2(1)==0)
        previous_evidence = [-1;trial_evidence(temp2(2:end))];
    else
        previous_evidence = trial_evidence(temp2);
    end
    for j = 1:length(tone_evidence)
        temp3 = find(previous_evidence == tone_evidence(j));
        t=temp2(temp3);
        t(t==0)=[];
        evi_stim_pre(i,j).matrix = t;
        evi_correct_pre(i,j).matrix = intersect(t, correct);
    end
end

time_long  = round((time_long_pre + time_long_post) ./ time_window);
frame_sound = frame_sound(use_long,:);
for i = 1:time_long
    frame_sound_use_on(:,i) = frame_sound - time_long_pre + (i-1)*time_window;
    %frame_sound_use_off(:,i) = frame_sound - time_long_pre + i*time_window - 1;
end

[stim1_evi,stim1_correct,stim1_evi_pre,stim1_correct_pre] = ...
    get_basic_activity(max_tif, time_long, frame_sound_use_on, tone_evidence,...
    frame_sound,use_long,time_window,...
    trial_H_correct,trial_L_correct,high,low,evi_stim,evi_correct,evi_stim_pre,evi_correct_pre);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stim1_evi,stim1_correct,stim1_evi_pre,stim1_correct_pre] = ...
    get_basic_activity(max_tif, time_long, frame_sound_use_on, tone_evidence,...
    frame_sound, use_long,time_window,...
    trial_H_correct,trial_L_correct,high,low,evi_stim,evi_correct,evi_stim_pre,evi_correct_pre)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stim1_stim_correct1 = nan(max_tif,time_long);
stim1_stim_correct2 = nan(max_tif,time_long);
stim1_stim1 = nan(max_tif,time_long);
stim1_stim2 = nan(max_tif,time_long);

stim1_evi = nan(max_tif,length(tone_evidence),time_long);
stim1_correct = nan(max_tif,length(tone_evidence),time_long);
stim1_evi_pre = nan(max_tif,length(tone_evidence),length(tone_evidence),time_long);
stim1_correct_pre = nan(max_tif,length(tone_evidence),length(tone_evidence),time_long);

parfor file_count = 1:max_tif
    temp_file = sprintf('task_spike_stripe20210520_%d',file_count);
    
    %clear data
    data = load(temp_file); %spike_mark
    
    spike_mark = data.spike_mark;
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
        stim1_stim1(file_count,j) = mean(spike_frame1(high,j));
        stim1_stim2(file_count,j) = mean(spike_frame1(low,j));
        
        for k = 1:6
            stim1_evi(file_count,k,j) = mean(spike_frame1(evi_stim(k).matrix,j));
            stim1_correct(file_count,k,j) = mean(spike_frame1(evi_correct(k).matrix,j));
            for s = 1:6
                stim1_evi_pre(file_count,k,s,j) = mean(spike_frame1(evi_stim_pre(k,s).matrix,j));
                stim1_correct_pre(file_count,k,s,j) = mean(spike_frame1(evi_correct_pre(k,s).matrix,j));
            end
        end
    end
end
% stim_correct(1).matrix = stim1_stim_correct1;
% stim_correct(2).matrix = stim1_stim_correct2;
% stim1(1).matrix = stim1_stim1;
% stim1(2).matrix = stim1_stim2;

end

