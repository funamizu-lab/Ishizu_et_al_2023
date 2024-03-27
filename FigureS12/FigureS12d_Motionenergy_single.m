function FigureS12d_Motionenergy_single

workpath='G:\Ishizu_data\movie\i34\20211119_AC';
cd(workpath);

close all;
n=4000;%ms / plot point (long-tone correct/error trial)

%% Get task parameter(Bpod)
temp = dir('Bpod*.mat');
Bpod_file=temp.name;
load(temp.name);

left  = find(Chosen_side == 0);
right = find(Chosen_side == 1);
stim_length = unique(StimDuration);
Long  = find(StimDuration == stim_length(2));
Short = find(StimDuration == stim_length(1));

Choice_trial = find(Outcome == 1 | Outcome == 2);
useBlock_trial= find(TrialBlock>1);
use_trial = intersect(Choice_trial,useBlock_trial);
use_trial(end) = [];

%%% Get chosen data
output.Choice = Chosen_side(use_trial);

%%% Get tone evidence in all trials form the true tone cloud value
binary_tone=zeros(length(Tone_cloud),1);
for i = 1:length(Tone_cloud)
    temp_tone = Tone_cloud(i).matrix;
    temp1 = find(temp_tone >= 9);  %Get the data in all sound
    binary_tone(i) = length(temp1) ./ length(temp_tone);
end

%Based on the correct trial, flip the tone cloud
if mean(binary_tone(Correct_side == 1)) < 0.5 %low for right correct
    binary_tone = 1 - binary_tone;    % flipping tones %
else
end
output.ToneES = binary_tone(use_trial);

LongCorrect= intersect(Long,find(Outcome==2));
LongError  = intersect(Long,find(Outcome==1));
longcor_usetrial = ismember(use_trial,LongCorrect);
longerr_usetrial = ismember(use_trial,LongError);

%% Get prior value form RL model
temp = dir('RL_20220818*'); 
load(temp.name,'para_max','N_trial');
Prior = Dual_RL_model_block1_20220314_para_determined(Bpod_file,para_max(3,:),N_trial); 
relativePrior=Prior(:,1)./(Prior(:,1)+Prior(:,2));
output.Prior=relativePrior(1:length(use_trial),:);

%% Get wheel speed(NI daq) 
temp = dir('task_frame_tokyo_ephys_20220210*');
load(temp.name,'frame_spout','frame_sound','frame_choice','frame_end','ave_velocity');

%frame_sound_off
frame_sound_off = frame_sound;
frame_sound_off(Long) = frame_sound_off(Long) + 1000;% Add 1000 ms
frame_sound_off(Short)= frame_sound_off(Short)+ round(stim_length(1)*10)*100; % Add short stim 200 or 400 ms

frame_sound = frame_sound(use_trial);
frame_sound_off = frame_sound_off(use_trial);

%frame_choice_select
frame_choice_select = nan(length(frame_choice),1);
frame_choice_select(left) = frame_choice(left,1);
frame_choice_select(right)= frame_choice(right,2);

frame_choice_select =frame_choice_select(use_trial);
frame_spout=frame_spout(use_trial,:);

% frame end
frame_end = frame_end(use_trial);
    
if(sum(ave_velocity)==0)
else
    %normalize ave_velocity%
    ave_velocity = rescale(ave_velocity);
    %Get wheel speed value
    speed_prior = nan(length(use_trial),1);
    speed_sound = nan(length(use_trial),1);
    speed_choice= nan(length(use_trial),1);
    speedtime_longcor_trial  = cell(length(find(longcor_usetrial)),1); j=1;
    speedtime_longerr_trial  = cell(length(find(longerr_usetrial)),1); k=1;
    for i=1:length(use_trial)
        time_trial = frame_spout(i,1):frame_end(i);
        time_prior = frame_spout(i,1):frame_sound(i);
        time_sound = frame_sound(i):frame_sound_off(i);
        time_choice= frame_sound_off(i):frame_choice_select(i);
        
        speed_prior(i) = mean(ave_velocity(time_prior));
        speed_sound(i) = mean(ave_velocity(time_sound));
        speed_choice(i)= mean(ave_velocity(time_choice));
        tmp=ave_velocity(time_trial);
        if(longcor_usetrial(i))
            speedtime_longcor_trial{j} = tmp(1:n);
            j=j+1;
        end
        if(longerr_usetrial(i))
            speedtime_longerr_trial{k} = tmp(1:n);
            k=k+1;
        end
    end
end


%% Get motion energy (DLC data) 
temp = dir('Motion*');
load(temp.name,'video_t','trialID','spout','mouth_MEraw','nose_MEraw','tongue_MEraw');
allME = mouth_MEraw+nose_MEraw+tongue_MEraw;
allME = rescale(allME);

% sychronize data between taskframe & DLC 
FirstUseTrial = use_trial(1);
FisrtSpoutAwayTime_taskframe = frame_spout(1,1);

tmp=find(spout>0.5);
tmp2=find(diff(tmp)>10);% the timing that spout go across 0.5
if(tmp(1)==1)% spout away
    spoutAwayID = tmp(tmp2); 
else
    spoutAwayID = tmp(tmp2+1); 
end
tmp=find(trialID==FirstUseTrial-1);
FisrtSpoutAwayTime_video = video_t(tmp(ismember(tmp,spoutAwayID)));
video_t = video_t-FisrtSpoutAwayTime_video;

% ajust task frame time
frame_sound = frame_sound - FisrtSpoutAwayTime_taskframe;
frame_sound_off = frame_sound_off - FisrtSpoutAwayTime_taskframe;
frame_choice_select =frame_choice_select - FisrtSpoutAwayTime_taskframe;
frame_spout = frame_spout - FisrtSpoutAwayTime_taskframe;
frame_end = frame_end - FisrtSpoutAwayTime_taskframe;

%Get motion value
ME_prior  = nan(length(use_trial),1);
ME_sound  = nan(length(use_trial),1);
ME_choice = nan(length(use_trial),1);
ME_longcor_trial = cell(length(find(longcor_usetrial)),1);
ME_longerr_trial = cell(length(find(longerr_usetrial)),1);
j=1; k=1;
for i=1:length(use_trial)
    time_trial = frame_spout(i,1):frame_end(i);
    time_prior = frame_spout(i,1):frame_sound(i);
    time_sound = frame_sound(i):frame_sound_off(i);
    time_choice= frame_sound_off(i):frame_choice_select(i);
    
    Frame_trial = ismember(video_t,time_trial);
    Frame_prior = ismember(video_t,time_prior);
    Frame_sound = ismember(video_t,time_sound);
    Frame_choice= ismember(video_t,time_choice);   

    vTime = video_t(Frame_trial);
    
    ME_prior(i) = mean(allME(Frame_prior));
    ME_sound(i) = mean(allME(Frame_sound));
    ME_choice(i)= mean(allME(Frame_choice));
    if(longcor_usetrial(i))
        t=vTime -vTime(1);
        ME_longcor_trial{j} = spCorr(n,t,allME(Frame_trial));
        j=j+1;
    end
    if(longerr_usetrial(i))
        t=vTime -vTime(1);
        ME_longerr_trial{k} = spCorr(n,t,allME(Frame_trial));
        k=k+1;
    end
end
%% Each
allME_longcor= cell2mat(ME_longcor_trial);
allME_longerr= cell2mat(ME_longerr_trial);

%%% Fig S12d %%%
h=figure('Position',[10 10 1000 1000]);
subplot(1,2,1); hold on;
plot(1:n,allME_longcor(1,:),'k');
title('single correct trial')
xticks(0:500:n)

subplot(1,2,2); hold on;
p1=errorplot(1:n,mean(allME_longcor,1),std(allME_longcor,1),std(allME_longcor,1),'r',.5,1);
p2=errorplot(1:n,mean(allME_longerr,1),std(allME_longerr,1),std(allME_longerr,1),'k',.5,1);
legend([p1,p2],{'cor','err'});
title('trials')
xticks(0:500:n)

set(h,'PaperPositionMode','auto');
print(h,'-r0','move','-dpng');
print(h,'-r0','move','-dsvg');


%%% source data %%% 
cd('G:\upload_code\FigureS12\FigS12d');
sdata = struct();
sdata.xtime_msec = (1:n)'-1500;
sdata.y = allME_longcor(1,:)';
T = struct2table(sdata);
writetable(T, 'source fig S12d left.csv');


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  output= spCorr(n,t,v)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F = griddedInterpolant(t,v);
t_correct = 0:1:t(end);
v_correct = F(t_correct);
output= v_correct(1:n);
end
